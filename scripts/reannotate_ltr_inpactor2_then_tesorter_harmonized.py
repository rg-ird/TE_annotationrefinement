#!/usr/bin/env python3

import argparse
import subprocess
import re
from pathlib import Path
from collections import defaultdict


def clean_repeat_id(seqid):
    """Remove RepeatMasker/FASTA classification suffix after #."""
    return seqid.split("#")[0]


def parse_fasta_headers(fasta):
    headers = {}
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                raw = line[1:].strip()
                seq_id = raw.split()[0]
                headers[seq_id] = raw
    return headers


def get_class_from_inpactor2_header(header):
    """
    Inpactor2 header examples:
      #LTR/RLG/TAT   -> LTR/Gypsy/TAT
      #LTR/RLC/SIRE  -> LTR/Copia/SIRE
    """
    m = re.search(r"#?LTR/(RLG|RLC)/([A-Za-z0-9_.-]+)", header)
    if not m:
        return None

    code = m.group(1)
    family = normalize_family_name(m.group(2))

    if code == "RLG":
        return f"LTR/Gypsy/{family}"
    if code == "RLC":
        return f"LTR/Copia/{family}"
    return None


def normalize_family_name(family):
    """
    Harmonize family names from Inpactor2 / REXdb / TEsorter.

    Target nomenclature used in the final RepeatMasker .out:
      Gypsy: TEKAY-DEL, GALADRIEL, REINA, CRM, ATHILA, TAT
      Copia: TAR-TORK, BIANCA, SIRE, ALE-RETROFIT, ANGELA, IVANA-ORYCO

    Compatibility rules requested:
      Tekay -> TEKAY-DEL
      Galadriel -> GALADRIEL
      Reina -> REINA
      Athila -> ATHILA
      TatI/TatII/TatIII/Ogre/Retand -> TAT
      Ikeros/Tork/TAR -> TAR-TORK
      Ale -> ALE-RETROFIT
      Ivana -> IVANA-ORYCO
    """
    if family is None:
        return None

    fam = str(family).strip()
    fam = fam.replace(" ", "_")
    fam_clean = fam.strip(";:,[](){}")
    key = fam_clean.lower().replace("_", "-")

    aliases = {
        # Gypsy
        "tekay": "TEKAY-DEL",
        "tekay-del": "TEKAY-DEL",
        "galadriel": "GALADRIEL",
        "reina": "REINA",
        "crm": "CRM",
        "athila": "ATHILA",

        # TAT group harmonization
        "tat": "TAT",
        "tati": "TAT",
        "tat-i": "TAT",
        "tatii": "TAT",
        "tat-ii": "TAT",
        "tatiii": "TAT",
        "tat-iii": "TAT",
        "ogre": "TAT",
        "retand": "TAT",

        # Copia
        "tar-tork": "TAR-TORK",
        "tar": "TAR-TORK",
        "tork": "TAR-TORK",
        "ikeros": "TAR-TORK",
        "bianca": "BIANCA",
        "sire": "SIRE",
        "ale": "ALE-RETROFIT",
        "ale-retrofit": "ALE-RETROFIT",
        "angela": "ANGELA",
        "ivana": "IVANA-ORYCO",
        "ivana-oryco": "IVANA-ORYCO",

        # Keep mixture as a special non-family category if present
        "mixture": "mixture",
    }

    return aliases.get(key, fam_clean)


def is_ltr_family_annotation(annotation):
    if not annotation:
        return False
    low = annotation.lower()
    if low.endswith("/unknown") or low.endswith("/unclassified"):
        return False
    return annotation.startswith("LTR/Gypsy/") or annotation.startswith("LTR/Copia/")


def is_ltr_to_reannotate(annotation):
    return annotation in {
        "LTR/Gypsy",
        "LTR/Copia",
        "LTR/Gypsy/Unknown",
        "LTR/Copia/Unknown",
    }


def unknown_ltr_annotation(annotation):
    if annotation in {"LTR/Gypsy", "LTR/Gypsy/Unknown"}:
        return "LTR/Gypsy/Unknown"
    if annotation in {"LTR/Copia", "LTR/Copia/Unknown"}:
        return "LTR/Copia/Unknown"
    return annotation


def order_from_annotation(annotation):
    if annotation.startswith("LTR/Gypsy"):
        return "Gypsy"
    if annotation.startswith("LTR/Copia"):
        return "Copia"
    return None


def run_blast(hite_fa, inpactor_fa, out_tsv, threads, reuse_blast=False):
    if reuse_blast and Path(out_tsv).exists() and Path(out_tsv).stat().st_size > 0:
        print(f"BLAST TSV existant réutilisé : {out_tsv}")
        return

    db_prefix = str(Path(inpactor_fa).with_suffix("")) + "_blastdb"

    subprocess.run(
        ["makeblastdb", "-in", inpactor_fa, "-dbtype", "nucl", "-out", db_prefix],
        check=True
    )

    subprocess.run(
        [
            "blastn",
            "-query", hite_fa,
            "-db", db_prefix,
            "-out", out_tsv,
            "-num_threads", str(threads),
            "-outfmt",
            "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
        ],
        check=True
    )


def parse_blast_best_hits(blast_tsv, inpactor_classes, min_identity, min_qcov):
    best = {}

    with open(blast_tsv) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue

            qseqid_raw = cols[0]
            sseqid = cols[1]
            qseqid = clean_repeat_id(qseqid_raw)

            try:
                pident = float(cols[2])
                aln_len = int(cols[3])
                qlen = int(cols[4])
                evalue = float(cols[10])
                bitscore = float(cols[11])
            except ValueError:
                continue

            if qlen <= 0:
                continue

            qcov = aln_len / qlen * 100

            if pident < min_identity or qcov < min_qcov:
                continue
            if sseqid not in inpactor_classes:
                continue

            new_class = inpactor_classes[sseqid]

            if qseqid not in best or bitscore > best[qseqid][3]:
                best[qseqid] = (
                    new_class,
                    pident,
                    qcov,
                    bitscore,
                    evalue,
                    sseqid,
                    qseqid_raw
                )

    return best


def parse_tesorter_cls(tsv):
    """
    Parse TEsorter *rexdb-plant.cls.tsv.

    The parser is intentionally permissive because TEsorter headers may differ.
    It searches each line for an ID and an LTR lineage such as:
      RLG/Athila, LTR/Gypsy/Athila, Gypsy/Athila
      RLC/SIRE,   LTR/Copia/SIRE,  Copia/SIRE

    Returns:
      dict[clean_repeat_id] = "LTR/Gypsy/FAMILY" or "LTR/Copia/FAMILY"
    """
    if not tsv:
        return {}

    result = {}
    skipped_unknown = 0

    with open(tsv) as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if header is None and any(x.lower() in cols[0].lower() for x in ["seq", "id", "name"]):
                header = [c.strip().lower() for c in cols]
                continue

            raw_id = cols[0].strip()
            seq_id = clean_repeat_id(raw_id)
            text = "\t".join(cols)

            superfamily = None
            family = None

            # Full class/family patterns
            m = re.search(r"LTR/(Gypsy|Copia)/([A-Za-z0-9_.-]+)", text, re.IGNORECASE)
            if m:
                superfamily = m.group(1).capitalize()
                family = m.group(2)
            else:
                # REXdb codes
                m = re.search(r"\b(RLG|RLC)/([A-Za-z0-9_.-]+)", text, re.IGNORECASE)
                if m:
                    code = m.group(1).upper()
                    superfamily = "Gypsy" if code == "RLG" else "Copia"
                    family = m.group(2)
                else:
                    # Separate columns often contain order and lineage
                    tokens = [c.strip() for c in cols]
                    for i, tok in enumerate(tokens):
                        tlow = tok.lower()
                        if tlow in {"gypsy", "rlg"} and i + 1 < len(tokens):
                            superfamily = "Gypsy"
                            family = tokens[i + 1]
                            break
                        if tlow in {"copia", "rlc"} and i + 1 < len(tokens):
                            superfamily = "Copia"
                            family = tokens[i + 1]
                            break

            if not superfamily or not family:
                continue

            family = normalize_family_name(family)
            if not family or family.lower() in {"unknown", "unclassified", "na", "none", "-"}:
                skipped_unknown += 1
                continue

            result[seq_id] = f"LTR/{superfamily}/{family}"

    print(f"TEsorter classifications chargées : {len(result)}")
    if skipped_unknown:
        print(f"TEsorter lignes ignorées car famille inconnue/non classée : {skipped_unknown}")
    return result


def reannotate_rmout(rmout, output, hite_to_inpactor, hite_to_tesorter, report_file):
    total_lines = 0
    candidate_lines = 0
    modified_lines = 0
    inpactor_lines = 0
    tesorter_lines = 0
    unknown_candidate_lines = 0
    tesorter_available_not_used = 0

    counts = defaultdict(int)
    status_counts = defaultdict(int)
    mapping_records = {}

    with open(rmout) as fin, open(output, "w") as fout:
        for line in fin:
            stripped = line.strip()

            if not stripped:
                fout.write(line)
                continue

            if stripped.startswith(("SW", "score", "perc", "Kimura", "There", "Matrix")):
                fout.write(line)
                continue

            cols = line.rstrip("\n").split()

            if len(cols) < 11:
                fout.write(line)
                continue

            total_lines += 1

            repeat_name_raw = cols[9]
            repeat_name = clean_repeat_id(repeat_name_raw)
            repeat_class = cols[10]

            if not is_ltr_to_reannotate(repeat_class):
                fout.write(line)
                continue

            candidate_lines += 1
            original_order = order_from_annotation(repeat_class)

            new_class = None
            status = None
            source = None
            identity = qcov = bitscore = evalue = best_hit = blast_raw_id = "NA"

            # Priority 1: Inpactor2
            if repeat_name in hite_to_inpactor:
                (
                    new_class,
                    identity,
                    qcov,
                    bitscore,
                    evalue,
                    best_hit,
                    blast_raw_id,
                ) = hite_to_inpactor[repeat_name]
                status = "Inpactor2_hit"
                source = "Inpactor2"
                inpactor_lines += 1

                if repeat_name in hite_to_tesorter:
                    tesorter_available_not_used += 1

            else:
                # First keep/add Unknown
                candidate_unknown = unknown_ltr_annotation(repeat_class)

                # Priority 2: TEsorter only for remaining LTR/Gypsy/Unknown or LTR/Copia/Unknown
                if repeat_name in hite_to_tesorter:
                    ts_class = hite_to_tesorter[repeat_name]
                    ts_order = order_from_annotation(ts_class)

                    # Avoid turning Gypsy/Unknown into Copia/family or reverse unless explicitly allowed later.
                    if original_order == ts_order:
                        new_class = ts_class
                        status = "TEsorter_hit_after_no_Inpactor2"
                        source = "TEsorter"
                        tesorter_lines += 1
                    else:
                        new_class = candidate_unknown
                        status = "TEsorter_conflict_order_kept_unknown"
                        source = "None"
                        unknown_candidate_lines += 1
                else:
                    new_class = candidate_unknown
                    status = "No_valid_Inpactor2_or_TEsorter_hit"
                    source = "None"
                    unknown_candidate_lines += 1

            cols[10] = new_class
            fout.write(" ".join(cols) + "\n")

            modified_lines += 1
            counts[new_class] += 1
            status_counts[status] += 1

            if repeat_name not in mapping_records:
                mapping_records[repeat_name] = {
                    "hite_id": repeat_name,
                    "rmout_raw_id": repeat_name_raw,
                    "blast_query_raw_id": blast_raw_id,
                    "original_annotation": repeat_class,
                    "new_annotation": new_class,
                    "status": status,
                    "source": source,
                    "identity": identity,
                    "query_coverage": qcov,
                    "bitscore": bitscore,
                    "evalue": evalue,
                    "best_inpactor2_hit": best_hit,
                    "tesorter_annotation": hite_to_tesorter.get(repeat_name, "NA"),
                    "rmout_candidate_line_count": 1,
                }
            else:
                mapping_records[repeat_name]["rmout_candidate_line_count"] += 1

    with open(report_file, "w") as rep:
        rep.write(f"Total RepeatMasker annotation lines\t{total_lines}\n")
        rep.write(f"Candidate LTR/Gypsy, LTR/Copia or Unknown LTR lines\t{candidate_lines}\n")
        rep.write(f"Lines rewritten in output\t{modified_lines}\n")
        rep.write(f"Lines reannotated by Inpactor2\t{inpactor_lines}\n")
        rep.write(f"Lines reannotated by TEsorter after no Inpactor2 hit\t{tesorter_lines}\n")
        rep.write(f"Candidate lines kept/set to Unknown\t{unknown_candidate_lines}\n")
        rep.write(f"Lines with TEsorter classification but not used because Inpactor2 had priority\t{tesorter_available_not_used}\n")
        rep.write("\nStatus\tCount\n")
        for status, count in sorted(status_counts.items()):
            rep.write(f"{status}\t{count}\n")
        rep.write("\nNew_family\tCount\n")
        for family, count in sorted(counts.items()):
            rep.write(f"{family}\t{count}\n")

    return mapping_records


def write_mapping_table(mapping_records, output):
    with open(output, "w") as f:
        f.write(
            "HiTE_ID_in_RMout\tHiTE_ID_raw_in_RMout\tBLAST_query_raw_ID\t"
            "Original_annotation\tNew_annotation\tStatus\tSource_used\tIdentity\t"
            "Query_coverage\tBitscore\tEvalue\tBest_Inpactor2_hit\t"
            "TEsorter_annotation\tRMout_candidate_line_count\n"
        )

        for hite_id, rec in sorted(mapping_records.items()):
            identity = rec["identity"]
            qcov = rec["query_coverage"]
            bitscore = rec["bitscore"]
            evalue = rec["evalue"]

            if isinstance(identity, float):
                identity = f"{identity:.2f}"
            if isinstance(qcov, float):
                qcov = f"{qcov:.2f}"
            if isinstance(bitscore, float):
                bitscore = f"{bitscore:.1f}"
            if isinstance(evalue, float):
                evalue = f"{evalue:.2e}"

            f.write(
                f"{rec['hite_id']}\t{rec['rmout_raw_id']}\t{rec['blast_query_raw_id']}\t"
                f"{rec['original_annotation']}\t{rec['new_annotation']}\t{rec['status']}\t"
                f"{rec['source']}\t{identity}\t{qcov}\t{bitscore}\t{evalue}\t"
                f"{rec['best_inpactor2_hit']}\t{rec['tesorter_annotation']}\t"
                f"{rec['rmout_candidate_line_count']}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Réannote les lignes LTR/Gypsy, LTR/Copia, LTR/Gypsy/Unknown "
            "et LTR/Copia/Unknown d'un fichier RepeatMasker .out. "
            "Priorité : 1) Inpactor2 par BLAST, 2) TEsorter optionnel seulement "
            "pour les séquences restant Unknown."
        )
    )

    parser.add_argument("--rmout", required=True, help="Fichier RepeatMasker .out à réannoter")
    parser.add_argument("--hite-lib", required=True, help="Librairie consensus utilisée par RepeatMasker")
    parser.add_argument("--inpactor-lib", required=True, help="Librairie Inpactor2 FASTA avec headers #LTR/RLG/famille ou #LTR/RLC/famille")
    parser.add_argument("--out", required=True, help="Nouveau fichier RepeatMasker .out réannoté")

    parser.add_argument("--tesorter-tsv", default=None, help="Optionnel : fichier TEsorter *rexdb-plant.cls.tsv")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-identity", type=float, default=80.0)
    parser.add_argument("--min-qcov", type=float, default=70.0)
    parser.add_argument("--blast-tsv", default="hite_vs_inpactor2.blast.tsv")
    parser.add_argument("--reuse-blast", action="store_true", help="Réutilise --blast-tsv s'il existe déjà")

    args = parser.parse_args()

    inpactor_headers = parse_fasta_headers(args.inpactor_lib)
    inpactor_classes = {}

    for seq_id, header in inpactor_headers.items():
        cls = get_class_from_inpactor2_header(header)
        if is_ltr_family_annotation(cls):
            inpactor_classes[seq_id] = cls

    if not inpactor_classes:
        raise ValueError(
            "Aucune annotation Inpactor2 de type #LTR/RLG/famille ou #LTR/RLC/famille trouvée."
        )

    run_blast(
        args.hite_lib,
        args.inpactor_lib,
        args.blast_tsv,
        args.threads,
        reuse_blast=args.reuse_blast,
    )

    hite_to_inpactor = parse_blast_best_hits(
        args.blast_tsv,
        inpactor_classes,
        args.min_identity,
        args.min_qcov,
    )

    hite_to_tesorter = parse_tesorter_cls(args.tesorter_tsv) if args.tesorter_tsv else {}

    mapping_file = args.out + ".mapping.tsv"
    report_file = args.out + ".summary.txt"

    mapping_records = reannotate_rmout(
        args.rmout,
        args.out,
        hite_to_inpactor,
        hite_to_tesorter,
        report_file,
    )

    write_mapping_table(mapping_records, mapping_file)

    print("Réannotation terminée")
    print(f"Fichier .out réannoté : {args.out}")
    print(f"Table de correspondance : {mapping_file}")
    print(f"Résumé : {report_file}")
    print(f"Hits Inpactor2 retenus : {len(hite_to_inpactor)}")
    if args.tesorter_tsv:
        print(f"Classifications TEsorter chargées : {len(hite_to_tesorter)}")


if __name__ == "__main__":
    main()

