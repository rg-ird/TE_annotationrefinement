#!/usr/bin/env python3

import argparse
import subprocess
import re
from pathlib import Path
from collections import defaultdict


def clean_repeat_id(seqid):
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
    Inpactor2 :
    #LTR/RLG/TAT -> LTR/Gypsy/TAT
    #LTR/RLC/SIRE -> LTR/Copia/SIRE
    """

    m = re.search(r"#?LTR/(RLG|RLC)/([A-Za-z0-9_.-]+)", header)

    if not m:
        return None

    code = m.group(1)
    family = m.group(2)

    if code == "RLG":
        return f"LTR/Gypsy/{family}"

    if code == "RLC":
        return f"LTR/Copia/{family}"

    return None


def is_ltr_family_annotation(annotation):
    return annotation and (
        annotation.startswith("LTR/Gypsy/") or
        annotation.startswith("LTR/Copia/")
    ) and not annotation.lower().endswith("/unclassified")


def is_ltr_to_reannotate(annotation):
    return annotation in {
        "LTR/Gypsy",
        "LTR/Copia",
    }


def unknown_ltr_annotation(annotation):
    if annotation == "LTR/Gypsy":
        return "LTR/Gypsy/Unknown"

    if annotation == "LTR/Copia":
        return "LTR/Copia/Unknown"

    return annotation


def run_blast(hite_fa, inpactor_fa, out_tsv, threads):
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

            pident = float(cols[2])
            aln_len = int(cols[3])
            qlen = int(cols[4])
            evalue = float(cols[10])
            bitscore = float(cols[11])

            qcov = aln_len / qlen * 100

            if pident < min_identity:
                continue

            if qcov < min_qcov:
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


def reannotate_rmout(rmout, output, hite_to_family, report_file):
    total_lines = 0
    candidate_lines = 0
    modified_lines = 0
    unchanged_candidate_lines = 0

    counts = defaultdict(int)

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

            if is_ltr_to_reannotate(repeat_class):
                candidate_lines += 1

                if repeat_name in hite_to_family:
                    new_class = hite_to_family[repeat_name][0]
                    cols[10] = new_class

                    modified_lines += 1
                    counts[new_class] += 1

                    fout.write(" ".join(cols) + "\n")
                else:
                    new_class = unknown_ltr_annotation(repeat_class)
                    cols[10] = new_class

                    modified_lines += 1
                    unchanged_candidate_lines += 1
                    counts[new_class] += 1

                    fout.write(" ".join(cols) + "\n")
            else:
                fout.write(line)

    with open(report_file, "w") as rep:
        rep.write(f"Total RepeatMasker annotation lines\t{total_lines}\n")
        rep.write(f"Candidate LTR/Gypsy or LTR/Copia lines\t{candidate_lines}\n")
        rep.write(f"Reannotated lines\t{modified_lines}\n")
        rep.write(f"Candidate lines set to Unknown because no valid Inpactor2 hit\t{unchanged_candidate_lines}\n")
        rep.write("\nNew_family\tCount\n")

        for family, count in sorted(counts.items()):
            rep.write(f"{family}\t{count}\n")


def write_mapping_table(hite_to_family, output):
    with open(output, "w") as f:
        f.write(
            "HiTE_ID_in_RMout\tHiTE_ID_raw\tNew_annotation\tIdentity\t"
            "Query_coverage\tBitscore\tEvalue\tBest_Inpactor2_hit\n"
        )

        for hite_id, values in sorted(hite_to_family.items()):
            new_class, pident, qcov, bitscore, evalue, subject, raw_id = values
            f.write(
                f"{hite_id}\t{raw_id}\t{new_class}\t{pident:.2f}\t"
                f"{qcov:.2f}\t{bitscore:.1f}\t{evalue:.2e}\t{subject}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Réannote uniquement les lignes LTR/Gypsy et LTR/Copia "
            "d'un fichier RepeatMasker .out à partir d'une librairie Inpactor2. "
            "Sans correspondance Inpactor2 valide, ajoute /Unknown."
        )
    )

    parser.add_argument("--rmout", required=True)
    parser.add_argument("--hite-lib", required=True)
    parser.add_argument("--inpactor-lib", required=True)
    parser.add_argument("--out", required=True)

    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-identity", type=float, default=80.0)
    parser.add_argument("--min-qcov", type=float, default=70.0)
    parser.add_argument("--blast-tsv", default="hite_vs_inpactor2.blast.tsv")

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
        args.threads
    )

    hite_to_family = parse_blast_best_hits(
        args.blast_tsv,
        inpactor_classes,
        args.min_identity,
        args.min_qcov
    )

    mapping_file = args.out + ".mapping.tsv"
    report_file = args.out + ".summary.txt"

    write_mapping_table(hite_to_family, mapping_file)

    reannotate_rmout(
        args.rmout,
        args.out,
        hite_to_family,
        report_file
    )

    print("Réannotation terminée")
    print(f"Fichier .out réannoté : {args.out}")
    print(f"Correspondances HiTE-Inpactor2 : {mapping_file}")
    print(f"Résumé : {report_file}")


if __name__ == "__main__":
    main()
