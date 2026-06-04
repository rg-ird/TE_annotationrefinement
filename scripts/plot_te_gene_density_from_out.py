#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict, Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Classification harmonisée pour le .out final Inpactor2 + TEsorter
# -----------------------------------------------------------------------------

def normalize_class(cls):
    cls = cls.strip()
    mapping = {
        "RC/Helitron": "Helitron",
        "LINE": "LINE/L1",
        "LTR/Gypsy": "LTR/Gypsy/Unknown",
        "LTR/Gypsy/unclassified": "LTR/Gypsy/Unknown",
        "LTR/Gypsy/Tekay": "LTR/Gypsy/TEKAY-DEL",
        "LTR/Gypsy/Galadriel": "LTR/Gypsy/GALADRIEL",
        "LTR/Gypsy/Reina": "LTR/Gypsy/REINA",
        "LTR/Gypsy/Athila": "LTR/Gypsy/ATHILA",
        "LTR/Gypsy/Tat": "LTR/Gypsy/TAT",
        "LTR/Gypsy/TatI": "LTR/Gypsy/TAT",
        "LTR/Gypsy/TatII": "LTR/Gypsy/TAT",
        "LTR/Gypsy/TatIII": "LTR/Gypsy/TAT",
        "LTR/Gypsy/Ogre": "LTR/Gypsy/TAT",
        "LTR/Gypsy/Retand": "LTR/Gypsy/TAT",
        "LTR/Copia": "LTR/Copia/Unknown",
        "LTR/Copia/": "LTR/Copia/Unknown",
        "LTR/Copia/unclassified": "LTR/Copia/Unknown",
        "LTR/Copia/Ale": "LTR/Copia/ALE-RETROFIT",
        "LTR/Copia/Angela": "LTR/Copia/ANGELA",
        "LTR/Copia/Bianca": "LTR/Copia/BIANCA",
        "LTR/Copia/Ivana": "LTR/Copia/IVANA-ORYCO",
        "LTR/Copia/Ikeros": "LTR/Copia/TAR-TORK",
        "LTR/Copia/Tork": "LTR/Copia/TAR-TORK",
        "LTR/Copia/TAR": "LTR/Copia/TAR-TORK",
    }
    return mapping.get(cls, cls)


def chrom_key(seqid, mode="auto"):
    s = str(seqid).strip()
    if mode == "exact":
        return s

    m = re.search(r"(?i)(chr(?:omosome)?[_\.-]*0*\d+|chr[_\.-]*[A-Za-z]+)", s)
    if m:
        k = m.group(1)
        k = re.sub(r"(?i)^chromosome", "Chr", k)
        k = re.sub(r"(?i)^chr", "Chr", k)
        k = k.replace("_", "").replace(".", "").replace("-", "")
        return k

    return re.split(r"[|\s]", s)[0].split(".")[-1]


def choose_category(repeat_class, mode):
    cls = normalize_class(repeat_class)

    if mode == "level1":
        if cls.startswith("SINE/"):
            return "SINE"
        if cls.startswith("LINE/"):
            return "LINE"
        if cls.startswith("LTR/"):
            return "LTR"
        if cls.startswith("DNA/"):
            return "DNA"
        if cls == "Helitron" or cls.startswith("Helitron/"):
            return "Helitron"
        if cls == "Unknown" or cls.startswith("Unknown/"):
            return "Unknown"
        return None

    if mode == "ltr_superfamily":
        if cls.startswith("LTR/Gypsy/"):
            return "Gypsy"
        if cls.startswith("LTR/Copia/"):
            return "Copia"
        return None

    if mode == "gypsy_family":
        if cls.startswith("LTR/Gypsy/"):
            return cls.split("/", 2)[2]
        return None

    if mode == "copia_family":
        if cls.startswith("LTR/Copia/"):
            return cls.split("/", 2)[2]
        return None

    raise ValueError(f"Unknown mode: {mode}")


# -----------------------------------------------------------------------------
# Lectures fichiers, sans pandas
# -----------------------------------------------------------------------------

def read_repeatmasker_out(path, mode, chrom_match="auto"):
    rows = []
    skipped = 0
    categories = Counter()
    chroms = set()
    keys = set()

    with open(path) as f:
        for line in f:
            if not line.strip():
                continue

            fields = line.split()
            try:
                int(fields[0])
            except Exception:
                skipped += 1
                continue

            if len(fields) < 11:
                skipped += 1
                continue

            try:
                seqid = fields[4]
                start = int(fields[5]) - 1   # 0-based half-open pour calculs
                end = int(fields[6])
                repeat_class = fields[10]
            except Exception:
                skipped += 1
                continue

            if start > end:
                start, end = end, start

            category = choose_category(repeat_class, mode)
            if category is None:
                continue

            key = chrom_key(seqid, chrom_match)
            row = {
                "chrom": seqid,
                "chrom_key": key,
                "start": start,
                "end": end,
                "category": category,
                "class": normalize_class(repeat_class),
            }
            rows.append(row)
            categories[category] += 1
            chroms.add(seqid)
            keys.add(key)

    if not rows:
        raise ValueError("No valid TE annotations retained. Check --mode and input .out file.")

    print(f"RepeatMasker lines skipped/header: {skipped}")
    print("TE categories retained:")
    for cat, n in categories.most_common():
        print(f"  {cat}: {n}")
    print(f"TE chromosomes: {len(chroms)} original names; {len(keys)} matching keys")
    return rows


def read_gene_gff(path, feature_type="gene", chrom_match="auto"):
    rows = []
    feature_counts = Counter()
    chroms = set()
    keys = set()

    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            feature_counts[parts[2]] += 1
            if parts[2] != feature_type:
                continue

            try:
                start = int(parts[3]) - 1
                end = int(parts[4])
            except Exception:
                continue

            key = chrom_key(parts[0], chrom_match)
            rows.append({
                "chrom": parts[0],
                "chrom_key": key,
                "start": start,
                "end": end,
            })
            chroms.add(parts[0])
            keys.add(key)

    print(f"Gene feature requested: {feature_type}")
    print(f"Gene annotations retained: {len(rows)}")

    if not rows:
        print("WARNING: no gene annotations retained.")
        if feature_counts:
            print("Feature types found in GFF3:")
            for k, v in feature_counts.most_common(20):
                print(f"  {k}: {v}")
    else:
        print(f"Gene chromosomes: {len(chroms)} original names; {len(keys)} matching keys")

    return rows


def read_sizes(path, chrom_match="auto"):
    sizes = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                sizes.append({
                    "chrom": parts[0],
                    "chrom_key": chrom_key(parts[0], chrom_match),
                    "length": int(parts[1]),
                })
            except Exception:
                continue

    if not sizes:
        raise ValueError("No chromosome sizes read. Provide a .fai or two-column sizes file.")
    return sizes


def infer_sizes_from_te(te_rows):
    max_end = {}
    chrom_for_key = {}
    for r in te_rows:
        key = r["chrom_key"]
        chrom_for_key.setdefault(key, r["chrom"])
        max_end[key] = max(max_end.get(key, 0), int(r["end"]))
    return [{"chrom": chrom_for_key[k], "chrom_key": k, "length": v} for k, v in max_end.items()]


# -----------------------------------------------------------------------------
# Indexation et densités
# -----------------------------------------------------------------------------

def index_by_key(rows):
    idx = defaultdict(list)
    for r in rows:
        idx[r["chrom_key"]].append(r)
    for key in idx:
        idx[key].sort(key=lambda x: x["start"])
    return idx


def overlap_sum(rows, win_start, win_end, category=None):
    total = 0
    if not rows:
        return 0
    for r in rows:
        # Optimisation simple grâce au tri par start
        if r["start"] >= win_end:
            break
        if r["end"] <= win_start:
            continue
        if category is not None and r.get("category") != category:
            continue
        ov_start = max(r["start"], win_start)
        ov_end = min(r["end"], win_end)
        if ov_end > ov_start:
            total += ov_end - ov_start
    return total


def build_windows(chrom_sizes, window, step, min_size):
    windows = []
    for row in chrom_sizes:
        chrom = row["chrom"]
        key = row["chrom_key"]
        chrom_len = int(row["length"])
        if chrom_len < min_size:
            continue

        for start in range(0, chrom_len, step):
            end = min(start + window, chrom_len)
            if end <= start:
                continue
            windows.append({
                "chrom": chrom,
                "chrom_key": key,
                "start": start,
                "end": end,
                "window_length": end - start,
            })

    if not windows:
        raise ValueError("No windows created. Check --sizes and --min-size.")
    return windows


def compute_density(te_rows, gene_rows, windows):
    categories = sorted(set(r["category"] for r in te_rows))
    te_idx = index_by_key(te_rows)
    gene_idx = index_by_key(gene_rows) if gene_rows else {}
    result = []

    for win in windows:
        key = win["chrom_key"]
        start = int(win["start"])
        end = int(win["end"])
        window_length = int(win["window_length"])

        row = {
            "chrom": win["chrom"],
            "chrom_key": key,
            "start": start,
            "end": end,
            "window_length": window_length,
        }

        te_sub = te_idx.get(key, [])
        for cat in categories:
            bp = overlap_sum(te_sub, start, end, category=cat)
            row[f"{cat}_bp"] = bp
            row[f"{cat}_density"] = bp / window_length if window_length else 0

        if gene_rows is not None:
            gene_bp = overlap_sum(gene_idx.get(key, []), start, end, category=None)
            row["GENE_bp"] = gene_bp
            row["GENE_density"] = gene_bp / window_length if window_length else 0

        result.append(row)

    if gene_rows is not None:
        total_gene_bp = sum(r.get("GENE_bp", 0) for r in result)
        if total_gene_bp == 0:
            print("WARNING: GFF genes were read but no overlap was found with plot windows.")
            print("Most likely cause: chromosome names differ between .out/.fai and GFF3.")
            print("Try default --chrom-match auto, or check names with: cut -f1 genes.gff3 | sort -u")
        else:
            print(f"Total gene bp counted across windows: {total_gene_bp}")

    return result, categories


def write_table(path, rows, categories, with_genes):
    header = ["chrom", "chrom_key", "start", "end", "window_length"]
    for cat in categories:
        header.extend([f"{cat}_bp", f"{cat}_density"])
    if with_genes:
        header.extend(["GENE_bp", "GENE_density"])

    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            vals = []
            for h in header:
                v = r.get(h, 0)
                vals.append(str(v))
            out.write("\t".join(vals) + "\n")


# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

def group_by_chrom(rows):
    grouped = defaultdict(list)
    chrom_order = []
    seen = set()
    for r in rows:
        chrom = r["chrom"]
        if chrom not in seen:
            chrom_order.append(chrom)
            seen.add(chrom)
        grouped[chrom].append(r)
    return chrom_order, grouped


def plot_by_chromosome(rows, categories, prefix, gene_label=None, dpi=300):
    chrom_order, grouped = group_by_chrom(rows)

    for chrom in chrom_order:
        subset = sorted(grouped[chrom], key=lambda r: r["start"])
        x = [(r["start"] + r["end"]) / 2 / 1e6 for r in subset]

        fig, ax = plt.subplots(figsize=(12, 4), dpi=dpi)

        for cat in categories:
            y = [r.get(f"{cat}_density", 0) for r in subset]
            ax.plot(x, y, linewidth=1.2, label=cat)

        if subset and "GENE_density" in subset[0]:
            y_gene = [r.get("GENE_density", 0) for r in subset]
            ax.plot(x, y_gene, linewidth=1.8, linestyle="-", color="black", label=gene_label or "Genes")

        ax.set_xlabel("Position (Mb)")
        ax.set_ylabel("Density / bp fraction")
        ax.set_title(f"{chrom} - TE density")
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=8, ncol=2, frameon=False)
        fig.tight_layout()

        safe_chrom = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(chrom))
        fig.savefig(f"{prefix}_{safe_chrom}.pdf")
        plt.close(fig)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot TE density from RepeatMasker .out with optional gene density. No pandas version."
    )

    parser.add_argument("--out-file", required=True, help="Input RepeatMasker .out file")
    parser.add_argument("--sizes", help="Chromosome sizes file: FASTA .fai or two-column file")
    parser.add_argument("--gff", help="Optional gene GFF3 file")
    parser.add_argument("--gene-feature", default="gene", help="GFF3 feature type to use, default: gene")
    parser.add_argument("--window", type=int, default=200000, help="Window size, default: 200000")
    parser.add_argument("--step", type=int, help="Step size, default = window")
    parser.add_argument("--min-size", type=int, default=1000000, help="Minimum chromosome/contig size")
    parser.add_argument("--prefix", default="te_gene_density", help="Output prefix")
    parser.add_argument("--dpi", type=int, default=300, help="PDF dpi metadata, default: 300")
    parser.add_argument(
        "--chrom-match",
        choices=["auto", "exact"],
        default="auto",
        help="auto matches names like Chr01 and CC1.8.Chr01; exact requires identical names. Default: auto"
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=["level1", "ltr_superfamily", "gypsy_family", "copia_family"],
        help=(
            "level1 = SINE/LINE/LTR/DNA/Helitron/Unknown; "
            "ltr_superfamily = Gypsy vs Copia; "
            "gypsy_family = families of LTR/Gypsy; "
            "copia_family = families of LTR/Copia"
        )
    )

    args = parser.parse_args()
    step = args.step if args.step else args.window

    te_rows = read_repeatmasker_out(args.out_file, args.mode, args.chrom_match)

    gene_rows = None
    gene_label = None
    if args.gff:
        gene_rows = read_gene_gff(args.gff, args.gene_feature, args.chrom_match)
        gene_label = f"Genes ({args.gene_feature})"

    if args.sizes:
        chrom_sizes = read_sizes(args.sizes, args.chrom_match)
    else:
        print("WARNING: no --sizes provided. Sizes inferred from max RepeatMasker end position.")
        chrom_sizes = infer_sizes_from_te(te_rows)

    windows = build_windows(chrom_sizes, args.window, step, args.min_size)
    result, categories = compute_density(te_rows, gene_rows, windows)

    out_table = f"{args.prefix}_table.tsv"
    write_table(out_table, result, categories, with_genes=(gene_rows is not None))
    print(f"Output table: {out_table}")

    plot_by_chromosome(result, categories, args.prefix, gene_label=gene_label, dpi=args.dpi)
    print(f"PDF files written with prefix: {args.prefix}_<chromosome>.pdf")


if __name__ == "__main__":
    main()

