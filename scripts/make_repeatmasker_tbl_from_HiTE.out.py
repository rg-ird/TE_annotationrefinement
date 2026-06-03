#!/usr/bin/env python3

import argparse
from collections import defaultdict

CLASS_ORDER = [
    "SINE/tRNA",
    "Unknown",
    "RC/Helitron",
    "LTR/Copia",
    "LTR/ERV",
    "LTR/Gypsy",
    "DNA/CMC-EnSpm",
    "DNA/hAT",
    "DNA/Merlin",
    "DNA/MULE",
    "DNA/PIF-Harbinger",
    "DNA/TcMar",
    "LINE/L1",
]

SUMMARY_GROUPS = {
    "All SINE": lambda cls: cls.startswith("SINE/"),
    "All LINE": lambda cls: cls.startswith("LINE/"),
    "All LTR": lambda cls: cls.startswith("LTR/"),
    "All DNA": lambda cls: cls.startswith("DNA/"),
    "All RC": lambda cls: cls.startswith("RC/"),
    "All Unknown": lambda cls: cls == "Unknown" or cls.startswith("Unknown/"),
}


def parse_out(path):
    counts = defaultdict(int)
    lengths = defaultdict(int)
    total_annotated_bp = 0
    seq_lengths = {}

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            fields = line.split()

            try:
                int(fields[0])
            except Exception:
                continue

            if len(fields) < 11:
                continue

            seqid = fields[4]
            start = int(fields[5])
            end = int(fields[6])
            left = fields[7].replace("(", "").replace(")", "")
            repeat_class = fields[10]

            if start > end:
                start, end = end, start

            length = end - start + 1

            counts[repeat_class] += 1
            lengths[repeat_class] += length
            total_annotated_bp += length

            try:
                left = int(left)
                seq_len = end + left
                seq_lengths[seqid] = max(seq_lengths.get(seqid, 0), seq_len)
            except Exception:
                pass

    genome_size = sum(seq_lengths.values())

    return counts, lengths, total_annotated_bp, genome_size


def pct_text(bp, genome_size):
    if genome_size:
        return f"{100 * bp / genome_size:8.3f}"
    return "      NA"


def write_tbl(output, counts, lengths, total_annotated_bp, genome_size):
    total_elements = sum(counts.values())

    with open(output, "w") as out:
        out.write("==================================================\n")
        out.write("RepeatMasker-like summary rebuilt from .out file\n")
        out.write("==================================================\n\n")

        out.write(f"Genome size estimated:       {genome_size if genome_size else 'NA'} bp\n")
        out.write(f"Total interspersed repeats:  {total_annotated_bp} bp\n")
        out.write(f"Number of elements:          {total_elements}\n\n")

        out.write("Detailed classification\n")
        out.write("--------------------------------------------------------------------\n")
        out.write(f"{'Class':<28} {'Count':>8} {'bp masked':>16} {'% genome':>10}\n")
        out.write("--------------------------------------------------------------------\n")

        for cls in CLASS_ORDER:
            c = counts.get(cls, 0)
            bp = lengths.get(cls, 0)
            out.write(f"{cls:<28} {c:>8} {bp:>16} {pct_text(bp, genome_size):>10}\n")

        other_count = 0
        other_bp = 0

        for cls in sorted(counts):
            if cls not in CLASS_ORDER:
                other_count += counts[cls]
                other_bp += lengths[cls]

        if other_count > 0:
            out.write(f"{'Other':<28} {other_count:>8} {other_bp:>16} {pct_text(other_bp, genome_size):>10}\n")

        out.write("--------------------------------------------------------------------\n\n")

        out.write("Summary groups\n")
        out.write("--------------------------------------------------------------------\n")
        out.write(f"{'Group':<28} {'Count':>8} {'bp masked':>16} {'% genome':>10}\n")
        out.write("--------------------------------------------------------------------\n")

        for group_name, rule in SUMMARY_GROUPS.items():
            group_count = sum(counts[cls] for cls in counts if rule(cls))
            group_bp = sum(lengths[cls] for cls in lengths if rule(cls))

            out.write(
                f"{group_name:<28} {group_count:>8} "
                f"{group_bp:>16} {pct_text(group_bp, genome_size):>10}\n"
            )

        out.write("--------------------------------------------------------------------\n")

        out.write(
            f"{'Total':<28} {total_elements:>8} "
            f"{total_annotated_bp:>16} {pct_text(total_annotated_bp, genome_size):>10}\n"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Rebuild a RepeatMasker-like .tbl summary from a RepeatMasker .out file."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input RepeatMasker .out file"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output .tbl file"
    )

    args = parser.parse_args()

    counts, lengths, total_annotated_bp, genome_size = parse_out(args.input)

    write_tbl(
        args.output,
        counts,
        lengths,
        total_annotated_bp,
        genome_size
    )

    print(f"Input file: {args.input}")
    print(f"Output tbl: {args.output}")
    print(f"Genome size estimated: {genome_size if genome_size else 'NA'} bp")
    print(f"Total annotated bp: {total_annotated_bp}")
    print(f"Total elements: {sum(counts.values())}")


if __name__ == "__main__":
    main()
