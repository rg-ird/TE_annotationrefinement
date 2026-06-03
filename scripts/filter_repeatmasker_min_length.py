#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict


def is_rm_hit(line):
    fields = line.split()
    if len(fields) < 11:
        return False
    try:
        int(fields[0])
        int(fields[5])
        int(fields[6])
        return True
    except ValueError:
        return False


def get_hit_length(line):
    fields = line.split()
    start = int(fields[5])
    end = int(fields[6])
    if start > end:
        start, end = end, start
    return end - start + 1


def get_class(line):
    fields = line.split()
    return fields[10]


def main():
    parser = argparse.ArgumentParser(
        description="Filtre un fichier RepeatMasker .out en retirant les annotations plus courtes qu'une longueur minimale."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Fichier RepeatMasker .out d'entrée"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Fichier RepeatMasker .out filtré"
    )

    parser.add_argument(
        "-m", "--min-length",
        type=int,
        default=100,
        help="Longueur minimale à conserver. Default: 100 bp"
    )

    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Ne pas conserver les lignes d'en-tête"
    )

    args = parser.parse_args()

    total_hits = 0
    kept_hits = 0
    removed_hits = 0

    total_bp = 0
    kept_bp = 0
    removed_bp = 0

    removed_by_class = defaultdict(lambda: [0, 0])
    kept_by_class = defaultdict(lambda: [0, 0])

    with open(args.input) as inp, open(args.output, "w") as out:
        for line in inp:
            if not is_rm_hit(line):
                if not args.no_header:
                    out.write(line)
                continue

            length = get_hit_length(line)
            repeat_class = get_class(line)

            total_hits += 1
            total_bp += length

            if length >= args.min_length:
                out.write(line)
                kept_hits += 1
                kept_bp += length
                kept_by_class[repeat_class][0] += 1
                kept_by_class[repeat_class][1] += length
            else:
                removed_hits += 1
                removed_bp += length
                removed_by_class[repeat_class][0] += 1
                removed_by_class[repeat_class][1] += length

    print("===== Résumé du filtrage =====")
    print(f"Fichier entrée          : {args.input}")
    print(f"Fichier sortie          : {args.output}")
    print(f"Longueur minimale       : {args.min_length} bp")
    print()
    print(f"Annotations totales     : {total_hits}")
    print(f"Annotations conservées  : {kept_hits}")
    print(f"Annotations retirées    : {removed_hits}")
    print(f"% annotations retirées  : {100 * removed_hits / total_hits:.4f}%" if total_hits else "% annotations retirées  : NA")
    print()
    print(f"Bases annotées totales  : {total_bp}")
    print(f"Bases conservées        : {kept_bp}")
    print(f"Bases retirées          : {removed_bp}")
    print(f"% bases retirées        : {100 * removed_bp / total_bp:.4f}%" if total_bp else "% bases retirées        : NA")
    print()
    print("Annotations retirées par classe")
    print("class\tremoved_count\tremoved_bp")
    for cls in sorted(removed_by_class):
        count, bp = removed_by_class[cls]
        print(f"{cls}\t{count}\t{bp}")


if __name__ == "__main__":
    main()
