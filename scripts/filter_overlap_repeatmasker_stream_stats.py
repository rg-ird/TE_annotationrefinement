#!/usr/bin/env python3

import argparse
import sys
import tempfile
import subprocess
import os

PRIORITY_BONUS = {
    "Unknown": -10000,
    "mixture": -5000,
    "LTR/mixture": -5000,
    "LTR/Gypsy/mixture": -4000,
    "pararetrovirus": -3000,

    "LTR/Gypsy": 1000,
    "LTR/Copia": 1000,

    "LTR/Gypsy/CRM": 3000,
    "LTR/Gypsy/Athila": 3000,
    "LTR/Gypsy/Ogre": 3000,
    "LTR/Gypsy/Reina": 3000,
    "LTR/Gypsy/Retand": 3000,
    "LTR/Gypsy/Tekay": 3000,
    "LTR/Gypsy/Galadriel": 3000,

    "LTR/Copia/Ale": 3000,
    "LTR/Copia/Alesia": 3000,
    "LTR/Copia/Angela": 3000,
    "LTR/Copia/Bianca": 3000,
    "LTR/Copia/Ikeros": 3000,
    "LTR/Copia/Ivana": 3000,
    "LTR/Copia/SIRE": 3000,
    "LTR/Copia/TAR": 3000,
    "LTR/Copia/Tork": 3000,

    "LINE": 1500,
    "Helitron": 1500,
    "TIR/EnSpm_CACTA": 1500,
    "TIR/hAT": 1500,
    "TIR/MuDR_Mutator": 1500,
    "TIR/PIF_Harbinger": 1500,
    "TIR/Tc1_Mariner": 1500,
}


class Hit:
    def __init__(self, line, idx):
        self.line = line.rstrip("\n")
        self.idx = idx
        f = self.line.split()

        self.score = int(f[0])
        self.seqid = f[4]
        self.start = int(f[5])
        self.end = int(f[6])
        self.name = f[9]
        self.cls = f[10]

        if self.start > self.end:
            self.start, self.end = self.end, self.start

        self.length = self.end - self.start + 1

    def rank(self):
        bonus = PRIORITY_BONUS.get(self.cls, 0)
        specificity = self.cls.count("/") * 100
        return (
            self.score + bonus + specificity,
            self.score,
            specificity,
            self.length,
            -self.idx,
        )


def is_hit_line(line):
    f = line.split()
    if len(f) < 11:
        return False
    try:
        int(f[0])
        int(f[5])
        int(f[6])
        return True
    except Exception:
        return False


def overlap_fraction(a, b):
    ov = max(0, min(a.end, b.end) - max(a.start, b.start) + 1)
    if ov == 0:
        return 0.0
    return ov / min(a.length, b.length)


def removal_reason(kept_hit, removed_hit):
    if kept_hit.score > removed_hit.score:
        return "removed_by_raw_score"

    if kept_hit.score <= removed_hit.score and kept_hit.rank() > removed_hit.rank():
        return "removed_by_priority"

    return "removed_by_tie_or_length"


def filter_component(component, min_overlap, stats):
    kept = []

    for h in sorted(component, key=lambda x: x.rank(), reverse=True):
        redundant = False

        for k in kept:
            if overlap_fraction(h, k) >= min_overlap:
                redundant = True
                reason = removal_reason(k, h)
                stats[reason] += 1
                break

        if not redundant:
            kept.append(h)

    return sorted(kept, key=lambda x: (x.seqid, x.start, x.end))


def sort_hits_to_temp(input_file):
    tmp = tempfile.NamedTemporaryFile(delete=False, mode="w")
    tmp.close()

    cmd = (
        f"awk 'NF>=11 && $1 ~ /^[0-9]+$/ {{print}}' {input_file} "
        f"| sort -k5,5 -k6,6n -k7,7n > {tmp.name}"
    )

    subprocess.run(cmd, shell=True, check=True)
    return tmp.name


def main():
    parser = argparse.ArgumentParser(
        description="Fast streaming overlap filter for RepeatMasker .out with removal statistics."
    )

    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)

    parser.add_argument(
        "--min-overlap-fraction",
        type=float,
        default=0.8,
        help="Overlap fraction relative to the shortest annotation. Default: 0.8",
    )

    parser.add_argument(
        "--sort",
        action="store_true",
        help="Sort RepeatMasker hits before filtering. Recommended.",
    )

    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Do not write RepeatMasker header lines.",
    )

    args = parser.parse_args()

    headers = []

    if not args.no_header:
        with open(args.input) as f:
            for line in f:
                if is_hit_line(line):
                    break
                headers.append(line.rstrip("\n"))

    input_to_read = args.input
    tmp_sorted = None

    if args.sort:
        print("Sorting input hits...", file=sys.stderr)
        tmp_sorted = sort_hits_to_temp(args.input)
        input_to_read = tmp_sorted

    stats = {
        "removed_by_raw_score": 0,
        "removed_by_priority": 0,
        "removed_by_tie_or_length": 0,
    }

    n_in = 0
    n_out = 0
    idx = 0

    current_seq = None
    component = []
    component_end = -1

    with open(input_to_read) as inp, open(args.output, "w") as out:

        if not args.no_header:
            for h in headers:
                out.write(h + "\n")

        for line in inp:
            if not is_hit_line(line):
                continue

            h = Hit(line, idx)
            idx += 1
            n_in += 1

            if n_in % 100000 == 0:
                print(f"Processed {n_in} annotations...", file=sys.stderr)

            if current_seq is None:
                current_seq = h.seqid
                component = [h]
                component_end = h.end
                continue

            if h.seqid != current_seq or h.start > component_end:
                filtered = filter_component(
                    component,
                    args.min_overlap_fraction,
                    stats
                )

                for x in filtered:
                    out.write(x.line + "\n")

                n_out += len(filtered)
                out.flush()

                current_seq = h.seqid
                component = [h]
                component_end = h.end

            else:
                component.append(h)
                if h.end > component_end:
                    component_end = h.end

        if component:
            filtered = filter_component(
                component,
                args.min_overlap_fraction,
                stats
            )

            for x in filtered:
                out.write(x.line + "\n")

            n_out += len(filtered)

    if tmp_sorted:
        os.remove(tmp_sorted)

    removed_total = n_in - n_out

    print("\nFiltering summary", file=sys.stderr)
    print("-----------------", file=sys.stderr)
    print(f"Input annotations:  {n_in}", file=sys.stderr)
    print(f"Output annotations: {n_out}", file=sys.stderr)
    print(f"Removed:            {removed_total}", file=sys.stderr)

    if removed_total > 0:
        print("\nRemoval causes", file=sys.stderr)
        print("--------------", file=sys.stderr)

        for key, label in [
            ("removed_by_raw_score", "Removed by raw RepeatMasker score"),
            ("removed_by_priority", "Removed by priority/class bonus"),
            ("removed_by_tie_or_length", "Removed by tie/length/specificity"),
        ]:
            value = stats[key]
            pct = 100 * value / removed_total
            print(f"{label}: {value} ({pct:.2f}%)", file=sys.stderr)

    else:
        print("No annotation removed.", file=sys.stderr)


if __name__ == "__main__":
    main()
