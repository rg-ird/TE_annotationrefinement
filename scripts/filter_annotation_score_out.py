#!/usr/bin/env python3

import argparse
import sys
import tempfile
import subprocess
import os


class Hit:
    def __init__(self, line, idx):
        self.line = line.rstrip("\n")
        self.idx = idx
        f = self.line.split()

        self.score = int(f[0])
        self.seqid = f[4]
        self.start = int(f[5])
        self.end = int(f[6])

        if self.start > self.end:
            self.start, self.end = self.end, self.start

        self.length = self.end - self.start + 1

    def rank(self):
        return (
            self.score,
            self.length,
            -self.idx
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


def filter_component(component, min_overlap, stats):
    kept = []

    for h in sorted(component, key=lambda x: x.rank(), reverse=True):
        redundant = False

        for k in kept:
            if overlap_fraction(h, k) >= min_overlap:
                redundant = True
                stats["removed_by_score"] += 1
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
        description="Filter overlapping RepeatMasker .out annotations by raw score only."
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

    stats = {"removed_by_score": 0}

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

                current_seq = h.seqid
                component = [h]
                component_end = h.end

            else:
                component.append(h)
                component_end = max(component_end, h.end)

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
    print(f"Removed by score:   {stats['removed_by_score']}", file=sys.stderr)


if __name__ == "__main__":
    main()
