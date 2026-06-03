#!/usr/bin/env python3

import argparse
from collections import defaultdict


def parse_repeatmasker_out(out_file):
    ann = defaultdict(list)

    with open(out_file) as f:
        for line in f:
            if not line.strip():
                continue

            fields = line.split()

            try:
                score = int(fields[0])
            except Exception:
                continue

            if len(fields) < 11:
                continue

            seqid = fields[4]
            start = int(fields[5])
            end = int(fields[6])
            strand = fields[8]
            repeat_name = fields[9]
            repeat_class = fields[10]

            if start > end:
                start, end = end, start

            ann[seqid].append({
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand,
                "score": score,
                "repeat_name": repeat_name,
                "repeat_class": repeat_class,
                "line": line.rstrip("\n")
            })

    return ann


def overlap(a, b):
    s = max(a["start"], b["start"])
    e = min(a["end"], b["end"])

    if s <= e:
        return s, e, e - s + 1

    return None


def main():
    parser = argparse.ArgumentParser(
        description="Compte les chevauchements entre annotations TE dans un fichier RepeatMasker .out."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Fichier RepeatMasker .out d'entrée"
    )

    parser.add_argument(
        "-s", "--summary",
        required=True,
        help="Fichier TSV de résumé par séquence"
    )

    parser.add_argument(
        "-d", "--details",
        required=True,
        help="Fichier TSV détaillant chaque chevauchement"
    )

    args = parser.parse_args()

    annotations = parse_repeatmasker_out(args.input)

    total_annotations = 0
    total_annotated_bp = 0
    total_overlap_events = 0
    total_overlap_bp = 0

    with open(args.summary, "w") as summary, open(args.details, "w") as details:

        summary.write(
            "seqid\tn_annotations\tannotated_bp\toverlap_events\toverlap_bp\toverlap_percent\n"
        )

        details.write(
            "seqid\tann1_start\tann1_end\tann1_score\tann1_name\tann1_class\t"
            "ann2_start\tann2_end\tann2_score\tann2_name\tann2_class\t"
            "overlap_start\toverlap_end\toverlap_length\n"
        )

        for seqid, feats in annotations.items():
            feats.sort(key=lambda x: x["start"])

            n_annotations = len(feats)
            annotated_bp = sum(f["end"] - f["start"] + 1 for f in feats)
            overlap_events = 0
            overlap_bp = 0

            total_annotations += n_annotations
            total_annotated_bp += annotated_bp

            for i in range(len(feats)):
                a = feats[i]

                for j in range(i + 1, len(feats)):
                    b = feats[j]

                    if b["start"] > a["end"]:
                        break

                    ov = overlap(a, b)

                    if ov:
                        ov_start, ov_end, ov_len = ov

                        overlap_events += 1
                        overlap_bp += ov_len

                        details.write(
                            f"{seqid}\t"
                            f"{a['start']}\t{a['end']}\t{a['score']}\t{a['repeat_name']}\t{a['repeat_class']}\t"
                            f"{b['start']}\t{b['end']}\t{b['score']}\t{b['repeat_name']}\t{b['repeat_class']}\t"
                            f"{ov_start}\t{ov_end}\t{ov_len}\n"
                        )

            total_overlap_events += overlap_events
            total_overlap_bp += overlap_bp

            overlap_percent = (
                overlap_bp / annotated_bp * 100 if annotated_bp > 0 else 0
            )

            summary.write(
                f"{seqid}\t{n_annotations}\t{annotated_bp}\t"
                f"{overlap_events}\t{overlap_bp}\t{overlap_percent:.4f}\n"
            )

    global_percent = (
        total_overlap_bp / total_annotated_bp * 100
        if total_annotated_bp > 0 else 0
    )

    print("===== Résumé global =====")
    print(f"Annotations totales        : {total_annotations}")
    print(f"Bases annotées totales     : {total_annotated_bp}")
    print(f"Événements de chevauchement: {total_overlap_events}")
    print(f"Bases chevauchantes        : {total_overlap_bp}")
    print(f"Pourcentage chevauchement  : {global_percent:.4f}%")
    print(f"Résumé par séquence        : {args.summary}")
    print(f"Détails des chevauchements : {args.details}")


if __name__ == "__main__":
    main()
