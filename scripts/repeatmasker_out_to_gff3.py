#!/usr/bin/env python3

import argparse
import re


def sanitize_attr(value):
    value = str(value)
    value = value.replace(";", "_")
    value = value.replace("=", "_")
    value = value.replace(",", "_")
    value = value.replace(" ", "_")
    return value


def normalize_class(cls):
    """
    Harmonisation optionnelle de la classification.
    Adaptée à ton fichier .out harmonisé.
    """
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


def repeatmasker_out_to_gff3(input_out, output_gff, normalize=True):
    n = 0
    skipped = 0

    with open(input_out) as inp, open(output_gff, "w") as out:
        out.write("##gff-version 3\n")

        for line in inp:
            raw = line.rstrip("\n")

            if not raw.strip():
                continue

            fields = raw.split()

            try:
                sw_score = int(fields[0])
            except Exception:
                skipped += 1
                continue

            if len(fields) < 14:
                skipped += 1
                continue

            perc_div = fields[1]
            perc_del = fields[2]
            perc_ins = fields[3]

            seqid = fields[4]
            start = int(fields[5])
            end = int(fields[6])
            left = fields[7]

            strand_raw = fields[8]
            repeat_name = fields[9]
            repeat_class = fields[10]

            if normalize:
                repeat_class = normalize_class(repeat_class)

            # RepeatMasker utilise parfois C pour le brin reverse
            if strand_raw == "C":
                strand = "-"
            elif strand_raw == "+":
                strand = "+"
            else:
                strand = "."

            if start > end:
                start, end = end, start

            n += 1

            feature_id = f"RM_{n:09d}"

            attrs = {
                "ID": feature_id,
                "Name": sanitize_attr(repeat_name),
                "Class": sanitize_attr(repeat_class),
                "Target": sanitize_attr(repeat_name),
                "SW_score": sw_score,
                "PercDiv": perc_div,
                "PercDel": perc_del,
                "PercIns": perc_ins,
                "QueryLeft": sanitize_attr(left),
            }

            attr_text = ";".join(f"{k}={v}" for k, v in attrs.items())

            gff_fields = [
                seqid,
                "RepeatMasker",
                "dispersed_repeat",
                str(start),
                str(end),
                str(sw_score),
                strand,
                ".",
                attr_text,
            ]

            out.write("\t".join(gff_fields) + "\n")

    print(f"Input: {input_out}")
    print(f"Output: {output_gff}")
    print(f"Annotations written: {n}")
    print(f"Lines skipped: {skipped}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert RepeatMasker .out file to GFF3."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input RepeatMasker .out file"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output GFF3 file"
    )

    parser.add_argument(
        "--no-normalize",
        action="store_true",
        help="Do not harmonize repeat class/family names"
    )

    args = parser.parse_args()

    repeatmasker_out_to_gff3(
        input_out=args.input,
        output_gff=args.output,
        normalize=not args.no_normalize
    )


if __name__ == "__main__":
    main()
