#!/usr/bin/env python3

import argparse

def filter_repeatmasker_out(input_file, output_file, remove_classes):
    removed = 0
    kept = 0
    header_lines = 0

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            stripped = line.strip()

            # Garder les lignes vides et l'en-tête
            if not stripped or stripped.startswith("SW") or stripped.startswith("score"):
                outfile.write(line)
                header_lines += 1
                continue

            cols = stripped.split()

            # Les vraies lignes RepeatMasker ont au moins 14 colonnes
            if len(cols) < 14:
                outfile.write(line)
                continue

            repeat_class = cols[10]

            if repeat_class in remove_classes:
                removed += 1
                continue
            else:
                outfile.write(line)
                kept += 1

    print("Résumé du filtrage")
    print(f"Fichier d'entrée : {input_file}")
    print(f"Fichier de sortie : {output_file}")
    print(f"Annotations conservées : {kept}")
    print(f"Annotations retirées : {removed}")
    print(f"Classes retirées : {', '.join(remove_classes)}")


def main():
    parser = argparse.ArgumentParser(
        description="Retire les annotations Simple_repeat et Low_complexity d'un fichier RepeatMasker .out"
    )

    parser.add_argument("-i", "--input", required=True, help="Fichier RepeatMasker .out en entrée")
    parser.add_argument("-o", "--output", required=True, help="Fichier .out filtré en sortie")

    parser.add_argument(
        "--remove",
        nargs="+",
        default=["Simple_repeat", "Low_complexity"],
        help="Classes à retirer. Défaut: Simple_repeat Low_complexity"
    )

    args = parser.parse_args()

    filter_repeatmasker_out(args.input, args.output, set(args.remove))


if __name__ == "__main__":
    main()
