#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np


def identify_and_remove_outliers(
    fasta_file,
    output_file,
    iqr_factor=1.5,
    pdf_prefix="LTR_Sequence_Sizes"
):
    classification_sizes = defaultdict(list)
    sequences_by_classification = defaultdict(list)

    for record in SeqIO.parse(fasta_file, "fasta"):
        classification = record.id.split("#")[1] if "#" in record.id else "Unknown"
        classification_sizes[classification].append(len(record.seq))
        sequences_by_classification[classification].append(record)

    outliers = defaultdict(list)

    for classification, sizes in classification_sizes.items():
        sizes_np = np.array(sizes)

        q1 = np.percentile(sizes_np, 25)
        q3 = np.percentile(sizes_np, 75)
        iqr = q3 - q1

        lower_bound = q1 - iqr_factor * iqr
        upper_bound = q3 + iqr_factor * iqr

        for i, size in enumerate(sizes):
            if size < lower_bound or size > upper_bound:
                outliers[classification].append(
                    sequences_by_classification[classification][i].id
                )

    filtered_sequences = []
    filtered_classification_sizes = defaultdict(list)

    removed_count = 0

    for classification, records in sequences_by_classification.items():
        for record in records:
            if record.id not in outliers[classification]:
                filtered_sequences.append(record)
                filtered_classification_sizes[classification].append(len(record.seq))
            else:
                removed_count += 1

    SeqIO.write(filtered_sequences, output_file, "fasta")

    # Boxplot avant filtrage
    plt.figure(figsize=(12, 8))
    plt.boxplot(
        classification_sizes.values(),
        labels=classification_sizes.keys()
    )
    plt.xlabel("Classification")
    plt.ylabel("Taille des séquences")
    plt.title("Distribution des tailles des séquences par classification avant filtrage")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(f"{pdf_prefix}_Before_Filtering.pdf")
    plt.close()

    # Boxplot après filtrage
    plt.figure(figsize=(12, 8))
    plt.boxplot(
        filtered_classification_sizes.values(),
        labels=filtered_classification_sizes.keys()
    )
    plt.xlabel("Classification")
    plt.ylabel("Taille des séquences")
    plt.title("Distribution des tailles des séquences par classification après filtrage")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(f"{pdf_prefix}_After_Filtering.pdf")
    plt.close()

    total_count = sum(len(v) for v in sequences_by_classification.values())
    kept_count = len(filtered_sequences)

    print("===== Résumé nettoyage librairie Inpactor2 =====")
    print(f"Fichier entrée     : {fasta_file}")
    print(f"Fichier sortie     : {output_file}")
    print(f"IQR factor         : {iqr_factor}")
    print(f"Séquences totales  : {total_count}")
    print(f"Séquences gardées  : {kept_count}")
    print(f"Séquences retirées : {removed_count}")
    print(f"PDF avant          : {pdf_prefix}_Before_Filtering.pdf")
    print(f"PDF après          : {pdf_prefix}_After_Filtering.pdf")


def main():
    parser = argparse.ArgumentParser(
        description="Remove extreme-size LTR sequences from a classified FASTA library."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input classified FASTA file"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output filtered FASTA file"
    )

    parser.add_argument(
        "--iqr-factor",
        type=float,
        default=1.5,
        help="IQR multiplier for outlier detection. Default: 1.5"
    )

    parser.add_argument(
        "--pdf-prefix",
        default="LTR_Sequence_Sizes",
        help="Prefix for output PDF boxplots"
    )

    args = parser.parse_args()

    identify_and_remove_outliers(
        fasta_file=args.input,
        output_file=args.output,
        iqr_factor=args.iqr_factor,
        pdf_prefix=args.pdf_prefix
    )


if __name__ == "__main__":
    main()