#!/usr/bin/env python3

import argparse
import subprocess
import shutil
import re
from pathlib import Path
from datetime import datetime


HELP_TEXT = """
TEAR pipeline

Pipeline automatique pour traiter un fichier RepeatMasker .out issu de HiTE,
nettoyer/filtrer les annotations, réannoter les LTR avec Inpactor2 puis TEsorter,
et produire les sorties finales .out, .tbl, .gff3 et éventuellement les graphiques.

Étapes :
  1. filter_repeatmasker_min_length.py
  2. filter_repeatmasker_simple_repeat.py
  3. count_TE_overlaps_out.py avant filtrage score
  4. filter_annotation_score_out.py
  5. count_TE_overlaps_out.py après filtrage score
  6A. nettoyage optionnel de la librairie Inpactor2
  6B. reannotate_ltr_inpactor2_then_tesorter_harmonized.py
  7. rebuild_tbl_harmonized.py
  8. repeatmasker_out_to_gff3.py
  9. graphiques optionnels

Utilisation :
  python tear_pipeline.py --param config.txt
"""


def read_config(path):
    config = {}

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            if "=" not in line:
                continue

            key, value = line.split("=", 1)
            config[key.strip()] = value.strip()

    return config


def get_bool(config, key, default=False):
    value = config.get(key, str(default)).strip().lower()
    return value in {"true", "yes", "1", "y"}


def get_int(config, key, default):
    return int(config.get(key, default))


def get_float(config, key, default):
    return float(config.get(key, default))


def require(config, key):
    if key not in config or not config[key]:
        raise ValueError(f"Paramètre obligatoire manquant dans config.txt : {key}")
    return config[key]


def check_file(path, label):
    if not Path(path).exists():
        raise FileNotFoundError(f"{label} introuvable : {path}")


def log_message(message, log_file):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    text = f"[{timestamp}] {message}"

    print(text)

    with open(log_file, "a") as log:
        log.write(text + "\n")


def run_cmd(cmd, log_file, step_name):
    cmd = [str(x) for x in cmd]

    log_message("", log_file)
    log_message(f"{step_name}", log_file)
    log_message("Commande : " + " ".join(cmd), log_file)

    with open(log_file, "a") as log:
        result = subprocess.run(
            cmd,
            stdout=log,
            stderr=log,
            text=True
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"Erreur pendant {step_name}. Voir le log : {log_file}"
        )

    log_message(f"{step_name} terminé", log_file)


def run_cmd_capture(cmd, log_file, step_name):
    cmd = [str(x) for x in cmd]

    log_message("", log_file)
    log_message(f"{step_name}", log_file)
    log_message("Commande : " + " ".join(cmd), log_file)

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    with open(log_file, "a") as log:
        log.write(result.stdout)
        log.write(result.stderr)

    if result.returncode != 0:
        raise RuntimeError(
            f"Erreur pendant {step_name}. Voir le log : {log_file}"
        )

    log_message(f"{step_name} terminé", log_file)

    return result.stdout + "\n" + result.stderr


def extract_overlap_percentage(text):
    patterns = [
        r"Pourcentage\s+chevauchement\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)",
        r"overlap\s+percentage\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)",
        r"percentage\s+overlap\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)",
        r"([0-9]+(?:\.[0-9]+)?)\s*%"
    ]

    for pattern in patterns:
        m = re.search(pattern, text, re.IGNORECASE)
        if m:
            return float(m.group(1))

    return None


def check_blast_available(log_file):
    blastn = shutil.which("blastn")

    if blastn is None:
        raise RuntimeError(
            "blastn n'est pas accessible. Charge le module BLAST avant de lancer TEAR.\n"
            "Exemple : module load blast"
        )

    log_message(f"blastn accessible : {blastn}", log_file)


def main():
    parser = argparse.ArgumentParser(
        description="TEAR pipeline: HiTE / RepeatMasker / Inpactor2 / TEsorter",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=HELP_TEXT
    )

    parser.add_argument(
        "--param",
        required=True,
        help="Fichier de configuration TEAR config.txt"
    )

    args = parser.parse_args()

    config = read_config(args.param)

    outdir = Path(config.get("outdir", "TEAR_RESULTS"))
    outdir.mkdir(parents=True, exist_ok=True)

    log_file = outdir / "tear_pipeline.log"

    with open(log_file, "w") as log:
        log.write("TEAR pipeline log\n")
        log.write("=================\n")
        log.write(f"Start: {datetime.now()}\n")
        log.write(f"Config file: {args.param}\n\n")

    scripts_dir = Path(config.get("scripts_dir", "."))

    rmout = Path(require(config, "rmout"))
    hite_lib = Path(require(config, "hite_lib"))
    inpactor_lib = Path(require(config, "inpactor_lib"))
    tesorter_tsv = Path(require(config, "tesorter_tsv"))

    check_file(rmout, "RepeatMasker .out initial")
    check_file(hite_lib, "Librairie HiTE")
    check_file(inpactor_lib, "Librairie Inpactor2")
    check_file(tesorter_tsv, "Fichier TEsorter TSV")

    min_len = get_int(config, "min_len", 100)
    min_overlap_fraction = get_float(config, "min_overlap_fraction", 0.8)
    min_identity = get_float(config, "min_identity", 50)
    min_qcov = get_float(config, "min_qcov", 40)
    threads = get_int(config, "threads", 8)

    script_min_length = scripts_dir / config.get(
        "script_min_length",
        "filter_repeatmasker_min_length.py"
    )
    script_simple_repeat = scripts_dir / config.get(
        "script_simple_repeat",
        "filter_repeatmasker_simple_repeat.py"
    )
    script_count_overlaps = scripts_dir / config.get(
        "script_count_overlaps",
        "count_TE_overlaps_out.py"
    )
    script_score_filter = scripts_dir / config.get(
        "script_score_filter",
        "filter_annotation_score_out.py"
    )
    script_reannot = scripts_dir / config.get(
        "script_reannot",
        "reannotate_ltr_inpactor2_then_tesorter_harmonized.py"
    )
    script_tbl = scripts_dir / config.get(
        "script_tbl",
        "rebuild_tbl_harmonized.py"
    )
    script_gff = scripts_dir / config.get(
        "script_gff",
        "repeatmasker_out_to_gff3.py"
    )
    script_plot = scripts_dir / config.get(
        "script_plot",
        "plot_te_gene_density_from_out.py"
    )

    for script in [
        script_min_length,
        script_simple_repeat,
        script_count_overlaps,
        script_score_filter,
        script_reannot,
        script_tbl,
        script_gff,
    ]:
        check_file(script, f"Script {script.name}")

    step1_out = outdir / "01.min_length_filtered.out"
    step2_out = outdir / "02.simple_repeats_removed.out"
    step4_out = outdir / "04.score_overlap_filtered.out"
    step6_out = outdir / "06.reannotated_inpactor2_tesorter.out"

    final_prefix = config.get("final_prefix", "tear")
    final_out = outdir / f"{final_prefix}.final.out"
    final_tbl = outdir / f"{final_prefix}.final.tbl"
    final_gff = outdir / f"{final_prefix}.final.gff3"

    log_message("TEAR pipeline started", log_file)
    log_message(f"Output directory: {outdir}", log_file)

    # 1. Filter minimum length
    run_cmd(
        [
            "python",
            script_min_length,
            "-i",
            rmout,
            "-o",
            step1_out,
            "--min-len",
            min_len,
        ],
        log_file,
        "STEP 1 - Filter RepeatMasker annotations by minimum length"
    )
    log_message(f"STEP 1 done: minimum length filter >= {min_len} bp", log_file)

    # 2. Remove simple repeats / SSR
    run_cmd(
        [
            "python",
            script_simple_repeat,
            "-i",
            step1_out,
            "-o",
            step2_out,
        ],
        log_file,
        "STEP 2 - Remove simple repeats / SSR"
    )
    log_message("STEP 2 done: simple repeats removed", log_file)

    # 3. Count overlaps before score filtering
    before_text = run_cmd_capture(
        [
            "python",
            script_count_overlaps,
            "-i",
            step2_out,
        ],
        log_file,
        "STEP 3 - Count TE overlaps before score filtering"
    )

    before_pct = extract_overlap_percentage(before_text)

    if before_pct is not None:
        log_message(
            f"Pourcentage chevauchement avant filtrage score : {before_pct:.3f} %",
            log_file
        )
    else:
        log_message(
            "Pourcentage chevauchement avant filtrage score : NA",
            log_file
        )

    # 4. Filter annotation by score and overlap fraction
    run_cmd(
        [
            "python",
            script_score_filter,
            "-i",
            step2_out,
            "-o",
            step4_out,
            "--min-overlap-fraction",
            min_overlap_fraction,
            "--sort",
        ],
        log_file,
        "STEP 4 - Filter overlapping annotations by score"
    )
    log_message(
        f"STEP 4 done: score filtering with min-overlap-fraction = {min_overlap_fraction}",
        log_file
    )

    # 5. Count overlaps after score filtering
    after_text = run_cmd_capture(
        [
            "python",
            script_count_overlaps,
            "-i",
            step4_out,
        ],
        log_file,
        "STEP 5 - Count TE overlaps after score filtering"
    )

    after_pct = extract_overlap_percentage(after_text)

    if after_pct is not None:
        log_message(
            f"Pourcentage chevauchement après filtrage score : {after_pct:.3f} %",
            log_file
        )
    else:
        log_message(
            "Pourcentage chevauchement après filtrage score : NA",
            log_file
        )

    # 6A. Optional Inpactor2 library cleaning
    inpactor_lib_for_reannotation = inpactor_lib

    clean_inpactor_lib = get_bool(config, "clean_inpactor_lib", False)

    if clean_inpactor_lib:
        clean_script = Path(
            config.get(
                "clean_inpactor_script",
                "count_ltr_classif_valeurextreme_cleaning.py"
            )
        )

        if not clean_script.is_absolute():
            clean_script = scripts_dir / clean_script

        check_file(clean_script, "Script de nettoyage Inpactor2")

        clean_output = Path(
            config.get(
                "clean_inpactor_output",
                outdir / "inpactor2_library.cleaned.fa"
            )
        )

        if not clean_output.is_absolute():
            clean_output = outdir / clean_output

        run_cmd(
            [
                "python",
                clean_script,
                "-i",
                inpactor_lib,
                "-o",
                clean_output,
            ],
            log_file,
            "STEP 6A - Cleaning Inpactor2 library by removing extreme-size LTR sequences"
        )

        inpactor_lib_for_reannotation = clean_output

        log_message(
            f"STEP 6A done: cleaned Inpactor2 library: {inpactor_lib_for_reannotation}",
            log_file
        )
    else:
        log_message(
            "STEP 6A skipped: clean_inpactor_lib = false",
            log_file
        )

    # 6B. Reannotation with Inpactor2 then TEsorter
    check_blast_available(log_file)

    run_cmd(
        [
            "python",
            script_reannot,
            "--rmout",
            step4_out,
            "--hite-lib",
            hite_lib,
            "--inpactor-lib",
            inpactor_lib_for_reannotation,
            "--tesorter-tsv",
            tesorter_tsv,
            "--out",
            step6_out,
            "--threads",
            threads,
            "--min-identity",
            min_identity,
            "--min-qcov",
            min_qcov,
        ],
        log_file,
        "STEP 6B - Reannotate LTR with Inpactor2 then TEsorter"
    )
    log_message(
        f"STEP 6B done: reannotation min-identity={min_identity}, min-qcov={min_qcov}",
        log_file
    )

    shutil.copyfile(step6_out, final_out)
    log_message(f"Final .out written: {final_out}", log_file)

    # 7. Build RepeatMasker-like tbl
    run_cmd(
        [
            "python",
            script_tbl,
            "-i",
            final_out,
            "-o",
            final_tbl,
        ],
        log_file,
        "STEP 7 - Build RepeatMasker-like .tbl"
    )
    log_message(f"STEP 7 done: .tbl written: {final_tbl}", log_file)

    # 8. Convert to GFF3
    run_cmd(
        [
            "python",
            script_gff,
            "-i",
            final_out,
            "-o",
            final_gff,
        ],
        log_file,
        "STEP 8 - Convert final .out to GFF3"
    )
    log_message(f"STEP 8 done: GFF3 written: {final_gff}", log_file)

    # 9. Optional plots
    make_plots = get_bool(config, "plot", False)

    if make_plots:
        check_file(script_plot, "Script graphique")

        plot_sizes = Path(require(config, "plot_sizes"))
        check_file(plot_sizes, "FAI pour graphiques")

        plot_prefix = config.get("plot_prefix", "te_density")
        plot_modes_raw = config.get("plot_modes", "all")

        if plot_modes_raw == "all":
            plot_modes = [
                "level1",
                "ltr_superfamily",
                "gypsy_family",
                "copia_family",
            ]
        else:
            plot_modes = [
                x.strip()
                for x in plot_modes_raw.split(",")
                if x.strip()
            ]

        plot_window = get_int(config, "plot_window", 200000)
        plot_step = get_int(config, "plot_step", plot_window)
        plot_gff = config.get("plot_gff", "").strip()
        plot_gene_feature = config.get("plot_gene_feature", "gene")

        for mode in plot_modes:
            mode_prefix = outdir / f"{plot_prefix}_{mode}"

            cmd = [
                "python",
                script_plot,
                "--out-file",
                final_out,
                "--sizes",
                plot_sizes,
                "--mode",
                mode,
                "--window",
                plot_window,
                "--step",
                plot_step,
                "--prefix",
                mode_prefix,
            ]

            if plot_gff:
                check_file(plot_gff, "GFF gènes pour graphiques")
                cmd.extend(
                    [
                        "--gff",
                        plot_gff,
                        "--gene-feature",
                        plot_gene_feature,
                    ]
                )

            run_cmd(
                cmd,
                log_file,
                f"STEP 9 - Plot TE/gene density mode={mode}"
            )

        log_message("STEP 9 done: plots generated", log_file)
    else:
        log_message("STEP 9 skipped: plot = false", log_file)

    log_message("", log_file)
    log_message("TEAR pipeline finished successfully", log_file)
    log_message(f"Final OUT : {final_out}", log_file)
    log_message(f"Final TBL : {final_tbl}", log_file)
    log_message(f"Final GFF : {final_gff}", log_file)
    log_message(f"Log file  : {log_file}", log_file)


if __name__ == "__main__":
    main()