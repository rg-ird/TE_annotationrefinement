#!/usr/bin/env python3
"""
hite_pipeline.py

Master pipeline for HiTE / RepeatMasker .out post-processing:
1) Filter annotations by minimum length
2) Remove Simple_repeat / Low_complexity
3) Evaluate overlaps before score filtering
4) Filter overlapping annotations by score
5) Evaluate overlaps after score filtering
6) Reannotate LTR with Inpactor2 first, then optional TEsorter
7) Build RepeatMasker-like .tbl
8) Convert final .out to GFF3
9) Optional graphical outputs: TE/gene density plots

Recommended use:
    python hite_pipeline.py --param config.txt

A config template can be generated with:
    python hite_pipeline.py --write-template config.txt
"""

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


TEMPLATE = """# HiTE TE annotation refinement pipeline parameters
# Lines starting with # are ignored.
# Format: key = value

# Required inputs
rmout = canephora.fna.out
hite_lib = confident_TE.cons.fa
inpactor_lib = inpactor2_library.fa

# Optional, but recommended for step 6
tesorter_tsv = confident_TE.cons.fa.rexdb-plant.cls.tsv

# Output directory
outdir = HITE_PIPELINE_RESULTS

# Directory containing helper scripts
scripts_dir = scripts

# Filtering parameters
min_length = 100
min_overlap_fraction = 0.8
sort = true

# Reannotation parameters
threads = 8
min_identity = 50
min_qcov = 40
blast_tsv =
reuse_blast = false

# Final GFF option
no_normalize_gff = false

# Optional graphical outputs using plot_te_gene_density_from_out.py
# plot = false disables all graphical outputs
plot = false
# Comma-separated list among: level1,ltr_superfamily,gypsy_family,copia_family,all
plot_modes = all
# Prefix for graphical files. If relative, it is written inside outdir.
plot_prefix = te_density
# Required if plot = true, unless plot script can infer sizes. Recommended: genome.fasta.fai
plot_sizes = genome.fna.fai
# Optional gene GFF3. Leave empty to plot only TE densities.
plot_gff = genes.gff3
plot_gene_feature = gene
plot_window = 200000
plot_step = 200000
plot_min_size = 1000000
plot_dpi = 300
plot_chrom_match = auto

# Optional switches
skip_tesorter = false
dry_run = false
"""


def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def str_to_bool(value):
    if isinstance(value, bool):
        return value
    value = str(value).strip().lower()
    if value in {"1", "true", "yes", "y", "on"}:
        return True
    if value in {"0", "false", "no", "n", "off", ""}:
        return False
    raise ValueError(f"Cannot parse boolean value: {value}")


def read_param_file(path):
    params = {}
    with open(path) as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if "=" not in line:
                raise ValueError(
                    f"Invalid parameter line {line_no} in {path}: {line}\n"
                    "Expected format: key = value"
                )

            key, value = line.split("=", 1)
            key = key.strip().replace("-", "_")
            value = value.strip()

            # Remove optional surrounding quotes
            if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
                value = value[1:-1]

            params[key] = value

    return params


def get_param(params, key, default=None, required=False, cast=str):
    if key not in params or params[key] == "":
        if required:
            raise ValueError(f"Missing required parameter in config: {key}")
        return default

    value = params[key]

    if cast is bool:
        return str_to_bool(value)

    return cast(value)


def log_message(message, log_file=None, also_print=True):
    text = f"[{timestamp()}] {message}"
    if also_print:
        print(text)
    if log_file:
        with open(log_file, "a") as log:
            log.write(text + "\n")


def run_cmd(cmd, log_file, dry_run=False, description=None):
    cmd = [str(x) for x in cmd]

    if description:
        log_message(description, log_file)

    log_message("RUN: " + " ".join(shlex.quote(x) for x in cmd), log_file)

    if dry_run:
        log_message("DRY-RUN: command not executed.", log_file)
        return ""

    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    with open(log_file, "a") as log:
        log.write(proc.stdout)
        if proc.stdout and not proc.stdout.endswith("\n"):
            log.write("\n")

    if proc.returncode != 0:
        print(proc.stdout)
        raise RuntimeError(
            f"Command failed with exit code {proc.returncode}:\n"
            + " ".join(shlex.quote(x) for x in cmd)
            + f"\nSee log: {log_file}"
        )

    return proc.stdout


def check_file(path, label):
    if path is None or str(path).strip() == "":
        return
    if not Path(path).exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def check_script(path):
    if not Path(path).exists():
        raise FileNotFoundError(f"Required helper script not found: {path}")


def extract_overlap_percent(text):
    """
    Extracts the line:
    Pourcentage chevauchement  : XX.XXXX%
    from count_TE_overlaps_out.py output.
    """
    patterns = [
        r"Pourcentage\s+chevauchement\s*:\s*([0-9]+(?:\.[0-9]+)?)\s*%",
        r"overlap_percent\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)"
    ]

    for pat in patterns:
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            return float(m.group(1))

    return None


def make_absolute(path):
    if path is None or str(path).strip() == "":
        return None
    return str(Path(path).expanduser().resolve())


def copy_text_file(src, dst):
    with open(src) as inp, open(dst, "w") as out:
        for line in inp:
            out.write(line)


def main():
    parser = argparse.ArgumentParser(
        description="Master pipeline for HiTE / RepeatMasker / Inpactor2 / TEsorter .out refinement."
    )

    parser.add_argument(
        "--param",
        help="Parameter file in key = value format. Recommended."
    )

    parser.add_argument(
        "--write-template",
        help="Write a template parameter file and exit."
    )

    # Minimal command-line overrides
    parser.add_argument("--rmout", help="Input RepeatMasker .out file")
    parser.add_argument("--hite-lib", dest="hite_lib", help="HiTE consensus library FASTA")
    parser.add_argument("--inpactor-lib", dest="inpactor_lib", help="Inpactor2 classified library FASTA")
    parser.add_argument("--tesorter-tsv", dest="tesorter_tsv", help="Optional TEsorter *rexdb-plant.cls.tsv")
    parser.add_argument("--outdir", help="Output directory")
    parser.add_argument("--scripts-dir", dest="scripts_dir", help="Directory containing helper scripts")
    parser.add_argument("--threads", type=int, help="Number of threads")
    parser.add_argument("--min-length", type=int, help="Minimum annotation length. Default: 100")
    parser.add_argument("--min-overlap-fraction", type=float, help="Overlap fraction for score filtering. Default: 0.8")
    parser.add_argument("--min-identity", type=float, help="BLAST minimum identity for reannotation. Default: 50")
    parser.add_argument("--min-qcov", type=float, help="BLAST minimum query coverage for reannotation. Default: 40")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them")
    parser.add_argument("--plot", action="store_true", help="Enable graphical outputs")
    parser.add_argument("--plot-modes", help="Comma-separated plot modes: level1,ltr_superfamily,gypsy_family,copia_family,all")
    parser.add_argument("--plot-prefix", help="Prefix for graphical outputs")
    parser.add_argument("--plot-sizes", help="Genome .fai or two-column sizes file for plots")
    parser.add_argument("--plot-gff", help="Optional GFF3 file for gene density plots")
    parser.add_argument("--plot-gene-feature", help="GFF3 feature type for genes, default: gene")

    args = parser.parse_args()

    if args.write_template:
        Path(args.write_template).write_text(TEMPLATE)
        print(f"Template written: {args.write_template}")
        return

    params = {}
    if args.param:
        params.update(read_param_file(args.param))

    # CLI overrides config file
    for key in [
        "rmout", "hite_lib", "inpactor_lib", "tesorter_tsv", "outdir",
        "scripts_dir", "threads", "min_length", "min_overlap_fraction",
        "min_identity", "min_qcov", "plot_modes", "plot_prefix",
        "plot_sizes", "plot_gff", "plot_gene_feature"
    ]:
        value = getattr(args, key, None)
        if value is not None:
            params[key] = str(value)

    if args.dry_run:
        params["dry_run"] = "true"
    if args.plot:
        params["plot"] = "true"

    rmout = make_absolute(get_param(params, "rmout", required=True))
    hite_lib = make_absolute(get_param(params, "hite_lib", required=True))
    inpactor_lib = make_absolute(get_param(params, "inpactor_lib", required=True))
    tesorter_tsv = make_absolute(get_param(params, "tesorter_tsv", default=None))
    outdir = Path(get_param(params, "outdir", default="HITE_PIPELINE_RESULTS")).expanduser().resolve()
    scripts_dir = Path(get_param(params, "scripts_dir", default="scripts")).expanduser().resolve()

    min_length = get_param(params, "min_length", default=100, cast=int)
    min_overlap_fraction = get_param(params, "min_overlap_fraction", default=0.8, cast=float)
    do_sort = get_param(params, "sort", default=True, cast=bool)

    threads = get_param(params, "threads", default=8, cast=int)
    min_identity = get_param(params, "min_identity", default=50.0, cast=float)
    min_qcov = get_param(params, "min_qcov", default=40.0, cast=float)
    blast_tsv = get_param(params, "blast_tsv", default="")
    reuse_blast = get_param(params, "reuse_blast", default=False, cast=bool)
    skip_tesorter = get_param(params, "skip_tesorter", default=False, cast=bool)

    no_normalize_gff = get_param(params, "no_normalize_gff", default=False, cast=bool)
    dry_run = get_param(params, "dry_run", default=False, cast=bool)

    do_plot = get_param(params, "plot", default=False, cast=bool)
    plot_modes_raw = get_param(params, "plot_modes", default="all")
    plot_prefix = get_param(params, "plot_prefix", default="te_density")
    plot_sizes = make_absolute(get_param(params, "plot_sizes", default=None))
    plot_gff = make_absolute(get_param(params, "plot_gff", default=None))
    plot_gene_feature = get_param(params, "plot_gene_feature", default="gene")
    plot_window = get_param(params, "plot_window", default=200000, cast=int)
    plot_step = get_param(params, "plot_step", default=200000, cast=int)
    plot_min_size = get_param(params, "plot_min_size", default=1000000, cast=int)
    plot_dpi = get_param(params, "plot_dpi", default=300, cast=int)
    plot_chrom_match = get_param(params, "plot_chrom_match", default="auto")

    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / "hite_pipeline.log"

    log_message("Starting HiTE TE annotation refinement pipeline", log_file)
    log_message(f"Output directory: {outdir}", log_file)

    check_file(rmout, "RepeatMasker .out")
    check_file(hite_lib, "HiTE library")
    check_file(inpactor_lib, "Inpactor2 library")
    if tesorter_tsv and not skip_tesorter:
        check_file(tesorter_tsv, "TEsorter TSV")
    if do_plot:
        if plot_sizes:
            check_file(plot_sizes, "Plot sizes / genome FAI")
        if plot_gff:
            check_file(plot_gff, "Plot gene GFF3")

    # Helper scripts
    script_length = scripts_dir / "filter_repeatmasker_min_length.py"
    script_simple = scripts_dir / "filter_repeatmasker_simple_repeat.py"
    script_overlap = scripts_dir / "count_TE_overlaps_out.py"
    script_score = scripts_dir / "filter_annotation_score_out.py"
    script_reannot = scripts_dir / "reannotate_ltr_inpactor2_then_tesorter_harmonized.py"
    script_tbl = scripts_dir / "rebuild_tbl_harmonized.py"
    script_gff = scripts_dir / "repeatmasker_out_to_gff3.py"
    script_plot = scripts_dir / "plot_te_gene_density_from_out.py"

    required_scripts = [
        script_length, script_simple, script_overlap, script_score,
        script_reannot, script_tbl, script_gff
    ]
    if do_plot:
        required_scripts.append(script_plot)

    for s in required_scripts:
        check_script(s)

    # BLAST availability for reannotation script
    if shutil.which("blastn") is None:
        raise EnvironmentError(
            "blastn not found in PATH. Load/install BLAST before running step 6 "
            "(reannotate_ltr_inpactor2_then_tesorter_harmonized.py)."
        )

    # Output files
    step1 = outdir / "01.min_length_filtered.out"
    step2 = outdir / "02.no_simple_repeat.out"
    overlap_before_summary = outdir / "03.overlaps_before_filter.summary.tsv"
    overlap_before_details = outdir / "03.overlaps_before_filter.details.tsv"
    step4 = outdir / "04.score_overlap_filtered.out"
    overlap_after_summary = outdir / "05.overlaps_after_filter.summary.tsv"
    overlap_after_details = outdir / "05.overlaps_after_filter.details.tsv"
    step6 = outdir / "06.reannotated_inpactor2_tesorter_harmonized.out"

    final_out = outdir / "final.filtered.reannotated.out"
    final_tbl = outdir / "final.repeatmasker_like.tbl"
    final_gff = outdir / "final.repeatmasker.gff3"

    # Step 1
    run_cmd([
        sys.executable, script_length,
        "-i", rmout,
        "-o", step1,
        "--min-length", min_length
    ], log_file, dry_run, "STEP 1/8 - Filtering annotations by minimum length")
    log_message(f"STEP 1 done: annotations < {min_length} bp removed.", log_file)

    # Step 2
    run_cmd([
        sys.executable, script_simple,
        "-i", step1,
        "-o", step2
    ], log_file, dry_run, "STEP 2/8 - Removing Simple_repeat / Low_complexity")
    log_message("STEP 2 done: Simple_repeat / Low_complexity removed.", log_file)

    # Step 3
    stdout_before = run_cmd([
        sys.executable, script_overlap,
        "-i", step2,
        "-s", overlap_before_summary,
        "-d", overlap_before_details
    ], log_file, dry_run, "STEP 3/8 - Evaluating overlaps before score filtering")

    overlap_before = extract_overlap_percent(stdout_before) if not dry_run else None
    if overlap_before is not None:
        log_message(f"Pourcentage chevauchement avant filtrage score : {overlap_before:.4f} %", log_file)
    else:
        log_message("Pourcentage chevauchement avant filtrage score : NA", log_file)

    # Step 4
    cmd_step4 = [
        sys.executable, script_score,
        "-i", step2,
        "-o", step4,
        "--min-overlap-fraction", min_overlap_fraction
    ]
    if do_sort:
        cmd_step4.append("--sort")

    run_cmd(cmd_step4, log_file, dry_run, "STEP 4/8 - Filtering overlapping annotations by score")
    log_message(
        f"STEP 4 done: score filtering applied with min-overlap-fraction={min_overlap_fraction}; sort={do_sort}.",
        log_file
    )

    # Step 5
    stdout_after = run_cmd([
        sys.executable, script_overlap,
        "-i", step4,
        "-s", overlap_after_summary,
        "-d", overlap_after_details
    ], log_file, dry_run, "STEP 5/8 - Evaluating overlaps after score filtering")

    overlap_after = extract_overlap_percent(stdout_after) if not dry_run else None
    if overlap_after is not None:
        log_message(f"Pourcentage chevauchement après filtrage score : {overlap_after:.4f} %", log_file)
    else:
        log_message("Pourcentage chevauchement après filtrage score : NA", log_file)

    # Step 6
    reannot_cmd = [
        sys.executable, script_reannot,
        "--rmout", step4,
        "--hite-lib", hite_lib,
        "--inpactor-lib", inpactor_lib,
        "--out", step6,
        "--threads", threads,
        "--min-identity", min_identity,
        "--min-qcov", min_qcov
    ]

    if blast_tsv:
        blast_tsv_path = Path(blast_tsv)
        if not blast_tsv_path.is_absolute():
            blast_tsv_path = outdir / blast_tsv_path
        reannot_cmd.extend(["--blast-tsv", blast_tsv_path])

    if reuse_blast:
        reannot_cmd.append("--reuse-blast")

    if tesorter_tsv and not skip_tesorter:
        reannot_cmd.extend(["--tesorter-tsv", tesorter_tsv])

    run_cmd(
        reannot_cmd,
        log_file,
        dry_run,
        "STEP 6/8 - Reannotating LTR with Inpactor2 first, then optional TEsorter"
    )
    log_message(
        f"STEP 6 done: reannotation completed with min_identity={min_identity}, min_qcov={min_qcov}.",
        log_file
    )

    # Final out copy
    if not dry_run:
        copy_text_file(step6, final_out)
    log_message(f"Final .out written: {final_out}", log_file)

    # Step 7
    run_cmd([
        sys.executable, script_tbl,
        "-i", final_out,
        "-o", final_tbl
    ], log_file, dry_run, "STEP 7/8 - Building RepeatMasker-like .tbl summary")
    log_message(f"STEP 7 done: .tbl produced: {final_tbl}", log_file)

    # Step 8
    cmd_gff = [
        sys.executable, script_gff,
        "-i", final_out,
        "-o", final_gff
    ]
    if no_normalize_gff:
        cmd_gff.append("--no-normalize")

    run_cmd(cmd_gff, log_file, dry_run, "STEP 8/8 - Converting final .out to GFF3")
    log_message(f"STEP 8 done: GFF3 produced: {final_gff}", log_file)

    # Optional graphical outputs
    if do_plot:
        valid_modes = ["level1", "ltr_superfamily", "gypsy_family", "copia_family"]
        if str(plot_modes_raw).strip().lower() == "all":
            plot_modes = valid_modes
        else:
            plot_modes = [m.strip() for m in str(plot_modes_raw).split(",") if m.strip()]
            invalid = [m for m in plot_modes if m not in valid_modes]
            if invalid:
                raise ValueError(
                    "Invalid plot mode(s): " + ", ".join(invalid) +
                    ". Valid modes: " + ", ".join(valid_modes) + ", all"
                )

        base_prefix = Path(plot_prefix)
        if not base_prefix.is_absolute():
            base_prefix = outdir / base_prefix

        log_message("STEP 9 - Producing optional graphical outputs", log_file)
        for mode in plot_modes:
            mode_prefix = f"{base_prefix}_{mode}"
            plot_cmd = [
                sys.executable, script_plot,
                "--out-file", final_out,
                "--mode", mode,
                "--window", plot_window,
                "--step", plot_step,
                "--min-size", plot_min_size,
                "--prefix", mode_prefix,
                "--dpi", plot_dpi,
                "--chrom-match", plot_chrom_match
            ]
            if plot_sizes:
                plot_cmd.extend(["--sizes", plot_sizes])
            if plot_gff:
                plot_cmd.extend(["--gff", plot_gff, "--gene-feature", plot_gene_feature])

            run_cmd(plot_cmd, log_file, dry_run, f"STEP 9 - Plot mode: {mode}")

        log_message(f"STEP 9 done: graphical outputs produced with prefix: {base_prefix}_<mode>", log_file)

    log_message("Pipeline finished successfully.", log_file)

    print("\n========== FINAL OUTPUTS ==========")
    print(f"Final .out : {final_out}")
    print(f"Final .tbl : {final_tbl}")
    print(f"Final GFF3 : {final_gff}")
    if do_plot:
        final_plot_prefix = Path(plot_prefix)
        if not final_plot_prefix.is_absolute():
            final_plot_prefix = outdir / final_plot_prefix
        print(f"Plots      : {final_plot_prefix}_<mode>_*.pdf")
    print(f"Log file   : {log_file}")
    print("===================================")


if __name__ == "__main__":
    main()
