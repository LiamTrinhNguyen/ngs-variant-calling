#!/usr/bin/env bash
# =============================================================================
# Run the Nextflow pipeline and keep ALL of Nextflow's own files inside the
# output folder, so the folder you launch from stays clean:
#
#   <outdir>/                 00_raw … 07_dashboard, logs/   (results)
#   <outdir>/pipeline_info/   nextflow.log, report.html, timeline.html, trace.tsv
#   <outdir>/nextflow/work/   per-task work folders (workDir in nextflow.config)
#   <outdir>/nextflow/cache/  run history used by -resume
#
# WHY a script: Nextflow reads the log and cache locations from environment
#     variables (NXF_LOG_FILE, NXF_CACHE_DIR) before it reads nextflow.config,
#     so they can't be set in the config file. Without them, a plain
#     `nextflow run` writes .nextflow.log and .nextflow/ where you launch it.
#
# Usage (any folder; arguments go straight to `nextflow run`):
#   bash workflow/run_nextflow.sh
#   bash workflow/run_nextflow.sh --outdir /tmp/run --threads 2 -resume
# =============================================================================
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Same output folder that main.nf will use: --outdir if given, else the default
outdir="results_nf"
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do
    case "${args[i]}" in
        --outdir)   outdir="${args[i+1]:-$outdir}" ;;
        --outdir=*) outdir="${args[i]#--outdir=}" ;;
    esac
done

mkdir -p "$outdir/pipeline_info" "$outdir/nextflow"
outdir="$(cd "$outdir" && pwd)"

export NXF_LOG_FILE="$outdir/pipeline_info/nextflow.log"
export NXF_CACHE_DIR="$outdir/nextflow/cache"

exec nextflow run "$here/main.nf" --outdir "$outdir" "$@"
