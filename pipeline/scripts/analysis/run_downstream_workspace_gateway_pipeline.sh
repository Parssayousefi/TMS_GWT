#!/bin/bash
#
# Run downstream workspace and gateway/broadcaster analyses after PhiID.
# 
# Included:
#   - 05_group_average_workspace.py for HCP, all control conditions, and all drug conditions
#   - 06_identify_gateways_broadcasters.py for the same set of conditions
#
# Explicitly excluded:
#   - dataset-specific control-defined workspaces via --workspace_dataset <dataset>
#   - cross-dataset pooled placebo/baseline comparison via 05 --compare
#   - pooled placebo/baseline gateway/broadcaster analysis via 06 --pooled
#
# Usage:
#   bash scripts/analysis/run_downstream_workspace_gateway_pipeline.sh
#   bash scripts/analysis/run_downstream_workspace_gateway_pipeline.sh --dry-run
#   bash scripts/analysis/run_downstream_workspace_gateway_pipeline.sh --max-jobs 8
#
# Environment overrides:
#   MAX_JOBS=8
#   PYTHON_BIN=python
#   OVERWRITE_06=0

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash scripts/analysis/run_downstream_workspace_gateway_pipeline.sh [options]

Options:
  --dry-run          Print commands without running them
  --max-jobs N       Maximum number of parallel jobs per phase
  --python CMD       Python executable to use
  --no-overwrite-06  Do not pass --overwrite to script 06
  -h, --help         Show this help

Notes:
  - Script 05 is always recomputed because it has no overwrite flag.
  - Script 06 uses --overwrite by default so its outputs stay in sync with fresh 05 outputs.
  - basel_mmm is intentionally not included because its dataset config file is missing in this repo.
EOF
}

DRY_RUN=0
MAX_JOBS="${MAX_JOBS:-4}"
PYTHON_BIN="${PYTHON_BIN:-python}"
OVERWRITE_06="${OVERWRITE_06:-1}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --max-jobs)
            MAX_JOBS="$2"
            shift 2
            ;;
        --python)
            PYTHON_BIN="$2"
            shift 2
            ;;
        --no-overwrite-06)
            OVERWRITE_06=0
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if ! [[ "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: --max-jobs must be a positive integer" >&2
    exit 1
fi

resolve_python() {
    if python_works "$PYTHON_BIN"; then
        return 0
    fi

    if [[ "$PYTHON_BIN" == "python" ]]; then
        local candidate
        for candidate in \
            /mnt/c/Miniforge3/envs/gwt/python.exe \
            /mnt/c/Users/YousefiP/.conda/envs/gwt/python.exe \
            python.exe \
            /mnt/c/Miniforge3/python.exe \
            /c/Miniforge3/envs/gwt/python.exe \
            /c/Users/YousefiP/.conda/envs/gwt/python.exe \
            /c/Miniforge3/python.exe \
            python3; do
            if python_works "$candidate"; then
                PYTHON_BIN="$candidate"
                return 0
            fi
        done
    fi

    return 1
}

python_works() {
    local candidate="$1"

    if command -v "$candidate" >/dev/null 2>&1; then
        :
    elif [[ "$candidate" == /* ]] && [[ -e "$candidate" ]]; then
        :
    else
        return 1
    fi

    "$candidate" -c "import numpy, pandas, yaml, scipy" >/dev/null 2>&1
}

if ! resolve_python; then
    echo "Error: Python executable not found: $PYTHON_BIN" >&2
    echo "Try: --python /mnt/c/Miniforge3/envs/gwt/python.exe" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="$PROJECT_DIR/logs/downstream_workspace_gateway/$TIMESTAMP"

mkdir -p "$LOG_DIR"
cd "$PROJECT_DIR"

build_cmd() {
    local quoted=()
    local arg
    for arg in "$@"; do
        quoted+=("$(printf '%q' "$arg")")
    done
    local joined
    joined="${quoted[*]}"
    printf '%s' "$joined"
}

append_job() {
    local -n target_array="$1"
    local job_name="$2"
    shift 2
    local cmd
    cmd="$(build_cmd "$@")"
    target_array+=("${job_name}|${cmd}")
}

append_05_job() {
    local array_name="$1"
    local job_name="$2"
    shift 2
    append_job "$array_name" "$job_name" \
        "$PYTHON_BIN" "scripts/analysis/05_group_average_workspace.py" "$@"
}

append_06_job() {
    local array_name="$1"
    local job_name="$2"
    shift 2

    if [[ "$OVERWRITE_06" == "1" ]]; then
        append_job "$array_name" "$job_name" \
            "$PYTHON_BIN" "scripts/analysis/06_identify_gateways_broadcasters.py" "$@" "--overwrite"
    else
        append_job "$array_name" "$job_name" \
            "$PYTHON_BIN" "scripts/analysis/06_identify_gateways_broadcasters.py" "$@"
    fi
}

run_parallel_phase() {
    local phase_name="$1"
    local -n jobs_ref="$2"

    echo "================================================================================"
    echo "$phase_name"
    echo "================================================================================"
    echo "Jobs: ${#jobs_ref[@]}"
    echo "Max parallel jobs: $MAX_JOBS"
    echo ""

    if [[ ${#jobs_ref[@]} -eq 0 ]]; then
        echo "No jobs queued."
        echo ""
        return 0
    fi

    if [[ "$DRY_RUN" == "1" ]]; then
        local spec
        for spec in "${jobs_ref[@]}"; do
            echo "${spec#*|}"
        done
        echo ""
        return 0
    fi

    printf '%s\0' "${jobs_ref[@]}" | xargs -0 -I{} -P "$MAX_JOBS" bash -lc '
        set -euo pipefail
        spec="$1"
        project_dir="$2"
        log_dir="$3"

        name="${spec%%|*}"
        cmd="${spec#*|}"
        log_file="$log_dir/${name}.log"

        echo "[START] $name"
        if (
            cd "$project_dir"
            echo "[COMMAND] $cmd"
            eval "$cmd"
        ) >"$log_file" 2>&1; then
            echo "[DONE]  $name"
            echo "        log: $log_file"
        else
            status=$?
            echo "[FAIL]  $name (exit $status)" >&2
            echo "        log: $log_file" >&2
            exit "$status"
        fi
    ' _ "{}" "$PROJECT_DIR" "$LOG_DIR"

    echo ""
}

declare -a PHASE_05_JOBS=()
declare -a PHASE_06_JOBS=()

# ---------------------------------------------------------------------------
# Phase 05: group-average workspace jobs
# ---------------------------------------------------------------------------

# HCP already done but included for completeness (05 will overwrite)
append_05_job PHASE_05_JOBS "05_hcp_rest" --dataset hcp

# Berlin Ketamine (2 conditions)
append_05_job PHASE_05_JOBS "05_berlin_ketamine_placebo" --dataset berlin_ketamine
append_05_job PHASE_05_JOBS "05_berlin_ketamine_ketamine" --dataset berlin_ketamine --task ketamine

# Basel LAM (4 conditions)
append_05_job PHASE_05_JOBS "05_basel_lam_placebo" --dataset basel_lam
append_05_job PHASE_05_JOBS "05_basel_lam_lsd" --dataset basel_lam --task lsd
append_05_job PHASE_05_JOBS "05_basel_lam_amph" --dataset basel_lam --task amph
append_05_job PHASE_05_JOBS "05_basel_lam_mdma" --dataset basel_lam --task mdma

# Basel LPM (4 conditions)
append_05_job PHASE_05_JOBS "05_basel_lpm_placebo" --dataset basel_lpm
append_05_job PHASE_05_JOBS "05_basel_lpm_lsd" --dataset basel_lpm --task lsd
append_05_job PHASE_05_JOBS "05_basel_lpm_psilocybin" --dataset basel_lpm --task psilocybin
append_05_job PHASE_05_JOBS "05_basel_lpm_mescaline" --dataset basel_lpm --task mescaline

# london_lsd: excluded — no PhiID data (not yet through CONN pipeline)
# missouri_psilo: excluded — not in current rollout

# ---------------------------------------------------------------------------
# Phase 06: gateway/broadcaster jobs using default script behavior
# - no dataset-specific workspace_dataset override
# - no pooled analysis
# ---------------------------------------------------------------------------

# HCP (baseline workspace reference)
append_06_job PHASE_06_JOBS "06_hcp_rest" --dataset hcp

# Berlin Ketamine — drug vs placebo (HCP workspace is default)
append_06_job PHASE_06_JOBS "06_berlin_ketamine_placebo" --dataset berlin_ketamine
append_06_job PHASE_06_JOBS "06_berlin_ketamine_ketamine" --dataset berlin_ketamine --task ketamine

# Basel LAM — 3 drugs vs placebo
append_06_job PHASE_06_JOBS "06_basel_lam_placebo" --dataset basel_lam
append_06_job PHASE_06_JOBS "06_basel_lam_lsd" --dataset basel_lam --task lsd
append_06_job PHASE_06_JOBS "06_basel_lam_amph" --dataset basel_lam --task amph
append_06_job PHASE_06_JOBS "06_basel_lam_mdma" --dataset basel_lam --task mdma

# Basel LPM — 3 drugs vs placebo
append_06_job PHASE_06_JOBS "06_basel_lpm_placebo" --dataset basel_lpm
append_06_job PHASE_06_JOBS "06_basel_lpm_lsd" --dataset basel_lpm --task lsd
append_06_job PHASE_06_JOBS "06_basel_lpm_psilocybin" --dataset basel_lpm --task psilocybin
append_06_job PHASE_06_JOBS "06_basel_lpm_mescaline" --dataset basel_lpm --task mescaline

# london_lsd: excluded — no PhiID data (not yet through CONN pipeline)
# missouri_psilo: excluded — not in current rollout

echo "================================================================================"
echo "DOWNSTREAM WORKSPACE + GATEWAY/BROADCASTER PIPELINE"
echo "================================================================================"
echo "Project directory: $PROJECT_DIR"
echo "Python: $PYTHON_BIN"
echo "Max parallel jobs: $MAX_JOBS"
echo "Dry run: $DRY_RUN"
echo "06 overwrite enabled: $OVERWRITE_06"
echo "Logs: $LOG_DIR"
echo ""
echo "This script runs:"
echo "  - Phase 1: all non-pooled script 05 jobs"
echo "  - Phase 2: all non-pooled script 06 jobs using default HCP workspace behavior"
echo ""

run_parallel_phase "PHASE 1/2: GROUP-AVERAGE WORKSPACE ANALYSES (SCRIPT 05)" PHASE_05_JOBS
run_parallel_phase "PHASE 2/2: GATEWAY/BROADCASTER ANALYSES (SCRIPT 06)" PHASE_06_JOBS

if [[ "$DRY_RUN" == "1" ]]; then
    echo "Dry run complete."
    exit 0
fi

echo "================================================================================"
echo "PIPELINE COMPLETE"
echo "================================================================================"
echo "Logs saved to: $LOG_DIR"
echo ""
