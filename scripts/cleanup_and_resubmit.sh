#!/bin/bash
# =============================================================================
# cleanup_and_resubmit.sh
# Fix corrupt VCFs / .tbi files caused by the 04_make_vcf_var.sh race condition,
# then identify which (rep × cov × freq) combos need resubmission.
#
# Run on the cluster login node (not via Slurm) after loading samtools/tabix.
# Usage: bash scripts/cleanup_and_resubmit.sh [--dry-run]
# =============================================================================
set -euo pipefail

# WORK: prefer explicit env var, then the hardcoded config value (matches generated configs)
WORK="${WORK:-/home/tbellagio/scratch/visor_freqk}"
# VCFs may live on a different mount (carnegie nobackup); detect from script location
_SCRIPT_WORK="$(cd "$(dirname "$0")/.." && pwd)"
VCF_ROOT="${_SCRIPT_WORK}/data/vcf/del"
RESULTS_ROOT="${WORK}/results/del"
SCRIPTS="${WORK}/scripts"
DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

SIZES=(100bp 500bp 1kb 5kb 10kb)

echo "============================================================"
echo " Phase 1: Validate and repair VCFs in ${VCF_ROOT}"
echo "============================================================"

N_CORRUPT=0
for rep_dir in "${VCF_ROOT}"/rep*/var; do
    [[ -d "${rep_dir}" ]] || continue
    rep=$(basename "$(dirname "${rep_dir}")")   # e.g. rep9

    for sz in "${SIZES[@]}"; do
        vcf="${rep_dir}/del_${sz}.vcf.gz"
        tbi="${vcf}.tbi"

        [[ -f "${vcf}" ]] || continue

        # Count non-header lines
        n_data=$(zcat "${vcf}" 2>/dev/null | grep -cv "^#" || true)

        if [[ "${n_data}" -ne 1 ]]; then
            echo "[CORRUPT] ${vcf}  (data lines=${n_data}, expected=1)"
            N_CORRUPT=$((N_CORRUPT + 1))
            if [[ "${DRY_RUN}" -eq 0 ]]; then
                rm -f "${vcf}" "${tbi}" "${rep_dir}/del_${sz}.lock"
                echo "  → Deleted corrupt VCF + tbi"
            else
                echo "  → [DRY-RUN] would delete"
            fi
        else
            # Data line count OK — still check .tbi exists and is non-empty
            if [[ ! -s "${tbi}" ]]; then
                echo "[MISSING TBI] ${vcf}"
                N_CORRUPT=$((N_CORRUPT + 1))
                if [[ "${DRY_RUN}" -eq 0 ]]; then
                    tabix -f "${vcf}"
                    echo "  → Regenerated .tbi"
                else
                    echo "  → [DRY-RUN] would regenerate .tbi"
                fi
            fi
        fi
    done
done

echo ""
if [[ "${N_CORRUPT}" -eq 0 ]]; then
    echo "No corrupt VCFs found."
else
    echo "Fixed ${N_CORRUPT} corrupt/missing file(s)."
fi

echo ""
echo "============================================================"
echo " Phase 2: Find missing results (need resubmission)"
echo "============================================================"

COVERAGES=(10 20 50)
SV_FREQS=(0.10 0.30 0.50 0.70 0.90)
N_SAMPLES=231
ERR_LABEL=001   # from ERROR_RATE=0.001 → "0.001".split('.')[-1] = "001"
K=31
# Optional: restrict to a subset of rep IDs (space-separated), e.g. "rep1 rep2 rep3"
# Leave empty to scan all reps found in the VCF dir.
REP_FILTER="${REP_FILTER:-}"

# Build cov/err label mappings (must match config_sv_var_deletions.sh logic)
# ERROR_RATE=0.001 → "0.001".split('.')[-1] → "001"
declare -A COV_LABEL_MAP=([10]="cov10_err001" [20]="cov20_err001" [50]="cov50_err001")

MISSING_JOBS=()

# Scan all reps that have a VCF dir (pipeline got at least to step 4)
for vcf_rep_dir in "${VCF_ROOT}"/rep*/var; do
    [[ -d "${vcf_rep_dir}" ]] || continue
    rep=$(basename "$(dirname "${vcf_rep_dir}")")
    # If REP_FILTER is set, skip reps not in the filter list
    if [[ -n "${REP_FILTER}" ]]; then
        [[ " ${REP_FILTER} " == *" ${rep} "* ]] || continue
    fi
    rep_dir="${RESULTS_ROOT}/${rep}/var"

    for cov in "${COVERAGES[@]}"; do
        cov_label="${COV_LABEL_MAP[$cov]}"

        for freq in "${SV_FREQS[@]}"; do
            # freq label matches 05_freqk_var.sh: f$(python3 -c "print(round(SV_FREQ*100))")
            # 0.10 → f10, 0.30 → f30, 0.50 → f50, 0.70 → f70, 0.90 → f90
            freq_int=$(awk -v f="${freq}" 'BEGIN{printf "%.0f", f*100}')
            freq_label="f${freq_int}"

            all_done=1
            missing_sizes=()

            for sz in "${SIZES[@]}"; do
                result_dir="${rep_dir}/${cov_label}/${sz}/n${N_SAMPLES}/${freq_label}/k31"
                # Output file: var_del_<sz>_n<N>_f<F>_err<E>.allele_frequencies.k<K>.tsv
                af_file="${result_dir}/var_del_${sz}_n${N_SAMPLES}_${freq_label}_err${ERR_LABEL}.allele_frequencies.k${K}.tsv"
                if [[ ! -s "${af_file}" ]]; then
                    all_done=0
                    missing_sizes+=("${sz}")
                fi
            done

            if [[ "${all_done}" -eq 0 ]]; then
                MISSING_JOBS+=("${rep}|cov${cov}|freq${freq}|missing:${missing_sizes[*]}")
            fi
        done
    done
done  # rep loop

if [[ "${#MISSING_JOBS[@]}" -eq 0 ]]; then
    echo "All results present — nothing to resubmit!"
else
    echo "Found ${#MISSING_JOBS[@]} (rep × cov × freq) combos with missing results:"
    echo ""
    printf "  %-10s %-8s %-8s %s\n" "REP" "COV" "FREQ" "MISSING_SIZES"
    printf "  %-10s %-8s %-8s %s\n" "---" "---" "----" "-------------"
    for job in "${MISSING_JOBS[@]}"; do
        IFS='|' read -r rep cov freq miss <<< "${job}"
        printf "  %-10s %-8s %-8s %s\n" "${rep}" "${cov}" "${freq}" "${miss}"
    done

    echo ""
    echo "============================================================"
    echo " Phase 3: Resubmit missing jobs"
    echo "============================================================"
    echo ""
    echo "To resubmit ALL missing jobs, run:"
    echo ""

    # Group by rep and collect unique (cov, freq) pairs
    declare -A REP_COVS
    declare -A REP_FREQS
    for job in "${MISSING_JOBS[@]}"; do
        IFS='|' read -r rep cov freq miss <<< "${job}"
        cov_num="${cov#cov}"
        freq_num="${freq#freq}"
        REP_COVS["${rep}"]+="${cov_num} "
        REP_FREQS["${rep}"]+="${freq_num} "
    done

    for rep in $(echo "${!REP_COVS[@]}" | tr ' ' '\n' | sort -V); do
        covs=$(echo "${REP_COVS[$rep]}" | tr ' ' '\n' | sort -nu | tr '\n' ' ')
        freqs=$(echo "${REP_FREQS[$rep]}" | tr ' ' '\n' | sort -nu | tr '\n' ' ')
        echo "  python ${SCRIPTS}/launch_experiment_var.py \\"
        echo "      --sv-type DEL \\"
        echo "      --coverage ${covs}\\"
        echo "      --sv-freq ${freqs}\\"
        echo "      --n-samples ${N_SAMPLES} --sizes 100bp 500bp 1kb 5kb 10kb \\"
        echo "      --rep-ids ${rep}"
        echo ""
    done

    if [[ "${DRY_RUN}" -eq 0 ]]; then
        read -r -p "Resubmit all missing jobs now? [y/N] " answer
        if [[ "${answer}" =~ ^[Yy]$ ]]; then
            for rep in $(echo "${!REP_COVS[@]}" | tr ' ' '\n' | sort -V); do
                covs=$(echo "${REP_COVS[$rep]}" | tr ' ' '\n' | sort -nu | tr '\n' ' ')
                freqs=$(echo "${REP_FREQS[$rep]}" | tr ' ' '\n' | sort -nu | tr '\n' ' ')
                echo "Submitting ${rep} …"
                # shellcheck disable=SC2086
                python3 "${SCRIPTS}/launch_experiment_var.py" \
                    --sv-type DEL \
                    --coverage ${covs} \
                    --sv-freq ${freqs} \
                    --n-samples "${N_SAMPLES}" \
                    --sizes 100bp 500bp 1kb 5kb 10kb \
                    --rep-ids "${rep}"
            done
        else
            echo "Skipped resubmission."
        fi
    fi
fi
