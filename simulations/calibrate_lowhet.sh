#!/bin/bash
# Measure the keep fraction categories O and P need. RUNS ON O2.
#
#   bash calibrate_lowhet.sh              # submit the measurement job
#   DRY=1 bash calibrate_lowhet.sh        # resolve every path, submit nothing
#   LOWHET_INLINE=1 bash calibrate_lowhet.sh   # run here instead of via sbatch
#
# WHY THIS EXISTS
#
# O and P are A and E with pi halved. That target fixes the keep fraction rather
# than leaving it to taste, because pi is additive over sites and thinning does
# not move any surviving site's frequency:
#
#     E[pi(k)] = k * pi_neutral + pi_selected
#     pi(k) = 0.5 * pi(1)   =>   k = 0.5 * (1 - pi_selected / pi_neutral)
#
# So there is no search and no iteration -- one measurement of the two components
# per replicate fixes the parameter for the whole arm. This runs that measurement
# over A1-A5 and E1-E5 and prints the mean k per species, which is what goes into
# config/{human,cattle}_lowhet_2Mb_g5t20_r3.yaml.
#
# k DIFFERS BETWEEN THE SPECIES and that is expected, not an error: cattle's
# selected sites carry a different share of its pi than human's do. Pinning one
# shared constant would halve pi in one species and not the other.
#
# WHY IT REPLAYS STAGE 2 RATHER THAN READING THE FILE. On the cattle branch the
# neutral class is laid down IN STAGE 2 by add_neutral, so the raw .full.ts holds
# only the 5.6e-9 DFE mutations and would report a selected share near 1.0 -- an
# "infeasible" verdict that is an artefact of measuring the wrong object.
# helpers/measure_pi_components.py imports the pipeline's own add_neutral for
# exactly this reason.
#
# IF IT SAYS INFEASIBLE. k comes out negative when the selected class alone
# already carries half of that species' pi, i.e. no amount of neutral thinning
# halves it. The script reports the numbers rather than clamping; that is a
# result about the simulation, and the arm has to be re-specified, not rounded.

set -uo pipefail
REPO="${REPO:-/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake}"
SCRATCH="${SCRATCH:-/n/scratch/users/n/njc12/snakemake}"
OUT_TSV="${OUT_TSV:-$REPO/helpers/pi_components_lowhet.tsv}"
PYTHON="${PYTHON:-/home/njc12/miniconda3/envs/coloc_sims/bin/python}"
REPS="${REPS:-A1 A2 A3 A4 A5 E1 E2 E3 E4 E5}"
DRY="${DRY:-}"

# The same two roots the submit scripts search, most-specific first: A2-A5 and
# E2-E5 ran their own stage 1 under the cmaf arm; A1 and E1 reuse the x35 pilot's.
STAGE1_ROOTS=(
    "$SCRATCH/simulations_round_3_2Mb_g5t20_cmaf_0.001"
    "$SCRATCH/simulations_round_3_2Mb"
)

cd "$REPO" || { echo "ERROR: no repo at $REPO" >&2; exit 1; }

# ------------------------------------------------------------------ pre-flight
[[ -f helpers/measure_pi_components.py ]] || {
    echo "ERROR: helpers/measure_pi_components.py missing -- rsync the repo to O2." >&2
    exit 1; }

seed_of() {  # the parent category's seed; O and P run at these, which is the point
    local n="${1:1}"
    case "${1:0:1}" in
        A) echo "1${n}" ;;
        E) echo "5${n}" ;;
        *) echo "ERROR: calibrate only A and E, got '$1'" >&2; return 1 ;;
    esac
}

declare -A TS_OF SEED_OF
for ID in $REPS; do
    SEED=$(seed_of "$ID") || exit 1
    TS=""
    for ROOT in "${STAGE1_ROOTS[@]}"; do
        for CAND in "$ROOT/${ID}/stage1"/*.ts; do
            [[ -e "$CAND" ]] && { TS="$CAND"; break 2; }
        done
    done
    [[ -n "$TS" ]] || {
        echo "ERROR: no stage-1 tree sequence for ${ID} in any of:" >&2
        for ROOT in "${STAGE1_ROOTS[@]}"; do echo "       $ROOT/${ID}/stage1" >&2; done
        echo "Nothing measured." >&2; exit 1; }
    # The seed is embedded in the cattle filename, so a mismatch here means the
    # donor is not the replicate it claims to be -- which would calibrate O
    # against a genome it will not actually run on.
    if [[ "${ID:0:1}" == "E" && "$TS" != *"seed_${SEED}."* ]]; then
        echo "ERROR: ${ID} resolves to $TS, which is not seed ${SEED}." >&2
        exit 1
    fi
    TS_OF["$ID"]="$TS"; SEED_OF["$ID"]="$SEED"
done

echo "Pre-flight OK: $(echo $REPS | wc -w) replicates."
for ID in $REPS; do printf "  %-3s seed %-3s %s\n" "$ID" "${SEED_OF[$ID]}" "${TS_OF[$ID]}"; done
echo "  target: pi halved     output: $OUT_TSV"
echo
[[ -n "$DRY" ]] && { echo "DRY=1 -- nothing submitted."; exit 0; }

# ---------------------------------------------------------------- the measuring
measure_all() {
    : > "$OUT_TSV"
    for ID in $REPS; do
        local sp; sp=$([[ "${ID:0:1}" == "E" ]] && echo cattle || echo human)
        echo "=== $ID ($sp) ==="
        if [[ "$sp" == "cattle" ]]; then
            "$PYTHON" helpers/measure_pi_components.py \
                --ts_file "${TS_OF[$ID]}" --species cattle \
                --Q_scaling 0.01 --handoff_ticks 2400 --deep_Q_scaling 1 \
                --seed "${SEED_OF[$ID]}" --label "$ID" --out_tsv "$OUT_TSV" || return 1
        else
            "$PYTHON" helpers/measure_pi_components.py \
                --ts_file "${TS_OF[$ID]}" --species human \
                --Q_scaling 10 --already_includes_neutral True \
                --seed "${SEED_OF[$ID]}" --label "$ID" --out_tsv "$OUT_TSV" || return 1
        fi
    done

    echo
    echo "================================================================"
    awk -F'\t' 'NR == 1 {
            for (i = 1; i <= NF; i++) col[$i] = i; next
        }
        $col["keep_fraction_for_half_pi"] == "NA" { bad[$col["species"]]++; next }
        {
            sp = $col["species"]
            k[sp] += $col["keep_fraction_for_half_pi"]; n[sp]++
            share[sp] += $col["selected_share"]
        }
        END {
            printf "%-8s %3s  %-22s %-22s %s\n", "species", "n", \
                   "mean selected pi share", "mean keep fraction", "pin as"
            for (sp in n) {
                mk = k[sp] / n[sp]
                printf "%-8s %3d  %-22.4f %-22.6f neutral_keep_fraction: %.3f\n", \
                       sp, n[sp], share[sp] / n[sp], mk, mk
            }
            for (sp in bad)
                printf "%-8s %3d replicate(s) INFEASIBLE -- see stderr above; do NOT pin a value\n", \
                       sp, bad[sp]
        }' "$OUT_TSV"
    echo "================================================================"
    echo "Per-replicate table: $OUT_TSV"
    echo "Pin the rounded value in:"
    echo "  config/human_lowhet_2Mb_g5t20_r3.yaml    (category O)"
    echo "  config/cattle_lowhet_2Mb_g5t20_r3.yaml   (category P)"
    echo "Both submit scripts refuse to launch O/P until the key is set, so an"
    echo "un-pinned run cannot silently reproduce A and E."
}

if [[ -n "${LOWHET_INLINE:-}" || -n "${SLURM_JOB_ID:-}" ]]; then
    measure_all
    exit $?
fi

# Sized against stage2_split_pheno, which loads the same tree sequences: that rule
# measured 369 MB peak at L = 2 Mb. diversity() over ~18k samples is the cost here,
# and 10 replicates run serially, so the walltime is the generous part.
JID=$(sbatch --parsable --job-name=lowhet_calibration \
    --partition=short --time=4:00:00 --mem=32G --cpus-per-task=4 \
    --output="$REPO/lowhet_calibration_%j.out" \
    --error="$REPO/lowhet_calibration_%j.err" \
    --export=ALL,LOWHET_INLINE=1,REPS="$REPS",OUT_TSV="$OUT_TSV" \
    "$0")
echo "Submitted calibration job $JID"
echo "  tail -f $REPO/lowhet_calibration_${JID}.out"
