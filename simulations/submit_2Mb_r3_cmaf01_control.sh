#!/bin/bash
# The five-replicate CONTROL arm: all three MAF floors at 0.01.
#
#   causal_min_maf 0.01   causative variants drawn from MAF >= 0.01
#   fm_min_maf     0.01   fine-mapped from the same floor
#   min_maf        0.01   and entered into the plink2 GLM from the same floor
#   ld_ctrl        0.75   dap-g signal-cluster r2 threshold
#   gwas_scaling 5, gtex_scaling 20 (from config/*_2Mb_g5t20_r3.yaml)
#
# This is the baseline the two causal-0.001 arms are measured against:
#   ..._g5t20_cmaf_0.001              causal 0.001, fm 0.01,  glm 0.01
#   ..._g5t20_cmaf_0.001_fm_0.001     causal 0.001, fm 0.001, glm 0.001
#   ..._g5t20_cmaf_0.01_fm_0.01       <- this one, all three at 0.01
# Each at A1-A5 + E1-E5, so every arm is 250 GWAS traits rather than 50.
#
# WHY A1 AND E1 ARE RERUN even though simulations_round_3_2Mb_g5t20_maf01 already
# holds them at exactly these settings: that pair predates the params file, so
# summarize_coloc.py reverse-engineers their settings and records
# params_source=inferred -- and params_source is part of its grouping key. Left
# as they are, they would form a separate one-replicate row beside a
# four-replicate row describing the SAME experimental condition. Rerunning them
# here puts all five replicates on one row. Their numbers should reproduce.
#
# STAGE 2 RUNS FOR ALL TEN, INCLUDING A1/E1. Stage 2 at causal 0.01 does exist
# for A1/E1 (in simulations_round_3_2Mb_g5t20) and could have been adopted, but
# it was written before causal_min_maf existed. The historical 0.01 is meant to
# be a no-op through the phenotype layer; rerunning settles that rather than
# assuming it, and stage 2 measured 2:17 / 369 MB, so uniformity is nearly free.
#
# STAGE 1 IS REUSED for all ten -- these are the same genetic simulations as the
# cmaf arms, which is what makes the three arms comparable. The tree sequences
# live in two places: A1/E1 in the x35 pilot root, A2-A5/E2-E5 in the
# cmaf_0.001 root where they were first generated. Both are searched.
#
# NOT dependent on the _fm_thins_candidates patch in rules/common.smk: at
# fm_min_maf=0.01 that helper returns true exactly as the old `> 0` test did, so
# this arm gets the same DAP-G resource tiers it always did.
#
# Run from a login node:
#   bash submit_2Mb_r3_cmaf01_control.sh
#   DRY=1 bash submit_2Mb_r3_cmaf01_control.sh      # print, submit nothing
#   REPS="A2 A3" bash submit_2Mb_r3_cmaf01_control.sh

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SCRATCH=/n/scratch/users/n/njc12/snakemake
PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake
cd "$REPO"

CELL="${CELL:-g5t20}"
CAUSAL_FLOOR="${CAUSAL_FLOOR:-0.01}"
FM_FLOOR="${FM_FLOOR:-0.01}"
GLM_FLOOR="${GLM_FLOOR:-0.01}"
LD_CTRL="${LD_CTRL:-0.75}"
REPS="${REPS:-A1 A2 A3 A4 A5 E1 E2 E3 E4 E5}"
DRY="${DRY:-}"
BEGIN="${BEGIN:-}"

OUT_SCRATCH="$SCRATCH/simulations_round_3_2Mb_${CELL}_cmaf_${CAUSAL_FLOOR}_fm_${FM_FLOOR}"
OUT_PUBLISH="$PUBLISH/simulations_round_3_2Mb_${CELL}_cmaf_${CAUSAL_FLOOR}_fm_${FM_FLOOR}"

# Where the stage-1 tree sequences actually are, most-specific first.
#
# DELIBERATELY NOT parameterised by CELL. Stage 1 is the GENETIC simulation and
# does not depend on the phenotype multipliers at all, so a replicate has exactly
# one tree sequence no matter which cell is being run -- it lives wherever it
# happened to be generated first. Interpolating $CELL here would make every new
# cell fail pre-flight on A2-A5/E2-E5, whose tree sequences only ever existed in
# the cmaf_0.001 root, and would invite a redundant stage-1 rerun per cell.
STAGE1_ROOTS=(
    "$SCRATCH/simulations_round_3_2Mb_g5t20_cmaf_0.001"   # A2-A5, E2-E5
    "$SCRATCH/simulations_round_3_2Mb"                    # A1, E1 (x35 pilot)
)

# Seed convention from the configs' own headers: seed = 10*letter_index +
# replicate, so A{N} -> 1{N}, E{N} -> 5{N}, J (the tenth letter) -> 10{N},
# M (the thirteenth) -> 13{N} and N (the fourteenth) -> 14{N}.
seed_of() {
    local n="${1:1}"
    case "${1:0:1}" in
        A) echo "1${n}" ;;
        E) echo "5${n}" ;;
        H) echo "8${n}" ;; I) echo "9${n}" ;;
        J) echo "10${n}" ;;
        M) echo "13${n}" ;; N) echo "14${n}" ;;
        *) echo "ERROR: no seed convention for '$1'" >&2; return 1 ;;
    esac
}

#
# M and N are the FINNISH FOUNDER PAIR, and unlike B/K/L they need real config
# files: `demographic_model: FinnishWang2014` selects a demes graph that is not
# in stdpopsim's catalog, and it forces `Q_scaling: 3` because stdpopsim's SLiM
# engine will not sample 9,000 individuals out of a FIN deme that holds 2,266 at
# Q=10 or 7,491 at Q=4. At Q=3 it holds 9,988 -- a margin of 988. M samples FIN, N samples NFE, and the two files are otherwise identical:
#   M - N  isolates the Finnish founder event
#   N - A  isolates the model swap (Wang NFE vs Tennessen EUR)
# Only the g5t20 cell exists so far; any other CELL fails the pre-flight
# `[[ -f "$REPO/$CFG" ]]` check, which is the intended way to find that out.
config_of() {
    case "${1:0:1}" in
        A) echo "config/human_2Mb_${CELL}_r3.yaml" ;;
        J) echo "config/human_afr_2Mb_${CELL}_r3.yaml" ;;
        M) echo "config/human_fin_2Mb_${CELL}_r3.yaml" ;;
        N) echo "config/human_nfe_2Mb_${CELL}_r3.yaml" ;;
        H) echo "config/human_neutral_2Mb_${CELL}_r3.yaml" ;;
        I) echo "config/cattle_neutral_2Mb_${CELL}_r3.yaml" ;;
        E) echo "config/cattle_baseline_from_midpoint_2Mb_${CELL}_r3.yaml" ;;
        *) return 1 ;;
    esac
}

# ---------------------------------------------------------------- pre-flight
# Resolve every reused stage-1 path before submitting anything. A stage-1 path
# that resolves to nothing does not fail loudly: Snakemake simply decides it has
# to BUILD the tree sequence, and a multi-hour genetic simulation starts instead
# of a reuse. The lookup is also seed-strict, so a directory holding some other
# replicate's .ts is silently no help -- hence checking per replicate.
declare -A STAGE1_SRC
for ID in $REPS; do
    CFG=$(config_of "$ID") || { echo "ERROR: unknown replicate '$ID'" >&2; exit 1; }
    [[ -f "$REPO/$CFG" ]] || { echo "ERROR: missing config: $REPO/$CFG" >&2; exit 1; }
    seed_of "$ID" > /dev/null || exit 1

    # Categories H, I and J build their own stage 1 rather than reusing one.
    # H and I are distinct pipelines (human_neutral / cattle_neutral, pure
    # coalescents -- seconds to minutes, not hours), so there is no x35
    # predecessor to point at. J runs the SAME pipeline as A but samples the AFR
    # branch, so its tree sequence is hts_AFR_{seed}.ts and no EUR stage1 dir can
    # satisfy it -- pointing at one would find nothing, and find nothing
    # silently. An empty STAGE1_SRC tells the controller to omit
    # stage1_search_dirs entirely.
    if [[ "${ID:0:1}" =~ ^(H|I|J|M|N)$ ]]; then
        STAGE1_SRC["$ID"]=""
    else
        S1=""
        for ROOT in "${STAGE1_ROOTS[@]}"; do
            if compgen -G "$ROOT/${ID}/stage1/*.ts" > /dev/null; then S1="$ROOT/${ID}/stage1"; break; fi
        done
        [[ -n "$S1" ]] || {
            echo "ERROR: no stage-1 tree sequence for ${ID} in any of:" >&2
            for ROOT in "${STAGE1_ROOTS[@]}"; do echo "       $ROOT/${ID}/stage1" >&2; done
            echo "Aborting -- no jobs launched." >&2; exit 1; }
        STAGE1_SRC["$ID"]="$S1"
    fi

    # Refuse to relaunch into a workdir that already holds a finished stage 2:
    # the trait draw would be regenerated over results that may already be
    # summarized.
    if compgen -G "$OUT_SCRATCH/${ID}/stage2/*/*/.stage2.done" > /dev/null; then
        echo "ERROR: $OUT_SCRATCH/${ID} already has a completed stage 2." >&2
        echo "       Remove it deliberately or pick another root. Nothing launched." >&2
        exit 1
    fi
done

echo "Pre-flight OK: $(echo $REPS | wc -w) replicates."
echo "  causal_min_maf=$CAUSAL_FLOOR  fm_min_maf=$FM_FLOOR  min_maf=$GLM_FLOOR  ld_ctrl=$LD_CTRL"
echo "  target root: $OUT_SCRATCH"
echo "  stage 1 REUSED; stage 2 RUNS for every replicate (including A1/E1)."
echo

# ------------------------------------------------------------------- submit
N=0
for ID in $REPS; do
    CFG=$(config_of "$ID")
    SEED=$(seed_of "$ID")
    WD="$OUT_SCRATCH/${ID}"
    PD="$OUT_PUBLISH/${ID}"
    # EXTRA_CONFIG is expanded unquoted by the controller so each space-separated
    # token becomes its own --config arg. No token may contain a space. All three
    # floors are passed explicitly even though 0.01 is what the configs and the
    # causal default already say, so the params file records what was INTENDED
    # rather than whatever the config happened to hold on the day.
    EXTRA="causal_min_maf=${CAUSAL_FLOOR} fm_min_maf=${FM_FLOOR} min_maf=${GLM_FLOOR} ld_ctrl=${LD_CTRL}"
    if [[ -n "$DRY" ]]; then
        printf "  DRY %-3s seed %-3s %s\n        STAGE1_SRC=%s\n        EXTRA_CONFIG=%s\n" \
            "$ID" "$SEED" "$CFG" "${STAGE1_SRC[$ID]}" "$EXTRA"
        N=$((N + 1)); continue
    fi
    mkdir -p "$WD"
    JID=$(sbatch --parsable --job-name="${ID}_${CELL}_cmaf${CAUSAL_FLOOR}_fm${FM_FLOOR}" \
        ${BEGIN:+--begin="$BEGIN"} \
        --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
        --export=ALL,CONFIGFILE="${CFG}",SEED="${SEED}",WD="${WD}",PD="${PD}",STAGE1_SRC="${STAGE1_SRC[$ID]}",EXTRA_CONFIG="${EXTRA}" \
        controller_2Mb_x10.sbatch)
    printf "  %-3s seed %-3s controller %s\n" "$ID" "$SEED" "$JID"
    N=$((N + 1))
done
echo
echo "${DRY:+Would submit }${N} controllers${DRY:+ (DRY=1, nothing launched)}."
[[ -n "$DRY" ]] && exit 0
echo "Watch with: squeue -u \$USER -o '%.10i %.24j %.9P %.2t %.10M %R'"
echo
echo "Check once stage 2 lands -- the causative pool must be SMALLER than the"
echo "cmaf=0.001 arm's, or the floor did not take effect:"
echo "  grep -E 'causative|neutral|flank' $OUT_SCRATCH/A1/logs/stage2_split_pheno.log"
echo "and the fine-mapping floor actually applied should read $FM_FLOOR:"
echo "  cat $OUT_SCRATCH/A1/stage3/*/hgwas/geno.sbams.gz.fmmaf"
