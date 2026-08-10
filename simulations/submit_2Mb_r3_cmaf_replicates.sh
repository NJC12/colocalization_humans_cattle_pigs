#!/bin/bash
# Additional REPLICATES for the causal-MAF-floor arm of the g5t20 cell.
#
# The cmaf=0.001 experiment currently rests on one stage-1 simulation per arm
# (A1, E1), so every rate in it is a count out of 50 loci with a binomial SE of
# roughly +/-3 counts -- too coarse to read the human/cattle comparison off. This
# adds A2-A5 and E2-E5 at the same settings so the summary table pools 5
# replicates per arm (n_gwas_traits 250 instead of 50).
#
# Settings, all of which must match the existing A1/E1 cmaf run exactly or the
# summary table will split them into separate rows instead of pooling:
#   causal_min_maf = 0.001    causative variants drawn from MAF >= 0.001
#   fm_min_maf     = 0.01     but only MAF >= 0.01 variants are fine-mapped
#   min_maf        = 0.01     and only MAF >= 0.01 variants enter the plink GLM
#   ld_ctrl        = 0.75     dap-g signal-cluster r2 threshold
#   gwas_scaling   = 5        from config/*_2Mb_g5t20_r3.yaml
#   gtex_scaling   = 20       from config/*_2Mb_g5t20_r3.yaml
# The last four already sit in the configs; they are passed explicitly anyway so
# the run is pinned to what was intended rather than to whatever the config says
# on the day, and so the params file records them unambiguously.
#
# THE IMPORTANT DIFFERENCE FROM ITS SIBLING SCRIPTS: stage 1 RUNS HERE.
#
# submit_2Mb_r3_{multipliers,finemap_variants,causal_maf}.sh all reuse A1/E1's
# stage-1 tree sequence via stage1_search_dirs, because they vary only downstream
# settings. These are NEW replicates -- new stage1_seed -- so there is nothing to
# reuse and the genetic simulation must actually execute. Hence:
#   - stage1_search_dirs is NOT set (and the lookup is seed-strict anyway, so
#     pointing it at A1 would find nothing and silently fall through);
#   - the controller is controller_2Mb.sbatch, which does not require STAGE1_SRC,
#     rather than controller_2Mb_x10.sbatch, which does.
#
# Cost is modest at L=2 Mb: stage1_human measured 1:05 / 232 MB, and the cattle
# arm only reruns epochs 8-12 from the SHARED deep-history checkpoint
# (cattle_baseline_seed 20250303, measured 4:29 for E1) -- the 30,000-generation
# burn-in is not repeated. Both stage-1 rules declare `short` themselves at this
# L, so nothing lands on priority or long.
#
# CATEGORIES. Originally A/E only; B, F and G were added when the category axis
# was extended to this arm. Seed bands come from the configs' own headers:
#   A{N} -> 1{N}   human, directional-negative selection on trait variants
#   B{N} -> 2{N}   human, NEUTRAL trait variants -- A's config plus the single
#                  override neutral_trait_vars=True. B has no config of its own.
#   E{N} -> 5{N}   cattle baseline from midpoint
#   F{N} -> 6{N}   cattle + positive selection, WITH continued bottlenecking
#   G{N} -> 7{N}   cattle + positive selection, WITHOUT it
# F and G have no round-3 tree sequences at all, so stage 1 runs for them here
# exactly as it does for a new A/E replicate. Both load the same shared ep7
# handoff that E does, which the pre-flight now checks for all three.
#
# Run from a login node:
#   bash submit_2Mb_r3_cmaf_replicates.sh
#   REPS="A2 A3" bash submit_2Mb_r3_cmaf_replicates.sh     # subset
#   DRY=1 bash submit_2Mb_r3_cmaf_replicates.sh            # print, do not submit
#   REPS="B1 F1 G1" CTRL_TIME=2-00:00:00 bash submit_2Mb_r3_cmaf_replicates.sh

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SCRATCH=/n/scratch/users/n/njc12/snakemake
PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake
cd "$REPO"

CELL="${CELL:-g5t20}"
FLOOR="${FLOOR:-0.001}"
FM_FLOOR="${FM_FLOOR:-0.01}"
GLM_FLOOR="${GLM_FLOOR:-0.01}"
LD_CTRL="${LD_CTRL:-0.75}"
REPS="${REPS:-A2 A3 A4 A5 E2 E3 E4 E5}"
DRY="${DRY:-}"
# Controller walltime. 1 day held for every A/E replicate, but G climbs the
# DAP-G retry ladder (32 -> 64 -> 96 GB; see _dapg_mem_mb_base) and has never run
# under round-3 code, so give B/F/G 2 days. `medium` allows up to 5.
CTRL_TIME="${CTRL_TIME:-1-00:00:00}"

# Same root as the existing A1/E1 cmaf run, so summarize_coloc.py pools them into
# one row per arm. A different root would still work but would report 1-replicate
# and 4-replicate rows side by side.
OUT_SCRATCH="$SCRATCH/simulations_round_3_2Mb_${CELL}_cmaf_${FLOOR}"
OUT_PUBLISH="$PUBLISH/simulations_round_3_2Mb_${CELL}_cmaf_${FLOOR}"

# Seed convention is fixed by the configs' own header comments:
#   A{N} -> 1{N}   (A1=11 ... A10=20)
#   B{N} -> 2{N}   (B1=21 ... B10=30)
#   E{N} -> 5{N}   (E1=51 ... E10=60)
#   F{N} -> 6{N}   (F1=61 ... F10=70)
#   G{N} -> 7{N}   (G1=71 ... G10=80)
#   H{N} -> 8{N}   (H1=81 ... H10=90)
#   I{N} -> 9{N}   (I1=91 ... I9=99; I is the last letter the
#                   tens-digit rule has room for)
# Deriving it rather than hard-coding a table keeps it from drifting.
seed_of() {
    local id="$1" letter="${1:0:1}" n="${1:1}"
    case "$letter" in
        A) echo "1${n}" ;;
        B) echo "2${n}" ;;
        E) echo "5${n}" ;;
        F) echo "6${n}" ;;
        G) echo "7${n}" ;;
        H) echo "8${n}" ;; I) echo "9${n}" ;;
        *) echo "ERROR: no seed convention for '$id'" >&2; return 1 ;;
    esac
}

# B shares A's config: the two categories differ only in neutral_trait_vars,
# which category_extra() supplies as a --config override.
#
# H and I do NOT share another category's config, unlike B. They are different
# stage-1 pipelines (human_neutral / cattle_neutral: pure coalescents, no SLiM),
# so each needs its own file rather than an override -- and having one means the
# seed convention below can stay derived from the letter. There is deliberately
# no config/human_2Mb_*_B.yaml -- one file, one genetic model.
#
# I also needs no cattle_baseline_search_dirs: unlike E/F/G it simulates all
# twelve epochs itself instead of resuming from the shared ep7 checkpoint, which
# is why the pre-flight handoff guard below stays scoped to ^(E|F|G)$.
config_of() {
    case "${1:0:1}" in
        A|B) echo "config/human_2Mb_${CELL}_r3.yaml" ;;
        H)   echo "config/human_neutral_2Mb_${CELL}_r3.yaml" ;;
        I)   echo "config/cattle_neutral_2Mb_${CELL}_r3.yaml" ;;
        E)   echo "config/cattle_baseline_from_midpoint_2Mb_${CELL}_r3.yaml" ;;
        F)   echo "config/cattle_sel_bottlenecked_2Mb_${CELL}_r3.yaml" ;;
        G)   echo "config/cattle_sel_not_bottlenecked_2Mb_${CELL}_r3.yaml" ;;
        *)   return 1 ;;
    esac
}

# Per-category --config tokens appended to EXTRA. NO TOKEN MAY CONTAIN A SPACE:
# EXTRA_CONFIG is expanded unquoted by the controller so each space-separated
# token becomes its own --config arg.
category_extra() {
    case "${1:0:1}" in
        B) echo "neutral_trait_vars=True" ;;
        *) echo "" ;;
    esac
}

# ---------------------------------------------------------------- pre-flight
for ID in $REPS; do
    CFG=$(config_of "$ID") || { echo "ERROR: unknown replicate '$ID'" >&2; exit 1; }
    [[ -f "$REPO/$CFG" ]] || { echo "ERROR: missing config: $REPO/$CFG" >&2; exit 1; }
    seed_of "$ID" > /dev/null || exit 1

    # Cattle replicates load the SHARED deep-history handoff. Without this guard a
    # missing checkpoint surfaces only when the stage-1 job runs and fails.
    # F and G (pipeline cattle_sel) read the same ep7 checkpoint E does -- they
    # differ only in applying positive selection on top of it -- so all three are
    # checked, not just E.
    if [[ "${ID:0:1}" =~ ^(E|F|G)$ ]]; then
        CB_DIR=$(awk '/^cattle_baseline_search_dirs:/{f=1;next} f&&/^ *- /{gsub(/^ *- */,"");print;exit}' "$REPO/$CFG")
        CB_SEED=$(awk '/^cattle_baseline_seed:/{print $2;exit}' "$REPO/$CFG")
        L_VAL=$(awk '/^L:/{print $2;exit}' "$REPO/$CFG")
        Q_VAL=$(awk '/^Q_scaling:/{print $2;exit}' "$REPO/$CFG")
        CKPT="${CB_DIR}/farm_selection_Q_${Q_VAL}.L_${L_VAL}.seed_${CB_SEED}.ep7.ts"
        [[ -f "$CKPT" ]] || {
            echo "ERROR: missing shared cattle handoff checkpoint for ${ID}:" >&2
            echo "       $CKPT" >&2; exit 1; }
    fi

    # Refuse to launch into a workdir that already holds a finished stage 2: the
    # trait draw would be regenerated over the top of results already summarized.
    if compgen -G "$OUT_SCRATCH/${ID}/stage2/*/*/.stage2.done" > /dev/null; then
        echo "ERROR: $OUT_SCRATCH/${ID} already has a completed stage 2." >&2
        echo "       Remove it deliberately or pick another root. Nothing launched." >&2
        exit 1
    fi
done
echo "Pre-flight OK: $(echo $REPS | wc -w) replicates."
echo "  causal_min_maf=$FLOOR  fm_min_maf=$FM_FLOOR  min_maf=$GLM_FLOOR  ld_ctrl=$LD_CTRL"
echo "  root: $OUT_SCRATCH"
echo "  replicates: $REPS"
for ID in $REPS; do
    CAT_EXTRA=$(category_extra "$ID")
    [[ -n "$CAT_EXTRA" ]] && echo "    $ID adds: $CAT_EXTRA"
done
echo "  STAGE 1 WILL RUN for every replicate (new seeds, nothing to reuse)."
echo "  controller walltime: $CTRL_TIME"
echo

# ------------------------------------------------------------------- submit
# controller_2Mb.sbatch declares long/30-day because it predates the resource
# rightsizing; these DAGs finish in a couple of hours, so the CLI flags below
# override it down to medium/1-day. sbatch options beat in-script #SBATCH.
N=0
for ID in $REPS; do
    CFG=$(config_of "$ID")
    SEED=$(seed_of "$ID")
    WD="$OUT_SCRATCH/${ID}"
    PD="$OUT_PUBLISH/${ID}"
    EXTRA="causal_min_maf=${FLOOR} fm_min_maf=${FM_FLOOR} min_maf=${GLM_FLOOR} ld_ctrl=${LD_CTRL}"
    CAT_EXTRA=$(category_extra "$ID")
    [[ -n "$CAT_EXTRA" ]] && EXTRA="${EXTRA} ${CAT_EXTRA}"
    if [[ -n "$DRY" ]]; then
        printf "  DRY %-3s seed %-3s %s\n        EXTRA_CONFIG=%s\n" "$ID" "$SEED" "$CFG" "$EXTRA"
        N=$((N + 1)); continue
    fi
    mkdir -p "$WD"
    JID=$(sbatch --parsable --job-name="${ID}_${CELL}_cmaf${FLOOR}" \
        --partition=medium --time="${CTRL_TIME}" \
        --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
        --export=ALL,CONFIGFILE="${CFG}",SEED="${SEED}",WD="${WD}",PD="${PD}",EXTRA_CONFIG="${EXTRA}" \
        controller_2Mb.sbatch)
    printf "  %-3s seed %-3s controller %s\n" "$ID" "$SEED" "$JID"
    N=$((N + 1))
done
echo
echo "${DRY:+Would submit }${N} controllers${DRY:+ (DRY=1, nothing launched)}."
[[ -n "$DRY" ]] && exit 0
echo "Watch with: squeue -u \$USER -o '%.10i %.24j %.9P %.2t %.10M %R'"
echo
echo "Check once stage 2 lands, per replicate -- the causative pool must be larger"
echo "than the cmaf=0.01 cell's, or the floor did not take effect:"
echo "  grep -E 'causative|neutral|flank' $OUT_SCRATCH/A2/logs/stage2_split_pheno.log"
