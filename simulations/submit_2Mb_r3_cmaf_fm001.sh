#!/bin/bash
# Lower the FINE-MAPPING and GLM floors to 0.001 on the existing causal-0.001 arm.
#
# The five-replicate cmaf arm answered its question with an unambiguous no: at
# causal_min_maf=0.001 with fm_min_maf=0.01, colocalization power at loci whose
# causal variant fell below the fine-mapping floor was 1 detection out of 341.
# The tagging route does not carry fastEnloc. 43% of human loci sat in that
# untestable group (vs 11% of cattle), which is the whole of the human/cattle gap
# in that arm.
#
# This run removes the floor that created the untestable group:
#
#   causal_min_maf = 0.001   causative variants drawn from MAF >= 0.001  (UNCHANGED)
#   fm_min_maf     = 0.001   and now fine-mapped down to the same floor  (was 0.01)
#   min_maf        = 0.001   and entered into the GWAS at the same floor (was 0.01)
#
# so nearly every causal variant becomes a fine-mapping candidate and the power
# question is asked without the composition effect in the way.
#
# WHAT 0.001 ACTUALLY MEANS PER PANEL. The floor is a frequency, so what it
# removes depends on panel size:
#     GWAS   n=8000  ->  MAC >= 16
#     GTEx   n=1000  ->  MAC >= 2
#     GTEx   n=500   ->  MAC >= 1, i.e. no filtering at all
# The GTEx-500 arm is therefore effectively unfiltered. That is intended -- it is
# the same variant set the existing fm-unfiltered cells used -- but it is worth
# knowing before reading the 500-panel numbers as "filtered at 0.001".
#
# PAIRED WITH THE fm=0.01 ARM. Stage 2 is ADOPTED wholesale from
# simulations_round_3_2Mb_g5t20_cmaf_0.001 rather than redrawn. causal_min_maf,
# both scalings, both seeds and the trait counts are all unchanged, so the
# phenotypes are identical by construction -- and only fm_min_maf and min_maf
# differ between the two arms. Every locus in this run is the SAME locus as in
# the fm=0.01 run, which is the pairing the first cmaf-vs-control comparison
# lacked. Redrawing stage 2 would confound the fine-mapping floor with the draw.
#
# Snakefile enforces the reuse is legitimate: stage2_params.txt in the adopted
# directory records causal_min_maf, and the provenance guard refuses a directory
# whose stage-2 keys disagree. min_maf and fm_min_maf are ANALYSIS keys, not
# stage-2 keys, so changing them does not (and must not) invalidate the adoption.
#
# DEPENDS ON A common.smk PATCH. rules/common.smk used to pick DAP-G's resource
# tier on `fm_min_maf > 0`, which was calibrated at 0.01 where human drops
# 2847 -> 327 candidates per window. At 0.001 human keeps essentially all 2847,
# so that test handed it a 20 min / 4 GB allocation for work measured at 23:50
# and 2580 MB peak -- a walltime kill on the upper tail of every human GWAS
# locus. The test is now _fm_thins_candidates(), which compares the floor against
# _FM_THINNING_FLOOR = 0.005. Verify before running:
#     grep -n "_FM_THINNING_FLOOR\|_fm_thins_candidates" rules/common.smk
# If that returns nothing, the O2 checkout predates the patch: rsync first, or
# every human hgwas locus burns two failed attempts before the retry ladder
# (20 -> 40 -> 60 min) gets it through.
#
# CATEGORIES. Originally A/E; B, F and G were added when the category axis was
# extended to this arm. Both the stage-1 AND the stage-2 source for B/F/G is the
# cmaf_0.001 root, so submit_2Mb_r3_cmaf_replicates.sh must have COMPLETED for
# the same replicate before this script can run it -- stage-2 adoption needs
# .stage2.done, which is stricter than the psamp arm's stage-1-only dependency.
# The per-replicate pre-flight below is what catches a missing one.
#
# Run from a login node:
#   bash submit_2Mb_r3_cmaf_fm001.sh
#   DRY=1 bash submit_2Mb_r3_cmaf_fm001.sh          # print, submit nothing
#   REPS="A1 E1" bash submit_2Mb_r3_cmaf_fm001.sh   # subset
#   REPS="B1 F1 G1" bash submit_2Mb_r3_cmaf_fm001.sh  # the B/F/G pilot

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SCRATCH=/n/scratch/users/n/njc12/snakemake
PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake
cd "$REPO"

CELL="${CELL:-g5t20}"
CAUSAL_FLOOR="${CAUSAL_FLOOR:-0.001}"   # must match the source arm, or stage 2 will not adopt
FM_FLOOR="${FM_FLOOR:-0.001}"
GLM_FLOOR="${GLM_FLOOR:-0.001}"
LD_CTRL="${LD_CTRL:-0.75}"
REPS="${REPS:-A1 A2 A3 A4 A5 E1 E2 E3 E4 E5}"
DRY="${DRY:-}"
BEGIN="${BEGIN:-}"

# Source of stage 1 and stage 2: the existing fm=0.01 arm.
SRC_ROOT="$SCRATCH/simulations_round_3_2Mb_${CELL}_cmaf_${CAUSAL_FLOOR}"
# A1 and E1 predate the replicate expansion and reused the x35 pilot's stage-1
# tree sequence, so their own root has no stage1/ directory. A2-A5 and E2-E5 ran
# stage 1 themselves and have one. Both are searched, most-specific first.
STAGE1_FALLBACK="$SCRATCH/simulations_round_3_2Mb"

OUT_SCRATCH="$SCRATCH/simulations_round_3_2Mb_${CELL}_cmaf_${CAUSAL_FLOOR}_fm_${FM_FLOOR}"
OUT_PUBLISH="$PUBLISH/simulations_round_3_2Mb_${CELL}_cmaf_${CAUSAL_FLOOR}_fm_${FM_FLOOR}"

# Seed convention from the configs' own headers:
#   A{N} -> 1{N}   human, directional-negative selection on trait variants
#   B{N} -> 2{N}   human, NEUTRAL trait variants (A's config + one override)
#   E{N} -> 5{N}   cattle baseline from midpoint
#   F{N} -> 6{N}   cattle + positive selection, WITH continued bottlenecking
#   G{N} -> 7{N}   cattle + positive selection, WITHOUT it
seed_of() {
    local n="${1:1}"
    case "${1:0:1}" in
        A) echo "1${n}" ;;
        B) echo "2${n}" ;;
        E) echo "5${n}" ;;
        F) echo "6${n}" ;;
        G) echo "7${n}" ;;
        H) echo "8${n}" ;; I) echo "9${n}" ;;
        *) echo "ERROR: no seed convention for '$1'" >&2; return 1 ;;
    esac
}

# B shares A's config: the two differ only in neutral_trait_vars, which
# category_extra() supplies as a --config override.
#
# H and I do NOT share another category's config, unlike B. They are different
# stage-1 pipelines (human_neutral / cattle_neutral: pure coalescents, no SLiM),
# so each needs its own file rather than an override -- and having one means the
# seed convention below can stay derived from the letter. One file, one genetic
# model. I is the cattle counterpart of H: same drawn-DFE effect sizes, cattle
# demography, and unlike E/F/G it simulates all twelve epochs itself rather than
# resuming from the shared ep7 checkpoint.
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
# Resolve and validate EVERY reused path before submitting anything. A stage-1
# path that resolves to nothing does not fail loudly -- Snakemake would simply
# decide it has to build the tree sequence, and a real multi-hour genetic
# simulation starts instead of a reuse. A stage-2 path that misses is milder
# (stage 2 reruns, correctly) but silently breaks the pairing this arm exists
# for, so both are hard errors.
if ! grep -q "_fm_thins_candidates" "$REPO/rules/common.smk" 2>/dev/null; then
    echo "ERROR: $REPO/rules/common.smk predates the DAP-G resource patch." >&2
    echo "       At fm_min_maf=$FM_FLOOR human would get the MAF>=0.01 tier" >&2
    echo "       (20 min / 4 GB) for an essentially unfiltered candidate set." >&2
    echo "       rsync the repo to O2 first. Nothing launched." >&2
    exit 1
fi

declare -A STAGE1_SRC STAGE2_SRC
for ID in $REPS; do
    CFG=$(config_of "$ID") || { echo "ERROR: unknown replicate '$ID'" >&2; exit 1; }
    [[ -f "$REPO/$CFG" ]] || { echo "ERROR: missing config: $REPO/$CFG" >&2; exit 1; }
    seed_of "$ID" > /dev/null || exit 1

    # Categories H and I build their own stage 1 rather than reusing one. They are
    # distinct pipelines (human_neutral / cattle_neutral, pure coalescents --
    # seconds to minutes, not hours), so there
    # is no predecessor to point at and nothing to save by looking. An empty
    # STAGE1_SRC tells the controller to omit stage1_search_dirs entirely.
    if [[ "${ID:0:1}" =~ ^(H|I)$ ]]; then
        STAGE1_SRC["$ID"]=""
    else
        S1=""
        for CAND in "$SRC_ROOT/${ID}/stage1" "$STAGE1_FALLBACK/${ID}/stage1"; do
            if compgen -G "$CAND/*.ts" > /dev/null; then S1="$CAND"; break; fi
        done
        [[ -n "$S1" ]] || {
            echo "ERROR: no stage-1 tree sequence for ${ID} in either" >&2
            echo "       $SRC_ROOT/${ID}/stage1" >&2
            echo "       $STAGE1_FALLBACK/${ID}/stage1" >&2
            echo "Aborting -- no jobs launched." >&2; exit 1; }
        STAGE1_SRC["$ID"]="$S1"
    fi

    # stage2_search_dirs wants the directory CONTAINING gwas_*_gtex_*_maf_*/,
    # i.e. <root>/<id>/stage2/<stage1 name>.cmaf_<floor>/.
    S2=$(ls -d "$SRC_ROOT/${ID}/stage2"/*/ 2>/dev/null | head -1 || true)
    [[ -n "$S2" && -d "$S2" ]] || {
        echo "ERROR: no stage-2 output for ${ID} under $SRC_ROOT/${ID}/stage2" >&2
        echo "Run the fm=0.01 cmaf arm first. Nothing launched." >&2; exit 1; }
    S2="${S2%/}"
    compgen -G "$S2/*/.stage2.done" > /dev/null || {
        echo "ERROR: stage 2 for ${ID} in $S2 is not marked complete." >&2
        echo "Adoption is all-or-nothing; a partial dir would silently rerun stage 2." >&2
        exit 1; }
    # The adopted directory must carry the causal floor this run claims, or the
    # arms are not paired. Snakefile's guard checks this too; checking here means
    # finding out now rather than after ten controllers have queued.
    if compgen -G "$S2/*/stage2_params.txt" > /dev/null; then
        SEEN=$(awk -F': *' '/^causal_min_maf:/{print $2; exit}' "$S2"/*/stage2_params.txt)
        awk -v a="$SEEN" -v b="$CAUSAL_FLOOR" 'BEGIN{exit !(a+0==b+0)}' </dev/null || {
            echo "ERROR: ${ID} stage 2 records causal_min_maf=$SEEN, not $CAUSAL_FLOOR." >&2
            exit 1; }
    fi
    STAGE2_SRC["$ID"]="$S2"

    # Refuse to write over a finished stage 4 in the TARGET root.
    if compgen -G "$OUT_SCRATCH/${ID}/stage4/*/*.enloc.sig.out" > /dev/null; then
        echo "ERROR: $OUT_SCRATCH/${ID} already holds stage-4 output." >&2
        echo "       Remove it deliberately or pick another root. Nothing launched." >&2
        exit 1
    fi
done

echo "Pre-flight OK: $(echo $REPS | wc -w) replicates."
echo "  causal_min_maf=$CAUSAL_FLOOR (unchanged)  fm_min_maf=$FM_FLOOR  min_maf=$GLM_FLOOR  ld_ctrl=$LD_CTRL"
echo "  stage 1 + stage 2 REUSED from $SRC_ROOT"
echo "  target root: $OUT_SCRATCH"
echo "  Only stages 3-5 will run; the trait draw is identical to the fm=0.01 arm."
echo

# ------------------------------------------------------------------- submit
N=0
for ID in $REPS; do
    CFG=$(config_of "$ID")
    SEED=$(seed_of "$ID")
    WD="$OUT_SCRATCH/${ID}"
    PD="$OUT_PUBLISH/${ID}"
    # EXTRA_CONFIG is expanded unquoted by the controller so each space-separated
    # token becomes its own --config arg. No token may contain a space -- the
    # list literal below deliberately has none.
    EXTRA="causal_min_maf=${CAUSAL_FLOOR} fm_min_maf=${FM_FLOOR} min_maf=${GLM_FLOOR}"
    EXTRA="${EXTRA} ld_ctrl=${LD_CTRL} stage2_search_dirs=['${STAGE2_SRC[$ID]}']"
    # Category tokens are NOT optional here. neutral_trait_vars is in
    # params_record.STAGE2_KEYS, so a B replicate that omitted it would be
    # refused by the stage-2 provenance guard rather than silently adopting the
    # source arm's stage 2 under the wrong label.
    CAT_EXTRA=$(category_extra "$ID")
    [[ -n "$CAT_EXTRA" ]] && EXTRA="${EXTRA} ${CAT_EXTRA}"
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
echo "First check once stage 3 starts, per replicate -- stage 2 must have been"
echo "ADOPTED, not rerun. The stage-2 dir should be a symlink into the source arm:"
echo "  ls -l $OUT_SCRATCH/A1/stage2/"
echo "and the fine-mapping floor actually applied should read $FM_FLOOR:"
echo "  cat $OUT_SCRATCH/A1/stage3/*/hgwas/geno.sbams.gz.fmmaf"
