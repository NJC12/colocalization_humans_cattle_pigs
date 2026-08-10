#!/bin/bash
# Power-weighted causal sampling at the x20 cell, fine-mapping and GLM floors 0.001.
#
#   causal_sampling  power    causal loci drawn with inclusion probability
#                             proportional to the power to detect them at n=8000
#   sampling_gwas_n  8000     the GWAS size that power is computed for
#   sampling_sig_p   5e-8     what "detected" means (= GWAS_SIG_P in figure 2)
#   gwas_scaling     20       from config/*_2Mb_x20_r3.yaml
#   gtex_scaling     20       from config/*_2Mb_x20_r3.yaml
#   fm_min_maf       0.001    fine-mapping floor
#   min_maf          0.001    plink2 --glm floor
#   ld_ctrl          0.75     dap-g signal-cluster r2 threshold
#   causal_min_maf   0        no floor on the causal pool. This USED to read 0.01
#                             with the note "not a floor on the causal pool under
#                             power sampling" -- the floor was scheme-dependent and
#                             power ignored it. It no longer is: a floor applies
#                             wherever it is set, in both schemes. 0 is therefore
#                             what reproduces the runs already on disk, which is
#                             also why their directories were relabelled maf_0.
#                             Set CAUSAL_FLOOR=0.01 to get an actually-floored arm.
#
# A1-A5 + E1-E5, so the arm is 500 GWAS traits rather than 50.
#
# WHY. Uniform sampling draws causal loci from the whole MAF >= causal_min_maf
# pool, most of which no GWAS of realistic size could find -- 52-64% of the human
# pool sits below MAF 0.01. Power weighting picks loci the way a real study does,
# so the colocalization rate is measured on a locus set an actual analysis would
# assemble. The fm/GLM floors are dropped to 0.001 for the same reason the
# cmaf_fm001 arm dropped them: at 0.01 the untestable group (causal variant below
# the fine-mapping floor) colocalized 1 time in 341, and it was 43% of human loci
# against 11% of cattle -- the whole of the human/cattle gap in that arm.
#
# READ THIS BEFORE TRUSTING ANY NUMBER FROM THIS RUN
#
# Power sampling lets the GWAS and shared-GTEx causal sets DIVERGE. Uniform
# sampling intersects the eligible pool with the GTEx panel, so the two sets are
# the same 50 positions by construction; power sampling does not, so a drawn locus
# the GTEx panel does not carry has no partner, and the GTEx central set is topped
# up uniformly to 50. Both scorers -- helpers/summarize_coloc.py and
# figures_and_tables/figure2_revision2.ipynb -- currently define a true positive as
# trait-NAME equality and keep every GWAS trait in the denominator. A partnerless
# locus therefore reads as a colocalization failure, and as a FALSE POSITIVE if it
# colocalizes with a neighbour. That fix is deliberately not part of this arm.
# Stage 2 prints the partner count and writes the explicit pairing table; read it
# before comparing power or FDR against any uniform arm.
#
# CATEGORY B HAS NEVER COMPLETED UNDER POWER SAMPLING. RUN B1 ALONE FIRST.
#
# B sets neutral_trait_vars=True, which redistributes each donor's effect onto a
# random NEUTRAL recipient and names the trait for the recipient. Two things about
# combining it with power sampling, one settled and one still open:
#
# SETTLED (2026-08-04). The combination used to crash outright, in stage 2, with
# `KeyError: <donor position>` from combine_phenos_to_df: the donor's selection
# coefficient was resolved against whichever panel was being phenotyped, and for
# the shared central loci that is the GTEx frame, which under power sampling does
# not carry ~76% of donors. The donor's selco now travels on its own row.
# test_causal_selection.py covers the loop (the "phenotype loop" section); before
# the fix every case there exercising the power draw passed redistribution=None.
#
# OPEN. Under UNIFORM sampling the recipient pool is intersected with the GTEx
# panel, so every GWAS trait is guaranteed a GTEx partner. Under POWER sampling
# that intersection is dropped for EVERY category
# (create_gwas_files_and_phenotypes.py, `topup_gtex`), so a
# GWAS trait can sit where the GTEx panel carries nothing. Both scorers count such
# a locus as a colocalization FAILURE and, if it hits a neighbour, as a FALSE
# POSITIVE. Measured in this arm at x30: A1 48/50, A2 47/50, and B1 48/50 before
# it died -- so B is NOT worse off than A, and power weighting does mitigate it as
# predicted, by favouring high-MAF recipients the GTEx panel also carries.
#
# Do not "fix" this by intersecting B's recipient pool with the GTEx panel: that
# would give B a partner rate A does not have, and the A-B difference would then
# be partly a partner-rate artefact rather than genetics.
#
# So: run B1, confirm it now COMPLETES and that its partner count is still 48/50
# (the draw is deterministic in stage2_seed, so a different number means the fix
# moved something it should not have), and only then launch B2-B5.
#
# STAGE 2 CAN REFUSE TO RUN. If fewer than 2 x n_central_traits (= 100) pool
# variants reach power 0.05, the draw would be a handful of certainties plus an
# arbitrary tail, and the script exits with the counts instead of producing it.
# At x20 that is not expected -- the multiplier is large -- but a refusal is a
# real answer about the cell, not a bug.
#
# THE CAUSAL SET IS CELL-SPECIFIC. beta = sqrt(|selco|) * gwas_scaling, so the
# weights depend on the multiplier and this arm is NOT paired variant-for-variant
# with the g5t20 / g10t20 / g5t30 cells the way the uniform grid is.
#
# COMPARISON ARM. Run submit_2Mb_r3_cmaf01_control.sh with CELL=x20 FM_FLOOR=0.001
# GLM_FLOOR=0.001 for the uniform-sampling counterpart: identical genetics,
# identical multipliers, identical analysis floors, differing only in which 50
# central variants were made causative.
#
# STAGE 1 IS REUSED for every replicate. The tree sequences live in two places:
# A1/E1 in the x35 pilot root, everything else (A2-A5, E2-E5, and all of B/F/G)
# in the cmaf_0.001 root where they were first generated. Both are searched.
# Stage 2 runs fresh for every replicate -- it must, because the causal draw is
# what this arm changes.
#
# B/F/G THEREFORE REQUIRE THAT submit_2Mb_r3_cmaf_replicates.sh HAS ALREADY RUN
# for the same replicate: that is the arm where their stage 1 is generated. The
# per-replicate pre-flight below is what catches a missing one.
#
# DEPENDS ON THE causal_sampling PATCH. If the O2 checkout predates it the run
# would silently fall back to a uniform draw and produce a mislabeled arm, so
# pre-flight refuses. rsync the repo first.
#
# OTHER MULTIPLIER CELLS. The x20 in the filename is this script's DEFAULT, not a
# restriction: CELL= selects any cell that has a config pair, and each lands in
# its own root (..._<CELL>_psamp_<n>_fm_<floor>), so the cells cannot collide.
# x10 and x30 were run this way. Note the pool-adequacy guard bites hardest at low
# multipliers -- at x10 the human pool clears the 100-variant floor with about 13%
# to spare, so a lower cell than x10 should be expected to refuse.
#
# Run from a login node:
#   bash submit_2Mb_r3_x20_psamp_fm001.sh
#   DRY=1 bash submit_2Mb_r3_x20_psamp_fm001.sh          # print, submit nothing
#   REPS="A1 E1" bash submit_2Mb_r3_x20_psamp_fm001.sh   # pilot pair first
#   CELL=x10 bash submit_2Mb_r3_x20_psamp_fm001.sh       # the 10x cell
#   CELL=x30 bash submit_2Mb_r3_x20_psamp_fm001.sh       # the 30x cell
#   CELL=x30 REPS="B1 F1 G1" bash submit_2Mb_r3_x20_psamp_fm001.sh   # B/F/G pilot

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SCRATCH=/n/scratch/users/n/njc12/snakemake
PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake
cd "$REPO"

CELL="${CELL:-x20}"
SAMPLING_N="${SAMPLING_N:-8000}"
SIG_P="${SIG_P:-5e-8}"
MIN_POWER="${MIN_POWER:-0.05}"
MIN_POOL_MULTIPLE="${MIN_POOL_MULTIPLE:-2}"
CAUSAL_FLOOR="${CAUSAL_FLOOR:-0}"
FM_FLOOR="${FM_FLOOR:-0.001}"
GLM_FLOOR="${GLM_FLOOR:-0.001}"
LD_CTRL="${LD_CTRL:-0.75}"
REPS="${REPS:-A1 A2 A3 A4 A5 E1 E2 E3 E4 E5}"
DRY="${DRY:-}"
BEGIN="${BEGIN:-}"

OUT_SCRATCH="$SCRATCH/simulations_round_3_2Mb_${CELL}_psamp_${SAMPLING_N}_fm_${FM_FLOOR}"
OUT_PUBLISH="$PUBLISH/simulations_round_3_2Mb_${CELL}_psamp_${SAMPLING_N}_fm_${FM_FLOOR}"

# Stage 1 is the GENETIC simulation. It does not depend on the phenotype
# multipliers or on how the causal variants are drawn, so a replicate has exactly
# one tree sequence, wherever it was first generated. DELIBERATELY NOT
# parameterised by CELL -- interpolating it would make A2-A5/E2-E5 fail pre-flight
# and invite a redundant multi-hour stage-1 rerun per cell.
STAGE1_ROOTS=(
    "$SCRATCH/simulations_round_3_2Mb_g5t20_cmaf_0.001"   # A2-A5, E2-E5
    "$SCRATCH/simulations_round_3_2Mb"                    # A1, E1 (x35 pilot)
)

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

# The whole point of the arm. Without the patch Snakemake would reject the unknown
# --config key, or worse, an older create_gwas_files_and_phenotypes.py would ignore
# the flag and draw uniformly into a directory named for power sampling.
for MARKER in 'causal_sampling:Snakefile' 'causal_sampling:helpers/paths.py' \
              'sampling_flag:rules/common.smk' \
              'select_central_power:create_gwas_files_and_phenotypes.py' \
              'synthetic_dfe_effects:create_gwas_files_and_phenotypes.py' \
              'synthetic_flag:rules/common.smk' \
              'DEFAULT_CAUSAL_MIN_MAF:helpers/paths.py' \
              'cattle_neutral:Snakefile' \
              'stage1_cattle_neutral_ts:helpers/paths.py'; do
    PAT="${MARKER%%:*}"; FILE="${MARKER#*:}"
    grep -q "$PAT" "$REPO/$FILE" 2>/dev/null || {
        echo "ERROR: $REPO/$FILE predates the causal_sampling patch (no '$PAT')." >&2
        echo "       rsync the repo to O2 before running this. Nothing launched." >&2
        exit 1; }
done

# fm_min_maf=0.001 leaves human with essentially all ~2847 candidates per window,
# so DAP-G needs the un-thinned resource tier. The old `fm_min_maf > 0` test would
# hand it a 20 min / 4 GB allocation for work measured at 23:50 and 2580 MB.
grep -q "_fm_thins_candidates" "$REPO/rules/common.smk" 2>/dev/null || {
    echo "ERROR: $REPO/rules/common.smk predates the DAP-G resource patch." >&2
    echo "       At fm_min_maf=$FM_FLOOR every human hgwas locus would burn two" >&2
    echo "       failed attempts before the retry ladder gets it through." >&2
    exit 1; }

# Resolve every reused stage-1 path before submitting anything. A stage-1 path that
# resolves to nothing does not fail loudly: Snakemake simply decides it has to BUILD
# the tree sequence, and a multi-hour genetic simulation starts instead of a reuse.
# The lookup is also seed-strict, so a directory holding some other replicate's .ts
# is silently no help -- hence checking per replicate.
declare -A STAGE1_SRC
for ID in $REPS; do
    CFG=$(config_of "$ID") || { echo "ERROR: unknown replicate '$ID'" >&2; exit 1; }
    [[ -f "$REPO/$CFG" ]] || { echo "ERROR: missing config: $REPO/$CFG" >&2; exit 1; }
    seed_of "$ID" > /dev/null || exit 1

    # causal_sampling=power needs an explicit n_central_traits: its pool is every
    # central selco != 0 variant with no MAF floor, so "use all eligible" would make
    # the whole region causative. The Snakefile refuses at parse time; catch it here
    # instead, before ten controllers queue only to die identically.
    grep -qE '^n_central_traits:[[:space:]]*[0-9]+' "$REPO/$CFG" || {
        echo "ERROR: $CFG sets no n_central_traits, which causal_sampling=power requires." >&2
        echo "       Nothing launched." >&2; exit 1; }

    # Categories H and I build their own stage 1 rather than reusing one. They are
    # distinct pipelines (human_neutral / cattle_neutral, pure coalescents --
    # seconds to minutes, not hours), so there is no x35 predecessor to point at
    # and nothing to save by looking. An empty STAGE1_SRC tells the controller to
    # omit stage1_search_dirs entirely.
    if [[ "${ID:0:1}" =~ ^(H|I)$ ]]; then
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

    # Refuse to relaunch into a workdir that already holds a finished stage 2: the
    # trait draw would be regenerated over results that may already be summarized.
    if compgen -G "$OUT_SCRATCH/${ID}/stage2/*/*/.stage2.done" > /dev/null; then
        echo "ERROR: $OUT_SCRATCH/${ID} already has a completed stage 2." >&2
        echo "       Remove it deliberately or pick another root. Nothing launched." >&2
        exit 1
    fi
done

echo "Pre-flight OK: $(echo $REPS | wc -w) replicates."
echo "  causal_sampling=power  sampling_gwas_n=$SAMPLING_N  sig_p=$SIG_P"
echo "  min_power=$MIN_POWER  min_pool_multiple=$MIN_POOL_MULTIPLE"
echo "  cell=$CELL (gwas_scaling=gtex_scaling=20)  ld_ctrl=$LD_CTRL"
echo "  fm_min_maf=$FM_FLOOR  min_maf=$GLM_FLOOR  causal_min_maf=$CAUSAL_FLOOR (central pool AND flanks)"
echo "  target root: $OUT_SCRATCH"
echo "  stage 1 REUSED; stage 2 RUNS for every replicate."
echo

# ------------------------------------------------------------------- submit
N=0
for ID in $REPS; do
    CFG=$(config_of "$ID")
    SEED=$(seed_of "$ID")
    WD="$OUT_SCRATCH/${ID}"
    PD="$OUT_PUBLISH/${ID}"
    # EXTRA_CONFIG is expanded unquoted by the controller so each space-separated
    # token becomes its own --config arg. NO TOKEN MAY CONTAIN A SPACE. Every floor
    # is passed explicitly even where it equals the default, so the params file
    # records what was INTENDED rather than whatever the config held on the day.
    EXTRA="causal_sampling=power sampling_gwas_n=${SAMPLING_N} sampling_sig_p=${SIG_P}"
    EXTRA="$EXTRA sampling_min_power=${MIN_POWER} sampling_min_pool_multiple=${MIN_POOL_MULTIPLE}"
    EXTRA="$EXTRA causal_min_maf=${CAUSAL_FLOOR} fm_min_maf=${FM_FLOOR} min_maf=${GLM_FLOOR}"
    EXTRA="$EXTRA ld_ctrl=${LD_CTRL}"
    CAT_EXTRA=$(category_extra "$ID")
    [[ -n "$CAT_EXTRA" ]] && EXTRA="${EXTRA} ${CAT_EXTRA}"
    if [[ -n "$DRY" ]]; then
        printf "  DRY %-3s seed %-3s %s\n        STAGE1_SRC=%s\n        EXTRA_CONFIG=%s\n" \
            "$ID" "$SEED" "$CFG" "${STAGE1_SRC[$ID]}" "$EXTRA"
        N=$((N + 1)); continue
    fi
    mkdir -p "$WD"
    JID=$(sbatch --parsable --job-name="${ID}_${CELL}_psamp${SAMPLING_N}_fm${FM_FLOOR}" \
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
FIRST_REP=$(echo $REPS | awk '{print $1}')
echo "CHECK FIRST, before anything downstream is believed -- how many GWAS causal"
echo "loci actually have a GTEx partner. Anything short of 50 means the scorers are"
echo "counting structural gaps as colocalization failures:"
echo "  grep -E 'GTEx partner|power:' $OUT_SCRATCH/*/logs/stage2_split_pheno.log"
echo
echo "Then confirm the draw did what it is for -- selected median power should be"
echo "far above the pool median:"
echo "  grep 'power:' $OUT_SCRATCH/*/logs/stage2_split_pheno.log"
echo
echo "The per-variant draw is auditable in (h* for A/B, c* for E/F/G):"
echo "  $OUT_SCRATCH/${FIRST_REP}/stage2/*.psamp_${SAMPLING_N}/*/[hc]gwas_causal_power_*.tsv"
echo "  $OUT_SCRATCH/${FIRST_REP}/stage2/*.psamp_${SAMPLING_N}/*/[hc]gwas_trait_partners_*.tsv"
if [[ " $REPS " == *" B"* ]]; then
    echo
    echo "*** B IS IN THIS BATCH AND HAS NEVER COMPLETED UNDER POWER SAMPLING. ***"
    echo "The stage-2 crash that blocked it was fixed 2026-08-04; confirm this run"
    echo "reaches 'Finished jobid: 0' and that its partner count above is still 48/50"
    echo "BEFORE launching further B replicates. See this script's header."
fi
