#!/bin/bash
# Power-weighted causal sampling at the x20 cell, fine-mapping and GLM floors 0.001.
#
#   causal_sampling  power    causal loci drawn with inclusion probability
#                             proportional to the power to detect them at n=8000
#   PLATEAU          1        weight = power. Set PLATEAU=0.5 for SATURATING
#                             power sampling: weight = min(power, 0.5), so once a
#                             variant is effectively certain to be found, being
#                             found more surely stops increasing its chance of
#                             being drawn. Changes the root name, so a saturated
#                             arm never lands on top of its unsaturated twin.
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
#
# CATEGORIES O AND P (low heterozygosity). O is A's genome and P is E's, with a
# measured fraction of the NEUTRAL (selco == 0) sites deleted in stage 2, sized so
# pi halves (`neutral_keep_fraction`, helpers/neutral_thinning.py). They are the
# ONLY categories that share another category's seed band -- O{N} runs at A's seed
# 1{N} and P{N} at E's 5{N} -- because they ADOPT their parent's stage-1 tree
# sequence rather than simulating a new one. That is the point: O_i's variant set
# is a strict subset of A_i's, the 50 central causal loci are the same 50
# positions, and A - O isolates variant density with no replicate noise in it.
# stage1_donor_of() below is what redirects the stage-1 lookup. They must always
# rebuild stage 2 (neutral_keep_fraction is in params_record.STAGE2_KEYS), and
# they must never be launched from submit_2Mb_r3_cmaf_replicates.sh, which builds
# stage 1. See "Categories O and P" in README_snakemake.md.
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
# SATURATING POWER SAMPLING. 1 = off (the historical weight, exactly proportional
# to power). At 0.5 every variant with power >= 0.5 carries the same weight, so
# once detection is effectively assured, relative power stops deciding which loci
# become causal. The value is passed straight through to
# sampling_power_plateau; helpers/paths.py turns it into the ".psamp_N_sat05"
# path segment, and SAT_SUFFIX below keeps the root name telling the same story.
PLATEAU="${PLATEAU:-1}"
# HOW MANY TRAIT LOCI. 50/50 is every arm through round 3: 50 causal loci in the
# central 1 Mb (the GWAS traits, and the same 50 shared with GTEx) plus 50
# GTEx-only loci drawn from the 500 kb flanks. Halving both gives an arm with
# half the loci per replicate. The count appears in NO filename -- stage-2 outputs
# are named gwas_G_gtex_T_maf_M only -- so it goes in the root name here and in
# helpers/paths.trait_count_segment, or two arms would collide file-for-file in
# the flat archive.
N_CENTRAL="${N_CENTRAL:-50}"
N_FLANK="${N_FLANK:-50}"
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

# A saturated run and its unsaturated twin differ in WHICH loci become causal, so
# they get different roots. At PLATEAU=1 the suffix is empty and the root is
# byte-identical to what this script has always produced.
if [[ "$PLATEAU" == "1" || "$PLATEAU" == "1.0" ]]; then
    SAT_SUFFIX=""
else
    SAT_SUFFIX="_sat_${PLATEAU}"
fi
if [[ "$N_CENTRAL" == "50" ]]; then NTR_SUFFIX=""; else NTR_SUFFIX="_ntr_${N_CENTRAL}"; fi
OUT_SCRATCH="$SCRATCH/simulations_round_3_2Mb_${CELL}_psamp_${SAMPLING_N}${SAT_SUFFIX}${NTR_SUFFIX}_fm_${FM_FLOOR}"
OUT_PUBLISH="$PUBLISH/simulations_round_3_2Mb_${CELL}_psamp_${SAMPLING_N}${SAT_SUFFIX}${NTR_SUFFIX}_fm_${FM_FLOOR}"

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
#   J{N} -> 10{N}  human, AFRICAN ancestry (A's config + population: AFR)
#   M{N} -> 13{N}  human, FINNISH founder demography (Wang 2014, FIN deme)
#   N{N} -> 14{N}  human, NON-FINNISH European (same model, NFE deme)
# The rule is seed = 10*letter_index + replicate (A=1 ... I=9, J=10), which
# stops reading off as a tens digit at J but keeps every band disjoint.
seed_of() {
    local n="${1:1}"
    case "${1:0:1}" in
        A) echo "1${n}" ;;
        B) echo "2${n}" ;;
        E) echo "5${n}" ;;
        F) echo "6${n}" ;;
        G) echo "7${n}" ;;
        H) echo "8${n}" ;; I) echo "9${n}" ;;
        J) echo "10${n}" ;;
        M) echo "13${n}" ;; N) echo "14${n}" ;;
        # O and P SHARE A's and E's seed bands, and that is not an oversight.
        # Every other category has a private band because it is a different
        # genetic simulation; O and P are the SAME simulation with a fraction of
        # the neutral sites deleted in stage 2, so O{N} adopts A{N}'s tree
        # sequence and P{N} adopts E{N}'s, at those seeds, and O1's variant set
        # is a strict subset of A1's. That pairing is why A - O carries no
        # replicate noise. Consequence: O and P must never run their own stage 1
        # (see stage1_donor_of below), and they must never adopt A's or E's
        # stage 2 (neutral_keep_fraction is in params_record.STAGE2_KEYS, which
        # makes that a loud refusal rather than a silent wrong answer).
        O) echo "1${n}" ;; P) echo "5${n}" ;;
        *) echo "ERROR: no seed convention for '$1'" >&2; return 1 ;;
    esac
}

# Which replicate's stage-1 tree sequence a run adopts. The identity for every
# category except O and P, which reuse A's and E's rather than simulating their
# own -- the pairing described in seed_of(). Their seeds already match, so the
# seed-strict lookup succeeds; only the DIRECTORY has to be redirected, because
# it is named for the replicate rather than for the seed.
stage1_donor_of() {
    case "${1:0:1}" in
        O) echo "A${1:1}" ;;
        P) echo "E${1:1}" ;;
        *) echo "$1" ;;
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
#
# J is category A with `population: AFR` -- the same OutOfAfrica_2T12 run, the
# other branch sampled. It gets its own file rather than a --config override
# (which is how B is done) because it also needs its own basename and seed band,
# and because one file per genetic model is the rule everywhere else here.
# Only the g5t20 cell exists so far; any other CELL fails the pre-flight
# `[[ -f "$REPO/$CFG" ]]` check, which is the intended way to find that out.
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
        A|B) echo "config/human_2Mb_${CELL}_r3.yaml" ;;
        J)   echo "config/human_afr_2Mb_${CELL}_r3.yaml" ;;
        M)   echo "config/human_fin_2Mb_${CELL}_r3.yaml" ;;
        N)   echo "config/human_nfe_2Mb_${CELL}_r3.yaml" ;;
        H)   echo "config/human_neutral_2Mb_${CELL}_r3.yaml" ;;
        I)   echo "config/cattle_neutral_2Mb_${CELL}_r3.yaml" ;;
        E)   echo "config/cattle_baseline_from_midpoint_2Mb_${CELL}_r3.yaml" ;;
        F)   echo "config/cattle_sel_bottlenecked_2Mb_${CELL}_r3.yaml" ;;
        G)   echo "config/cattle_sel_not_bottlenecked_2Mb_${CELL}_r3.yaml" ;;
        O)   echo "config/human_lowhet_2Mb_${CELL}_r3.yaml" ;;
        P)   echo "config/cattle_lowhet_2Mb_${CELL}_r3.yaml" ;;
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
              'stage1_cattle_neutral_ts:helpers/paths.py' \
              'neutral_keep_fraction:helpers/paths.py' \
              'thin_flag:rules/common.smk'; do
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

    # O and P with no neutral_keep_fraction would be byte-identical copies of A
    # and E under a different name. The key is only meaningful once
    # helpers/measure_pi_components.py has measured it, so refuse until it is set.
    if [[ "${ID:0:1}" =~ ^(O|P)$ ]]; then
        grep -qE '^neutral_keep_fraction:[[:space:]]*[0-9.]+' "$REPO/$CFG" || {
            echo "ERROR: $CFG does not set neutral_keep_fraction, so ${ID} would" >&2
            echo "       reproduce its parent category exactly. Measure it with" >&2
            echo "       helpers/measure_pi_components.py and pin it. Nothing launched." >&2
            exit 1; }
    fi

    # Categories H and I build their own stage 1 rather than reusing one. They are
    # distinct pipelines (human_neutral / cattle_neutral, pure coalescents --
    # seconds to minutes, not hours), so there is no x35 predecessor to point at.
    # J runs the SAME pipeline as A but samples the AFR branch, so its tree
    # sequence is hts_AFR_{seed}.ts and no EUR stage1 dir can satisfy it --
    # pointing at one would find nothing, and find nothing silently.
    # An empty STAGE1_SRC tells the controller to omit stage1_search_dirs
    # entirely.
    if [[ "${ID:0:1}" =~ ^(H|I|J|M|N)$ ]]; then
        STAGE1_SRC["$ID"]=""
    else
        S1=""
        DONOR=$(stage1_donor_of "$ID")   # itself, except O->A and P->E
        for ROOT in "${STAGE1_ROOTS[@]}"; do
            if compgen -G "$ROOT/${DONOR}/stage1/*.ts" > /dev/null; then S1="$ROOT/${DONOR}/stage1"; break; fi
        done
        [[ -n "$S1" ]] || {
            echo "ERROR: no stage-1 tree sequence for ${ID} (donor ${DONOR}) in any of:" >&2
            for ROOT in "${STAGE1_ROOTS[@]}"; do echo "       $ROOT/${DONOR}/stage1" >&2; done
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

# The trailing `:` is load-bearing under `set -e`. The loop's exit status is
# its LAST iteration's, and `[[ ... ]] && printf` returns 1 whenever that
# replicate does not match -- so `VAR=$(...)` inherits a non-zero status and
# kills the script before the submit loop, silently, whenever the letter is
# absent from REPS. That is exactly the common case for this line.
THINNED=$(for r in $REPS; do [[ "${r:0:1}" =~ ^(O|P)$ ]] && printf ' %s' "$r"; done; :)
if [[ -n "$THINNED" ]]; then
    echo "Low-heterozygosity categories in this launch:${THINNED}"
    for r in $THINNED; do
        K=$(awk -F': *' '/^neutral_keep_fraction:/{print $2; exit}' "$REPO/$(config_of "$r")")
        printf "  %-3s stage 1 adopted from %s; keeps %s of the selco == 0 sites\n" \
            "$r" "$(stage1_donor_of "$r")" "$K"
    done
    echo "  The causal draw is UNAFFECTED: power weights are built from selco != 0"
    echo "  variants, which thinning never touches, so the drawn loci match the parent."
fi
echo "Pre-flight OK: $(echo $REPS | wc -w) replicates."
echo "  causal_sampling=power  sampling_gwas_n=$SAMPLING_N  sig_p=$SIG_P"
echo "  min_power=$MIN_POWER  min_pool_multiple=$MIN_POOL_MULTIPLE"
CFG_PROBE=$(config_of "${REPS%% *}")
GW_S=$(awk '/^gwas_scaling:/{print $2; exit}' "$REPO/$CFG_PROBE")
GT_S=$(awk '/^gtex_scaling:/{print $2; exit}' "$REPO/$CFG_PROBE")
echo "  cell=$CELL (gwas_scaling=$GW_S  gtex_scaling=$GT_S)  ld_ctrl=$LD_CTRL"
if [[ "$PLATEAU" == "1" || "$PLATEAU" == "1.0" ]]; then
    echo "  weight = power (no plateau)"
else
    echo "  weight = min(power, $PLATEAU)  SATURATING -- past $PLATEAU, extra power buys no extra representation"
fi
echo "  fm_min_maf=$FM_FLOOR  min_maf=$GLM_FLOOR  causal_min_maf=$CAUSAL_FLOOR (central pool AND flanks)"
echo "  trait loci: $N_CENTRAL central (GWAS + shared GTEx) + $N_FLANK GTEx-only flank"
echo "  pool guard needs $((MIN_POOL_MULTIPLE * N_CENTRAL)) variants at power >= $MIN_POWER"
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
    EXTRA="$EXTRA sampling_power_plateau=${PLATEAU}"
    EXTRA="$EXTRA n_central_traits=${N_CENTRAL} n_flank_gtex_traits=${N_FLANK}"
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
echo "loci actually have a GTEx partner. Anything short of $N_CENTRAL means the scorers are"
echo "counting structural gaps as colocalization failures:"
echo "  grep -E 'GTEx partner|power:' $OUT_SCRATCH/*/logs/stage2_split_pheno.log"
echo
echo "Then confirm the draw did what it is for -- selected median power should be"
echo "far above the pool median:"
echo "  grep 'power:' $OUT_SCRATCH/*/logs/stage2_split_pheno.log"
echo
echo "The per-variant draw is auditable in (h* for A/B, c* for E/F/G):"
echo "  $OUT_SCRATCH/${FIRST_REP}/stage2/*.psamp_${SAMPLING_N}*/*/[hc]gwas_causal_power_*.tsv"
echo "  $OUT_SCRATCH/${FIRST_REP}/stage2/*.psamp_${SAMPLING_N}*/*/[hc]gwas_trait_partners_*.tsv"
if [[ " $REPS " == *" B"* ]]; then
    echo
    echo "*** B IS IN THIS BATCH AND HAS NEVER COMPLETED UNDER POWER SAMPLING. ***"
    echo "The stage-2 crash that blocked it was fixed 2026-08-04; confirm this run"
    echo "reaches 'Finished jobid: 0' and that its partner count above is still 48/50"
    echo "BEFORE launching further B replicates. See this script's header."
fi
