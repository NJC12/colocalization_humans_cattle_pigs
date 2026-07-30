#!/bin/bash
# Fine-mapping variants of the round-3 A1 + E1 multiplier grid.
#
# The multiplier cells that already exist (x5 x10 x15 x20 x35 g5t20) were all
# fine-mapped one way: DAP-G clustering at r2 >= 0.75, over the UNFILTERED
# variant set. This script reruns those same cells under the other three
# corners of that 2x2:
#
#   r25        ld_ctrl=0.25                        (looser clustering)
#   maf01      fm_min_maf=0.01                     (only common variants fine-mapped)
#   r25maf01   ld_ctrl=0.25 fm_min_maf=0.01        (both)
#
# Why it matters: in a +/-250 kb window human carries 2847 candidates of which
# 2516 (88%) are MAF < 0.01, against cattle's 366 of which 34 (9%). At
# MAF >= 0.01 the pools are 327 vs 330 -- essentially identical. Under DAP-G's
# exchangeable prior the unfiltered runs therefore hand cattle a ~1 log10 unit
# prior-odds advantage that has nothing to do with demography. The maf01 arms
# test whether the cattle > human colocalization ordering survives once that is
# removed.
#
# STAGES 1 AND 2 ARE NOT RERUN. Each variant reuses the source cell's own
# stage-1 tree sequence (stage1_search_dirs) and its whole stage-2 output
# (stage2_search_dirs, all-or-nothing symlink adoption in Snakefile). Only
# stages 3-5 execute. That is both a large saving -- stage 2 writes ~600 MB of
# SBAMS per cell -- and what keeps the comparison clean: every variant of a
# given cell shares byte-identical genetics and phenotypes, so any difference
# in the results is attributable to the fine-mapping settings alone.
#
# The MAF filter is applied in stage3_index_geno when building geno.sbams.gz,
# which is exactly the file DAP-G reads, so stage 2 stays reusable. Resource
# requests drop automatically for the human arm (see _dapg_runtime and
# _dapg_mem_mb_base in rules/common.smk); cattle's are left alone because
# filtering removes only ~10% of its candidates.
#
# Run from a login node:
#   bash submit_2Mb_r3_finemap_variants.sh                 # all 3 variants, all 6 cells
#   VARIANTS="r25" bash submit_2Mb_r3_finemap_variants.sh  # batch 1 only
#   CELLS="x35" VARIANTS="maf01" bash submit_2Mb_r3_finemap_variants.sh

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SCRATCH=/n/scratch/users/n/njc12/snakemake
PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake
cd "$REPO"

CELLS="${CELLS:-x5 x10 x15 x20 x35 g5t20}"
VARIANTS="${VARIANTS:-r25 maf01 r25maf01}"

# Every cell's stage-1 tree sequence lives in the x35 root and nowhere else:
# the multiplier runs were launched with stage1_search_dirs pointing here, so
# their own roots have no stage1/ directory at all. Stage 1 is therefore always
# sourced from here regardless of which cell is being rerun. (Stage 2 is
# per-cell -- the phenotypes depend on the multipliers -- and is resolved
# separately below.)
STAGE1_ROOT="$SCRATCH/simulations_round_3_2Mb"

# Optional SLURM --begin for the controllers, e.g. BEGIN=now+90minutes.
#
# Each controller runs Snakemake with `jobs: 200` (profiles/o2/config.yaml), so
# 12 controllers can hold ~2400 jobs in the queue and all 36 would approach the
# 10000 MaxSubmit association limit -- at which point sbatch starts refusing
# and a controller dies mid-DAG. Launching one variant per wave, each deferred
# ~90 min, keeps the peak near one wave's worth. The waves overlap harmlessly
# because DAP-G loci are minutes long, not hours.
BEGIN="${BEGIN:-}"

# id | seed | config stem; the full path is config/<stem><cell suffix>_r3.yaml
REPLICATES=(
  "A1|11|human_2Mb"
  "E1|51|cattle_baseline_from_midpoint_2Mb"
)

# x35 is the original pilot: its configs and output root carry no cell label.
# Every other cell is suffixed. Keeping this in one function stops the config
# path and the source-root path from drifting apart.
cell_suffix() {
    case "$1" in
        x35) echo "" ;;
        *)   echo "_$1" ;;
    esac
}

variant_config() {
    case "$1" in
        r25)      echo "ld_ctrl=0.25" ;;
        maf01)    echo "fm_min_maf=0.01" ;;
        r25maf01) echo "ld_ctrl=0.25 fm_min_maf=0.01" ;;
        *) echo "ERROR: unknown variant '$1'" >&2; return 1 ;;
    esac
}

# ---------------------------------------------------------------- pre-flight
# Validate every config, every reused stage-1 tree sequence and every reused
# stage-2 directory BEFORE submitting anything. Without the stage-1 guard a
# wrong path silently launches a real genetic simulation instead of reusing
# one; without the stage-2 guard the adoption quietly falls through and stage 2
# reruns, which would still be correct but wastes hours and ~600 MB per cell.
declare -A STAGE2_SRC
for CELL in $CELLS; do
    SFX=$(cell_suffix "$CELL")
    SRC_ROOT="$SCRATCH/simulations_round_3_2Mb${SFX}"
    for row in "${REPLICATES[@]}"; do
        IFS='|' read -r ID SEED CFGSTEM <<< "$row"
        CFG="config/${CFGSTEM}${SFX}_r3.yaml"
        [[ -f "$REPO/$CFG" ]] || { echo "ERROR: missing config: $REPO/$CFG" >&2; exit 1; }

        S1_DIR="$STAGE1_ROOT/${ID}/stage1"
        compgen -G "$S1_DIR/*.ts" > /dev/null || {
            echo "ERROR: no stage-1 tree sequence for ${ID} in $S1_DIR" >&2
            echo "Run submit_2Mb_r3_pilot.sh first. Aborting -- no jobs launched." >&2
            exit 1; }

        # stage2_search_dirs wants the dir CONTAINING gwas_*_gtex_*_maf_*/,
        # i.e. <root>/<id>/stage2/<stage1 name>/.
        S2_DIR=$(ls -d "$SRC_ROOT/${ID}/stage2"/*/ 2>/dev/null | head -1 || true)
        [[ -n "$S2_DIR" && -d "$S2_DIR" ]] || {
            echo "ERROR: no stage-2 output for ${ID} under $SRC_ROOT/${ID}/stage2" >&2
            exit 1; }
        S2_DIR="${S2_DIR%/}"
        compgen -G "$S2_DIR/*/.stage2.done" > /dev/null || {
            echo "ERROR: stage 2 for ${ID} in $S2_DIR is not marked complete" >&2
            echo "Adoption is all-or-nothing; a partial dir would silently rerun stage 2." >&2
            exit 1; }
        STAGE2_SRC["${CELL}|${ID}"]="$S2_DIR"
    done
done
for V in $VARIANTS; do variant_config "$V" > /dev/null; done
echo "Pre-flight OK: $(echo $CELLS | wc -w) cells x $(echo $VARIANTS | wc -w) variants x ${#REPLICATES[@]} arms."
echo

# ------------------------------------------------------------------- submit
N=0
for V in $VARIANTS; do
    VCFG=$(variant_config "$V")
    echo "=== variant ${V}  (${VCFG}) ==="
    for CELL in $CELLS; do
        SFX=$(cell_suffix "$CELL")
        SRC_ROOT="$SCRATCH/simulations_round_3_2Mb${SFX}"
        OUT_SCRATCH="$SCRATCH/simulations_round_3_2Mb${SFX}_${V}"
        OUT_PUBLISH="$PUBLISH/simulations_round_3_2Mb${SFX}_${V}"
        for row in "${REPLICATES[@]}"; do
            IFS='|' read -r ID SEED CFGSTEM <<< "$row"
            CFG="config/${CFGSTEM}${SFX}_r3.yaml"
            S2_DIR="${STAGE2_SRC["${CELL}|${ID}"]}"
            WD="$OUT_SCRATCH/${ID}"
            PD="$OUT_PUBLISH/${ID}"
            mkdir -p "$WD"
            # EXTRA_CONFIG is expanded unquoted by the controller so each
            # space-separated token becomes its own --config arg. No token may
            # contain a space -- the list literal below deliberately has none.
            EXTRA="${VCFG} stage2_search_dirs=['${S2_DIR}']"
            JID=$(sbatch --parsable --job-name="${ID}_${CELL}_${V}" \
                ${BEGIN:+--begin="$BEGIN"} \
                --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
                --export=ALL,CONFIGFILE="${CFG}",SEED="${SEED}",WD="${WD}",PD="${PD}",STAGE1_SRC="${STAGE1_ROOT}/${ID}/stage1",EXTRA_CONFIG="${EXTRA}" \
                controller_2Mb_x10.sbatch)
            printf "  %-3s %-6s %-9s controller %s\n" "$ID" "$CELL" "$V" "$JID"
            N=$((N + 1))
        done
    done
    echo
done
echo "Submitted ${N} controllers."
echo "Watch with: squeue -u \$USER -o '%.10i %.24j %.9P %.2t %.10M %R'"
