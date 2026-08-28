#!/bin/bash
#
# Launch one ARM of the publication simulation set: every category, every
# replicate, one controller each.
#
#   ARM=causal_maf001_paired bash submit_publication.sh
#   ARM=causal_power_n8000 REPS="1" UNTIL=stage2 JOBS=4 bash submit_publication.sh
#   DRY=1 ARM=... bash submit_publication.sh          # pre-flight only, no sbatch
#
# Env:
#   ARM      required; a row of helpers/publication_arms.tsv
#   REPS     replicates, default "1 2 3 4 5"
#   CATS     category letters, default every row of publication_categories.tsv
#   JOBS     max concurrent SLURM jobs PER CONTROLLER, default 8
#   UNTIL    optional snakemake --until target
#   DRY      1 = pre-flight and print, submit nothing
#   CELL     multiplier cell, default g5t20
#   MAX_QUEUE  refuse to launch if the account queue plus this wave exceeds it
#              (default 6000; account MaxSubmit is 10000 and is SHARED)
#
# WHY ONE ARM PER INVOCATION. 120 controllers at the profile's `jobs: 200` would
# queue 24,000 against an account MaxSubmit of 10,000, on a shared account. A
# refused sbatch is a failed job to Snakemake's SLURM executor; it retries twice,
# gives up, and with keep-going: True the replicate ends with no .enloc.sig.out --
# which reads downstream as "this arm colocalized nothing" rather than as a gap.
# One arm at a time, at a modest -j, is how that does not happen.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${REPO:-$HERE}"
source "$REPO/lib/publication_common.sh"

: "${ARM:?set ARM to one of: $(tsv_keys "$ARMS_TSV" | tr '\n' ' ')}"
REPS="${REPS:-1 2 3 4 5}"
CATS="${CATS:-$(tsv_keys "$CATEGORIES_TSV" | tr '\n' ' ')}"
JOBS="${JOBS:-8}"
UNTIL="${UNTIL:-}"
DRY="${DRY:-0}"
MAX_QUEUE="${MAX_QUEUE:-6000}"

STAGE1_INPUTS="${STAGE1_INPUTS:-$PUBLISH_ROOT/stage1_inputs}"

arm_field "$ARM" arm >/dev/null || {
    echo "ERROR: unknown arm '$ARM'. Known: $(tsv_keys "$ARMS_TSV" | tr '\n' ' ')" >&2
    exit 1
}

# Keys the arm table SETS, passed as --config. Everything the arm varies is here.
PASS_KEYS="causal_min_maf causal_sampling sampling_gwas_n sampling_sig_p
           sampling_min_power sampling_power_plateau sampling_min_pool_multiple
           require_gtex_partner fm_min_maf min_maf n_central_traits
           n_flank_gtex_traits ld_ctrl dapg_window fastenloc_prior params_strict"

# Keys the arm table RECORDS but does not set: they must already be right in the
# config YAML. Asserted rather than passed, so there is one source of truth and a
# YAML that drifted is caught at pre-flight instead of producing a run whose
# params file disagrees with ARMS.tsv.
ASSERT_KEYS="gwas_scaling gtex_scaling gtex_size gwas_size L"

echo "=============================================================="
echo " arm         $ARM"
echo "             $(arm_field "$ARM" description)"
echo " categories  $CATS"
echo " replicates  $REPS"
echo " cell        $CELL"
echo " repo        $REPO"
echo " workdirs    $SCRATCH_ROOT/$ARM/<category>_rep<N>"
echo " publish to  $PUBLISH_ROOT/$ARM/<category>_rep<N>"
echo " stage 1     $STAGE1_INPUTS (adopted; nothing re-simulated)"
echo " -j          $JOBS per controller${UNTIL:+   --until $UNTIL}"
echo "=============================================================="

FAIL=0

echo
echo "-- pre-flight: deployed code --"
preflight_repo || FAIL=1

echo
echo "-- pre-flight: per-run --"
declare -a RUN_IDS=() RUN_CFG=() RUN_SEED=() RUN_WD=() RUN_PD=() RUN_EXTRA=() RUN_S1=() RUN_DIR_NAMES=()

for L in $CATS; do
    CFG=$(config_of "$L") || { echo "  ERROR: unknown category '$L'" >&2; FAIL=1; continue; }
    if [[ ! -f "$REPO/$CFG" ]]; then
        echo "  ERROR: $L: config not found: $REPO/$CFG" >&2; FAIL=1; continue
    fi
    preflight_adoption "$CFG" || FAIL=1

    # The recorded-but-not-passed keys must already match the YAML.
    for k in $ASSERT_KEYS; do
        want=$(arm_field "$ARM" "$k")
        got=$(awk -v k="$k" -F': *' '$1 == k { print $2; exit }' "$REPO/$CFG" | sed 's/ *#.*//; s/[[:space:]]*$//')
        # Numeric compare so 2000000 and 2.0e6 do not read as different.
        if ! awk -v a="$want" -v b="$got" 'BEGIN { exit !(a+0 == b+0 && a != "" && b != "") }'; then
            echo "  ERROR: $L: $CFG has $k=$got but ARMS.tsv records $want" >&2
            FAIL=1
        fi
    done

    for R in $REPS; do
        ID="${L}${R}"
        SEED=$(seed_of "$ID") || { FAIL=1; continue; }
        DIR="$(category_dir_of "$L")_rep${R}"
        EXTRA=""
        for k in $PASS_KEYS; do
            v=$(arm_field "$ARM" "$k")
            [[ -n "$v" ]] && EXTRA+="${k}=${v} "
        done
        CE=$(category_extra "$L")
        [[ -n "$CE" ]] && EXTRA+="$CE "
        # NO TOKEN MAY CONTAIN A SPACE: EXTRA is expanded unquoted on the
        # snakemake command line so each key=value becomes its own --config arg.
        if grep -q "= \|=$" <<<"$EXTRA"; then
            echo "  ERROR: $ID: empty or space-containing --config token in: $EXTRA" >&2
            FAIL=1
        fi

        RUN_IDS+=("$ID"); RUN_CFG+=("$CFG"); RUN_SEED+=("$SEED")
        RUN_WD+=("$SCRATCH_ROOT/$ARM/$DIR"); RUN_PD+=("$PUBLISH_ROOT/$ARM/$DIR")
        RUN_EXTRA+=("$EXTRA"); RUN_S1+=("$STAGE1_INPUTS")
        RUN_DIR_NAMES+=("$DIR")
    done
done

# Stage 1 must resolve BY EXACT FILENAME. This is the only guard for the human
# categories: search_dirs' seed pattern matches sd(\d+) and seed_(\d+), and
# `hts_11.ts` matches neither, so _extract_seed returns None and the
# seed-collision check never fires. A wrong human tree sequence adopts silently.
echo
echo "-- pre-flight: stage-1 inputs --"
MANIFEST="$PUBLISH_ROOT/RUNS.tsv"
if [[ -f "$MANIFEST" ]]; then
    missing=0
    # Only the runs THIS invocation submits. Checking the whole arm made a
    # single-replicate wave refuse because replicates 2-5 were not built yet,
    # which is the wrong answer for `REPS=1` -- and a pre-flight that refuses
    # correct requests gets bypassed.
    while IFS=$'\t' read -r m_arm m_dir m_s1; do
        [[ "$m_arm" == "$ARM" ]] || continue
        local_wanted=0
        for d in "${RUN_DIR_NAMES[@]}"; do [[ "$d" == "$m_dir" ]] && local_wanted=1 && break; done
        (( local_wanted )) || continue
        if [[ ! -f "$STAGE1_INPUTS/$m_s1" ]]; then
            echo "  MISSING $m_s1  (needed by $m_dir)" >&2; missing=1
        fi
    done < <(awk -F'\t' 'NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
                         { print $h["arm"] "\t" $h["run_dir"] "\t" $h["stage1_file"] }' "$MANIFEST")
    if (( missing )); then
        echo "  ERROR: stage-1 inputs incomplete. An unresolved stage-1 path does" >&2
        echo "         not fail loudly -- Snakemake decides it must BUILD the tree" >&2
        echo "         sequence and a multi-hour genetic simulation starts." >&2
        FAIL=1
    else
        echo "  all stage-1 files for the ${#RUN_IDS[@]} runs in this wave are present"
    fi
else
    echo "  ERROR: $MANIFEST not found. Run helpers/write_run_manifests.py first --" >&2
    echo "         it is what turns paths.py's own filename builders into the" >&2
    echo "         expected-name list this check uses." >&2
    FAIL=1
fi

echo
echo "-- pre-flight: no controller already on these workdirs --"
# A second controller on a live workdir is not a no-op: controller_publication
# runs `snakemake --unlock || true` first, which would release the RUNNING
# controller's lock and let two DAGs write the same outputs. Detected from the
# filesystem, not squeue, because O2's slurmctld goes unresponsive under load and
# a guard that depends on it is a guard that disappears exactly when the cluster
# is busy enough to need it.
locked=0
for i in "${!RUN_IDS[@]}"; do
    if compgen -G "${RUN_WD[$i]}/.snakemake/locks/*" >/dev/null 2>&1; then
        echo "  ${RUN_IDS[$i]}: ${RUN_WD[$i]} holds a snakemake lock" >&2
        locked=1
    fi
done
if (( locked )); then
    echo "  A controller is already running on at least one of these workdirs." >&2
    echo "  Cancel it, or narrow REPS/CATS to exclude it." >&2
    FAIL=1
else
    echo "  clear"
fi

echo
echo "-- pre-flight: nothing to overwrite --"
existing=0
for i in "${!RUN_IDS[@]}"; do
    if compgen -G "${RUN_WD[$i]}/stage4/*/*.enloc.sig.out" >/dev/null 2>&1; then
        echo "  ${RUN_IDS[$i]}: ${RUN_WD[$i]} already holds stage-4 output" >&2
        existing=1
    fi
done
(( existing )) && { echo "  refusing to overwrite finished output; move it aside or pick another ARM" >&2; FAIL=1; }
(( existing )) || echo "  clear"

echo
echo "-- pre-flight: queue headroom --"
if command -v squeue >/dev/null 2>&1; then
    # Bounded. O2's login-node slurmctld goes unresponsive under load -- squeue
    # can hang for minutes -- and an unbounded call here hung the launcher itself,
    # which is strictly worse than having no headroom check: the operator sees
    # nothing happen and cannot tell a refusal from a stall.
    ACCT=$(timeout 20 sacctmgr -n -P show assoc where user="$USER" format=Account 2>/dev/null | head -1)
    QUEUED=$(timeout 20 squeue -h -u "$USER" 2>/dev/null | wc -l | tr -d ' ')
    if [[ -z "$QUEUED" || "$QUEUED" == "0" ]] && ! timeout 20 squeue -h -u "$USER" >/dev/null 2>&1; then
        echo "  squeue did not answer within 20s -- scheduler is loaded." >&2
        echo "  SKIPPING the headroom check. The account GrpTRES cap bounds this" >&2
        echo "  regardless of what we submit; SLURM queues the excess." >&2
        QUEUED=""
    fi
    WANT=$(( ${#RUN_IDS[@]} * JOBS + ${#RUN_IDS[@]} ))
    echo "  account=${ACCT:-?}  currently queued/running (this user): ${QUEUED:-unknown}"
    echo "  this wave would add at most $WANT (${#RUN_IDS[@]} controllers x $JOBS + controllers)"
    if [[ -n "$QUEUED" ]] && (( QUEUED + WANT > MAX_QUEUE )); then
        echo "  ERROR: $((QUEUED + WANT)) exceeds MAX_QUEUE=$MAX_QUEUE. The account cap is" >&2
        echo "         10000 and shared; lower JOBS, split REPS, or wait." >&2
        FAIL=1
    fi
else
    echo "  (no squeue here -- not on O2; skipping)"
fi

if (( FAIL )); then
    echo
    echo "PRE-FLIGHT FAILED -- nothing submitted." >&2
    exit 1
fi

echo
echo "-- pre-flight passed: ${#RUN_IDS[@]} runs --"

if [[ "$DRY" == "1" ]]; then
    for i in "${!RUN_IDS[@]}"; do
        printf "  DRY %-3s seed %-3s %s\n        WD=%s\n        EXTRA=%s\n" \
            "${RUN_IDS[$i]}" "${RUN_SEED[$i]}" "${RUN_CFG[$i]}" "${RUN_WD[$i]}" "${RUN_EXTRA[$i]}"
    done
    echo
    echo "DRY=1: nothing submitted."
    exit 0
fi

echo
for i in "${!RUN_IDS[@]}"; do
    JID=$(sbatch --parsable \
        --job-name="${RUN_IDS[$i]}_${ARM}" \
        --output="${RUN_WD[$i]}/controller_%j.out" \
        --error="${RUN_WD[$i]}/controller_%j.err" \
        --export=ALL,REPO="${REPO}",CONFIGFILE="${RUN_CFG[$i]}",SEED="${RUN_SEED[$i]}",WD="${RUN_WD[$i]}",PD="${RUN_PD[$i]}",STAGE1_SRC="${RUN_S1[$i]}",EXTRA_CONFIG="${RUN_EXTRA[$i]}",JOBS="${JOBS}",UNTIL="${UNTIL}" \
        "$REPO/controller_publication.sbatch")
    printf "  %-3s seed %-3s controller %s\n" "${RUN_IDS[$i]}" "${RUN_SEED[$i]}" "$JID"
done

echo
echo "Submitted ${#RUN_IDS[@]} controllers for arm $ARM."
echo "Watch:   squeue -u $USER -o '%.10i %.30j %.8T %.10M'"
echo "Publish: ARM=$ARM bash publish_to_data2.sh   (gated; refuses partial runs)"
