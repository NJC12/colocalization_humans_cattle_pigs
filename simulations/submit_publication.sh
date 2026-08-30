#!/bin/bash
#
# Launch one ARM of the publication simulation set: every category, every
# replicate, one controller each.
#
#   ARM=causal_maf001_paired bash submit_publication.sh
#   ARM=causal_power_n30000 REPS="1" UNTIL=stage2 JOBS=4 bash submit_publication.sh
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
# From the deep-history table, not a literal: the manifest's --reps default comes
# from the same place, so a wave and a manifest cannot silently disagree about how
# many replicates exist.
REPS="${REPS:-$(replicate_list)}"
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

# One launcher per arm at a time. Two concurrent launchers for the same arm is
# the proximate cause of the controller collision: each passed the live-lock
# check minutes apart, and both submitted the same runs. mkdir is atomic; the
# claim is released on exit however this script ends.
LAUNCH_CLAIM="${TMPDIR:-/tmp}/.publaunch_${ARM}.claim"
if ! mkdir "$LAUNCH_CLAIM" 2>/dev/null; then
    echo "ERROR: another launcher is already running for arm $ARM" >&2
    echo "       ($LAUNCH_CLAIM exists). Wait for it, or remove that directory" >&2
    echo "       if you are certain no launcher is alive." >&2
    exit 1
fi
trap 'rmdir "$LAUNCH_CLAIM" 2>/dev/null || true' EXIT

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
declare -a RUN_CB=() RUN_SPECIES=()

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

        # The deep history is per-REPLICATE, so it cannot live in EXTRA (which is
        # the arm, by contract) or in category_extra (which is the category).
        # Human never gets it: nothing on the human path reads it, and a config
        # key that is silently ignored is how a run ends up not being the run
        # someone thinks they ran.
        SP=$(species_of "$L") || { FAIL=1; continue; }
        CB=""
        if [[ "$SP" == "cattle" ]]; then
            CB=$(cattle_baseline_seed_of "$R") || { FAIL=1; continue; }
        fi
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
        RUN_CB+=("$CB"); RUN_SPECIES+=("$SP")
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
    # Iterate the WAVE, not the manifest. The previous version walked manifest
    # rows and filtered them to the wave, so a run with NO row was never checked
    # -- and it still printed a pass counting ${#RUN_IDS[@]}. With the manifest at
    # its old 5-replicate default and a 20-replicate wave, that reported all 120
    # runs healthy having verified 30.
    #
    # One awk pass builds "run_dir -> stage1_file, marks_file, cattle_baseline_seed"
    # for this arm; bash 3.2 has no associative arrays, and the test suite sources
    # this library, so a temp file is the portable way to carry it.
    IDX=$(mktemp)
    awk -F'\t' -v arm="$ARM" '
        NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
        $h["arm"] == arm {
            print $h["run_dir"] "\t" $h["stage1_file"] "\t" \
                  $h["stage1_marks_file"] "\t" $h["cattle_baseline_seed"] }
    ' "$MANIFEST" > "$IDX"

    missing=0; checked=0
    for i in "${!RUN_DIR_NAMES[@]}"; do
        d="${RUN_DIR_NAMES[$i]}"
        row=$(awk -F'\t' -v d="$d" '$1 == d { print; exit }' "$IDX")
        if [[ -z "$row" ]]; then
            echo "  NO MANIFEST ROW for $d" >&2; missing=1; continue
        fi
        m_s1=$(cut -f2 <<<"$row"); m_marks=$(cut -f3 <<<"$row"); m_cb=$(cut -f4 <<<"$row")

        [[ -f "$STAGE1_INPUTS/$m_s1" ]] || { echo "  MISSING $m_s1  (needed by $d)" >&2; missing=1; }

        # cattle_sel needs BOTH; rules/stage1_cattle_sel.smk requires full AND
        # marks, and a present .full.ts with a missing marks file does not fail --
        # it falls through to a fresh SLiM run inside a 4 GB controller.
        if [[ "$m_marks" != "-" && ! -f "$STAGE1_INPUTS/$m_marks" ]]; then
            echo "  MISSING $m_marks  (needed by $d)" >&2; missing=1
        fi

        # The launcher and the manifest must name the SAME deep history, or the
        # run adopts a stage-1 file built from a different ancestry.
        want_cb="${RUN_CB[$i]}"; [[ -z "$want_cb" ]] && want_cb="-"
        if [[ "$m_cb" != "$want_cb" ]]; then
            echo "  DEEP HISTORY MISMATCH for $d: manifest says $m_cb, launcher says $want_cb" >&2
            echo "    regenerate RUNS.tsv -- one of them is stale" >&2
            missing=1
        fi
        checked=$((checked + 1))
    done
    rm -f "$IDX"

    if (( missing )); then
        echo "  ERROR: stage-1 inputs incomplete. An unresolved stage-1 path does" >&2
        echo "         not fail loudly -- Snakemake decides it must BUILD the tree" >&2
        echo "         sequence and a multi-hour genetic simulation starts." >&2
        FAIL=1
    else
        echo "  all stage-1 files for the $checked runs in this wave are present"
        echo "  (checked $checked of ${#RUN_IDS[@]} wave runs against the manifest)"
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
# SKIP the locked runs rather than refusing the whole arm. Refusing was the
# first version and it is the wrong shape: a launch truncated partway through
# (a client timeout, a dropped connection) leaves some runs live and some
# unstarted, and the only way to finish the arm is then to hand-pick the
# remainder -- which this launcher cannot even express, since CATS x REPS is a
# cross product, not an arbitrary set. Skipping makes a relaunch idempotent,
# which is what you actually need after a truncated one.
#
# Only a LIVE lock is skipped. A run that stopped at stage 2 (Wave 0) holds no
# lock and must be relaunched to continue, which is the common case.
declare -a KEEP_IDX=()
skipped_locked=0
for i in "${!RUN_IDS[@]}"; do
    if compgen -G "${RUN_WD[$i]}/.snakemake/locks/*" >/dev/null 2>&1; then
        printf "  SKIP %-3s a controller already holds this workdir\n" "${RUN_IDS[$i]}"
        skipped_locked=$((skipped_locked+1))
    else
        KEEP_IDX+=("$i")
    fi
done
if (( skipped_locked )); then
    echo "  skipping $skipped_locked run(s) already in progress; will submit the rest"
else
    echo "  clear"
fi

echo
echo "-- pre-flight: nothing to overwrite --"
# SKIP the finished, RESUME the partial. The predicate used to be "holds ANY
# stage-4 output -> refuse the whole wave", which conflated two opposite cases.
#
# `keep-going: True` means a controller can exit COMPLETED having lost one panel
# to a walltime kill, leaving 1 of 2 .enloc.sig.out. That run is exactly the one
# that must be relaunched -- Snakemake rebuilds only what is missing, so a
# relaunch is a cheap resume -- and the old check refused it, which is how a
# truncated wave became unfinishable without hand-editing the launcher. Measured:
# block 2 (replicates 6-10) came back with 30 of 120 runs in that state.
#
# A run with the FULL set is skipped rather than refused, matching the live-lock
# skip above: refusing the whole arm because part of it is already done is the
# wrong answer to a relaunch.
ENLOC_PER_RUN=2      # gwas x {gtex, gtex_smaller}; publish_to_data2.sh wants the same 2
declare -a RESUME_IDX=()
skipped_done=0; resuming=0
for i in "${KEEP_IDX[@]}"; do
    # Counted with a glob array, NOT `ls | wc -l`. This script runs under
    # `set -euo pipefail`: when a run has ZERO enloc files `ls` exits non-zero,
    # pipefail propagates it, the assignment fails and -e exits the script
    # SILENTLY mid-loop, having submitted nothing. That is the same silent-exit
    # trap the queue-headroom check below documents, and it happened here too --
    # the log ended after the last matching run with no error anywhere. An
    # unmatched glob stays literal, so `-f` is false and n stays 0.
    ENLOC_HITS=( "${RUN_WD[$i]}"/stage4/*/*.enloc.sig.out )
    n=0
    for f in "${ENLOC_HITS[@]}"; do [[ -f "$f" ]] && n=$((n+1)); done
    if [[ "$n" -ge "$ENLOC_PER_RUN" ]]; then
        printf "  SKIP %-4s already has the full stage-4 output (%s/%s)\n" \
            "${RUN_IDS[$i]}" "$n" "$ENLOC_PER_RUN"
        skipped_done=$((skipped_done+1))
    else
        [[ "$n" -gt 0 ]] && {
            printf "  RESUME %-4s partial stage-4 output (%s/%s); Snakemake will finish it\n" \
                "${RUN_IDS[$i]}" "$n" "$ENLOC_PER_RUN"
            resuming=$((resuming+1))
        }
        RESUME_IDX+=("$i")
    fi
done
KEEP_IDX=("${RESUME_IDX[@]+"${RESUME_IDX[@]}"}")
if (( skipped_done || resuming )); then
    echo "  skipping $skipped_done finished run(s), resuming $resuming partial one(s)"
else
    echo "  clear"
fi

echo
echo "-- pre-flight: queue headroom --"
if command -v squeue >/dev/null 2>&1; then
    # Bounded. O2's login-node slurmctld goes unresponsive under load -- squeue
    # can hang for minutes -- and an unbounded call here hung the launcher itself,
    # which is strictly worse than having no headroom check: the operator sees
    # nothing happen and cannot tell a refusal from a stall.
    # `|| true` is load-bearing, not defensive noise. This script runs under
    # `set -euo pipefail`, and `timeout` exits 124 when it fires. With pipefail
    # that makes the whole pipeline fail, which makes the assignment fail, which
    # with -e exits the script -- SILENTLY, mid-pre-flight, having submitted
    # nothing. That is exactly what happened: both arm launchers died here with
    # the log ending at this section header and no error anywhere. The fix for an
    # unbounded hang introduced a silent exit in its place.
    ACCT=$(timeout 20 sacctmgr -n -P show assoc where user="$USER" format=Account 2>/dev/null | head -1 || true)
    QUEUED=$(timeout 20 squeue -h -u "$USER" 2>/dev/null | wc -l | tr -d ' ' || true)
    if [[ -z "$QUEUED" || "$QUEUED" == "0" ]] && ! timeout 20 squeue -h -u "$USER" >/dev/null 2>&1; then
        echo "  squeue did not answer within 20s -- scheduler is loaded." >&2
        echo "  SKIPPING the headroom check. The account GrpTRES cap bounds this" >&2
        echo "  regardless of what we submit; SLURM queues the excess." >&2
        QUEUED=""
    fi
    WANT=$(( ${#KEEP_IDX[@]} * JOBS + ${#KEEP_IDX[@]} ))
    echo "  account=${ACCT:-?}  currently queued/running (this user): ${QUEUED:-unknown}"
    echo "  this wave would add at most $WANT (${#KEEP_IDX[@]} controllers x $JOBS + controllers)"
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
echo "-- pre-flight passed: ${#KEEP_IDX[@]} runs to submit --"
if [[ ${#KEEP_IDX[@]} -eq 0 ]]; then echo "Nothing to submit; every run is already in progress."; exit 0; fi

if [[ "$DRY" == "1" ]]; then
    for i in "${KEEP_IDX[@]}"; do
        printf "  DRY %-4s seed %-4s deep_history %-9s %s\n        WD=%s\n        EXTRA=%s\n" \
            "${RUN_IDS[$i]}" "${RUN_SEED[$i]}" "${RUN_CB[$i]:--}" \
            "${RUN_CFG[$i]}" "${RUN_WD[$i]}" "${RUN_EXTRA[$i]}"
    done
    echo
    echo "DRY=1: nothing submitted."
    exit 0
fi

echo
for i in "${KEEP_IDX[@]}"; do
    # sbatch --output/--error cannot create the directory. The controller does
    # mkdir -p "$WD" itself, but that runs AFTER slurmstepd has tried to open the
    # log files -- fine while every workdir already existed from a previous wave,
    # and not fine for the 360 replicate-6-and-up runs whose directories have
    # never been created.
    mkdir -p "${RUN_WD[$i]}"
    JID=$(sbatch --parsable \
        --job-name="${RUN_IDS[$i]}_${ARM}" \
        --output="${RUN_WD[$i]}/controller_%j.out" \
        --error="${RUN_WD[$i]}/controller_%j.err" \
        --export=ALL,REPO="${REPO}",CONFIGFILE="${RUN_CFG[$i]}",SEED="${RUN_SEED[$i]}",WD="${RUN_WD[$i]}",PD="${RUN_PD[$i]}",STAGE1_SRC="${RUN_S1[$i]}",EXTRA_CONFIG="${RUN_EXTRA[$i]}",JOBS="${JOBS}",UNTIL="${UNTIL}",CB_SEED="${RUN_CB[$i]}",SPECIES="${RUN_SPECIES[$i]}" \
        "$REPO/controller_publication.sbatch")
    printf "  %-4s seed %-4s cb %-9s controller %s\n" \
        "${RUN_IDS[$i]}" "${RUN_SEED[$i]}" "${RUN_CB[$i]:--}" "$JID"
done

echo
echo "Submitted ${#KEEP_IDX[@]} controllers for arm $ARM."
echo "Watch:   squeue -u $USER -o '%.10i %.30j %.8T %.10M'"
echo "Publish: ARM=$ARM bash publish_to_data2.sh   (gated; refuses partial runs)"
