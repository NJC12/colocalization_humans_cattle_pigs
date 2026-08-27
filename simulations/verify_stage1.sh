#!/bin/bash
#
# Re-run one replicate's stage 1 under the frozen tag and byte-compare it against
# the copy already on scratch.
#
#   bash verify_stage1.sh                       # A1 E1 F1 G1
#   IDS="A1" bash verify_stage1.sh
#   bash verify_stage1.sh --check               # compare results of a previous run
#
# WHY. The publication run adopts 20 existing stage-1 tree sequences rather than
# re-simulating them. They were built between 2026-07-27 and 2026-08-04 by code
# that was not committed at the time, and every params record on disk says
# `git_commit: null`, so "which code made these" has no answer from the files
# themselves.
#
# Reading the diff says they should be reproducible: the only stage-1 code touched
# since is the human path (2026-08-24, for the Finnish M/N pair), and along the A
# path it is additive -- sample_dict('EUR', n, 'OutOfAfrica_2T12') returns the old
# {"AFR": 0, "EUR": n}, stage1_ts_name returns the old hts_{seed}.ts, and the
# catalog branch is the old two lines. Cattle stage-1 code was already committed
# and is unchanged.
#
# This turns that reading into a measurement. One replicate per category is enough
# to catch a code or environment change, because a change that alters stage 1 at
# all would alter every replicate of that pipeline.
#
# NOT a claim about seed-reproducibility of the cattle deep history: F, G, E and L
# all resume from one archived farm_selection_*.ep7.ts, and the ancestor of that
# checkpoint was built with an unseeded np.random.permutation
# (farm_create_orig_pop_e2.py). What this verifies is the epochs downstream of the
# checkpoint, which is the part this run actually executes.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${REPO:-$HERE}"
source "$REPO/lib/publication_common.sh"

SANDBOX="${SANDBOX:-/n/scratch/users/n/njc12/pubfreeze_sandbox/stage1_verify}"
SNAKEMAKE="${SNAKEMAKE:-/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake}"
PYTHON="${PYTHON:-/home/njc12/miniconda3/envs/coloc_sims/bin/python}"
IDS="${IDS:-A1 E1 F1 G1}"
REPORT="$SANDBOX/VERIFY_REPORT.tsv"

# Where each replicate's existing tree sequence lives. Resolved by searching the
# round-3 roots rather than hardcoded, so a moved arm does not silently make the
# comparison vacuous (a missing reference must FAIL, not pass).
EXISTING_ROOTS="${EXISTING_ROOTS:-/n/scratch/users/n/njc12/snakemake}"

find_existing() {
    local fname="$1"
    find "$EXISTING_ROOTS" -maxdepth 4 -path "*/stage1/$fname" -print -quit 2>/dev/null
}

if [[ "${1:-}" == "--check" ]]; then
    echo "id	file	status	detail" > "$REPORT.tmp"
    ALLOK=1
    for ID in $IDS; do
        L="${ID:0:1}"; R="${ID:1}"
        WD="$SANDBOX/$ID"
        FRESH=$(find "$WD/stage1" -maxdepth 1 -name '*.ts' 2>/dev/null | head -1)
        if [[ -z "$FRESH" ]]; then
            printf "%s\t-\tNO_OUTPUT\tstage 1 produced nothing in %s\n" "$ID" "$WD" >> "$REPORT.tmp"
            ALLOK=0; continue
        fi
        FN=$(basename "$FRESH")
        REF=$(find_existing "$FN")
        if [[ -z "$REF" ]]; then
            printf "%s\t%s\tNO_REFERENCE\tnothing to compare against under %s\n" "$ID" "$FN" "$EXISTING_ROOTS" >> "$REPORT.tmp"
            ALLOK=0; continue
        fi
        # Content comparison, not `cmp`. A .trees file embeds a provenance
        # record -- timestamp, software versions, resolved command line -- so two
        # runs of identical code from identical seeds differ in BYTES while being
        # the same simulation. Measured: a re-run E1 came out 24 bytes larger than
        # the copy on disk, entirely in that record. `cmp` would have said "the
        # code changed, re-run all 20", which is both costly and wrong.
        VERDICT=$("$PYTHON" "$REPO/helpers/compare_tree_sequences.py" "$FRESH" "$REF" 2>&1 | tail -1)
        STATUS=$(cut -f1 <<<"$VERDICT")
        DETAIL=$(cut -f4 <<<"$VERDICT")
        if [[ "$STATUS" == "SAME" ]]; then
            printf "%s\t%s\tSAME\t%s (%s)\n" "$ID" "$FN" "$DETAIL" "$REF" >> "$REPORT.tmp"
        else
            printf "%s\t%s\t%s\t%s (%s)\n" "$ID" "$FN" "${STATUS:-ERROR}" "${DETAIL:-$VERDICT}" "$REF" >> "$REPORT.tmp"
            ALLOK=0
        fi
        # cattle_sel also emits the m4 marks table, which stage 2 reads.
        for extra in "${FRESH%.full.ts}.m4_marks.tsv"; do
            [[ -f "$extra" ]] || continue
            EN=$(basename "$extra"); ER=$(find_existing "$EN")
            # A plain TSV, so a byte comparison IS the content comparison here.
            if [[ -n "$ER" ]] && cmp -s "$extra" "$ER"; then
                printf "%s\t%s\tIDENTICAL\t%s\n" "$ID" "$EN" "$ER" >> "$REPORT.tmp"
            else
                printf "%s\t%s\t%s\t%s\n" "$ID" "$EN" \
                    "$([[ -n "$ER" ]] && echo DIFFERS || echo NO_REFERENCE)" "${ER:-none}" >> "$REPORT.tmp"
                ALLOK=0
            fi
        done
    done
    mv "$REPORT.tmp" "$REPORT"
    column -t -s$'\t' "$REPORT"
    echo
    if (( ALLOK )); then
        echo "ALL SAME -- the existing tree sequences are, genetically, what this tag"
        echo "produces. Reuse them."
        echo "Report: $REPORT"
    else
        echo "NOT all the same. Do NOT reuse: re-run stage 1 for all 20 replicates" >&2
        echo "under the tag, or explain the difference first. Report: $REPORT" >&2
        exit 1
    fi
    exit 0
fi

mkdir -p "$SANDBOX"
echo "sandbox: $SANDBOX"
echo "repo:    $REPO  ($(git -C "$REPO" describe --exact-match --tags 2>/dev/null || git -C "$REPO" rev-parse --short HEAD))"
echo

for ID in $IDS; do
    L="${ID:0:1}"
    CFG=$(config_of "$L") || exit 1
    SEED=$(seed_of "$ID") || exit 1
    WD="$SANDBOX/$ID"
    mkdir -p "$WD"
    # stage1_search_dirs=[] is the whole point: BUILD it, do not adopt it.
    # Trait counts and floors are irrelevant to stage 1 and are left at the
    # config's own values so nothing about the comparison depends on the arm.
    # --unlock first: a re-launch (or a launch interrupted partway through this
    # loop, which is easy over a slow ssh) otherwise dies with LockException
    # against its own predecessor's stale lock.
    JID=$(sbatch --parsable \
        --job-name="verify_stage1_$ID" \
        --partition=short --time=8:00:00 --mem=16G --cpus-per-task=4 \
        --output="$WD/verify_%j.out" --error="$WD/verify_%j.err" \
        --wrap="cd '$WD' && '$SNAKEMAKE' --snakefile '$REPO/Snakefile' \
            --configfile '$REPO/$CFG' --unlock \
            --config stage1_seed=$SEED stage2_seed=$SEED workdir='$WD' \
                     publishdir='$WD' stage1_search_dirs=[] || true; \
            '$SNAKEMAKE' --snakefile '$REPO/Snakefile' \
            --configfile '$REPO/$CFG' --profile '$REPO/profiles/o2' \
            --rerun-triggers mtime -j 4 --until stage1 \
            --config stage1_seed=$SEED stage2_seed=$SEED workdir='$WD' \
                     publishdir='$WD' stage1_search_dirs=[]")
    printf "  %-3s seed %-3s job %s   %s\n" "$ID" "$SEED" "$JID" "$CFG"
done

echo
echo "When they finish:  bash verify_stage1.sh --check"
