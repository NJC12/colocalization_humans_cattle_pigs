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

# Which round-3 roots may supply the reference. Restricted on purpose: an
# unrestricted search finds simulations_round_2_2Mb/A1/stage1/hts_11.ts, a
# 25 MB file from a different round at different parameters, and comparing
# against it would report DIFFERS for a run that is perfectly fine. The first bug
# this script had.
REF_ROOTS="${REF_ROOTS:-/n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb /n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb_g5t20_cmaf_0.001 /n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb_g5t20_cmaf_0.001_fm_0.001}"

# The expected stage-1 filename, from paths.py -- the same builder the pipeline
# uses. Globbing the fresh directory instead picked up `.ep7.ts`, an INTERMEDIATE
# checkpoint, rather than the `.full.ts` stage 2 actually reads. The second bug
# this script had.
expected_stage1_file() {
    local letter="${1:0:1}" rep="${1:1}" cfg seed
    cfg=$(config_of "$letter") || return 1
    seed=$(seed_of "$1") || return 1
    "$PYTHON" - "$REPO" "$REPO/$cfg" "$seed" <<'PY'
import sys, yaml
sys.path.insert(0, sys.argv[1])
from helpers import paths
cfg = yaml.safe_load(open(sys.argv[2]))
cfg["stage1_seed"] = int(sys.argv[3])
print(paths.stage1_full_filename(cfg))
PY
}

find_reference() {
    local fname="$1" root
    for root in $REF_ROOTS; do
        local hit
        hit=$(find "$root" -maxdepth 3 -path "*/stage1/$fname" -print -quit 2>/dev/null)
        [[ -n "$hit" ]] && { echo "$hit"; return 0; }
    done
    return 1
}

if [[ "${1:-}" == "--check" ]]; then
    printf 'id\tfile\tstatus\tdetail\n' > "$REPORT.tmp"
    ALLOK=1
    for ID in $IDS; do
        WD="$SANDBOX/$ID"
        FN=$(expected_stage1_file "$ID") || { ALLOK=0; continue; }
        FRESH="$WD/stage1/$FN"
        if [[ ! -f "$FRESH" ]]; then
            printf '%s\t%s\tNO_OUTPUT\tstage 1 did not produce it in %s\n' \
                "$ID" "$FN" "$WD" >> "$REPORT.tmp"
            ALLOK=0; continue
        fi
        if ! REF=$(find_reference "$FN"); then
            printf '%s\t%s\tNO_REFERENCE\tnot found under: %s\n' \
                "$ID" "$FN" "$REF_ROOTS" >> "$REPORT.tmp"
            ALLOK=0; continue
        fi
        # Content, not bytes. A .trees file embeds a provenance record -- an ISO
        # timestamp, software versions, the command line -- so identical code from
        # identical seeds differs in BYTES while being the same simulation.
        # Measured: a re-run E1 came out 24 bytes larger, entirely in that record.
        VERDICT=$("$PYTHON" "$REPO/helpers/compare_tree_sequences.py" "$FRESH" "$REF" 2>&1 | tail -1)
        STATUS=$(cut -f1 <<<"$VERDICT"); DETAIL=$(cut -f4 <<<"$VERDICT")
        printf '%s\t%s\t%s\t%s (%s)\n' "$ID" "$FN" "${STATUS:-ERROR}" \
            "${DETAIL:-$VERDICT}" "$REF" >> "$REPORT.tmp"
        # COSMETIC counts as reproduced: same simulation, spatial-coordinate noise.
        [[ "$STATUS" == "SAME" || "$STATUS" == "COSMETIC" ]] || ALLOK=0

        # cattle_sel also writes the m4 marks table, which stage 2 reads. A plain
        # TSV, so there a byte comparison IS the content comparison.
        MARKS="${FRESH%.full.ts}.m4_marks.tsv"
        if [[ -f "$MARKS" ]]; then
            MN=$(basename "$MARKS")
            if MR=$(find_reference "$MN") && cmp -s "$MARKS" "$MR"; then
                printf '%s\t%s\tSAME\tbyte-identical (%s)\n' "$ID" "$MN" "$MR" >> "$REPORT.tmp"
            else
                printf '%s\t%s\t%s\t%s\n' "$ID" "$MN" \
                    "$([[ -n "${MR:-}" ]] && echo DIFFERS || echo NO_REFERENCE)" "${MR:-none}" >> "$REPORT.tmp"
                ALLOK=0
            fi
        fi
    done
    mv "$REPORT.tmp" "$REPORT"
    column -t -s$'\t' "$REPORT"
    echo
    if (( ALLOK )); then
        echo "ALL SAME -- the existing tree sequences are, genetically, what this tag"
        echo "produces. Report: $REPORT"
    else
        echo "NOT all the same. Read the detail column before deciding: a DIFFERENT" >&2
        echo "verdict names which tables moved, and sites+mutations alone is a" >&2
        echo "different problem from nodes+edges. Report: $REPORT" >&2
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
