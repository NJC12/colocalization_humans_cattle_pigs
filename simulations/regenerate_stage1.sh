#!/bin/bash
#
# Regenerate every replicate's stage-1 tree sequence under the frozen tag, then
# collect them into the publication's stage1_inputs/.
#
#   bash regenerate_stage1.sh              # submit all 20
#   IDS="A1 A2" bash regenerate_stage1.sh
#   bash regenerate_stage1.sh --collect    # once they finish
#   bash regenerate_stage1.sh --status
#
# Twenty, not thirty: K and L adopt A's and E's tree sequences at A's and E's
# seeds, so they never run stage 1.
#
# WHY REGENERATE RATHER THAN ADOPT. Two RNG calls were unseeded until the freeze
# commit -- stdpopsim's SLiM engine (human) and the selection-coefficient
# permutation (cattle deep history) -- so nothing on disk could be regenerated
# from the seed it was named for. Everything here is built by the tagged code
# from its seed, and can be rebuilt from them.
#
# THE OVERRIDE THAT MATTERS. Every cattle config's `cattle_baseline_search_dirs`
# still points at the round-3 deep history, whose ancestor was built with the
# UNSEEDED permutation. Left alone, E/F/G would resume from that checkpoint and
# the regeneration would be pointless -- silently, because resuming from a valid
# checkpoint looks exactly like success. This script points them at
# stage1_inputs/ instead, and refuses to run if the new handoff is not there.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${REPO:-$HERE}"
source "$REPO/lib/publication_common.sh"

SNAKEMAKE="${SNAKEMAKE:-/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake}"
PYTHON="${PYTHON:-/home/njc12/miniconda3/envs/coloc_sims/bin/python}"
BUILD="${BUILD:-$SCRATCH_ROOT/stage1_build}"
STAGE1_INPUTS="${STAGE1_INPUTS:-$PUBLISH_ROOT/stage1_inputs}"
# The self-donor categories x every replicate the deep-history table defines.
# K and L adopt A's and E's tree sequences, so they never build one.
_default_ids() {
    local l r
    for l in A E F G; do
        for r in $(replicate_list); do printf "%s%s " "$l" "$r"; done
    done
}
IDS="${IDS:-$(_default_ids)}"

expected_stage1_file() {
    local letter="${1:0:1}" cfg seed
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

if [[ "${1:-}" == "--status" || "${1:-}" == "--collect" ]]; then
    MISSING=0
    declare -a READY=()

    # Two staleness guards, both learned the hard way.
    #
    # A cancelled run leaves its outputs behind. When E/F/G were cancelled and
    # relaunched against a corrected cattle handoff, the PREVIOUS round's
    # .full.ts files were still on disk, so --status reported "20 of 20 ready"
    # while the new runs had not produced anything yet. Snakemake removed them
    # itself a minute later, as incomplete -- but in that window --collect would
    # have published tree sequences built from the WRONG checkpoint, and nothing
    # downstream would have noticed: they are valid tree sequences at the right
    # seeds with the right filenames.
    #
    # 1. Refuse while any controller is still running. A file that exists while
    #    its producer is still going is either stale or half-written.
    ACTIVE=$(squeue -h -u "$USER" -o "%j" 2>/dev/null | grep -c "^regen_stage1_" || true)
    if [[ "${ACTIVE:-0}" -gt 0 && "${1:-}" == "--collect" ]]; then
        echo "ERROR: $ACTIVE regen_stage1 controllers are still active." >&2
        echo "       A stage-1 file present while its producer is still running is" >&2
        echo "       either a leftover from a previous attempt or half-written." >&2
        exit 1
    fi
    # 2. No CATTLE stage-1 file may predate the handoff, because it resumes from
    #    it and an older file was therefore built from a different one. Human is
    #    exempt: it does not touch the handoff, so an A tree sequence built before
    #    the handoff was installed is perfectly current.
    # Per-replicate: with four deep histories a single reference would compare
    # every replicate against block 1's handoff and call blocks 2-4 stale.
    for ID in $IDS; do
        FN=$(expected_stage1_file "$ID") || { MISSING=1; continue; }
        SRC="$BUILD/$ID/stage1/$FN"
        HANDOFF_REF=""
        if [[ "$(species_of "${ID:0:1}")" == "cattle" ]]; then
            HANDOFF_REF="$DEEP_HISTORIES/$(deep_history_handoff_of "${ID:1}")"
        fi
        if [[ ! -f "$SRC" ]]; then
            printf "  pending  %-3s %s\n" "$ID" "$FN"
            MISSING=1
        elif [[ -n "$HANDOFF_REF" && -f "$HANDOFF_REF" && "$SRC" -ot "$HANDOFF_REF" ]]; then
            printf "  STALE    %-3s %s  (older than the handoff it should derive from)\n" "$ID" "$FN"
            MISSING=1
        else
            READY+=("$ID:$SRC")
            printf "  ready    %-3s %s\n" "$ID" "$FN"
        fi
    done
    echo
    if [[ "${1:-}" == "--status" ]]; then
        echo "$(( ${#READY[@]} )) of $(wc -w <<<"$IDS") ready"
        exit 0
    fi
    if (( MISSING )); then
        echo "Not collecting: some replicates have not finished. An incomplete" >&2
        echo "stage1_inputs/ does not fail loudly -- Snakemake decides it must BUILD" >&2
        echo "the missing tree sequence and starts a genetic simulation inside a" >&2
        echo "controller sized for symlinks." >&2
        exit 1
    fi
    mkdir -p "$STAGE1_INPUTS"
    for entry in "${READY[@]}"; do
        SRC="${entry#*:}"
        cp -pu "$SRC" "$STAGE1_INPUTS/"
        # cattle_sel also writes the m4 marks table, which stage 2 reads.
        MARKS="${SRC%.full.ts}.m4_marks.tsv"
        [[ -f "$MARKS" ]] && cp -pu "$MARKS" "$STAGE1_INPUTS/"
    done
    echo "collected $(( ${#READY[@]} )) tree sequences into $STAGE1_INPUTS"
    ls -1 "$STAGE1_INPUTS" | wc -l | xargs echo "  files now in stage1_inputs:"
    exit 0
fi

# The new handoff must be in place first, or the cattle categories resume from
# the archived one and the regeneration silently achieves nothing.
# Every handoff this wave needs, one per deep history it touches. A missing one
# does not fail loudly downstream: Snakemake would decide it must BUILD the deep
# history, which is a multi-hour job inside a stage-1 sbatch.
MISSING_HANDOFF=0
for ID in $IDS; do
    [[ "$(species_of "${ID:0:1}")" == "cattle" ]] || continue
    H="$DEEP_HISTORIES/$(deep_history_handoff_of "${ID:1}")" || { MISSING_HANDOFF=1; continue; }
    if [[ ! -f "$H" ]]; then
        echo "ERROR: missing handoff for $ID: $H" >&2
        MISSING_HANDOFF=1
    fi
done
if (( MISSING_HANDOFF )); then
    echo "       Run rebuild_cattle_deep_history.sh (and --rescale) for the missing" >&2
    echo "       seeds, then collect them into $DEEP_HISTORIES. Without the right" >&2
    echo "       handoff a replicate resumes from a DIFFERENT ancestry than the one" >&2
    echo "       publication_deep_histories.tsv assigns it, and the filename would" >&2
    echo "       still say otherwise." >&2
    exit 1
fi

echo "repo:          $REPO ($(git -C "$REPO" describe --exact-match --tags 2>/dev/null || git -C "$REPO" rev-parse --short HEAD))"
echo "build dir:     $BUILD"
echo "deep histories: $DEEP_HISTORIES ($(awk -F'\t' 'NR>1{print $2}' "$DEEP_HISTORIES_TSV" | sort -u | wc -l | tr -d ' ') distinct)"
echo

for ID in $IDS; do
    L="${ID:0:1}"
    CFG=$(config_of "$L") || exit 1
    SEED=$(seed_of "$ID") || exit 1
    WD="$BUILD/$ID"
    mkdir -p "$WD"
    # stage1_search_dirs=[] so it BUILDS. cattle_baseline_search_dirs only for the
    # cattle pipelines -- the Snakefile is strict about keys reaching the wrong
    # pipeline, and the human configs have no such key.
    EXTRA="stage1_search_dirs=[]"
    CB=""
    if [[ "$(species_of "$L")" == "cattle" ]]; then
        # The deep history this replicate resumes from. Without the explicit seed
        # every replicate would take the config YAML's 20250303 and the four
        # ancestries would collapse into one -- while the F/G filename, which now
        # carries _cb{seed}, would claim otherwise.
        CB=$(cattle_baseline_seed_of "${ID:1}") || exit 1
        EXTRA+=" cattle_baseline_seed=$CB cattle_baseline_search_dirs=['$DEEP_HISTORIES']"
    fi
    JID=$(sbatch --parsable \
        --job-name="regen_stage1_$ID" \
        --partition=short --time=8:00:00 --mem=16G --cpus-per-task=4 \
        --output="$WD/regen_%j.out" --error="$WD/regen_%j.err" \
        --wrap="cd '$WD' && '$SNAKEMAKE' --snakefile '$REPO/Snakefile' \
            --configfile '$REPO/$CFG' --unlock \
            --config stage1_seed=$SEED stage2_seed=$SEED workdir='$WD' \
                     publishdir='$WD' $EXTRA || true; \
            '$SNAKEMAKE' --snakefile '$REPO/Snakefile' --configfile '$REPO/$CFG' \
            --profile '$REPO/profiles/o2' --rerun-triggers mtime -j 4 --until stage1 \
            --config stage1_seed=$SEED stage2_seed=$SEED workdir='$WD' \
                     publishdir='$WD' $EXTRA")
    printf "  %-4s seed %-4s cb %-9s job %s   %s\n" "$ID" "$SEED" "${CB:--}" "$JID" "$CFG"
done

echo
echo "Then:  bash regenerate_stage1.sh --status"
echo "       bash regenerate_stage1.sh --collect"
