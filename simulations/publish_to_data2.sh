#!/bin/bash
#
# Copy finished publication runs from the scratch workdirs to
# /n/data2/.../farm_sims_for_publication/, one arm at a time, and refuse to
# publish anything incomplete.
#
#   ARM=causal_maf001_paired bash publish_to_data2.sh
#   ARM=... DRY=1 bash publish_to_data2.sh          # report only
#   ARM=... WITH_GLM=1 bash publish_to_data2.sh     # also the full GLM stats
#
# A COPY, not a move: scratch is left untouched, and re-running is idempotent
# (cp -u).
#
# WHY THE GATE IS THE POINT. profiles/o2/config.yaml sets keep-going: True, so a
# failed job does not stop the rest of the DAG. That is the right setting -- one
# flaky trait should not waste a controller -- but it means a replicate can finish
# with a truncated or missing .enloc.sig.out, and downstream that does not read as
# "unfinished". It reads as "this arm colocalized nothing", which is a plausible
# wrong number rather than an obvious gap. Nothing else in the pipeline checks
# completeness, so this script is the only thing standing between a partial run
# and a figure.
#
# Expected counts are derived from RUNS.tsv and ARMS.tsv rather than hardcoded,
# because they depend on the arm (power arms write a causal-power sidecar) and on
# the category (K and L write a synthetic-DFE table).

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${REPO:-$HERE}"
source "$REPO/lib/publication_common.sh"

: "${ARM:?set ARM to one of: $(tsv_keys "$ARMS_TSV" | tr '\n' ' ')}"
DRY="${DRY:-0}"
WITH_GLM="${WITH_GLM:-0}"
MANIFEST="${MANIFEST:-$PUBLISH_ROOT/RUNS.tsv}"
LOG="$PUBLISH_ROOT/PUBLISH_LOG.tsv"

[[ -f "$MANIFEST" ]] || { echo "ERROR: $MANIFEST not found" >&2; exit 1; }

SAMPLING=$(arm_field "$ARM" causal_sampling)

echo "=============================================================="
echo " publishing arm  $ARM"
echo " from            $SCRATCH_ROOT/$ARM/"
echo " to              $PUBLISH_ROOT/$ARM/"
echo " full GLM stats  $([[ "$WITH_GLM" == 1 ]] && echo yes || echo "no (lead SNPs only)")"
echo "=============================================================="

# ------------------------------------------------------------------ the gate

# n_expected <glob> -> how many files must be there
# Returns "<found>/<wanted>" and 0 when they match.
check_count() {
    local wanted="$1"; shift
    local found
    found=$(ls -1 "$@" 2>/dev/null | wc -l | tr -d ' ')
    echo "$found"
    [[ "$found" -eq "$wanted" ]]
}

# complete_p <workdir> <letter> <run_tag>
# Sets GATE_DETAIL. Returns 0 only when every expected output is present.
complete_p() {
    local wd="$1" letter="$2" tag="$3"
    GATE_DETAIL=""
    local s2 fails=0 n

    s2=$(ls -d "$wd"/stage2/"$tag"/gwas_*_gtex_*_maf_* 2>/dev/null | head -1)
    if [[ -z "$s2" ]]; then
        GATE_DETAIL="no stage-2 directory for run tag $tag"
        return 1
    fi

    _need() {  # _need <label> <wanted> <glob...>
        local label="$1" wanted="$2"; shift 2
        n=$(ls -1 "$@" 2>/dev/null | wc -l | tr -d ' ')
        if [[ "$n" -ne "$wanted" ]]; then
            GATE_DETAIL+="${GATE_DETAIL:+; }$label=$n/want$wanted"
            fails=1
        fi
    }

    # Four panels: gwas, gtex, gtex_smaller, gtex_smallest. The smallest is
    # written by stage 2 but not analysed, and is still expected here -- its
    # absence means stage 2 did not finish, not that it was skipped.
    _need vars      4 "$s2"/*_vars_*.tsv
    _need pheno     4 "$s2"/*_pheno.sbams
    # Written for every publication run, because an ABSENT partner table is
    # ambiguous between "intersected, so all-True" and "topped up, written by a
    # version that did not emit it". A scorer must never have to guess.
    _need partners  1 "$s2"/*_trait_partners_*.tsv
    _need s2params  1 "$s2"/stage2_params.txt
    case "$letter" in
        K|L) _need synthdfe 1 "$s2"/*_synthetic_dfe_*.tsv ;;
    esac
    [[ "$SAMPLING" == "power" ]] && _need causalpower 1 "$s2"/*_causal_power_*.tsv

    _need stage3done 3 "$wd"/stage3/"$tag"/*/.stage3.done
    _need enlocsig   2 "$wd"/stage4/"$tag"/*.enloc.sig.out
    _need s4params   1 "$wd"/stage4/"$tag"/*.params.txt
    _need glmleads   3 "$wd"/stage5/"$tag"/*_glm_lead_snps.tsv
    _need stage5done 3 "$wd"/stage5/"$tag"/.*.stage5.done

    # The provenance file must describe THIS run, not one adopted from elsewhere.
    local wdparams
    wdparams=$(ls -1 "$wd"/params/"$tag".params.txt 2>/dev/null | head -1)
    if [[ -z "$wdparams" ]]; then
        GATE_DETAIL+="${GATE_DETAIL:+; }workdir params file missing for $tag"
        fails=1
    else
        # git_commit: null is what every run before the freeze recorded. A
        # publication run that records it means the code was not deployed as a
        # git checkout, which is the one thing this freeze exists to fix.
        if grep -Eq '^\s*git_commit:\s*null\s*$' "$wdparams"; then
            GATE_DETAIL+="${GATE_DETAIL:+; }git_commit is null -- code was not a git checkout"
            fails=1
        fi
        if grep -Eq '^\s*git_dirty:\s*true\s*$' "$wdparams"; then
            GATE_DETAIL+="${GATE_DETAIL:+; }git_dirty is true -- deployed tree had local changes"
            fails=1
        fi
    fi

    # A provenance write that failed leaves the summary reading params_source=
    # inferred, with the MAF floors reverse-engineered from directory names.
    if compgen -G "$wd/controller_*.out" >/dev/null 2>&1; then
        if grep -qs '\[params\] WARNING' "$wd"/controller_*.out; then
            GATE_DETAIL+="${GATE_DETAIL:+; }controller logged a [params] WARNING"
            fails=1
        fi
    fi

    return $fails
}

# --------------------------------------------------------------- the copying

copy_run() {
    local wd="$1" tag="$2" dest="$3"
    local s2; s2=$(ls -d "$wd"/stage2/"$tag"/gwas_*_gtex_*_maf_* | head -1)
    mkdir -p "$dest"
    # -L so an adopted symlink is published as the file it points at: scratch will
    # be purged, and a dangling link is worse than no link.
    local pat
    for pat in "$s2"/*_pheno.sbams "$s2"/*_vars_*.tsv "$s2"/*_trait_partners_*.tsv \
               "$s2"/*_synthetic_dfe_*.tsv "$s2"/*_causal_power_*.tsv \
               "$s2"/stage2_params.txt \
               "$wd"/stage4/"$tag"/*.enloc.sig.out "$wd"/stage4/"$tag"/*.params.txt \
               "$wd"/stage5/"$tag"/*_glm_lead_snps.tsv \
               "$wd"/params/"$tag".params.txt; do
        compgen -G "$pat" >/dev/null 2>&1 && cp -pLu $pat "$dest"/ 2>/dev/null
    done
    # The only record of the pool-adequacy guard's counts and the realised partner
    # counts. A few KB, and the thing you want when a number looks wrong.
    if [[ -f "$wd/logs/stage2_split_pheno.log" ]]; then
        cp -pLu "$wd/logs/stage2_split_pheno.log" "$dest/stage2_split_pheno.log"
    fi
    if [[ "$WITH_GLM" == "1" ]]; then
        mkdir -p "$dest/glm"
        for f in "$wd"/stage5/"$tag"/*.glm.linear; do
            [[ -f "$f" ]] || continue
            local out="$dest/glm/$(basename "$f").gz"
            [[ -f "$out" ]] || gzip -c "$f" > "$out"
        done
    fi
}

# ------------------------------------------------------------------ main loop

mkdir -p "$PUBLISH_ROOT/$ARM"
PUBLISHED=0; SKIPPED=0
declare -a FAILED=()

while IFS=$'\t' read -r m_arm m_dir m_letter m_tag; do
    [[ "$m_arm" == "$ARM" ]] || continue
    WD="$SCRATCH_ROOT/$ARM/$m_dir"
    DEST="$PUBLISH_ROOT/$ARM/$m_dir"
    if [[ ! -d "$WD" ]]; then
        FAILED+=("$m_dir	no workdir at $WD"); SKIPPED=$((SKIPPED+1)); continue
    fi
    if complete_p "$WD" "$m_letter" "$m_tag"; then
        if [[ "$DRY" == "1" ]]; then
            printf "  WOULD PUBLISH  %s\n" "$m_dir"
        else
            copy_run "$WD" "$m_tag" "$DEST"
            printf "  published      %s\n" "$m_dir"
        fi
        PUBLISHED=$((PUBLISHED+1))
    else
        printf "  INCOMPLETE     %-44s %s\n" "$m_dir" "$GATE_DETAIL"
        FAILED+=("$m_dir	$GATE_DETAIL")
        SKIPPED=$((SKIPPED+1))
    fi
done < <(awk -F'\t' 'NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
                     { print $h["arm"] "\t" $h["run_dir"] "\t" $h["letter"] "\t" $h["stage2_run_tag"] }' "$MANIFEST")

# stage1_inputs and the manifests live once at the publication root, not per
# replicate: 20 tree sequences serve all 120 runs (A and K share one, and every
# arm reuses it), so a per-replicate copy would be 8x redundant. RUNS.tsv's
# stage1_file column is the pointer.
if [[ "$DRY" != "1" ]]; then
    for f in ARMS.tsv RUNS.tsv; do
        [[ -f "$REPO/../$f" ]] && cp -pu "$REPO/../$f" "$PUBLISH_ROOT/$f" 2>/dev/null
    done
    {
        echo "tag        $(git -C "$REPO" describe --exact-match --tags 2>/dev/null || echo '(none)')"
        echo "commit     $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
        echo "clean      $([[ -z "$(git -C "$REPO" status --porcelain 2>/dev/null)" ]] && echo yes || echo NO)"
        echo "repo       $REPO"
    } > "$PUBLISH_ROOT/CODE_VERSION"
fi

echo
echo "--------------------------------------------------------------"
echo " arm $ARM: published=$PUBLISHED  incomplete_skipped=$SKIPPED"
if (( SKIPPED )); then
    echo
    echo " NOT PUBLISHED:"
    printf '   %s\n' "${FAILED[@]}" | column -t -s$'\t'
    if [[ "$DRY" != "1" ]]; then
        { printf '%s\t%s\n' "$ARM" "$(date -Iseconds)"
          printf '%s\t%s\n' "$ARM" "${FAILED[@]}" ; } >> "$LOG" 2>/dev/null || true
    fi
    echo
    echo " A skipped replicate is NOT a null result. Fix or explain each one" >&2
    echo " before any summary table is read." >&2
    exit 1
fi
echo " all runs for this arm are complete and published."
