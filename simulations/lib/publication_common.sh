# Shared mechanics for the publication simulation set. Sourced, not executed.
#
# The four round-3 launchers each carry their own copy of the category dispatch
# tables, and they have already drifted -- only submit_2Mb_r3_cmaf_fm001.sh knows
# categories K and L, and submit_2Mb_r3_cmaf01_control.sh has no category_extra at
# all, which means a launcher can dispatch a category whose config belongs to a
# DIFFERENT category and write that category's simulation under the wrong name.
# helpers/tests/test_launcher_tables.py asserts they agree where they overlap.
#
# This file exists so the publication launcher does not become a fifth copy. The
# tables live in helpers/publication_categories.tsv and
# helpers/publication_arms.tsv; nothing here hardcodes a letter.
#
# Roots are ${VAR:-default} so the whole thing can be pointed at a sandbox, which
# is how the stage-1 verification pass runs. The round-3 launchers hardcode them.

REPO="${REPO:-/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/publication_repo/simulations}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/n/scratch/users/n/njc12/snakemake/farm_sims_for_publication}"
PUBLISH_ROOT="${PUBLISH_ROOT:-/n/data2/hms/dbmi/sunyaev/lab/nconnally/farm_sims_for_publication}"
CELL="${CELL:-g5t20}"

CATEGORIES_TSV="$REPO/helpers/publication_categories.tsv"
ARMS_TSV="$REPO/helpers/publication_arms.tsv"
DEEP_HISTORIES_TSV="$REPO/helpers/publication_deep_histories.tsv"

# The rescaled end-of-epoch-7 checkpoints the cattle replicates resume from.
# DELIBERATELY NOT stage1_inputs/: if they shared a directory then
# stage1_search_dirs and cattle_baseline_search_dirs would be one path holding
# four handoffs, and search_dirs._seed_collision_check cannot catch a wrong pick
# among them -- its glob rewrites seed_<N> to sd*, which matches no seed_ file.
DEEP_HISTORIES="${DEEP_HISTORIES:-$PUBLISH_ROOT/deep_histories}"

# ---------------------------------------------------------------- table access

# tsv_field <file> <key-column-value> <wanted-column-name>
# Looks a value up by the first column. Returns 1 and prints nothing if absent,
# so callers can test rather than parse an error message.
tsv_field() {
    local file="$1" key="$2" col="$3"
    awk -F'\t' -v key="$key" -v col="$col" '
        /^#/ || /^[ \t]*$/ { next }
        !seen_header { for (i = 1; i <= NF; i++) idx[$i] = i; seen_header = 1; next }
        $1 == key { if (!(col in idx)) exit 3; print $idx[col]; found = 1; exit 0 }
        END { if (!found) exit 1 }
    ' "$file"
}

tsv_keys() {
    awk -F'\t' '/^#/ || /^[ \t]*$/ { next } !h { h = 1; next } { print $1 }' "$1"
}

cat_field()  { tsv_field "$CATEGORIES_TSV" "$1" "$2"; }
arm_field()  { tsv_field "$ARMS_TSV" "$1" "$2"; }
dh_field()   { tsv_field "$DEEP_HISTORIES_TSV" "$1" "$2"; }

# Every replicate the set defines, in table order. The launcher's REPS default
# and write_run_manifests' --reps default both come from here, so a wave and a
# manifest can never silently disagree about how many replicates exist.
replicate_list() { tsv_keys "$DEEP_HISTORIES_TSV" | tr "\n" " "; }

# The deep history replicate N resumes from. Cattle only -- nothing on the human
# path reads it, and a key that is silently ignored is how a run ends up not
# being the run someone thinks they ran.
cattle_baseline_seed_of() {
    local cb
    if ! cb=$(dh_field "$1" cattle_baseline_seed); then
        echo "ERROR: replicate $1 has no row in publication_deep_histories.tsv." >&2
        echo "       Without one the run would fall back to the config YAML's" >&2
        echo "       20250303 and the independent deep histories would silently" >&2
        echo "       collapse back to one." >&2
        return 1
    fi
    echo "$cb"
}
deep_history_handoff_of() { dh_field "$1" handoff_file; }

# The published directory name for a category letter, and its inverse.
category_dir_of() { cat_field "$1" dir_name; }
species_of()      { cat_field "$1" species; }
pipeline_of()     { cat_field "$1" pipeline; }
config_of()       { local c; c=$(cat_field "$1" config) || return 1
                   echo "${c//\$\{CELL\}/$CELL}"; }
category_extra()  { local e; e=$(cat_field "$1" extra_config) || return 1
                   [[ "$e" == "-" ]] && echo "" || echo "$e"; }

# Which category's stage-1 tree sequence this one uses. K -> A, L -> E; every
# other category is its own donor.
stage1_donor_of() {
    local d; d=$(cat_field "$1" stage1_donor) || return 1
    [[ "$d" == "self" ]] && echo "$1" || echo "$d"
}

# seed = <donor's prefix><replicate>. Read off the DONOR, not the category, so a
# paired category cannot drift out of its donor's seed band -- stage-1 lookup is
# seed-strict, and a mismatch does not fail, it starts a fresh genetic simulation.
seed_of() {
    local letter="${1:0:1}" rep="${1:1}" donor prefix donor_prefix
    donor=$(stage1_donor_of "$letter") || return 1
    prefix=$(cat_field "$letter" seed_prefix) || return 1
    donor_prefix=$(cat_field "$donor" seed_prefix) || return 1
    if [[ "$prefix" != "$donor_prefix" ]]; then
        echo "ERROR: category $letter adopts $donor's stage 1 but sits in seed band ${prefix}N against ${donor_prefix}N" >&2
        return 1
    fi
    echo "${prefix}${rep}"
}

# ------------------------------------------------------------------ pre-flight

# require_marker <pattern> <file> -- the deployed repo must contain <pattern>.
#
# Snakemake accepts UNKNOWN --config keys silently: the Snakefile validates
# REQUIRED_KEYS and rejects nothing else. So a deployment that predates the
# require_gtex_partner work would take `require_gtex_partner=False`, ignore it,
# and run the PAIRED arm under the unpaired arm's name -- and since the two arms
# differ in nothing else, no output would look wrong. These markers are the only
# thing standing between that and a published figure.
require_marker() {
    local pattern="$1" file="$2"
    if ! grep -q -- "$pattern" "$REPO/$file" 2>/dev/null; then
        echo "ERROR: $REPO/$file does not contain '$pattern'." >&2
        echo "       The deployed code predates a change this run depends on." >&2
        return 1
    fi
}

preflight_repo() {
    local rc=0
    require_marker "require_gtex_partner" Snakefile              || rc=1
    require_marker "require_gtex_partner_segment" helpers/paths.py || rc=1
    require_marker "require_gtex_partner" helpers/params_record.py || rc=1
    require_marker "partner_flag" rules/common.smk               || rc=1
    require_marker "str2bool_or_auto" create_gwas_files_and_phenotypes.py || rc=1
    # params_strict has to actually gate. It used to raise ValueError into a
    # blanket `except Exception` in onstart, so the flag could be set, the
    # disagreement found, and the run carry on regardless.
    require_marker "except ValueError:" Snakefile                || rc=1
    require_marker "_fm_thins_candidates" rules/common.smk       || rc=1
    # The cattle deep history must be visible in the F/G stage-1 filename. Without
    # the segment, two histories at one stage1_seed write the same .full.ts, the
    # same .m4_marks.tsv and the same stage2_run_tag, and cattle_baseline_seed is
    # not in params_record.STAGE2_KEYS either -- so nothing anywhere separates them.
    require_marker "cattle_baseline_segment" helpers/paths.py    || rc=1
    # And it must reach the command line, or every replicate silently resumes from
    # the config YAML's single 20250303 history.
    require_marker "CB_SEED" controller_publication.sbatch       || rc=1
    for f in helpers/publication_categories.tsv helpers/publication_arms.tsv \
             helpers/publication_deep_histories.tsv \
             controller_publication.sbatch helpers/write_run_manifests.py; do
        [[ -f "$REPO/$f" ]] || { echo "ERROR: missing $REPO/$f" >&2; rc=1; }
    done
    # The tag is the publication's code provenance. The O2 checkout used to be an
    # rsync copy rather than a git repo, which is why every existing run on disk
    # records `git_commit: null`.
    if ! git -C "$REPO" rev-parse --git-dir >/dev/null 2>&1; then
        echo "ERROR: $REPO is not a git repository, so params_record._git() returns" >&2
        echo "       None and every run records git_commit: null -- exactly the gap" >&2
        echo "       this freeze exists to close. Deploy with git clone, not rsync." >&2
        rc=1
    else
        local dirty; dirty=$(git -C "$REPO" status --porcelain)
        if [[ -n "$dirty" ]]; then
            echo "ERROR: $REPO has uncommitted changes; every run would record" >&2
            echo "       git_dirty: true. Commit or reset before launching." >&2
            git -C "$REPO" status --short >&2
            rc=1
        fi
        echo "  repo $REPO at $(git -C "$REPO" rev-parse --short HEAD)$(git -C "$REPO" describe --exact-match --tags 2>/dev/null | sed 's/^/ tag /')"
    fi
    return $rc
}

# Refuse to reuse a stage-2 directory across arms, and refuse to adopt stages 3-5
# from another workdir. Both are unsafe HERE specifically: K adopts A's tree
# sequence at A's seed, so in the two power arms A and K have byte-identical
# stage2_run_tags and are separated only by their workdir. previous_workdirs
# resolves <prev>/<stage>/<run_tag>, so pointed anywhere near an arm root it would
# cross-adopt A's stage 3/4/5 into K's workdir without a word.
preflight_adoption() {
    local cfg="$1" rc=0
    if grep -Eq '^\s*previous_workdirs\s*:\s*[^[]|^\s*previous_workdirs\s*:\s*\[[^]]' "$REPO/$cfg"; then
        echo "ERROR: $cfg sets previous_workdirs. In this run A and K share a" >&2
        echo "       stage2_run_tag in the power arms, so subtree adoption would" >&2
        echo "       silently copy one category's fine-mapping into the other's." >&2
        rc=1
    fi
    if grep -Eq '^\s*stage2_search_dirs\s*:\s*\[[^]]' "$REPO/$cfg"; then
        echo "ERROR: $cfg sets a non-empty stage2_search_dirs. Stage-2 adoption" >&2
        echo "       keys on the run tag, which A and K share; and the arms differ" >&2
        echo "       only in require_gtex_partner." >&2
        rc=1
    fi
    return $rc
}
