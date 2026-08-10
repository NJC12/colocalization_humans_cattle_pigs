#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#
# Replay one replicate's stage 2 into a throwaway workdir and prove the outputs
# are byte-identical to the ones already on disk.
#
# WHY THIS EXISTS. Stage 2 is the only stage whose outputs are not reproducible by
# inspection: the trait draw, the redistribution map and the phenotypes all come
# out of seeded RNG streams. So a change to create_gwas_files_and_phenotypes.py
# that is *supposed* to be behaviour-preserving cannot be verified by reading it.
# Worse, nothing downstream would notice if it were not: stage2_run_tag carries
# only the stage-1 stem, .cmaf_* and .psamp_*, and the provenance guard in the
# Snakefile compares PARAMETER VALUES, not code -- git_commit is recorded in the
# sidecar but never compared. A behaviour change at fixed parameters is therefore
# invisible to every guard in the pipeline, and any arm adopting that stage 2
# would silently inherit it. This script is the missing check.
#
# It re-runs stage 2 ONLY (--until stage2_split_pheno), reusing the replicate's
# own stage-1 tree sequence, with exactly the invocation recorded in the
# replicate's params sidecar -- so the only difference between the two runs is the
# code.
#
# Env:
#   SRC   replicate workdir to replay. Required. e.g.
#         /n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb_x30_psamp_8000_fm_0.001/A1
#   DEST  workdir for the replay. Default $(dirname SRC)/.replay/$(basename SRC).
#
# Usage:
#   sbatch --export=ALL,SRC=/n/scratch/.../A1 \
#          --output=/n/scratch/.../replay_A1.log --error=/n/scratch/.../replay_A1.log \
#          replay_stage2_check.sbatch

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SNAKEMAKE=/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake
PYTHON=/home/njc12/miniconda3/envs/coloc_sims/bin/python

: "${SRC:?set SRC to the replicate workdir to replay}"
SRC="${SRC%/}"
DEST="${DEST:-$(dirname "$SRC")/.replay/$(basename "$SRC")}"

PARAMS=$(ls "$SRC"/params/*.params.txt 2>/dev/null | head -1)
[[ -n "$PARAMS" ]] || { echo "ERROR: no params sidecar under $SRC/params/" >&2; exit 1; }

# Recover the original invocation. Only workdir and publishdir are replaced;
# EVERY other --config token is replayed verbatim, which is the whole point -- a
# hand-retyped flag list would be testing the retyping.
#
# stage1_search_dirs is KEPT and merely extended with $SRC/stage1. It must not be
# replaced: an arm that adopted its stage 1 from another root has no stage1/ dir
# of its own (the psamp arms adopt from simulations_round_3_2Mb/), and pointing
# the search somewhere empty does not fail -- the stage-1 rule quietly REGENERATES
# the tree sequence. That is not byte-reproducible even at the same seed (measured
# on A1: 13009284 vs 13042540 bytes, different md5), so every stage-2 output would
# differ for reasons having nothing to do with the code under test.
eval "$("$PYTHON" - "$PARAMS" "$SRC" <<'PY'
import ast, shlex, sys, yaml
d = yaml.safe_load(open(sys.argv[1]))
src = sys.argv[2].rstrip('/')
argv = d['_meta']['argv']
cfg = argv[argv.index('--configfile') + 1]

dirs, tokens = [], []
for t in argv[argv.index('--config') + 1:]:
    if t.startswith('-') or t.startswith(('workdir=', 'publishdir=')):
        continue
    if t.startswith('stage1_search_dirs='):
        dirs = list(ast.literal_eval(t.split('=', 1)[1]))
        continue
    tokens.append(t)
for extra in (f'{src}/stage1',):
    if extra not in dirs:
        dirs.append(extra)

print(f'CFG={shlex.quote(cfg)}')
print(f'TOKENS={shlex.quote(" ".join(tokens))}')
print("STAGE1_CFG=" + shlex.quote('stage1_search_dirs=' + repr(dirs)))
print(f'STAGE1_DIRS={shlex.quote(" ".join(dirs))}')
print(f"STAGE1_FILE={shlex.quote(d['_derived']['stage1_file'])}")
print(f"OLD_INNER={shlex.quote(d['_derived']['stage2_inner'])}")
PY
)"

echo "SRC        $SRC"
echo "DEST       $DEST"
echo "configfile $CFG"
echo "config     $TOKENS"
echo "stage1     $STAGE1_FILE in [$STAGE1_DIRS]"
echo "old inner  $OLD_INNER"
[[ -d "$OLD_INNER" ]] || { echo "ERROR: $OLD_INNER is gone; nothing to compare against." >&2; exit 1; }

# Pin the tree sequence BEFORE running. A replay that regenerates stage 1 compares
# two different simulations and its verdict is meaningless, so refuse rather than
# report a difference the code change did not cause.
STAGE1_SRC=""
for d in $STAGE1_DIRS; do
    [[ -f "$d/$STAGE1_FILE" ]] && { STAGE1_SRC="$d/$STAGE1_FILE"; break; }
done
[[ -n "$STAGE1_SRC" ]] || {
    echo "ERROR: $STAGE1_FILE is in none of the search dirs:" >&2
    for d in $STAGE1_DIRS; do echo "       $d" >&2; done
    echo "       Refusing to replay -- stage 1 would be regenerated and the" >&2
    echo "       comparison would be against a different simulation." >&2; exit 1; }
STAGE1_MD5=$(md5sum "$STAGE1_SRC" | awk '{print $1}')
echo "stage1 src $STAGE1_SRC  ($STAGE1_MD5)"

mkdir -p "$DEST"
cd "$DEST"

"$SNAKEMAKE" --snakefile "$REPO/Snakefile" --configfile "$CFG" --unlock \
    --config workdir="$DEST" publishdir="$DEST/publish" "$STAGE1_CFG" $TOKENS || true

"$SNAKEMAKE" --snakefile "$REPO/Snakefile" --configfile "$CFG" \
    --profile "$REPO/profiles/o2" --rerun-triggers mtime \
    --until stage2_split_pheno \
    --config workdir="$DEST" publishdir="$DEST/publish" "$STAGE1_CFG" $TOKENS

# Belt and braces: even with the search dirs right, confirm the tree sequence the
# replay actually consumed is the one the original run did.
if [[ -f "$DEST/$STAGE1_FILE" || -f "$DEST/stage1/$STAGE1_FILE" ]]; then
    USED=$([[ -f "$DEST/stage1/$STAGE1_FILE" ]] && echo "$DEST/stage1/$STAGE1_FILE" \
                                                || echo "$DEST/$STAGE1_FILE")
    USED_MD5=$(md5sum "$USED" | awk '{print $1}')
    if [[ "$USED_MD5" != "$STAGE1_MD5" ]]; then
        echo "ERROR: the replay consumed a DIFFERENT stage-1 tree sequence" >&2
        echo "       expected $STAGE1_MD5 ($STAGE1_SRC)" >&2
        echo "       got      $USED_MD5 ($USED)" >&2
        echo "       Stage 1 was regenerated; the comparison would be meaningless." >&2
        exit 1
    fi
    echo "stage1 ok  replay used the same tree sequence ($USED_MD5)"
fi

NEW_INNER=$(ls -d "$DEST"/stage2/*/*/ 2>/dev/null | head -1)
[[ -n "$NEW_INNER" ]] || { echo "ERROR: replay produced no stage-2 inner dir." >&2; exit 1; }
NEW_INNER="${NEW_INNER%/}"
echo
echo "new inner  $NEW_INNER"
if [[ "$(basename "$NEW_INNER")" != "$(basename "$OLD_INNER")" ]]; then
    echo "WARNING: inner dir name changed ($(basename "$OLD_INNER") -> $(basename "$NEW_INNER"))."
    echo "         That means a path-determining parameter moved, which is a bigger"
    echo "         change than a behaviour-preserving fix. Comparing anyway."
fi

# stage2_params.txt records host, timestamp and argv, so it differs by design.
# Compressed outputs carry an mtime in the gzip header and are not written by
# stage 2 itself, so they are skipped too -- everything stage 2 actually produces
# is plain text or sbams.
echo
echo "=== byte comparison ==="
SAME=0; DIFF=0; ONLY=0
for f in "$OLD_INNER"/*; do
    b=$(basename "$f")
    case "$b" in stage2_params.txt|*.gz|.stage2.done) continue ;; esac
    if [[ ! -e "$NEW_INNER/$b" ]]; then
        printf '  MISSING   %s\n' "$b"; ONLY=$((ONLY + 1)); continue
    fi
    if cmp -s "$f" "$NEW_INNER/$b"; then
        SAME=$((SAME + 1))
    else
        printf '  DIFFERS   %s  (%s vs %s bytes)\n' "$b" \
            "$(stat -c %s "$f")" "$(stat -c %s "$NEW_INNER/$b")"
        DIFF=$((DIFF + 1))
    fi
done
for f in "$NEW_INNER"/*; do
    b=$(basename "$f")
    case "$b" in stage2_params.txt|*.gz|.stage2.done) continue ;; esac
    [[ -e "$OLD_INNER/$b" ]] || { printf '  NEW ONLY  %s\n' "$b"; ONLY=$((ONLY + 1)); }
done

echo
echo "identical: $SAME   differing: $DIFF   present in only one: $ONLY"
if [[ "$DIFF" -eq 0 && "$ONLY" -eq 0 ]]; then
    echo "REPLAY CLEAN -- the change is behaviour-preserving for $(basename "$SRC")."
else
    echo "REPLAY DIRTY -- the change moved stage-2 output. Do not submit anything."
    exit 1
fi
