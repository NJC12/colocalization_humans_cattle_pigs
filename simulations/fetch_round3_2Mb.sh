#!/bin/bash
# Pull round-3 2 Mb replicate results from O2 into simulation_data/, fetching only
# what is COMPLETE upstream and MISSING locally.
#
# WHY GATE ON COMPLETENESS
#
# A replicate is only usable once stage 4 has written its `.enloc.sig.out`: the R
# notebook keys colocalization off that file, and a replicate directory holding
# stage-2 and stage-5 output but no sig.out reads as "this arm colocalized
# nothing" rather than "this arm is still running". Copying a half-finished
# replicate down therefore does not just waste a transfer, it produces a
# plausible wrong number. So a replicate is fetched only when all four pieces are
# present upstream, and a partially-finished one is REPORTED and skipped.
#
# Re-run it as jobs land; rsync skips what is already here, so it converges.
#
#   bash simulations/fetch_round3_2Mb.sh              # check + fetch
#   DRY=1 bash simulations/fetch_round3_2Mb.sh        # report only, transfer nothing
#   ARMS=cmaf0_fm01_g5t20_2Mb bash simulations/fetch_round3_2Mb.sh   # one arm
#
# HOST defaults to o2.hms.harvard.edu because ~/.ssh/config already keeps a
# ControlMaster open for it (ControlPersist 8h), so repeated runs need no further
# Duo push. For a bulk first pull, HOST=transfer.rc.hms.harvard.edu is the
# polite choice -- it opens its own master and costs one push.

set -uo pipefail

HOST="${HOST:-o2.hms.harvard.edu}"
SCRATCH=/n/scratch/users/n/njc12/snakemake
LOCAL="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/simulation_data"
DRY="${DRY:-}"

# local arm directory  <-  O2 root suffix (simulations_round_3_2Mb_<suffix>)
#
# The local names are what figures_and_tables/figure2_revision2.ipynb's parse_arm()
# decodes: `cmaf<digits>` and `fm<digits>` both mean "0." + digits, so cmaf01 is
# 0.01, cmaf001 is 0.001 and cmaf0 is 0. They are NOT the upstream names and the
# mapping cannot be derived mechanically -- hence this table.
ARM_MAP="
cmaf01_fm01_g5t20_2Mb:g5t20_cmaf_0.01_fm_0.01
cmaf01_fm01_g10t20_2Mb:g10t20_cmaf_0.01_fm_0.01
cmaf01_fm01_g5t30_2Mb:g5t30_cmaf_0.01_fm_0.01
cmaf001_fm01_g5t20_2Mb:g5t20_cmaf_0.001
cmaf001_fm001_g5t20_2Mb:g5t20_cmaf_0.001_fm_0.001
cmaf0_fm01_g5t20_2Mb:g5t20_cmaf_0
cmaf0_fm001_10x_psamp8000_2Mb:x10_psamp_8000_fm_0.001
cmaf0_fm001_20x_psamp8000_2Mb:x20_psamp_8000_fm_0.001
cmaf0_fm001_30x_psamp8000_2Mb:x30_psamp_8000_fm_0.001
"

WANT_ARMS="${ARMS:-}"

# One remote pass to classify every replicate, rather than a round trip each.
# `complete` needs all four: >=4 vars tables, >=4 pheno sbams, >=2 enloc.sig.out
# (one per analyzed GTEx panel) and >=3 glm lead-SNP tables.
remote_status() {
    # $(echo ...) collapses the map onto ONE line before it is interpolated: a
    # multi-line value would land mid-`for` on the remote side and be a syntax error.
    ssh "$HOST" 'bash -s' <<REMOTE
SCRATCH=$SCRATCH
for pair in $(echo $ARM_MAP); do
    arm=\${pair%%:*}; sfx=\${pair#*:}
    ROOT="\$SCRATCH/simulations_round_3_2Mb_\$sfx"
    [ -d "\$ROOT" ] || { echo "\$arm MISSINGROOT -"; continue; }
    for D in "\$ROOT"/[A-I][0-9]*; do
        [ -d "\$D" ] || continue
        ID=\$(basename "\$D")
        nv=\$(compgen -G "\$D/stage2/*/*/*_vars_*.tsv" | wc -l)
        np=\$(compgen -G "\$D/stage2/*/*/*_pheno.sbams" | wc -l)
        n4=\$(compgen -G "\$D/stage4/*/*.enloc.sig.out" | wc -l)
        n5=\$(compgen -G "\$D/stage5/*/*_glm_lead_snps.tsv" | wc -l)
        if [ "\$nv" -ge 4 ] && [ "\$np" -ge 4 ] && [ "\$n4" -ge 2 ] && [ "\$n5" -ge 3 ]; then
            echo "\$arm \$ID complete"
        else
            echo "\$arm \$ID partial v\$nv/p\$np/e\$n4/g\$n5"
        fi
    done
done
REMOTE
}

echo "Querying $HOST ..."
STATUS=$(remote_status) || { echo "ERROR: could not query $HOST" >&2; exit 1; }

n_fetch=0; n_have=0; n_wait=0; n_new_dir=0
for pair in $ARM_MAP; do
    arm="${pair%%:*}"; sfx="${pair#*:}"
    [ -n "$WANT_ARMS" ] && [[ " $WANT_ARMS " != *" $arm "* ]] && continue

    todo=""; waiting=""; have=""
    while read -r a id state extra; do
        [ "$a" = "$arm" ] || continue
        [ "$id" = "MISSINGROOT" ] && continue
        if [ "$state" != "complete" ]; then
            waiting="$waiting ${id}(${extra})"
        elif [ -d "$LOCAL/$arm/$id" ] && compgen -G "$LOCAL/$arm/$id/*.enloc.sig.out" > /dev/null 2>&1; then
            have="$have $id"
        else
            todo="$todo $id"
        fi
    done <<< "$STATUS"

    [ -z "$todo$waiting$have" ] && continue
    echo
    echo "### $arm  <-  simulations_round_3_2Mb_$sfx"
    [ -n "$have" ]    && { echo "  already local:$have";  n_have=$(($n_have + $(echo $have | wc -w))); }
    # Incomplete means running OR failed -- this cannot tell them apart, and a
    # six-day-old empty workdir looks identical to one queued a minute ago.
    # Check squeue, then the workdir's controller_*.err, before assuming patience.
    [ -n "$waiting" ] && { echo "  INCOMPLETE   :$waiting"; n_wait=$(($n_wait + $(echo $waiting | wc -w))); }
    [ -z "$todo" ]    && continue
    echo "  TO FETCH:$todo"

    if [ ! -d "$LOCAL/$arm" ]; then
        echo "  creating $LOCAL/$arm"
        [ -z "$DRY" ] && mkdir -p "$LOCAL/$arm"
        n_new_dir=$((n_new_dir + 1))
    fi

    for ID in $todo; do
        DEST="$LOCAL/$arm/$ID"
        R="$SCRATCH/simulations_round_3_2Mb_$sfx/$ID"
        # Flat into the replicate directory, matching the layout every existing arm
        # uses. The two sidecars are new: _trait_partners_ is the GWAS/GTEx pairing
        # answer key (only written where the GTEx set is topped up, i.e. power
        # sampling and the drawn-DFE arms H/I), and _synthetic_dfe_ carries the drawn
        # coefficients, which the *_vars_*.tsv files deliberately leave at zero.
        SRC="$R/stage2/*/*/*_pheno.sbams \
             $R/stage2/*/*/*_vars_*.tsv \
             $R/stage2/*/*/*_trait_partners_*.tsv \
             $R/stage2/*/*/*_synthetic_dfe_*.tsv \
             $R/stage4/*/*.enloc.sig.out \
             $R/stage5/*/*_glm_lead_snps.tsv"
        if [ -n "$DRY" ]; then
            echo "    DRY $ID"
            n_fetch=$((n_fetch + 1)); continue
        fi
        mkdir -p "$DEST"
        # --ignore-missing-args so the optional sidecars do not fail the whole
        # replicate on an arm that never writes them.
        if rsync -a --ignore-missing-args -e ssh "$HOST:$SRC" "$DEST/" 2>/dev/null; then
            printf "    %-4s %2d files\n" "$ID" "$(ls "$DEST" | wc -l | tr -d ' ')"
            n_fetch=$((n_fetch + 1))
        else
            echo "    $ID: WARN rsync returned non-zero" >&2
        fi
    done
done

echo
echo "================================================================"
echo "fetched=$n_fetch  already_local=$n_have  incomplete_upstream=$n_wait  new_dirs=$n_new_dir"
[ -n "$DRY" ] && echo "DRY=1 -- nothing transferred."
[ "$n_wait" -gt 0 ] && echo "Re-run when those finish; it only pulls what is missing. If one is not in
squeue, it FAILED rather than being slow -- read its controller_*.err."
exit 0
