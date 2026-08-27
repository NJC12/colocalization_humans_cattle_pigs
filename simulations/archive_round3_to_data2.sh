#!/bin/bash
# Collect the round-3 2 Mb analysis files out of the scratch workdirs into a
# tidy, durable archive on /n/data2. RUNS ON O2.
#
#   bash archive_round3_to_data2.sh          # copy what is complete and missing
#   DRY=1 bash archive_round3_to_data2.sh    # report only
#   ARMS="cmaf0_fm01_g5t20_2Mb" bash ...     # one arm
#
# WHY THIS EXISTS
#
# The Snakemake workdirs under /n/scratch hold the whole pipeline -- tree
# sequences, DAP-G per-locus outputs, .snakemake metadata, SLURM logs -- most of
# which is intermediate, none of which is backed up, and all of which scratch
# will eventually purge. What the analysis actually consumes is six small file
# patterns per replicate. This copies exactly those onto /n/data2 in the same
# flat `<arm>/<replicate>/` layout the local simulation_data/ tree uses, so the
# two are interchangeable and either can be re-derived from the other.
#
# It is a COPY, not a move: scratch is left untouched.
#
# WHAT COUNTS AS COMPLETE
#
# All four stages must have landed: >=4 vars tables and >=4 pheno sbams from
# stage 2, >=2 .enloc.sig.out from stage 4, >=3 glm lead-SNP tables from stage 5.
# A replicate missing any of them is REPORTED and skipped rather than
# half-copied, because a replicate directory with stage-2 output but no
# .enloc.sig.out does not read as "unfinished" downstream -- it reads as "this
# arm colocalized nothing", which is a plausible wrong number rather than an
# obvious gap.

set -uo pipefail

SRC="${SRC:-/n/scratch/users/n/njc12/snakemake}"
DEST="${DEST:-/n/data2/hms/dbmi/sunyaev/lab/nconnally/new_simulation_set}"
REPO="${REPO:-/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake}"
ARM_TSV="${ARM_TSV:-$REPO/helpers/round3_arms.tsv}"
DRY="${DRY:-}"
WANT="${ARMS:-}"

[ -f "$ARM_TSV" ] || { echo "ERROR: arm map not found: $ARM_TSV" >&2; exit 1; }

# The six patterns the analysis consumes, plus the run's parameter record.
#
# The two sidecars are recent and only some arms have them: _trait_partners_ is
# the GWAS-locus/GTEx-partner answer key, written wherever the GTEx causal set is
# topped up rather than intersected (power sampling, and the drawn-DFE arms H and
# I, and the background-selection arms K and L); _synthetic_dfe_ carries the drawn
# selection coefficients for H, I, K and L, which the *_vars_*.tsv files
# deliberately leave at zero. Absent elsewhere, and that
# absence is meaningful rather than a copy failure.
copy_replicate() {   # $1 = source replicate dir, $2 = destination dir
    local R="$1" D="$2" n=0 f
    for f in "$R"/stage2/*/*/*_pheno.sbams \
             "$R"/stage2/*/*/*_vars_*.tsv \
             "$R"/stage2/*/*/*_trait_partners_*.tsv \
             "$R"/stage2/*/*/*_synthetic_dfe_*.tsv \
             "$R"/stage4/*/*.enloc.sig.out \
             "$R"/stage5/*/*_glm_lead_snps.tsv; do
        [ -e "$f" ] || continue          # unmatched glob stays literal; skip it
        cp -p -u "$f" "$D/" && n=$((n + 1))
    done
    # Provenance. Not part of the local file set, but an archive without the
    # parameters that produced it is much harder to trust later; it is ~2 KB.
    for f in "$R"/stage2/*/*/stage2_params.txt; do
        [ -e "$f" ] && cp -p -u "$f" "$D/stage2_params.txt" && n=$((n + 1))
    done
    echo "$n"
}

complete_p() {       # $1 = source replicate dir -> 0 if all four stages present
    local R="$1"
    local nv np n4 n5
    nv=$(compgen -G "$R/stage2/*/*/*_vars_*.tsv" | wc -l)
    np=$(compgen -G "$R/stage2/*/*/*_pheno.sbams" | wc -l)
    n4=$(compgen -G "$R/stage4/*/*.enloc.sig.out" | wc -l)
    n5=$(compgen -G "$R/stage5/*/*_glm_lead_snps.tsv" | wc -l)
    STATUS_DETAIL="v$nv/p$np/e$n4/g$n5"
    [ "$nv" -ge 4 ] && [ "$np" -ge 4 ] && [ "$n4" -ge 2 ] && [ "$n5" -ge 3 ]
}

species_of() {  case "$1" in [ABCDHJKMNO]) echo human ;; [EFGILP]) echo cattle ;; *) echo "?" ;; esac; }

[ -z "$DRY" ] && mkdir -p "$DEST"
MANIFEST="$DEST/MANIFEST.tsv"
[ -z "$DRY" ] && printf 'arm\treplicate\tcategory\tspecies\tn_files\tbytes\to2_source\n' > "$MANIFEST"

n_arm=0; n_rep=0; n_skip=0
while IFS=$'\t' read -r arm sfx desc; do
    case "$arm" in ''|\#*) continue ;; esac
    [ -n "$WANT" ] && [[ " $WANT " != *" $arm "* ]] && continue
    ROOT="$SRC/simulations_round_3_2Mb_$sfx"
    [ -d "$ROOT" ] || { echo "### $arm -- source root absent, skipped"; continue; }

    echo "### $arm  <-  simulations_round_3_2Mb_$sfx"
    n_arm=$((n_arm + 1))
    got=""; skipped=""
    for R in "$ROOT"/[A-Z][0-9]*; do
        [ -d "$R" ] || continue
        ID=$(basename "$R")
        if ! complete_p "$R"; then
            skipped="$skipped ${ID}($STATUS_DETAIL)"; n_skip=$((n_skip + 1)); continue
        fi
        D="$DEST/$arm/$ID"
        if [ -n "$DRY" ]; then
            got="$got $ID"; n_rep=$((n_rep + 1)); continue
        fi
        mkdir -p "$D"
        cnt=$(copy_replicate "$R" "$D")
        bytes=$(du -sb "$D" 2>/dev/null | cut -f1)
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$arm" "$ID" "${ID:0:1}" "$(species_of "${ID:0:1}")" \
            "$cnt" "${bytes:-0}" "$ROOT/$ID" >> "$MANIFEST"
        got="$got $ID"; n_rep=$((n_rep + 1))
    done
    [ -n "$got" ]     && echo "  archived:$got"
    [ -n "$skipped" ] && echo "  INCOMPLETE, skipped:$skipped"
done < "$ARM_TSV"

# ------------------------------------------------------------------ README
if [ -z "$DRY" ]; then
cat > "$DEST/README.md" <<'DOC'
# Round-3 2 Mb colocalization simulations -- analysis file set

Generated by `simulations/archive_round3_to_data2.sh`. This is a COPY of the
analysis-relevant files from the Snakemake workdirs under `/n/scratch`; scratch
is not backed up and is periodically purged, this is not.

## Layout

    <arm>/<replicate>/<files, flat>

identical to the local `simulation_data/` tree, so the two are interchangeable.
`MANIFEST.tsv` lists every replicate with its category, species, file count,
size, and the scratch path it came from.

## Arm names

Two settings are encoded in the directory name, as digits after an implied "0.":

    cmaf<digits>   causal_min_maf -- the MAF floor causative variants were drawn from
    fm<digits>     fm_min_maf     -- the MAF floor for the variants handed to DAP-G

so `cmaf001_fm01` is causal 0.001 / fine-map 0.01, and `cmaf0` is no causal floor
at all. The rest of the name gives the phenotype-multiplier cell (`g5t20` = GWAS
x5 / GTEx x20, `10x` = both x10), the causal-sampling scheme where it is not
uniform (`psamp8000` = power-weighted, sampling n=8000), and the region size.

These are NOT the workdir names on scratch, which spell the same settings out in
full. `helpers/round3_arm.tsv` in the pipeline repo is the mapping, and it is the
only place the two dialects are tied together.

## Categories

    A  human,  directional (negative) selection on trait variants
    B  human,  NEUTRAL trait variants (effects moved onto neutral variants)
    C  human,  100k GWAS / 50k GTEx
    D  human,  as C with neutral trait variants
    E  cattle, baseline, bottlenecked
    F  cattle, positive selection, WITH continued bottlenecking
    G  cattle, positive selection, WITHOUT bottlenecking
    H  human,  NEUTRAL genealogy (pure coalescent), drawn DFE effect sizes
    I  cattle, NEUTRAL genealogy (pure coalescent), drawn DFE effect sizes
    J  human,  AFRICAN ancestry -- A with the AFR branch of the demographic
               model sampled instead of the EUR one, nothing else changed
    K  human,  BACKGROUND SELECTION -- A's genome (shaped by selection) but the
               causal loci are the NEUTRAL variants, betas drawn from the DFE
    L  cattle, BACKGROUND SELECTION -- E's genome, same recipe as K
    M  human,  FINNISH founder demography -- Wang et al. 2014 (a demes graph,
               not an stdpopsim catalog model), FIN deme sampled, Q_scaling 3
    N  human,  NON-FINNISH EUROPEAN -- the same Wang model with the NFE deme
               sampled instead. M's companion arm.
    O  human,  LOW HETEROZYGOSITY -- A's genome with a measured fraction of its
               NEUTRAL (selco == 0) sites deleted, sized so pi halves. Runs at
               A's seeds and adopts A's tree sequence, so O{N}'s variant set is a
               strict subset of A{N}'s and A - O isolates variant density
    P  cattle, LOW HETEROZYGOSITY -- E's genome, same recipe as O

H and I draw effect sizes from the SAME distribution, so H vs I isolates the
demography; H vs A and I vs E isolate the genealogy. A vs J isolates ancestry:
same species, same model, same selection, the other side of the out-of-Africa
split.

M and N are a pair for the same reason. M vs A alone would swap the founder
event and the whole demographic model at once, so N holds the model fixed:
M - N is the Finnish founder event, N - A is Wang-NFE vs Tennessen-EUR. Both
run Q_scaling 3 rather than 10 -- forced, because the FIN deme holds only 2,266
individuals at Q=10 against the 9,000 stage 1 samples, and 9,988 at Q=3.

K and L separate the two things A vs H confounds. K has A's genealogy and H's
effect model, so K vs H is the genealogy alone -- background selection -- and
A vs K is the effect assignment alone. L does the same for E vs I.

## Files per replicate

    *_vars_*.tsv               every variant in that panel: position, MAF, beta.
                               beta != 0 marks a causal variant. NOTE selco is 0
                               for H and I by design -- those variants really have
                               no selection coefficient; the drawn one is in the
                               synthetic_dfe sidecar.
    *_pheno.sbams              simulated phenotypes, DAP-G's sbams format
    *.enloc.sig.out            fastEnloc signal-level colocalization output;
                               RCP/LCP per GWAS-trait/eQTL-trait pair
    *_glm_lead_snps.tsv        plink2 GLM lead SNP per trait
    *_trait_partners_*.tsv     which GWAS causal loci have a GTEx trait at the
                               SAME causal variant. Present wherever the GTEx
                               causal set is topped up rather than intersected
                               (power sampling; the drawn-DFE arms H, I, K, L in
                               both schemes; any arm at require_gtex_partner
                               False), and also wherever the pairing rule was
                               stated explicitly -- there it is all-True, which
                               is the point: an ABSENT file is ambiguous between
                               "intersected" and "written before this existed".
    *_synthetic_dfe_*.tsv      the drawn selection coefficients (the drawn-DFE
                               arms: H, I, K, L)
    stage2_params.txt          the full parameter record for that run

Panel prefixes: `h` human / `c` cattle, then `gwas`, `gtex` (1000 individuals),
`gtex_smaller` (500), `gtex_smallest` (250, written but not analyzed).

## Completeness

Only replicates with all four stages present were copied. Anything skipped is
listed in the run log, not silently omitted.
DOC
sed -i 's/round3_arm\.tsv/round3_arms.tsv/' "$DEST/README.md" 2>/dev/null
fi

echo
echo "================================================================"
echo "arms=$n_arm  replicates=$n_rep  incomplete_skipped=$n_skip"
if [ -z "$DRY" ]; then
    echo "destination: $DEST"
    du -sh "$DEST" 2>/dev/null
    echo "manifest: $MANIFEST ($(($(wc -l < "$MANIFEST") - 1)) rows)"
else
    echo "DRY=1 -- nothing copied."
fi
exit 0
