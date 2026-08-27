# Snakemake-driven simulation pipelines

Five pipelines, one Snakefile.

## Pipelines

| Config                                   | Pipeline               | Stage 1 producer                                                  |
|------------------------------------------|------------------------|-------------------------------------------------------------------|
| `config/human.yaml`                      | human                  | `human_simulation_o2.py` (or `_larger.py`)                        |
| `config/cattle_baseline.yaml`            | cattle_baseline        | `farm_create_orig_pop_e2.py` -> `farm_burn_in_e2.slim` -> `farm_selection.slim` |
| `config/cattle_sel_bottlenecked.yaml`    | cattle_sel             | `farm_selection_from_ep8.slim` (`num_muts_selected>0`, `continue_bottlenecking=1`) |
| `config/cattle_sel_not_bottlenecked.yaml`| cattle_sel             | `farm_selection_from_ep8.slim` (`num_muts_selected>0`, `continue_bottlenecking=0`) |
| `config/human_neutral_2Mb_*_r3.yaml`     | human_neutral          | `human_neutral_simulation.py` (msprime coalescent, no SLiM)       |
| `config/cattle_neutral_2Mb_*_r3.yaml`    | cattle_neutral         | `cattle_neutral_simulation.py` (msprime coalescent, all 12 epochs, no SLiM) |

All share stages 2 (split + traits), 3 (DAP-G), 4 (fastEnloc), 5 (plink GLM).

## Quick start

```bash
# DAG dry-run (every config). Works on Mac or O2 login node.
for cfg in config/*.yaml; do
  snakemake --snakefile Snakefile --configfile $cfg -n
done

# Stop early at a specific stage
snakemake --configfile config/cattle_sel_bottlenecked.yaml --until stage2 -n
snakemake --configfile config/cattle_sel_bottlenecked.yaml --until stage3 -n
```

## Running on O2

Submit the **Snakemake controller itself** as a SLURM batch job so it lives
on a compute node, not the login node. The controller is light (~4 GB / 1 CPU)
but needs to outlast the longest stage; the `priority` partition with 7-day
walltime works.

```bash
ssh o2.hms.harvard.edu '
sbatch --parsable \
  --job-name=snake_cattle_sel_bot_sdYYYYMMDD \
  --partition=priority \
  --time=7-00:00:00 \
  --mem=4G \
  --cpus-per-task=1 \
  --output=/n/scratch/users/n/njc12/snakemake/cattle_sel_bot/logs/controller_%j.out \
  --error=/n/scratch/users/n/njc12/snakemake/cattle_sel_bot/logs/controller_%j.err \
  --chdir=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake \
  --wrap="/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake \
    --configfile config/cattle_sel_bottlenecked.yaml \
    --config stage1_seed=YYYYMMDD stage2_seed=YYYYMMDD stage3_seed=YYYYMMDD stage4_seed=YYYYMMDD stage5_seed=YYYYMMDD \
    --profile profiles/o2 -j 200 \
    --keep-going --rerun-incomplete"
'
```

Notes:
- Set the date / seed to whatever's appropriate. `stage{1,2}_seed` are the only
  ones currently consumed (see "Per-stage seeds" below).
- `--keep-going` continues other independent branches when a single rule fails.
- `--rerun-incomplete` reruns any rule whose outputs were partially produced
  (useful after a controller crash).

### When a previous controller died via scancel/SIGKILL

Snakemake leaves a directory lock that blocks new controllers. Run
`snakemake --configfile <cfg> --unlock` once from the snakemake dir on O2,
then resubmit.

### Submissions log

Every submission has been recorded in `~/comparative_colocalization/snakemake_submissions.txt`
with the job id, command, status, and (for failed runs) the post-mortem.

## On-cluster prerequisites

Issues in setting up on o2 (HMS compute cluster)

| Concern | Resolution |
|---|---|
| Python env for stages 1/2 | `coloc_sims` (with `tskit`, `pyslim`, `msprime`, `tstrait`, `matplotlib`, `seaborn`, `snakemake>=9`, `snakemake-executor-plugin-slurm`). All installed via `mamba install -n coloc_sims …`. |
| Snakemake version | 9.x in `coloc_sims`. The older `smilenfer` env's 6.x snakemake is not compatible with the SLURM executor plugin we use. |
| `dap-g` shared libs | `dap-g` (and `fastenloc`) need `libgsl.so.28` from GSL 2.8 + a recent `libstdc++`. A separate conda env `fastenloc_env` holds those (`mamba create -n fastenloc_env -c conda-forge gsl=2.8 libstdcxx-ng`). `coloc_sims` is pinned at `gsl=2.7` because `msprime`/`pyslim` depend on it -- don't mix the two envs. |
| `dap2enloc` | Perl script at `/home/njc12/bin/fastenloc/utility/dap2enloc` (not `summarize_dap2enloc.pl`; this naming is enforced as a config default in the Snakefile). |
| `plink2` | `/home/njc12/bin/plink2`. Built for AVX2 — compute nodes are fine; login nodes (SSE4.2 only) will refuse to even print `--version`. The rules submit it to a compute node so this is fine in practice. |
| `slim`, `tabix`, `bgzip` | All at known cluster paths; defaults are in `Snakefile`'s `config.setdefault(...)` block. |

The Snakefile sets these as `config.setdefault(...)`-style defaults so YAMLs
don't have to repeat them. Override per-pipeline in the YAML if necessary.

## Where existing data is reused

Each config lists per-stage **search dirs**. This is because I added this snakemake
on top of my existing (not streamlined) simulation setup.
At parse time the Snakefile walks them and:

1. Runs `helpers/alias.py` to symlink legacy filenames into the canonical
   names that the rules expect (e.g.
   `revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24.full.ts`
   → `cattle_sel_bot_mult100_gen23_muts26_sd24.full.ts`).
2. For stage-2 search dirs: if a directory matching the expected
   `gwas_{S}_gtex_{S}_maf_{M}/` shape is found and **all required files**
   are present (counting both the top-level and `plink_analysis/` sub-dir
   where legacy `.pheno` files live), symlinks them into the canonical
   workdir/stage2/ location and touches the marker. Otherwise leaves stage 2
   to re-run.

If the workdir is read-only (e.g. doing a dry-run on a non-O2 host while
the YAML points at `/n/scratch`), parse-time adoption is a no-op and prints
a notice -- nothing breaks.

### Resuming from a previous workdir (stages 3-5)

For the cheap-to-redo stages 1+2 the rules look up individual files via search
dirs. For stages 3-5 (where the work to redo is *expensive*) the Snakefile
supports a different mechanism — adopting whole subtrees from one or more
**previous workdirs**:

```yaml
# In any pipeline YAML:
previous_workdirs:
  - /n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake/cattle_sel_bot   # e.g., rsynced from /n/scratch
  - /n/scratch/users/n/njc12/snakemake/cattle_sel_bot_v1                # an older run
```

At parse time the Snakefile walks each entry's `stage3/<run_tag>/`,
`stage4/<run_tag>/`, and `stage5/<run_tag>/` subtrees and symlinks every
file found into the current workdir, skipping any path that already exists.
`<run_tag>` includes the basename, selection params, and seed, so adoption
is automatically scoped to runs with the **same seed configuration** as the
current run.

This is per-file granularity. Practical implications:

- **Stage 3 dap-g** is per-trait, so partial adoption helps proportionally —
  if you have 100/323 cgwas dap-g outputs from a prior run, those 100 are
  reused and only the remaining 223 get re-submitted.
- **Stage 4 fastenloc** produces all 5 `.enloc.*.out` files from a single
  rule invocation. Partial adoption (some `.enloc.*.out` files present,
  others not) means Snakemake re-runs fastenloc anyway. In practice you
  either have the whole set or none.
- **Stage 5 plink GLM** runs `plink2 --glm` once per category to produce
  many `.glm.linear` files plus a `.stage5.done` sentinel; only the sentinel
  is a tracked output. If the sentinel is in the previous workdir, stage 5
  is skipped wholesale; otherwise re-runs.

Concrete use cases this unlocks:
- Controller died mid-run: same workdir still works (Snakemake's native
  in-workdir existence check). No need for `previous_workdirs` -- this is
  belt-and-suspenders.
- `/n/scratch` purged: rsync your `<workdir>` to `/n/data2` before purge,
  set `previous_workdirs: [/n/data2/.../snakemake/cattle_sel_bot]`,
  submit. Stage 1+2 outputs are re-found via the regular search dirs and
  search-dir lookups; stage 3-5 outputs are re-found via this mechanism.
- Diff'd workdir: you want a clean workdir but reuse expensive stage 3
  outputs from an old workdir. Same as above.

### Permanent override: cattle_sel pipelines always start from `epoch_8.ts`

Both `cattle_sel_bottlenecked.yaml` and `cattle_sel_not_bottlenecked.yaml`
have `stage1_search_dirs: []` so the rule never short-circuits to a
pre-existing `revision_*.full.ts`. The `epoch_8.ts` checkpoint at
`/n/data2/.../comparative_colocalization/simulations/selection/farm_selection_Q_0.01.L_10000000.seed_20250303.epoch_8.ts`
is still found via `cattle_baseline_search_dirs`, and the SLiM job runs
fresh from there for every submission.

## The three MAF floors

There are three independent MAF cutoffs. They are set separately, none is derived
from another, and each governs exactly one thing:

| Config key | Governs | Where it is applied | Default |
|---|---|---|---|
| `causal_min_maf` | Which variants may be **selected as causative** — central pool (in **both** sampling schemes), central neutral recipients (categories B/D), and flanking GTEx alike | `create_gwas_files_and_phenotypes.py`, stage 2 | `0` = no floor |
| `fm_min_maf` | Which variants are **fine-mapped** (the SBAMS handed to DAP-G) | `stage3_index_geno`'s awk filter | `0` = no filter |
| `min_maf` | Which variants enter the **GWAS** | `plink2 --maf`, stage 5 | required key |

The interesting combination is a causal floor *below* the tested floors:

```
--config causal_min_maf=0.001 fm_min_maf=0.01     # min_maf stays 0.01
```

`causal_min_maf` is a **floor with no ceiling**, so this does not make the run
unfindable — most of the MAF ≥ 0.001 pool is still above 0.01. It splits the loci:

- causal MAF ≥ 0.01 → in the SBAMS and in the GLM, behaves exactly as at baseline;
- causal MAF in [0.001, 0.01) → in neither tested set, findable only through a
  common variant that tags it.

`helpers/summarize_coloc.py` splits its rows on `causal_tested` accordingly, so the
second group does not drag down the first. For a shared locus `causal_tested=True`
means MAF ≥ `fm_min_maf` in **both** samples — the same test the baseline pool
applies — so that subset is directly comparable to a baseline run.

**Size the trait count for it.** `n_central_traits` is drawn from the *lowered* pool,
so the tested subset is only `n_central_traits × P(MAF ≥ 0.01 | MAF ≥ 0.001)` loci.
That fraction is far smaller for human than for bottlenecked cattle (88% of human
variants in a ±250 kb window are below MAF 0.01, against cattle's 9% — see the
comment in `rules/common.smk`), so at the usual `n_central_traits: 50` the human
arm's `causal_tested=True` row rests on the fewest loci in exactly the comparison
that matters. Raise `n_central_traits` for these arms to keep it near 50, or read
the `replicates` and pool counts before trusting a small-n row.
`submit_2Mb_r3_causal_maf.sh` launches the arms.

`causal_min_maf` is the only one of the three that changes what stage 2 produces,
so it is the one that appears in paths — and it appears in **two** places:

```
stage2/hts_11.cmaf_0.001/gwas_35_gtex_35_maf_0.001/   <- dir named for the causal floor
stage3/hts_11.cmaf_0.001/hgwas/...                       (so search-dir adoption
stage4/hts_11.cmaf_0.001/A1_human.hgtex.enloc.sig.out     cannot cross floors)
stage5/hts_11.cmaf_0.001/hgwas_glm.*
```

The `.cmaf_<v>` segment comes from `stage2_run_tag()`, the single namespace for
stages 2–5, and is **empty at 0.01** (`paths.LEGACY_CAUSAL_MIN_MAF`) so every path
produced before this key existed is unchanged. It composes with `output_tag`
rather than replacing it: `output_tag` still means "an analysis variant of one
fixed stage-2 output" (r2, fine-mapping floor), which is why stage 5 remains
tag-agnostic but does carry the causal segment.

Two guards, both at parse time, before any job is submitted:

- A stage-2 directory whose `stage2_params.txt` records a different
  content-determining parameter is refused, naming the key and both values.
- A stage-2 directory that predates `stage2_params.txt` is refused for any run
  whose causal floor is not `0.01`, since its provenance cannot be verified.
  At `0.01` it is silently accepted — that is what every existing run is.

Unlike the fine-mapping variants, these arms **cannot** reuse stage 2
(`stage2_search_dirs`), because the causal floor decides which variants become
traits at all. Stage 1 is still shared, and must be, or the arms are not
comparable.

### Semantics change: the floor defaults to 0 and always applies

`causal_min_maf` used to default to `0.01` **and** to be skipped by the central
pool under `causal_sampling: power`. Both changed together:

- **Default is now `0`** (`paths.DEFAULT_CAUSAL_MIN_MAF`), i.e. no floor.
- **A floor that is set applies in both sampling schemes**, to the central pool,
  the GTEx top-up candidates, the `neutral_trait_vars` recipient pool, and the
  flanks.
- `0` is not "no filter": `maf > 0` well-formedness survives it. A site
  monomorphic *within a panel* is not a variant there, and an effect assigned to
  it would produce a phenotype of pure noise. That test used to be implicit in a
  non-zero floor.

**Every config in `config/` pins `causal_min_maf: 0.01` explicitly**, so the
default change is inert for them and every stage-2..5 path they produce is
byte-identical to before (verified across all 55). A hand-rolled `--config`
invocation that omits the key now gets `0` rather than `0.01`.

`paths.LEGACY_CAUSAL_MIN_MAF = 0.01` survives as the **path-suppression sentinel
only**, decoupled from the default — which is what lets the default move without
renaming every directory ever written. A run at the new default therefore *does*
emit `.cmaf_0`. Do not re-couple the two.

#### The `psamp` arms were relabelled `maf_0`

`submit_2Mb_r3_x20_psamp_fm001.sh` never passed `causal_min_maf`; it inherited
`0.01` and its own header said the value was "NOT a floor on the causal pool under
power sampling". Under the new rule that content is a `causal_min_maf: 0` run, so
the directories were relabelled to say so (`helpers/migrate_cmaf_psamp.py`), and
the script now defaults `CAUSAL_FLOOR=0` to keep reproducing them.

What makes the relabel safe rather than a guess: the new predicate at floor 0 is
*exactly* the old no-floor predicate. Checked on a 24k-site fixture — central pool
4024 vs 4024, top-up candidates 2778 vs 2778, identical position sets. The only
rows that move are 4 flank sites monomorphic in the GTEx panel, which the new
`maf > 0` term drops and which could never have been eQTLs.


## Causal-variant sampling (`causal_sampling`)

How the `n_central_traits` causal loci are picked out of the eligible pool. Two
schemes; `uniform` is the default and is what every run through round 3 used.

| Value | The draw |
|---|---|
| `uniform` | Uniform sample of the pool that clears `causal_min_maf` in **both** panels. The pool is intersected with the GTEx panel, so the GWAS and shared-GTEx causal sets are the same positions by construction. |
| `power` | Every eligible central variant is weighted by its probability of being **detected** in a GWAS of `sampling_gwas_n` individuals, and the draw has inclusion probability proportional to that power. |

A uniform draw is dominated by variants no realistic GWAS could find — 52–64% of
the human pool sits below MAF 0.01 — so the simulated "GWAS loci" are not the kind
of loci a real study reports. `power` picks loci the way a GWAS does.

The model (`helpers/causal_power.py`). The phenotype is `y = ±β·g + N(0,1)` with
`β = sqrt(|selco|) × gwas_scaling` and `g ∈ {0,1,2}`, never standardized, so

```
vexp  = 2f(1-f)β²          ncp = N · vexp / (1 + vexp)
power = P(χ²₁(ncp) > qchisq(sampling_sig_p, 1, lower = FALSE))
```

the same identity as `add_gwas_significance()` in `figure2_revision2.ipynb` and
`helpers/pval_rescale.py`. The draw is **πPS** (exact first-order inclusion
probabilities), not `rng.choice(p=weights, replace=False)`, which is successive
sampling and over-represents the high-power tail relative to its own weights.

| Config key | Meaning | Default |
|---|---|---|
| `causal_sampling` | `uniform` or `power` | `uniform` |
| `sampling_gwas_n` | GWAS size the power is computed for | `gwas_size` (8000 for A/E) |
| `sampling_sig_p` | What "detected" means | `5e-8` |
| `sampling_min_power` | Power floor used by the guard below | `0.05` |
| `sampling_min_pool_multiple` | Refuse unless this × `n_central_traits` variants clear that floor | `2` |

`sampling_gwas_n` is a knob in its own right: raise it to ask "which loci would a
larger study find?" and the causal set shifts toward rarer, smaller-effect
variants **without** changing the simulated GWAS.

Four things that behave differently under `power`, all deliberate:

- **`power` requires an explicit `n_central_traits`.** At the default
  `causal_min_maf: 0` the pool is every polymorphic central variant, so "use all
  eligible" would make the whole region causative.

  `causal_min_maf` **does** gate the central pool here. It used to be skipped
  under `power`, on the grounds that the detection-power weight subsumed a
  frequency floor. It does not — a weight is a soft preference, a floor is a hard
  cut — and having one knob mean two things depending on another knob made "the
  pool at `causal_min_maf=0.01`" ambiguous. See the semantics note below.
- **The two panels' causal sets can diverge.** The GWAS draw is not intersected
  with the GTEx panel, so the shared GTEx set is only whichever drawn loci that
  panel carries; the rest of the central GTEx slots are filled uniformly from
  central `selco != 0` GTEx variants. **`helpers/summarize_coloc.py` and
  `figure2_revision2.ipynb` do not yet know this** — both define a true positive
  as trait-name equality and keep every GWAS trait in the denominator, so a
  partnerless locus reads as a colocalization failure, and as a *false positive*
  if it colocalizes with a neighbour. Read the partner count out of the stage-2
  log before trusting any power or FDR number from a `power` run.
- **The causal set differs between multiplier cells.** `β` scales with
  `gwas_scaling`, so g5t20 / g10t20 / g5t30 are no longer paired
  variant-for-variant.
- **Stage 2 can refuse to run**, with the counts, when too little of the pool is
  detectable. A draw from a pool where only a handful of variants carry real
  weight is a few certainties plus an arbitrary tail, not a weighted sample.

Two sidecars are written next to the phenotypes, only in this mode:

```
{h,c}gwas_causal_power_gwas_<g>_gtex_<t>_maf_<c>.tsv    whole pool: power, pi, selected
{h,c}gwas_trait_partners_gwas_<g>_gtex_<t>_maf_<c>.tsv  gwas_trait -> gtex_trait, shared
```

The partner table is the explicit pairing a corrected scorer should read instead
of inferring it from trait names.

Paths carry a `.psamp_<sampling_gwas_n>` segment from `stage2_run_tag()`, composing
with `.cmaf_<v>` and empty under `uniform`, so every pre-existing path is
unchanged:

```
stage2/hts_11.psamp_8000/gwas_10_gtex_20_maf_0.01/
stage3/hts_11.psamp_8000/hgwas/...
```

Only `sampling_gwas_n` is in the path. Two `power` runs differing solely in
`sampling_sig_p` land on the same path and are caught by the stage-2 provenance
guard instead — loudly, rather than silently reusing the wrong phenotypes. All
five keys are in `STAGE2_KEYS`, so the guard covers them; records written before
these keys existed are skipped rather than flagged, and are unaffected.

`submit_2Mb_r3_x20_psamp_fm001.sh` launches the x20 arm (fine-mapping and GLM
floors at 0.001). Its uniform counterpart is `submit_2Mb_r3_cmaf01_control.sh`
run at the same cell and floors — `CELL=x20 FM_FLOOR=0.001 GLM_FLOOR=0.001` —
which differs only in which 50 central variants were made causative.


## Per-simulation parameters file

Every real invocation writes the resolved config — after defaults and every
`--config` override — to two identical files:

```
{workdir}/params/{run_tag}[.{output_tag}].params.txt      # from onstart, survives a dead DAG
{stage4_dir}/{basename}[.{output_tag}].params.txt         # beside the outputs it describes
```

Mapping an output back to its parameters is textual: for any
`{prefix}.{gtex_cat}.enloc.{kind}.out`, strip `.{gtex_cat}.enloc.{kind}.out` and
append `.params.txt`. Uniqueness comes from the directory — `stage4_dir` embeds
`run_tag` — not from `basename`, which is `human` in every round-3 config.

Stage 2 additionally writes `{stage2_inner}/stage2_params.txt`, which records the
resulting **pool sizes** as well as the parameters, and which travels with the
directory when it is symlink-adopted into another output root. That is what the
provenance guards above read.

The format is a flat YAML mapping under a `.txt` extension: readable, `yaml.safe_load`-able,
and re-feedable with `snakemake --configfile <params.txt>`. `_meta` carries host,
timestamp, Snakemake version and the repo's git commit; `_derived` carries the
computed paths, the file list this record describes, and two hashes — `stage2_uid`
(over the stage-2-determining keys, so it can be checked against the stage-2
sidecar after adoption) and `run_uid` (over those plus the analysis keys).

Written at `onstart`, not as a rule output: it adds no job to the DAG, changes
nothing for runs already complete, and correctly writes nothing under `-n`. Same
reasoning as the `.fmmaf` sidecar. Two runs that share an `output_tag` but differ
in an analysis parameter are caught here too — the older record is archived to
`.params.<timestamp>.txt` with a warning, or raises with `params_strict: true`.


## Per-stage seeds

Each config has `stage{1..5}_seed`.

- **`stage1_seed`** — passed to SLiM (`slim -s {seed}`) or to
  `human_simulation_o2.py --seed`. The search-dir lookup uses it to
  fail-loudly on seed mismatch (see below).
- **`stage2_seed`** — passed to `create_gwas_files_and_phenotypes.py --seed`.
  Currently only affects the `add_neutral` msprime overlay; the per-variant
  trait simulation in `tstrait.sim_phenotype(...)` still uses
  `random_seed=position`, so trait sim is deterministic per causative variant
  regardless of `stage2_seed`. To make the trait sim seed-driven too, patch
  line 106 of `create_gwas_files_and_phenotypes.py`.
- `stage3_seed`, `stage4_seed`, `stage5_seed` — reserved YAML keys; the
  underlying tools (dap-g, fastenloc, plink2 GLM) are effectively deterministic
  for these inputs, so they have no effect. Kept in the YAMLs for parallelism
  and future use.

## Seed-mismatch behavior

If a search dir contains a file with the right shape (e.g. `*_sd24.full.ts`)
but the config requests a different seed (`stage1_seed: 99`), the pipeline
fails with `WorkflowError` rather than silently regenerating.
Either update the config or remove the search dir.

## Trait selection cap (optional)

Bounds the number of traits fine-mapped (and colocalized via fastEnloc) per
replicate, sharing the budget across a category of replicate runs
(e.g. A1, A2, A3, A4). Add three optional fields to each replicate's config:

```yaml
selection_cap: 1500        # total cap across replicates in the category
selection_seed: 42         # RNG seed (shared across replicates)
num_replicates: 4          # category size (e.g. 4 for A1-A4)
```

When set, the `stage3_manifest` checkpoint emits a per-replicate quota of
`ceil(selection_cap / num_replicates)` trait IDs, drawn from a
`(selection_seed, stage1_seed)`-seeded permutation of the pheno file's traits.
Quota = `ceil(1500/4) = 375` per replicate; total ≈ `selection_cap`, slightly
over-shooting when the division is not exact.

The same trait subset is used across `hgwas`, `hgtex`, and `hgtex_smaller`
within a replicate, because all three tiers derive from the same
`ccausative_maf01` position set in stage 2.

**Extending the cap later.** Raise `selection_cap` (e.g. 1500 → 1700) and rerun
snakemake. The permutation is fixed by the seed, so the new quota (`ceil(1700/4) = 425`)
extends the previous prefix — already-completed DAP-G/fastEnloc outputs are
reused, and only the newly-selected traits run.

**Caveats:**
- Replicates with fewer available traits than their quota take all available
  and do not redistribute. Actual total fine-mapped can fall below
  `selection_cap` if individual replicates are small.
- `selection_cap` not divisible by `num_replicates` over-shoots by up to
  `num_replicates − 1` traits in aggregate.
- If `selection_cap` is unset, the checkpoint falls back to the previous
  behavior (fine-map every trait in the pheno file).

## Workdir vs publishdir

- `workdir`: where intermediate jobs run; defaults to `/n/scratch/...`. Fast,
  but temporary storage
- `publishdir`: persistent destination on `/n/data2/...`. **No `publish:` rule
  exists yet** — the pipeline writes everything to workdir, and rsync to
  `/n/data2` is a manual step at the moment.

## Files / directories

```
simulations/
  Snakefile                         # top-level dispatch + parse-time symlinking
                                    # + binary-path defaults
  config/                           # one YAML per pipeline
  rules/
    common.smk                      # stages 2, 3, 4, 5 (shared)
    stage1_human.smk
    stage1_cattle_baseline.smk
    stage1_cattle_sel.smk
  helpers/
    paths.py                        # canonical filename builders
    search_dirs.py                  # multi-directory lookup with seed-mismatch check
    alias.py                        # legacy-filename -> canonical symlinks
    params_record.py                # per-simulation parameters files + provenance
    summarize_coloc.py              # cross-cell colocalization summary table
  scripts/
    dapg_one_locus.sh               # per-trait DAP-G runner (stage 3 fan-out)
    run_fastenloc.sh                # stage 4 wrapper
    run_plink_glm.sh                # stage 5 wrapper
  envs/
    coloc_sims.yaml                 # informational spec; pipeline does NOT
                                    # use snakemake's --use-conda
  profiles/o2/
    config.yaml                     # SLURM executor profile
```

## Resources (per-rule)

Stages can be split across partitions because their cost varies wildly. The
current settings come from real observation, not guesses:

| Rule | Partition | Walltime | Memory | CPUs |
|---|---|---|---|---|
| `stage1_human` | priority | 30 d | 200 GB | 4 |
| `stage1_cattle_create_pop` | short | 4 h | 50 GB | 4 |
| `stage1_cattle_burn_in` | priority | 30 d | 150 GB | 4 |
| `stage1_cattle_selection` | priority | 30 d | 150 GB | 4 |
| `stage1_cattle_sel` (cattle+sel from `epoch_8`) | short | 8 h | 32 GB | 4 |
| `stage2_split_pheno` | short | 8 h | 64 GB | 4 |
| `stage3_index_geno` | short | 4 h | 32 GB | 8 |
| `stage3_dapg_locus` — **gwas** category | **medium** | **48 h** | **200 GB** | **16** |
| `stage3_dapg_locus` — **gtex** category | short | 4 h | 8 GB | 4 |
| `stage4_annot_vcf` | default (short) | 4 h | 8 GB | 2 |
| `stage4_fastenloc` | short | 2 h | 16 GB | 4 |
| `stage5_plink_glm` | short | 4 h | 16 GB | 4 |

Per-category branching on stage-3 dap-g uses Snakemake callable resources:

```python
slurm_partition = lambda wc: "medium" if wc.cat.endswith("gwas") else "short",
runtime         = lambda wc: 48 * 60  if wc.cat.endswith("gwas") else 4 * 60,
mem_mb          = lambda wc: 200000   if wc.cat.endswith("gwas") else 8000,
cpus_per_task   = lambda wc: 16       if wc.cat.endswith("gwas") else 4,
```

This split is necessary because dap-g on a ~99k-sample × ~5k-SNP locus
window comfortably exceeds 32 GB during matrix construction and runs for
roughly 30-40 wall-clock hours; the gtex side (~1k samples) finishes in
seconds on 8 GB.

## Partition policy reminders (O2)

- `short` accepts jobs ≤ 12 h.
- `medium` accepts jobs > 12 h.
- A 12-hour job on `medium` will be rejected with `Job not submitted: please
  submit jobs that are less than 12 hours long to the short (or priority)
  partition.` All rule walltimes are tuned to avoid this trap.
- `priority` accepts very long jobs (used here for the controllers).

## Stage 4 fastEnloc `-t` (total SNPs)

Defaulted to `"auto"` in every YAML. The wrapper computes:

```
comm -12 <(ids in cgwas.dap.vcf.gz) <(ids in cgtex.dap.vcf.gz) | wc -l
```

via a header-flexible awk (skips lines starting with `#` and lines whose
position column isn't purely numeric, so it works on both VCFs with and
without `##fileformat=…` headers). Override by setting
`fastenloc_total_snps: <int>` in the YAML.

## What's not yet wired up

- **Publish to /n/data2**: no `rule publish:`. Outputs live on `/n/scratch`
  until manually rsynced. Add a final rsync rule if you want persistence.
- **Stage-3/4/5 search-dir reuse**: only stage 1 and stage 2 read search dirs.
  Stages 3-5 always regenerate (cgwas dap-g is by far the most expensive part,
  so this matters if you ever want to reuse it -- not currently supported).
- **Crediting a tagging variant**: `summarize_coloc.py` scores a true positive by
  the *identity* of the causal variant, so an individual locus whose causal variant
  fell below `fm_min_maf` can only ever report 0% TP. Those loci are isolated on
  their own `causal_tested=False` rows, which makes the table legible, but nothing
  yet credits a credible set for containing a variant that merely tags the truth.
  A window match (against `dapg_window`) or an LD match (r2 from the stage-2 VCF,
  which still holds the causal variant) would need new code.
- **`fetch_big_results_2Mb.sh` and causal-MAF arms**: it pins
  `SUB=stage2/*/gwas_35_gtex_35_maf_0.01` and rsyncs flat into `$LOCAL/$CAT/$ID/`,
  so it will not match a `maf_0.001` directory -- and widening the glob would
  flatten two arms into one local directory. Parameterise it per output root
  before fetching a causal-MAF run.


## Submitting on O2 issues

1. **The controller must be a `sbatch` job, not a foreground process.** A
   foreground `snakemake … --profile profiles/o2 -j 200` on the login node
   dies when ssh disconnects, leaving the SLURM children orphaned.
2. **`--use-conda` was a trap.** With `use-conda: True`, the controller's
   first action is `mamba env create`, which exceeded 4 GB on the controller
   sbatch. Removed in this repo.
3. **Stage 2's partition was wrong.** Asked for `medium` with 12-h walltime;
   was rejected because `medium` is for jobs *longer* than 12 h. Changed to
   `short` with 8 h.
4. **`matplotlib`/`seaborn` were missing from `coloc_sims`.** Imported at
   module level in `create_gwas_files_and_phenotypes.py` but never used.
   Installed via `mamba install -n coloc_sims -c conda-forge matplotlib seaborn`.
5. **`libgsl.so.28` missing for both `dap-g` and `fastenloc`.** Fixed by
   prepending `LD_LIBRARY_PATH=/home/njc12/miniconda3/envs/fastenloc_env/lib`
   in `dapg_one_locus.sh` (via `--dapg-libpath`) and `run_fastenloc.sh`
   (via `--fastenloc-libpath`).
6. **OOM on cgwas dap-g at 32 GB.** Fixed by per-cat branching (cgwas/hgwas
   → 200 GB/16 CPU/medium/48h). The legacy `submit_revision_dapg_o2.sh`
   profiles were correct; I'd started with the small-profile values.
7. **`LockException` after `scancel`** of a previous controller. Fix:
   `snakemake --unlock` on the workdir, then resubmit.

# List of simulations

| Population | GWAS size | GTEx size       | Selection on trait associated variants | Bottlenecking  | ID | Pipeline                      | Replicates |
| ---------- | --------- | --------------- | -------------------------------------- | -------------  | -- | ----------------------------- | ---------- |
| Human      |     8,000 | 1,000, 500, 250 | Directional (negative)                 | NA             | A  |                         human |          4 |
| Human      |     8,000 | 1,000, 500, 250 | Neutral *(legacy — see below)*         | NA             | B  |                         human |          4 |
| Human      |   100,000 |          50,000 | Directional (negative)                 | NA             | C  |                         human |          1 |
| Human      |   100,000 |          50,000 | Neutral *(legacy — see below)*         | NA             | D  |                         human |          1 |
| Cattle     |     8,000 | 1,000, 500, 250 | Directional (negative)                 | yes            | E  | cattle_baseline_from_midpoint |          4 |
| Cattle     |     8,000 | 1,000, 500, 250 | Directional (positive)                 | yes            | F  |       cattle_sel_bottlenecked |          4 |
| Cattle     |     8,000 | 1,000, 500, 250 | Directional (positive)                 | no             | G  |   cattle_sel_not_bottlenecked |          4 |
| Human      |     8,000 | 1,000, 500, 250 | None anywhere (neutral genealogy)      | NA             | H  |                 human_neutral |          5 |
| Cattle     |     8,000 | 1,000, 500, 250 | None anywhere (neutral genealogy)      | yes            | I  |                cattle_neutral |          5 |
| Human (AFR)|     8,000 | 1,000, 500, 250 | Directional (negative)                 | NA             | J  |                         human | 5 *(planned)* |
| Human      |     8,000 | 1,000, 500, 250 | In the GENOME, not on trait variants   | NA             | K  |                         human | 5 *(planned)* |
| Cattle     |     8,000 | 1,000, 500, 250 | In the GENOME, not on trait variants   | yes            | L  | cattle_baseline_from_midpoint | 5 *(planned)* |
| Human (FIN)|     8,000 | 1,000, 500, 250 | Directional (negative)                 | founder event  | M  |                         human | 5 *(planned)* |
| Human (NFE)|     8,000 | 1,000, 500, 250 | Directional (negative)                 | NA             | N  |                         human | 5 *(planned)* |
| Human      |     8,000 | 1,000, 500, 250 | Directional (negative), pi HALVED       | NA             | O  |                         human | 5 *(planned)* |
| Cattle     |     8,000 | 1,000, 500, 250 | Directional (negative), pi HALVED       | yes            | P  | cattle_baseline_from_midpoint | 5 *(planned)* |

## Categories K and L: background selection

A and H differ in two things at once. A's genome comes from a forward run under
the DFE and its causal loci are the selected variants; H's genome is a pure
coalescent and its causal loci are neutral variants carrying drawn effect sizes.
So A - H is a genealogy difference and an effect-assignment difference summed
together, and nothing in the category set separated them.

K and L are the missing cell. They take **A's and E's genomes** -- forward runs
under the DFE, so the genealogy carries background selection -- and put **H's and
I's effect model** on top: the causal loci are drawn from the strictly NEUTRAL
variants (`selco == 0`) and their betas come from the truncated DFE rather than
from the variant's own coefficient. That gives a 2x2:

|                    | causal = selected variants | causal = neutral, drawn betas |
| ------------------ | -------------------------- | ----------------------------- |
| genome under selection | A / E                  | **K / L**                     |
| genome neutral (coalescent) | --                | H / I                         |

so **K - H isolates the genealogy** (background selection alone, effect model
held fixed) and **A - K isolates the effect assignment** (genealogy held fixed).
L - I and E - L do the same on the cattle side.

Mechanically they are A's and E's configs plus one override,
`synthetic_dfe_effects=True`, exactly as B is A's config plus
`neutral_trait_vars=True`. The stage-1 genetic model is byte-identical to A/E.

One subtlety in the implementation. `synthetic_dfe_effects` used to mean "drop
the `selco != 0` requirement on the causal pool", which is fine for H and I --
their coalescent genomes have no selected variants at all, so "all variants" and
"the neutral variants" are the same set. On A's and E's genomes those differ, and
dropping the filter would have let K and L draw *selected* variants and then
assign them a second, unrelated coefficient. The predicate is now `selco == 0`
when synthetic, which is provably inert for H and I (0 of 39,842 H1 variants and
0 of 4,196 I1 variants carry a coefficient) and load-bearing for K and L.


## Categories M and N: the Finnish founder pair

M and N are the first arms to run a demography that is **not in stdpopsim's
catalog**. `demographic_model: FinnishWang2014` names an entry in
`helpers/human_demography.MODELS`, which points at the demes graph
`demography/finnish_demo_wang_2014.yaml` — Wang, Agarwala, Flannick, Chiang,
Altshuler & Hirschhorn (2014), *Am J Hum Genet* 94:710–720, Table S2 class 3.
`human_simulation_o2.py` turns it into a `stdpopsim.DemographicModel` via
`msprime.Demography.from_demes` and hands it to the same SLiM engine every other
human arm uses; the contig, the DFE, the mutation and recombination rates and
every one of stages 2–5 are unchanged.

**They are a pair, and only make sense read as one.**

| contrast | what it isolates |
| -------- | ---------------- |
| M − N    | the Finnish founder event. Same model, same Q, same seed rule, same everything — only which deme is retained differs. |
| N − A    | the model swap: Wang-2014 NFE against Tennessen-2012 EUR. |

M against A on its own would swap the founder event *and* the whole model at
once (root Ne 8,100 vs 7,310, present-day 111,394 vs 512,000, different growth
schedules), which is why N exists at all. The two config files are identical
except for `basename`, `population`, the two seeds and the workdir/publishdir
leaf, and `helpers/tests/test_finnish_demography.py` asserts exactly that.

### Why `Q_scaling` is 3 and not 10

This is forced, not a preference, and it is the one number in these arms that
cannot be changed casually. stdpopsim's SLiM engine refuses to sample more
individuals than the **Q-rescaled** deme holds at the sampling tick
(`slim_engine.py`: *"Request to sample N individuals from pX … but only M
individuals will be alive"*) — and it finds out at the very end, after the whole
forward simulation has run. Stage 1 samples `gwas_size + 1000 = 9,000`:

    individuals alive at the sampling tick, MEASURED from SLiM:

      Q      FIN      NFE     verdict for a 9,000-individual draw
     10    2,266   10,653     FIN fails
      5    4,531   21,306     FIN fails
      4    7,491   27,437     FIN fails
      3    9,988   36,583     both OK -- FIN by a margin of 988

Note these are NOT `present_size / Q`. stdpopsim rounds epoch boundaries to
whole ticks, which truncates the tail of the exponential growth, so the naive
figure overshoots by 13–34% (it promises 11,486 at Q=3 where SLiM delivers
9,988). `helpers/human_demography.SLIM_CAPACITY` holds the measured numbers and
the recipe for measuring more; `test_finnish_demography.py` refuses a config
whose (model, population, Q) has not been measured rather than falling back to
the arithmetic.

NFE would in fact fit at Q=10 (10,653 alive). It runs at 3 anyway, because M − N
is only the founder event if the two arms are Q-matched.

SLiM cost scales as 1/Q², and the burn-in dominates: stdpopsim's default
`slim_burn_in=10` gives 10 × 8,100 / 3 = 27,000 ticks over 2,700 diploids,
against A's 7,310 ticks over 731. That is roughly 11× A's stage 1, so ~12 h
rather than ~1 h. `rules/stage1_human.smk` routes any small-`L` run below Q=10
to `medium` with 4 days and 16 GB for exactly this reason.

Two consequences worth stating rather than discovering:

* `Q_scaling` is in `params_record.STAGE2_KEYS`, so M and N can never silently
  adopt A's stage 2 — they build their own, which is what the `^(K|L|M|N)$`
  branch in `submit_2Mb_r3_cmaf_fm001.sh` arranges.
* M and N run at a **different Q than A**, so an M-vs-A or N-vs-A difference is
  not entirely free of rescaling artifact (Q=3 is the more accurate of the two).
  M − N is internally Q-matched, which is the other reason the pair exists.

Two smaller things the scaling does. Sampling 9,000 of FIN's 9,988 is a 90%
sample where A takes 17.6% of EUR — relatedness is not the problem that looks
like (the expected number of full-sib pairs in a Wright-Fisher deme that size is
under one), but the sample's site-frequency spectrum is essentially the
population's, which matters when MAF distributions are compared against A's. And
the model's 13-generation fast-growth epoch is 4.33 ticks at Q=3, since stdpopsim
rounds epoch boundaries to whole ticks. That rounding is exactly what costs the
1,498 individuals above, and it makes the shape of the final expansion coarse.

### Filenames

M's stage-1 tree sequence is `hts_finwang_FIN_{seed}.ts` and N's is
`hts_finwang_NFE_{seed}.ts`. Both halves of that name are load-bearing: the tag
because M and N share a model, the deme because that is all they differ by. For
human filenames the seed half of stage-1 adoption never fires — `search_dirs`
matches `sd\d+|seed_\d+`, which no `hts_*` name contains — so the filename is
the whole guard, and `params_record.STAGE2_KEYS` omits every demography key on
exactly that basis.

### The three CHOSEN parameters

The demes file marks three values as chosen within ranges the paper reports:
slow-phase growth 2%/generation (range 0.5–5%), fast-phase growth 15% (range
8–30%) and NFE→FIN migration 2% (range 0.5–7%). Together they leave present-day
Finnish Ne very loosely determined — the corners span ~4,200 to ~2.1 million.
Treat them as a sensitivity axis. The registry is keyed by model id precisely so
a corner run can be added as a **second entry** with its own tag and its own
`hts_*` namespace, rather than by editing the file underneath results that
already exist. Note that changing the growth rates moves the present-day size
that `Q_scaling: 3` was chosen against; `test_finnish_demography.py` pins both
numbers and checks the cap for every config.

## Categories O and P: low heterozygosity

Every other contrast on the category axis changes *which* variants are causal, or
*which* demography produced them. None of them changes how much neutral variation
the genome carries -- and that is not a nuisance parameter here. Variant density
sets LD-tagging, the size of the DAP-G candidate set and therefore fine-mapping
difficulty, and the two species differ in it enormously: E carries 1,365
segregating sites at 2 Mb where A carries ~36,000. Any A-vs-E number is partly a
density difference, and nothing separated it out.

**O is A with pi halved, and P is E with pi halved.** Nothing else differs -- not
the demography, not the DFE, not the seeds, not the causal loci.

### Why it is a deletion and not a lower mutation rate

Neutral mutations are not part of either species' forward simulation. They are a
Poisson overlay on branches the forward run has already fixed: 8.4e-9 in
`human_simulation_o2.py` (stage 1, so they arrive inside `hts_*.ts`) and 8.4e-9
in `create_gwas_files_and_phenotypes.add_neutral` (stage 2, in two passes across
the cattle Q handoff). Dropping each neutral site independently with probability
`1 - k` is the thinning of a Poisson process, and a thinned Poisson process is a
Poisson process at the thinned rate.

So `neutral_keep_fraction: k` is not an approximation of "re-simulate at
`k x 8.4e-9`" -- conditional on the genealogy it has exactly that distribution.
Three things follow, and they are the whole reason the arm is built this way:

* **One mechanism for both species.** The human overlay is a stage-1 literal and
  the cattle one a stage-2 literal, but the thinning happens in stage 2 after
  `add_neutral` and after `remove_fixed`, so a single code path
  (`helpers/neutral_thinning.py`) serves both.
* **Stage 1 is reused, not rerun.** O_i adopts A_i's `hts_1i.ts` and P_i adopts
  E_i's `farm_selection_from_ep8...seed_5i.full.ts`.
* **The causal architecture is untouched, provably.** `causal_eligible`,
  `flank_eligible`, `select_central_power` and `select_gtex_topup` all select on
  `selco != 0`; the thinning only ever deletes `selco == 0`. So the pools, the
  pi-PS weights, the 50 drawn central positions, the flank loci, the effect sizes
  and the `*_trait_partners_*.tsv` tables come out **identical** to the parent
  category's at the same seed. `helpers/tests/test_lowhet_categories.py` pins
  that the two modules read the same metadata field by the same route, which is
  the one way it could quietly stop being true.

### They share A's and E's seed bands, deliberately

This is the only place the `seed = 10 * letter_index + replicate` rule is broken,
and it is broken on purpose. O{N} runs at seed 1{N} and P{N} at 5{N} -- their
parents' seeds -- because that is what makes O_i's variant set a **strict subset**
of A_i's: the same genealogy, the same individuals in the same panels, the same
causal loci, with some of the background removed. A - O therefore carries no
replicate noise at all, which no other pair on this axis manages.

Two consequences, both enforced:

* **O and P must never run their own stage 1.** A private seed band would make
  them do exactly that and silently turn a paired contrast into an unpaired one.
  `stage1_donor_of()` in the two submit scripts redirects the stage-1 lookup
  (`O<i> -> A<i>`, `P<i> -> E<i>`), the stage-1 filename builders read
  `stage1_seed` and not `basename` so the adoption still resolves, and
  `submit_2Mb_r3_cmaf_replicates.sh` -- the script that *does* build stage 1 --
  has no O or P entry at all.
* **O and P must never adopt A's or E's stage 2.** `neutral_keep_fraction` is in
  `params_record.STAGE2_KEYS`, so the provenance guard refuses a mismatched
  directory loudly, and `stage2_run_tag` carries a `.nkeep_<token>` segment so
  they cannot land on the same path in the first place. In
  `submit_2Mb_r3_cmaf_fm001.sh` they sit in the `^(K|L|M|N|O|P)$` fresh-stage-2
  branch.

### The keep fraction is measured, not chosen

pi is additive over sites and thinning does not move any surviving site's
frequency, so

```
E[pi(k)] = k * pi_neutral + pi_selected
```

and halving pi is solved, not searched:

```
k = 0.5 * (1 - pi_selected / pi_neutral)
```

`helpers/measure_pi_components.py` replays the stage-2 preamble on a stage-1 tree
sequence -- including the cattle `add_neutral` overlay, without which the cattle
selected share would read as ~1.0 -- and reports both components and the implied
k. `calibrate_lowhet.sh` runs it over all ten donors and prints the mean per
species, which is the line to paste into the config:

```bash
DRY=1 bash calibrate_lowhet.sh     # resolve every donor path, submit nothing
bash calibrate_lowhet.sh           # one short job; writes helpers/pi_components_lowhet.tsv
```

**k differs between the species** (cattle's selected sites carry a different
share of its pi), which is why the two configs carry different numbers rather
than a shared constant.

Measured 2026-08-27 (job 51523508, `helpers/pi_components_lowhet.tsv`):

| | selected share of pi | k, mean over 5 | k range | pinned | EXPECTED pi ratio |
|---|---:|---:|---|---:|---|
| A1-A5 (-> O) | 0.2794 | 0.3061 | 0.2991-0.3163 | **0.306** | 0.4924-0.5049 |
| E1-E5 (-> P) | 0.2948 | 0.2909 | 0.2777-0.3056 | **0.291** | 0.4895-0.5092 |

The last column is `E[pi(k)]/pi(1)` evaluated at the pinned k -- an expectation,
not an outcome. The deletion is a random subset, so the realized ratio scatters
around it: the O1 and P1 pilots both landed at **0.508**, slightly above the
expected 0.502 and 0.500 and outside the tabulated range. That scatter is the
draw, not a mis-set k, and `thin_pi_ratio` in each replicate's
`stage2_params.txt` is the number to read for what actually happened.

The two species land close together, which is a coincidence of this cell and not
a reason to share one constant. Note also how much *more* than half the neutral
class this removes -- ~69-71% of it -- which is the composition effect below,
quantified. Do this BEFORE launching anything: both submit scripts
refuse an O or P whose config has no `neutral_keep_fraction`, because the key
defaults to keep-everything and the run would otherwise be a byte-identical copy
of its parent under a different name -- a null result that looks like a
measurement.

If `pi_selected >= 0.5 * pi_total` for a species, k comes out negative: no amount
of neutral thinning halves that species' pi. `keep_fraction_for_pi_target` raises
with the numbers rather than clamping, because a clamped k would deliver a
different experiment under the same name.

### What the target costs in composition -- quote this with any O - A number

Halving pi removes **more** than half the neutral sites, because pi is dominated
by the neutral class, which is also the common half of the frequency spectrum. So
the total site count falls by more than the pi target, and the selected:neutral
ratio shifts toward selected: O and P are **enriched for rare, low-MAF selected
variants** relative to A and E.

That is a consequence of sizing the cut against pi rather than against the neutral
count, not a defect -- but it means O is not "A with fewer markers" in every
respect, and a difference in, say, the MAF distribution of the fine-mapping
candidate set is expected rather than surprising.

At the measured fractions it is a large effect, so here it is in numbers rather
than in the abstract. A1 carries 38,174 sites, 23,164 of them neutral (60.7%);
keeping 0.306 of those leaves ~22,100 sites, a **42% cut in site count** for a
50% cut in pi, and moves the selected share from 39% to **68%**. E1 goes from
1,366 sites (925 neutral, 67.7%) to ~710, with the selected share moving from
32% to **62%**. Both arms are therefore majority-selected where their parents
were majority-neutral. `stage2_params.txt` records
`thin_sites_before` / `thin_sites_after` / `thin_neutral_sites_before` /
`thin_neutral_sites_kept` / `thin_pi_before` / `thin_pi_after` / `thin_pi_ratio`
per replicate, so the realized numbers are always beside the outputs.

### Where they run

They are new **categories**, not new **arms**: both parameter sets are existing
rows of `helpers/round3_arms.tsv` that already hold A1-A5 and E1-E5, so that file,
`fetch_round3_2Mb.sh` and the notebook's `parse_arm()` need no changes.

```bash
# cmaf 0.001 / fm 0.001 / GLM 0.001, g5t20   -> cmaf001_fm001_g5t20_2Mb
REPS="O1 O2 O3 O4 O5 P1 P2 P3 P4 P5" CELL=g5t20 \
    bash submit_2Mb_r3_cmaf_fm001.sh

# power-weighted n=100000, no plateau        -> cmaf0_fm001_g5t20_psamp100000_2Mb
REPS="O1 O2 O3 O4 O5 P1 P2 P3 P4 P5" CELL=g5t20 SAMPLING_N=100000 \
    PLATEAU=1 CAUSAL_FLOOR=0 FM_FLOOR=0.001 GLM_FLOOR=0.001 \
    bash submit_2Mb_r3_x20_psamp_fm001.sh
```

The pool guard on the psamp arm behaves identically to A's and E's, incidentally
-- it counts central `selco != 0` variants clearing power 0.05, and thinning does
not touch them.

## Category J: the African-ancestry arm

J is category A with one parameter changed. stdpopsim's `HomSap` /
`OutOfAfrica_2T12` (Tennessen et al. 2012) has two populations and simulates
**both** forward regardless of which is sampled; every human arm through I took
`{"AFR": 0, "EUR": n}` and J takes `{"AFR": n, "EUR": 0}`. Same SLiM run, same
DFE, same recombination and mutation rates, same sample size, same trait
architecture, same stages 2-5, same cost. Only the retained individuals differ.

That makes A − J the ancestry contrast and nothing else. It sits alongside the
other single-parameter contrasts the category axis is built out of: B − A
isolates the effect assignment, H − A the genealogy, H − I the demography, and
J − A the LD and allele-frequency structure that colocalization power is
actually a function of. The EUR branch carries the out-of-Africa bottleneck, so
J should show more diversity, more common variation and shorter-range LD.

Two things to be precise about:

* This is the African **branch of a two-population model**, not an isolated
  African population — 2T12 has migration between AFR and EUR. `Africa_1T12`
  would be the isolated model, and would make the *model* differ from A rather
  than only the sampling.
* J's stage-1 tree sequence is `hts_AFR_{seed}.ts`, not `hts_{seed}.ts`.
  Stage-1 adoption matches on filename plus embedded seed, and a J run is
  otherwise indistinguishable from an A run, so the population is in the name
  for the same reason `nts_`/`cnts_` exist. EUR keeps the bare `hts_{seed}.ts`
  so every tree sequence written before the key existed stays adoptable. The
  key is `population:` (helpers/human_demography.py); it is read only by the
  `human` pipeline and the Snakefile refuses it anywhere else rather than
  ignoring it.

## The neutral models: B, H and I

They answer different questions and are not interchangeable.

**B (and D) keep a genome produced by a forward run UNDER selection.** Only the
*assignment* of effects is neutral: each `selco != 0` donor's coefficient is moved
onto a random `selco == 0` recipient and the trait is named for the recipient
(`build_redistribution_map`). The genealogy still carries background selection, and
the effect size still comes from a real coefficient — so a locus's effect size is
tied to the donor's frequency history rather than to the recipient's.

**H removes selection from the genome itself.** Stage 1 is a pure msprime
coalescent under the same `OutOfAfrica_2T12` demography A uses, with no SLiM phase
and no selected mutations at all (`human_neutral_simulation.py`). The causal
variants' effect-size parameters are *drawn* in stage 2 from the DFE with its
neutral class removed (`helpers/synthetic_dfe.py`, proportions 0.975/0.024/0.001
over m1/m2/m3 — the neutral class was never a DFE component, it is the separate
8.4e-9 overlay, so dropping it and renormalizing lands on the g1 vector verbatim).
`beta = sqrt(|s|) * multiplier` exactly as in A, so effect sizes stay on A's scale.

So: **B − A isolates the effect assignment; H − A isolates the genealogy.**

B is retained for comparison and is fully runnable. It is marked legacy only in the
sense that H is the model we expect to use going forward; nothing about it has been
moved or disabled.

**I is H with the cattle demography**, and that is the whole difference. The DFE
they draw from is the same object — `human_simulation_o2.py:67-92` and
`farm_selection.slim:40-45` define identical components with identical proportions,
so `helpers/synthetic_dfe.py` needs no species variant. So **H − I isolates the
demography**, and **I − E isolates the genealogy on the cattle side** exactly as
H − A does on the human side.

Two things about I differ from E/F/G rather than from H:

*It simulates all twelve epochs.* E/F/G resume from the shared
`farm_selection_*.ep7.ts` checkpoint because 29,800 ticks of burn-in plus epochs 2–7
is unaffordable to re-run per replicate in SLiM. A coalescent has no burn-in — epoch
1 is simply the ancestral state — so I runs the whole bottleneck from N=17,000 down
to N=90 itself, in seconds. Its configs therefore carry neither `cattle_baseline_seed`
nor `cattle_baseline_search_dirs`, and the pre-flight handoff guard in
`submit_2Mb_r3_cmaf_replicates.sh` stays scoped to `^(E|F|G)$`.

*Its `mutation_rate` is 5.6e-9, not 1.4e-8.* This is the opposite of H and it is not
a typo — see below.

### Why H's rate is the total and I's is the 5.6e-9 component

The two families put the neutral overlay in different places. The human arms apply
both halves in stage 1, so H must ask for the 1.4e-8 total (see the A1 measurement
below). The cattle arms apply 5.6e-9 in SLiM and let
`create_gwas_files_and_phenotypes.py:add_neutral` overlay 8.4e-9 in stage 2 — and
that overlay runs **unconditionally** on the `--cattle_ts_file` branch, ignoring
`already_includes_neutral` (which every cattle config sets to `True`, so making the
branch honour it would silently switch E/F/G's overlay off). So I reproduces the same
split: 5.6e-9 × Q in `cattle_neutral_simulation.py`, 8.4e-9 × Q in stage 2, total
1.4e-8, and **no change to shared stage-2 code**. Both classes are neutral in I, so
the split carries no meaning beyond landing on E/F/G's total.

### What I measures against E

Measured at L = 2 Mb, seed 91, `causal_min_maf: 0.01`, uniform sampling, both
`cgwas` panels at 8,000 individuals:

| | segregating sites | implied Ne | MAF < 0.01 | MAF ≥ 0.05 | mean MAF |
| --- | ---: | ---: | ---: | ---: | ---: |
| E1 (forward, under selection) | 1,365 | 1,188 | 11.4% | 64.3% | 0.146 |
| I1 (coalescent, no selection) | 4,196 | 3,652 | 1.6% | 83.7% | 0.212 |

**I carries 3.07× E's variation, and E is also skewed toward rare variants.** Both
are what background selection does, and the direction is the same one H shows
against A — but the magnitude is far larger here, which is the interesting part.
A's site count is dominated by young singletons thrown up by the out-of-Africa
expansion, and BGS barely touches those, so A − H came to only 1.10× in site count
(against 2.17× in π). Cattle has no such expansion: the deep history is 29,800
generations at a flat N = 17,000 with the whole 2 Mb under the DFE, so the
reduction lands on the site count directly. Read the 3.07× as an upper bound on
BGS rather than a clean estimate of it — E's deep phase also ran at `Q_scaling: 1`
and only epochs 8–12 at 0.01, whereas I is single-Q throughout.

Note θ_W implies Ne = 3,652 while π implies ~10,980. That gap (positive Tajima's D)
is the bottleneck, and it is present in I with no selection at all, which is what
makes it usable as the null.

Partner coverage is not a problem for I the way it is for H at a zero floor: at
`causal_min_maf: 0.01` the run gives 50 central GWAS loci, 50 shared GTEx, 0 top-up
and 50 flank — **50/50 partners**, no top-up needed. Still read
`cgwas_trait_partners_*.tsv` before trusting a coloc rate at any other cell.

### Why I runs at `Q_scaling: 0.01` when H runs at 1

Under the Hudson coalescent the two are the same distribution: scaling every size and
duration by 1/Q and every rate by Q leaves ∫dt/2N, µt and rt unchanged, and the tree
sequences are identical (verified — the implied real Ne agrees to four decimal places
between Q=0.01 and Q=1). What the value states is which process the run *claims* to
approximate. E/F/G are Wright-Fisher populations of 9,000 whose entire membership is
sampled; Hudson is the standard approximation to that. At `Q_scaling: 1` the same
9,000-individual sample would be drawn from a population of 90 — an approximation to
a process nobody simulates. `Q_scaling: 0.01` also keeps stage 2 correct without an
exemption: `get_vars_df` multiplies the tick times by it, and `add_neutral` divides
its overlay rate by `1/it`. There is no `handoff_ticks` — that exists for E/F/G,
whose tree sequences span the Q=1 → Q=0.01 deep-history boundary and so carry a
piecewise time scale. I is a single coalescent at a single Q.

### H's mutation rate is the TOTAL, not the neutral component

A mutates at 1.4e-8: SLiM applies the DFE at 5.6e-9 and msprime overlays neutral
variation at 8.4e-9, and both classes end up in the tree sequence stage 2 reads.
(Confirmed on A1: 60.6% of its 35,929 `hgwas` variants have `selco == 0`, against
the predicted 8.4/14 = 60%.) H therefore runs at `mutation_rate: 1.4e-8` — matching
*mutational input*, which is the thing that is not itself a consequence of
selection. Simulating only the 8.4e-9 neutral component would give ~40% fewer
variants and change fine-mapping difficulty alongside the genetic model.

Measured against A1 at L = 2 Mb, H then carries **~10% more segregating sites and
~2.2× the nucleotide diversity**. That gap is the background selection A has and H
does not, and it is a result rather than a miscalibration: a mis-set rate would move
site count and π by the *same* factor, and these move by 1.10 and 2.17. Under a
recent expansion most segregating sites are young singletons, whose count barely
depends on the deep coalescent rate BGS suppresses, while π depends on it strongly.

### Run H with a causal MAF floor, or with power sampling — not uniform-at-zero

H does not intersect its causal pool with the GTEx panel (it uses `select_gtex_topup`
in both sampling schemes, per the model as specified: draw 50 in the GWAS panel,
whichever GTEx carries become causal there, top up the rest). A GWAS locus the GTEx
panel lacks has no partner, and both scorers count that as a colocalization failure.
Measured on H1 at L = 2 Mb, x35:

| sampling | `causal_min_maf` | central pool | GTEx shared | top-up | partners |
| -------- | ---------------- | -----------: | ----------: | -----: | -------: |
| uniform  | 0.01             |        1,299 |          50 |      0 |    50/50 |
| power    | 0.01             |        1,299 |          50 |      0 |    50/50 |
| power    | 0                |       20,097 |          46 |      4 |    46/50 |
| uniform  | 0                |       20,097 |          10 |     40 |  **10/50** |

At a 0.01 floor a common GWAS variant almost always segregates in a
1,000-individual GTEx panel, so the issue does not arise. At a zero floor the pool
is dominated by very rare variants; power weighting still lands on common ones, but
a uniform draw does not, and four fifths of the loci end up with nothing to find.
The configs pin `causal_min_maf: 0.01` for this reason. Stage 2 prints the partner
count and writes `<panel>_trait_partners_*.tsv`; neither scorer reads it yet.

### Round-3 notes on B, F and G

**B has no config of its own.** It is A's `config/human_2Mb_<cell>_r3.yaml` plus a
single `--config neutral_trait_vars=True` override, which redistributes each
donor's effect onto a random neutral recipient and names the trait for the
recipient. One file, one genetic model.

**B under `causal_sampling: power` loses a GTEx partner or two, and so does A.**
Under uniform sampling the neutral recipient pool is intersected with the GTEx
panel and floored at `causal_min_maf`, so every GWAS trait is guaranteed a partner.
Power sampling drops that intersection for *every* category, so a GWAS trait can
sit where the GTEx panel carries nothing — which both scorers still count as a
colocalization failure rather than as "nothing to find". Measured in
`…_2Mb_x30_psamp_8000_fm_0.001`: A1 48/50, A2 47/50, B1 48/50; cattle F1 and G1
both 50/50, their panels being far less sparse. So this is a property of the
sampling scheme, not of B, and it does not bias an A-vs-B comparison at one cell.
Stage 2 prints the count and writes the answer key to
`<panel>_trait_partners_*.tsv`; neither scorer reads that file yet.

Do not re-introduce a GTEx intersection on B's recipient pool to "fix" this: it
would give B a partner rate A does not have, and the A−B difference would then be
partly a partner-rate artefact.

**B × power crashed before 2026-08-04 and the outputs of that combination do not
exist anywhere.** `combine_phenos_to_df` resolved the *donor's* selection
coefficient against the panel it was phenotyping, which is the GTEx frame for the
shared central loci — and under power sampling the donor pool is not intersected
with that panel, so ~76% of donors are absent from it and the lookup raised
`KeyError`. The donor's selco now travels on its own row (it is mutation metadata,
identical in every panel that carries it); only the recipient, whose site id must
be a row of the tree sequence being phenotyped, is still resolved against `vars`.

**Round-3 F and G run at `num_muts_selected: 5`, not the 26 rounds 1 and 2 used.**
26 was sized for the 10 Mb region and was carried onto 2 Mb unchanged, ~5× the
intended density of positively-selected variants per Mb; `26 × (2/10) = 5.2 → 5`
restores it. This deliberately breaks comparability with round-1/2 F/G. The value
is embedded in the stage-1 filename (`..._muts5_...`), so a `muts5` run can never
adopt or collide with a `muts26` tree sequence.

**F and G have no round-3 stage 1 anywhere**, unlike A and E. Whichever arm runs
them first must use `controller_2Mb.sbatch` (which does not require `STAGE1_SRC`);
later arms reuse that tree sequence via `stage1_search_dirs` on the command line.
Both read the same shared `…ep7.ts` handoff that E does, via
`cattle_baseline_search_dirs`.

## List of snakemake submission commands

The block below submits all 22 controller jobs in one shell session.
SLURM dependencies (`--dependency=afterany:...`) hold C and D until
every other controller has exited, so the original table order is
preserved here even though the two heavyweight runs actually execute
last. Each run gets its own flat publishdir
(`/n/data2/.../simulations_for_revision/<ID><rep>/`) and matching
scratch workdir, so the 22 runs never collide and existing in-flight
controllers (e.g. `snakemake/cattle_sel_bot/`) keep running undisturbed.

Seeds encode the run identity:

    seed = 10 * letter_index + replicate
    (A=1, B=2, …, I=9, J=10, K=11, L=12, M=13, N=14)

So A1 → 11; D1 → 41; F2 → 62; G4 → 74; H1 → 81; I1 → 91; **J1 → 101**;
K1 → 111; L1 → 121; M1 → 131; N1 → 141. Through
I this reads off as "tens digit = letter, ones digit = replicate", and that
shorthand is what stops working at J — the arithmetic does not. J's band is
101–109, which is disjoint from the 11–99 above it; K is 111–119, L 121–129,
M 131–139 and N 141–149, and the rule keeps producing disjoint bands
indefinitely.

(The `seed_of()` implementations concatenate rather than add, so the two agree
for one-digit replicates — which is every replicate that has ever been run —
and diverge at replicate 10. Nothing has ever gone past 5.)

```bash
ssh o2.hms.harvard.edu << 'EOF'
set -e
SIMS_PUB=/n/data2/hms/dbmi/sunyaev/lab/nconnally/simulations_for_revision
SIMS_SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_for_revision

mkdir -p \
  $SIMS_SCRATCH/{A1,A2,A3,A4,B1,B2,B3,B4,C1,D1,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4}/logs \
  $SIMS_PUB/{A1,A2,A3,A4,B1,B2,B3,B4,C1,D1,E1,E2,E3,E4,F1,F2,F3,F4,G1,G2,G3,G4}

SUBMIT_DIR=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SNAKEMAKE=/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake

# Generic submit: ID CFG BASE SEED GWAS_SIZE GTEX_SIZE QV NEUTRAL EXTRA DEP
submit() {
  local ID=$1 CFG=$2 BASE=$3 SEED=$4 GWAS=$5 GTEX=$6 QV=$7 NEUTRAL=$8 EXTRA=$9 DEP=${10}
  sbatch --parsable \
    --job-name=snake_${ID}_sd${SEED} \
    --partition=long --time=7-00:00:00 --mem=4G --cpus-per-task=1 \
    --output=$SIMS_SCRATCH/$ID/logs/controller_%j.out \
    --error=$SIMS_SCRATCH/$ID/logs/controller_%j.err \
    --chdir=$SUBMIT_DIR \
    $DEP \
    --wrap="$SNAKEMAKE --snakefile $SUBMIT_DIR/Snakefile \
      --directory $SIMS_SCRATCH/$ID \
      --configfile $SUBMIT_DIR/$CFG \
      --config basename=${ID}_${BASE} stage1_seed=$SEED stage2_seed=$SEED \
               gwas_size=$GWAS gtex_size=$GTEX Q_scaling=$QV neutral_trait_vars=$NEUTRAL \
               workdir=$SIMS_SCRATCH/$ID publishdir=$SIMS_PUB/$ID $EXTRA \
      --profile $SUBMIT_DIR/profiles/o2 -j 200 --keep-going --rerun-incomplete"
}

# ---- A: human, 8k GWAS, multi-GTEx, directional negative selection ----
JOB_A1=$(submit A1 config/human.yaml human 11 8000 -1 10 False ""); echo "A1 $JOB_A1"
JOB_A2=$(submit A2 config/human.yaml human 12 8000 -1 10 False ""); echo "A2 $JOB_A2"
JOB_A3=$(submit A3 config/human.yaml human 13 8000 -1 10 False ""); echo "A3 $JOB_A3"
JOB_A4=$(submit A4 config/human.yaml human 14 8000 -1 10 False ""); echo "A4 $JOB_A4"

# ---- B: human, 8k GWAS, multi-GTEx, neutral trait variants ----
JOB_B1=$(submit B1 config/human.yaml human 21 8000 -1 10 True ""); echo "B1 $JOB_B1"
JOB_B2=$(submit B2 config/human.yaml human 22 8000 -1 10 True ""); echo "B2 $JOB_B2"
JOB_B3=$(submit B3 config/human.yaml human 23 8000 -1 10 True ""); echo "B3 $JOB_B3"
JOB_B4=$(submit B4 config/human.yaml human 24 8000 -1 10 True ""); echo "B4 $JOB_B4"

# C and D are submitted at the bottom with a dependency on all other JOB_*.

# ---- E: cattle baseline-from-midpoint, multi-GTEx, directional negative ----
JOB_E1=$(submit E1 config/cattle_baseline_from_midpoint.yaml cattle_baseline_from_midpoint 51 8000 -1 0.01 False ""); echo "E1 $JOB_E1"
JOB_E2=$(submit E2 config/cattle_baseline_from_midpoint.yaml cattle_baseline_from_midpoint 52 8000 -1 0.01 False ""); echo "E2 $JOB_E2"
JOB_E3=$(submit E3 config/cattle_baseline_from_midpoint.yaml cattle_baseline_from_midpoint 53 8000 -1 0.01 False ""); echo "E3 $JOB_E3"
JOB_E4=$(submit E4 config/cattle_baseline_from_midpoint.yaml cattle_baseline_from_midpoint 54 8000 -1 0.01 False ""); echo "E4 $JOB_E4"

# ---- F: cattle, bottlenecked, directional positive selection ----
JOB_F1=$(submit F1 config/cattle_sel_bottlenecked.yaml cattle_sel_bot 61 8000 -1 0.01 False ""); echo "F1 $JOB_F1"
JOB_F2=$(submit F2 config/cattle_sel_bottlenecked.yaml cattle_sel_bot 62 8000 -1 0.01 False ""); echo "F2 $JOB_F2"
JOB_F3=$(submit F3 config/cattle_sel_bottlenecked.yaml cattle_sel_bot 63 8000 -1 0.01 False ""); echo "F3 $JOB_F3"
JOB_F4=$(submit F4 config/cattle_sel_bottlenecked.yaml cattle_sel_bot 64 8000 -1 0.01 False ""); echo "F4 $JOB_F4"

# ---- G: cattle, non-bottlenecked, directional positive selection ----
JOB_G1=$(submit G1 config/cattle_sel_not_bottlenecked.yaml cattle_sel_nobot 71 8000 -1 0.01 False ""); echo "G1 $JOB_G1"
JOB_G2=$(submit G2 config/cattle_sel_not_bottlenecked.yaml cattle_sel_nobot 72 8000 -1 0.01 False ""); echo "G2 $JOB_G2"
JOB_G3=$(submit G3 config/cattle_sel_not_bottlenecked.yaml cattle_sel_nobot 73 8000 -1 0.01 False ""); echo "G3 $JOB_G3"
JOB_G4=$(submit G4 config/cattle_sel_not_bottlenecked.yaml cattle_sel_nobot 74 8000 -1 0.01 False ""); echo "G4 $JOB_G4"

# ---- C and D: held last via --dependency=afterany on every prior controller ----
DEP="--dependency=afterany:$JOB_A1:$JOB_A2:$JOB_A3:$JOB_A4:$JOB_B1:$JOB_B2:$JOB_B3:$JOB_B4:$JOB_E1:$JOB_E2:$JOB_E3:$JOB_E4:$JOB_F1:$JOB_F2:$JOB_F3:$JOB_F4:$JOB_G1:$JOB_G2:$JOB_G3:$JOB_G4"

# C: human, 50k GWAS, 50k GTEx, directional negative selection, L=2.5 Mb
JOB_C1=$(submit C1 config/human.yaml human 31 50000 50000 1 False "L=2500000" "$DEP"); echo "C1 $JOB_C1"

# D: human, 50k GWAS, 50k GTEx, neutral trait variants, L=2.5 Mb
JOB_D1=$(submit D1 config/human.yaml human 41 50000 50000 1 True  "L=2500000" "$DEP"); echo "D1 $JOB_D1"
EOF
```

Monitor progress with:

```bash
ssh o2 'squeue -u $USER --format="%.10i %.30j %.2t %.10M %.5D %R"'
ssh o2 'sacct -u $USER --format=JobID,JobName,State,Elapsed,MaxRSS --starttime=$(date -d "1 day ago" +%Y-%m-%d)'
```

# Re-run stages 3-4 at r2 = 0.25 (parallel result set)

Goal: produce a second set of stage-3 (DAP-G) + stage-4 (fastEnloc) outputs at
`ld_ctrl=0.25` next to the existing `ld_ctrl=0.75` results, for A/B/E/F/G (no C
or D). Stages 1, 2, and 5 are reused from the original runs.

The new outputs are namespaced by `output_tag=r2_0_25`:

| File                          | r2 = 0.75 (default)                                | r2 = 0.25 (new)                                                |
| ----------------------------- | -------------------------------------------------- | -------------------------------------------------------------- |
| DAP-G loci                    | `stage3/<run>/<cat>/outputs/{trait}.dapg.out`      | `stage3/<run>/<cat>/outputs_r2_0_25/{trait}.dapg.out`          |
| DAP-G logs                    | `stage3/<run>/<cat>/logs/{trait}.dapg.out`         | `stage3/<run>/<cat>/logs_r2_0_25/{trait}.dapg.out`             |
| stage-3 done sentinel         | `stage3/<run>/<cat>/.stage3.done`                  | `stage3/<run>/<cat>/.stage3.r2_0_25.done`                      |
| fastEnloc per gtex_cat        | `stage4/<run>/{basename}.{gtex_cat}.enloc.*.out`   | `stage4/<run>/{basename}.r2_0_25.{gtex_cat}.enloc.*.out`       |

The original (r2=0.75) files are untouched. Each rep reuses its existing
`$SIMS_SCRATCH/<ID>/` workdir so stages 1/2 and the stage-3 genotype index +
manifest are detected as up-to-date.

## Resource-estimation phase (10 jobs)

Submit one snakemake controller per rep in A1/B1/E1/F1/G1, each targeting just
two stage-3 DAP-G outputs (the first trait of each category's manifest). After
they complete, read `MaxRSS` / `Elapsed` from `sacct` and scale by 1.5x to size
the production submission.

```bash
ssh o2.hms.harvard.edu << 'EOF'
set -e
SIMS_PUB=/n/data2/hms/dbmi/sunyaev/lab/nconnally/simulations_for_revision
SIMS_SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_for_revision
SUBMIT_DIR=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SNAKEMAKE=/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake

# Test: ID CFG BASE SEED Q NEUTRAL GWAS_CAT GTEX_CAT
submit_r2_test() {
  local ID=$1 CFG=$2 BASE=$3 SEED=$4 QV=$5 NEUTRAL=$6 GWAS_CAT=$7 GTEX_CAT=$8
  # Resolve stage2_run_tag (single subdir of stage3/ in this workdir) and the
  # first trait of each category's manifest.
  local RUN_TAG TRAIT_G TRAIT_T
  RUN_TAG=$(basename $(ls -1d $SIMS_SCRATCH/$ID/stage3/*/ | head -1))
  TRAIT_G=$(head -1 $SIMS_SCRATCH/$ID/stage3/$RUN_TAG/${GWAS_CAT}/manifest.txt)
  TRAIT_T=$(head -1 $SIMS_SCRATCH/$ID/stage3/$RUN_TAG/${GTEX_CAT}/manifest.txt)
  sbatch --parsable \
    --job-name=snake_${ID}_r2_025_test \
    --partition=short --time=8:00:00 --mem=4G --cpus-per-task=1 \
    --output=$SIMS_SCRATCH/$ID/logs/controller_r2_025_test_%j.out \
    --error=$SIMS_SCRATCH/$ID/logs/controller_r2_025_test_%j.err \
    --chdir=$SUBMIT_DIR \
    --wrap="$SNAKEMAKE --snakefile $SUBMIT_DIR/Snakefile \
      --directory $SIMS_SCRATCH/$ID \
      --configfile $SUBMIT_DIR/$CFG \
      --config basename=${ID}_${BASE} stage1_seed=$SEED stage2_seed=$SEED \
               gwas_size=8000 gtex_size=-1 Q_scaling=$QV neutral_trait_vars=$NEUTRAL \
               ld_ctrl=0.25 output_tag=r2_0_25 \
               workdir=$SIMS_SCRATCH/$ID publishdir=$SIMS_PUB/$ID \
      --profile $SUBMIT_DIR/profiles/o2 -j 2 --rerun-incomplete \
      $SIMS_SCRATCH/$ID/stage3/$RUN_TAG/${GWAS_CAT}/outputs_r2_0_25/${TRAIT_G}.dapg.out \
      $SIMS_SCRATCH/$ID/stage3/$RUN_TAG/${GTEX_CAT}/outputs_r2_0_25/${TRAIT_T}.dapg.out"
}

submit_r2_test A1 config/human.yaml                          human                          11 10   False hgwas hgtex
submit_r2_test B1 config/human.yaml                          human                          21 10   True  hgwas hgtex
submit_r2_test E1 config/cattle_baseline_from_midpoint.yaml  cattle_baseline_from_midpoint  51 0.01 False cgwas cgtex
submit_r2_test F1 config/cattle_sel_bottlenecked.yaml        cattle_sel_bot                 61 0.01 False cgwas cgtex
submit_r2_test G1 config/cattle_sel_not_bottlenecked.yaml    cattle_sel_nobot               71 0.01 False cgwas cgtex
EOF
```

After all 10 DAP-G jobs complete, inspect:

```bash
ssh o2 'sacct -u $USER --starttime=$(date -d "1 day ago" +%Y-%m-%d) \
        --format=JobID,JobName%40,State,Elapsed,MaxRSS,ReqMem' | grep stage3_dapg_locus
```

Then bump `_dapg_mem_mb_base` (in `simulations/rules/common.smk` lines ~152-174)
and the `runtime` lambda (line ~209) to `ceil(1.5 * MaxRSS)` /
`ceil(1.5 * Elapsed)` per category. If any test exceeds ~9h (1.5 x 6h
current ceiling), escalate the gwas partition from `short` to `medium` for
the production run.

## Production submission (20 reps: A1-4, B1-4, E1-4, F1-4, G1-4)

After resources are tuned, submit the full r2=0.25 set. Targeting `stage4`
(not `all`) prevents the stage-5 plink GLM rules from being scheduled.

```bash
ssh o2.hms.harvard.edu << 'EOF'
set -e
SIMS_PUB=/n/data2/hms/dbmi/sunyaev/lab/nconnally/simulations_for_revision
SIMS_SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_for_revision
SUBMIT_DIR=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SNAKEMAKE=/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake

# ID CFG BASE SEED GWAS_SIZE GTEX_SIZE QV NEUTRAL EXTRA
submit_r2() {
  local ID=$1 CFG=$2 BASE=$3 SEED=$4 GWAS=$5 GTEX=$6 QV=$7 NEUTRAL=$8 EXTRA=$9
  sbatch --parsable \
    --job-name=snake_${ID}_r2_025_sd${SEED} \
    --partition=long --time=7-00:00:00 --mem=4G --cpus-per-task=1 \
    --output=$SIMS_SCRATCH/$ID/logs/controller_r2_025_%j.out \
    --error=$SIMS_SCRATCH/$ID/logs/controller_r2_025_%j.err \
    --chdir=$SUBMIT_DIR \
    --wrap="$SNAKEMAKE --snakefile $SUBMIT_DIR/Snakefile \
      --directory $SIMS_SCRATCH/$ID \
      --configfile $SUBMIT_DIR/$CFG \
      --config basename=${ID}_${BASE} stage1_seed=$SEED stage2_seed=$SEED \
               gwas_size=$GWAS gtex_size=$GTEX Q_scaling=$QV neutral_trait_vars=$NEUTRAL \
               ld_ctrl=0.25 output_tag=r2_0_25 \
               workdir=$SIMS_SCRATCH/$ID publishdir=$SIMS_PUB/$ID $EXTRA \
      --profile $SUBMIT_DIR/profiles/o2 -j 200 --keep-going --rerun-incomplete \
      stage4"
}

# A: human, 8k GWAS, multi-GTEx, directional negative selection
for i in 1 2 3 4; do submit_r2 A$i config/human.yaml                         human                         $((10+i)) 8000 -1 10   False ""; done
# B: human, 8k GWAS, multi-GTEx, neutral trait variants
for i in 1 2 3 4; do submit_r2 B$i config/human.yaml                         human                         $((20+i)) 8000 -1 10   True  ""; done
# E: cattle baseline-from-midpoint, multi-GTEx, directional negative
for i in 1 2 3 4; do submit_r2 E$i config/cattle_baseline_from_midpoint.yaml cattle_baseline_from_midpoint $((50+i)) 8000 -1 0.01 False ""; done
# F: cattle, bottlenecked, directional positive selection
for i in 1 2 3 4; do submit_r2 F$i config/cattle_sel_bottlenecked.yaml       cattle_sel_bot                $((60+i)) 8000 -1 0.01 False ""; done
# G: cattle, non-bottlenecked, directional positive selection
for i in 1 2 3 4; do submit_r2 G$i config/cattle_sel_not_bottlenecked.yaml   cattle_sel_nobot              $((70+i)) 8000 -1 0.01 False ""; done
EOF
```

Sanity checks after a rep completes:

```bash
# Outputs match the manifest, one per trait.
ls $SIMS_SCRATCH/A1/stage3/$RUN_TAG/hgwas/outputs_r2_0_25/ | wc -l
wc -l $SIMS_SCRATCH/A1/stage3/$RUN_TAG/hgwas/manifest.txt

# Stage-4 result set exists.
ls $SIMS_SCRATCH/A1/stage4/$RUN_TAG/A1_human.r2_0_25.*.enloc.*.out

# Existing r2=0.75 outputs are unmodified (mtime check).
stat -c '%y %n' $SIMS_SCRATCH/A1/stage3/$RUN_TAG/hgwas/outputs/*.dapg.out | head
```


