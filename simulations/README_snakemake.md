# Snakemake-driven simulation pipelines

Four pipelines, one Snakefile, four config YAMLs.

## Pipelines

| Config                                   | Pipeline               | Stage 1 producer                                                  |
|------------------------------------------|------------------------|-------------------------------------------------------------------|
| `config/human.yaml`                      | human                  | `human_simulation_o2.py` (or `_larger.py`)                        |
| `config/cattle_baseline.yaml`            | cattle_baseline        | `farm_create_orig_pop_e2.py` -> `farm_burn_in_e2.slim` -> `farm_selection.slim` |
| `config/cattle_sel_bottlenecked.yaml`    | cattle_sel             | `farm_selection_from_ep8.slim` (`num_muts_selected>0`, `continue_bottlenecking=1`) |
| `config/cattle_sel_not_bottlenecked.yaml`| cattle_sel             | `farm_selection_from_ep8.slim` (`num_muts_selected>0`, `continue_bottlenecking=0`) |

All four share stages 2 (split + traits), 3 (DAP-G), 4 (fastEnloc), 5 (plink GLM).

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
| Human      |     8,000 | 1,000, 500, 250 | Neutral                                | NA             | B  |                         human |          4 |
| Human      |   100,000 |          50,000 | Directional (negative)                 | NA             | C  |                         human |          1 |
| Human      |   100,000 |          50,000 | Neutral                                | NA             | D  |                         human |          1 |
| Cattle     |     8,000 | 1,000, 500, 250 | Directional (negative)                 | yes            | E  | cattle_baseline_from_midpoint |          4 |
| Cattle     |     8,000 | 1,000, 500, 250 | Directional (positive)                 | yes            | F  |       cattle_sel_bottlenecked |          4 |
| Cattle     |     8,000 | 1,000, 500, 250 | Directional (positive)                 | no             | G  |   cattle_sel_not_bottlenecked |          4 |

## List of snakemake submission commands

The block below submits all 22 controller jobs in one shell session.
SLURM dependencies (`--dependency=afterany:...`) hold C and D until
every other controller has exited, so the original table order is
preserved here even though the two heavyweight runs actually execute
last. Each run gets its own flat publishdir
(`/n/data2/.../simulations_for_revision/<ID><rep>/`) and matching
scratch workdir, so the 22 runs never collide and existing in-flight
controllers (e.g. `snakemake/cattle_sel_bot/`) keep running undisturbed.

Seeds encode the run identity: tens digit = letter (A=1, B=2, …, G=7),
ones digit = replicate. So A1 → 11; D1 → 41; F2 → 62; G4 → 74.

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


