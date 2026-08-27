# Simulation comparison dashboard

Side-by-side comparison of the simulation arms in `simulation_data/`, in one
static HTML page. Open `simulation_dashboard.html` in any browser — no server,
no build step.

```
open figures_and_tables/dashboard/simulation_dashboard.html
```

## What is in it

- **Left pane**: every simulation with its causative-variant MAF
  cutoff, its fine-mapping/GLM cutoff, and TPR / FPR / FDR at RCP > 0.5 and
  > 0.9. Sortable by any column, filterable, and clicking a row loads that
  simulation into the highlighted comparison box. The **simulation table**
  toggle in the header (or the pane's own *hide* button) collapses it and gives
  the width to the comparison boxes.
- **Comparison boxes**: pick how many (1–8). Each box selects a simulation by
  category, causative cutoff, fine-mapping/GLM cutoff and parameter arm.
- **Plot selector**: one universal set of toggles chooses which plots appear in
  every box — the plots the report PDFs carry, plus the parameter and results
  tables.
- **Table-row selector**: a second universal set of toggles chooses which rows
  the results table shows in every box — power, false positives and FDR at each
  cutoff, each of those again restricted to GWAS-significant loci, and the two
  denominators.
- **Outcome comparison plot** (toggle at the end of the Plots row): one chart
  above the boxes with a group per selected percentage row and, inside each
  group, one bar per compared simulation, coloured and numbered to match the
  boxes. It plots exactly the rows ticked under *Table rows*, minus the two
  denominator rows, which have no meaning on a percent axis. Where two boxes
  share a category colour the duplicate is stepped lighter or darker within the
  same hue until it separates, so the category family still reads correctly.

Rate definitions: TPR and FPR are out of the GWAS traits scored (power and
false-positive rate as the report tables define them), while **FDR is
FP / (TP + FP)** — out of the loci actually called at that cutoff, i.e.
1 − precision. Rates are shown as `–` when the denominator is zero. Because the
reduction keeps one call per GWAS trait and prefers a correct pairing over a
higher-RCP wrong one, FDR is a lower bound on the rate of spurious
above-cutoff pairings.

**Significant loci.** The `sig. loci` rows, and the left table's `all loci` /
`sig. loci` switch, restrict the whole tally to the loci a GWAS of the selected
size would find: a true positive counts only if its GWAS trait is significant,
a false positive likewise, and the TPR/FPR denominator becomes the
significant-locus count (the `n` column). Significance is the expected-power
definition the report PDFs use — the causal variant's non-centrality parameter
at the selected study size clears χ² = 29.72 — not the observed lead-SNP
p-value of the simulated study.

A *simulation* here is one simulation category (A–N) inside one parameter arm,
pooled over that arm's replicates — the same unit the category report PDFs
compare.

Global controls (eQTL panel size, RCP cutoff, GWAS size for "significant
loci", allele-frequency MAF floor and y-axis) scope every box at once, so the
boxes are always showing the same thing about different simulations. The
settings persist in `localStorage`, and a view can be linked:

```
simulation_dashboard.html?sims=cmaf001_fm01_g5t20_2Mb::A,cmaf001_fm01_g5t20_2Mb::E&plots=all&rcp=0.9
```

Recognised parameters: `sims`, `n`, `plots` (`all`, `none`, or ids), `rows`
(`all`, `none`, or ids), `tloci` (`all`/`sig`), `table` (`0` hides the left
pane), `compare` (`0` hides the outcome comparison plot), `panel`, `rcp`,
`gwasn`, `floor`, `afy`, `sharedy`.

## Rebuilding the data

`dashboard_data.js` is generated. The full refresh cycle, after new replicates
finish on O2:

```
o2-ssh 'bash /n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake/archive_round3_to_data2.sh'
bash simulations/fetch_round3_2Mb.sh                     # mirror into simulation_data/
Rscript figures_and_tables/dashboard/build_dashboard_data.R
```

The archive step copies only replicates whose four stages have all landed, and
the mirror is a single rsync of the whole archive, so both are safe to re-run.
The builder walks every loadable replicate under `simulation_data/` (~10 s for
207 replicates across 22 arms → 55 simulations) and writes ~0.8 MB of JSON.
Nothing else needs editing when a new arm or category appears: the arm list,
the category dropdown, the colours and the labels are all read from the data
and from the notebook's own `CATEGORY_LABEL` / `color_key`.

A new arm's `min_maf` (the plink2 GLM floor) is the one setting no filename
encodes — it comes only from a hand-written `simulation_data/<arm>/PROVENANCE.txt`,
which is also where the arm's purpose should be recorded. Without one the
dashboard shows `n/r` in the GLM column.

The builder **sources the helper cells of `figure2_revision2.ipynb` directly**
(setup, filename grammar, `select_runs`, the readers, `load_sims`, the
GWAS-significance calculation, the outcome definitions), so the dashboard and
the report PDFs cannot drift apart. What it exports is deliberately upstream of
any plot:

- **traits** — the `enloc_best_rows()` reduction: one best row per GWAS trait
  per replicate per eQTL panel, carrying the MAF bin, the RCP, the three-valued
  `correct`, and `vexp/(1+vexp)` so the browser can recompute GWAS significance
  at any study size. Nothing here depends on an RCP cutoff, which is why the
  cutoff can be a live control.
- **af** — allele-frequency histograms (50 bins of width 0.01, so the 1% floor
  lands on a bin edge) and MAF-by-selection-coefficient box statistics, pooled
  over replicates.

Two notes on what the builder had to handle:

- The two psamp arms relabelled `cmaf0` (10x and 30x) carry both namings of
  every variant table (`_maf_0.01.tsv` and `_maf_0.tsv`) because the relabel
  copied rather than moved. All 88 pairs are byte-identical, so `find_file()`
  is overridden to resolve identical duplicates to one file instead of
  stopping. **The notebook itself still stops on these**, so those 22
  replicates cannot currently be loaded with `load_sims()`.
- Array fields are wrapped in `I()`: `toJSON(auto_unbox = TRUE)` would
  otherwise write a single-replicate arm's `reps` as a bare string.
