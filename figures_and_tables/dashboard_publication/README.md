# Publication coloc matrix

Comparison dashboard for the frozen publication simulation set in
`simulation_data/farm_sims_for_publication/`. One static HTML file, no server and
no build step:

```
open figures_and_tables/dashboard_publication/publication_dashboard.html
```

This is a **separate** dashboard from `figures_and_tables/dashboard/`, not a
replacement. That one covers the round-3 exploratory arms in `simulation_data/`,
whose directory layout (`<arm>/<LETTER><N>/`) and naming the publication tree does
not share. Keeping both means the old figures stay reproducible.

## What it shows

A *cell* is one simulation category inside one causal-sampling strategy, pooled
over that strategy's finished replicates — the same unit the older dashboard
compares.

- **Filters** — checkboxes for each of the six categories and each of the four
  sampling strategies, with all/none shortcuts. They scope every view below.
- **Category × strategy matrix** — categories across, sampling strategies down.
  Each cell is two-sided about a shared zero: **power** grows upward to 100%,
  **FDR** downward to 50%, both drawn at the same pixels per percentage point
  (0.88 px/%), so a bar of a given length means the same on either side. The FDR
  track is therefore half the height of the power track. 50% is the ceiling
  because the largest FDR anywhere in the set is 35.3% (cattle baseline,
  500-person panel, RCP > 0.5), so nothing clips.
- **Allele frequency by selection strength** — vertical box-and-whisker of MAF
  against selection-coefficient bin, in the same format as
  `figures_and_tables/dashboard`: box is the IQR on the surface fill, heavy line
  is the median, whiskers reach the furthest point within 1.5×IQR and the tooltip
  counts what lies beyond. One panel per category × sampling strategy.
  A **dashed line** per selected set marks that set's panel-wide median, pooled
  across selection bins (stored under `sb = -1`; a median cannot be recovered
  from per-bin medians, so it needs its own histogram).
  **Variant-set checkboxes** choose which series are drawn:
  `all` · `non-neutral` (selco ≠ 0, the pool causal loci are drawn from) ·
  `neutral` (selco = 0) · `causative` (actually assigned a trait effect).
  Under the synthetic DFE the background-selection categories draw their
  causative variants from the *neutral* pool, so that series lands in the `0`
  bin there — visible directly in the plot.
- **Outcome by causative-variant frequency** — vertical stacked counts of GWAS
  traits against the MAF bin of the causative variant (the variant assigned that
  trait's effect), stacked bottom-to-top as true positive, false positive,
  underpowered colocalization, underpowered fine-mapping. One panel per category × sampling strategy, with power / FP / FDR in
  the corner. Four outcome classes, matching the old dashboard's `outcomeOf`:
  `tp` correct pairing above the cutoff · `fp` only a wrong pairing above it ·
  `uc` underpowered colocalization, a signal at the right pairing but under the
  cutoff · `ufm` underpowered fine-mapping, no signal at the right pairing at
  all. (The old dashboard defines `ufm` as "no fastEnloc row at all"; here it is
  "no signal at the correct pairing", which is what the trait dump can
  distinguish.)

Both grids give **each panel its own y-axis** by default. A shared scale
compressed most panels into the bottom of the frame — measured, the per-panel box
maxima span 0.085 to 0.5 and the stack maxima 51 to 307. The **shared y-axis**
control turns sharing back on when panels need to be compared directly.

### Colours

Taken from the project's own `color_key` in `figures_and_tables/figure2_revision2.ipynb`,
so the dashboard and the report figures agree: **human coloc** is `firebrick4`
(`#8B1A1A`) and **cattle coloc** is `purple4` (`#551A8B`). Power bars and
true-positive segments use the species colour; FDR and false positives use the
maximum-contrast slot, underpowered colocalization dark grey, underpowered
fine-mapping light grey.

That contrast slot is black in light mode and near-white in dark mode. Pure black
on a dark ground is invisible, and what the slot encodes is "maximum contrast
against the panel", not the literal colour.

### Terminology

Following the pipeline's own usage: a **causal variant** is the variant that
causes a trait; a **causative** locus or variant is one actually assigned an
effect (`summarize_coloc`'s `*_cs_causative` is "at the causative loci only");
**non-neutral** means `selco ≠ 0`, the pool from which causal loci are drawn.


Global controls (eQTL panel size, RCP cutoff) scope everything at once.

## Definitions

**Power** = GWAS loci colocalizing with the *correct* eQTL trait, over every
defined GWAS locus. Loci with no eQTL partner stay in the denominator and can
never contribute a hit — that is deliberate (`summarize_coloc.py:843`), and the
drag it produces is about 0.45 points on the unpaired arms.

**FDR** = FP / (TP + FP), where a false positive is a locus that colocalized only
with the wrong trait. TP and FP are mutually exclusive; a locus clearing nothing
is neither.

## Rebuilding

Data is embedded in the HTML, so the page is self-contained. To refresh it after
more replicates land:

1. On O2, filter `RUNS.tsv` to runs passing the publish gate (2 `enloc.sig.out` +
   3 `_glm_lead_snps.tsv`) and drop it at the scratch root. **This step is not
   optional**: `summarize_coloc.py` scores whatever files exist, so an in-flight
   run is not a smaller sample but a wrong one — its GWAS loci land in the
   denominator while its missing stage-4 output contributes no numerator, and the
   cell reads as low power rather than as incomplete.
2. Run `summarize_coloc.py` with `--dump-traits`.
3. Rebuild the embedded JSON from the trait dump (per arm × category × panel ×
   replicate: trait count, TP and FP at RCP 0.5 and 0.9) and reassemble.

The build was verified by re-deriving all 96 cell × threshold combinations from
the embedded per-replicate data and checking them against `summarize_coloc.py`'s
own summary — 0 mismatches. Worth repeating on any rebuild.

### MAF histogram resolution

The aggregate bins MAF at **0.0001**, not 0.001. With 9,000 GWAS individuals the
smallest non-zero MAF is ~5.6e-5 and the human panel is dominated by
near-singleton variants: at 0.001 resolution more than half of every human bin
fell into bucket 0 and every human median, quartile and whisker reported exactly
0.0000. Cattle was unaffected (median 0.119), which is why it went unnoticed —
it only became visible when the median line was added.

For scale, human: panel median 0.0001, causative median 0.0056. The causative
variants are ~56× commoner than the panel at large, which is `causal_min_maf =
0.001` gating them. A linear MAF axis will always compress the human panels;
reading them well would want a log axis.

## Checking a rebuild

`node --check` only validates syntax. It does **not** catch a function that is
called but never defined — which happened twice: a slice-and-replace edit removed
`renderSpread` and `renderBlocks`, `renderAll` threw partway, and two panels
rendered blank while everything above them looked fine.

    node render_smoke.js publication_dashboard.fragment.html

executes the page's real logic against a stub DOM and asserts every target
element actually received content. Run it after any edit to the script block.

## Files

- `publication_dashboard.html` — standalone, open directly in a browser.
- `render_smoke.js` — executes the dashboard logic against a stub DOM and fails
  if any panel is never painted.
- `build_maf_aggregates.py` — runs on O2; produces the MAF/selection and
  outcome/MAF aggregates the last two plots embed. Reads ~11M variant rows and
  emits quantiles from a 0.001-resolution histogram rather than holding every
  value.
- `publication_dashboard.fragment.html` — same page without the document
  wrapper, which is the form the Artifact publisher takes.
