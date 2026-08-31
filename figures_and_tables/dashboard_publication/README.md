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
- **Category × strategy matrix** — each cell is two-sided about a shared zero:
  **power** grows rightward, **FDR** grows leftward, both on a 0–100% scale, so
  a cell that is good on both counts reads as a short left arm and a long right
  one. The `n` column is the replicate count, highlighted when a row's
  strategies pool different numbers.
- **Replicate spread** — grouped by category, then by a labelled lane per
  sampling strategy. One mark per replicate with the mean, so a difference can
  be read against the scatter it sits in.
- **Cattle ancestry blocks** — power per deep history. Cattle replicates come in
  blocks of five that each resume from an independently simulated ancestry;
  human replicates are each their own population, so their blocks are an
  arbitrary split and act as a noise reference.

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

## Files

- `publication_dashboard.html` — standalone, open directly in a browser.
- `publication_dashboard.fragment.html` — same page without the document
  wrapper, which is the form the Artifact publisher takes.
