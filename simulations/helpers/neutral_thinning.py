"""Thin the NEUTRAL class of a tree sequence -- categories O and P.

Why this is a deletion and not a lower mutation rate
----------------------------------------------------
Neutral mutations are not part of either species' forward simulation. They are a
Poisson process overlaid on branches the forward run has already fixed:

    human  ``human_simulation_o2.py:203``  msprime.sim_mutations(rate=8.4e-9) in
                                           STAGE 1, so they arrive in ``hts_*.ts``
    cattle ``create_gwas_files_and_phenotypes.add_neutral``  the same 8.4e-9 (in
                                           two passes across the Q handoff) in
                                           STAGE 2

Dropping each neutral site independently with probability ``1 - k`` is the
thinning of a Poisson process, and a thinned Poisson process is a Poisson process
at the thinned rate. So this is not an approximation of "simulate at
``k * 8.4e-9``" -- it has exactly that distribution, conditional on the same
genealogy. Three things follow, and they are the reason categories O and P are
built this way:

* One code path serves both species, rather than editing a stage-1 literal for
  human and a stage-2 literal for cattle.
* Stage 1 is REUSED. O_i adopts A_i's tree sequence and P_i adopts E_i's, at A's
  and E's own seeds, so O_i's variant set is a strict SUBSET of A_i's and the
  contrast carries no replicate noise.
* The causal machinery is untouched by construction. ``causal_eligible``,
  ``flank_eligible``, ``select_central_power`` and ``select_gtex_topup`` all
  select on ``selco != 0``; this module only ever deletes ``selco == 0`` sites.
  Under ``synthetic_dfe_effects`` that predicate is INVERTED (the neutral class
  is the causal pool), which is why the Snakefile refuses to combine the two.
  The pools, the pi-PS weights, the drawn positions, the effect sizes and the
  trait-partner tables therefore come out IDENTICAL to the un-thinned arm.

What it is NOT safe to combine with
-----------------------------------
``synthetic_dfe_effects`` (categories H, I, K, L) draws the causal loci from the
strictly neutral variants, and ``neutral_trait_vars`` (B, D) draws its effect
RECIPIENTS from them. For those categories the neutral class IS the causal pool,
so thinning it would silently halve the one thing the category is defined by.
The Snakefile refuses the combination rather than running it.

Sizing the cut
--------------
pi is additive over sites, so the keep fraction that hits a target pi is solved,
not searched -- see ``keep_fraction_for_pi_target``.
"""

import numpy as np

__all__ = [
    "NO_NEUTRAL_THINNING",
    "neutral_site_ids",
    "thin_neutral_sites",
    "pi_components",
    "keep_fraction_for_pi_target",
]

#: Keep everything. Mirrors ``helpers.paths.NO_NEUTRAL_THINNING``; duplicated as a
#: plain float so this module has no import cycle back into the path builders.
NO_NEUTRAL_THINNING = 1.0

#: Second word of the thinning RNG's seed. The draw gets its own
#: ``default_rng([seed, _THIN_STREAM])`` rather than consuming from any existing
#: generator, so adding thinning to a run cannot shift the trait subsample, the
#: pi-PS draw or the nucleotide overlay -- all of which build their own
#: ``default_rng(seed)`` from the same stage-2 seed.
_THIN_STREAM = 0x6E7468  # "nth"


def _site_selection_coeff(site):
    """The selection coefficient of a site, read exactly as stage 2 reads it.

    ``create_gwas_files_and_phenotypes.get_vars_df`` takes
    ``v.site.mutations[0].metadata['mutation_list'][0]['selection_coeff']`` and
    ``causal_eligible`` tests the result against 0. Using the same accessor here
    is deliberate: the two can then never disagree about which sites are neutral,
    which is what makes "thinning never removes a causal candidate" a fact rather
    than a hope. (Sites are single-mutation here -- ``remove_fixed`` drops
    multiallelics upstream in both species.)

    The recorded coefficient is Q-scaled, but the test is against 0 and 0 is
    scale-free, so no ``Q_scaling`` is needed.
    """
    return site.mutations[0].metadata["mutation_list"][0]["selection_coeff"]


def neutral_site_ids(ts):
    """Ids of the sites carrying no selection coefficient, as an int array."""
    return np.array(
        [s.id for s in ts.sites() if _site_selection_coeff(s) == 0],
        dtype=np.int64,
    )


def selected_site_ids(ts):
    """Ids of the sites carrying a non-zero selection coefficient."""
    return np.array(
        [s.id for s in ts.sites() if _site_selection_coeff(s) != 0],
        dtype=np.int64,
    )


def thin_neutral_sites(ts, keep_fraction, seed, measure_pi=True):
    """Delete a seeded random ``1 - keep_fraction`` of the NEUTRAL sites.

    Returns ``(thinned_ts, counts)`` where ``counts`` records what happened, for
    the stage-2 parameters sidecar. ``keep_fraction == 1.0`` returns the input
    unchanged with an EMPTY ``counts``, and does no work at all -- no site scan,
    no pi -- so a run that does not ask for thinning is bit-identical to one built
    before this module existed, provenance sidecar included.

    The number kept is ``round(keep_fraction * n_neutral)`` exactly, not a
    per-site coin flip. Both are valid thinnings; the fixed count removes the
    binomial noise between replicates, which matters because the whole point of
    the arm is that O_i differs from A_i in one controlled way.

    ``measure_pi`` adds ``pi_before`` / ``pi_after`` / ``pi_neutral_before`` /
    ``pi_selected`` to ``counts``. That is three ``diversity()`` passes, which is
    why it only happens when thinning is actually on: it is the number the arm is
    calibrated against, so it belongs beside the outputs rather than in a
    separate audit.
    """
    keep_fraction = float(keep_fraction)
    if not (0.0 < keep_fraction <= 1.0):
        raise ValueError(
            f"neutral_keep_fraction must be in (0, 1], got {keep_fraction!r}"
        )

    if keep_fraction == NO_NEUTRAL_THINNING:
        # Return BEFORE the site scan, not after. An un-thinned arm must pay
        # nothing for this key existing -- neither the metadata decode over every
        # site nor a set of extra rows in its stage2_params.txt sidecar, which
        # would make its provenance record differ from the ones already on disk.
        return ts, {}

    neutral = neutral_site_ids(ts)
    counts = {
        "sites_before": int(ts.num_sites),
        "neutral_sites_before": int(neutral.size),
        "selected_sites": int(ts.num_sites - neutral.size),
        "neutral_keep_fraction": keep_fraction,
    }

    if measure_pi:
        pi_neutral, pi_selected = pi_components(ts)
        counts["pi_neutral_before"] = pi_neutral
        counts["pi_selected"] = pi_selected
        counts["pi_before"] = pi_neutral + pi_selected

    n_keep = int(round(keep_fraction * neutral.size))
    rng = np.random.default_rng([int(seed), _THIN_STREAM])
    drop = rng.choice(neutral, size=neutral.size - n_keep, replace=False)
    thinned = ts.delete_sites(np.sort(drop))

    counts["neutral_sites_kept"] = n_keep
    counts["sites_after"] = int(thinned.num_sites)
    if measure_pi:
        counts["pi_after"] = float(thinned.diversity())
        counts["pi_ratio"] = (counts["pi_after"] / counts["pi_before"]
                              if counts["pi_before"] else float("nan"))
    return thinned, counts


def pi_components(ts):
    """Nucleotide diversity split into its neutral and selected halves.

    Returns ``(pi_neutral, pi_selected)``, both per base pair over the whole
    sequence, so ``pi_neutral + pi_selected`` is ``ts.diversity()``. Computed by
    deleting the other class and asking tskit, which keeps the estimator (and its
    sample-size correction) identical to the one every other pi number in this
    project is quoted from.
    """
    neutral = neutral_site_ids(ts)
    selected = selected_site_ids(ts)
    pi_neutral = float(ts.delete_sites(selected).diversity())
    pi_selected = float(ts.delete_sites(neutral).diversity())
    return pi_neutral, pi_selected


def keep_fraction_for_pi_target(pi_neutral, pi_selected, target=0.5):
    """The keep fraction that scales pi to ``target`` times its current value.

    pi is a sum over sites, and thinning removes neutral sites without touching
    any surviving site's frequency, so

        E[pi(k)] = k * pi_neutral + pi_selected

    and setting that equal to ``target * (pi_neutral + pi_selected)`` gives

        k = target + (target - 1) * pi_selected / pi_neutral

    which at ``target = 0.5`` is ``0.5 * (1 - pi_selected / pi_neutral)``. No
    search, no calibration loop: one measurement of the two components fixes k.

    Raises when the target is unreachable -- if the selected class already
    carries more than ``target`` of the total, no amount of neutral thinning can
    get there, and clamping to a positive k would quietly deliver a different
    experiment than the one asked for.
    """
    target = float(target)
    if not (0.0 < target < 1.0):
        raise ValueError(f"target must be in (0, 1), got {target!r}")
    if pi_neutral <= 0:
        raise ValueError(
            f"pi_neutral must be positive to thin against it, got {pi_neutral!r}"
        )

    k = target + (target - 1.0) * (pi_selected / pi_neutral)
    if k <= 0:
        total = pi_neutral + pi_selected
        raise ValueError(
            f"pi cannot be scaled to {target:g}x by thinning the neutral class: "
            f"the SELECTED sites alone carry pi={pi_selected:.6g} of a total "
            f"{total:.6g} ({pi_selected / total:.1%}), which is already at or "
            f"above the target. Deleting every neutral site would still leave "
            f"pi above {target:g}x."
        )
    return k
