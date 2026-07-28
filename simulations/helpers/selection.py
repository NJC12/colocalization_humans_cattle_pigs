"""Trait selection for category-wide fine-mapping caps.

A category spans several replicate Snakemake runs that differ only in
`stage1_seed`. To bound the cost of fine-mapping + fastEnloc across the
category, each replicate independently selects a quota of traits from its
own pheno file. The quota and per-replicate RNG are derived purely from
config fields, so no cross-replicate coordination is needed.
"""

import math

import numpy as np


def select_traits_for_replicate(all_traits, cap, seed, num_replicates, stage1_seed):
    """Deterministic per-replicate trait selection.

    The per-replicate RNG is seeded from (seed, stage1_seed) so different
    replicates draw distinct permutations and the choice is reproducible.
    Cap extension is automatic: raising `cap` lengthens the prefix without
    reshuffling, so previously-selected traits remain selected.

    If quota exceeds the number of available traits, returns the full
    permutation without error.
    """
    quota = math.ceil(cap / num_replicates)
    quota = min(quota, len(all_traits))
    ss = np.random.SeedSequence([int(seed), int(stage1_seed)])
    rng = np.random.default_rng(ss)
    perm = list(all_traits)
    rng.shuffle(perm)
    return perm[:quota]
