"""Categories O and P: A's and E's genomes with pi halved, and nothing else.

Three things here can go wrong quietly, and each has a test.

1. **The thinning could touch a selected site.** The whole claim of the arm is
   that the causal architecture is untouched -- same 50 central loci, same flank
   loci, same effect sizes -- and that rests entirely on the deletion set being a
   subset of ``selco == 0``. If ``neutral_site_ids`` ever read a different field
   than ``get_vars_df`` does, thinning would start removing causal candidates and
   the arm would still run, still produce plausible numbers, and no longer be A
   minus background.

2. **The keep fraction could go missing.** An unset ``neutral_keep_fraction``
   means keep-everything, so O would be a byte-identical copy of A under another
   name -- a null result that looks like a measurement.

3. **The seeds could be "fixed".** Every other category has a private seed band;
   O and P deliberately share A's and E's, because that is what makes O_i's
   variant set a strict subset of A_i's. A well-meaning edit to give them their
   own band would silently turn a paired contrast into an unpaired one.

The tree-sequence tests use a fake rather than tskit, following
test_causal_selection.py: this machine may not carry the SLiM Python stack, and
what can be wrong here is which sites get chosen, not whether tskit can delete
them.
"""

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from helpers import neutral_thinning as nt
from helpers import params_record, paths


SIM_DIR = Path(paths.__file__).parent.parent
CONFIG_DIR = SIM_DIR / "config"
SRC = (SIM_DIR / "create_gwas_files_and_phenotypes.py").read_text()

HUMAN_LOWHET = sorted(CONFIG_DIR.glob("human_lowhet_*.yaml"))
CATTLE_LOWHET = sorted(CONFIG_DIR.glob("cattle_lowhet_*.yaml"))


def load_yaml(path):
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh)


# --------------------------------------------------------------- the fake ts

class _Mutation:
    def __init__(self, selco, mutation_type, time=0.0):
        self.metadata = {"mutation_list": [{"selection_coeff": float(selco),
                                            "mutation_type": mutation_type}]}
        self.time = time


class _Site:
    def __init__(self, sid, position, selco):
        self.id = sid
        self.position = float(position)
        # type 0 is msprime's neutral overlay, 1-3 the DFE components. The class
        # test is on the COEFFICIENT, never on the type: the human overlay is
        # written with SLiMMutationModel(type=2), which collides with DFE
        # component m2, so the type does not separate the two classes and the
        # fake reproduces that collision on purpose.
        self.mutation_type = 0 if selco == 0 else 2
        self.mutations = [_Mutation(selco, self.mutation_type)]


class _Variant:
    def __init__(self, site, genotypes):
        self.site = site
        self.genotypes = np.asarray(genotypes)


class FakeTS:
    """Enough of a tskit.TreeSequence for the thinning code path."""

    def __init__(self, selcos, dafs=None, sequence_length=2e6, n_hap=100):
        self.num_samples = n_hap
        self.sequence_length = sequence_length
        selcos = list(selcos)
        if dafs is None:
            dafs = [0.25] * len(selcos)
        self._sites = [_Site(i, 1000 * (i + 1), s) for i, s in enumerate(selcos)]
        self._dafs = list(dafs)

    @property
    def num_sites(self):
        return len(self._sites)

    def sites(self):
        return iter(self._sites)

    def variants(self):
        for site, daf in zip(self._sites, self._dafs):
            n1 = int(round(daf * self.num_samples))
            g = np.zeros(self.num_samples, dtype=int)
            g[:n1] = 1
            yield _Variant(site, g)

    def delete_sites(self, site_ids):
        drop = set(int(i) for i in np.asarray(site_ids, dtype=np.int64).ravel())
        keep = [(s, d) for s, d in zip(self._sites, self._dafs) if s.id not in drop]
        out = FakeTS([s.mutations[0].metadata["mutation_list"][0]["selection_coeff"]
                      for s, _ in keep],
                     dafs=[d for _, d in keep],
                     sequence_length=self.sequence_length,
                     n_hap=self.num_samples)
        # Site ids are re-numbered by delete_sites, exactly as tskit does; keep
        # the positions so a caller can still identify what survived.
        for new, (old, _) in zip(out._sites, keep):
            new.position = old.position
        return out

    def diversity(self):
        return sum(2.0 * d * (1.0 - d) for d in self._dafs) / self.sequence_length

    def positions(self):
        return [s.position for s in self._sites]


def mixed_ts(n_neutral=200, n_selected=100, seed=3):
    """A forward-genome shape: a common neutral majority, a rarer selected tail."""
    rng = np.random.default_rng(seed)
    selcos = [0.0] * n_neutral + list(-rng.uniform(1e-4, 1e-2, n_selected))
    dafs = list(rng.uniform(0.05, 0.5, n_neutral)) + list(rng.uniform(0.001, 0.05, n_selected))
    order = rng.permutation(len(selcos))
    return FakeTS([selcos[i] for i in order], dafs=[dafs[i] for i in order])


# ------------------------------------------------ the deletion set is neutral

def test_only_neutral_sites_are_ever_deleted():
    """The claim the whole arm rests on. Selected sites are the causal pool."""
    ts = mixed_ts()
    before = {s.position for s in ts.sites()
              if s.mutations[0].metadata["mutation_list"][0]["selection_coeff"] != 0}
    thinned, _ = nt.thin_neutral_sites(ts, 0.4, seed=11, measure_pi=False)
    after = {s.position for s in thinned.sites()
             if s.mutations[0].metadata["mutation_list"][0]["selection_coeff"] != 0}
    assert before == after, "thinning removed a selected site"


def test_the_survivors_are_a_strict_subset():
    """O_i's variant set must be a SUBSET of A_i's, not merely the same size."""
    ts = mixed_ts()
    thinned, _ = nt.thin_neutral_sites(ts, 0.4, seed=11, measure_pi=False)
    assert set(thinned.positions()) < set(ts.positions())


def test_the_kept_count_is_exact_and_not_binomial():
    """A fixed count, so O_i differs from A_i in one controlled way rather than
    in one controlled way plus a draw."""
    ts = mixed_ts(n_neutral=200, n_selected=100)
    thinned, counts = nt.thin_neutral_sites(ts, 0.375, seed=11, measure_pi=False)
    assert counts["neutral_sites_before"] == 200
    assert counts["neutral_sites_kept"] == 75          # round(0.375 * 200)
    assert counts["selected_sites"] == 100
    assert thinned.num_sites == 175


def test_the_draw_is_seeded():
    ts = mixed_ts()
    a, _ = nt.thin_neutral_sites(ts, 0.4, seed=11, measure_pi=False)
    b, _ = nt.thin_neutral_sites(ts, 0.4, seed=11, measure_pi=False)
    c, _ = nt.thin_neutral_sites(ts, 0.4, seed=12, measure_pi=False)
    assert a.positions() == b.positions()
    assert a.positions() != c.positions()


def test_keeping_everything_does_nothing_at_all():
    """An arm that does not ask for thinning must be bit-identical to one built
    before this key existed -- including its stage2_params.txt sidecar, which is
    why counts comes back EMPTY rather than merely unchanged."""
    ts = mixed_ts()
    out, counts = nt.thin_neutral_sites(ts, 1.0, seed=11)
    assert out is ts
    assert counts == {}


@pytest.mark.parametrize("bad", [0.0, -0.1, 1.5])
def test_an_out_of_range_fraction_is_refused(bad):
    with pytest.raises(ValueError):
        nt.thin_neutral_sites(mixed_ts(), bad, seed=11)


# ------------------------------- the neutral predicate matches the pipeline's

def _load_from_source(fn):
    """Exec one top-level function out of the stage-2 script, as
    test_background_selection.py does -- importing it would pull in tskit."""
    m = re.search(rf'^def {fn}\(.*?\n(?=^def |\Z)', SRC, re.S | re.M)
    assert m, f"{fn} not found"
    ns = {'np': np, 'pd': pd}
    exec(compile(m.group(0), fn, 'exec'), ns)
    return ns[fn]


def test_neutral_site_ids_agrees_with_get_vars_df():
    """The two must read the same metadata field by the same route. If they ever
    diverge, thinning starts deleting variants that stage 2 still considers
    causal candidates -- and nothing downstream would notice."""
    get_vars_df = _load_from_source('get_vars_df')
    ts = mixed_ts()

    from_thinning = set(nt.neutral_site_ids(ts).tolist())
    vars_df = get_vars_df(ts, Q_scaling=1.0, times_already_unscaled=True)
    from_stage2 = set(vars_df.loc[vars_df['selco'] == 0, 'id'].astype(int).tolist())

    assert from_thinning == from_stage2
    assert from_thinning, "the fixture must actually contain neutral sites"


def test_the_causal_predicate_never_sees_a_thinned_site():
    """causal_eligible selects selco != 0, and thinning deletes selco == 0, so
    the pool is identical before and after. This is the property that makes the
    drawn trait positions identical to the parent category's."""
    causal_eligible = _load_from_source('causal_eligible')
    get_vars_df = _load_from_source('get_vars_df')
    ts = mixed_ts()
    thinned, _ = nt.thin_neutral_sites(ts, 0.4, seed=11, measure_pi=False)

    def pool(t):
        df = get_vars_df(t, Q_scaling=1.0, times_already_unscaled=True)
        return set(causal_eligible(df, 0.0, 1e9, 0.0, False)['position'])

    assert pool(ts) == pool(thinned)


# ------------------------------------------------------------- sizing the cut

def test_pi_splits_into_its_two_classes():
    ts = mixed_ts()
    pi_neutral, pi_selected = nt.pi_components(ts)
    assert pi_neutral + pi_selected == pytest.approx(ts.diversity())
    assert pi_neutral > pi_selected     # the fixture's neutral class is common


@pytest.mark.parametrize("pi_neu,pi_sel,want", [
    (1.0, 0.0, 0.5),        # no selected pi: halve the neutral class
    (0.8, 0.2, 0.375),      # 0.5 * (1 - 0.25)
    (0.9, 0.1, 0.5 - 0.5 / 9),
])
def test_the_keep_fraction_solves_the_pi_target(pi_neu, pi_sel, want):
    assert nt.keep_fraction_for_pi_target(pi_neu, pi_sel) == pytest.approx(want)


def test_an_unreachable_target_raises_rather_than_clamping():
    """If the selected class alone already carries half of pi, no amount of
    neutral thinning gets there. Clamping to a small positive k would deliver a
    different experiment than the one asked for, silently."""
    with pytest.raises(ValueError, match="cannot be scaled"):
        nt.keep_fraction_for_pi_target(0.4, 0.6)


def test_the_solved_fraction_actually_lands_on_the_target():
    """End to end on the fake: solve for k, thin, and check pi really halved."""
    ts = mixed_ts(n_neutral=4000, n_selected=1000, seed=5)
    pi_neutral, pi_selected = nt.pi_components(ts)
    k = nt.keep_fraction_for_pi_target(pi_neutral, pi_selected, target=0.5)
    thinned, counts = nt.thin_neutral_sites(ts, k, seed=11)
    # The deletion is a random subset, so the realized ratio scatters around the
    # target; with 4000 neutral sites that scatter is small.
    assert counts["pi_ratio"] == pytest.approx(0.5, abs=0.02)
    assert thinned.diversity() == pytest.approx(0.5 * ts.diversity(), rel=0.05)


# ------------------------------------------------------------- the path segment

def test_the_segment_is_empty_when_nothing_is_thinned():
    """Every path ever written must stay reachable, so the historical value emits
    nothing -- the same rule LEGACY_CAUSAL_MIN_MAF and LEGACY_N_CENTRAL_TRAITS
    follow."""
    assert paths.neutral_thin_segment({}) == ""
    assert paths.neutral_thin_segment({"neutral_keep_fraction": 1.0}) == ""


def test_the_segment_names_the_fraction():
    assert paths.neutral_thin_segment({"neutral_keep_fraction": 0.375}) == ".nkeep_0375"
    assert paths.neutral_thin_segment({"neutral_keep_fraction": 0.42}) == ".nkeep_042"


def test_an_unthinned_run_tag_is_unchanged():
    cfg = {"pipeline": "human", "stage1_seed": 11, "causal_min_maf": 0.001}
    assert paths.stage2_run_tag(cfg) == "hts_11.cmaf_0.001"


def test_a_thinned_run_tag_cannot_collide_with_its_parent():
    """O1 and A1 share a stage-1 filename by design, so the segment is what keeps
    a mis-aimed stage2_search_dirs from adopting A's phenotypes."""
    parent = {"pipeline": "human", "stage1_seed": 11, "causal_min_maf": 0.001}
    child = dict(parent, neutral_keep_fraction=0.375)
    assert paths.stage2_run_tag(child) == "hts_11.cmaf_0.001.nkeep_0375"
    assert paths.stage2_run_tag(child) != paths.stage2_run_tag(parent)


def test_the_segment_composes_with_power_sampling():
    """Parameter set 2 is the psamp arm, so both segments appear at once."""
    cfg = {"pipeline": "human", "stage1_seed": 11, "causal_min_maf": 0,
           "causal_sampling": "power", "sampling_gwas_n": 100000,
           "neutral_keep_fraction": 0.375}
    assert paths.stage2_run_tag(cfg) == "hts_11.cmaf_0.psamp_100000.nkeep_0375"


def test_the_stage1_filename_carries_no_segment():
    """Stage 1 is SHARED with the parent category -- that is the pairing -- so the
    thinning must not reach the filename the adoption lookup matches on."""
    cfg = {"pipeline": "human", "stage1_seed": 11, "neutral_keep_fraction": 0.375}
    assert paths.stage1_full_filename(cfg) == "hts_11.ts"


# ----------------------------------------------------------- the provenance guard

def test_the_key_determines_stage2_and_is_recorded_as_such():
    """In STAGE2_KEYS, so pointing an O run at A's stage 2 is refused rather than
    silently reusing the un-thinned phenotypes. Not in ANALYSIS_KEYS: it is not a
    re-analysis of a fixed stage 2, it is a different stage 2."""
    assert "neutral_keep_fraction" in params_record.STAGE2_KEYS
    assert "neutral_keep_fraction" not in params_record.ANALYSIS_KEYS


def test_two_fractions_hash_to_different_stage2_uids():
    base = {k: None for k in params_record.STAGE2_KEYS}
    a = dict(base, neutral_keep_fraction=1.0)
    b = dict(base, neutral_keep_fraction=0.375)
    assert (params_record.uid(a, params_record.STAGE2_KEYS)
            != params_record.uid(b, params_record.STAGE2_KEYS))


# -------------------------------------------- incompatible with the other models

@pytest.mark.parametrize("key", ["synthetic_dfe_effects", "neutral_trait_vars"])
def test_the_snakefile_refuses_the_neutral_effect_models(key):
    """H/I/K/L draw their causal loci from selco == 0 and B/D draw their effect
    recipients from it. Thinning there would halve the causal pool while looking
    like a change in background density."""
    sf = (SIM_DIR / "Snakefile").read_text()
    guard = sf[sf.index('config.setdefault("neutral_keep_fraction"'):]
    assert key in guard
    assert "cannot be" in guard


@pytest.mark.parametrize("key", ["synthetic_dfe_effects", "neutral_trait_vars"])
def test_the_stage2_script_refuses_them_too(key):
    """A hand invocation does not go through the Snakefile."""
    guard = SRC[SRC.index("neutral_keep_fraction = args.neutral_keep_fraction"):]
    guard = guard[:guard.index("# Trait-associated")]
    assert key in guard
    assert "raise SystemExit" in guard


# ------------------------------------------------------------------ the configs

def test_the_configs_exist():
    assert HUMAN_LOWHET, "expected config/human_lowhet_2Mb_g5t20_r3.yaml"
    assert CATTLE_LOWHET, "expected config/cattle_lowhet_2Mb_g5t20_r3.yaml"


@pytest.mark.parametrize("path", HUMAN_LOWHET + CATTLE_LOWHET, ids=lambda p: p.name)
def test_the_lowhet_config_is_its_parent_with_one_key_changed(path):
    """The identity IS the experiment, as it is for J against A. Note what is NOT
    in `differ`: the SEEDS. O runs at A's seed and P at E's, which is what makes
    the variant sets nest and the contrast paired."""
    child = load_yaml(path)
    parent_name = (path.name.replace("human_lowhet_", "human_")
                   if path.name.startswith("human_lowhet_")
                   else path.name.replace("cattle_lowhet_",
                                          "cattle_baseline_from_midpoint_"))
    parent = load_yaml(CONFIG_DIR / parent_name)

    differ = {"basename", "neutral_keep_fraction", "workdir", "publishdir"}
    for k in sorted(set(parent) | set(child)):
        if k in differ:
            continue
        assert parent.get(k) == child.get(k), \
            f"{path.name} diverges from {parent_name} on {k!r}"


@pytest.mark.parametrize("path,want", [(p, 11) for p in HUMAN_LOWHET]
                                      + [(p, 51) for p in CATTLE_LOWHET],
                         ids=lambda x: getattr(x, "name", str(x)))
def test_the_seeds_are_the_parents_and_that_is_deliberate(path, want):
    """Every other category has a private band (B=2N, H=8N, K=11N, M=13N ...).
    O and P do not, because O1 must adopt A1's tree sequence and P1 must adopt
    E1's. Giving them a private band would make them run their own stage 1 and
    turn a paired contrast into an unpaired one."""
    cfg = load_yaml(path)
    assert cfg["stage1_seed"] == want
    assert cfg["stage2_seed"] == want


@pytest.mark.parametrize("path", HUMAN_LOWHET + CATTLE_LOWHET, ids=lambda p: p.name)
def test_the_config_does_not_build_its_own_stage1(path):
    assert load_yaml(path)["stage1_search_dirs"] == []


@pytest.mark.parametrize("path", HUMAN_LOWHET + CATTLE_LOWHET, ids=lambda p: p.name)
def test_the_workdir_leaf_is_not_the_parents(path):
    """A shared leaf would put O's stage 1 symlink in A's directory."""
    cfg = load_yaml(path)
    for key in ("workdir", "publishdir"):
        assert cfg[key].rstrip("/").endswith("lowhet"), f"{path.name}: {key}"


# ----------------------------------------------------------- downstream wiring

@pytest.mark.parametrize("letter,species", [("O", "human"), ("P", "cattle")])
def test_the_categories_are_known_to_the_summariser(letter, species):
    from helpers import summarize_coloc
    assert summarize_coloc.DEMOGRAPHY[letter] == species


@pytest.mark.parametrize("letter,species", [("O", "human"), ("P", "cattle")])
def test_the_archiver_knows_the_species(letter, species):
    """archive_round3_to_data2.sh globs [A-Z][0-9]* so it picks the replicates up
    on its own, but species_of() is a hand-written character class and a miss
    there writes '?' into MANIFEST.tsv rather than failing."""
    sh = (SIM_DIR / "archive_round3_to_data2.sh").read_text()
    m = re.search(r"species_of\(\).*", sh)
    assert m
    classes = dict(re.findall(r"\[([A-Z]+)\]\) echo (\w+)", m.group(0)))
    assert any(letter in cls and sp == species for cls, sp in classes.items()), \
        f"{letter} missing from species_of()"


@pytest.mark.parametrize("letter,species", [("O", "human"), ("P", "cattle")])
def test_the_notebook_knows_the_species(letter, species):
    """CATEGORY_SPECIES drives build_manifest()'s cross-check; an unknown letter
    makes every replicate of it warn as a species mismatch."""
    import json
    nb = json.loads((SIM_DIR.parent / "figures_and_tables"
                     / "figure2_revision2.ipynb").read_text())
    src = "".join("".join(c["source"]) for c in nb["cells"])
    block = src[src.index("CATEGORY_SPECIES = c("):]
    block = block[:block.index(")")]
    assert re.search(rf"\b{letter} = '{species}'", block)
    assert re.search(rf"^\s*{letter} = '{letter}\s", src, re.M), \
        f"{letter} missing from CATEGORY_LABEL"


# --------------------------------------------------------- the submit scripts

SUBMIT = ["submit_2Mb_r3_cmaf_fm001.sh", "submit_2Mb_r3_x20_psamp_fm001.sh"]


@pytest.mark.parametrize("name", SUBMIT)
def test_the_submit_script_maps_the_letters(name):
    sh = (SIM_DIR / name).read_text()
    assert 'O) echo "1${n}" ;; P) echo "5${n}" ;;' in sh, "O/P must run at A's and E's seeds"
    assert 'config/human_lowhet_2Mb_${CELL}_r3.yaml' in sh
    assert 'config/cattle_lowhet_2Mb_${CELL}_r3.yaml' in sh
    assert 'O) echo "A${1:1}" ;;' in sh
    assert 'P) echo "E${1:1}" ;;' in sh


@pytest.mark.parametrize("name", SUBMIT)
def test_they_adopt_stage1_rather_than_building_it(name, guard_letters):
    """The pairing depends on it, and building their own would silently give them
    a different genome at the same seed."""
    sh = (SIM_DIR / name).read_text()
    assert not ({"O", "P"} & guard_letters(sh, 'STAGE1_SRC["$ID"]=""'))
    assert "$(stage1_donor_of " in sh


def test_they_never_adopt_the_parents_stage2(guard_letters):
    """Only the cmaf script adopts stage 2 at all; the psamp one always rebuilds
    it, power sampling being a stage-2 change."""
    sh = (SIM_DIR / "submit_2Mb_r3_cmaf_fm001.sh").read_text()
    assert {"O", "P"} <= guard_letters(sh, 'STAGE2_SRC["$ID"]=""')


@pytest.mark.parametrize("name", SUBMIT)
def test_an_unset_keep_fraction_is_refused_before_launch(name):
    """Without the key O is a byte-identical copy of A -- a null result that looks
    like a measurement. Catch it at submit time, not in the summary table."""
    sh = (SIM_DIR / name).read_text()
    assert "neutral_keep_fraction:[[:space:]]*[0-9.]+" in sh


@pytest.mark.parametrize("name", ["submit_2Mb_r3_cmaf_replicates.sh",
                                  "submit_2Mb_r3_cmaf01_control.sh"])
def test_the_stage1_building_scripts_have_no_entry_for_them(name):
    """submit_2Mb_r3_cmaf_replicates.sh RUNS stage 1. An O entry there would
    simulate a fresh genome at seed 11 and overwrite the pairing. seed_of()'s
    `*) ERROR` is the guard; this test is what keeps it a guard."""
    sh = (SIM_DIR / name).read_text()
    assert 'O) echo "1${n}"' not in sh
    assert "human_lowhet" not in sh


# ------------------------------------------------------------ the calibration

CALIBRATE = SIM_DIR / "calibrate_lowhet.sh"


def test_the_calibration_script_exists_and_is_runnable():
    assert CALIBRATE.exists()
    assert os.access(CALIBRATE, os.X_OK)


def test_it_measures_the_donors_not_the_thinned_categories():
    """pi is measured on A's and E's tree sequences -- O and P do not have their
    own, and measuring an already-thinned one would solve for a second halving."""
    sh = CALIBRATE.read_text()
    assert "REPS:-A1 A2 A3 A4 A5 E1 E2 E3 E4 E5" in sh
    assert "calibrate only A and E" in sh


def test_the_cattle_call_carries_the_split_q_handoff():
    """Without handoff_ticks/deep_Q_scaling, add_neutral applies the recent
    (0.01-scaled) rate to the whole tree and loses ~99% of the neutral variation.
    The measurement would still succeed and report a selected share near 1 --
    a plausible number, and completely wrong."""
    sh = CALIBRATE.read_text()
    cattle = sh[sh.index("--species cattle"):]
    cattle = cattle[:cattle.index("else")]
    assert "--handoff_ticks 2400" in cattle
    assert "--deep_Q_scaling 1" in cattle
    assert "--Q_scaling 0.01" in cattle


def test_it_searches_the_same_stage1_roots_the_submit_scripts_do():
    """A donor resolved from a different root than the one the run will adopt
    would calibrate against a genome the arm never sees."""
    sh = CALIBRATE.read_text()
    submit = (SIM_DIR / "submit_2Mb_r3_x20_psamp_fm001.sh").read_text()
    for root in ("simulations_round_3_2Mb_g5t20_cmaf_0.001", "simulations_round_3_2Mb"):
        assert root in sh and root in submit


def test_an_infeasible_replicate_is_reported_and_not_averaged_in():
    """k is NA when the selected class alone already carries half of pi. Folding
    an NA into the mean as a zero would pin a fraction that halves nothing."""
    sh = CALIBRATE.read_text()
    assert 'keep_fraction_for_half_pi"] == "NA"' in sh
    assert "INFEASIBLE" in sh
    assert "do NOT pin a value" in sh


@pytest.mark.parametrize("name", SUBMIT)
def test_the_summary_lines_survive_set_e(name):
    """`VAR=$(for r in $REPS; do [[ ... ]] && printf ...; done)` inherits the LAST
    iteration's status, which is 1 whenever that replicate does not match. Under
    `set -euo pipefail` the assignment then kills the script before the submit
    loop -- silently, with nothing launched and no error printed.

    This bit for real: the pre-existing FRESH line aborted every launch whose
    REPS ended in a non-K/L/M/N letter, which is the common case. Each such
    assignment must end its subshell with `:`."""
    sh = (SIM_DIR / name).read_text()
    hazards = re.findall(r"^[A-Z_]+=\$\(for r in \$REPS;.*$", sh, re.M)
    assert hazards, f"{name}: expected at least one REPS summary line"
    for line in hazards:
        assert line.rstrip().endswith("; :)"), \
            f"{name}: `{line.strip()}` aborts under set -e when no replicate matches"
