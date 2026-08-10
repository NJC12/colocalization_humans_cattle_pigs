"""Tests for the GWAS-locus/GTEx-partner join the scorer uses.

WHY THIS MATTERS MORE THAN IT LOOKS

`summarize_coloc.py` pairs a GWAS trait with a GTEx trait by name equality, which
cannot distinguish "this locus has no eQTL to find" from "colocalization failed".
Those are opposite results. Where the GTEx causal set is topped up rather than
intersected -- under `causal_sampling: power`, and under the drawn-DFE arms H and
I in both schemes -- a substantial share of GWAS loci have no partner at all, and
every RCP those loci earn against some OTHER GTEx trait is a false positive.

So `shared` is the flag that splits the trait dump into a power denominator and a
false-positive denominator. Reading it wrong is silent: a parse that returned
False for everything would report a false-positive rate over loci that DID have a
partner, which is a different and much smaller number.
"""

import os

import pytest

from helpers.summarize_coloc import shared_partner_by_trait


HEADER = "gwas_trait\tgtex_trait\tposition\tshared\n"


def write_sidecar(tmp_path, cat="cgwas", suffix="gwas_5_gtex_20_maf_0.tsv", rows=()):
    """The real shape: pandas.to_csv of trait_partner_table's DataFrame, so the
    booleans arrive as the strings 'True'/'False'."""
    p = tmp_path / f"{cat}_trait_partners_{suffix}"
    with open(p, "w") as fh:
        fh.write(HEADER)
        for trait, partner, pos, shared in rows:
            fh.write(f"{trait}\t{partner}\t{pos}\t{shared}\n")
    return str(tmp_path)


def test_it_reads_the_booleans_pandas_actually_writes():
    """`to_csv` writes Python bools as 'True'/'False', not '1'/'0'. A truthiness
    test on the raw string would call 'False' shared."""
    import tempfile, pathlib
    with tempfile.TemporaryDirectory() as d:
        d = pathlib.Path(d)
        write_sidecar(d, rows=[("tr100", "tr100", 100, "True"),
                               ("tr200", "", 200, "False")])
        got = shared_partner_by_trait(str(d), "cgwas")
    assert got == {"tr100": True, "tr200": False}


def test_the_unpartnered_loci_are_the_false_positive_denominator():
    """The H1-at-floor-0 shape: 10 of 50 partnered. The other 40 cannot produce a
    true colocalization, so they are exactly the loci a false-positive rate is
    measured over."""
    import tempfile, pathlib
    rows = [(f"tr{i}", f"tr{i}" if i < 10 else "", i, str(i < 10)) for i in range(50)]
    with tempfile.TemporaryDirectory() as d:
        d = pathlib.Path(d)
        write_sidecar(d, rows=rows)
        got = shared_partner_by_trait(str(d), "cgwas")
    assert sum(got.values()) == 10
    assert sum(1 for v in got.values() if not v) == 40


def test_a_missing_sidecar_gives_no_answer_rather_than_a_wrong_one():
    """Every run predating the sidecar. `.get(tr)` then yields None, and the dump
    column says 'unknown' instead of claiming every locus is unpartnered."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        assert shared_partner_by_trait(d, "cgwas") == {}
        assert shared_partner_by_trait(d, "cgwas").get("tr100") is None


def test_a_sidecar_without_the_column_is_refused_not_guessed():
    import tempfile, pathlib
    with tempfile.TemporaryDirectory() as d:
        d = pathlib.Path(d)
        p = d / "cgwas_trait_partners_x.tsv"
        p.write_text("gwas_trait\tgtex_trait\tposition\ntr100\ttr100\t100\n")
        assert shared_partner_by_trait(str(d), "cgwas") == {}


def test_it_picks_the_panel_it_was_asked_for():
    """hgwas and cgwas sidecars can sit in one stage-2 dir only by accident, but
    a glob that matched both would silently take whichever sorted first."""
    import tempfile, pathlib
    with tempfile.TemporaryDirectory() as d:
        d = pathlib.Path(d)
        write_sidecar(d, cat="cgwas", rows=[("tr1", "tr1", 1, "True")])
        write_sidecar(d, cat="hgwas", rows=[("tr2", "", 2, "False")])
        assert shared_partner_by_trait(str(d), "cgwas") == {"tr1": True}
        assert shared_partner_by_trait(str(d), "hgwas") == {"tr2": False}


@pytest.mark.parametrize("value,expected", [
    ("True", True), ("true", True), ("TRUE", True), ("1", True),
    ("False", False), ("false", False), ("0", False), ("", False),
])
def test_boolean_spellings(value, expected):
    import tempfile, pathlib
    with tempfile.TemporaryDirectory() as d:
        d = pathlib.Path(d)
        write_sidecar(d, rows=[("tr1", "", 1, value)])
        assert shared_partner_by_trait(str(d), "cgwas")["tr1"] is expected
