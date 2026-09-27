"""Data-driven NMR model and GFN2-xTB vibrations."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import pytest
from rdkit import Chem
from core import hose as H, spectra as X, vibrations as V

model = H.load()


@pytest.mark.skipif(model is None, reason="HOSE model not built")
def test_hose_model_validation_numbers_are_recorded_and_sane():
    v = model.validation
    assert v["13C"]["MAE"] < 3.5 and v["13C"]["coverage"] > 0.95 and v["13C"]["by_depth"][4]["MAE"] < 1.6
    assert v["1H"]["MAE"] < 0.5


@pytest.mark.skipif(model is None, reason="HOSE model not built")
@pytest.mark.parametrize("smiles,exp", [
    ("CC(=O)c1ccccc1", {198.1, 137.1, 133.1, 128.6, 128.3, 26.6}),   # acetophenone
    ("CCOC(C)=O", {171.1, 60.4, 21.0, 14.2}),                       # ethyl acetate
    ("c1ccc2ccccc2c1", {133.5, 127.9, 125.8}),                       # naphthalene
    ("OC(=O)c1ccccc1", {133.8, 130.3, 129.4, 128.5}),                # benzoic acid (COOH carbon is 172.6 CDCl3 / 167.4 DMSO — solvent-dependent, excluded)
])
def test_c13_matches_experiment_within_3ppm(smiles, exp):
    pred = [e["shift"] for e in X.c_environments(Chem.MolFromSmiles(smiles))]
    for e in exp:
        assert min(abs(p - e) for p in pred) < 3.0, (smiles, e, pred)
    assert any(e["source"] == "nmrshiftdb2" for e in X.c_environments(Chem.MolFromSmiles(smiles)))


def test_hose_code_is_canonical_and_depth_sensitive():
    m = Chem.MolFromSmiles("CCOC(C)=O")
    a = H.hose(m, 0, 2); b = H.hose(Chem.MolFromSmiles("O=C(C)OCC"), 5, 2)   # same methyl written differently
    assert a == b and H.hose(m, 0, 1) != H.hose(m, 0, 3)


def test_gfn2_carbonyl_and_ring_modes():
    g = V.normal_modes_gfn2(Chem.MolFromSmiles("CC(=O)C"))
    assert g.get("level") == "GFN2-xTB" and len(g["modes"]) >= 3 * 10 - 7      # a methyl torsion may come out imaginary and is dropped
    co = max((m for m in g["modes"] if 1650 <= m["freq"] <= 1800), key=lambda m: m["intensity"])
    assert abs(co["freq"] - 1731) < 25 and co["intensity"] > 0.9          # the C=O is the strongest IR band
    b = V.normal_modes_gfn2(Chem.MolFromSmiles("c1ccccc1"))
    assert any(660 <= m["freq"] <= 690 and m["intensity"] > 0.15 for m in b["modes"])


def test_best_level_reports_itself():
    r = V.normal_modes_best(Chem.MolFromSmiles("CCO"))
    assert r["level"] in ("GFN2-xTB", "MMFF94s") and r["modes"]
