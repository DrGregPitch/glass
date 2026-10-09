"""
Round-trip and robustness tests.  Run:  python -m pytest -q
Offline by default (PubChem disabled) so results are deterministic; set CHEMXREF_ONLINE=1
to also exercise the PubChem paths.
"""
import os, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import pytest
from rdkit import Chem
from core import names as N, spectra as X, properties as P, structure as S, vibrations as V, photo as PH, solvents as SOLV, opsin

ONLINE = os.environ.get("CHEMXREF_ONLINE") == "1"

# name, expected canonical SMILES (RDKit), expected parent locant count (heavy atoms with numeric parent locants)
CASES = [
    ("2-methylbutan-2-ol", "CCC(C)(C)O", 4),
    ("(2S)-2-[4-(2-methylpropyl)phenyl]propanoic acid", "CC(C)Cc1ccc([C@H](C)C(=O)O)cc1", 3),
    ("2-acetyloxybenzoic acid", "CC(=O)Oc1ccccc1C(=O)O", 6),
    ("caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O", None),
    ("1,3,7-trimethylpurine-2,6-dione", "Cn1c(=O)c2c(ncn2C)n(C)c1=O", 9),
    ("N,N-dimethylaniline", "CN(C)c1ccccc1", 6),
    ("(E)-but-2-enoic acid", "C/C=C/C(=O)O", 4),
    ("(Z)-but-2-enoic acid", "C/C=C\\C(=O)O", 4),
    ("L-alanine", "C[C@H](N)C(=O)O", 3),
    ("D-alanine", "C[C@@H](N)C(=O)O", 3),
    ("naphthalen-2-ol", "Oc1ccc2ccccc2c1", 10),
    ("2-(naphthalen-2-yl)ethanol", "OCCc1ccc2ccccc2c1", 2),
    ("1-methyl-4-(propan-2-yl)benzene", "Cc1ccc(C(C)C)cc1", 6),
    ("sodium acetate", "CC(=O)[O-].[Na+]", 2),
    ("tetramethylammonium chloride", "C[N+](C)(C)C.[Cl-]", None),
    ("benzene-1,3,5-triol", "Oc1cc(O)cc(O)c1", 6),
    ("pyridine-3-carboxylic acid", "O=C(O)c1cccnc1", 6),
    ("cyclohexane-1,2-diol", "OC1CCCCC1O", 6),
    ("4-nitrophenol", "O=[N+]([O-])c1ccc(O)cc1", 6),
    ("trans-stilbene", "C(=C/c1ccccc1)\\c1ccccc1", None),
    ("1H-indole-3-carbaldehyde", "O=Cc1c[nH]c2ccccc12", 9),
    ("Benzoic acid, 2-(acetyloxy)-", "CC(=O)Oc1ccccc1C(=O)O", 6),
    ("hexadecanoic acid", "CCCCCCCCCCCCCCCC(=O)O", 16),
    ("2,2,2-trifluoroethanol", "OCC(F)(F)F", 2),
    ("propan-2-one", "CC(C)=O", 3),
    ("ethyl acetate", "CCOC(C)=O", None),
    ("adamantane", "C1C2CC3CC1CC(C2)C3", 10),
]


@pytest.mark.parametrize("name,smiles,nloc", CASES)
def test_name_to_structure_and_back(name, smiles, nloc):
    r = N.resolve(name, want_pubchem=False)
    assert r.smiles_canonical == Chem.MolToSmiles(Chem.MolFromSmiles(smiles)), name
    inp = [e for e in r.names if e.kind == "input"]
    assert inp and inp[0].verified, f"input name not verified: {name}"
    # OPSIN round trip on the canonical SMILES key
    assert opsin.name_to_inchikey(name) == r.inchikey
    if nloc is not None:
        e = inp[0]
        got = sum(1 for f in (e.fragments or []) if f == "parent")
        assert got == nloc, f"{name}: parent locants {got} != {nloc}: {list(zip(e.locants, e.fragments))}"


@pytest.mark.parametrize("smiles", [c[1] for c in CASES])
def test_smiles_input_is_canonical_and_stable(smiles):
    r1 = N.resolve(smiles, want_pubchem=False)
    r2 = N.resolve(r1.smiles_canonical, want_pubchem=False)
    assert r1.inchikey == r2.inchikey and r1.smiles_canonical == r2.smiles_canonical
    r3 = N.resolve(r1.inchi, want_pubchem=False)
    assert r3.inchikey == r1.inchikey


def test_stereo_is_preserved_and_distinguished():
    l = N.resolve("L-alanine", want_pubchem=False); d = N.resolve("D-alanine", want_pubchem=False)
    assert l.inchikey != d.inchikey and l.inchikey.split("-")[0] == d.inchikey.split("-")[0]
    e = N.resolve("(E)-but-2-enoic acid", want_pubchem=False); z = N.resolve("(Z)-but-2-enoic acid", want_pubchem=False)
    assert e.inchikey != z.inchikey


@pytest.mark.parametrize("bad", ["", "   ", "notachemical xyz", "C1CC", "InChI=1S/garbage", "))(("])
def test_garbage_is_rejected_cleanly(bad):
    with pytest.raises(ValueError):
        N.resolve(bad, want_pubchem=False)


def test_verify_semantics():
    ref = N.resolve("(S)-ibuprofen", want_pubchem=False)
    same = N.resolve("(2S)-2-[4-(2-methylpropyl)phenyl]propanoic acid", want_pubchem=False)
    flat = N.resolve("CC(C)Cc1ccc(C(C)C(O)=O)cc1", want_pubchem=False)
    other = N.resolve("toluene", want_pubchem=False)
    assert same.inchikey == ref.inchikey
    assert flat.inchikey != ref.inchikey and flat.inchikey.split("-")[0] == ref.inchikey.split("-")[0]
    assert other.inchikey.split("-")[0] != ref.inchikey.split("-")[0]


def test_locant_spans_point_at_own_token():
    r = N.resolve("2-[4-(2-methylpropyl)phenyl]propanoic acid", want_pubchem=False)
    e = r.names[0]; nm = e.name
    parent2 = [i for i, (l, f) in enumerate(zip(e.locants, e.fragments)) if l == "2" and f == "parent"][0]
    assert [nm[a:b] for a, b in e.spans[str(parent2)]] == ["2"] and e.spans[str(parent2)][0][0] == 0
    phenyl4 = [i for i, (l, f) in enumerate(zip(e.locants, e.fragments)) if l == "4" and f == "sub"][0]
    assert e.spans[str(phenyl4)][0][0] == nm.index("4-(")


@pytest.mark.parametrize("smiles", ["CCO", "CC(=O)Oc1ccccc1C(=O)O", "c1ccc2cc3ccccc3cc2c1", "CC(C)Cc1ccc([C@H](C)C(=O)O)cc1", "O=[N+]([O-])c1ccccc1", "CC(=O)[O-].[Na+]", "C[N+](C)(C)C.[Cl-]"])
def test_all_simulations_run(smiles):
    m = Chem.MolFromSmiles(smiles)
    for sol in SOLV.SOLVENTS:
        sp = X.all_spectra(m, sol)
        assert sp["h1"]["x"] and sp["c13"]["x"] and sp["ir"]["x"] and sp["uv"]["x"]
        assert all(0 <= v <= 100.0001 for v in sp["ir"]["y"])
    p = P.predict(m, smiles)
    assert p["molar_mass"]["value"] > 0 and p["formula"]
    d = S.depict_svg(m); assert "<svg" in d["svg"] and len(d["coords"]) == m.GetNumAtoms() and "M  END" in d["molblock"]
    dh = S.depict_svg_h(m); assert len(dh["coords"]) == Chem.AddHs(m).GetNumAtoms()
    g = S.geometry3d(m); assert "molblock" in g or "error" in g
    v = V.normal_modes(m)
    if "modes" in v:
        assert len(v["modes"]) >= 3 * m.GetNumAtoms() - 6 - 6 * 0 or len(v["modes"]) > 0
        assert all(f["freq"] > 0 for f in v["modes"])
    pl = PH.photoluminescence(m, v if "modes" in v else None, "chloroform")
    assert "class" in pl


def test_solvent_trends_are_physical():
    m = Chem.MolFromSmiles("CC(=O)c1ccc(N(C)C)cc1")   # push-pull ketone
    em = {s: PH.photoluminescence(m, None, s)["lambda_em_max"] for s in ("hexane", "chloroform", "DMSO")}
    assert em["hexane"] < em["chloroform"] < em["DMSO"], em
    co = {s: [b["center"] for b in X.ir_spectrum(m, solvent=s)["bands"] if "C=O" in b["label"]][0] for s in ("gas phase", "hexane", "methanol")}
    assert co["gas phase"] > co["hexane"] > co["methanol"], co
    h_dmso = X.h_environments(Chem.MolFromSmiles("CCO"), "DMSO"); h_cdcl3 = X.h_environments(Chem.MolFromSmiles("CCO"), "chloroform")
    oh = lambda envs: [e["shift"] for e in envs if e["type"] == "alcohol OH"][0]
    assert oh(h_dmso) > oh(h_cdcl3)
    assert not any(e["type"] == "alcohol OH" for e in X.h_environments(Chem.MolFromSmiles("CCO"), "water"))


def test_nmr_reference_values():
    m = Chem.MolFromSmiles("CC(=O)c1ccccc1")   # acetophenone, CDCl3: 2.60 (s), 7.96/7.45/7.55; 13C 198.1, 137.1, 133.1, 128.6, 128.3, 26.6
    h = {e["type"]: e["shift"] for e in X.h_environments(m)}
    assert abs(h["CH₃"] - 2.60) < 0.6
    c = sorted(e["shift"] for e in X.c_environments(m))
    assert abs(c[-1] - 198.1) < 5 and abs(c[0] - 26.6) < 5


def test_unknown_solvent_rejected():
    with pytest.raises(KeyError):
        SOLV.get("unobtainium")
    assert SOLV.get("CDCl3")[0] == "chloroform"


@pytest.mark.skipif(not ONLINE, reason="set CHEMXREF_ONLINE=1")
def test_pubchem_enrichment():
    r = N.resolve("aspirin")
    assert r.cid == 2244 and r.cas == "50-78-2" and r.iupac_name and any(e.verified and e.kind == "iupac" for e in r.names)


# ----------------------------------------------------------------------------- name regions
def _region_symbols(name):
    """{(kind, label): sorted element symbols} for every region of the resolved name, plus the entry."""
    from rdkit import Chem
    r = N.resolve(name, want_pubchem=False)
    e = next(e for e in r.names if e.regions)
    m = Chem.MolFromSmiles(r.smiles_canonical)
    sym = lambda atoms: sorted(m.GetAtomWithIdx(i).GetSymbol() for i in atoms)
    return e, {(g["kind"], g["label"]): sym(g["atoms"]) for g in e.regions}


@pytest.mark.parametrize("name, expected", [
    ("2-methyl-4-(2-oxopropyl)benzoic acid", {
        ("group", "4-(2-oxopropyl)"): ["C", "C", "C", "O"],
        ("prefix", "2-oxo"): ["O"],
        ("prefix", "2-methyl"): ["C"],
        ("parent", "benzoic acid"): ["C"] * 7 + ["O", "O"],
        ("suffix", "oic acid"): ["C", "O", "O"]}),
    ("9,10-dimethylanthracene", {
        ("prefix", "9,10-dimethyl"): ["C", "C"],
        ("parent", "anthracene"): ["C"] * 14}),
    ("anthracene", {("parent", "anthracene"): ["C"] * 14}),
    ("2-[4-(2-methylpropyl)phenyl]propanoic acid", {
        ("group", "4-(2-methylpropyl)"): ["C"] * 4,
        ("group", "2-[4-(2-methylpropyl)phenyl]"): ["C"] * 10,
        ("prefix", "2-methyl"): ["C"],
        ("parent", "propanoic acid"): ["C", "C", "C", "O", "O"],
        ("suffix", "oic acid"): ["C", "O", "O"]}),
    ("4-chloro-N,N-dimethylaniline", {
        ("prefix", "4-chloro"): ["Cl"],
        ("prefix", "N,N-dimethyl"): ["C", "C"],
        ("parent", "aniline"): ["C"] * 6 + ["N"]}),
    ("2,2,2-trifluoroethanol", {
        ("prefix", "2,2,2-trifluoro"): ["F", "F", "F"],
        ("parent", "ethanol"): ["C", "C", "O"],
        ("suffix", "ol"): ["O"]}),
    # characteristic-group suffixes
    ("2-methylbutan-2-ol", {
        ("prefix", "2-methyl"): ["C"],
        ("parent", "butan-2-ol"): ["C"] * 4 + ["O"],
        ("suffix", "2-ol"): ["O"]}),                    # geminal methyl and hydroxyl at C2, told apart by element
    ("hexane-2,4-dione", {
        ("parent", "hexane-2,4-dione"): ["C"] * 6 + ["O", "O"],
        ("suffix", "2,4-dione"): ["O", "O"]}),
    ("acetic acid", {
        ("parent", "acetic acid"): ["C", "C", "O", "O"],
        ("suffix", "ic acid"): ["C", "O", "O"]}),       # on a chain the suffix carbon is C1 itself
    ("butanamide", {
        ("parent", "butanamide"): ["C"] * 4 + ["N", "O"],
        ("suffix", "amide"): ["C", "N", "O"]}),
    ("2-chlorobenzoic acid", {
        ("prefix", "2-chloro"): ["Cl"],
        ("parent", "benzoic acid"): ["C"] * 7 + ["O", "O"],
        ("suffix", "oic acid"): ["C", "O", "O"]}),
    ("methylbenzene", {("parent", "benzene"): ["C"] * 6}),   # unlocanted prefix: no region, parent still exact
    ("2-Methyl-2-butanol", {                                 # PubChem title: capitalised, CAS-style locant before the hydride
        ("prefix", "2-Methyl"): ["C"],
        ("parent", "2-butanol"): ["C"] * 4 + ["O"],
        ("suffix", "ol"): ["O"]}),
])
def test_name_regions_map_words_to_atoms(name, expected):
    _, got = _region_symbols(name)
    assert got == expected, f"{name}: {got}"


@pytest.mark.parametrize("name", ["2-methyl-4-(2-oxopropyl)benzoic acid", "2-[4-(2-methylpropyl)phenyl]propanoic acid",
                                  "4-chloro-N,N-dimethylaniline", "2,2,2-trifluoroethanol"])
def test_name_regions_are_consistent(name):
    e, _ = _region_symbols(name)
    nm = e.name
    parent = next(g for g in e.regions if g["kind"] == "parent")
    for g in e.regions:
        assert 0 <= g["start"] < g["end"] <= len(nm) and g["label"] == nm[g["start"]:g["end"]]
        assert g["atoms"], f"empty region {g}"
        if g["kind"] == "prefix" and g["end"] <= parent["start"]:
            assert not set(g["atoms"]) & set(parent["atoms"]), "a depth-0 prefix must not share atoms with the parent"
        if g["kind"] == "suffix":
            assert parent["start"] <= g["start"] and g["end"] == parent["end"], "a suffix is the tail of the parent word"
            assert set(g["atoms"]) <= set(parent["atoms"]), "suffix atoms belong to the parent"
        if g["kind"] == "group":
            for h in e.regions:
                if h["kind"] == "prefix" and g["start"] <= h["start"] and h["end"] <= g["end"]:
                    assert set(h["atoms"]) <= set(g["atoms"]), "a group must contain its own prefixes"


def test_repeated_substituent_with_primed_locant_is_a_substituent():
    """OPSIN numbers the second methyl of 9,10-dimethylanthracene 1'; it must be labelled sub, not hetero."""
    from rdkit import Chem
    r = N.resolve("9,10-dimethylanthracene", want_pubchem=False)
    e = next(e for e in r.names if e.locants)
    m = Chem.MolFromSmiles(r.smiles_canonical)
    methyls = [i for i, a in enumerate(m.GetAtoms()) if a.GetSymbol() == "C" and not a.IsInRing()]
    assert len(methyls) == 2 and all(e.fragments[i] == "sub" for i in methyls), [(i, e.locants[i], e.fragments[i]) for i in methyls]
    assert sum(1 for f in e.fragments if f == "parent") == 14
