"""API contract tests (offline: pubchem=false)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from fastapi.testclient import TestClient
import app as appmod

client = TestClient(appmod.app)
SMI = "CC(=O)Oc1ccccc1C(=O)O"


def test_resolve_ok_and_shape():
    r = client.get("/api/resolve", params={"q": "2-acetyloxybenzoic acid", "pubchem": "false"})
    assert r.status_code == 200
    j = r.json()
    for k in ("smiles_canonical", "smiles_kekule", "inchi", "inchikey", "formula", "names", "warnings", "input_kind"):
        assert k in j
    assert j["smiles_canonical"] == SMI and j["names"][0]["verified"]


def test_resolve_rejects_garbage_with_422():
    r = client.get("/api/resolve", params={"q": "definitely not a molecule", "pubchem": "false"})
    assert r.status_code == 422 and "OPSIN" in r.json()["detail"]


def test_structure_properties_spectra_modes_photo():
    s = client.get("/api/structure", params={"smiles": SMI}).json()
    assert "<svg" in s["depiction"]["svg"] and "M  END" in s["depiction"]["molblock"] and s["geometry"]["bond_lengths"]
    p = client.get("/api/properties", params={"smiles": SMI}).json()
    assert p["formula"] == "C9H8O4" and abs(p["molar_mass"]["value"] - 180.16) < 0.01
    for sol in ("chloroform", "CDCl3", "DMSO-d6", "water", "gas phase"):
        sp = client.get("/api/spectra", params={"smiles": SMI, "solvent": sol, "mhz": 500}).json()
        assert sp["h1"]["mhz"] == 500 and sp["solvent"]["name"]
    assert client.get("/api/spectra", params={"smiles": SMI, "solvent": "unobtainium"}).status_code == 400
    inc = client.get("/api/spectra", params={"smiles": SMI, "nmr": "increments"}).json()
    assert inc["nmr_level"] == "increments" and all(p["source"] == "increments" for p in inc["c13"]["peaks"])
    assert client.get("/api/spectra", params={"smiles": SMI, "nmr": "bogus"}).status_code == 400
    mm = client.get("/api/modes", params={"smiles": "CCO", "level": "mmff"}).json(); assert mm["level"] == "MMFF94s"
    m = client.get("/api/modes", params={"smiles": SMI}).json()
    assert len(m["modes"]) == 3 * 21 - 6 and all(x["freq"] > 0 for x in m["modes"])
    ph = client.get("/api/photo", params={"smiles": SMI, "solvent": "ethanol"}).json()
    assert ph["class"] in ("bright", "moderate", "weak", "non-emissive")


def test_verify_four_way():
    key = client.get("/api/resolve", params={"q": "(S)-ibuprofen", "pubchem": "false"}).json()["inchikey"]
    v = lambda x: client.get("/api/verify", params={"name": x, "inchikey": key}).json()["status"]
    assert v("(2S)-2-[4-(2-methylpropyl)phenyl]propanoic acid") == "identical"
    assert v("CC(C)Cc1ccc(C(C)C(O)=O)cc1") == "same_connectivity"
    assert v("toluene") == "different"
    assert v("qwertyuiop") == "uninterpretable"


def test_electronic():
    j = client.get("/api/electronic", params={"smiles": "c1ccccc1"}).json()
    assert j["gap_ev"] and 1.0 < j["gap_ev"] < 8.0 and abs(j["dipole_debye"]) < 0.05
    assert any(l["label"] == "HOMO" for l in j["levels"]) and len(j["charges"]) == 12
    j2 = client.get("/api/electronic", params={"smiles": "c1ccc2cc3ccccc3cc2c1"}).json()
    assert j2["gap_ev"] < j["gap_ev"]                     # anthracene gap < benzene gap


def test_giao_level():
    j = client.get("/api/spectra", params={"smiles": "c1ccccc1", "nmr": "giao"}).json()
    c = [p for p in j["c13"]["peaks"]]
    assert all(p["source"] == "giao" for p in c) and abs(c[0]["shift"] - 128.4) < 4
    h = [p for p in j["h1"]["peaks"] if not p["exchangeable"]]
    assert abs(h[0]["shift"] - 7.26) < 0.5 and h[0]["source"] == "giao"
    big = client.get("/api/spectra", params={"smiles": "CCCCCCCCCCCCCCCCCC", "nmr": "giao"}).json()
    assert big["giao_note"] and "limited" in big["giao_note"]


def test_orbital_cube():
    j = client.get("/api/orbital", params={"smiles": "C=C", "which": "homo"}).json()
    assert int(j["cube"].splitlines()[2].split()[0]) == 6 and j["label"] == "HOMO" and j["occ"] == 2
    j2 = client.get("/api/orbital", params={"smiles": "C=C", "which": "lumo", "offset": 1}).json()
    assert j2["label"] == "LUMO+1" and j2["energy_ev"] > j["energy_ev"]
    # symmetric molecule: HOMO centroid sits on the nuclear charge centre
    import math
    d = math.dist(j["centroid"], j["nuc_center"]); assert d < 0.05
    assert client.get("/api/orbital", params={"smiles": "C=C", "which": "sigma"}).status_code == 400


def test_read_endpoint_roundtrip():
    txt = client.get("/api/cdxml", params={"smiles": SMI}).text
    j = client.post("/api/read", params={"fmt": "cdxml"}, content=txt).json()
    from rdkit import Chem
    assert Chem.MolToInchiKey(Chem.MolFromSmiles(j["smiles"])) == Chem.MolToInchiKey(Chem.MolFromSmiles(SMI))
    assert client.post("/api/read", params={"fmt": "cdxml"}, content="not xml").status_code == 422


def test_cdxml_export():
    import xml.dom.minidom as MD
    r = client.get("/api/cdxml", params={"smiles": SMI})
    assert r.status_code == 200 and "CDXML" in r.text
    MD.parseString(r.text)
    assert r.text.count("<n ") == 13 and r.text.count("<b ") == 13
    rn = client.get("/api/cdxml", params={"smiles": SMI, "name": "2-acetyloxybenzoic acid", "numbers": "true"})
    assert rn.text.count('size="7"') >= 8            # locant annotations present
    MD.parseString(rn.text)


def test_cdxml_graph_roundtrip():
    """Parse the CDXML back and compare the molecular graph to the source structure."""
    import xml.etree.ElementTree as ET
    from rdkit import Chem
    for smi in [SMI, "C[C@H](N)C(=O)O", "O=[N+]([O-])c1ccc(O)cc1", "C/C=C/C(=O)O"]:
        mol = Chem.MolFromSmiles(smi)
        txt = client.get("/api/cdxml", params={"smiles": smi}).text
        root = ET.fromstring(txt.split("?>", 1)[1].split(">", 1)[1].rsplit("</CDXML", 1)[0].join(["<CDXML", "</CDXML>"])) if False else ET.fromstring(txt[txt.index("<CDXML"):])
        frag = root.find(".//fragment")
        atoms = frag.findall("n"); bonds = frag.findall("b")
        assert len(atoms) == mol.GetNumAtoms() and len(bonds) == mol.GetNumBonds(), smi
        z = {a.get("id"): int(a.get("Element", "6")) for a in atoms}
        rw = Chem.RWMol()
        ids = {}
        for a in atoms:
            at = Chem.Atom(z[a.get("id")])
            if a.get("Charge"):
                at.SetFormalCharge(int(a.get("Charge")))
            ids[a.get("id")] = rw.AddAtom(at)
        for b in bonds:
            rw.AddBond(ids[b.get("B")], ids[b.get("E")], {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[int(b.get("Order", "1"))])
        back = rw.GetMol()
        Chem.SanitizeMol(back)
        # graph identity: same connectivity-layer InChIKey (stereo lives in wedges/coords, not the graph)
        ik = lambda m: Chem.MolToInchiKey(m).split("-")[0]
        assert ik(back) == ik(mol), smi
        # stereo present as wedges when the molecule has stereocentres
        if Chem.FindMolChiralCenters(mol):
            assert "Wedge" in txt or "WedgedHash" in txt, smi


def test_methods_and_solvents_and_about():
    m = client.get("/api/methods").json()
    assert {"nmr", "uv", "photoluminescence", "name_generation"} <= set(m) and all("fails_at" in v for v in m.values())
    s = client.get("/api/solvents").json()
    assert any(x["name"] == "chloroform" and x["nmr"] == "CDCl3" for x in s)
    assert "OPSIN" in client.get("/api/about").text
    assert client.get("/").headers["cache-control"] == "no-cache"
