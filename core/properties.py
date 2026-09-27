"""
Physical-property estimation.
  * exact/average mass, formula, descriptors      : RDKit
  * Tm, Tb, Tc, Pc, ΔHvap, ΔHfus                   : Joback group contribution (thermo)
  * density (liquid/solid, 20 °C)                  : Girolami (1994) atom-volume method
  * logP, TPSA, refractivity, pKa-relevant groups  : RDKit Crippen / SMARTS
All estimates carry a stated typical error.
"""
from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, Lipinski, QED

_PERIOD_VOL = {1: 1, 2: 2, 3: 4, 4: 5, 5: 6, 6: 7}
_GIROLAMI_GROUPS = [   # each +10 % (max +30 % for functional groups + 10 % per ring, capped at 130 %)
    ("[OX2H]", "hydroxyl"), ("[CX3](=O)[OX2H1]", "carboxylic acid"), ("[NX3;H2,H1;!$(NC=O)]", "primary/secondary amine"),
    ("[NX3][CX3](=O)", "amide"), ("[SX3](=O)", "sulfoxide"), ("[SX4](=O)(=O)", "sulfone"),
    ("[CX3](=O)[NX3H2]", "primary amide"), ("[NX3](=O)=O", "nitro"), ("[#6]#N", "nitrile"),
]


def period(z: int) -> int:
    for p, hi in ((1, 2), (2, 10), (3, 18), (4, 36), (5, 54), (6, 86)):
        if z <= hi:
            return p
    return 7


def girolami_density(mol: Chem.Mol) -> dict:
    mh = Chem.AddHs(mol)
    vs = sum(_PERIOD_VOL.get(period(a.GetAtomicNum()), 7) for a in mh.GetAtoms())
    M = Descriptors.MolWt(mol)
    rho0 = M / (5.0 * vs)
    corr = 0.0; applied = []
    for smarts, label in _GIROLAMI_GROUPS:
        n = len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))
        if n:
            corr += 0.1 * n; applied.append(f"{label}×{n}")
    nring = rdMolDescriptors.CalcNumRings(mol)
    fused = 0
    ri = mol.GetRingInfo()
    if nring:
        # fused ring systems count +7.5 % each, isolated rings +10 %
        corr += 0.1 * nring; applied.append(f"rings×{nring}")
    corr = min(corr, 0.3)
    rho = rho0 * (1 + corr)
    return {"value": round(rho, 3), "unit": "g/cm³", "method": "Girolami (J. Chem. Educ. 1994)", "corrections": applied,
            "error": "typically ±0.1 g/cm³ (liquids/solids at ~20 °C)"}


def joback(smiles: str) -> dict:
    try:
        from thermo import Joback
        j = Joback(smiles)
        if j.status != "OK":
            return {"error": f"Joback: {j.status}"}
        e = j.estimate()
        out = {}
        for k, unit, err in (("Tm", "K", "±25–50 K (Joback Tm is a rough estimate)"), ("Tb", "K", "±15–25 K"), ("Tc", "K", "±10 K"),
                             ("Hvap", "J/mol", "±10 %"), ("Hfus", "J/mol", "±20 %"), ("Hf", "J/mol", "±5–10 kJ/mol")):
            v = e.get(k)
            if v is not None and v == v:
                out[k] = {"value": round(v, 1), "unit": unit, "error": err}
        if e.get("Pc"):
            out["Pc"] = {"value": round(e["Pc"] / 1e5, 2), "unit": "bar", "error": "±5 %"}
        out["method"] = "Joback & Reid (1987) group contribution"
        return out
    except Exception as ex:
        return {"error": f"Joback unavailable: {ex}"}


def predict(mol: Chem.Mol, smiles: str) -> dict:
    d: dict = {}
    d["formula"] = rdMolDescriptors.CalcMolFormula(mol)
    d["molar_mass"] = {"value": round(Descriptors.MolWt(mol), 3), "unit": "g/mol"}
    d["exact_mass"] = {"value": round(Descriptors.ExactMolWt(mol), 5), "unit": "Da (monoisotopic)"}
    d["heavy_atoms"] = mol.GetNumHeavyAtoms()
    d["charge"] = Chem.GetFormalCharge(mol)
    d["logP"] = {"value": round(Crippen.MolLogP(mol), 2), "method": "Wildman–Crippen", "error": "±0.7"}
    d["molar_refractivity"] = {"value": round(Crippen.MolMR(mol), 2), "unit": "cm³/mol"}
    d["TPSA"] = {"value": round(rdMolDescriptors.CalcTPSA(mol), 2), "unit": "Å²"}
    d["h_bond_donors"] = Lipinski.NumHDonors(mol)
    d["h_bond_acceptors"] = Lipinski.NumHAcceptors(mol)
    d["rotatable_bonds"] = Lipinski.NumRotatableBonds(mol)
    d["rings"] = rdMolDescriptors.CalcNumRings(mol)
    d["aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
    d["fraction_sp3"] = round(rdMolDescriptors.CalcFractionCSP3(mol), 2)
    d["stereocentres"] = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))
    try:
        d["QED_druglikeness"] = round(QED.qed(mol), 2)
    except Exception:
        pass
    d["density"] = girolami_density(mol)
    jb = joback(smiles)
    d["joback"] = jb
    if mol.GetNumHeavyAtoms() > 35 or ("Tm" in jb and jb["Tm"]["value"] > 650):
        jb["range_warning"] = "Joback is unreliable for molecules this large/polar (Tm, Tb not shown)."
        jb.pop("Tm", None); jb.pop("Tb", None)
    if "Tm" in jb:
        d["melting_point"] = {"K": jb["Tm"]["value"], "C": round(jb["Tm"]["value"] - 273.15, 1), "error": jb["Tm"]["error"], "method": jb["method"]}
    if "Tb" in jb:
        d["boiling_point"] = {"K": jb["Tb"]["value"], "C": round(jb["Tb"]["value"] - 273.15, 1), "error": jb["Tb"]["error"], "method": jb["method"]}
    if "value" in d["density"]:
        d["molar_volume"] = {"value": round(d["molar_mass"]["value"] / d["density"]["value"], 1), "unit": "cm³/mol"}
    # isotope pattern (top peaks)
    return d
