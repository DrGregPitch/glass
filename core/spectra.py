"""
Spectrum simulation from structure: fully offline, rule based.

 1H NMR  : Pretsch/Shoolery additivity (sp3), Pascual (alkene), Hammett-type ring
           increments (aromatic); first-order J splitting; solvent-dependent exchangeable
           protons; spectrometer frequency adjustable (J in Hz -> ppm scale).
 13C NMR : Grant–Paul (alkanes) + heteroatom / functional group increments, aromatic
           ipso/ortho/meta/para increments, carbonyl classes.  Proton decoupled.
 IR      : group-frequency table (SMARTS -> bands), rendered as transmittance.
 UV-Vis  : Woodward–Fieser (dienes/enones), Scott (aryl carbonyls), aromatic bases.

Accuracy is that of textbook increment schemes (~±0.3 ppm 1H, ±3–5 ppm 13C), which is
adequate for cross-checking assignments but not a substitute for an experiment.  For a
data-driven prediction the UI links to nmrdb.org with the same SMILES.
"""
from __future__ import annotations
import math, itertools
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# =============================================================================== utilities
def _S(s):  # compiled SMARTS cache
    if s not in _S.cache:
        _S.cache[s] = Chem.MolFromSmarts(s)
    return _S.cache[s]
_S.cache = {}


def _neighbors_heavy(a):
    return [n for n in a.GetNeighbors() if n.GetAtomicNum() > 1]


def substituent_class(mol: Chem.Mol, attach: Chem.Atom, via: Chem.Atom) -> str:
    """
    Classify the group `via` seen from `attach` (via is directly bonded to attach).
    Returns a key used by the increment tables.
    """
    z = via.GetAtomicNum(); sym = via.GetSymbol()
    if z == 6:
        if via.GetIsAromatic():
            return "aryl"
        # carbonyl-type?
        dbl_o = [n for n in via.GetNeighbors() if n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(via.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2]
        if dbl_o:
            others = [n for n in via.GetNeighbors() if n.GetIdx() not in (attach.GetIdx(), dbl_o[0].GetIdx())]
            if not others:
                return "CHO"
            o = others[0]
            if o.GetAtomicNum() == 8:
                return "COOH" if o.GetTotalNumHs() > 0 or o.GetFormalCharge() < 0 else "COOR"
            if o.GetAtomicNum() == 7:
                return "CONR"
            if o.GetAtomicNum() in (17, 35):
                return "COCl"
            return "COR"
        for b in via.GetBonds():
            o = b.GetOtherAtom(via)
            if o.GetIdx() == attach.GetIdx():
                continue
            t = b.GetBondTypeAsDouble()
            if t == 3:
                return "CN" if o.GetAtomicNum() == 7 else "alkynyl"
            if t == 2:
                if o.GetAtomicNum() == 6:
                    return "vinyl"
                if o.GetAtomicNum() == 7:
                    return "C=N"
        # CF3?
        if sum(1 for n in via.GetNeighbors() if n.GetAtomicNum() == 9) == 3:
            return "CF3"
        # alkyl bearing heteroatom (CH2X)
        het = [n for n in via.GetNeighbors() if n.GetAtomicNum() in (7, 8, 9, 16, 17, 35, 53) and n.GetIdx() != attach.GetIdx()]
        if het:
            return "CH2X"
        return "alkyl"
    if z in (7, 8, 16) and mol.GetBondBetweenAtoms(attach.GetIdx(), via.GetIdx()).GetBondTypeAsDouble() >= 2:
        return "=X"
    if z == 8:
        others = [n for n in via.GetNeighbors() if n.GetIdx() != attach.GetIdx()]
        if via.GetFormalCharge() < 0:
            return "O-"
        if not others:
            return "OH"
        o = others[0]
        if o.GetAtomicNum() == 6:
            if any(n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(o.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2 for n in o.GetNeighbors()):
                return "OCOR"
            if o.GetIsAromatic():
                return "OAr"
            return "OR"
        if o.GetAtomicNum() == 7:
            return "ONO2" if via.GetFormalCharge() == 0 and o.GetFormalCharge() == 1 else "OR"
        return "OR"
    if z == 7:
        others = [n for n in via.GetNeighbors() if n.GetIdx() != attach.GetIdx()]
        if via.GetFormalCharge() == 1 and any(n.GetAtomicNum() == 8 for n in others):
            return "NO2"
        if via.GetIsAromatic():
            return "aryl"          # part of heteroaromatic ring
        if any(n.GetAtomicNum() == 6 and any(m.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(n.GetIdx(), m.GetIdx()).GetBondTypeAsDouble() == 2 for m in n.GetNeighbors()) for n in others):
            return "NHCOR"
        if via.GetFormalCharge() == 1:
            return "NR3+"
        if any(mol.GetBondBetweenAtoms(via.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2 for n in others):
            return "N=C"
        if via.GetTotalNumHs() == 2:
            return "NH2"
        return "NR2"
    if z == 16:
        if any(n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(via.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2 for n in via.GetNeighbors()):
            return "SO2"
        return "SR"
    return {9: "F", 17: "Cl", 35: "Br", 53: "I", 15: "P", 14: "Si", 5: "B"}.get(z, "X")


# =============================================================================== 1H NMR
# alpha increments (relative to an sp3-carbon substituent)
_H_ALPHA = {"alkyl": 0.0, "CH2X": 0.15, "aryl": 1.40, "vinyl": 0.80, "alkynyl": 0.95, "CHO": 1.35, "COR": 1.25, "COOH": 1.20, "COOR": 1.15,
            "CONR": 1.05, "COCl": 1.85, "CN": 1.15, "C=N": 0.9, "OH": 2.35, "OR": 2.35, "OAr": 2.95, "OCOR": 2.85, "O-": 2.0, "ONO2": 3.6,
            "NH2": 1.45, "NR2": 1.35, "NHCOR": 2.05, "NR3+": 2.4, "NO2": 3.25, "N=C": 2.2, "F": 3.05, "Cl": 2.25, "Br": 2.25, "I": 2.00,
            "SR": 1.35, "SO2": 1.85, "CF3": 0.9, "P": 0.9, "Si": -0.8, "B": 0.1, "X": 1.0}
_H_BETA = {"aryl": 0.25, "vinyl": 0.15, "alkynyl": 0.2, "CHO": 0.25, "COR": 0.25, "COOH": 0.25, "COOR": 0.25, "CONR": 0.2, "CN": 0.35,
           "OH": 0.30, "OR": 0.30, "OAr": 0.40, "OCOR": 0.40, "NH2": 0.15, "NR2": 0.15, "NHCOR": 0.3, "NR3+": 0.5, "NO2": 0.7, "F": 0.45,
           "Cl": 0.45, "Br": 0.55, "I": 0.6, "SR": 0.3, "SO2": 0.45, "CF3": 0.35, "alkyl": 0.0, "CH2X": 0.05}
_H_BASE = {3: 0.86, 2: 1.25, 1: 1.50, 0: 0}     # CH3-C, CH2(C)2, CH(C)3

# aromatic ring increments (ortho, meta, para) relative to 7.26
_H_AROM = {"alkyl": (-0.18, -0.10, -0.20), "CH2X": (-0.05, -0.05, -0.10), "aryl": (0.37, 0.20, 0.10), "vinyl": (0.06, -0.03, -0.10), "alkynyl": (0.15, -0.02, -0.01),
           "CHO": (0.56, 0.22, 0.29), "COR": (0.62, 0.14, 0.21), "COOH": (0.85, 0.18, 0.27), "COOR": (0.71, 0.10, 0.21), "CONR": (0.61, 0.10, 0.17), "COCl": (0.84, 0.22, 0.36),
           "CN": (0.36, 0.18, 0.28), "OH": (-0.56, -0.12, -0.45), "OR": (-0.48, -0.09, -0.44), "OAr": (-0.29, -0.05, -0.23), "OCOR": (-0.25, 0.03, -0.13), "O-": (-0.8, -0.2, -0.6),
           "NH2": (-0.75, -0.25, -0.65), "NR2": (-0.66, -0.18, -0.67), "NHCOR": (0.12, -0.07, -0.28), "NR3+": (0.69, 0.36, 0.31), "NO2": (0.95, 0.26, 0.38), "N=C": (0.2, 0.1, 0.1),
           "F": (-0.26, 0.00, -0.20), "Cl": (0.03, -0.02, -0.09), "Br": (0.18, -0.08, -0.04), "I": (0.39, -0.21, 0.00), "SR": (0.01, -0.15, -0.17), "SO2": (0.76, 0.35, 0.45),
           "CF3": (0.32, 0.14, 0.20), "P": (0.3, 0.1, 0.1), "Si": (0.2, 0.0, 0.0), "B": (0.3, 0.1, 0.1), "X": (0.1, 0.0, 0.0)}
# alkene increments (gem, cis, trans) relative to 5.25
_H_VINYL = {"alkyl": (0.45, -0.22, -0.28), "CH2X": (0.64, -0.01, -0.02), "aryl": (1.38, 0.36, -0.07), "vinyl": (1.24, 0.02, -0.05), "alkynyl": (0.47, 0.38, 0.12),
            "CHO": (1.02, 0.95, 1.17), "COR": (1.10, 1.12, 0.87), "COOH": (0.97, 1.41, 0.71), "COOR": (0.80, 1.18, 0.55), "CONR": (1.37, 0.98, 0.46), "COCl": (1.11, 1.46, 1.01),
            "CN": (0.27, 0.75, 0.55), "OH": (1.22, -1.07, -1.21), "OR": (1.22, -1.07, -1.21), "OAr": (1.22, -1.07, -1.21), "OCOR": (2.11, -0.35, -0.64),
            "NH2": (0.80, -1.26, -1.21), "NR2": (0.80, -1.26, -1.21), "NHCOR": (2.08, -0.57, -0.72), "NO2": (1.87, 1.30, 0.62), "F": (1.54, -0.40, -1.02),
            "Cl": (1.08, 0.18, 0.13), "Br": (1.07, 0.45, 0.55), "I": (1.14, 0.81, 0.88), "SR": (1.11, -0.29, -0.13), "SO2": (1.55, 1.16, 0.93), "CF3": (0.66, 0.61, 0.32),
            "X": (0.5, 0.2, 0.2)}

from . import solvents as SOLV
from . import hose as HOSE

def _nmr_solvent(name):
    """Map any solvent name (deuterated or not) to the parameters the NMR model needs."""
    key, sv = SOLV.get(name)
    return key, {"residual": sv["residual"], "c13": sv["c13"], "OH": sv["OH"], "phenol": sv["phenol"], "NH": sv["NH"], "amide": sv["amide"], "COOH": sv["COOH"], "exchange": sv["exchange"], "nmr": sv["nmr"], "et30": sv["et30"], "hbd": sv["hbd"], "hba": sv["hba"]}

SOLVENTS = {k: v["nmr"] for k, v in SOLV.SOLVENTS.items()}   # kept for API listing


def _ring_position_offsets(mol, ring_atoms, target_idx):
    """For aromatic 6-ring: map each other ring atom -> 'ortho'/'meta'/'para' relative to target."""
    n = len(ring_atoms); pos = ring_atoms.index(target_idx); out = {}
    for k, idx in enumerate(ring_atoms):
        d = min((k - pos) % n, (pos - k) % n)
        if d == 0:
            continue
        out[idx] = {1: 0, 2: 1, 3: 2}.get(d, 2)  # index into (o,m,p)
    return out


def _aromatic_H_shift(mol, atom):
    ri = mol.GetRingInfo()
    rings = [r for r in ri.AtomRings() if atom.GetIdx() in r and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    if not rings:
        return 7.26
    ring = min(rings, key=len)
    hetero = [mol.GetAtomWithIdx(i) for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() != 6]
    # base by heteroaromatic position
    base = 7.26
    if hetero:
        # distance to nearest heteroatom along the ring
        n = len(ring); pos = ring.index(atom.GetIdx())
        dmin = min(min((ring.index(h.GetIdx()) - pos) % n, (pos - ring.index(h.GetIdx())) % n) for h in hetero)
        h0 = hetero[0]
        if n == 6:   # pyridine-like
            base = {1: 8.55, 2: 7.30, 3: 7.70}.get(dmin, 7.4)
            if len(hetero) >= 2:
                base += 0.4
        else:        # 5-membered
            if h0.GetAtomicNum() == 8:
                base = 7.40 if dmin == 1 else 6.35
            elif h0.GetAtomicNum() == 16:
                base = 7.30 if dmin == 1 else 7.10
            else:
                base = 6.70 if dmin == 1 else 6.20
                if h0.GetTotalNumHs() == 0 and h0.GetDegree() == 2:  # pyridine-type N in azole
                    base += 1.0
    # fused ring: naphthalene-like alpha/beta positions
    nfused = sum(1 for i in ring if ri.NumAtomRings(i) > 1)
    if nfused >= 2:
        nb_fused = any(ri.NumAtomRings(n.GetIdx()) > 1 for n in atom.GetNeighbors())
        base += 0.55 if nb_fused else 0.20
    # substituent increments
    offs = _ring_position_offsets(mol, list(ring), atom.GetIdx())
    delta = 0.0
    for idx, k in offs.items():
        ra = mol.GetAtomWithIdx(idx)
        for nb in ra.GetNeighbors():
            if nb.GetIdx() in ring or nb.GetAtomicNum() == 1:
                continue
            if ri.NumAtomRings(idx) > 1 and nb.GetIsAromatic():
                continue   # fused ring atom handled above
            cls = substituent_class(mol, ra, nb)
            delta += _H_AROM.get(cls, (0.1, 0, 0))[k]
    return base + delta


def _vinyl_H_shift(mol, atom):
    # find C=C partner
    partner = None
    for b in atom.GetBonds():
        if b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(atom).GetAtomicNum() == 6 and not b.GetIsAromatic():
            partner = b.GetOtherAtom(atom)
    if partner is None:
        return 5.25
    delta = 0.0
    for nb in atom.GetNeighbors():
        if nb.GetIdx() == partner.GetIdx() or nb.GetAtomicNum() == 1:
            continue
        delta += _H_VINYL.get(substituent_class(mol, atom, nb), _H_VINYL["X"])[0]
    # cis / trans from bond stereo (unknown -> average)
    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), partner.GetIdx())
    for nb in partner.GetNeighbors():
        if nb.GetIdx() == atom.GetIdx() or nb.GetAtomicNum() == 1:
            continue
        inc = _H_VINYL.get(substituent_class(mol, partner, nb), _H_VINYL["X"])
        rel = _cis_trans(mol, bond, atom, nb)
        delta += inc[1] if rel == "cis" else inc[2] if rel == "trans" else (inc[1] + inc[2]) / 2
    return 5.25 + delta


def _cis_trans(mol, bond, h_carbon, sub_on_partner):
    """Relationship between the H on h_carbon and substituent on partner across bond (uses E/Z)."""
    st = bond.GetStereo()
    if st in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY):
        return None
    sa = list(bond.GetStereoAtoms())
    if len(sa) != 2:
        return None
    begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
    # reference atoms: sa[0] on begin, sa[1] on end
    ref_h_side = sa[0] if h_carbon.GetIdx() == begin.GetIdx() else sa[1]
    ref_p_side = sa[1] if h_carbon.GetIdx() == begin.GetIdx() else sa[0]
    # heavy substituent on h_carbon (other than partner): is it the ref?
    hc_subs = [n.GetIdx() for n in h_carbon.GetNeighbors() if n.GetIdx() not in (begin.GetIdx(), end.GetIdx())]
    if not hc_subs:
        return None
    same_side_refs = st in (Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOZ)
    hc_sub_is_ref = hc_subs[0] == ref_h_side
    p_sub_is_ref = sub_on_partner.GetIdx() == ref_p_side
    subs_same_side = same_side_refs == (hc_sub_is_ref == p_sub_is_ref)
    # H is opposite to the heavy substituent on its carbon
    return "trans" if subs_same_side else "cis"


def _sp3_H_shift(mol, atom):
    """delta = base(CH3 0.86 / CH2 1.25 / CH 1.50, i.e. all-alkyl) + sum(alpha_X) + sum(beta_X)."""
    nH = atom.GetTotalNumHs()
    heavy = _neighbors_heavy(atom)
    base = {3: 0.86, 2: 1.25, 1: 1.50}.get(nH, 1.5)
    if not heavy:
        return 0.23
    delta = 0.0
    for n in heavy:
        cls = substituent_class(mol, atom, n)
        if cls in ("alkyl", "CH2X"):
            pass                                   # alkyl neighbour: already in base
        else:
            delta += _H_ALPHA.get(cls, 1.0)
        # beta substituents seen through a carbon neighbour
        if n.GetAtomicNum() == 6 and not n.GetIsAromatic() and cls in ("alkyl", "CH2X"):
            for nn in _neighbors_heavy(n):
                if nn.GetIdx() == atom.GetIdx():
                    continue
                c2 = substituent_class(mol, n, nn)
                if c2 in ("alkyl", "CH2X"):
                    continue
                delta += _H_BETA.get(c2, 0.1) * (0.7 if nH == 3 else 1.0)
    ri = mol.GetRingInfo()
    if ri.NumAtomRings(atom.GetIdx()):
        rs = min(len(r) for r in ri.AtomRings() if atom.GetIdx() in r)
        delta += {3: -0.9, 4: 0.7, 5: 0.25, 6: 0.15}.get(rs, 0.2)
    return base + delta


_GIAO_CACHE: dict = {}


def _giao_for(mol):
    from rdkit import Chem as _C
    key = _C.MolToSmiles(mol)
    if key not in _GIAO_CACHE:
        from . import giao as GIAO
        _GIAO_CACHE[key] = GIAO.giao_shifts(mol)
    return _GIAO_CACHE[key]


def h_environments(mol: Chem.Mol, solvent="chloroform", level="hose"):
    """
    Returns list of dicts: {atoms:[heavy idx...], nH, shift, type, label, couplings:[(J, n)]}
    grouped by topological symmetry.  Solvent: exchangeable-proton positions from the table,
    plus a small polarity/H-bond-acceptor shift on protons alpha to O/N (DMSO, DMF ~ +0.1 ppm).
    """
    _, sol = _nmr_solvent(solvent)
    polar_shift = 0.12 * max(0.0, sol["hba"] - 0.4)      # 0 in CDCl3, ~0.04 DMSO
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    groups = defaultdict(list)
    for a in mol.GetAtoms():
        if a.GetTotalNumHs() == 0:
            continue
        groups[ranks[a.GetIdx()]].append(a.GetIdx())
    envs = []
    for rk, idxs in groups.items():
        a = mol.GetAtomWithIdx(idxs[0]); z = a.GetAtomicNum(); nH = a.GetTotalNumHs()
        typ = ""; shift = None; exch = False
        if z == 6:
            if a.GetIsAromatic():
                shift = _aromatic_H_shift(mol, a); typ = "aromatic C–H"
            elif any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(a).GetAtomicNum() == 6 for b in a.GetBonds()):
                shift = _vinyl_H_shift(mol, a); typ = "vinylic C–H"
            elif any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(a).GetAtomicNum() == 8 for b in a.GetBonds()):
                arom = any(n.GetIsAromatic() for n in a.GetNeighbors()); shift = 10.0 if arom else 9.75; typ = "aldehyde C–H"
                if any(n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1 for n in a.GetNeighbors()):
                    shift = 8.05; typ = "formate C–H"
                if any(n.GetAtomicNum() == 7 for n in a.GetNeighbors()):
                    shift = 8.1; typ = "formamide C–H"
            elif any(b.GetBondTypeAsDouble() == 3 for b in a.GetBonds()):
                shift = 2.0 + (0.9 if any(n.GetIsAromatic() for nb in a.GetNeighbors() for n in nb.GetNeighbors()) else 0); typ = "alkyne C–H"
            elif any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(a).GetAtomicNum() == 7 for b in a.GetBonds()):
                shift = 8.0; typ = "imine C–H"
            else:
                shift = _sp3_H_shift(mol, a); typ = {3: "CH₃", 2: "CH₂", 1: "CH"}[nH]
                if any(n.GetAtomicNum() in (7, 8) for n in a.GetNeighbors()):
                    shift += polar_shift
        elif z == 8:
            exch = True
            nb = _neighbors_heavy(a)[0] if _neighbors_heavy(a) else None
            if nb is not None and nb.GetAtomicNum() == 6 and any(m.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(nb.GetIdx(), m.GetIdx()).GetBondTypeAsDouble() == 2 for m in nb.GetNeighbors()):
                shift = sol["COOH"]; typ = "COOH"
            elif nb is not None and nb.GetIsAromatic():
                shift = sol["phenol"]; typ = "phenol OH"
            else:
                shift = sol["OH"]; typ = "alcohol OH"
        elif z == 7:
            exch = True
            nbs = _neighbors_heavy(a)
            if any(n.GetAtomicNum() == 6 and any(m.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(n.GetIdx(), m.GetIdx()).GetBondTypeAsDouble() == 2 for m in n.GetNeighbors()) for n in nbs):
                shift = sol["amide"]; typ = "amide N–H"
            elif a.GetIsAromatic():
                shift = (8.0 if sol["NH"] is not None else None); typ = "aromatic N–H"
            elif a.GetFormalCharge() == 1:
                shift = (8.5 if sol["NH"] is not None else None); typ = "N⁺–H"
            elif any(n.GetIsAromatic() for n in nbs):
                shift = (3.7 if sol["NH"] is not None else None); typ = "aniline N–H"
            else:
                shift = sol["NH"]; typ = "amine N–H"
        elif z == 16:
            exch = True; shift = 1.6 if sol["NH"] is not None else None; typ = "S–H"
        else:
            shift = 4.0; typ = f"{a.GetSymbol()}–H"
        if shift is None:
            continue   # exchanged with deuterated solvent
        src = {"source": "increments", "sigma": None, "n_ref": None, "depth": None}
        if z == 6 and level == "hose":
            hm = HOSE.load()
            if hm is not None:
                pr = hm.predict(mol, a.GetIdx(), "1H")
                if pr and pr[3] >= 3 and pr[2] >= 2:
                    shift = pr[0]; src = {"source": "nmrshiftdb2", "sigma": pr[1], "n_ref": pr[2], "depth": pr[3]}
        elif z == 6 and level == "giao":
            g = _giao_for(mol)
            if "h1" in g and a.GetIdx() in g["h1"]:
                shift = round(g["h1"][a.GetIdx()], 2); src = {"source": "giao", "sigma": g["mae_h"], "n_ref": None, "depth": None}
        envs.append({"atoms": idxs, "nH": nH * len(idxs), "nH_each": nH, "shift": round(shift, 2), "type": typ, "exchangeable": exch, **src})
    # couplings (first order, 3J H–H; 2J for alkene gem; 4J meta aromatic)
    idx2env = {i: e for e in envs for i in e["atoms"]}
    for e in envs:
        a = mol.GetAtomWithIdx(e["atoms"][0]); cps = []
        if e["exchangeable"]:
            e["couplings"] = []; continue
        seen_env = set()
        for nb in a.GetNeighbors():
            if nb.GetIdx() in e["atoms"] or nb.GetIdx() not in idx2env:
                continue
            o = idx2env[nb.GetIdx()]
            if o["exchangeable"] or id(o) in seen_env:
                continue
            seen_env.add(id(o))
            bond = mol.GetBondBetweenAtoms(a.GetIdx(), nb.GetIdx())
            if a.GetIsAromatic() and nb.GetIsAromatic():
                J = 8.0 if len([r for r in mol.GetRingInfo().AtomRings() if a.GetIdx() in r and nb.GetIdx() in r and len(r) == 6]) else 3.5
            elif bond.GetBondTypeAsDouble() == 2:
                # vinyl: cis 10 / trans 16 depending on stereo, unknown -> 12
                rel = None
                if nb.GetTotalNumHs() and a.GetTotalNumHs():
                    subs = [n for n in nb.GetNeighbors() if n.GetIdx() != a.GetIdx() and n.GetAtomicNum() > 1]
                    if subs:
                        r = _cis_trans(mol, bond, a, subs[0])
                        rel = {"cis": "trans", "trans": "cis"}.get(r)
                    else:
                        rel = "both"
                J = {"cis": 10.5, "trans": 16.0}.get(rel, 12.0)
                if rel == "both":
                    cps.append((10.5, 1)); cps.append((16.0, 1)); continue
            elif a.GetIsAromatic() or nb.GetIsAromatic():
                J = 7.0 if (nb.GetIsAromatic() and not a.GetIsAromatic()) else 0.0
                J = 0.0
            elif any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(a).GetAtomicNum() == 8 for b in a.GetBonds()):
                J = 2.0   # aldehyde CHO–CH
            elif any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(nb).GetAtomicNum() == 8 for b in nb.GetBonds()):
                J = 2.0
            elif bond.GetBondTypeAsDouble() == 1 and (any(b.GetBondTypeAsDouble() == 2 for b in a.GetBonds()) or any(b.GetBondTypeAsDouble() == 2 for b in nb.GetBonds())):
                J = 6.5   # allylic vicinal
            else:
                J = 7.0
            if J > 0:
                cps.append((J, o["nH_each"] * sum(1 for x in a.GetNeighbors() if x.GetIdx() in o["atoms"])))
        # meta coupling in aromatic rings (4J ~ 2 Hz) -> small, ignore for clarity
        # merge identical J's
        merged = defaultdict(int)
        for J, n in cps:
            merged[J] += n
        e["couplings"] = sorted(merged.items(), reverse=True)
        e["multiplicity"] = _mult_label(e["couplings"])
    for e in envs:
        if e["exchangeable"]:
            e["multiplicity"] = "br s"
    return envs


def _mult_label(cps):
    names = {1: "s", 2: "d", 3: "t", 4: "q", 5: "quint", 6: "sext", 7: "sept"}
    if not cps:
        return "s"
    parts = [names.get(n + 1, "m") for J, n in cps]
    if len(parts) > 2 or "m" in parts:
        return "m"
    return "".join(parts)


def _multiplet(center_hz, cps):
    lines = [(center_hz, 1.0)]
    for J, n in cps:
        new = []
        coeffs = [math.comb(n, k) for k in range(n + 1)]
        for pos, amp in lines:
            for k, c in enumerate(coeffs):
                new.append((pos + (k - n / 2) * J, amp * c / sum(coeffs)))
        lines = new
    return lines


def h_spectrum(mol, solvent="chloroform", mhz=400.0, npoints=6000, lw_hz=0.9, level="hose"):
    envs = h_environments(mol, solvent, level)
    key, sol = _nmr_solvent(solvent)
    lo, hi = -0.5, 13.0
    xs = [lo + (hi - lo) * i / (npoints - 1) for i in range(npoints)]
    ys = [0.0] * npoints
    peaks = []
    for e in envs:
        lines = _multiplet(e["shift"] * mhz, e["couplings"])
        lw = lw_hz if not e["exchangeable"] else 12.0
        for pos, amp in lines:
            ppm = pos / mhz; A = amp * e["nH"]
            hw = lw / mhz
            for i, x in enumerate(xs):
                dx = x - ppm
                if abs(dx) < 40 * hw:
                    ys[i] += A * hw * hw / (dx * dx + hw * hw)
        peaks.append({**e, "J": [round(J, 1) for J, n in e["couplings"]]})
    return {"x": [round(v, 4) for v in xs], "y": [round(v, 4) for v in ys], "peaks": peaks, "solvent": key, "nmr_solvent": sol["nmr"], "mhz": mhz,
            "residual_solvent_ppm": sol["residual"], "unit": "ppm",
            "solvent_note": ("exchangeable protons replaced by deuterium (not observed)" if sol["exchange"] else "exchangeable protons shown at typical positions; they broaden and shift with concentration/water") + ("" if sol["nmr"] else "; no common deuterated form of this solvent; shifts shown as in a non-polar reference")}


# =============================================================================== 13C NMR
_C_INC = {  # alpha, beta, gamma for substituent replacing H on sp3 carbon
    "OH": (49, 10, -6), "OR": (58, 8, -4), "OAr": (56, 6, -4), "OCOR": (51, 6, -3), "O-": (60, 8, -4), "ONO2": (62, 4, -4),
    "NH2": (28.3, 11.3, -5.1), "NR2": (40, 5, -5), "NHCOR": (28, 7, -4), "NR3+": (42, 4, -4), "NO2": (61.6, 3, -4.6), "N=C": (30, 6, -4),
    "F": (70, 8, -7), "Cl": (31, 11, -4), "Br": (20, 11, -3), "I": (-6, 11, -1), "SR": (11, 12, -3), "SO2": (35, 3, -3),
    "CN": (4, 3, -3), "COR": (30, 1, -2), "CHO": (31, 0, -2), "COOH": (21, 3, -2), "COOR": (20, 3, -2), "CONR": (22, 2.5, -3), "COCl": (33, 2, -2),
    "aryl": (23, 9, -2), "vinyl": (21.5, 6.9, -2.1), "alkynyl": (4.5, 5.4, -3.5), "C=N": (20, 5, -2), "CF3": (25, 0, -2), "P": (10, 5, -2), "Si": (-5, 5, -2), "B": (10, 3, -2), "X": (10, 5, -2),
}
_C_AROM = {  # ipso, ortho, meta, para (benzene 128.5)
    "alkyl": (9.3, 0.7, -0.1, -2.9), "CH2X": (12.0, -1.0, 0.0, -2.5), "aryl": (13.1, -1.1, 0.4, -1.2), "vinyl": (9.1, -2.4, 0.2, -0.5), "alkynyl": (-6.1, 3.8, 0.4, -0.2),
    "F": (34.8, -12.9, 1.4, -4.5), "Cl": (6.2, 0.4, 1.3, -1.9), "Br": (-5.5, 3.4, 1.7, -1.6), "I": (-34.1, 8.7, 1.6, -1.1),
    "OH": (26.9, -12.7, 1.4, -7.3), "OR": (31.4, -14.4, 1.0, -7.7), "OAr": (29.1, -9.5, 0.3, -5.3), "OCOR": (22.4, -7.1, 0.4, -3.2), "O-": (39.6, -8.2, 1.9, -13.6),
    "NH2": (18.2, -13.4, 0.8, -10.0), "NR2": (22.4, -15.7, 0.8, -11.8), "NHCOR": (9.7, -8.1, 0.2, -4.4), "NR3+": (19.5, -6.8, 2.5, 4.0), "NO2": (19.6, -5.3, 0.9, 6.0), "N=C": (10, -2, 1, 2),
    "CN": (-15.4, 3.6, 0.6, 3.9), "CHO": (8.2, 1.2, 0.6, 5.8), "COR": (9.1, 0.1, 0.0, 4.2), "COOH": (2.1, 1.5, 0.0, 5.1), "COOR": (2.0, 1.2, -0.1, 4.8), "CONR": (5.0, -1.2, 0.1, 3.4), "COCl": (4.6, 2.9, 0.6, 7.0),
    "CF3": (2.6, -3.1, 0.4, 3.4), "SR": (10.2, -1.8, 0.2, -3.6), "SO2": (12.0, -1.0, 1.0, 5.0), "P": (5, 3, 0, 1), "Si": (12, 5, 0, 0), "B": (10, 5, 0, 3), "X": (5, 0, 0, 0),
}
_C_VINYL = {  # (alpha, beta) increments on C=C relative to 123.3
    "alkyl": (10.6, -7.9), "CH2X": (10, -7), "aryl": (12.5, -11.0), "vinyl": (13.6, -7.0), "alkynyl": (5, 3), "OH": (29, -39), "OR": (29, -39), "OAr": (28, -37), "OCOR": (18, -27),
    "NH2": (28, -32), "NR2": (28, -32), "NHCOR": (20, -25), "NO2": (22, -1), "F": (24.9, -34.3), "Cl": (2.6, -6.1), "Br": (-7.9, -1.4), "I": (-38, 7), "SR": (9, -13),
    "SO2": (14, 7), "CN": (-15.1, 14.2), "CHO": (13.1, 12.7), "COR": (15, 5.8), "COOH": (4.2, 8.9), "COOR": (6.3, 7.0), "CONR": (5, 7), "COCl": (8, 14), "CF3": (-3, 8), "X": (5, 0),
}
_GP_STERIC = {("1", "3"): -1.1, ("1", "4"): -3.4, ("2", "3"): -2.5, ("2", "4"): -7.2, ("3", "2"): -3.7, ("3", "3"): -9.5, ("3", "4"): -15.0,
              ("4", "1"): -1.5, ("4", "2"): -8.4, ("4", "3"): -15.0, ("4", "4"): -25.0}


def _carbon_degree(a):
    return str(min(4, sum(1 for n in a.GetNeighbors() if n.GetAtomicNum() == 6)))


def _sp3_C_shift(mol, atom):
    """
    Grant–Paul alkane term over sp3 carbons reached through alkyl chains, plus one
    alpha/beta/gamma increment per functional substituent (not traversed through).
    """
    i = atom.GetIdx()
    counts = defaultdict(int)          # distance -> number of alkyl carbons
    shift = -2.3
    visited = {i}
    frontier = [(atom, 0)]
    while frontier:
        cur, d = frontier.pop()
        for nb in _neighbors_heavy(cur):
            j = nb.GetIdx()
            if j in visited:
                continue
            cls = substituent_class(mol, cur, nb)
            if cls in ("alkyl", "CH2X") and nb.GetAtomicNum() == 6:
                visited.add(j)
                if d + 1 <= 4:
                    counts[d + 1] += 1
                    frontier.append((nb, d + 1))
            else:
                if cls == "=X":
                    continue
                visited.add(j)
                if d + 1 <= 3:
                    shift += _C_INC.get(cls, (5, 3, -2))[d]
    shift += 9.1 * counts[1] + 9.4 * counts[2] - 2.5 * counts[3] + 0.3 * counts[4]
    mine = _carbon_degree(atom)
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() == 6 and not n.GetIsAromatic() and substituent_class(mol, atom, n) in ("alkyl", "CH2X"):
            shift += _GP_STERIC.get((mine, _carbon_degree(n)), 0.0)
    ri = mol.GetRingInfo()
    if ri.NumAtomRings(i):
        rs = min(len(r) for r in ri.AtomRings() if i in r)
        shift += {3: -19.0, 4: -12.0, 5: -9.0, 6: -5.0, 7: -4.0}.get(rs, -3.0)
    return shift


def _aromatic_C_shift(mol, atom):
    ri = mol.GetRingInfo()
    rings = [r for r in ri.AtomRings() if atom.GetIdx() in r and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    ring = min(rings, key=len)
    n = len(ring); pos = ring.index(atom.GetIdx())
    hetero = [mol.GetAtomWithIdx(i) for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() != 6]
    base = 128.5
    if hetero:
        h0 = hetero[0]; dmin = min(min((ring.index(h.GetIdx()) - pos) % n, (pos - ring.index(h.GetIdx())) % n) for h in hetero)
        if n == 6:
            base = {1: 150.0, 2: 124.0, 3: 136.0}.get(dmin, 130)
            if len(hetero) >= 2:
                base += 6
        else:
            if h0.GetAtomicNum() == 8:
                base = 143.0 if dmin == 1 else 110.0
            elif h0.GetAtomicNum() == 16:
                base = 125.0 if dmin == 1 else 127.0
            else:
                base = 118.0 if dmin == 1 else 108.0
                if len(hetero) >= 2:
                    base += 15 if dmin == 1 else 10
    if ri.NumAtomRings(atom.GetIdx()) > 1:
        base += 5.0     # ring-fusion carbon (naphthalene C4a ~133.5)
    elif sum(1 for i in ring if ri.NumAtomRings(i) > 1) >= 2:
        nb_fused = any(ri.NumAtomRings(x.GetIdx()) > 1 for x in atom.GetNeighbors())
        base += -0.5 if nb_fused else -2.5
    delta = 0.0
    for k, idx in enumerate(ring):
        d = min((k - pos) % n, (pos - k) % n)
        if d > 3:
            continue
        ra = mol.GetAtomWithIdx(idx)
        for nb in ra.GetNeighbors():
            if nb.GetIdx() in ring or nb.GetAtomicNum() == 1:
                continue
            if ri.NumAtomRings(idx) > 1 and nb.GetIsAromatic():
                continue
            cls = substituent_class(mol, ra, nb)
            delta += _C_AROM.get(cls, (5, 0, 0, 0))[d]
    return base + delta


def _vinyl_C_shift(mol, atom):
    partner = None
    for b in atom.GetBonds():
        if b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(atom).GetAtomicNum() == 6:
            partner = b.GetOtherAtom(atom)
    shift = 123.3
    for nb in atom.GetNeighbors():
        if partner is not None and nb.GetIdx() == partner.GetIdx() or nb.GetAtomicNum() == 1:
            continue
        shift += _C_VINYL.get(substituent_class(mol, atom, nb), (5, 0))[0]
    if partner is not None:
        for nb in partner.GetNeighbors():
            if nb.GetIdx() == atom.GetIdx() or nb.GetAtomicNum() == 1:
                continue
            shift += _C_VINYL.get(substituent_class(mol, partner, nb), (5, 0))[1]
    return shift


def _carbonyl_C_shift(mol, atom):
    o = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2]
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != o[0].GetIdx()]
    conj = any(n.GetIsAromatic() or any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(n).GetAtomicNum() == 6 for b in n.GetBonds()) for n in others if n.GetAtomicNum() == 6)
    hetero = [n for n in others if n.GetAtomicNum() != 6]
    nH = atom.GetTotalNumHs()
    if not hetero:
        if nH:
            return 192.0 if conj else 201.0
        return 198.0 if conj else 208.0
    kinds = sorted(n.GetAtomicNum() for n in hetero)
    if kinds == [8]:
        oa = hetero[0]
        if oa.GetTotalNumHs() or oa.GetFormalCharge() < 0:
            return 172.0 if conj else 178.5
        # ester / anhydride / carbonate
        other_c = [n for n in oa.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
        if other_c and other_c[0].GetAtomicNum() == 6 and any(m.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(other_c[0].GetIdx(), m.GetIdx()).GetBondTypeAsDouble() == 2 for m in other_c[0].GetNeighbors()):
            return 168.0
        if nH:
            return 161.0
        return 166.5 if conj else 173.5
    if kinds == [7]:
        return 167.5 if conj else 172.5
    if kinds == [7, 7]:
        return 158.0
    if kinds == [7, 8]:
        return 156.0
    if kinds == [8, 8]:
        return 155.0
    if 17 in kinds or 35 in kinds:
        return 168.0 if conj else 173.0
    if 16 in kinds:
        return 196.0
    return 170.0


def c_environments(mol: Chem.Mol, level="hose"):
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    groups = defaultdict(list)
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 6:
            groups[ranks[a.GetIdx()]].append(a.GetIdx())
    envs = []
    for rk, idxs in groups.items():
        a = mol.GetAtomWithIdx(idxs[0])
        if a.GetIsAromatic():
            s = _aromatic_C_shift(mol, a); t = "aromatic"
        elif any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(a).GetAtomicNum() == 8 for b in a.GetBonds()):
            s = _carbonyl_C_shift(mol, a); t = "C=O"
        elif any(b.GetBondTypeAsDouble() == 3 for b in a.GetBonds()):
            tb = [b for b in a.GetBonds() if b.GetBondTypeAsDouble() == 3][0]
            if tb.GetOtherAtom(a).GetAtomicNum() == 7:
                s = 118.0 + (1.0 if any(n.GetIsAromatic() for n in a.GetNeighbors()) else 0); t = "C≡N"
            else:
                s = 68.0 if a.GetTotalNumHs() else 82.0; t = "C≡C"
                if any(n.GetIsAromatic() for n in a.GetNeighbors()):
                    s += 4
        elif any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(a).GetAtomicNum() == 6 for b in a.GetBonds()):
            s = _vinyl_C_shift(mol, a); t = "C=C"
        elif any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(a).GetAtomicNum() == 7 for b in a.GetBonds()):
            s = 160.0; t = "C=N"
        else:
            s = _sp3_C_shift(mol, a); t = "sp³"
        nH = a.GetTotalNumHs()
        src = {"source": "increments", "sigma": None, "n_ref": None, "depth": None}
        hm = HOSE.load() if level == "hose" else None
        if hm is not None:
            pr = hm.predict(mol, a.GetIdx(), "13C")
            if pr and pr[3] >= 2:
                s = pr[0]; src = {"source": "nmrshiftdb2", "sigma": pr[1], "n_ref": pr[2], "depth": pr[3]}
        if level == "giao":
            g = _giao_for(mol)
            if "c13" in g and a.GetIdx() in g["c13"]:
                s = round(g["c13"][a.GetIdx()], 1); src = {"source": "giao", "sigma": g["mae_c"], "n_ref": None, "depth": None}
            elif "error" in g:
                src = {"source": "increments", "sigma": None, "n_ref": None, "depth": None, "giao_error": g["error"]}
        envs.append({"atoms": idxs, "shift": round(s, 1), "type": t, "nH": nH, "dept": {0: "C", 1: "CH", 2: "CH₂", 3: "CH₃"}[min(nH, 3)], "n": len(idxs), **src})
    envs.sort(key=lambda e: -e["shift"])
    return envs


def c_spectrum(mol, solvent="chloroform", mhz=100.0, npoints=6000, lw_hz=1.5, level="hose"):
    envs = c_environments(mol, level)
    lo, hi = -10.0, 230.0
    xs = [lo + (hi - lo) * i / (npoints - 1) for i in range(npoints)]
    ys = [0.0] * npoints
    for e in envs:
        # NOE-ish intensity: quaternary carbons weaker
        A = e["n"] * (0.35 if e["nH"] == 0 else 1.0)
        hw = lw_hz / mhz
        for i, x in enumerate(xs):
            dx = x - e["shift"]
            if abs(dx) < 60 * hw:
                ys[i] += A * hw * hw / (dx * dx + hw * hw)
    key, sol = _nmr_solvent(solvent)
    return {"x": [round(v, 3) for v in xs], "y": [round(v, 4) for v in ys], "peaks": envs, "solvent": key, "nmr_solvent": sol["nmr"], "mhz": mhz,
            "solvent_c13_ppm": sol["c13"], "unit": "ppm"}


# =============================================================================== IR
_IR_BANDS = [  # (SMARTS, [(center, width, intensity, label)], per-match)
    ("[CX4;!$(C=*)][#1]", [(2925, 40, 0.75, "C–H stretch (sp³)"), (2870, 30, 0.5, "C–H stretch (sp³, sym)"), (1460, 25, 0.4, "C–H bend (CH₂/CH₃)")]),
    ("[CH3]", [(1375, 15, 0.35, "CH₃ symmetric bend (umbrella)")]),
    ("[c][#1]", [(3050, 30, 0.3, "C–H stretch (aromatic)")]),
    ("[C;!c]=[C;!c][#1]", [(3080, 30, 0.3, "=C–H stretch (vinylic)")]),
    ("[C]#[CH]", [(3300, 20, 0.6, "≡C–H stretch"), (2120, 15, 0.3, "C≡C stretch (terminal)"), (640, 30, 0.6, "≡C–H bend")]),
    ("[C;!$(C#[CH])]#[C;!$(C#[CH])]", [(2220, 15, 0.15, "C≡C stretch (internal, weak)")]),
    ("[C]#N", [(2245, 15, 0.55, "C≡N stretch")]),
    ("[c]C#N", [(2228, 15, 0.6, "Ar–C≡N stretch")]),
    ("[OX2H][CX4]", [(3350, 150, 0.85, "O–H stretch (alcohol, H-bonded, broad)"), (1060, 40, 0.7, "C–O stretch (alcohol)")]),
    ("[OX2H]c", [(3350, 160, 0.8, "O–H stretch (phenol, broad)"), (1220, 40, 0.75, "C–O stretch (phenol)")]),
    ("[CX3](=O)[OX2H1]", [(3000, 350, 0.7, "O–H stretch (carboxylic acid, very broad)"), (1290, 40, 0.6, "C–O stretch (acid)"), (930, 60, 0.4, "O–H out-of-plane bend (acid dimer)")]),
    ("[CX3;!$(Cc);!$(CC=C)](=O)[OX2H1]", [(1710, 20, 0.95, "C=O stretch (carboxylic acid dimer)")]),
    ("[c,$(C=C)][CX3](=O)[OX2H1]", [(1690, 20, 0.95, "C=O stretch (aryl acid, conjugated)")]),
    ("[CX3;!$(C[N,O,S,Cl,Br]);!$(Cc);!$(CC=C);!$([CH1])](=O)[#6]", [(1715, 18, 0.95, "C=O stretch (ketone)")]),
    ("[CX3;!$(C[N,O,S,Cl,Br]);!$([CH1])](=O)[c,$(C=C)]", [(1685, 18, 0.95, "C=O stretch (aryl ketone)")]),
    ("[CX3H1;!$(Cc);!$(CC=C)](=O)[#6]", [(1725, 18, 0.95, "C=O stretch (aldehyde)"), (2820, 20, 0.35, "C–H stretch (aldehyde)"), (2720, 20, 0.35, "C–H stretch (aldehyde, Fermi doublet)")]),
    ("[CX3H1](=O)[c,$(C=C)]", [(1700, 18, 0.95, "C=O stretch (aryl aldehyde)")]),
    ("[CX3;!$(Cc);!$(CC=C)](=O)[OX2][#6;!c]", [(1740, 18, 0.95, "C=O stretch (ester)"), (1240, 40, 0.8, "C–O stretch (ester, asym)"), (1100, 40, 0.5, "C–O stretch (ester)")]),
    ("[c,$(C=C)][CX3](=O)[OX2][#6;!c]", [(1720, 18, 0.95, "C=O stretch (aryl ester, conjugated)")]),
    ("[CX3](=O)[OX2]c", [(1765, 18, 0.95, "C=O stretch (phenyl ester)")]),
    ("[CX3;!$(C(=O)(N)N);!$(C(=O)(N)O)](=O)[NX3]", [(1655, 25, 0.95, "C=O stretch (amide I)"), (1550, 30, 0.6, "N–H bend / C–N (amide II)")]),
    ("[CX3](=O)[NX3H2]", [(3350, 40, 0.5, "N–H stretch (1° amide, asym)"), (3180, 40, 0.5, "N–H stretch (1° amide, sym)")]),
    ("[CX3](=O)[NX3H1]", [(3300, 50, 0.5, "N–H stretch (2° amide)")]),
    ("[CX3](=O)Cl", [(1800, 18, 0.95, "C=O stretch (acid chloride)")]),
    ("[CX3](=O)O[CX3]=O", [(1820, 18, 0.8, "C=O stretch (anhydride, asym)"), (1760, 18, 0.9, "C=O stretch (anhydride, sym)")]),
    ("[NX3;H2;!$(NC=O)]", [(3400, 40, 0.4, "N–H stretch (1° amine, asym)"), (3320, 40, 0.4, "N–H stretch (1° amine, sym)"), (1610, 30, 0.4, "N–H scissor"), (800, 80, 0.5, "N–H wag (broad)")]),
    ("[NX3;H1;!$(NC=O)]", [(3350, 40, 0.3, "N–H stretch (2° amine)")]),
    ("[CX4][NX3;!$(NC=O)]", [(1120, 40, 0.5, "C–N stretch (aliphatic amine)")]),
    ("c[NX3;!$(NC=O)]", [(1300, 40, 0.7, "C–N stretch (aryl amine)")]),
    ("[C;!c]=[C;!c]", [(1650, 15, 0.3, "C=C stretch")]),
    ("c", [(1600, 15, 0.4, "C=C ring stretch (aromatic)"), (1495, 15, 0.5, "C=C ring stretch (aromatic)"), (1450, 15, 0.35, "C=C ring stretch (aromatic)")]),
    ("c1ccccc1[!#1]", [(750, 20, 0.7, "C–H out-of-plane bend (monosubst.)"), (700, 20, 0.7, "ring bend (monosubst.)")]),
    ("[!#1]c1ccc([!#1])cc1", [(830, 20, 0.7, "C–H out-of-plane bend (para)")]),
    ("[!#1]c1ccccc1[!#1]", [(755, 20, 0.7, "C–H out-of-plane bend (ortho)")]),
    ("[!#1]c1cccc([!#1])c1", [(780, 20, 0.6, "C–H out-of-plane bend (meta)"), (880, 20, 0.4, "C–H out-of-plane bend (meta, isolated H)")]),
    ("[CX4][OX2][CX4]", [(1110, 40, 0.8, "C–O–C stretch (ether)")]),
    ("c[OX2][CX4]", [(1250, 30, 0.8, "C–O–C stretch (aryl ether, asym)"), (1040, 30, 0.6, "C–O–C stretch (aryl ether, sym)")]),
    ("[NX3+](=O)[O-]", [(1530, 20, 0.9, "N=O stretch (nitro, asym)"), (1350, 20, 0.9, "N=O stretch (nitro, sym)")]),
    ("[SX3](=O)", [(1050, 30, 0.8, "S=O stretch (sulfoxide)")]),
    ("[SX4](=O)(=O)", [(1330, 25, 0.9, "S=O stretch (sulfone/sulfonyl, asym)"), (1150, 25, 0.9, "S=O stretch (sulfone/sulfonyl, sym)")]),
    ("[SX2H]", [(2560, 20, 0.15, "S–H stretch")]),
    ("[#6]F", [(1150, 60, 0.85, "C–F stretch")]),
    ("[#6]Cl", [(720, 50, 0.6, "C–Cl stretch")]),
    ("[#6]Br", [(600, 40, 0.5, "C–Br stretch")]),
    ("[#6]I", [(520, 40, 0.4, "C–I stretch")]),
    ("C=N", [(1650, 20, 0.5, "C=N stretch (imine)")]),
    ("[NX3][CX3]=[NX2]", [(1640, 20, 0.6, "C=N stretch (amidine)")]),
    ("C1CC1", [(3080, 20, 0.2, "C–H stretch (cyclopropane)"), (1020, 20, 0.4, "ring breathing (cyclopropane)")]),
]


def ir_spectrum(mol: Chem.Mol, npoints=1800, solvent="chloroform"):
    """
    Group-frequency simulation.  Solvent handling (empirical, per Bellamy / Socrates):
      C=O stretch   : gas → non-polar → polar aprotic → protic shifts down ~5–25 cm⁻¹
      O–H / N–H     : H-bond acceptor solvents broaden and lower free X–H stretches
    """
    key, sv = SOLV.get(solvent)
    co_shift = -(6.0 * (sv["et30"] - 31.0) / 10.0 + 12.0 * sv["hbd"])          # ~0 hexane, -5 CHCl3, -8 DMSO, -20 MeOH
    xh_shift = -(60.0 * sv["hba"] + 40.0 * sv["hbd"])                          # H-bond acceptors pull X–H down
    xh_broad = 1.0 + 1.2 * sv["hba"] + 0.8 * sv["hbd"]
    if key == "gas phase":
        co_shift, xh_shift, xh_broad = +15.0, +120.0, 0.35                       # free molecules: sharp, higher X–H
    mh = Chem.AddHs(mol)
    bands = []
    for smarts, blist in _IR_BANDS:
        patt = _S(smarts)
        if patt is None:
            continue
        n = len(mh.GetSubstructMatches(patt))
        if n == 0:
            continue
        for center, width, inten, label in blist:
            c, w = center, width
            if "C=O" in label:
                c += co_shift
            elif ("O–H stretch" in label or "N–H stretch" in label) and center > 3000:
                c += xh_shift; w *= xh_broad
            bands.append({"center": round(c), "width": round(w), "intensity": min(1.0, inten * (1 + 0.15 * math.log(n))), "label": label, "count": n})
    # merge duplicate labels (e.g. aromatic ring bands counted once)
    seen = {}
    for b in bands:
        lab = b["label"]
        if lab not in seen or b["intensity"] > seen[lab]["intensity"]:
            seen[lab] = b
    bands = sorted(seen.values(), key=lambda b: -b["center"])
    xs = [4000 - (4000 - 400) * i / (npoints - 1) for i in range(npoints)]
    absorb = [0.0] * npoints
    for b in bands:
        for i, x in enumerate(xs):
            dx = x - b["center"]
            absorb[i] += b["intensity"] * math.exp(-(dx * dx) / (2 * b["width"] * b["width"]))
    trans = [round(100 * math.exp(-1.6 * a), 2) for a in absorb]
    # fingerprint noise-free baseline: add generic weak skeletal bands
    return {"x": [round(v, 1) for v in xs], "y": trans, "bands": bands, "unit": "cm⁻¹", "y_label": "% Transmittance", "solvent": key,
            "solvent_note": f"{key}: C=O shifted {co_shift:+.0f} cm⁻¹, free X–H stretches {xh_shift:+.0f} cm⁻¹ and ×{xh_broad:.1f} width relative to the gas-phase/non-polar table values.",
            "note": "Group-frequency simulation; band positions ±20 cm⁻¹, intensities qualitative."}


# =============================================================================== UV-Vis
def _longest_conjugated_path(mol):
    """Count C=C/C=O/aromatic units in the largest conjugated system (crude)."""
    conj_bonds = [b for b in mol.GetBonds() if b.GetIsConjugated()]
    if not conj_bonds:
        return 0, set()
    # union-find over atoms in conjugated bonds
    parent = {}
    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for b in conj_bonds:
        a, c = find(b.GetBeginAtomIdx()), find(b.GetEndAtomIdx())
        if a != c:
            parent[a] = c
    systems = defaultdict(set)
    for b in conj_bonds:
        systems[find(b.GetBeginAtomIdx())].update((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
    best = max(systems.values(), key=len)
    ndouble = sum(1 for b in mol.GetBonds() if b.GetBeginAtomIdx() in best and b.GetEndAtomIdx() in best and (b.GetBondTypeAsDouble() == 2 or b.GetIsAromatic()))
    return ndouble, best


def uv_vis(mol: Chem.Mol, npoints=700, solvent="chloroform"):
    """
    Empirical λmax rules plus a linear solvatochromic correction in Reichardt E_T(30):
      π→π* bands red-shift ~ +0.35 nm per kcal/mol of E_T(30) above hexane (31.0)
      n→π* bands blue-shift with polarity and strongly with H-bond donors (≈ −0.6 nm/kcal − 12 nm·α)
    (Reichardt, Chem. Rev. 1994, 94, 2319; Kosower).  Tabulated rule values refer to ethanol/methanol
    for Woodward–Fieser (dienes/enones) and are re-based accordingly.
    """
    key, sv = SOLV.get(solvent)
    dpol = sv["et30"] - 31.0
    shift_pipi = 0.35 * dpol
    shift_npi = -0.6 * dpol - 12.0 * sv["hbd"]
    wf_rebase = -0.35 * (51.9 - 31.0)     # rules tabulated in ethanol -> remove that solvent's contribution first
    peaks = []; notes = []
    ri = mol.GetRingInfo()
    n_arom = rdMolDescriptors.CalcNumAromaticRings(mol)
    # --- aromatic chromophores
    if n_arom:
        # fused systems
        fused = 0
        for r in ri.AtomRings():
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r) and sum(1 for i in r if ri.NumAtomRings(i) > 1) >= 2:
                fused += 1
        if fused >= 3:
            peaks.append((375, 8000, "π→π* (anthracene-type fused aromatic)")); peaks.append((255, 150000, "π→π* (fused aromatic, intense)"))
        elif fused >= 1:
            peaks.append((312, 250, "π→π* (naphthalene-type, ¹Lb)")); peaks.append((275, 5600, "π→π* (naphthalene-type, ¹La)")); peaks.append((221, 100000, "π→π* (¹Bb)"))
        else:
            # Scott rules for Ar–C=O, else generic benzene with auxochromes
            arcarbonyl = mol.GetSubstructMatches(_S("c[CX3]=O"))
            base = 254.0; eps = 200; kind = "benzenoid ¹Lb"
            sub_shift = 0.0
            if arcarbonyl:
                ar_c = mol.GetAtomWithIdx(arcarbonyl[0][0]); co = mol.GetAtomWithIdx(arcarbonyl[0][1])
                others = [n for n in co.GetNeighbors() if n.GetIdx() != ar_c.GetIdx() and n.GetAtomicNum() != 8 or (n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(co.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1)]
                if co.GetTotalNumHs():
                    base = 250; kind = "Ar–CHO (Scott)"
                elif any(n.GetAtomicNum() in (7, 8) for n in others):
                    base = 230; kind = "Ar–COOR/COOH/CONR (Scott)"
                else:
                    base = 246; kind = "Ar–COR (Scott)"
                eps = 12000
                ring = [r for r in ri.AtomRings() if ar_c.GetIdx() in r][0]
                offs = _ring_position_offsets(mol, list(ring), ar_c.GetIdx())
                scott = {"alkyl": (3, 3, 10), "CH2X": (3, 3, 10), "OH": (7, 7, 25), "OR": (7, 7, 25), "OAr": (7, 7, 25), "O-": (11, 20, 78), "Cl": (0, 0, 10), "Br": (2, 2, 15),
                         "NH2": (13, 13, 58), "NHCOR": (20, 20, 45), "NR2": (20, 20, 85), "OCOR": (2, 2, 5)}
                for idx, k in offs.items():
                    ra = mol.GetAtomWithIdx(idx)
                    for nb in ra.GetNeighbors():
                        if nb.GetIdx() in ring or nb.GetAtomicNum() == 1:
                            continue
                        sub_shift += scott.get(substituent_class(mol, ra, nb), (0, 0, 0))[k]
                peaks.append((base + sub_shift, eps, f"{kind} + substituents ({sub_shift:+.0f} nm)"))
                peaks.append((base + sub_shift + 70, 60, "n→π* (aryl carbonyl)"))
            else:
                # generic benzene with auxochromes (bathochromic shifts)
                aux = {"OH": 16, "OR": 16, "O-": 33, "NH2": 26, "NR2": 30, "NHCOR": 12, "Cl": 10, "Br": 9, "alkyl": 7, "CH2X": 7, "vinyl": 28, "aryl": 20, "CN": 20, "NO2": 60, "COOH": 18, "COOR": 18, "CONR": 14, "CHO": 30, "COR": 30, "SR": 30, "F": 6, "I": 8}
                shift = 0.0; found = []
                for r in ri.AtomRings():
                    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r):
                        continue
                    for i in r:
                        ra = mol.GetAtomWithIdx(i)
                        for nb in ra.GetNeighbors():
                            if nb.GetIdx() in r or nb.GetAtomicNum() == 1:
                                continue
                            cls = substituent_class(mol, ra, nb); v = aux.get(cls, 5)
                            shift += v; found.append(cls)
                    break
                shift = min(shift, 90)
                hetero_ring = any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for r in ri.AtomRings() for i in r if mol.GetAtomWithIdx(i).GetIsAromatic())
                if hetero_ring:
                    base = 257 if any(len(r) == 6 for r in ri.AtomRings()) else 240
                eps = 200 + 40 * shift
                peaks.append((base + shift, eps, f"π→π* ({'heteroaromatic' if hetero_ring else 'benzenoid'} secondary band; auxochromes {'+'.join(found) or 'none'} {shift:+.0f} nm)"))
                peaks.append((203 + shift * 0.6, 7000 + 100 * shift, "π→π* (primary band)"))
                if any(c in ("NO2", "CHO", "COR", "COOH", "COOR") for c in found):
                    peaks.append((base + shift + 25, 100, "n→π* (weak)"))
    # --- conjugated dienes / polyenes / enones (Woodward–Fieser)
    ndb, system = _longest_conjugated_path(mol)
    non_arom_double = [b for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 2 and not b.GetIsAromatic() and b.GetBeginAtomIdx() in system and b.GetEndAtomIdx() in system]
    cc = [b for b in non_arom_double if b.GetBeginAtom().GetAtomicNum() == 6 and b.GetEndAtom().GetAtomicNum() == 6]
    co = [b for b in non_arom_double if 8 in (b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum())]
    enone = mol.GetSubstructMatches(_S("[C;!c]=[C;!c][CX3]=O"))
    if enone:
        c_b, c_a, c_co, o = enone[0]
        co_atom = mol.GetAtomWithIdx(c_co)
        others = [n for n in co_atom.GetNeighbors() if n.GetIdx() not in (c_a, o)]
        if co_atom.GetTotalNumHs():
            base = 210; kind = "enal"
        elif others and others[0].GetAtomicNum() == 8:
            base = 195; kind = "enoic acid/ester"
        else:
            ring_sizes = [len(r) for r in ri.AtomRings() if c_co in r and c_a in r]
            base = 202 if 5 in ring_sizes else 215; kind = "enone"
        add = 0; parts = []
        for pos, idx, name in ((1, c_a, "α"), (2, c_b, "β")):
            at = mol.GetAtomWithIdx(idx)
            for nb in at.GetNeighbors():
                if nb.GetIdx() in (c_a, c_b, c_co) or nb.GetAtomicNum() == 1:
                    continue
                cls = substituent_class(mol, at, nb)
                inc = {"alkyl": {1: 10, 2: 12}, "CH2X": {1: 10, 2: 12}, "OH": {1: 35, 2: 30}, "OR": {1: 35, 2: 30}, "OCOR": {1: 6, 2: 6}, "Cl": {1: 15, 2: 12}, "Br": {1: 25, 2: 30},
                       "NR2": {1: 0, 2: 95}, "NH2": {1: 0, 2: 95}, "SR": {1: 0, 2: 85}, "aryl": {1: 0, 2: 30}, "vinyl": {1: 30, 2: 30}}.get(cls, {1: 5, 2: 5})[pos]
                add += inc; parts.append(f"{name}-{cls} +{inc}")
        # extended conjugation
        extra = max(0, len(cc) - 1)
        if extra:
            add += 30 * extra; parts.append(f"extended conjugation +{30 * extra}")
        # exocyclic double bond
        if ri.NumAtomRings(c_co) == 0 and ri.NumAtomRings(c_a) and not ri.NumAtomRings(c_b):
            pass
        lam = base + add
        peaks.append((lam, 10000, f"π→π* ({kind}, Woodward–Fieser: {base} {'+ '.join([''] + parts) if parts else ''})"))
        peaks.append((lam + 100, 50, "n→π* (enone, weak)"))
    elif len(cc) >= 2:
        # diene
        base = 217; parts = []
        homoannular = any(all(mol.GetBondBetweenAtoms(*sorted((r[i], r[(i + 1) % len(r)]))[:2]) is not None for i in range(len(r))) and sum(1 for b in cc if b.GetBeginAtomIdx() in r and b.GetEndAtomIdx() in r) >= 2 for r in ri.AtomRings())
        if homoannular:
            base = 253; parts.append("homoannular")
        elif any(ri.NumAtomRings(b.GetBeginAtomIdx()) for b in cc):
            base = 214; parts.append("heteroannular/ring")
        add = 30 * (len(cc) - 2)
        if add:
            parts.append(f"extra C=C +{add}")
        diene_atoms = set(i for b in cc for i in (b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        nalk = 0
        for i in diene_atoms:
            for nb in mol.GetAtomWithIdx(i).GetNeighbors():
                if nb.GetIdx() not in diene_atoms and nb.GetAtomicNum() > 1:
                    cls = substituent_class(mol, mol.GetAtomWithIdx(i), nb)
                    inc = {"alkyl": 5, "CH2X": 5, "OR": 6, "OH": 6, "OCOR": 0, "Cl": 5, "Br": 5, "SR": 30, "NR2": 60, "NH2": 60, "aryl": 5}.get(cls, 5)
                    add += inc; nalk += 1
        if nalk:
            parts.append(f"{nalk} substituent(s)")
        peaks.append((base + add, 10000 + 5000 * (len(cc) - 2), f"π→π* (polyene, Woodward–Fieser: {base} {'; '.join(parts)})"))
    elif len(cc) == 1 and not n_arom:
        peaks.append((180, 10000, "π→π* (isolated C=C; below solvent cutoff)"))
    # --- isolated carbonyl / other n→π*
    if not enone and not any("carbonyl" in p[2] for p in peaks):
        if mol.HasSubstructMatch(_S("[CX3;!$(Cc)](=O)[#6,#1]")):
            peaks.append((280, 15, "n→π* (isolated C=O, weak)"))
        elif mol.HasSubstructMatch(_S("[CX3;!$(Cc)](=O)[N,O]")):
            peaks.append((210, 50, "n→π* (ester/amide, weak, below 220 nm)"))
    if mol.HasSubstructMatch(_S("[NX3+](=O)[O-]")) and not n_arom:
        peaks.append((275, 20, "n→π* (nitro)"))
    if mol.HasSubstructMatch(_S("N=N")):
        peaks.append((350, 30, "n→π* (azo)"))
    if mol.HasSubstructMatch(_S("[SX2][SX2]")):
        peaks.append((250, 300, "n→σ* (disulfide)"))
    if not peaks:
        peaks.append((None, None, "No chromophore absorbing above ~200 nm (saturated compound)"))
        notes.append("Compound is expected to be transparent in the near-UV/visible region.")
    # solvatochromic correction
    corr = []
    for lam, eps, lab in peaks:
        if lam is None:
            corr.append((lam, eps, lab)); continue
        if "n→π*" in lab or "n→σ*" in lab:
            corr.append((lam + shift_npi, eps, lab))
        elif "Woodward" in lab:
            corr.append((lam + wf_rebase + shift_pipi, eps, lab))
        else:
            corr.append((lam + shift_pipi, eps, lab))
    peaks = corr
    # simulated spectrum
    xs = [190 + (700 - 190) * i / (npoints - 1) for i in range(npoints)]
    ys = [0.0] * npoints
    for lam, eps, lab in peaks:
        if lam is None:
            continue
        w = 18 + 0.05 * lam
        for i, x in enumerate(xs):
            ys[i] += eps * math.exp(-((x - lam) ** 2) / (2 * w * w))
    lam_max = max((p for p in peaks if p[0]), key=lambda p: p[0], default=None)
    colour = None
    if lam_max and lam_max[0] > 380:
        colour = "likely coloured (absorbs in visible)"
    elif lam_max and lam_max[0] > 340:
        colour = "possibly pale yellow (absorption tail into visible)"
    return {"x": [round(v, 1) for v in xs], "y": [round(v, 1) for v in ys], "peaks": [{"lambda_max": p[0] and round(p[0]), "epsilon": p[1], "assignment": p[2]} for p in peaks],
            "colour": colour, "notes": notes, "unit": "nm", "y_label": "ε (L mol⁻¹ cm⁻¹, approx.)", "solvent": key,
            "solvent_note": f"{key} (E_T(30) = {sv['et30']}): π→π* bands shifted {shift_pipi:+.0f} nm, n→π* bands {shift_npi:+.0f} nm relative to hexane.",
            "note": "Woodward–Fieser / Scott empirical rules; ±5–10 nm for the chromophores covered, qualitative otherwise."}


def all_spectra(mol, solvent="chloroform", h_mhz=400.0, nmr_level="hose"):
    """nmr_level: 'hose' (nmrshiftdb2 experimental-environment prediction, increments as fallback) or 'increments' (pure additivity)."""
    key, _ = SOLV.get(solvent)
    giao_note = None
    if nmr_level == "giao":
        g = _giao_for(mol)
        giao_note = g.get("error") or f"GIAO HF/6-31G* (scaled), SCF energy {g['e_scf']:.4f} Eh; measured fit MAE 13C 1.9 ppm / 1H 0.18 ppm"
    return {"h1": h_spectrum(mol, key, h_mhz, level=nmr_level), "c13": c_spectrum(mol, key, h_mhz / 3.976, level=nmr_level), "nmr_level": nmr_level, "giao_note": giao_note, "ir": ir_spectrum(mol, solvent=key), "uv": uv_vis(mol, solvent=key), "solvent": SOLV.summary(key)}
