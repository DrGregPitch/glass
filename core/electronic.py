"""
Electronic structure at the GFN2-xTB level (tblite): molecular-orbital energy levels,
HOMO/LUMO gap, dipole moment and per-atom Mulliken charges.

Used for the MO level diagram and for colouring molecular surfaces (an approximate
electrostatic view: vdW/SAS surface coloured by the partial charge of the nearest atom).
Limits: tight-binding orbital energies are systematically compressed vs hybrid DFT
(gaps too small by ~2-3 eV for organics); charges are Mulliken-type; no solvent.
"""
from __future__ import annotations
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

HARTREE_EV = 27.211386245988
BOHR = 0.529177210903


def electronic_structure(mol: Chem.Mol, max_atoms=80, seed=0xF00D) -> dict:
    try:
        from tblite.interface import Calculator
    except Exception:
        return {"error": "tblite (GFN2-xTB) not installed"}
    if mol.GetNumAtoms() > max_atoms:
        return {"error": f"limited to {max_atoms} heavy atoms"}
    mh = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed
    if AllChem.EmbedMolecule(mh, ps) < 0:
        return {"error": "3D embedding failed"}
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=2000)
    except Exception:
        pass
    Z = np.array([a.GetAtomicNum() for a in mh.GetAtoms()])
    X = np.array(mh.GetConformer().GetPositions()) / BOHR
    try:
        calc = Calculator("GFN2-xTB", Z, X, charge=Chem.GetFormalCharge(mh))
        calc.set("verbosity", 0)
        res = calc.singlepoint()
    except Exception as e:
        return {"error": f"GFN2-xTB failed ({type(e).__name__}: {e})"}
    e_orb = np.array(res.get("orbital-energies")) * HARTREE_EV
    occ = np.array(res.get("orbital-occupations"))
    charges = np.array(res.get("charges"))
    dip = np.array(res.get("dipole"))
    homo = int(np.max(np.where(occ > 0.5))) if np.any(occ > 0.5) else -1
    lumo = homo + 1 if homo + 1 < len(e_orb) else -1
    lo = max(0, homo - 7); hi = min(len(e_orb), (lumo if lumo >= 0 else homo) + 8)
    levels = [{"i": int(i), "e_ev": round(float(e_orb[i]), 3), "occ": round(float(occ[i]), 2),
               "label": ("HOMO" if i == homo else "LUMO" if i == lumo else (f"HOMO−{homo - i}" if i < homo else f"LUMO+{i - lumo}"))}
              for i in range(lo, hi)]
    masses = None
    nuc_center = (X * Z[:, None]).sum(axis=0) / Z.sum() * BOHR      # nuclear charge centre, Å
    return {
        "level": "GFN2-xTB", "n_orbitals": int(len(e_orb)), "n_electrons": int(round(float(occ.sum()))),
        "homo_ev": round(float(e_orb[homo]), 3) if homo >= 0 else None,
        "lumo_ev": round(float(e_orb[lumo]), 3) if lumo >= 0 else None,
        "gap_ev": round(float(e_orb[lumo] - e_orb[homo]), 3) if homo >= 0 and lumo >= 0 else None,
        "dipole_debye": round(float(np.linalg.norm(dip)) * 2.541746, 3),
        "dipole_vec_debye": [round(float(v) * 2.541746, 4) for v in dip],
        "nuc_center": [round(float(v), 4) for v in nuc_center],
        "levels": levels,
        "charges": [round(float(c), 4) for c in charges],
        "symbols": [a.GetSymbol() for a in mh.GetAtoms()],
        "xyz": np.round(np.array(mh.GetConformer().GetPositions()), 4).tolist(),
        "molblock": Chem.MolToMolBlock(mh),
        "note": "GFN2-xTB (tight-binding DFT): orbital-energy gaps are compressed vs hybrid DFT/experiment (typically 2–3 eV too small for organics); charges are Mulliken-type partial charges; gas phase, single conformer.",
    }
