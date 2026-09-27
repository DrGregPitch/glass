"""
HOMO/LUMO isosurface cubes at the HF/STO-3G level (PySCF).

Qualitative by design: minimal-basis Hartree–Fock reproduces the nodal structure and
symmetry of frontier orbitals well, which is what an isosurface shows; energies and
diffuse character are NOT reliable at this level (documented in the methods registry).
Geometry: ETKDG + MMFF94s (same convention as the rest of the app).
"""
from __future__ import annotations
import threading
from rdkit import Chem
from rdkit.Chem import AllChem

_SCF_CACHE: dict = {}          # smiles -> (pyscf mol, converged mf, molblock)  (SCF done once per molecule)
_lock = threading.Lock()


def orbital_cube(mol: Chem.Mol, which: str = "homo", offset: int = 0, max_heavy=30, seed=0xF00D, npts=44) -> dict:
    """which='homo'|'lumo', offset>=0 counts away from the gap (homo,1 -> HOMO-1; lumo,2 -> LUMO+2).
    Returns {cube, energy_ev, occ, label, centroid, nuc_center, orb_dipole_debye}; centroid is the
    electron-position expectation <phi|r|phi> (Å); orb_dipole is the orbital's electronic dipole
    contribution -n_e*e*(<r> - nuclear charge centre) in Debye components."""
    if mol.GetNumHeavyAtoms() > max_heavy:
        raise ValueError(f"orbital surfaces limited to {max_heavy} heavy atoms")
    from pyscf import gto, scf
    from pyscf.tools import cubegen
    import tempfile, os
    key = Chem.MolToSmiles(mol)
    with _lock:
        cached = _SCF_CACHE.get(key)
    if cached is None:
        mh = Chem.AddHs(mol)
        ps = AllChem.ETKDGv3(); ps.randomSeed = seed
        if AllChem.EmbedMolecule(mh, ps) < 0:
            raise ValueError("3D embedding failed")
        try:
            AllChem.MMFFOptimizeMolecule(mh, maxIters=2000)
        except Exception:
            pass
        conf = mh.GetConformer()
        atoms = [(a.GetSymbol(), tuple(conf.GetAtomPosition(i))) for i, a in enumerate(mh.GetAtoms())]
        charge = Chem.GetFormalCharge(mh)
        nelec = sum(a.GetAtomicNum() for a in mh.GetAtoms()) - charge
        if nelec % 2:
            raise ValueError("open-shell species not supported for orbital surfaces")
        m = gto.M(atom=atoms, basis="sto-3g", charge=charge, verbose=0)
        mf = scf.RHF(m).run()
        if not mf.converged:
            raise ValueError("HF did not converge")
        with _lock:
            _SCF_CACHE[key] = (m, mf, Chem.MolToMolBlock(mh))
            if len(_SCF_CACHE) > 24:
                _SCF_CACHE.pop(next(iter(_SCF_CACHE)))
        cached = _SCF_CACHE[key]
    m, mf, molblock = cached
    import numpy as np
    homo_idx = m.nelectron // 2 - 1
    idx = homo_idx - offset if which == "homo" else homo_idx + 1 + offset
    if idx < 0 or idx >= mf.mo_coeff.shape[1]:
        raise ValueError(f"orbital {which}{'-' if which == 'homo' else '+'}{offset} does not exist at this basis size")
    label = ("HOMO" if offset == 0 else f"HOMO−{offset}") if which == "homo" else ("LUMO" if offset == 0 else f"LUMO+{offset}")
    occ = 2 if which == "homo" else 0
    HARTREE_EV = 27.211386245988; BOHR = 0.529177210903; AU_D = 2.541746
    c = mf.mo_coeff[:, idx]
    with m.with_common_orig((0, 0, 0)):
        r_ints = m.intor("int1e_r")                        # (3, nao, nao), Bohr
    centroid_bohr = np.array([c @ r_ints[k] @ c for k in range(3)])
    zs = m.atom_charges(); coords = m.atom_coords()
    nuc_center_bohr = (coords * zs[:, None]).sum(axis=0) / zs.sum()
    orb_dip_D = [-occ * float(v) * AU_D for v in (centroid_bohr - nuc_center_bohr)] if occ else [0.0, 0.0, 0.0]
    with tempfile.NamedTemporaryFile(suffix=".cube", delete=False) as f:
        path = f.name
    try:
        cubegen.orbital(m, path, c, nx=npts, ny=npts, nz=npts)
        cube = open(path).read()
    finally:
        os.unlink(path)
    return {"cube": cube, "label": label, "energy_ev": round(float(mf.mo_energy[idx]) * HARTREE_EV, 3), "occ": occ,
            "centroid": [round(float(v) * BOHR, 4) for v in centroid_bohr],
            "nuc_center": [round(float(v) * BOHR, 4) for v in nuc_center_bohr],
            "orb_dipole_debye": [round(v, 3) for v in orb_dip_D],
            "molblock": molblock}
