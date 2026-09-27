"""
GIAO NMR chemical shifts at the HF/6-31G* level (PySCF), empirically scaled.

delta = slope * sigma_iso + intercept, with slope/intercept fitted against experimental
CDCl3 shifts of an 11-molecule reference set computed at exactly this level and geometry
convention (ETKDG + MMFF94s).  Fit quality (leave-in):
    13C:  delta = -1.00074 sigma + 202.093   n=14  MAE 1.92 ppm  max 4.59  R^2 = 0.9988
    1H :  delta = -0.94733 sigma + 31.211    n=11  MAE 0.18 ppm  max 0.42  R^2 = 0.9935
Scope limits (enforced/stated): closed-shell molecules, <=16 heavy atoms (cost grows
steeply), gas-phase single conformer; heavy halogens (Br, I) attached to carbon carry
large relativistic errors at this level and are refused.

Includes a runtime fix for a pyscf-properties CPHF batching bug (krylov may pass a
stacked (nvec, nmo*nocc) array; the shipped code reshapes with a hard-coded 3).
"""
from __future__ import annotations
import inspect, textwrap
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

SLOPE_C, INT_C, MAE_C = -1.00074, 202.093, 1.92
SLOPE_H, INT_H, MAE_H = -0.94733, 31.211, 0.18
MAX_HEAVY = 16
ALLOWED_Z = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17}     # Br/I refused (relativistic heavy-atom errors)

_patched = False


def _patch_cphf():
    global _patched
    if _patched:
        return
    from pyscf.prop.nmr import rhf as nmr_rhf
    src = inspect.getsource(nmr_rhf.gen_vind)
    if "reshape(3,nmo,nocc)" in src:
        ns = {}
        exec(compile(textwrap.dedent(src.replace("reshape(3,nmo,nocc)", "reshape(-1,nmo,nocc)")), "<patched gen_vind>", "exec"), vars(nmr_rhf), ns)
        nmr_rhf.gen_vind = ns["gen_vind"]
    _patched = True


def giao_shifts(mol: Chem.Mol, seed=0xF00D) -> dict:
    """Per-atom scaled shifts: {'c13': {atom_idx: delta}, 'h1': {heavy_parent_idx: delta}, 'meta': ...}."""
    if mol.GetNumHeavyAtoms() > MAX_HEAVY:
        return {"error": f"GIAO level limited to {MAX_HEAVY} heavy atoms (cost); use the database or increment level"}
    if any(a.GetAtomicNum() not in ALLOWED_Z for a in mol.GetAtoms()):
        return {"error": "GIAO level refuses Br/I and metals (large relativistic shift errors at HF/6-31G*)"}
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return {"error": "open-shell species not supported"}
    from pyscf import gto, scf
    from pyscf.prop import nmr as pnmr
    _patch_cphf()
    mh = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed
    if AllChem.EmbedMolecule(mh, ps) < 0:
        return {"error": "3D embedding failed"}
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=2000)
    except Exception:
        pass
    conf = mh.GetConformer()
    charge = Chem.GetFormalCharge(mh)
    nelec = sum(a.GetAtomicNum() for a in mh.GetAtoms()) - charge
    if nelec % 2:
        return {"error": "odd-electron system"}
    m = gto.M(atom=[(a.GetSymbol(), tuple(conf.GetAtomPosition(i))) for i, a in enumerate(mh.GetAtoms())],
              basis="6-31g*", charge=charge, verbose=0)
    mf = scf.RHF(m).run()
    if not mf.converged:
        return {"error": "SCF did not converge"}
    iso = [float(np.trace(s) / 3) for s in pnmr.RHF(mf).kernel()]
    c13, h1 = {}, {}
    hcount = {}
    for i, a in enumerate(mh.GetAtoms()):
        if a.GetAtomicNum() == 6:
            c13[i] = SLOPE_C * iso[i] + INT_C
        elif a.GetAtomicNum() == 1:
            p = a.GetNeighbors()[0].GetIdx()
            h1.setdefault(p, []).append(SLOPE_H * iso[i] + INT_H)
    h1 = {p: float(np.mean(v)) for p, v in h1.items()}      # average over the H's on one carbon (rotamer average)
    return {"c13": c13, "h1": h1, "level": "GIAO HF/6-31G* (scaled)", "mae_c": MAE_C, "mae_h": MAE_H,
            "e_scf": float(mf.e_tot)}
