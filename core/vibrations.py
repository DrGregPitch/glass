"""
Harmonic normal-mode analysis on the MMFF94s force field.

Hessian by central finite differences of the analytic MMFF gradient, mass-weighted,
diagonalised; the six (five for linear) lowest modes are translations/rotations and are
dropped.  IR intensities are approximated from MMFF partial charges (fixed-charge dipole
derivative), which gives qualitatively right relative intensities.

Accuracy: MMFF harmonic frequencies are typically within ~5–10 % of experiment for
stretches (C=O ~ +50 cm⁻¹ too high, X–H ~ +150 cm⁻¹ too high, as with any unscaled
harmonic model); a uniform scale factor of 0.95 is applied and reported.
"""
from __future__ import annotations
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

SCALE = 0.95
AMU = 1.66053906660e-27
C_CM = 2.99792458e10
KCAL_A2 = 4184.0 / 6.02214076e23 * 1e20     # kcal/mol/Å² -> J/m²


def _mode_label(mh: Chem.Mol, disp: np.ndarray, xyz: np.ndarray) -> str:
    """Describe a mode by the internal coordinate that changes most."""
    best, label = 0.0, "skeletal deformation"
    for b in mh.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        v = xyz[j] - xyz[i]; L = np.linalg.norm(v); u = v / L
        stretch = abs(np.dot(disp[j] - disp[i], u)) / L
        if stretch > best:
            best = stretch
            a, c = mh.GetAtomWithIdx(i), mh.GetAtomWithIdx(j)
            order = {1.0: "–", 1.5: "⋯", 2.0: "=", 3.0: "≡"}.get(b.GetBondTypeAsDouble(), "–")
            label = f"{a.GetSymbol()}{order}{c.GetSymbol()} stretch"
    # compare with bending: if total displacement is mostly perpendicular to bonds -> bend
    total = np.linalg.norm(disp)
    if best * 1.0 < 0.15 * total / max(1, len(xyz)) ** 0.5:
        heavy_h = any(mh.GetAtomWithIdx(int(k)).GetAtomicNum() == 1 for k in np.argsort(-np.linalg.norm(disp, axis=1))[:2])
        label = "X–H bend" if heavy_h else "skeletal bend / torsion"
    return label


def normal_modes(mol: Chem.Mol, max_atoms=60, seed=0xF00D) -> dict:
    if mol.GetNumAtoms() > max_atoms:
        return {"error": f"normal-mode analysis limited to {max_atoms} heavy atoms"}
    mh = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed
    if AllChem.EmbedMolecule(mh, ps) < 0:
        return {"error": "3D embedding failed"}
    props = AllChem.MMFFGetMoleculeProperties(mh, mmffVariant="MMFF94s")
    if props is None:
        return {"error": "MMFF94s parameters unavailable for this molecule"}
    ff = AllChem.MMFFGetMoleculeForceField(mh, props)
    ff.Minimize(maxIts=5000)
    ff.Minimize(maxIts=5000)
    n = mh.GetNumAtoms()
    x0 = np.array(ff.Positions(), dtype=float)
    h = 1e-3
    H = np.zeros((3 * n, 3 * n))
    for k in range(3 * n):
        xp = x0.copy(); xp[k] += h
        xm = x0.copy(); xm[k] -= h
        ff.CalcEnergy(xp.tolist()); gp = np.array(ff.CalcGrad(xp.tolist()))   # CalcEnergy refreshes cached geometry
        ff.CalcEnergy(xm.tolist()); gm = np.array(ff.CalcGrad(xm.tolist()))
        H[k] = (gp - gm) / (2 * h)
    H = 0.5 * (H + H.T) * KCAL_A2                         # J/m²
    masses = np.array([a.GetMass() for a in mh.GetAtoms()]) * AMU
    m3 = np.repeat(masses, 3)
    Hw = H / np.sqrt(np.outer(m3, m3))
    w2, V = np.linalg.eigh(Hw)
    freqs = np.sign(w2) * np.sqrt(np.abs(w2)) / (2 * np.pi * C_CM)   # cm⁻¹
    charges = np.array([props.GetMMFFPartialCharge(i) for i in range(n)])
    xyz = x0.reshape(n, 3)
    order = np.argsort(np.abs(freqs))
    ntr = 5 if _is_linear(xyz) else 6
    modes = []
    for idx in order[ntr:]:
        f = float(freqs[idx])
        if f < 0:
            continue
        cart = (V[:, idx] / np.sqrt(m3)).reshape(n, 3)
        cart /= np.linalg.norm(cart)
        dmu = np.sum(charges[:, None] * cart, axis=0)
        inten = float(np.dot(dmu, dmu))
        modes.append({"freq": round(f * SCALE, 1), "freq_unscaled": round(f, 1), "intensity": inten,
                      "label": _mode_label(mh, cart, xyz), "vec": np.round(cart, 4).tolist()})
    if modes:
        mx = max(m["intensity"] for m in modes) or 1.0
        for m in modes:
            m["intensity"] = round(m["intensity"] / mx, 3)
    modes.sort(key=lambda m: -m["freq"])
    return {"symbols": [a.GetSymbol() for a in mh.GetAtoms()], "xyz": np.round(xyz, 4).tolist(),
            "bonds": [[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in mh.GetBonds()],
            "modes": modes, "scale": SCALE, "method": "MMFF94s harmonic (finite-difference Hessian), intensities from fixed MMFF charges",
            "note": f"Frequencies scaled by {SCALE}. Typical error ±30–80 cm⁻¹; intensities qualitative."}


def _is_linear(xyz: np.ndarray) -> bool:
    if len(xyz) < 3:
        return True
    c = xyz - xyz.mean(axis=0)
    s = np.linalg.svd(c, compute_uv=False)
    return s[1] < 1e-3 * s[0]


# =============================================================================== GFN2-xTB (tblite)
BOHR = 0.529177210903
HARTREE_J = 4.3597447222071e-18
GFN2_SCALE = 0.97          # harmonic → fundamental scaling recommended for GFN2-xTB stretches


def _xtb_available() -> bool:
    try:
        import tblite.interface  # noqa: F401
        return True
    except Exception:
        return False


def normal_modes_gfn2(mol: Chem.Mol, max_atoms=45, seed=0xF00D, step=0.005) -> dict:
    """
    Harmonic analysis at the GFN2-xTB level (Bannwarth, Ehlert, Grimme 2019) via tblite.
      1. ETKDG + MMFF94s pre-optimisation, then L-BFGS optimisation on the GFN2 surface
      2. Hessian from central differences of the analytic GFN2 gradient
      3. IR intensities from the finite-difference dipole derivative  (∂μ/∂Q)²
    This is a genuine electronic-structure calculation (tight-binding DFT, D4 dispersion,
    self-consistent charges): typical harmonic-frequency MAE ≈ 25–35 cm⁻¹ vs experiment
    after scaling, intensities semi-quantitative.  Cost ≈ 6N gradient evaluations.
    """
    if mol.GetNumAtoms() > max_atoms:
        return {"error": f"GFN2-xTB normal modes limited to {max_atoms} heavy atoms (use MMFF level)"}
    if not _xtb_available():
        return {"error": "tblite (GFN2-xTB) not installed"}
    try:
        return _normal_modes_gfn2_inner(mol, max_atoms, seed, step)
    except Exception as e:                      # SCF non-convergence etc. → caller falls back to MMFF
        return {"error": f"GFN2-xTB failed ({type(e).__name__}: {e})"}


def _normal_modes_gfn2_inner(mol, max_atoms, seed, step):
    from tblite.interface import Calculator
    from scipy.optimize import minimize
    mh = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed
    if AllChem.EmbedMolecule(mh, ps) < 0:
        return {"error": "3D embedding failed"}
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=2000)
    except Exception:
        pass
    Z = np.array([a.GetAtomicNum() for a in mh.GetAtoms()])
    charge = Chem.GetFormalCharge(mh)
    x0 = np.array(mh.GetConformer().GetPositions()) / BOHR
    n = len(Z)
    calc = Calculator("GFN2-xTB", Z, x0, charge=charge)
    calc.set("verbosity", 0)
    calc.set("max-iter", 300)

    def eg(flat):
        calc.update(flat.reshape(n, 3))
        r = calc.singlepoint()
        return float(r.get("energy")), np.array(r.get("gradient")).ravel()

    res = minimize(eg, x0.ravel(), jac=True, method="L-BFGS-B", options={"maxiter": 400, "gtol": 1e-5})
    xopt = res.x.copy()
    e0, g0 = eg(xopt)
    if np.max(np.abs(g0)) > 2e-3:
        note_conv = f" (residual gradient {np.max(np.abs(g0)):.1e} Eh/a0, not fully converged)"
    else:
        note_conv = ""
    # Hessian and dipole derivatives
    H = np.zeros((3 * n, 3 * n)); dmu = np.zeros((3 * n, 3))
    for k in range(3 * n):
        xp = xopt.copy(); xp[k] += step
        calc.update(xp.reshape(n, 3)); rp = calc.singlepoint(); gp = np.array(rp.get("gradient")).ravel(); mp = np.array(rp.get("dipole"))
        xm = xopt.copy(); xm[k] -= step
        calc.update(xm.reshape(n, 3)); rm = calc.singlepoint(); gm = np.array(rm.get("gradient")).ravel(); mm = np.array(rm.get("dipole"))
        H[k] = (gp - gm) / (2 * step)
        dmu[k] = (mp - mm) / (2 * step)
    H = 0.5 * (H + H.T)                                # Eh / a0²
    masses = np.array([a.GetMass() for a in mh.GetAtoms()]) * AMU
    m3 = np.repeat(masses, 3)
    H_si = H * HARTREE_J / (BOHR * 1e-10) ** 2         # J / m²
    Hw = H_si / np.sqrt(np.outer(m3, m3))
    w2, V = np.linalg.eigh(Hw)
    freqs = np.sign(w2) * np.sqrt(np.abs(w2)) / (2 * np.pi * C_CM)
    xyz = xopt.reshape(n, 3) * BOHR
    order = np.argsort(np.abs(freqs))
    ntr = 5 if _is_linear(xyz) else 6
    modes = []
    for idx in order[ntr:]:
        f = float(freqs[idx])
        if f < 0:
            continue
        lvec = V[:, idx] / np.sqrt(m3)                 # Cartesian displacement per unit normal coordinate
        cart = lvec.reshape(n, 3) / np.linalg.norm(lvec)
        d = dmu.T @ lvec                                # ∂μ/∂Q (a.u.)
        inten = float(np.dot(d, d))
        modes.append({"freq": round(f * GFN2_SCALE, 1), "freq_unscaled": round(f, 1), "intensity": inten,
                      "label": _mode_label(mh, cart, xyz), "vec": np.round(cart, 4).tolist()})
    if modes:
        mx = max(m["intensity"] for m in modes) or 1.0
        for m in modes:
            m["intensity"] = round(m["intensity"] / mx, 3)
    modes.sort(key=lambda m: -m["freq"])
    return {"symbols": [a.GetSymbol() for a in mh.GetAtoms()], "xyz": np.round(xyz, 4).tolist(),
            "bonds": [[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in mh.GetBonds()],
            "modes": modes, "scale": GFN2_SCALE, "level": "GFN2-xTB", "energy_hartree": round(e0, 6),
            "method": "GFN2-xTB harmonic analysis (tblite): L-BFGS optimisation, finite-difference Hessian of the analytic gradient, IR intensities from dipole derivatives",
            "note": f"Frequencies scaled by {GFN2_SCALE}. Typical error ±25–35 cm⁻¹ vs experiment; intensities semi-quantitative{note_conv}."}


def normal_modes_best(mol: Chem.Mol) -> dict:
    """GFN2-xTB when available and affordable, otherwise MMFF94s; the result says which."""
    r = normal_modes_gfn2(mol)
    if "modes" in r:
        return r
    m = normal_modes(mol)
    m["level"] = "MMFF94s"
    m["fallback_reason"] = r.get("error")
    return m
