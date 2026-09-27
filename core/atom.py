"""
Single-atom visualisation data: electron shells (from the ground-state configuration) and
true hydrogenic atomic-orbital isosurfaces (analytic H-like wavefunctions on a grid → Gaussian
cube), scaled by an effective nuclear charge so the size trend across the table is right.

The orbital shapes (nodal structure, angular lobes) are exact for one-electron wavefunctions;
we render the highest-occupied subshell set. This is a teaching/である visual, not a
correlated many-electron density; stated in the methods registry.
"""
from __future__ import annotations
import re, math
import numpy as np

_SHELL_RE = re.compile(r"(\d)([spdf])(\d+)")
_NOBLE = {"[He]": "1s2", "[Ne]": "1s2 2s2 2p6", "[Ar]": "1s2 2s2 2p6 3s2 3p6",
          "[Kr]": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",
          "[Xe]": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5p6",
          "[Rn]": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5p6 6s2 4f14 5d10 6p6"}


def expand_config(cfg: str) -> list[tuple[int, str, int]]:
    for k, v in _NOBLE.items():
        cfg = cfg.replace(k, v + " ")
    return [(int(n), l, int(e)) for n, l, e in _SHELL_RE.findall(cfg)]


def shells(cfg: str) -> list[dict]:
    """Electrons per principal shell n (for a Bohr-style ring diagram)."""
    agg: dict[int, int] = {}
    for n, l, e in expand_config(cfg):
        agg[n] = agg.get(n, 0) + e
    return [{"n": n, "electrons": agg[n], "capacity": 2 * n * n} for n in sorted(agg)]


def valence_subshell(cfg: str) -> tuple[int, str, int]:
    exp = expand_config(cfg)
    return exp[-1] if exp else (1, "s", 0)


# ---- hydrogenic wavefunctions (real form) -------------------------------------------------
def _radial(n, l, r, Z):
    from scipy.special import genlaguerre, factorial
    a = 1.0
    rho = 2 * Z * r / (n * a)
    norm = math.sqrt((2 * Z / (n * a)) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-rho / 2) * rho ** l * genlaguerre(n - l - 1, 2 * l + 1)(rho)


def _angular(l, m, x, y, z, r):
    r = np.where(r == 0, 1e-9, r)
    if l == 0:
        return np.full_like(x, 0.28209479)                      # s
    if l == 1:
        return {("p", "z"): 0.4886025 * z / r, ("p", "x"): 0.4886025 * x / r, ("p", "y"): 0.4886025 * y / r}[("p", m)]
    if l == 2:
        d = {"z2": 0.3153916 * (3 * z * z - r * r) / (r * r), "xz": 1.0925484 * x * z / (r * r),
             "yz": 1.0925484 * y * z / (r * r), "xy": 1.0925484 * x * y / (r * r),
             "x2y2": 0.5462742 * (x * x - y * y) / (r * r)}
        return d[m]
    return 0.746353 * z * (5 * z * z - 3 * r * r) / (r ** 3)     # one representative f (fz3)


_ORBS = {"s": [("s", None)], "p": [("p", "z"), ("p", "x"), ("p", "y")],
         "d": [("d", "xy"), ("d", "x2y2"), ("d", "xz"), ("d", "yz"), ("d", "z2")], "f": [("f", "z3")]}
_LMAP = {"s": 0, "p": 1, "d": 2, "f": 3}


def orbital_cube(Z_eff: float, n: int, l_sym: str, which: int = 0, npts: int = 56) -> dict:
    l = _LMAP[l_sym]
    orb, m = _ORBS[l_sym][which % len(_ORBS[l_sym])]
    extent = 2.6 * n * n / max(Z_eff, 1.0) + 4.0    # Bohr; matches the orbital's spatial size
    ax = np.linspace(-extent, extent, npts)
    X, Y, Zc = np.meshgrid(ax, ax, ax, indexing="ij")
    R = np.sqrt(X * X + Y * Y + Zc * Zc)
    psi = _radial(n, l, R, Z_eff) * _angular(l, m, X, Y, Zc, R)
    step = (2 * extent) / (npts - 1)
    header = ["atomic orbital", f"{n}{l_sym} psi", f"    1    {-extent:.5f}  {-extent:.5f}  {-extent:.5f}",
              f"  {npts}   {step:.5f}  0.0  0.0", f"  {npts}   0.0  {step:.5f}  0.0", f"  {npts}   0.0  0.0  {step:.5f}",
              f"    1    1.0    0.0    0.0    0.0"]
    vals = psi.flatten()
    lines = [" ".join(f"{v: .5e}" for v in vals[i:i + 6]) for i in range(0, len(vals), 6)]
    iso = 0.30 * float(np.abs(psi).max())    # fraction-of-max isosurface → clean, consistent lobes
    return {"cube": "\n".join(header + lines), "iso": round(iso, 5), "orbital": f"{n}{l_sym}", "which": which, "n_lobes": len(_ORBS[l_sym])}


def atom_view(el: dict) -> dict:
    cfg = el["config"]
    sh = shells(cfg)
    n_v, l_v, e_v = valence_subshell(cfg)
    Z = el["z"]
    # Slater-like effective charge for the valence shell (rough, for orbital size only)
    inner = sum(s["electrons"] for s in sh if s["n"] < n_v)
    Z_eff = max(1.0, Z - 0.85 * (sum(s["electrons"] for s in sh if s["n"] == n_v) - 1) - inner)
    return {"z": Z, "symbol": el["symbol"], "name": el["name"], "config": cfg, "shells": sh,
            "valence": {"n": n_v, "l": l_v, "electrons": e_v}, "z_eff": round(Z_eff, 2),
            "note": "Shells from the ground-state configuration; 3-D shapes are exact hydrogenic (one-electron) orbitals scaled by an effective nuclear charge: correct nodal structure and symmetry, illustrative size."}
