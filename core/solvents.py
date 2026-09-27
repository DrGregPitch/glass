"""
Solvent table used consistently by every simulation.

Each solvent carries the physical constants the individual models need:
  eps   static dielectric constant
  n     refractive index (20 °C, D line)
  et30  Reichardt E_T(30) polarity (kcal/mol)
  hbd   hydrogen-bond-donor ability (protic) – Kamlet–Taft alpha (approx.)
  hba   hydrogen-bond-acceptor ability – Kamlet–Taft beta (approx.)
  nmr   deuterated analogue used for NMR (residual 1H shift, 13C shift, exchangeable-proton behaviour)
Lippert–Mataga orientation polarisability  Δf = (ε−1)/(2ε+1) − (n²−1)/(2n²+1).
"""
from __future__ import annotations

SOLVENTS: dict[str, dict] = {
    # name           eps     n      ET30   alpha  beta   deuterated  residual  c13    OH   phenol NH   amide COOH  exchange
    "chloroform":   dict(eps=4.81,  n=1.446, et30=39.1, hbd=0.20, hba=0.10, nmr="CDCl3",      residual=7.26, c13=77.16, OH=2.0, phenol=5.0, NH=1.5, amide=6.0, COOH=11.5, exchange=False),
    "DMSO":         dict(eps=46.7,  n=1.479, et30=45.1, hbd=0.00, hba=0.76, nmr="DMSO-d6",    residual=2.50, c13=39.52, OH=4.4, phenol=9.4, NH=2.0, amide=7.9, COOH=12.2, exchange=False),
    "water":        dict(eps=80.1,  n=1.333, et30=63.1, hbd=1.17, hba=0.47, nmr="D2O",        residual=4.79, c13=None,  OH=None, phenol=None, NH=None, amide=None, COOH=None, exchange=True),
    "methanol":     dict(eps=32.7,  n=1.329, et30=55.4, hbd=0.98, hba=0.66, nmr="CD3OD",      residual=3.31, c13=49.00, OH=None, phenol=None, NH=None, amide=None, COOH=None, exchange=True),
    "benzene":      dict(eps=2.27,  n=1.501, et30=34.3, hbd=0.00, hba=0.10, nmr="C6D6",       residual=7.16, c13=128.06, OH=1.0, phenol=4.5, NH=1.0, amide=5.5, COOH=11.0, exchange=False),
    "acetone":      dict(eps=20.7,  n=1.359, et30=42.2, hbd=0.08, hba=0.43, nmr="acetone-d6", residual=2.05, c13=29.84, OH=3.5, phenol=8.3, NH=1.8, amide=7.2, COOH=11.0, exchange=False),
    "acetonitrile": dict(eps=37.5,  n=1.344, et30=45.6, hbd=0.19, hba=0.40, nmr="CD3CN",      residual=1.94, c13=1.32,  OH=2.8, phenol=7.0, NH=1.6, amide=6.8, COOH=10.5, exchange=False),
    "hexane":       dict(eps=1.88,  n=1.375, et30=31.0, hbd=0.00, hba=0.00, nmr=None,         residual=None, c13=None,  OH=1.0, phenol=4.5, NH=0.9, amide=5.5, COOH=11.0, exchange=False),
    "cyclohexane":  dict(eps=2.02,  n=1.426, et30=30.9, hbd=0.00, hba=0.00, nmr="C6D12",      residual=1.38, c13=26.43, OH=1.0, phenol=4.5, NH=0.9, amide=5.5, COOH=11.0, exchange=False),
    "toluene":      dict(eps=2.38,  n=1.497, et30=33.9, hbd=0.00, hba=0.11, nmr="toluene-d8", residual=7.09, c13=137.86, OH=1.0, phenol=4.5, NH=1.0, amide=5.5, COOH=11.0, exchange=False),
    "THF":          dict(eps=7.58,  n=1.407, et30=37.4, hbd=0.00, hba=0.55, nmr="THF-d8",     residual=3.58, c13=67.57, OH=3.0, phenol=8.0, NH=1.5, amide=7.0, COOH=11.0, exchange=False),
    "dichloromethane": dict(eps=8.93, n=1.424, et30=40.7, hbd=0.13, hba=0.10, nmr="CD2Cl2",   residual=5.32, c13=53.84, OH=2.0, phenol=5.0, NH=1.5, amide=6.0, COOH=11.5, exchange=False),
    "ethanol":      dict(eps=24.5,  n=1.361, et30=51.9, hbd=0.86, hba=0.75, nmr="ethanol-d6", residual=3.56, c13=56.96, OH=None, phenol=None, NH=None, amide=None, COOH=None, exchange=True),
    "DMF":          dict(eps=36.7,  n=1.431, et30=43.2, hbd=0.00, hba=0.69, nmr="DMF-d7",     residual=8.03, c13=162.62, OH=4.0, phenol=9.0, NH=2.0, amide=7.8, COOH=12.0, exchange=False),
    "gas phase":    dict(eps=1.0,   n=1.0,   et30=27.0, hbd=0.00, hba=0.00, nmr=None,         residual=None, c13=None,  OH=0.5, phenol=4.0, NH=0.5, amide=5.0, COOH=10.0, exchange=False),
}
DEFAULT = "chloroform"

# accept deuterated names too
_ALIAS = {v["nmr"]: k for k, v in SOLVENTS.items() if v["nmr"]}


def get(name: str | None) -> tuple[str, dict]:
    if not name:
        return DEFAULT, SOLVENTS[DEFAULT]
    key = name if name in SOLVENTS else _ALIAS.get(name)
    if key is None:
        raise KeyError(f"unknown solvent {name!r}; choose one of {list(SOLVENTS)}")
    return key, SOLVENTS[key]


def delta_f(s: dict) -> float:
    """Lippert–Mataga orientation polarisability."""
    e, n = s["eps"], s["n"]
    if e <= 1.0:
        return 0.0
    return (e - 1) / (2 * e + 1) - (n * n - 1) / (2 * n * n + 1)


def summary(name: str) -> dict:
    k, s = get(name)
    return {"name": k, "nmr": s["nmr"], "eps": s["eps"], "n": s["n"], "et30": s["et30"], "delta_f": round(delta_f(s), 3), "protic": s["hbd"] > 0.5}
