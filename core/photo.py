"""
Photoluminescence estimate + Franck–Condon vibronic structure.

Model (deliberately simple and stated on the page):
  * S0→S1 0-0 energy from the empirical UV–Vis λmax (Woodward–Fieser / Scott / aromatic
    bases in spectra.uv_vis), shifted from the band maximum to the 0-0 origin.
  * One effective totally-symmetric mode ω (cm⁻¹), taken from the MMFF94s normal modes
    (highest-intensity skeletal C=C / C⋯C stretch between 1250 and 1700 cm⁻¹; default 1400).
  * Displaced-harmonic-oscillator Franck–Condon factors with a Huang–Rhys factor S
    estimated from rigidity (fused / planar conjugation → small S, flexible or twisted
    chromophores → larger S):  I(0→n) = e^{-S} S^n / n!
  * Mirror-image emission, Stokes shift ≈ 2 S ħω + solvent reorganisation term.
  * Emissivity class from known quenching motifs (heavy atoms, nitro, n→π* carbonyl,
    free rotors, aliphatic amines / azo, lack of chromophore).
Accuracy: λ ±30–50 nm; vibronic spacing ±100 cm⁻¹; S qualitative.  Good enough to
say "blue-emitting, structured band, weak because of the nitro group", not to
replace a fluorimeter.
"""
from __future__ import annotations
import math
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from . import spectra as X
from . import solvents as SOLV


def _S(s):
    return Chem.MolFromSmarts(s)


def photoluminescence(mol: Chem.Mol, modes: dict | None = None, solvent: str = "chloroform") -> dict:
    """
    Solvent enters twice: through the UV–Vis solvatochromic shift of the 0-0 energy and
    through a Lippert–Mataga reorganisation term  λ_solv ≈ 2 Δμ² Δf / (h c a³); here
    parametrised as  λ_solv = (400 + 4200·ct) · Δf/0.30 cm⁻¹  with ct ∈ [0,1] the estimated
    charge-transfer character (donor–acceptor substitution on the chromophore).
    """
    key, sv = SOLV.get(solvent)
    uv = X.uv_vis(mol, solvent=key)
    n_arom = rdMolDescriptors.CalcNumAromaticRings(mol)
    ndb, system = X._longest_conjugated_path(mol)
    reasons: list[str] = []
    real = [p for p in uv["peaks"] if p["lambda_max"]]
    if not real or (n_arom == 0 and ndb < 3):
        return {"emissive": False, "class": "non-emissive", "reasons": ["No extended chromophore: the lowest excited state lies above ~250 nm / is σ or n-σ* in character; no useful photoluminescence expected."],
                "uv": uv, "note": "Model only applies to π-conjugated chromophores."}

    # --- lowest π→π* band (largest λ among π→π* entries, else largest λ overall)
    pipi = [p for p in real if "n→π*" not in p["assignment"] and "n→σ*" not in p["assignment"]]
    band = max(pipi or real, key=lambda p: p["lambda_max"])
    lam_max = band["lambda_max"]
    npi = [p for p in real if "n→π*" in p["assignment"]]
    n_pi_lowest = bool(npi) and max(p["lambda_max"] for p in npi) > lam_max + 15

    # --- Huang–Rhys factor from rigidity
    ri = mol.GetRingInfo()
    fused = sum(1 for r in ri.AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r) and sum(1 for i in r if ri.NumAtomRings(i) > 1) >= 2)
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    conj_rot = 0
    for b in mol.GetBonds():
        if b.GetBondTypeAsDouble() == 1 and b.GetIsConjugated() and not b.IsInRing() and b.GetBeginAtom().GetDegree() > 1 and b.GetEndAtom().GetDegree() > 1:
            conj_rot += 1          # biaryl / styryl-type torsions inside the chromophore
    S = 0.9
    if fused >= 2:
        S = 0.55; reasons.append("rigid fused aromatic core → small geometry change on excitation (small Huang–Rhys factor, well-resolved vibronic structure)")
    elif fused == 1:
        S = 0.7
    if conj_rot:
        S += 0.45 * min(conj_rot, 3); reasons.append(f"{conj_rot} rotatable bond(s) inside the chromophore → larger excited-state relaxation (broader, more Stokes-shifted emission)")
    if any(a.GetIsAromatic() and a.GetAtomicNum() in (7, 8, 16) for a in mol.GetAtoms()):
        S += 0.1
    S = round(min(S, 2.8), 2)

    # --- effective vibrational mode from MMFF normal modes
    omega = 1400.0; omega_src = "default aromatic ring stretch"
    if modes and modes.get("modes"):
        cands = [m for m in modes["modes"] if 1250 <= m["freq"] <= 1700 and ("C⋯C" in m["label"] or "C=C" in m["label"] or "C–C" in m["label"])]
        if cands:
            best = max(cands, key=lambda m: m["intensity"])
            omega = round(best["freq"], 0); omega_src = f"MMFF94s mode {best['freq']:.0f} cm⁻¹ ({best['label']})"

    # --- charge-transfer character: donor and acceptor both on the chromophore
    donor = mol.HasSubstructMatch(_S("c[NX3;!$(NC=O)]")) or mol.HasSubstructMatch(_S("c[OX2H,OX2C]"))
    acceptor = mol.HasSubstructMatch(_S("c[CX3]=O")) or mol.HasSubstructMatch(_S("cC#N")) or mol.HasSubstructMatch(_S("c[NX3+](=O)[O-]")) or mol.HasSubstructMatch(_S("c[SX4](=O)(=O)"))
    ct = 1.0 if (donor and acceptor) else 0.35 if (donor or acceptor) else 0.0
    df = max(0.0, SOLV.delta_f(sv))
    solvent_reorg = (400.0 + 4200.0 * ct) * (df / 0.30)
    if ct:
        reasons.append(f"{'donor–acceptor' if ct == 1 else 'polar'} substitution → charge-transfer character; emission red-shifts with solvent polarity (Δf = {df:.2f} for {key})")
    # --- 0-0 origin: the band maximum sits ~S·ω above the origin for the absorption
    nu_max = 1e7 / lam_max
    nu00 = nu_max - S * omega
    lam00 = 1e7 / nu00
    stokes = 2 * S * omega + solvent_reorg
    nu_em_max = nu00 - S * omega - solvent_reorg
    lam_em = 1e7 / nu_em_max

    # --- vibronic progressions (Poisson)
    nmax = 6
    fc = [math.exp(-S) * S ** n / math.factorial(n) for n in range(nmax + 1)]
    absorb = [{"n": n, "nu": round(nu00 + n * omega), "nm": round(1e7 / (nu00 + n * omega), 1), "fc": round(fc[n], 4), "label": f"0→{n}′"} for n in range(nmax + 1)]
    emit = [{"n": n, "nu": round(nu00 - solvent_reorg - n * omega), "nm": round(1e7 / (nu00 - solvent_reorg - n * omega), 1), "fc": round(fc[n], 4), "label": f"0′→{n}"} for n in range(nmax + 1)]

    # --- simulated spectra on a nm axis (Gaussian in energy space; broader in polar / protic solvents)
    fwhm = 450 + 500 * df / 0.30 + 250 * sv["hbd"] + 400 * ct
    sigma = fwhm / 2.355
    lo, hi = max(200, lam00 - 130), min(900, lam_em + 220)
    xs = [lo + (hi - lo) * i / 599 for i in range(600)]
    ya, ye = [], []
    for x in xs:
        nu = 1e7 / x
        ya.append(sum(p["fc"] * math.exp(-((nu - p["nu"]) ** 2) / (2 * sigma ** 2)) for p in absorb))
        ye.append(sum(p["fc"] * math.exp(-((nu - p["nu"]) ** 2) / (2 * sigma ** 2)) for p in emit))
    ma, me = max(ya) or 1, max(ye) or 1
    ya = [round(v / ma, 4) for v in ya]; ye = [round(v / me, 4) for v in ye]

    # --- emissivity assessment
    cls, score = "moderate", 0
    if n_pi_lowest or (mol.HasSubstructMatch(_S("[CX3]=O")) and n_arom <= 1 and fused == 0):
        score -= 2; reasons.append("lowest singlet has n→π* character (carbonyl) → fast intersystem crossing, weak fluorescence")
    heavy = [a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() in (35, 53)]
    if heavy:
        score -= 3; reasons.append(f"heavy atom(s) {','.join(sorted(set(heavy)))} → spin–orbit coupling promotes intersystem crossing (fluorescence quenched, phosphorescence possible in rigid media)")
    if mol.HasSubstructMatch(_S("[NX3+](=O)[O-]")):
        score -= 3; reasons.append("nitro group → efficient non-radiative decay; typically non-fluorescent")
    if mol.HasSubstructMatch(_S("N=N")):
        score -= 3; reasons.append("azo group → rapid E/Z photoisomerisation quenches emission")
    if fused >= 2:
        score += 2; reasons.append("rigid polycyclic aromatic → high radiative rate")
    if mol.HasSubstructMatch(_S("c[NX3;H0;!$(NC=O)]")) or mol.HasSubstructMatch(_S("c[NX3;H1,H2;!$(NC=O)]")) or mol.HasSubstructMatch(_S("c[OX2H]")):
        score += 1; reasons.append("electron-donating amino / hydroxy substituent on the chromophore → charge-transfer character, solvatochromic emission")
    if conj_rot >= 2 and fused == 0:
        score -= 1; reasons.append("free rotors allow non-radiative twisting in solution (emission often recovers in viscous media / aggregates)")
    if n_arom == 1 and fused == 0 and score >= 0:
        score -= 1; reasons.append("single benzene-type chromophore → emission in the UV, modest quantum yield")
    cls = "bright" if score >= 2 else "moderate" if score >= 0 else "weak" if score >= -2 else "non-emissive"
    emissive = cls != "non-emissive"

    # --- phosphorescence (T1 ≈ 0.72 × S1 for aromatics; heavy atoms shift the balance)
    nu_T1 = nu00 * 0.72
    lam_phos = 1e7 / (nu_T1 - S * omega)

    return {
        "emissive": emissive, "class": cls, "reasons": reasons,
        "S": S, "omega": omega, "omega_source": omega_src,
        "lambda_abs_max": round(lam_max, 1), "lambda_00": round(lam00, 1), "nu_00": round(nu00), "lambda_em_max": round(lam_em, 1),
        "stokes_shift_cm": round(stokes), "stokes_shift_nm": round(lam_em - lam_max, 1), "solvent": key, "delta_f": round(df, 3), "ct_character": ct, "solvent_reorg_cm": round(solvent_reorg), "fwhm_cm": round(fwhm),
        "lambda_phos": round(lam_phos, 1), "colour": _colour(lam_em), "colour_abs": _complement(lam_max),
        "absorption": absorb, "emission": emit, "x": [round(v, 1) for v in xs], "y_abs": ya, "y_em": ye,
        "note": f"Displaced-harmonic-oscillator Franck–Condon model on empirical S₀→S₁ energies; λ ±30–50 nm, S qualitative. Fluorescence only. Solvent {key}: Δf = {df:.2f}, reorganisation {solvent_reorg:.0f} cm⁻¹, band FWHM {fwhm:.0f} cm⁻¹.",
    }


def _colour(nm: float) -> str:
    if nm < 400: return "ultraviolet (invisible)"
    if nm < 450: return "violet-blue"
    if nm < 495: return "blue"
    if nm < 570: return "green"
    if nm < 590: return "yellow"
    if nm < 620: return "orange"
    if nm < 750: return "red"
    return "near-infrared (invisible)"


def _complement(nm: float) -> str:
    return "colourless" if nm < 380 else "coloured (absorbs in the visible)"
