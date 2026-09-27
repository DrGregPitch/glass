"""
HOSE-code NMR shift prediction trained on nmrshiftdb2 (experimental, atom-assigned spectra).

A HOSE code (Bremser 1978) describes an atom by concentric spheres of its neighbours.  For a
query atom we look up the deepest sphere (4 → 1) for which the database has examples and
return the mean shift, its spread and the number of examples.  This is a nearest-environment
*experimental* prediction: no theory, but directly anchored to measured spectra.  Accuracy is
measured on a held-out split at build time and stored in the model (see `validation`).

Model file: data/hose_model.pkl.gz  (built by `python -m core.hose build`).
"""
from __future__ import annotations
import gzip, pickle, os, sys, math, random, re
from collections import defaultdict
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
ROOT = getattr(sys, "_MEIPASS", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MODEL_PATH = os.path.join(ROOT, "data", "hose_model.pkl.gz")
MAX_DEPTH = 4
_BOND = {1.0: "-", 1.5: ":", 2.0: "=", 3.0: "#"}


def atom_token(a: Chem.Atom) -> str:
    sym = a.GetSymbol()
    if a.GetIsAromatic():
        sym = sym.lower()
    tok = f"{sym}{a.GetTotalNumHs()}"
    if a.GetFormalCharge():
        tok += f"{a.GetFormalCharge():+d}"
    if a.IsInRing():
        tok += "r"
    return tok


def hose(mol: Chem.Mol, idx: int, depth: int) -> str:
    """Canonical sphere code of atom idx up to `depth` bonds (heavy atoms only, H as counts)."""
    def rec(i, parent, d):
        a = mol.GetAtomWithIdx(i)
        t = atom_token(a)
        if d == 0:
            return t
        parts = []
        for b in a.GetBonds():
            j = b.GetOtherAtomIdx(i)
            if j == parent or mol.GetAtomWithIdx(j).GetAtomicNum() == 1:
                continue
            parts.append(_BOND.get(b.GetBondTypeAsDouble(), "-") + rec(j, i, d - 1))
        return t + "(" + ",".join(sorted(parts)) + ")" if parts else t + "()"
    return rec(idx, -1, depth)


# ----------------------------------------------------------------------------- model
class HoseModel:
    def __init__(self, tables: dict, validation: dict):
        self.tables = tables            # nucleus -> depth -> code -> (sum, sumsq, n)
        self.validation = validation

    def predict(self, mol: Chem.Mol, idx: int, nucleus="13C"):
        """Return (shift, sigma, n, depth) or None."""
        t = self.tables.get(nucleus)
        if not t:
            return None
        for d in range(MAX_DEPTH, 0, -1):
            code = hose(mol, idx, d)
            e = t[d].get(code)
            if e and e[2] >= (1 if d >= 3 else 3):
                s, ss, n = e
                mean = s / n
                var = max(0.0, ss / n - mean * mean)
                return round(mean, 2), round(math.sqrt(var), 2), n, d
        return None


_model: HoseModel | None = None


def load() -> HoseModel | None:
    global _model
    if _model is None and os.path.exists(MODEL_PATH):
        with gzip.open(MODEL_PATH, "rb") as f:
            d = pickle.load(f)                       # plain dict: portable across module names
        _model = HoseModel(d["tables"], d["validation"])
    return _model


# ----------------------------------------------------------------------------- build
_SPEC = re.compile(r"^Spectrum (13C|1H) \d+$")


def _records(path):
    """Yield (mol_with_H, {nucleus: [(shift, atom_idx)]}) from the nmrshiftdb2 SD file."""
    suppl = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
    for m in suppl:
        if m is None:
            continue
        spectra = defaultdict(list)
        for k in m.GetPropNames():
            mm = _SPEC.match(k)
            if not mm:
                continue
            nuc = mm.group(1)
            for sig in m.GetProp(k).split("|"):
                bits = sig.split(";")
                if len(bits) < 3:
                    continue
                try:
                    shift = float(bits[0]); ai = int(bits[2])
                except ValueError:
                    continue
                spectra[nuc].append((shift, ai))
        if spectra:
            yield m, spectra


def build(path: str, out: str = MODEL_PATH, holdout=0.05, seed=7):
    random.seed(seed)
    tables = {"13C": {d: defaultdict(lambda: [0.0, 0.0, 0]) for d in range(1, MAX_DEPTH + 1)},
              "1H": {d: defaultdict(lambda: [0.0, 0.0, 0]) for d in range(1, MAX_DEPTH + 1)}}
    test = {"13C": [], "1H": []}
    nmol = 0
    for m, spectra in _records(path):
        nmol += 1
        is_test = random.random() < holdout
        for nuc, sigs in spectra.items():
            seen = set()
            for shift, ai in sigs:
                if ai >= m.GetNumAtoms():
                    continue
                a = m.GetAtomWithIdx(ai)
                if nuc == "13C":
                    if a.GetAtomicNum() != 6:
                        continue
                    target = ai
                else:
                    if a.GetAtomicNum() != 1 or a.GetDegree() != 1:
                        continue
                    target = a.GetNeighbors()[0].GetIdx()         # key protons by their parent heavy atom
                    if m.GetAtomWithIdx(target).GetAtomicNum() != 6:
                        continue                                   # only C–H (exchangeable H excluded)
                if (target, shift) in seen:
                    continue
                seen.add((target, shift))
                codes = {d: hose(m, target, d) for d in range(1, MAX_DEPTH + 1)}
                if is_test:
                    test[nuc].append((codes, shift))
                else:
                    for d, c in codes.items():
                        e = tables[nuc][d][c]; e[0] += shift; e[1] += shift * shift; e[2] += 1
        if nmol % 5000 == 0:
            print(f"  {nmol} molecules", file=sys.stderr)
    # freeze to plain dicts/tuples
    frozen = {nuc: {d: {c: tuple(v) for c, v in t.items()} for d, t in dd.items()} for nuc, dd in tables.items()}
    model = HoseModel(frozen, {})
    # validation on held-out molecules
    val = {}
    for nuc in ("13C", "1H"):
        errs = defaultdict(list); miss = 0
        for codes, shift in test[nuc]:
            pred = None
            for d in range(MAX_DEPTH, 0, -1):
                e = frozen[nuc][d].get(codes[d])
                if e and e[2] >= (1 if d >= 3 else 3):
                    pred = (e[0] / e[2], d); break
            if pred is None:
                miss += 1; continue
            errs[pred[1]].append(abs(pred[0] - shift))
        allerr = [x for v in errs.values() for x in v]
        val[nuc] = {"n_test": len(test[nuc]), "coverage": round(1 - miss / max(1, len(test[nuc])), 3),
                    "MAE": round(sum(allerr) / max(1, len(allerr)), 3),
                    "by_depth": {d: {"n": len(v), "MAE": round(sum(v) / len(v), 3), "RMSE": round(math.sqrt(sum(x * x for x in v) / len(v)), 3)} for d, v in sorted(errs.items())},
                    "n_train_entries": sum(e[2] for e in frozen[nuc][1].values())}
    model.validation = {"molecules": nmol, "holdout": holdout, **val, "source": "nmrshiftdb2 (https://nmrshiftdb.nmr.uni-koeln.de), nmrshiftdb2 Database Licence "
                                  "(derived from ODbL; section 4.5 requires dependent software to be open source)"}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with gzip.open(out, "wb", compresslevel=6) as f:
        pickle.dump({"tables": model.tables, "validation": model.validation}, f, protocol=4)
    print(model.validation)
    return model


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "build":
        build(sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "data", "nmrshiftdb2withsignals.sd"))
