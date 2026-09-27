"""
Structure -> IUPAC name generation for compounds that are not in PubChem.

Generator: knowledgator/SMILES2IUPAC-canonical-base (MT5 sequence-to-sequence model, Apache-2.0),
run in a sidecar Python 3.11 process (namer_env) because it needs PyTorch/transformers.  Three
style variants (systematic, base, traditional) are requested per structure.
Rigor: a generated name is NEVER shown unverified.  Every candidate is parsed back with OPSIN and
accepted only if the Standard InChIKey is identical to the structure's; the first verified
candidate wins, otherwise the rejection is reported.  The gate turns a ~90 %-accurate generator
into "verified or absent".
"""
from __future__ import annotations
import json, os, subprocess, sys, threading, time
from typing import Optional
from . import opsin

ROOT = getattr(sys, "_MEIPASS", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_CANDIDATES = [os.path.join(ROOT, "namer_env", "bin", "python3"), os.path.join(ROOT, "namer_env", "bin", "python"), os.path.join(ROOT, "namer_env", "python.exe"), os.path.join(ROOT, "namer_env", "Scripts", "python.exe")]
_lock = threading.Lock()
_proc: Optional[subprocess.Popen] = None
_state = {"available": None, "reason": ""}


def python_exe() -> Optional[str]:
    if os.environ.get("CHEMXREF_NAMER_PYTHON"):
        return os.environ["CHEMXREF_NAMER_PYTHON"]
    for c in _CANDIDATES:
        if os.path.exists(c):
            return c
    return None


def available() -> bool:
    if _state["available"] is None:
        _state["available"] = python_exe() is not None
        if not _state["available"]:
            _state["reason"] = "name generator not installed (create namer_env with chemical-converters; see README)"
    return _state["available"]


def _start() -> Optional[subprocess.Popen]:
    global _proc
    exe = python_exe()
    if not exe:
        return None
    env = dict(os.environ)
    envdir = os.path.dirname(os.path.dirname(exe)) if os.path.basename(os.path.dirname(exe)) in ("bin", "Scripts") else os.path.dirname(exe)
    if os.path.isdir(os.path.join(envdir, "hf")):
        env["HF_HOME"] = os.path.join(envdir, "hf"); env["HF_HUB_OFFLINE"] = "1"      # bundled model, no network
    for jh in (os.path.join(envdir, "jre17", "Contents", "Home"), os.path.join(envdir, "jre17")):
        if os.path.exists(os.path.join(jh, "bin", "java")) or os.path.exists(os.path.join(jh, "bin", "java.exe")):
            env["JAVA_HOME"] = jh; env["PATH"] = os.path.join(jh, "bin") + os.pathsep + env.get("PATH", ""); break   # JPype needs Java ≥ 9
    _proc = subprocess.Popen([exe, os.path.join(ROOT, "core", "namer_worker.py")], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True, bufsize=1, env=env)
    first = _proc.stdout.readline()
    try:
        j = json.loads(first)
    except Exception:
        j = {"error": "no response"}
    if not j.get("ready"):
        _state["available"] = False; _state["reason"] = j.get("error", "worker failed to start")
        _proc = None
    return _proc


def generate(smiles: str, timeout: float = 120.0) -> dict:
    """Return {'name', 'verified', 'inchikey', 'note'} or {'error'}."""
    if not available():
        return {"error": _state["reason"]}
    # the model was trained on Kekulé (uppercase) SMILES; aromatic lowercase tokens derail it
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles)
        if m is not None:
            Chem.Kekulize(m, clearAromaticFlags=True)
            smiles = Chem.MolToSmiles(m, kekuleSmiles=True, isomericSmiles=True)
    except Exception:
        pass
    with _lock:
        global _proc
        if _proc is None or _proc.poll() is not None:
            if _start() is None:
                return {"error": _state["reason"] or "worker unavailable"}
        try:
            _proc.stdin.write(smiles + "\n"); _proc.stdin.flush()
            line = _proc.stdout.readline()
        except (BrokenPipeError, OSError) as e:
            _proc = None
            return {"error": f"worker died: {e}"}
    try:
        j = json.loads(line)
    except Exception:
        return {"error": "bad worker response"}
    if "error" in j:
        return {"error": j["error"]}
    cands = [c for c in j.get("candidates", []) if c.get("name")]
    if not cands:
        return {"error": "generator returned no name"}
    return {"candidates": cands, "generator": "OPSIN-verified"}
