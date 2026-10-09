"""
Thin, robust wrapper around OPSIN (Open Parser for Systematic IUPAC Nomenclature).

OPSIN is the de-facto reference open-source IUPAC name parser.  We drive the
bundled CLI jar (shipped with py2opsin) in a persistent JVM process so that
each conversion costs milliseconds instead of a JVM start-up.

Outputs we use:
  * SMILES                     (-osmi)
  * Extended SMILES with $_AV$ (-oextendedsmi)  -> per-atom IUPAC locants
  * StdInChIKey                (-ostdinchikey)   -> identity checks
"""
from __future__ import annotations
import os, subprocess, threading, re
from typing import Optional

import glob, shutil, sys as _sys

_ROOT = getattr(_sys, "_MEIPASS", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _find_jar() -> str:
    if os.environ.get("OPSIN_JAR"):
        return os.environ["OPSIN_JAR"]
    for pat in (os.path.join(_ROOT, "vendor", "opsin*.jar"), os.path.join(_ROOT, "py2opsin", "opsin*.jar")):
        hits = glob.glob(pat)
        if hits:
            return hits[0]
    try:
        import py2opsin
        hits = glob.glob(os.path.join(os.path.dirname(py2opsin.__file__), "opsin*.jar"))
        if hits:
            return hits[0]
    except Exception:
        pass
    return ""


def _find_java() -> str:
    """Bundled JRE (vendor/jre) wins, then JAVA_HOME, then PATH."""
    for cand in (os.path.join(_ROOT, "vendor", "jre", "bin", "java"), os.path.join(_ROOT, "vendor", "jre", "bin", "java.exe"),
                 os.path.join(_ROOT, "vendor", "jre", "Contents", "Home", "bin", "java")):
        if os.path.exists(cand):
            return cand
    jh = os.environ.get("JAVA_HOME")
    if jh and os.path.exists(os.path.join(jh, "bin", "java")):
        return os.path.join(jh, "bin", "java")
    return shutil.which("java") or "java"


JAR = _find_jar()
JAVA = _find_java()


class OpsinError(RuntimeError):
    pass


class _Worker:
    """One persistent `java -jar opsin.jar -o<fmt>` process, line-oriented."""

    def __init__(self, fmt: str, extra: tuple[str, ...] = ()):
        self.fmt, self.extra = fmt, extra
        self.lock = threading.Lock()
        self.proc: Optional[subprocess.Popen] = None

    def _start(self):
        if not os.path.exists(JAR):
            raise OpsinError("OPSIN jar not found; run `pip install py2opsin` and make sure Java 8+ is installed.")
        if shutil.which(JAVA) is None and not os.path.exists(JAVA):
            raise OpsinError("Java runtime not found. Install a JRE (https://adoptium.net) or place one in vendor/jre.")
        self.proc = subprocess.Popen(
            [JAVA, "-Xss4m", "-jar", JAR, "-o" + self.fmt, *self.extra],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, bufsize=1,
        )

    def _read_line(self, timeout: float) -> str:
        # readline() on a hung JVM blocks forever (and holds the resolve lock, wedging the
        # whole app) -- select gives us a watchdog. POSIX-only, which is where we run the worker.
        import select
        r, _, _ = select.select([self.proc.stdout], [], [], timeout)
        if not r:
            raise TimeoutError
        return self.proc.stdout.readline()

    def query(self, name: str, _retried: bool = False) -> str:
        name = " ".join(name.split())
        with self.lock:
            if self.proc is None or self.proc.poll() is not None:
                self._start()
            try:
                self.proc.stdin.write(name + "\n")
                self.proc.stdin.flush()
                out = self._read_line(20.0)
            except TimeoutError:
                try:
                    self.proc.kill()
                except Exception:
                    pass
                self.proc = None
            except (BrokenPipeError, OSError):
                self.proc = None
                raise OpsinError("OPSIN process died; retry.")
            else:
                return out.rstrip("\n")
        if _retried:
            raise OpsinError("OPSIN worker timed out twice; giving up on this query.")
        return self.query(name, _retried=True)


_workers: dict[str, _Worker] = {}


def _w(fmt: str, allow_stereo_ignore=False) -> _Worker:
    key = fmt + ("s" if allow_stereo_ignore else "")
    if key not in _workers:
        _workers[key] = _Worker(fmt, ("-s",) if allow_stereo_ignore else ())
    return _workers[key]


def name_to_smiles(name: str) -> str | None:
    """Return SMILES for a systematic name, or None if OPSIN cannot parse it."""
    if not name.strip():
        return None
    out = _w("smi").query(name)
    return out or None


def name_to_extended_smiles(name: str) -> str | None:
    """CXSMILES with `$_AV:...$` atom values holding the IUPAC locant of each heavy atom."""
    out = _w("extendedsmi").query(name)
    return out or None


def name_to_inchikey(name: str) -> str | None:
    out = _w("stdinchikey").query(name)
    return out or None


def failure_reason(name: str) -> str:
    """Slow path: re-run OPSIN with detailed failure analysis to explain why a name failed."""
    try:
        r = subprocess.run([JAVA, "-jar", JAR, "-f"], input=name + "\n", capture_output=True, text=True, timeout=60)
        msg = (r.stderr or "").strip().splitlines()
        msg = [m for m in msg if "unparsable" in m or "uninterpretable" in m or "not parseable" in m]
        return msg[-1] if msg else "OPSIN could not interpret this as a systematic name."
    except Exception as e:  # pragma: no cover
        return f"OPSIN failure analysis unavailable ({e})"


_AV_RE = re.compile(r"\|\$_AV:([^$]*)\$\|")


def parse_locants(extended_smiles: str) -> tuple[str, list[str]]:
    """Split CXSMILES into (plain smiles, [locant per heavy atom in SMILES order])."""
    m = _AV_RE.search(extended_smiles)
    if not m:
        return extended_smiles.strip(), []
    plain = extended_smiles[: m.start()].strip()
    return plain, m.group(1).split(";")


# ----------------------------------------------------------------------------- CML batch
import xml.etree.ElementTree as ET

_CML_NS = "{http://www.xml-cml.org/schema}"


def names_to_cml(names: list[str]) -> dict[str, ET.Element]:
    """Batch-convert names -> CML <molecule> elements (keyed by name). One JVM start per call."""
    clean = [" ".join(n.split()) for n in names if n.strip()]
    if not clean:
        return {}
    r = subprocess.run([JAVA, "-Xss4m", "-jar", JAR, "-ocml"], input="\n".join(clean) + "\n",
                       capture_output=True, text=True, timeout=120)
    out: dict[str, ET.Element] = {}
    try:
        root = ET.fromstring(r.stdout)
    except ET.ParseError:
        return out
    for mol in root.iter(_CML_NS + "molecule"):
        nm = mol.find(_CML_NS + "name")
        if nm is not None and nm.text:
            out[nm.text] = mol
    return out


# ----------------------------------------------------------------------------- parse-tree sidecar
# vendor/opsin-shim/GlassOpsin.class runs OPSIN's own pipeline and reports the token tree with the
# atom ids each token created (see the Java source).  Optional: when the class is missing or the
# JVM fails, callers fall back to the CLI workers above.
_SHIM_DIR = os.path.join(_ROOT, "vendor", "opsin-shim")
_SHIM_CLASS = "uk.ac.cam.ch.wwmm.opsin.GlassOpsin"


class _ShimWorker:
    """One persistent `java -cp opsin.jar:shim GlassOpsin` process; JSON object per line."""

    def __init__(self):
        self.lock = threading.Lock()
        self.proc: Optional[subprocess.Popen] = None
        self.broken = False

    def available(self) -> bool:
        cls = os.path.join(_SHIM_DIR, *_SHIM_CLASS.split(".")) + ".class"
        return not self.broken and os.path.exists(JAR) and os.path.exists(cls)

    def _start(self):
        self.proc = subprocess.Popen(
            [JAVA, "-Xss4m", "-cp", JAR + os.pathsep + _SHIM_DIR, _SHIM_CLASS],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1,
        )
        import select
        r, _, _ = select.select([self.proc.stdout], [], [], 30.0)
        first = self.proc.stdout.readline() if r else ""
        if '"ready"' not in first:
            try:
                self.proc.kill()
            except Exception:
                pass
            self.proc = None
            self.broken = True
            raise OpsinError("OPSIN parse sidecar did not start")

    def query(self, name: str) -> Optional[dict]:
        import json, select
        name = " ".join(name.split())
        if not name or not self.available():
            return None
        with self.lock:
            try:
                if self.proc is None or self.proc.poll() is not None:
                    self._start()
                self.proc.stdin.write(name + "\n")
                self.proc.stdin.flush()
                r, _, _ = select.select([self.proc.stdout], [], [], 20.0)
                if not r:
                    raise TimeoutError
                line = self.proc.stdout.readline()
            except (TimeoutError, BrokenPipeError, OSError, OpsinError):
                try:
                    if self.proc:
                        self.proc.kill()
                except Exception:
                    pass
                self.proc = None
                return None
        try:
            d = json.loads(line)
        except ValueError:
            return None
        return d if "tokens" in d else None


_shim = _ShimWorker()


def shim_available() -> bool:
    return _shim.available()


def name_to_parse(name: str) -> Optional[dict]:
    """OPSIN parse tree for a name: {name, pre, smiles, cml, tokens[{i,p,el,v,g,atoms}], extras}, or None.
    Radicals are allowed, so substituent names (methyl, phenyl) parse to radical SMILES."""
    return _shim.query(name)


def names_to_parse(names: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for n in names:
        d = name_to_parse(n)
        if d:
            out[" ".join(n.split())] = d
    return out


def cml_to_mol(mol_el: ET.Element):
    """
    Build an RDKit molecule from an OPSIN CML <molecule>.  Returns
    (mol_without_H, locant_per_atom, creation_rank_per_atom); all indexed by RDKit atom idx.
    Locant preference: numeric (1, 2, 4a) > primed/lettered (1', N, O') > Greek (alpha).
    """
    from rdkit import Chem
    rw = Chem.RWMol()
    id2idx: dict[str, int] = {}
    locs: list[str | None] = []
    ranks: list[int] = []
    for a in mol_el.find(_CML_NS + "atomArray").findall(_CML_NS + "atom"):
        el = a.get("elementType")
        atom = Chem.Atom(el)
        if a.get("formalCharge"):
            atom.SetFormalCharge(int(a.get("formalCharge")))
        if a.get("isotopeNumber"):
            atom.SetIsotope(int(a.get("isotopeNumber")))
        if a.get("hydrogenCount") is not None and el != "H":
            atom.SetNumExplicitHs(int(a.get("hydrogenCount"))); atom.SetNoImplicit(True)
        idx = rw.AddAtom(atom)
        id2idx[a.get("id")] = idx
        labels = [l.get("value") for l in a.findall(_CML_NS + "label") if l.get("dictRef", "").endswith("locant")]
        num = [l for l in labels if re.fullmatch(r"\d+[a-z]{0,2}", l)]
        lettered = [l for l in labels if re.fullmatch(r"[A-Za-z0-9]{1,4}'*", l) and l not in num and not re.fullmatch(r"[a-z]{3,}", l)]
        greek = [l for l in labels if re.fullmatch(r"[a-z]{3,}", l)]
        locs.append((num or lettered or greek or [None])[0])
        ranks.append(int(re.sub(r"\D", "", a.get("id")) or 0))
    for b in mol_el.find(_CML_NS + "bondArray").findall(_CML_NS + "bond"):
        a1, a2 = b.get("atomRefs2").split()
        order = {"S": Chem.BondType.SINGLE, "D": Chem.BondType.DOUBLE, "T": Chem.BondType.TRIPLE, "A": Chem.BondType.AROMATIC}[b.get("order", "S")]
        rw.AddBond(id2idx[a1], id2idx[a2], order)
    m = rw.GetMol()
    for at in m.GetAtoms():                     # explicit H atoms carry the count; let RDKit recompute
        if at.GetSymbol() != "H" and not at.GetNoImplicit():
            at.SetNoImplicit(False)
    m.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(m, catchErrors=True)
    # drop explicit hydrogens, keep index bookkeeping
    keep = [i for i, at in enumerate(m.GetAtoms()) if at.GetAtomicNum() != 1]
    mh = Chem.RemoveHs(m, sanitize=False)
    Chem.SanitizeMol(mh, catchErrors=True)
    if mh.GetNumAtoms() != len(keep):
        raise OpsinError("H removal mismatch")
    return mh, [locs[i] for i in keep], [ranks[i] for i in keep]
