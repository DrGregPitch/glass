"""
Glass: local web app for rigorous name <-> SMILES cross-checking, journal-style
structure depiction with IUPAC locant numbering, and property / spectrum prediction.

    python app.py            (then open http://127.0.0.1:8765)
"""
from __future__ import annotations
import os, sys, json, functools, threading, time
from collections import defaultdict, deque
from pathlib import Path
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from rdkit import Chem

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from core import names as N, structure as S, properties as P, spectra as X, vibrations as V, photo as PH, solvents as SOLV, methods as M   # noqa: E402

app = FastAPI(title="Glass", version="1.2", docs_url=None, redoc_url=None, openapi_url=None)
_lock = threading.Lock()


MAX_INPUT_CHARS = 10_000     # longest sane identifier; blocks CPU-exhaustion parses
MAX_ATOMS = 600              # heavy+H atoms; far above anything the predictors accept


def _mol(smiles: str) -> Chem.Mol:
    if len(smiles) > MAX_INPUT_CHARS:
        raise HTTPException(413, f"input too long (> {MAX_INPUT_CHARS} characters)")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise HTTPException(400, "bad SMILES")
    if m.GetNumAtoms() > MAX_ATOMS:
        raise HTTPException(413, f"molecule too large (> {MAX_ATOMS} atoms)")
    return m


@app.get("/api/resolve")
def api_resolve(q: str = Query(..., min_length=1, max_length=10_000), pubchem: bool = True):
    try:
        with _lock:                    # OPSIN CML batch uses one JVM; serialise to keep memory sane
            r = N.resolve(q, want_pubchem=pubchem)
    except ValueError as e:
        raise HTTPException(422, str(e))
    return N.resolution_to_dict(r)


def _reject_fragment(smiles: str) -> None:
    """Substituent inputs carry attachment points (*); spectra and properties are meaningless for them."""
    if "*" in smiles:
        raise HTTPException(400, "substituent fragment (attachment point *): spectra and properties are not computed")


@functools.lru_cache(maxsize=512)
def _structure(smiles: str, three_d: bool):
    m = _mol(smiles)
    out = {"depiction": S.depict_svg(m), "depiction_h": S.depict_svg_h(m)}
    if three_d:
        out["geometry"] = S.geometry3d(m)
    return out


@app.get("/api/structure")
def api_structure(smiles: str, three_d: bool = True):
    return _structure(smiles, three_d)


@functools.lru_cache(maxsize=512)
def _props(smiles: str):
    _reject_fragment(smiles)
    return P.predict(_mol(smiles), smiles)


@app.get("/api/properties")
def api_properties(smiles: str):
    return _props(smiles)


@functools.lru_cache(maxsize=1024)
def _spectra(smiles: str, solvent: str, mhz: float, nmr_level: str = "hose"):
    _reject_fragment(smiles)
    return X.all_spectra(_mol(smiles), solvent, mhz, nmr_level)


@app.get("/api/spectra")
def api_spectra(smiles: str, solvent: str = "chloroform", mhz: float = 400.0, nmr: str = "hose"):
    """nmr: hose (experimental-database environments) | increments | giao (quantum GIAO, scaled; ≤16 heavy atoms)"""
    if nmr not in ("hose", "increments", "giao"):
        raise HTTPException(400, "nmr must be hose, increments or giao")
    try:
        key, _ = SOLV.get(solvent)
    except KeyError as e:
        raise HTTPException(400, str(e))
    return _spectra(smiles, key, float(mhz), nmr)


@functools.lru_cache(maxsize=256)
def _modes(smiles: str, level: str = "auto"):
    _reject_fragment(smiles)
    m = _mol(smiles)
    if level == "mmff":
        r = V.normal_modes(m); r["level"] = "MMFF94s"; return r
    if level == "gfn2":
        return V.normal_modes_gfn2(m)
    return V.normal_modes_best(m)


@app.get("/api/modes")
def api_modes(smiles: str, level: str = "auto"):
    """level: auto (GFN2-xTB when ≤25 heavy atoms, else MMFF94s) | gfn2 (up to 45, slow for large molecules) | mmff"""
    if level not in ("auto", "gfn2", "mmff"):
        raise HTTPException(400, "level must be auto, gfn2 or mmff")
    return _modes(smiles, level)


@app.get("/api/validation")
def api_validation():
    """Measured accuracy of the data-driven models (held-out test sets) and IR benchmarks."""
    out = {}
    from core import hose as H
    hm = H.load()
    if hm and getattr(hm, "validation", None):
        v = hm.validation
        def _nuc(k):
            d = v.get(k, {})
            bd = d.get("by_depth") or {}
            deep = bd.get(4) or bd.get("4") or {}
            return {"MAE": d.get("MAE"), "n_test": d.get("n_test"),
                    "coverage": d.get("coverage"), "units": "ppm",
                    "deep_match": {"MAE": deep.get("MAE"), "n": deep.get("n")}}
        out["nmr_nmrshiftdb2"] = {
            "13C": _nuc("13C"),
            "1H": _nuc("1H"),
            "molecules": v.get("molecules"),
            "holdout_fraction": v.get("holdout"),
            "validation": "held-out experimental spectra",
        }
    else:
        out["nmr_nmrshiftdb2"] = {"error": "model not built"}
    bench = ROOT / "data" / "benchmarks.json"
    if bench.exists():
        out["ir"] = json.loads(bench.read_text(encoding="utf-8"))
    return out


@functools.lru_cache(maxsize=512)
def _photo(smiles: str, solvent: str, level: str = "auto"):
    _reject_fragment(smiles)
    m = _mol(smiles)
    modes = _modes(smiles, level) if m.GetNumAtoms() <= 60 else None
    return PH.photoluminescence(m, modes if modes and not modes.get("error") else None, solvent)


@app.get("/api/photo")
def api_photo(smiles: str, solvent: str = "chloroform", level: str = "auto"):
    try:
        key, _ = SOLV.get(solvent)
    except KeyError as e:
        raise HTTPException(400, str(e))
    return _photo(smiles, key, level)


@functools.lru_cache(maxsize=1)
def _elements():
    return json.loads((ROOT / "data" / "periodic_table.json").read_text(encoding="utf-8"))


@functools.lru_cache(maxsize=1)
def _el_by_z():
    return {e["z"]: e for e in _elements()["elements"]}


@app.get("/api/atom")
def api_atom(z: int):
    from core import atom as AT
    e = _el_by_z().get(z)
    if not e:
        raise HTTPException(404, "element not found")
    return AT.atom_view(e)


@app.get("/api/atom_orbital")
def api_atom_orbital(z: int, n: int, l: str, which: int = 0):
    from fastapi.responses import PlainTextResponse
    from core import atom as AT
    e = _el_by_z().get(z)
    if not e or l not in ("s", "p", "d", "f"):
        raise HTTPException(400, "bad request")
    if not (1 <= n <= 7 and {"s": 0, "p": 1, "d": 2, "f": 3}[l] < n):
        raise HTTPException(400, f"orbital {n}{l} does not exist")
    v = AT.atom_view(e)
    return AT.orbital_cube(v["z_eff"], n, l, which)


@app.get("/api/image")
def api_image(smiles: str, w: int = 500, h: int = 380):
    """PNG depiction of a structure (for link previews / embeds)."""
    from fastapi.responses import Response
    from rdkit.Chem.Draw import rdMolDraw2D
    w = max(100, min(w, 2000)); h = max(100, min(h, 2000))
    m = _mol(smiles)
    S.rdDepictor.SetPreferCoordGen(True); S.rdDepictor.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DCairo(w, h)
    d.drawOptions().useBWAtomPalette(); d.DrawMolecule(m); d.FinishDrawing()
    return Response(d.GetDrawingText(), media_type="image/png")


@app.get("/api/elements")
def api_elements():
    """Periodic-table data (PubChem, public domain): 118 elements with physical properties and layout positions."""
    return _elements()


@functools.lru_cache(maxsize=256)
def _electronic(smiles: str):
    _reject_fragment(smiles)
    from core import electronic as EL
    return EL.electronic_structure(_mol(smiles))


@app.get("/api/electronic")
def api_electronic(smiles: str):
    """GFN2-xTB molecular-orbital levels, gap, dipole and per-atom charges (for surfaces)."""
    return _electronic(smiles)


@functools.lru_cache(maxsize=96)
def _orbital(smiles: str, which: str, offset: int = 0):
    _reject_fragment(smiles)
    from core import orbitals as ORB
    return ORB.orbital_cube(_mol(smiles), which, offset)


@app.get("/api/orbital")
def api_orbital(smiles: str, which: str = "homo", offset: int = 0):
    """Frontier-orbital isosurface (HF/STO-3G): cube data + orbital energy + charge-centroid vector.
    which=homo|lumo with offset 0–5 away from the gap (homo,2 → HOMO−2)."""
    if which not in ("homo", "lumo") or not (0 <= offset <= 5):
        raise HTTPException(400, "which must be homo or lumo with offset 0–5")
    try:
        return _orbital(smiles, which, offset)
    except ValueError as e:
        raise HTTPException(422, str(e))


@app.post("/api/read")
async def api_read(request: Request, fmt: str = "cdxml"):
    """Parse an uploaded CDXML or MOL/SDF text body → canonical SMILES (for drag-and-drop)."""
    from rdkit import Chem
    body = (await request.body()).decode("utf-8", errors="replace")
    if len(body) > 5_000_000:
        raise HTTPException(413, "file too large")
    try:
        if fmt == "cdxml":
            mols = [m for m in Chem.MolsFromCDXML(body) if m is not None and 0 < m.GetNumAtoms() <= MAX_ATOMS]
        else:
            m1 = Chem.MolFromMolBlock(body)
            mols = [m1] if m1 is not None else []
        if not mols:
            raise ValueError("no molecule found in file")
        combined = mols[0]
        for m in mols[1:]:
            combined = Chem.CombineMols(combined, m)
        return {"smiles": Chem.MolToSmiles(combined), "n_fragments": len(mols)}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(422, f"could not parse {fmt}: {e}")


@app.get("/api/cdxml")
def api_cdxml(smiles: str, name: str | None = None, numbers: bool = False):
    """ChemDraw CDXML export of the 2D structure; optional IUPAC locant annotations from `name`."""
    from fastapi.responses import PlainTextResponse
    m = _mol(smiles)
    locs = frags = None
    if numbers and name:
        e = N.NameEntry(name, "adhoc", "user")
        N.verify_name(e, N.std_inchikey(m))
        if e.verified:
            cml = N.opsin.names_to_cml([name])
            el = cml.get(" ".join(name.split()))
            if el is not None:
                try:
                    N.attach_locants(e, m, el)
                    locs, frags = e.locants, e.fragments
                except Exception:
                    pass
    from core import cdxml as CX
    txt = CX.mol_to_cdxml(m, locs, frags, include_numbers=numbers, title=name)
    return PlainTextResponse(txt, media_type="chemical/x-cdxml")


@app.get("/api/methods")
def api_methods():
    """Level of theory, scope and limitations of every simulation in the app (with measured accuracy where available)."""
    from core import hose as H
    methods = json.loads(json.dumps(M.METHODS))
    hm = H.load()
    if hm and hm.validation:
        v = hm.validation
        c, h = v.get("13C", {}), v.get("1H", {})
        methods["nmr"]["error"] = ("Database-level predictions are validated on a held-out set of experimental spectra; "
                                   "each peak reports its own match quality and uncertainty. Atoms without a close database "
                                   "match fall back to additivity increments.")
    return methods


@app.get("/api/experimental")
def api_experimental(cid: int):
    return N.pubchem_experimental(cid)


@app.get("/api/verify")
def api_verify(name: str, inchikey: str):
    """
    Compare ANY identifier (systematic name, SMILES, InChI, InChIKey) with a reference structure.
    status: identical | same_connectivity | different | uninterpretable
    Nothing is judged 'different' unless the input was actually interpreted.
    """
    text = name.strip()
    try:
        with _lock:
            r = N.resolve(text, want_pubchem=False)
    except ValueError as e:
        return {"input": text, "status": "uninterpretable", "kind": None,
                "note": "Not recognised as a systematic name, SMILES, InChI or InChIKey; no comparison was made.", "detail": str(e)}
    if r.inchikey == inchikey:
        status, note = "identical", "Same compound: identical Standard InChIKey (connectivity, stereo, isotopes, charge)."
    elif r.inchikey.split("-")[0] == inchikey.split("-")[0]:
        status, note = "same_connectivity", "Same skeleton but the stereo / isotope / protonation layer differs (InChIKey first block matches, second/third do not)."
    else:
        status, note = "different", "Different compound: InChIKeys do not match."
    return {"input": text, "status": status, "kind": r.input_kind, "smiles": r.smiles_canonical, "inchikey": r.inchikey,
            "formula": r.formula, "note": note}


@app.get("/api/suggest")
def api_suggest(q: str):
    """PubChem name autocomplete for 'did you mean' on failed inputs."""
    j = N._get_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/{N.requests.utils.quote(q)}/json?limit=8")
    try:
        return {"suggestions": j["dictionary_terms"]["compound"]}
    except Exception:
        return {"suggestions": []}


# AGPL-3.0 section 13: network users must be offered the corresponding source.
SOURCE_URL = os.environ.get("GLASS_SOURCE_URL", "https://github.com/DrGregPitch/glass")


@app.get("/api/about")
def api_about():
    """Provenance record: licence, source offer, third-party components and method references."""
    from fastapi.responses import PlainTextResponse
    header = ("Glass — Copyright (C) 2026 Gregory Pitch\n"
              "Free software under the GNU Affero General Public License v3 or later.\n"
              f"Corresponding source: {SOURCE_URL}\n"
              "Full licence: /api/license · data and model terms: data/LICENSE-DATA.md\n"
              "The NMR model is derived from nmrshiftdb2 and carries the nmrshiftdb2 Database\n"
              "Licence (ODbL-derived, attribution and share-alike): /api/model\n\n")
    return PlainTextResponse(header + (ROOT / "THIRD_PARTY.md").read_text(encoding="utf-8"))


@app.get("/api/license")
def api_license():
    """The project licence, verbatim."""
    from fastapi.responses import PlainTextResponse
    return PlainTextResponse((ROOT / "LICENSE").read_text(encoding="utf-8"))


@app.get("/api/model")
def api_model():
    """The NMR model file, offered as the nmrshiftdb2 share-alike terms require."""
    from fastapi.responses import FileResponse
    return FileResponse(ROOT / "data" / "hose_model.pkl.gz", media_type="application/gzip",
                        filename="hose_model.pkl.gz")


@app.get("/api/solvents")
def api_solvents():
    return [{"name": k, "nmr": v["nmr"], "eps": v["eps"], "et30": v["et30"], "delta_f": round(SOLV.delta_f(v), 3)} for k, v in SOLV.SOLVENTS.items()]


_BOT_UA = ("twitterbot", "facebookexternalhit", "linkedinbot", "discordbot",
           "telegrambot", "whatsapp", "embedly", "iframely", "opengraph")


@app.get("/")
def index(request: Request, q: str = ""):
    """Serve the app; for link-preview bots with ?q=, inject Open Graph tags (structure image +
    name). Interactive browsers get the HTML instantly; the resolve happens client-side."""
    from fastapi.responses import HTMLResponse
    import urllib.parse
    html = (ROOT / "static" / "index.html").read_text(encoding="utf-8")
    dash = chr(8212); dot = chr(183)
    ua = (request.headers.get("user-agent") or "").lower()
    if q.strip() and not any(b in ua for b in _BOT_UA):
        q = ""    # human browser: serve immediately with default tags
    og = ('<meta property="og:title" content="Glass ' + dash + ' chemistry console">'
          '<meta property="og:description" content="Names, structures, spectra and electronic structure.">'
          '<meta name="twitter:card" content="summary_large_image">')
    if q.strip():
        try:
            with _lock:
                r = N.resolve(q.strip(), want_pubchem=True)
            title = r.common_name or r.iupac_name or r.formula
            base = os.environ.get("GLASS_PUBLIC_URL", "")
            img = (base + "/api/image?smiles=" + urllib.parse.quote(r.smiles_canonical)) if base else ""
            desc = r.formula
            if r.iupac_name:
                desc = desc + "  " + dot + "  " + r.iupac_name
            if r.cas:
                desc = desc + "  " + dot + "  CAS " + r.cas
            og = ('<meta property="og:title" content="' + _hesc(title + " " + dash + " Glass") + '">'
                  '<meta property="og:description" content="' + _hesc(desc) + '">'
                  + ('<meta property="og:image" content="' + _hesc(img) + '">' if img else "")
                  + '<meta name="twitter:card" content="summary_large_image">')
        except Exception:
            pass
    html = html.replace("</head>", og + "</head>", 1)
    return HTMLResponse(html, headers={"Cache-Control": "no-cache"})


def _hesc(t: str) -> str:
    return (t or "").replace("&", "&amp;").replace(chr(34), "&quot;").replace("<", "&lt;").replace(">", "&gt;")


@app.middleware("http")
async def _host_guard(request: Request, call_next):
    # DNS-rebinding defense: when bound to localhost, only accept localhost Host headers.
    # A hosted deployment sets GLASS_ALLOWED_HOST (or GLASS_PUBLIC_URL) to its domain.
    host = (request.headers.get("host") or "").split(":")[0].lower()
    allowed = {"localhost", "127.0.0.1", "[::1]", "testserver", ""}
    extra = os.environ.get("GLASS_ALLOWED_HOST", "")
    if extra:
        allowed.add(extra.lower())
    pub = os.environ.get("GLASS_PUBLIC_URL", "")
    if pub:
        from urllib.parse import urlparse
        allowed.add((urlparse(pub).hostname or "").lower())
    if os.environ.get("GLASS_HOST", os.environ.get("CHEMXREF_HOST", "127.0.0.1")) != "127.0.0.1":
        return await call_next(request)   # LAN-shared mode: colleagues connect by IP/mDNS name
    if host not in allowed:
        from fastapi.responses import PlainTextResponse
        return PlainTextResponse("forbidden host", status_code=403)
    return await call_next(request)


@app.middleware("http")
async def _no_cache_for_app_assets(request, call_next):
    """Our own JS/HTML must never be served stale; vendored libraries may be cached."""
    resp = await call_next(request)
    if request.url.path.startswith("/static/js/") or request.url.path == "/":
        resp.headers["Cache-Control"] = "no-cache"
    return resp


# --- per-IP rate limiting: throttles a single source flooding the (public,
# unauthenticated) compute endpoints. A distributed flood still needs a CDN/WAF,
# but this stops the common single-origin abuse without adding any user friction. ---
_HEAVY_PATHS = ("/api/spectra", "/api/modes", "/api/orbital", "/api/electronic", "/api/photo")
# NOTE: limiter state is per uvicorn worker, so the effective per-IP allowance is these
# numbers times the worker count (see deploy/setup.sh). Values chosen accordingly.
# Requires GLASS_TRUSTED_PROXY=1 behind a proxy, or every client keys on the proxy IP.
_RL_GENERAL = (90, 30.0)       # ~180 requests / 30 s per IP effective (UI-generous)
_RL_HEAVY = (10, 60.0)         # ~20 heavy computations / 60 s per IP effective
_rl_general: dict[str, deque] = defaultdict(deque)
_rl_heavy: dict[str, deque] = defaultdict(deque)
_rl_lock = threading.Lock()


_TRUST_PROXY = os.environ.get("GLASS_TRUSTED_PROXY", "").strip().lower() in ("1", "true", "yes", "on")


def _client_ip(request: Request) -> str:
    # X-Forwarded-For is client-controlled and only trustworthy behind a proxy that appends the
    # real client IP as the last hop. Trust it only when GLASS_TRUSTED_PROXY says such a proxy is
    # in front; otherwise use the direct peer address, so a forged header cannot mint fresh
    # rate-limit buckets and bypass the throttle on the expensive compute endpoints.
    if _TRUST_PROXY:
        xff = request.headers.get("x-forwarded-for", "")
        if xff:
            return xff.split(",")[-1].strip()
    return request.client.host if request.client else "unknown"


def _rl_allow(bucket: dict, ip: str, limit: int, window: float, now: float) -> bool:
    dq = bucket[ip]
    cutoff = now - window
    while dq and dq[0] < cutoff:
        dq.popleft()
    if len(dq) >= limit:
        return False
    dq.append(now)
    return True


def _rl_prune(now: float) -> None:
    for bucket in (_rl_general, _rl_heavy):
        for ip in [ip for ip, dq in bucket.items() if not dq or dq[-1] < now - 120]:
            del bucket[ip]


@app.middleware("http")
async def _rate_limit(request: Request, call_next):
    path = request.url.path
    if path.startswith("/static") or path.startswith("/favicon"):
        return await call_next(request)
    peer = request.client.host if request.client else ""
    if not request.headers.get("x-forwarded-for") and peer in ("127.0.0.1", "::1", "localhost", "testclient"):
        return await call_next(request)          # loopback and tests: no proxy, no throttling
    ip = _client_ip(request)
    now = time.time()
    with _rl_lock:
        ok = _rl_allow(_rl_general, ip, *_RL_GENERAL, now)
        if ok and path in _HEAVY_PATHS:
            ok = _rl_allow(_rl_heavy, ip, *_RL_HEAVY, now)
        if len(_rl_general) > 20000:
            _rl_prune(now)
    if not ok:
        from fastapi.responses import PlainTextResponse
        return PlainTextResponse("rate limit exceeded \u2014 slow down and retry shortly",
                                 status_code=429, headers={"Retry-After": "10"})
    return await call_next(request)


# --- security response headers (applied to every response, incl. static and the Ketcher iframe) ---
_SECURITY_HEADERS = {
    "X-Content-Type-Options": "nosniff",
    "X-Frame-Options": "SAMEORIGIN",
    "Referrer-Policy": "strict-origin-when-cross-origin",
    "Permissions-Policy": "geolocation=(), microphone=(), camera=(), usb=(), interest-cohort=()",
    # Allows exactly what the app loads: its own assets, inline scripts/styles, the Ketcher
    # WebAssembly editor, 3Dmol from cdnjs, and Google Fonts. object/base/frame-ancestors locked down.
    "Content-Security-Policy": (
        "default-src 'self'; "
        "script-src 'self' 'unsafe-inline' 'unsafe-eval' 'wasm-unsafe-eval' https://cdnjs.cloudflare.com; "
        "style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; "
        "font-src 'self' https://fonts.gstatic.com data:; "
        "img-src 'self' data: blob:; "
        "connect-src 'self'; "
        "worker-src 'self' blob:; "
        "child-src 'self' blob:; "
        "frame-src 'self'; "
        "frame-ancestors 'self'; "
        "base-uri 'self'; "
        "object-src 'none'; "
        "form-action 'self'"
    ),
}


@app.middleware("http")
async def _security_headers(request: Request, call_next):
    resp = await call_next(request)
    for k, v in _SECURITY_HEADERS.items():
        resp.headers.setdefault(k, v)
    if request.headers.get("x-forwarded-proto", "").lower() == "https":   # only meaningful over TLS
        resp.headers.setdefault("Strict-Transport-Security", "max-age=31536000; includeSubDomains")
    return resp


app.mount("/static", StaticFiles(directory=ROOT / "static"), name="static")

if __name__ == "__main__":
    # Direct entry point for running the server; deployments serve `app:app` with any ASGI server.
    import uvicorn
    host = os.environ.get("GLASS_HOST", "127.0.0.1")
    port = int(os.environ.get("GLASS_PORT", "8765"))
    uvicorn.run(app, host=host, port=port, log_level="warning")
