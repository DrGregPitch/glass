# Glass

A web console for chemical identity, structure, spectra and electronic structure.
Hosted at **[chemistryconsole.glass](https://chemistryconsole.glass)**.

![CI](https://github.com/DrGregPitch/glass/actions/workflows/ci.yml/badge.svg)
&nbsp;·&nbsp; AGPL-3.0-or-later &nbsp;·&nbsp; Python 3.12

Glass cross-checks chemical names against SMILES, draws journal-style structures with IUPAC
locant numbering, and predicts properties and spectra. PubChem is used only to enrich results
with database names and experimental data, and can be switched off per request.

## What it does

| Input (any of) | Output |
|---|---|
| IUPAC / systematic name (any accepted style: preferred, CAS-inverted, trivial-systematic mixes), trivial or trade name, SMILES (any form), InChI, InChIKey, CAS RN | canonical SMILES (isomeric), canonical non-stereo SMILES, Kekulé SMILES, SMILES as entered / from OPSIN / from PubChem, InChI, StdInChIKey, formula, PubChem CID |
| | every name PubChem knows (IUPAC, title, synonyms), each re-parsed by OPSIN and marked `verified` only if the StdInChIKey is identical |
| | ACS-1996-style 2D depiction, per-name IUPAC atom numbering (parent / substituent / heteroatom locants toggleable), hover an atom → the locant is highlighted in the name |
| | 3D geometry (ETKDGv3 + MMFF94s): bond lengths, bond angles, interactive viewer |
| | formula, molar mass, exact mass, melting/boiling point (Joback), density (Girolami), logP, TPSA, etc., plus PubChem experimental values side-by-side |
| | ¹H NMR (solvent & MHz adjustable, first-order multiplets), ¹³C{¹H} NMR — **experimental-database prediction** (nmrshiftdb2 HOSE codes, measured held-out ¹³C MAE 2.75 ppm / 1.18 ppm on deep matches; per-peak basis, σ and n shown) with additivity fallback; IR, UV-Vis; zoomable, peak ↔ table ↔ atom linked |
| | **Solvent** selectable once for all techniques (15 solvents incl. gas phase): NMR exchangeables & reference, IR C=O / X–H shifts and broadening, UV–Vis solvatochromism (E_T(30)), emission Stokes shift (Lippert–Mataga Δf × CT character) |
| | **Editable structure** (Ketcher) — draw changes and everything re-resolves; exports: SVG/PNG/MOL/CDXML structure (numbering annotations optional), PNG/CSV for every spectrum, SVG Franck–Condon diagram, PNG of the 3D and vibration viewers, Markdown/JSON reports |
| | **Photoluminescence** (estimate): S₀→S₁ 0-0 energy, emission maximum and colour, Stokes shift, Huang–Rhys factor, vibronic progressions from a displaced-oscillator Franck–Condon model, emissivity class with reasons; interactive Franck–Condon diagram |
| | **Vibrational normal modes** at the **GFN2-xTB** level (tight-binding via tblite: optimisation, gradient Hessian, dipole-derivative intensities; ≤45 heavy atoms) with MMFF94s fallback — level shown, benchmark vs experiment in the methods modal; animated in 3D with displacement arrows |

## Using it

* **Single mode** — type anything; ⌘K / Ctrl-K focuses the search box. Sections: Identity → Names →
  Structure → Properties → Spectra (sticky navigation). *Export ▾* copies a Markdown report, downloads
  JSON, the numbered structure as SVG, or prints to PDF.
* **Check an identifier** (Names panel) — paste a name from a label, a SMILES from a supplier or an
  InChI from a paper: verdicts are *same compound* / *same skeleton, different stereo* /
  *different compound* / *not interpreted* (nothing is judged unless it was actually parsed).
* Two top-level tabs: **Molecules** (single / batch cross-checking, structures, spectra) and **Elements**.
* **Elements** tab — interactive periodic table (PubChem public-domain data): electron configuration on
  every cell; category, s/p/d/f block or viridis property heat maps, state-at-temperature slider, search,
  detail card, one click to open the element in Glass.
* **Level of theory toggles** in the Spectra panel: NMR basis (nmrshiftdb2 experimental environments vs
  pure increments) and vibrational level (GFN2-xTB vs MMFF94s vs auto), each with a note on what that
  level predicts well.
* Failed inputs get PubChem "did you mean" suggestions when online.

## Rigor model

* **Name → structure** is done by [OPSIN](https://opsin.ch.cam.ac.uk/) 2.9, the reference open-source
  IUPAC parser. If OPSIN cannot parse the input it falls back to a PubChem *name* lookup and **warns
  you** that the structure came from a database synonym, not from grammar.
* **Structure → names**: there is no complete open-source IUPAC *name generator*. Names are taken from
  PubChem and each candidate is **round-tripped through OPSIN**, accepted as `verified` only if the full
  27-character Standard InChIKey matches. Mismatches are shown with a red badge.
* All SMILES variants are regenerated from one RDKit molecule, so they are mutually consistent.
* Atom numbering comes from OPSIN's own locant assignment, correct by construction for the parent hydride.
* Property / spectrum predictions are **estimates** with stated typical errors (Joback ±25–50 K for Tm;
  ¹H ±0.3 ppm; ¹³C ±3–5 ppm; IR ±20 cm⁻¹; UV ±5–10 nm). Use them for sanity checks and assignment support,
  not as data.

## Name generation for novel structures

PubChem cannot name a structure it has never seen. Glass names those with a SMILES→IUPAC model
(knowledgator/SMILES2IUPAC-canonical-base, an MT5 sequence-to-sequence model, Apache-2.0), run
as a sidecar process because it needs PyTorch. Several candidates are generated per structure
in systematic, base and traditional styles, and every one is parsed back with OPSIN: a name is
shown only when its Standard InChIKey is identical to the structure's, and the others are
reported as rejected. An unverified name is never displayed.

## How it runs

Glass is a FastAPI application served by uvicorn behind Caddy, which terminates TLS, as a
systemd service on a Linux VM. OPSIN runs beside it as a persistent JVM worker, and requests
are rate-limited per client IP. The deployment is the one at
[chemistryconsole.glass](https://chemistryconsole.glass).

## Verifying the science

1. **Identity layer** (exact by construction): the InChIKey shown must match PubChem / ChemSpider /
   CAS Common Chemistry for the same structure.
2. **NMR**: `/api/validation` returns the held-out statistics
   (¹³C MAE 2.75 ppm on 21,896 test shifts). Per peak, the table shows basis, σ and n. Spot-check against
   SDBS (https://sdbs.db.aist.go.jp).
3. **IR / vibrations**: the ⓘ modal shows a 12-band benchmark vs gas-phase literature for both levels;
   compare any band with the NIST WebBook (https://webbook.nist.gov).
4. **UV–Vis / photoluminescence / Tm / density**: empirical rules — estimates with stated errors; the
   methods modal lists what each cannot resolve.
5. Every method carries literature references, and `/api/about` lists all third-party components.

## API

Every endpoint returns JSON and needs no key. Rate limits are per IP: 90 requests / 30 s
overall, and 10 / 60 s on the five compute endpoints (`/api/spectra`, `/api/modes`,
`/api/orbital`, `/api/electronic`, `/api/photo`), which answer `429` with `Retry-After`.
The interactive OpenAPI docs are disabled on the public deployment.

* **Identity** — `/api/resolve?q=…&pubchem=true`, `/api/verify?name=…&inchikey=…`,
  `/api/suggest?q=…`, `/api/experimental?cid=…`
* **Structure** — `/api/structure?smiles=…&three_d=true`, `/api/image?smiles=…&w=500&h=380`,
  `/api/cdxml?smiles=…&numbers=false`, `/api/read` (POST: CDXML / MOL upload)
* **Properties and spectra** — `/api/properties?smiles=…`,
  `/api/spectra?smiles=…&solvent=CDCl3&mhz=400`, `/api/solvents`
* **Quantum and vibrational** — `/api/modes?smiles=…&level=auto`,
  `/api/photo?smiles=…&solvent=…`, `/api/electronic?smiles=…`,
  `/api/orbital?smiles=…&which=homo`, `/api/atom_orbital?z=…&n=…&l=…`
* **Elements** — `/api/elements`, `/api/atom?z=…`
* **Provenance** — `/api/methods` (scope and limitations per method), `/api/validation`
  (held-out NMR statistics), `/api/about` (third-party components, source offer),
  `/api/license`, `/api/model` (the NMR model file, 9 MB, offered as nmrshiftdb2's
  share-alike terms require)

## Credits and licence

Glass is free software: Copyright (C) 2026 Gregory Pitch, under the **GNU Affero General Public License
v3 or later** (`LICENSE`, served at `/api/license`). Use it, study it, change it and share it, commercially
or not. If you distribute it, or run a modified copy as a network service, you must offer the corresponding
source under the same licence.

`THIRD_PARTY.md` records every third-party component, its licence and the published methods implemented
here (also served at `/api/about`).

The data files carry their own terms (`data/LICENSE-DATA.md`). In particular the NMR shift model is derived
from **nmrshiftdb2**, whose database licence requires attribution, share-alike on anything derived from it,
and an OSI-approved licence for software that depends on it — part of why Glass is AGPL. The model file is
offered for download at `/api/model`; the raw database is not redistributed here.

## Architecture

```
browser ──HTTP/JSON──▶ FastAPI (app.py) ──▶ core/  (pure functions, no web code)
   │                       │                  ├─ opsin.py       persistent OPSIN JVM worker (name→SMILES/InChIKey/CML)
   │  static/js/*.js       │  lru_cache per   ├─ names.py       resolution, PubChem enrichment, OPSIN verification, locant mapping
   │  (plain modules,      │  (smiles,args)   ├─ structure.py   ACS-style depiction, 3D geometry
   │   no build step)      │                  ├─ properties.py  Joback / Girolami / Crippen
   │                       │                  ├─ spectra.py     NMR / IR / UV simulators (solvent-aware)
   │                       │                  ├─ vibrations.py  GFN2-xTB / MMFF94s normal modes
   │                       │                  ├─ photo.py       photoluminescence + Franck–Condon
   │                       │                  ├─ solvents.py    one solvent table shared by every model
   │                       │                  └─ methods.py     level of theory / limitations registry (→ /api/methods)
   └── Ketcher (vendored)  └─ tests/          pytest round-trip + physics-trend tests
```

Design rules: every simulation is a pure function of (RDKit mol, parameters) so it is cacheable and
testable; the only mutable process state is the OPSIN worker (guarded by a lock) and the LRU caches. Every
method declares its scope and limitations in `core/methods.py`, which the UI shows next to each panel.

Rebuild the NMR model (optional, needs the nmrshiftdb2 SD file from sourceforge, not redistributed here):
`python -m core.hose build data/nmrshiftdb2withsignals.sd`. Measured accuracies are served at `/api/validation`.

Run the tests: `./run_tests.sh` (offline); `CHEMXREF_ONLINE=1 ./run_tests.sh` adds the PubChem checks.

## Layout

```
app.py              FastAPI server
core/opsin.py       persistent OPSIN worker (SMILES / InChIKey / CML with locants)
core/names.py       input classification, resolution, PubChem, verification, locant mapping
core/structure.py   ACS-style depiction, 3D geometry
core/properties.py  descriptors, Joback, Girolami density
core/spectra.py     1H/13C NMR, IR, UV-Vis simulators
core/vibrations.py  normal-mode analysis: GFN2-xTB (tblite) and MMFF94s
core/hose.py        nmrshiftdb2 HOSE-code NMR model
core/photo.py       photoluminescence estimate + Franck–Condon vibronic structure
static/index.html   page shell; static/js/*.js modules
static/vendor/ketcher vendored Ketcher structure editor (standalone, Apache-2.0)
tests/              pytest suite
```
