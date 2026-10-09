"""
Identity resolution: ANY input (IUPAC/systematic/trivial name, SMILES, InChI, InChIKey, CAS)
-> one canonical structure -> every name/SMILES representation, each cross-verified.

Rigor model
-----------
* A name is only labelled VERIFIED if OPSIN parses it and the resulting StdInChIKey is
  identical to the StdInChIKey of the resolved structure (full 27-char key => same
  connectivity, stereo, isotopes, and protonation).  Otherwise it is shown but flagged.
* Every SMILES variant is regenerated from one RDKit molecule, so they are mutually
  consistent by construction; the InChIKey of each is re-checked as a guard.
"""
from __future__ import annotations
import re, json, functools, threading, concurrent.futures
from dataclasses import dataclass, field, asdict
import xml.etree.ElementTree as ET
from typing import Optional
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi as rdinchi
from . import opsin

RDLogger.DisableLog("rdApp.*")
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound"
_TIMEOUT = 4
_session = requests.Session()
_session.headers["User-Agent"] = "Glass/1.0 (local productivity tool)"


# ----------------------------------------------------------------------------- helpers
def std_inchikey(mol: Chem.Mol) -> str:
    return rdinchi.MolToInchiKey(mol)


def looks_like_inchi(s: str) -> bool:
    return s.startswith("InChI=")


def looks_like_inchikey(s: str) -> bool:
    return re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", s) is not None


def looks_like_cas(s: str) -> bool:
    return re.fullmatch(r"\d{2,7}-\d{2}-\d", s) is not None


def mol_from_smiles(s: str) -> Optional[Chem.Mol]:
    try:
        m = Chem.MolFromSmiles(s)
        if m is not None and m.GetNumAtoms() > 0:
            return m
    except Exception:
        pass
    return None


@functools.lru_cache(maxsize=4096)
def _get_json(url: str):
    try:
        r = _session.get(url, timeout=_TIMEOUT)
        if r.status_code == 200:
            return r.json()
    except requests.RequestException:
        return None
    return None


# ----------------------------------------------------------------------------- data classes
@dataclass
class NameEntry:
    name: str
    kind: str                 # 'input' | 'iupac' | 'title' | 'synonym' | 'opsin-generated'
    source: str               # 'user' | 'PubChem' | ...
    verified: bool = False    # OPSIN round trip matched InChIKey
    opsin_smiles: Optional[str] = None
    note: str = ""
    locants: Optional[list] = None      # per resolved-mol atom index: locant string or None
    fragments: Optional[list] = None    # per atom: 'parent' | 'sub' | None
    spans: Optional[dict] = None        # atom index (str) -> [[start,end],...] char spans of ITS locant in the name
    fragment_ids: Optional[list] = None # per atom: 0 = parent hydride, 1.. = substituent groups in name order
    regions: Optional[list] = None      # name char ranges -> atoms: [{start,end,atoms,kind,label}], kind parent|group|prefix


@dataclass
class Resolution:
    input: str
    input_kind: str
    smiles_canonical: str
    smiles_canonical_nostereo: str
    smiles_kekule: str
    smiles_input: Optional[str]
    smiles_opsin: Optional[str]
    smiles_pubchem: Optional[str]
    inchi: str
    inchikey: str
    formula: str
    cid: Optional[int]
    generated_note: Optional[str] = None
    cas: Optional[str] = None
    iupac_name: Optional[str] = None     # systematic name (OPSIN-verified)
    common_name: Optional[str] = None    # PubChem title
    names: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    online: bool = True
    substituent: bool = False   # input was a substituent/radical name; structure drawn with attachment points (*)


# ----------------------------------------------------------------------------- locant mapping
_NUM_LOC = re.compile(r"^\d+[a-z]{0,2}$")   # 1, 2, 4a, 8b ...
_PRIMED_NUM = re.compile(r"^\d+[a-z]{0,2}['\u2032]+$")   # 1', 4'' : a numbered atom of a repeated substituent


def _grow(mol, seed, numbered, locants):
    frag = {seed}; used = {locants[seed]}; frontier = [seed]
    while frontier:
        a = frontier.pop()
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            j = nb.GetIdx()
            if j in frag or j not in numbered or locants[j] in used:
                continue
            frag.add(j); used.add(locants[j]); frontier.append(j)
    return frag


def _assign_fragments(mol: Chem.Mol, locants: list, ranks: list, n_groups=None):
    """
    Split numbered atoms into fragments (parent hydride + substituent groups).
    Facts used: locants are unique within a fragment (so growth never crosses into another
    fragment, whose attachment atom re-uses a locant the parent already owns), and OPSIN
    instantiates substituents before the parent, inner groups first, so the parent has the
    highest creation ids and substituent fragments are ordered by their creation ids in the
    same order in which their bracket groups close in the name.
    Returns (kind_per_atom, fragment_id_per_atom, ordered substituent fragment ids).
    kind: 'parent' | 'sub' | 'hetero' | None.  Parent fragment id = 0; substituents 1..n in
    name (closing-bracket) order.
    """
    n = mol.GetNumAtoms()
    remaining = set(i for i in range(n) if locants[i] and _NUM_LOC.match(locants[i]))
    frags = []
    while remaining:
        best, best_key = None, None
        for seed in remaining:
            f = _grow(mol, seed, remaining, locants)
            key = (max(ranks[i] for i in f), sum(ranks[i] for i in f) / len(f), len(f))
            if best_key is None or key > best_key:
                best, best_key = f, key
        frags.append(best); remaining -= best
    kind = [None] * n; fid = [-1] * n
    if not frags:
        return kind, fid, []
    parent = frags[0]
    subs = sorted(frags[1:], key=lambda f: min(ranks[i] for i in f))    # creation order
    # only bracketed substituents own locant tokens in the name; when there are more
    # fragments than bracket groups the extras are un-bracketed simple groups (methyl,
    # chloro-methyl carbons...) -> smallest fragments are dropped from the mapping
    if n_groups is not None:
        mapped = list(subs)
        while len(mapped) > n_groups and mapped:
            smallest = min(mapped, key=len)
            if len(smallest) > 1 and len(mapped) <= n_groups:
                break
            mapped.remove(smallest)
    else:
        mapped = list(subs)
    for i in parent:
        kind[i] = "parent"; fid[i] = 0
    for f in subs:
        for i in f:
            kind[i] = "sub"; fid[i] = -2
    for k, f in enumerate(mapped, 1):
        for i in f:
            fid[i] = k
    for i in range(n):
        if kind[i] is None and locants[i]:
            if mol.GetAtomWithIdx(i).GetSymbol() == "C" and _PRIMED_NUM.match(locants[i]):
                kind[i] = "sub"; fid[i] = -2       # second methyl of a dimethyl: numbered 1', not a heteroatom
                continue
            kind[i] = "hetero"
            nb = [x.GetIdx() for x in mol.GetAtomWithIdx(i).GetNeighbors() if fid[x.GetIdx()] >= 0]
            fid[i] = min((fid[j] for j in nb), default=0)
    return kind, fid, list(range(1, len(mapped) + 1))


_OPEN, _CLOSE = "([{", ")]}"


def _bracket_groups(name: str):
    """[(open_pos, close_pos, depth)] for every bracket group, sorted by closing position."""
    stack, groups = [], []
    for i, ch in enumerate(name):
        if ch in _OPEN:
            stack.append(i)
        elif ch in _CLOSE and stack:
            o = stack.pop(); groups.append((o, i, len(stack)))
    return sorted(groups, key=lambda g: g[1])


def _depth_at(name: str, pos: int) -> int:
    d = 0
    for ch in name[:pos]:
        if ch in _OPEN: d += 1
        elif ch in _CLOSE: d = max(0, d - 1)
    return d


def _spans(name: str, locs: list, fid: list, sub_ids: list) -> dict:
    """
    Per-atom character spans of *its own* locant token in the name.
    Parent-hydride atoms match tokens outside all brackets; substituent fragment k matches
    tokens directly inside its bracket group (k-th group in closing order); tokens nested in
    an inner group belong to that inner group.  Substituents that are not bracketed
    (e.g. 'methyl' in 2-methylbutan-2-ol) carry no own locant in the name -> no span.
    """
    groups = _bracket_groups(name)
    tok = re.compile(r"(?<![A-Za-z0-9'′])([0-9]+[a-z]{0,2}|[A-Z]|[a-z]{3,7})(['′]*)(?![A-Za-z0-9])")
    tokens = [(m.start(), m.end(), m.group(1) + m.group(2).replace("′", "'"), _depth_at(name, m.start())) for m in tok.finditer(name)]
    out: dict = {}
    for i, L in enumerate(locs):
        if not L or fid[i] < 0:
            continue
        if fid[i] == 0:
            region = [(0, len(name), 0)]
        else:
            k = fid[i] - 1
            if k >= len(groups):
                continue
            o, c, d = groups[k]; region = [(o + 1, c, d + 1)]
        hits = [[a, b] for (a, b, val, depth) in tokens if val == L and any(r0 <= a < r1 and depth == rd for r0, r1, rd in region)]
        if hits:
            out[str(i)] = hits
    return out


# ----------------------------------------------------------------------------- name regions
# Unbracketed substituent prefixes we can place on the structure.  Matching is by the prefix's
# own locant(s), so an unrecognised or unlocanted prefix simply yields no region -- it never
# highlights the wrong atoms.  Element sets disambiguate geminal substituents at one locant.
_PREFIX_ELEMENTS = {
    "fluoro": {"F"}, "chloro": {"Cl"}, "bromo": {"Br"}, "iodo": {"I"},
    "oxo": {"O"}, "hydroxy": {"O"}, "methoxy": {"O"}, "ethoxy": {"O"}, "phenoxy": {"O"}, "acetoxy": {"O"},
    "amino": {"N"}, "nitro": {"N"}, "cyano": {"C"}, "azido": {"N"}, "nitroso": {"N"},
    "sulfanyl": {"S"}, "mercapto": {"S"}, "methylsulfanyl": {"S"}, "methylthio": {"S"},
    "sulfo": {"S"}, "methylsulfonyl": {"S"}, "phosphono": {"P"},
}
_CARBON_PREFIXES = [
    "trifluoromethyl", "cyclopropyl", "cyclobutyl", "cyclopentyl", "cyclohexyl", "tert-butyl", "sec-butyl",
    "isopropyl", "isobutyl", "methyl", "ethyl", "propyl", "butyl", "pentyl", "hexyl", "heptyl", "octyl",
    "vinyl", "allyl", "ethynyl", "ethenyl", "propenyl", "phenyl", "benzyl", "naphthyl", "tolyl",
    "acetyl", "formyl", "benzoyl", "carboxy", "carbamoyl", "methoxycarbonyl", "ethoxycarbonyl",
]
_STEMS = sorted(list(_PREFIX_ELEMENTS) + _CARBON_PREFIXES, key=len, reverse=True)
# Characteristic-group suffixes we can place: element sets select the group's atoms among the
# parent's heteroatom-labelled atoms (those OPSIN gives O, O', N, alpha ... rather than a number).
_SUFFIX_ELEMENTS = {
    "carboxylic acid": {"C", "O"}, "carbothioic acid": {"C", "O", "S"}, "sulfonic acid": {"S", "O"},
    "carbaldehyde": {"C", "O"}, "carboxamide": {"C", "N", "O"}, "carbonitrile": {"C", "N"},
    "oic acid": {"C", "O"}, "ic acid": {"C", "O"}, "aldehyde": {"C", "O"}, "nitrile": {"C", "N"},
    "amide": {"C", "N", "O"}, "amine": {"N"}, "thiol": {"S"}, "one": {"O"}, "ol": {"O"}, "al": {"C", "O"},
}
_MULT = r"(?:tetrakis|pentakis|hexakis|tris|bis|tetra|penta|hexa|hepta|octa|tri|di)?"
_ONE_LOC = r"(?!tert|sec|iso|neo|cis|trans)(?:[0-9]+[a-z]{0,2}|[A-Z]|[a-z]{3,7})['\u2032]*"
_BRACKET_LOC_RE = re.compile(r"(?:" + _ONE_LOC + r"(?:," + _ONE_LOC + r")*)-(?= )", re.IGNORECASE)
_PREFIX_RE = re.compile(r"(?:(" + _ONE_LOC + r"(?:," + _ONE_LOC + r")*)-)?(" + _MULT + r")(" + "|".join(map(re.escape, _STEMS)) + r")-?", re.IGNORECASE)   # PubChem titles capitalise: 2-Methyl-
_SUFFIX_RE = re.compile(r"(?:-?(" + _ONE_LOC + r"(?:," + _ONE_LOC + r")*)-)?(" + _MULT + r")("
                        + "|".join(map(re.escape, sorted(_SUFFIX_ELEMENTS, key=len, reverse=True))) + r")$", re.IGNORECASE)


def _blank_nested(text: str) -> str:
    """Replace bracket groups with spaces (same length) so offsets are preserved."""
    out, depth = [], 0
    for ch in text:
        if ch in _OPEN:
            depth += 1; out.append(" ")
        elif ch in _CLOSE:
            depth = max(0, depth - 1); out.append(" ")
        else:
            out.append(ch if depth == 0 else " ")
    return "".join(out)


def _attached_group(mol, seed: int, k: int, fid: list, frags: list, elements) -> set:
    """Atoms of the substituent hanging off `seed` that is not part of numbered fragment k."""
    def outside(j):
        return fid[j] in (-2, -1) or (frags[j] == "hetero" and fid[j] == k)
    out = set()
    for nb in mol.GetAtomWithIdx(seed).GetNeighbors():
        j = nb.GetIdx()
        if not outside(j) or (elements and nb.GetSymbol() not in elements):
            continue
        grp, front = {j}, [j]
        while front:
            a = front.pop()
            for nb2 in mol.GetAtomWithIdx(a).GetNeighbors():
                b = nb2.GetIdx()
                if b != seed and b not in grp and outside(b):
                    grp.add(b); front.append(b)
        out |= grp
    return out


def _regions(name: str, mol, locs: list, fid: list, frags: list) -> list:
    """
    Map character ranges of the name onto atoms.  Three kinds:
      group  -- a bracketed substituent (k-th bracket in closing order <-> fragment id k),
                including anything nested inside it;
      prefix -- an unbracketed substituent prefix placed by its own locant(s);
      parent -- the parent hydride word (whatever is left at depth 0 once prefixes are peeled).
    """
    groups = _bracket_groups(name)
    regions = []
    n = len(locs)
    # bracketed groups, with nested groups folded in
    for k, (o, c, _d) in enumerate(groups, 1):
        atoms = {i for i in range(n) if fid[i] == k}
        for j, (oj, cj, _dj) in enumerate(groups, 1):
            if o < oj and cj < c:
                atoms |= {i for i in range(n) if fid[i] == j}
        regions.append({"start": o, "end": c + 1, "atoms": sorted(atoms), "kind": "group", "label": name[o:c + 1]})

    def peel(text_start: int, text_end: int, k: int):
        """Peel recognised prefixes from name[text_start:text_end] at fragment k; return end of last prefix."""
        seg = _blank_nested(name[text_start:text_end])
        pos, last = 0, 0
        while pos < len(seg):
            mb = _BRACKET_LOC_RE.match(seg, pos)
            if seg[pos] == " " or mb:         # a nested bracket group (with its locants): skip it and a following hyphen
                pos = mb.end() if mb else pos
                while pos < len(seg) and seg[pos] == " ":
                    pos += 1
                if pos < len(seg) and seg[pos] == "-":
                    pos += 1
                last = pos
                continue
            m = _PREFIX_RE.match(seg, pos)
            if not m:
                break
            loc_run, _mult, stem = m.group(1), m.group(2), m.group(3).lower()
            atoms = set()
            if loc_run:
                for L in loc_run.replace("\u2032", "'").split(","):
                    for a in range(n):
                        if fid[a] == k and locs[a] == L:
                            atoms |= _attached_group(mol, a, k, fid, frags, _PREFIX_ELEMENTS.get(stem, {"C"}))
            if atoms:
                regions.append({"start": text_start + m.start(), "end": text_start + m.end() - (1 if m.group(0).endswith("-") else 0),
                                "atoms": sorted(atoms), "kind": "prefix", "label": name[text_start + m.start():text_start + m.end()].rstrip("-")})
            pos = last = m.end()
        return text_start + last

    # prefixes inside each bracket group, at that group's fragment id
    for k, (o, c, _d) in enumerate(groups, 1):
        peel(o + 1, c, k)
    # depth-0 prefixes, then the parent word is what remains
    start = peel(0, len(name), 0)
    claimed = set()
    for g in regions:
        if g["kind"] == "prefix" and _depth_at(name, g["start"]) == 0:
            claimed |= set(g["atoms"])
    parent_atoms = sorted(i for i in range(n) if fid[i] == 0 and i not in claimed)
    if parent_atoms and start < len(name):
        regions.append({"start": start, "end": len(name), "atoms": parent_atoms, "kind": "parent", "label": name[start:]})
        # the characteristic-group suffix, as a region nested inside the parent word
        ms = _SUFFIX_RE.search(name, start)
        if ms:
            loc_run, _mult, suffix = ms.group(1), ms.group(2), ms.group(3).lower()
            elements = _SUFFIX_ELEMENTS[suffix]
            hetero0 = {i for i in range(n) if fid[i] == 0 and frags[i] == "hetero" and i not in claimed
                       and mol.GetAtomWithIdx(i).GetSymbol() in elements}
            if loc_run:
                atoms = set()
                for L in loc_run.replace("\u2032", "'").split(","):
                    for a in range(n):
                        if fid[a] == 0 and locs[a] == L and frags[a] != "hetero":
                            atoms |= _attached_group(mol, a, 0, fid, frags, elements) & hetero0
                            # grow through the group's other heteroatoms (the two O of a carboxyl)
                            atoms |= {j for j in hetero0 if any(nb.GetIdx() in atoms for nb in mol.GetAtomWithIdx(j).GetNeighbors())}
            else:
                atoms = hetero0
            if atoms and "C" in elements:
                # on a chain the suffix carbon is a numbered atom (ethanoic acid: C1); include it, but not a ring carbon
                atoms |= {nb.GetIdx() for i in list(atoms) for nb in mol.GetAtomWithIdx(i).GetNeighbors()
                          if fid[nb.GetIdx()] == 0 and frags[nb.GetIdx()] == "parent" and nb.GetSymbol() == "C" and not nb.IsInRing()}
            if atoms:
                s0 = ms.start() + (1 if name[ms.start()] == "-" else 0)
                regions.append({"start": s0, "end": len(name), "atoms": sorted(atoms), "kind": "suffix", "label": name[s0:]})
    # a bracketed group is the whole substituent: fold in its nested prefixes and take its own locants
    for g in regions:
        if g["kind"] != "group":
            continue
        inner = set(g["atoms"])
        for h in regions:
            if h["kind"] == "prefix" and g["start"] <= h["start"] and h["end"] <= g["end"]:
                inner |= set(h["atoms"])
        g["atoms"] = sorted(inner)
        mloc = re.search(r"(?:" + _ONE_LOC + r"(?:," + _ONE_LOC + r")*)-$", name[:g["start"]])
        if mloc and _depth_at(name, mloc.start()) == _depth_at(name, g["start"]):
            g["start"] = mloc.start(); g["label"] = name[g["start"]:g["end"]]
    return regions


def cap_radicals(m: Chem.Mol) -> Chem.Mol:
    """Replace each unpaired electron with a bond to a dummy atom (*): a substituent with explicit attachment points."""
    rw = Chem.RWMol(m)
    for a in list(rw.GetAtoms()):
        n = a.GetNumRadicalElectrons()
        if not n:
            continue
        for _ in range(n):
            d = rw.AddAtom(Chem.Atom(0))
            rw.AddBond(a.GetIdx(), d, Chem.BondType.SINGLE)
        a.SetNumRadicalElectrons(0)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


def _regions_from_parse(name: str, parse: dict, id2ti: dict, mol, locs: list) -> Optional[list]:
    """
    Name regions from OPSIN's own parse tree (core/opsin.name_to_parse).  Brackets are paired
    into `group` regions (with the locant that precedes the bracket), every `substituent`
    element becomes a `prefix` (its bracket tokens excluded), the `root` the `parent`, and a
    root suffix that created atoms of its own a nested `suffix`.  Offsets come from walking the
    raw tokens in order; a short gap between tokens (the elided 'o' of benz|o|ic acid) folds
    into the following token.  Returns None when the tree cannot be aligned to `name`, and the
    caller then falls back to the vocabulary-based _regions().
    """
    pre = parse.get("pre") or name
    if len(pre) != len(name):
        return None
    toks = parse["tokens"]
    by_i = {t["i"]: t for t in toks}
    children: dict = {}
    for t in toks:
        children.setdefault(t["p"], []).append(t["i"])
    leaves = [t for t in toks if t["v"] is not None]
    low = pre.lower()
    cursor, prev_start, span = 0, 0, {}
    for t in leaves:
        k = low.find(t["v"].lower(), cursor)
        if k < 0:                                   # a token may overlap the one before it (bi|biphenyl)
            k = low.find(t["v"].lower(), prev_start)
            if k < 0:
                return None
        # an elided letter belongs to the word that follows it (benz|o|ic acid), never to a locant or hyphen (purin|e|-2,6-dione)
        fold = 0 < k - cursor <= 2 and t["el"] not in ("locant", "multiplier", "hyphen", "openbracket", "closebracket", "structuralOpenBracket", "structuralCloseBracket")
        span[t["i"]] = (cursor if fold else k, k + len(t["v"]))
        prev_start, cursor = k, max(cursor, k + len(t["v"]))

    extra_by_v: dict = {}
    for ex in parse.get("extras", []):
        extra_by_v.setdefault(ex["v"], []).extend(ex["atoms"])

    def tok_atoms(t):
        ids = list(t["atoms"]) + (extra_by_v.get(t["v"], []) if t["el"] == "group" else [])
        return {id2ti[a] for a in ids if a in id2ti}

    def desc_leaves(i):
        t = by_i[i]
        if t["v"] is not None:
            return [t]
        out = []
        for c in children.get(i, []):
            out += desc_leaves(c)
        return out

    # brackets: pair them in document order; the locant just before an opening bracket belongs to it
    order = {t["i"]: n for n, t in enumerate(leaves)}
    bracket_locants = set()
    pairs, stack = [], []
    for n, t in enumerate(leaves):
        if t["el"] in ("openbracket", "structuralOpenBracket"):
            stack.append(n)
        elif t["el"] in ("closebracket", "structuralCloseBracket") and stack:
            pairs.append((stack.pop(), n))
    regions = []
    for o, c in pairs:
        if leaves[o]["el"] == "structuralOpenBracket":   # ring assemblies ([1,1'-biphenyl]) are part of the parent word
            continue
        s0 = span[leaves[o]["i"]][0]
        if o > 0 and leaves[o - 1]["el"] == "locant":
            s0 = span[leaves[o - 1]["i"]][0]; bracket_locants.add(leaves[o - 1]["i"])
        s1 = span[leaves[c]["i"]][1]
        atoms = set()
        for t in leaves[o + 1:c]:
            atoms |= tok_atoms(t)
        if atoms:
            regions.append({"start": s0, "end": s1, "atoms": sorted(atoms), "kind": "group", "label": name[s0:s1]})

    skip = {"openbracket", "closebracket", "structuralOpenBracket", "structuralCloseBracket"}
    for t in toks:
        if t["v"] is not None or t["el"] not in ("substituent", "root"):
            continue
        ls = desc_leaves(t["i"])
        if t["el"] == "substituent":
            ls = [l for l in ls if l["el"] not in skip and l["i"] not in bracket_locants]
            while ls and ls[-1]["el"] == "hyphen":
                ls.pop()
            if not ls:
                continue
            atoms = set()
            for l in ls:
                atoms |= tok_atoms(l)
            s0, s1 = span[ls[0]["i"]][0], span[ls[-1]["i"]][1]
            if atoms and s1 > s0:
                regions.append({"start": s0, "end": s1, "atoms": sorted(atoms), "kind": "prefix", "label": name[s0:s1]})
            continue
        # root: the parent word, through to the end of the name (a trailing elided 'e' has no token)
        atoms = set()
        for l in ls:
            atoms |= tok_atoms(l)
        if not atoms or not ls:
            continue
        s0, s1 = span[ls[0]["i"]][0], len(name)
        regions.append({"start": s0, "end": s1, "atoms": sorted(atoms), "kind": "parent", "label": name[s0:s1]})
        kids = children.get(t["i"], [])
        group = next((by_i[c] for c in kids if by_i[c]["el"] == "group"), None)
        sufs = [by_i[c] for c in kids if by_i[c]["el"] == "suffix"]
        if not group or not sufs:
            continue
        gm = Chem.MolFromSmiles(group["g"]) if group.get("g") else None
        n_own = gm.GetNumAtoms() if gm is not None else 0
        gat = sorted(a for a in group["atoms"] if a in id2ti)            # ids ascend in creation order
        suffix_ids = gat[n_own:] if n_own and len(gat) > n_own else []
        satoms = {id2ti[a] for a in suffix_ids}
        if not satoms:
            # no usable group template: the suffix's atoms are the group atoms OPSIN labelled with a
            # non-numeric locant (O, N, alpha) rather than a chain/ring number
            satoms = {id2ti[a] for a in gat if locs[id2ti[a]] and not _NUM_LOC.match(locs[id2ti[a]]) and not _PRIMED_NUM.match(locs[id2ti[a]])}
        if not satoms:
            continue
        suf = sufs[-1]
        sv = suf["v"].lower()
        if "C" in (_SUFFIX_ELEMENTS.get(sv) or _SUFFIX_ELEMENTS.get("o" + sv) or set()):
            # a chain carbon bearing the suffix heteroatoms is part of the group (ethanoic acid: C1); a ring carbon is not
            satoms |= {nb.GetIdx() for a in list(satoms) for nb in mol.GetAtomWithIdx(a).GetNeighbors()
                       if nb.GetSymbol() == "C" and not nb.IsInRing() and nb.GetIdx() in atoms}
        s_start = span[suf["i"]][0]
        j = kids.index(suf["i"]) - 1
        while j >= 0 and by_i[kids[j]]["el"] in ("locant", "multiplier", "hyphen") and kids[j] in span:
            s_start = span[kids[j]][0]; j -= 1
        lab = name[s_start:s1].lstrip("-")
        regions.append({"start": s1 - len(lab), "end": s1, "atoms": sorted(satoms), "kind": "suffix", "label": lab})
    return regions or None


def attach_locants(entry: NameEntry, target: Chem.Mol, cml_el) -> None:
    """Fill entry.locants / fragments / spans / regions for entry.name.
    Prefers OPSIN's parse sidecar (one persistent JVM, its own CML); `cml_el` from the CLI is the fallback."""
    parse = opsin.name_to_parse(entry.name)
    if parse and parse.get("cml"):
        try:
            root = ET.fromstring(parse["cml"])
            cml_el = root if root.tag.endswith("molecule") else next(root.iter(opsin._CML_NS + "molecule"))
        except Exception:
            parse = None
    if cml_el is None:
        return
    om, av, ranks = opsin.cml_to_mol(cml_el)
    if om.GetNumAtoms() != target.GetNumAtoms():
        return
    match = target.GetSubstructMatch(om, useChirality=True)
    if len(match) != target.GetNumAtoms():
        match = target.GetSubstructMatch(om)
    if len(match) != target.GetNumAtoms():
        return
    kind_o, fid_o, sub_ids = _assign_fragments(om, av, ranks, len(_bracket_groups(entry.name)))
    n = target.GetNumAtoms()
    locs = [None] * n; frags = [None] * n; fid = [-1] * n
    for oi, ti in enumerate(match):
        locs[ti] = av[oi]; frags[ti] = kind_o[oi]; fid[ti] = fid_o[oi]
    entry.locants, entry.fragments, entry.fragment_ids = locs, frags, fid
    entry.spans = _spans(entry.name, locs, fid, sub_ids)
    entry.regions = None
    if parse:
        try:
            entry.regions = _regions_from_parse(entry.name, parse, {ranks[oi]: ti for oi, ti in enumerate(match)}, target, locs)
        except Exception:
            entry.regions = None
    if not entry.regions:
        entry.regions = _regions(entry.name, target, locs, fid, frags)


# ----------------------------------------------------------------------------- style ranking
_STYLE_PENALTIES = [
    (r",\s", 40), (r"\b[omp]-", 25), (r"\b(iso|sec-|tert-|neo)", 12),
    (r"(oxidanyl|azanyl|sulfanyl|λ|lambda|anyl\b|chloranyl|fluoranyl|bromanyl|iodanyl)", 30),
    (r"cyclohexa-1,3,5-triene", 35), (r"\]methanoic|\]methanenitrile|\]methanal", 20),
    (r"\b(alpha|beta|gamma|α|β|γ)\b", 15), (r"\d[a-z]?-(ol|one|al|oic|amine|amide)\b", 8),
    (r"acetic|formic|propionic|butyric", 6),
]


_RETAINED = re.compile(r"^(methanol|ethanol|methanal|ethene|ethyne|phenol|aniline|toluene|acet(ic|one|aldehyde)|formic acid|form(aldehyde|amide)|benzoic acid|benzaldehyde|acetic acid|pyridine|furan|pyrrole|thiophene|naphthalene|anthracene|styrene|urea|oxalic acid|glycerol|catechol|resorcinol|hydroquinone|purine|indole)$", re.I)


def style_score(name: str) -> int:
    """Approximate IUPAC-2013 style preference for RANKING verified names (lower = preferred).
    A documented heuristic, not a P-rule engine: every ranked name is an exact synonym
    (identical InChIKey); the score only orders presentation."""
    n = name.strip(); score = 0
    if _RETAINED.match(n):
        return -50
    for pat, pen in _STYLE_PENALTIES:
        if re.search(pat, n, re.I):
            score += pen
    if not re.search(r"\d", n) and len(n) > 10:
        score += 10
    score += max(0, len(n) - 25) // 10
    return score


# ----------------------------------------------------------------------------- verification
def verify_name(entry: NameEntry, inchikey: str) -> None:
    smi = opsin.name_to_smiles(entry.name)
    entry.opsin_smiles = smi
    if not smi:
        entry.verified = False
        entry.note = entry.note or "not a systematic name (OPSIN cannot parse); database synonym only"
        return
    key = opsin.name_to_inchikey(entry.name)
    if key == inchikey:
        entry.verified = True
        entry.note = "OPSIN round-trip: identical StdInChIKey"
    elif key and key.split("-")[0] == inchikey.split("-")[0]:
        entry.verified = False
        entry.note = "same skeleton but stereo/isotope/charge layer differs (InChIKey block 2/3 mismatch)"
    else:
        entry.verified = False
        entry.note = "OPSIN parses this name to a DIFFERENT structure; do not trust"


# ----------------------------------------------------------------------------- PubChem
def pubchem_by_inchikey(key: str):
    j = _get_json(f"{PUBCHEM}/compound/inchikey/{key}/property/IUPACName,SMILES,ConnectivitySMILES,Title/JSON")
    if not j:
        return None
    return j["PropertyTable"]["Properties"][0]


def pubchem_by_name(name: str):
    j = _get_json(f"{PUBCHEM}/compound/name/{requests.utils.quote(name)}/property/IUPACName,SMILES,ConnectivitySMILES,InChIKey,Title/JSON")
    if not j:
        return None
    return j["PropertyTable"]["Properties"][0]


def pubchem_cas(cid: int) -> Optional[str]:
    j = _get_json(f"{PUBCHEM}/compound/cid/{cid}/synonyms/JSON")
    if not j:
        return None
    for s in j["InformationList"]["Information"][0].get("Synonym", []):
        if re.fullmatch(r"\d{2,7}-\d{2}-\d", s):
            return s
    return None


def pubchem_by_smiles(smiles: str):
    j = _get_json(f"{PUBCHEM}/compound/smiles/{requests.utils.quote(smiles, safe='')}/property/IUPACName,SMILES,ConnectivitySMILES,InChIKey,Title/JSON")
    if not j:
        return None
    return j["PropertyTable"]["Properties"][0]


def pubchem_synonyms(cid: int, limit=40) -> list[str]:
    j = _get_json(f"{PUBCHEM}/compound/cid/{cid}/synonyms/JSON")
    if not j:
        return []
    syns = j["InformationList"]["Information"][0].get("Synonym", [])
    out = []
    for s in syns:
        if re.fullmatch(r"[A-Z0-9\-]{6,}", s) and "-" in s and s.count("-") >= 2:
            continue      # registry numbers, InChIKeys, etc.
        if re.fullmatch(r"\d{2,7}-\d{2}-\d", s):
            continue
        if s.lower().startswith(("dtxsid", "chebi", "chembl", "unii", "schembl", "hsdb", "nsc", "ec ", "einecs", "brn ", "mfcd", "zinc", "akos", "bdbm", "cas-", "dtxcid", "gtpl")):
            continue
        out.append(s)
        if len(out) >= limit:
            break
    return out


def pubchem_experimental(cid: int) -> dict:
    """Experimental property strings from PubChem PUG-View (melting/boiling point, density...)."""
    out = {}
    for heading, key in (("Melting Point", "mp"), ("Boiling Point", "bp"), ("Density", "density"), ("Solubility", "solubility"), ("LogP", "logp")):
        j = _get_json(f"{PUBCHEM_VIEW}/{cid}/JSON?heading={requests.utils.quote(heading)}")
        if not j:
            continue
        vals = []
        def walk(sec):
            for s in sec.get("Section", []):
                walk(s)
            for info in sec.get("Information", []):
                v = info.get("Value", {})
                for sw in v.get("StringWithMarkup", []):
                    vals.append(sw["String"])
                if "Number" in v:
                    vals.append(" ".join(str(x) for x in v["Number"]) + " " + v.get("Unit", ""))
        walk(j.get("Record", {}))
        if vals:
            out[key] = vals[:6]
    return out


# ----------------------------------------------------------------------------- main entry
def resolve(text: str, want_pubchem: bool = True) -> Resolution:
    if len(text) > 10_000:
        raise ValueError("input too long (> 10,000 characters)")
    raw = text.strip()
    raw = re.sub(r"^(InChIKey|CAS)\s*[=:]\s*", "", raw, flags=re.I)
    if not raw:
        raise ValueError("empty input")
    warnings: list[str] = []
    mol = None; kind = ""; smiles_input = None; smiles_opsin = None; pc = None; online = True; substituent = False

    # 1. InChI / InChIKey / CAS
    if looks_like_inchi(raw):
        mol = Chem.MolFromInchi(raw); kind = "InChI"
        if mol is None:
            raise ValueError("Invalid InChI")
    elif looks_like_inchikey(raw) or looks_like_cas(raw):
        kind = "InChIKey" if looks_like_inchikey(raw) else "CAS RN"
        pc = pubchem_by_inchikey(raw) if kind == "InChIKey" else pubchem_by_name(raw)
        if not pc:
            raise ValueError(f"{kind} not found in PubChem (or offline)")
        mol = mol_from_smiles(pc.get("SMILES") or pc.get("ConnectivitySMILES"))
    else:
        # 2. systematic name via OPSIN
        smi = opsin.name_to_smiles(raw)
        if smi:
            mol = mol_from_smiles(smi); kind = "IUPAC / systematic name"; smiles_opsin = smi
        # 3. SMILES
        if mol is None:
            m = mol_from_smiles(raw)
            if m is not None and (re.search(r"[=#()\[\]@/\\0-9]", raw) or re.fullmatch(r"[A-Za-z]{1,12}", raw) is None or raw.lower() in ("c", "cc", "co", "ccc", "cco", "occo", "cn", "n", "o", "s", "cs")):
                mol = m; kind = "SMILES"; smiles_input = raw
        # 3b. substituent name (methyl, phenyl, acetyl): OPSIN with radicals allowed, drawn with attachment points
        if mol is None:
            parse = opsin.name_to_parse(raw)
            if parse and parse.get("smiles"):
                m = mol_from_smiles(parse["smiles"])
                if m is not None and any(a.GetNumRadicalElectrons() for a in m.GetAtoms()):
                    mol = cap_radicals(m); kind = "substituent name (OPSIN, radicals allowed)"; smiles_opsin = parse["smiles"]
                    substituent = True
                    warnings.append("Substituent name: OPSIN parsed it as a radical, shown here with attachment points (*). "
                                    "Properties are those of the radical fragment; spectra are not computed for fragments.")
        # 4. trivial / trade name via PubChem
        if mol is None and want_pubchem:
            pc = pubchem_by_name(raw)
            if pc:
                mol = mol_from_smiles(pc.get("SMILES") or pc.get("ConnectivitySMILES"))
                kind = "common / trade name (PubChem lookup)"
                warnings.append("Input was NOT a systematic name; structure came from a PubChem name lookup. Check the structure.")
        if mol is None:
            reason = opsin.failure_reason(raw)
            raise ValueError(f"Could not interpret input as a name, SMILES, InChI or InChIKey. OPSIN says: {reason}")

    # --- canonical structure -----------------------------------------------------------
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    can = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    mol = Chem.MolFromSmiles(can)          # normalise atom order to canonical
    nostereo = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
    kek = Chem.Mol(mol); Chem.Kekulize(kek, clearAromaticFlags=True)
    kekule = Chem.MolToSmiles(kek, kekuleSmiles=True, isomericSmiles=True, canonical=True)
    inchi = rdinchi.MolToInchi(mol)
    key = std_inchikey(mol)
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    formula = CalcMolFormula(mol)

    # --- names -------------------------------------------------------------------------
    names: list[NameEntry] = []
    seen: set[str] = set()

    def add(e: NameEntry):
        k = e.name.strip().lower()
        if not k or k in seen:
            return
        seen.add(k); names.append(e)

    if kind == "IUPAC / systematic name":
        add(NameEntry(raw, "input", "user"))

    cid = None; smiles_pubchem = None; cas = None
    if want_pubchem:
        if pc is None:
            pc = pubchem_by_inchikey(key) or pubchem_by_smiles(can)
            if pc is None and nostereo != can:
                pc = pubchem_by_smiles(nostereo)
                if pc:
                    warnings.append("PubChem only has this compound without the specified stereochemistry; database names below refer to the stereo-unspecified entry (they will show as 'unverified' unless stereo-complete).")
        if pc:
            cid = pc.get("CID"); smiles_pubchem = pc.get("SMILES") or pc.get("ConnectivitySMILES")
            # CAS and synonyms are independent PubChem calls; overlap them
            with concurrent.futures.ThreadPoolExecutor(max_workers=2) as _ex:
                _fcas = _ex.submit(pubchem_cas, cid)
                _fsyn = _ex.submit(pubchem_synonyms, cid)
                cas = _fcas.result()
                _syns = _fsyn.result()
            if pc.get("IUPACName"):
                add(NameEntry(pc["IUPACName"], "iupac", "PubChem systematic name"))
            if pc.get("Title"):
                add(NameEntry(pc["Title"], "title", "PubChem preferred/common name"))
            for s in _syns:
                add(NameEntry(s, "synonym", "PubChem synonym"))
        else:
            online = _get_json(f"{PUBCHEM}/compound/cid/2244/property/Title/JSON") is not None
            if not online:
                warnings.append("PubChem unreachable, working offline; only OPSIN/RDKit-derived data shown.")
            else:
                warnings.append("Structure not found in PubChem; no database names available. Names shown are what you entered / OPSIN can verify.")
    else:
        online = False

    # verify every name, attach locants for verified systematic names
    for e in names:
        verify_name(e, key)
    verified = [e for e in names if e.verified]
    cml = opsin.names_to_cml([e.name for e in verified]) if verified and not opsin.shim_available() else {}
    for e in verified:
        el = cml.get(" ".join(e.name.split()))
        if el is None and not opsin.shim_available():
            continue
        try:
            attach_locants(e, mol, el)
        except Exception as ex:  # never let numbering break resolution
            e.note += f" (locant mapping failed: {ex})"

    # order: verified first; among verified, approximate Blue-Book style preference
    names.sort(key=lambda e: (not e.verified, style_score(e.name) if e.verified else 0,
                              {"generated": 0, "iupac": 1, "input": 2, "title": 3, "synonym": 4}.get(e.kind, 5), len(e.name)))
    best_sys = next((e.name for e in names if e.verified), None)

    iupac = best_sys or next((e.name for e in names if e.kind == "iupac" and e.verified), None) or next((e.name for e in names if e.kind == "input" and e.verified), None)
    generated_note = None
    if iupac is None:
        from . import namer
        if namer.available():
            g = namer.generate(can)
            if "candidates" in g:
                rejected, accepted = [], []
                for c in g["candidates"]:
                    e = NameEntry(c["name"], "generated", f"{g['generator']} · {c['style'].lower()} style")
                    verify_name(e, key)
                    if e.verified:
                        if e.name.strip().lower() not in seen:
                            accepted.append(e); seen.add(e.name.strip().lower())
                    else:
                        rejected.append(f"“{c['name']}” ({c['style'].lower()}: {e.note.split('; ')[0]})")
                # all accepted names are correct (identical InChIKey); prefer conventional style:
                # penalise "oxidanyl/fluoranyl/…-anyl" and "cyclohexa-1,3,5-triene"-type spellings, then shorter names
                def badness(e):
                    n = e.name.lower()
                    return (sum(tok in n for tok in ("oxidanyl", "anyl", "cyclohexa-1,3,5-triene", "lambda", "azanyl", "sulfanyl", "]methanoic", "]methanenitrile", "]methanal", "ethanoic acid", "ethanoate")), 0 if "base" in e.source else 1, len(n))
                accepted.sort(key=badness)
                for e in reversed(accepted):
                    names.insert(0, e)
                if accepted:
                    e = accepted[0]; iupac = e.name
                    cml = {} if opsin.shim_available() else opsin.names_to_cml([e.name])
                    if cml.get(" ".join(e.name.split())) or opsin.shim_available():
                        try:
                            attach_locants(e, mol, cml.get(" ".join(e.name.split())))
                        except Exception:
                            pass
                if iupac:
                    generated_note = "generated systematic name, verified by OPSIN round-trip (identical Standard InChIKey)" + (f"; {len(rejected)} other candidate(s) rejected" if rejected else "")
                else:
                    generated_note = "no candidate name survived OPSIN verification, so none are shown: " + "; ".join(rejected)
            else:
                generated_note = "no verified systematic name available for this structure"
        else:
            generated_note = "no verified systematic name available for this structure"
    common = next((e.name for e in names if e.kind == "title"), None)
    if substituent:
        # a substituent name has no PubChem record and no systematic name to generate: the input is the name
        names = [NameEntry(name=raw, kind="input", source="user", verified=True,
                           note="OPSIN (radicals allowed): substituent; attachment points shown as *")]
        iupac, common = raw, raw
        generated_note = "substituent name; attachment points shown as *"
        warnings = [w for w in warnings if "not found in PubChem" not in w]
    return Resolution(
        input=raw, input_kind=kind, cas=cas, iupac_name=iupac, common_name=common, generated_note=generated_note,
        smiles_canonical=can, smiles_canonical_nostereo=nostereo, smiles_kekule=kekule,
        smiles_input=smiles_input, smiles_opsin=smiles_opsin, smiles_pubchem=smiles_pubchem,
        inchi=inchi, inchikey=key, formula=formula, cid=cid,
        names=names, warnings=warnings, online=online, substituent=substituent,
    )


def resolution_to_dict(r: Resolution) -> dict:
    d = asdict(r)
    return d
