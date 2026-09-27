"""
CDXML (ChemDraw XML) export.

Writes a single-fragment CDXML document from an RDKit molecule using its 2D depiction
coordinates: heteroatom labels, formal charges, isotopes, wedge/hash stereo bonds, and
optionally the IUPAC locant of each atom as an attached text annotation.

The geometry follows ChemDraw's ACS 1996 document settings (bond length 30 model units
= 0.3 in); CDXML y grows downwards, matching RDKit's SVG coordinate frame.
Reference: CambridgeSoft CDX/CDXML format specification (public web archive).
"""
from __future__ import annotations
from xml.sax.saxutils import escape
from rdkit import Chem
from rdkit.Chem import rdDepictor

BOND_LEN = 30.0          # ChemDraw model units per bond (ACS style: 0.3 in * 100 units/in)


def _wedges(mol: Chem.Mol):
    """RDKit wedge directions for chiral centres (uses the depiction conformer)."""
    Chem.rdDepictor.Compute2DCoords  # noqa - ensure module loaded
    m = Chem.Mol(mol)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    Chem.rdCoordGen  # noqa
    try:
        Chem.WedgeMolBonds(m, m.GetConformer())
    except Exception:
        pass
    out = {}
    for b in m.GetBonds():
        if b.GetBondDir() == Chem.BondDir.BEGINWEDGE:
            out[b.GetIdx()] = "WedgeBegin"
        elif b.GetBondDir() == Chem.BondDir.BEGINDASH:
            out[b.GetIdx()] = "WedgedHashBegin"
    return out


def mol_to_cdxml(mol: Chem.Mol, locants=None, fragments=None, include_numbers=False, title: str | None = None) -> str:
    """
    locants/fragments: per-atom lists as produced by the resolver (numbers shown as
    text annotations when include_numbers is True; substituent locants primed).
    """
    m = Chem.Mol(mol)
    if m.GetNumConformers() == 0:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(m)
    kek = Chem.Mol(m)
    Chem.Kekulize(kek, clearAromaticFlags=True)
    conf = m.GetConformer()
    # scale RDKit coords (median bond ~1.5) to ChemDraw units, flip y (CDXML y is down)
    import statistics
    blens = []
    for b in m.GetBonds():
        p, q = conf.GetAtomPosition(b.GetBeginAtomIdx()), conf.GetAtomPosition(b.GetEndAtomIdx())
        blens.append(((p.x - q.x) ** 2 + (p.y - q.y) ** 2) ** 0.5)
    scale = BOND_LEN / (statistics.median(blens) if blens else 1.5)
    xs = [conf.GetAtomPosition(i).x for i in range(m.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(m.GetNumAtoms())]
    x0, y0 = min(xs), max(ys)
    P = lambda i: (round((conf.GetAtomPosition(i).x - x0) * scale + 40, 2), round((y0 - conf.GetAtomPosition(i).y) * scale + 40, 2))

    wedge = _wedges(m)
    aid = lambda i: i + 1
    parts = []
    parts.append('<?xml version="1.0" encoding="UTF-8"?>')
    parts.append('<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">')
    parts.append('<CDXML CreationProgram="Glass" BondLength="30" LabelFont="3" LabelSize="10" LabelFace="96">'
                 '<fonttable><font id="3" charset="iso-8859-1" name="Arial"/></fonttable>'
                 f'<page WidthPages="1" HeightPages="1">')
    if title:
        parts.append(f'<t p="40 20"><s font="3" size="12">{escape(title)}</s></t>')
    parts.append('<fragment>')
    nid = m.GetNumAtoms() + m.GetNumBonds() + 10
    for a in m.GetAtoms():
        i = a.GetIdx(); x, y = P(i)
        attrs = [f'id="{aid(i)}"', f'p="{x} {y}"']
        if a.GetAtomicNum() != 6:
            attrs.append(f'Element="{a.GetAtomicNum()}"')
        if a.GetTotalNumHs() and a.GetAtomicNum() != 6:
            attrs.append(f'NumHydrogens="{a.GetTotalNumHs()}"')
        if a.GetFormalCharge():
            attrs.append(f'Charge="{a.GetFormalCharge()}"')
        if a.GetIsotope():
            attrs.append(f'Isotope="{a.GetIsotope()}"')
        label = None
        if a.GetAtomicNum() != 6:
            h = a.GetTotalNumHs()
            label = a.GetSymbol() + ("H" if h == 1 else f"H{h}" if h > 1 else "")
        if label:
            parts.append(f'<n {" ".join(attrs)}><t p="{x} {y}" LabelJustification="Auto"><s font="3" size="10" face="96">{escape(label)}</s></t></n>')
        else:
            parts.append(f'<n {" ".join(attrs)}/>')
    for b in kek.GetBonds():
        attrs = [f'id="{m.GetNumAtoms() + b.GetIdx() + 1}"', f'B="{aid(b.GetBeginAtomIdx())}"', f'E="{aid(b.GetEndAtomIdx())}"']
        order = int(b.GetBondTypeAsDouble())
        if order > 1:
            attrs.append(f'Order="{order}"')
        d = wedge.get(b.GetIdx())
        if d:
            attrs.append(f'Display="{d}"')
        parts.append(f'<b {" ".join(attrs)}/>')
    parts.append('</fragment>')
    if include_numbers and locants:
        for a in m.GetAtoms():
            i = a.GetIdx()
            loc = locants[i] if i < len(locants) else None
            frag = fragments[i] if fragments and i < len(fragments) else None
            if not loc or frag == "hetero":
                continue
            txt = str(loc) + ("'" if frag == "sub" else "")
            x, y = P(i)
            colr = ' color="4"' if frag == "sub" else ""
            parts.append(f'<t p="{x + 7} {y - 7}"><s font="3" size="7"{colr}>{escape(txt)}</s></t>')
    parts.append('</page></CDXML>')
    return "\n".join(parts)
