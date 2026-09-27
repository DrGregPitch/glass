"""
2D depiction (ACS 1996 journal style) and 3D geometry (bond lengths / angles).
"""
from __future__ import annotations
import re, math
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

rdDepictor.SetPreferCoordGen(True)   # CoordGen gives cleaner journal-like layouts


def depict_svg(mol: Chem.Mol, width=520, height=380, show_h_on_hetero=True, font=14, bond=42) -> dict:
    """
    Returns {'svg': str, 'coords': [[x,y],...], 'label_dirs': [[dx,dy],...]}.
    coords are in SVG pixel space so the UI can overlay locant numbers & hover targets.
    Style follows ACS 1996 conventions (bond length 14.4 pt, 0.6 pt lines, Arial 10 pt,
    120° geometry from the 2D coordinate generator).
    """
    m = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(m)
    rdDepictor.StraightenDepiction(m)
    d = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = d.drawOptions()
    try:
        rdMolDraw2D.SetACS1996Mode(opts, 14.4 * 1.9)   # scale up a bit for screen
    except Exception:
        pass
    opts.addStereoAnnotation = True
    opts.explicitMethyl = False
    opts.fixedBondLength = bond
    opts.padding = 0.12
    opts.useBWAtomPalette()
    opts.bondLineWidth = 1.6
    opts.minFontSize = font; opts.maxFontSize = font
    opts.fontFile = ""      # default (Arial-like) font
    d.DrawMolecule(m)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    svg = re.sub(r"<\?xml[^>]*>\s*", "", svg)
    svg = re.sub(r"<rect[^>]*fill:#FFFFFF[^>]*/>", "", svg, count=1)   # transparent background
    coords, dirs = [], []
    for a in m.GetAtoms():
        p = d.GetDrawCoords(a.GetIdx())
        coords.append([round(p.x, 1), round(p.y, 1)])
    # direction for label placement: away from neighbours
    for a in m.GetAtoms():
        i = a.GetIdx(); nx = ny = 0.0
        for nb in a.GetNeighbors():
            j = nb.GetIdx()
            nx += coords[j][0] - coords[i][0]; ny += coords[j][1] - coords[i][1]
        if abs(nx) + abs(ny) < 1e-6:
            # e.g. 2 collinear neighbours or 3 symmetric ones: use perpendicular of first bond
            nbrs = list(a.GetNeighbors())
            if nbrs:
                j = nbrs[0].GetIdx(); vx, vy = coords[j][0] - coords[i][0], coords[j][1] - coords[i][1]
                nx, ny = -vy, vx
            else:
                nx, ny = 0, -1
        else:
            nx, ny = -nx, -ny
        L = math.hypot(nx, ny) or 1
        dirs.append([round(nx / L, 3), round(ny / L, 3)])
    return {"svg": svg, "coords": coords, "dirs": dirs, "width": width, "height": height, "molblock": Chem.MolToMolBlock(m)}


def depict_svg_h(mol: Chem.Mol, width=680, height=500) -> dict:
    """Same as depict_svg but with all hydrogens explicit; heavy-atom indices are preserved
    (AddHs appends H atoms), and h_parent gives the heavy atom each H is bonded to."""
    mh = Chem.AddHs(mol)
    d = depict_svg(mh, width, height, font=15, bond=46)
    d["h_parent"] = [(a.GetNeighbors()[0].GetIdx() if a.GetAtomicNum() == 1 and a.GetDegree() else -1) for a in mh.GetAtoms()]
    d["is_h"] = [a.GetAtomicNum() == 1 for a in mh.GetAtoms()]
    return d


def geometry3d(mol: Chem.Mol, max_atoms=150) -> dict:
    """
    ETKDGv3 embedding + MMFF94s minimisation.  Returns molblock (with H) plus bond lengths
    (Å) and bond angles (°) for the heavy-atom skeleton, ready for tabulation.
    """
    if mol.GetNumAtoms() > max_atoms:
        return {"error": f"molecule too large for 3D ({mol.GetNumAtoms()} heavy atoms)"}
    mh = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3(); ps.randomSeed = 0xF00D; ps.useSmallRingTorsions = True
    cid = AllChem.EmbedMolecule(mh, ps)
    if cid < 0:
        ps.useRandomCoords = True
        cid = AllChem.EmbedMolecule(mh, ps)
    if cid < 0:
        return {"error": "3D embedding failed"}
    ff_name = "MMFF94s"
    try:
        props = AllChem.MMFFGetMoleculeProperties(mh, mmffVariant="MMFF94s")
        ff = AllChem.MMFFGetMoleculeForceField(mh, props)
        ff.Minimize(maxIts=2000)
        energy = ff.CalcEnergy()
    except Exception:
        AllChem.UFFOptimizeMolecule(mh, maxIters=2000); ff_name = "UFF"; energy = None
    conf = mh.GetConformer()

    def lbl(a):
        return f"{a.GetSymbol()}{a.GetIdx() + 1}"

    lengths = []
    for b in mh.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue
        d = (conf.GetAtomPosition(a1.GetIdx()) - conf.GetAtomPosition(a2.GetIdx())).Length()
        lengths.append({"a": a1.GetIdx(), "b": a2.GetIdx(), "label": f"{lbl(a1)}–{lbl(a2)}",
                        "order": b.GetBondTypeAsDouble(), "aromatic": b.GetIsAromatic(), "length": round(d, 3)})
    angles = []
    from rdkit.Chem import rdMolTransforms as T
    for a in mh.GetAtoms():
        if a.GetAtomicNum() == 1:
            continue
        nbrs = [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() != 1]
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                ang = T.GetAngleDeg(conf, nbrs[i], a.GetIdx(), nbrs[j])
                angles.append({"a": nbrs[i], "b": a.GetIdx(), "c": nbrs[j],
                               "label": f"{lbl(mh.GetAtomWithIdx(nbrs[i]))}–{lbl(a)}–{lbl(mh.GetAtomWithIdx(nbrs[j]))}",
                               "angle": round(ang, 1)})
    # C–H lengths summarised
    ch = [(conf.GetAtomPosition(b.GetBeginAtomIdx()) - conf.GetAtomPosition(b.GetEndAtomIdx())).Length()
          for b in mh.GetBonds() if {b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()} == {1, 6}]
    return {"molblock": Chem.MolToMolBlock(mh), "forcefield": ff_name, "symbols": [a.GetSymbol() for a in mol.GetAtoms()],
            "energy_kcal": round(energy, 2) if energy is not None else None,
            "bond_lengths": lengths, "bond_angles": angles,
            "ch_mean": round(sum(ch) / len(ch), 3) if ch else None,
            "note": "Gas-phase force-field geometry (ETKDGv3 + MMFF94s); typical accuracy ±0.02 Å / ±2° vs crystal structures."}
