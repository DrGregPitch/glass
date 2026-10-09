// 2D depiction overlay (IUPAC locants, hover → name highlighting), proton/carbon mirrors, 3D viewer.
// Load order matters only for the top-level statements; all files share the page's global scope.

// ---------------- structure
// Crop an RDKit SVG (fixed canvas) to the molecule's bounding box so drawings fill their frames.
function bbox(coords, pad) { let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9; for (const [x, y] of coords) { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); } return [x0 - pad, y0 - pad, Math.max(x1 - x0 + 2 * pad, 120), Math.max(y1 - y0 + 2 * pad, 90)]; }
function cropped(dep, extraInner, pad, maxH) {
  const [x, y, w, h] = bbox(dep.coords, pad);
  const inner = dep.svg.replace(/^[\s\S]*?<svg[^>]*>/, '').replace(/<\/svg>\s*$/, '');
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="${x} ${y} ${w} ${h}" style="width:100%;height:auto;max-height:${maxH}px;display:block">${inner}${extraInner || ''}</svg>`;
}
function renderStructure(d) {
  DEP = d.depiction; DEPH = d.depiction_h || null; GEOM = d.geometry || null;
  $('#thumb').innerHTML = cropped(DEP, '', 34, 300);
  const [bx, by, bw, bh] = bbox(DEP.coords, 40);
  $('#struct').innerHTML = cropped(DEP, '', 40, 460) + `<svg id="overlay" viewBox="${bx} ${by} ${bw} ${bh}" style="width:100%;height:auto;max-height:460px;position:absolute;left:0;top:0"></svg>`;
  drawOverlay();
  if (d.geometry && !d.geometry.error) render3D(d.geometry); else $('#geomnote').textContent = d.geometry ? d.geometry.error : '';
}
function drawOverlay() {
  if (!DEP) return;
  const ov = $('#overlay'); const showP = $('#t_loc').checked, showS = showP, showH = showP, showI = $('#t_idx').checked;   // one switch shows every locant kind
  let s = '';
  DEP.coords.forEach(([x, y], i) => {
    const loc = SEL && SEL.locants ? SEL.locants[i] : null; const frag = SEL && SEL.fragments ? SEL.fragments[i] : null;
    s += `<circle class="hit${(i === hoverAtom || nameHover.has(i)) ? ' on' : ''}" cx="${x}" cy="${y}" r="11" data-i="${i}"></circle>`;
    const [dx, dy] = DEP.dirs[i]; const lx = x + dx * 14, ly = y + dy * 14 + 4;
    if (loc && ((frag === 'parent' && showP) || (frag === 'sub' && showS) || (frag === 'hetero' && showH))) s += `<text class="loc ${frag}" x="${lx}" y="${ly}" text-anchor="middle">${esc(frag === 'sub' ? loc + '′' : loc)}</text>`;
    if (showI) s += `<text class="loc idx" x="${x - dx * 13}" y="${y - dy * 13 + 3}" text-anchor="middle">${i}</text>`;
  });
  ov.innerHTML = s;
  mirrorStructure(); mirrorH();
  ov.querySelectorAll('.hit').forEach(c => {
    c.onmouseenter = () => { const i = +c.dataset.i; hoverAtom = i; hoverLoc = SEL && SEL.locants ? SEL.locants[i] : null; hoverInfo(i); renderNameBox(); updateHover(); };
    c.onmouseleave = () => { hoverLoc = null; hoverAtom = -1; renderNameBox(); updateHover(); $('#hoverinfo').textContent = HOVER_HINT; };
  });
}
function mirrorStructure() {
  const m = $('#struct2'); if (!m || !DEP) return;
  const base = $('#struct').querySelector('svg:not(#overlay)'); if (!base) return;
  m.innerHTML = cropped(DEP, $('#overlay').innerHTML, 40, 420);
}
function hLabel(parent) {
  if (!SEL || !hasNum(SEL) || !SEL.locants) return '';
  const loc = SEL.locants[parent], f = SEL.fragments[parent]; if (!loc) return '';
  return f === 'sub' ? `H-${loc}′` : f === 'hetero' ? `${loc}H` : `H-${loc}`;
}
function mirrorH() {
  const m = $('#structH'); if (!m || !DEPH) return;
  let s = ''; const labelled = new Set();
  DEPH.coords.forEach(([x, y], i) => {
    if (!DEPH.is_h[i]) return;
    s += `<circle class="hit" cx="${x}" cy="${y}" r="11" data-i="${i}"></circle>`;
    const par = DEPH.h_parent[i]; if (labelled.has(par)) return; labelled.add(par);
    const lab = hLabel(par); if (!lab) return;
    const [dx, dy] = DEPH.dirs[i];
    s += `<text class="loc hlab" x="${x + dx * 15}" y="${y + dy * 15 + 4}" text-anchor="middle">${esc(lab)}</text>`;
  });
  m.innerHTML = cropped(DEPH, s, 34, 460);
}
function hoverInfo(i) {
  const loc = SEL && SEL.locants ? SEL.locants[i] : null; const frag = SEL && SEL.fragments ? SEL.fragments[i] : null;
  $('#hoverinfo').textContent = `atom ${i}` + (loc ? ` · locant ${loc} (${frag})` : ' · no locant in this name');
}
['t_loc', 't_idx'].forEach(id => $('#' + id).onchange = drawOverlay);
// Hover only flips classes on the circles drawOverlay() already drew: no innerHTML rebuild, no
// handler re-attachment, and the circle under the pointer is never destroyed mid-hover.
function updateHover() {
  $$('#overlay .hit').forEach(c => c.classList.toggle('on', +c.dataset.i === hoverAtom || nameHover.has(+c.dataset.i)));
  if (typeof mirrorStructure === 'function') mirrorStructure();
}
function structureSVG(withNumbers) {
  if (withNumbers === undefined) withNumbers = $('#expnum').checked;
  const [x, y, w, h] = bbox(DEP.coords, 40);
  let overlay = $('#overlay').innerHTML.replace(/<circle[^>]*><\/circle>/g, '');
  if (!withNumbers) overlay = '';
  return `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="${x} ${y} ${w} ${h}"><style>text.loc{font:bold 13px Arial,sans-serif}text.loc.parent{fill:#4f46e5}text.loc.sub{fill:#b45309}text.loc.hetero{fill:#0f766e}text.loc.idx{fill:#888;font-weight:400}</style>` + $('#struct').querySelector('svg').innerHTML + overlay + '</svg>';
}
function exportPNG() {
  if (!DEP) return; const svg = structureSVG(); const img = new Image(); const [x, y, w, h] = bbox(DEP.coords, 40);
  img.onload = () => { const c = document.createElement('canvas'); c.width = w * 3; c.height = h * 3; const g = c.getContext('2d'); g.fillStyle = '#fff'; g.fillRect(0, 0, c.width, c.height); g.drawImage(img, 0, 0, c.width, c.height); const a = document.createElement('a'); a.href = c.toDataURL('image/png'); a.download = ((R && R.common_name) || 'compound').replace(/[^\w\-]+/g, '_') + '.png'; a.click(); };
  img.src = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(svg)));
}
async function exportCDXML() {
  if (!R) return; const withN = $('#expnum').checked && SEL && hasNum(SEL);
  const u = new URL('/api/cdxml', location.origin);
  u.searchParams.set('smiles', R.smiles_canonical); u.searchParams.set('numbers', withN);
  if (withN) u.searchParams.set('name', SEL.name);
  const txt = await (await fetch(u)).text();
  download(((R.common_name || R.formula) + '').replace(/[^\w\-]+/g, '_') + '.cdxml', txt, 'chemical/x-cdxml');
  toast(withN ? 'CDXML with numbering annotations' : 'CDXML saved');
}
$('#dlsvg').onclick = () => download((SEL ? SEL.name : 'structure').replace(/[^\w\-]+/g, '_') + '.svg', structureSVG(), 'image/svg+xml');
$('#dlpng2').onclick = exportPNG;
$('#dlcdxml').onclick = exportCDXML;

let viewer = null, EL = null, SURF = null, CUBES = {};
// deep-saturation diverging colormap (blue = +, red = −, cream at 0), 9 anchor stops
const DIV_STOPS = [[103,0,31],[178,24,43],[214,96,77],[244,165,130],[247,247,245],[146,197,222],[67,147,195],[33,102,172],[5,48,97]];
function divColor(t) {   // t in [0,1]: 0 = most negative (deep red) … 1 = most positive (deep blue)
  const x = t * (DIV_STOPS.length - 1), i = Math.max(0, Math.min(DIV_STOPS.length - 2, Math.floor(x))), f = x - i;
  const c = DIV_STOPS[i].map((v, k) => Math.round(v + (DIV_STOPS[i + 1][k] - v) * f));
  return '#' + c.map(v => v.toString(16).padStart(2, '0')).join('');
}
async function loadElectronic() {
  EL = null; $('#modiag').innerHTML = ''; $('#mokv').innerHTML = '<span class="hint"><span class="spin"></span> GFN2-xTB…</span>';
  try { EL = await api('/api/electronic', { smiles: R.smiles_canonical }); } catch (e) { $('#mokv').innerHTML = `<span class="hint">${esc(e.message)}</span>`; return; }
  if (EL.error) { $('#mokv').innerHTML = `<span class="hint">${esc(EL.error)}</span>`; return; }
  buildMOOptions(); applySurface();
  // pre-warm the two most requested orbitals in the background (SCF is shared server-side)
  for (const key of [['homo', 0], ['lumo', 0]]) {
    const k = key[0] + key[1];
    if (!CUBES[k]) fetch(`/api/orbital?smiles=${encodeURIComponent(R.smiles_canonical)}&which=${key[0]}&offset=${key[1]}`).then(r => r.ok ? r.json() : null).then(j => { if (j) CUBES[k] = j; }).catch(() => {});
  }
}
function buildMOOptions() {
  const g = $('#mogroup'); if (!g || !EL || !EL.levels) return;
  const occ = EL.levels.filter(l => l.occ > 0.5).length;
  const virt = EL.levels.length - occ;
  let html = '';
  for (let k = Math.min(3, occ - 1); k >= 1; k--) html += `<option value="mo:homo:${k}">HOMO−${k}</option>`;
  html += '<option value="mo:homo:0">HOMO</option><option value="mo:lumo:0">LUMO</option>';
  for (let k = 1; k <= Math.min(3, virt - 1); k++) html += `<option value="mo:lumo:${k}">LUMO+${k}</option>`;
  const cur = $('#surface').value; g.innerHTML = html;
  if ([...$('#surface').options].some(o => o.value === cur)) $('#surface').value = cur;
}
$('#mobtn').onclick = () => { $('#moModal').classList.add('open'); if (EL && !EL.error) drawMODiagram(); };
function drawMODiagram() {
  const svg = $('#modiag'); const L = EL.levels; if (!L || !L.length) return;
  const OCC = '#17171a', VIRT = '#9a968a', GAPC = '#3d3d42';
  // degenerate grouping, then RANK-BASED vertical layout: constant pitch, no cramming
  const groups = []; L.forEach(l => { const g = groups.find(g => Math.abs(g.e - l.e_ev) < 0.06); if (g) g.ls.push(l); else groups.push({ e: l.e_ev, ls: [l] }); });
  groups.sort((p, q) => p.e - q.e);
  const occG = groups.filter(g => g.ls[0].occ > 0.5), virtG = groups.filter(g => g.ls[0].occ <= 0.5);
  const PITCH = 30, GAPBAND = 74, T = 34, B = 26;
  const W = 360, X0 = 96, X1 = 236, CHIPX = 246;
  const H = T + virtG.length * PITCH + GAPBAND + occG.length * PITCH + B;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const yOf = g => {
    if (g.ls[0].occ > 0.5) { const i = occG.indexOf(g); return T + virtG.length * PITCH + GAPBAND + (occG.length - 1 - i) * PITCH + PITCH / 2; }
    const i = virtG.indexOf(g); return T + (virtG.length - 1 - i) * PITCH + PITCH / 2;
  };
  const yHomoTop = T + virtG.length * PITCH + GAPBAND, yVirtBot = T + virtG.length * PITCH;
  let out = `<defs><linearGradient id="occg" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="${OCC}" stop-opacity="0.12"/><stop offset="1" stop-color="${OCC}" stop-opacity="0.03"/></linearGradient></defs><rect x="0" y="0" width="${W}" height="${H}" rx="14" fill="#fbfaf6"/>`;
  out += `<rect x="4" y="${yHomoTop - 6}" width="${W - 8}" height="${occG.length * PITCH + 10}" rx="10" fill="url(#occg)"/>`;
  out += `<text x="14" y="${yHomoTop + occG.length * PITCH + 0}" font-size="9.5" font-weight="800" fill="${OCC}" opacity="0.7" letter-spacing="2">FILLED</text>`;
  out += `<rect x="4" y="${T - 8}" width="${W - 8}" height="${virtG.length * PITCH + 10}" rx="10" fill="${VIRT}" opacity="0.04"/>`;
  out += `<text x="14" y="${T + 2}" font-size="9.5" font-weight="700" fill="${VIRT}" opacity="0.65" letter-spacing="2">EMPTY (VIRTUAL)</text>`;
  const arrows = (x, y, col, paired) => {
    let a = `<g stroke="${col}" stroke-width="1.6" stroke-linecap="round" fill="none">`;
    const up = dx => `<line x1="${x + dx}" y1="${y + 6}" x2="${x + dx}" y2="${y - 6}"/><path d="M${x + dx - 2.8} ${y - 2.8} L${x + dx} ${y - 6.3} L${x + dx + 2.8} ${y - 2.8}"/>`;
    const dn = dx => `<line x1="${x + dx}" y1="${y - 6}" x2="${x + dx}" y2="${y + 6}"/><path d="M${x + dx - 2.8} ${y + 2.8} L${x + dx} ${y + 6.3} L${x + dx + 2.8} ${y + 2.8}"/>`;
    return a + (paired ? up(-3.6) + dn(3.6) : up(0)) + '</g>';
  };
  for (const g of groups) {
    const y = yOf(g); const n = g.ls.length; const w = (X1 - X0 - (n - 1) * 12) / n;
    const isF = g.ls.some(l => l.label === 'HOMO' || l.label === 'LUMO');
    g.ls.forEach((l, k) => {
      const x = X0 + k * (w + 12); const occ = l.occ > 0.5; const col = occ ? OCC : VIRT;
      out += `<line x1="${x}" y1="${y}" x2="${x + w}" y2="${y}" stroke="${col}" stroke-width="${occ ? (isF ? 5 : 3.6) : (isF ? 3 : 2)}" stroke-linecap="round" ${occ ? '' : 'opacity="0.65"'}/>`;
      if (occ) out += arrows(x + w / 2, y, OCC, l.occ > 1.5);
    });
    const name = g.ls.length > 1 ? `${g.ls[0].label} ×${g.ls.length}` : g.ls[0].label;
    out += `<text x="${X0 - 10}" y="${y + 3}" font-size="9" text-anchor="end" fill="${isF ? (g.ls[0].occ > 0.5 ? OCC : VIRT) : '#9a968a'}" font-weight="${isF ? 800 : 400}">${name}</text>`;
    out += `<text x="${CHIPX}" y="${y + 3}" font-size="9.5" font-family="ui-monospace,Menlo,monospace" fill="${isF ? (g.ls[0].occ > 0.5 ? OCC : VIRT) : '#8a8a93'}" font-weight="${isF ? 700 : 400}">${g.e.toFixed(2)} eV</text>`;
  }
  // the gap band
  const ym = (yVirtBot + yHomoTop) / 2, xm = (X0 + X1) / 2;
  out += `<g stroke="${GAPC}" stroke-width="2.2" fill="${GAPC}" stroke-linecap="round">
    <line x1="${xm}" y1="${yVirtBot + 8}" x2="${xm}" y2="${ym - 14}"/><line x1="${xm}" y1="${ym + 14}" x2="${xm}" y2="${yHomoTop - 8}"/>
    <path d="M${xm - 4.5} ${yVirtBot + 13} L${xm} ${yVirtBot + 4} L${xm + 4.5} ${yVirtBot + 13} Z"/><path d="M${xm - 4.5} ${yHomoTop - 13} L${xm} ${yHomoTop - 4} L${xm + 4.5} ${yHomoTop - 13} Z"/></g>`;
  out += `<rect x="${xm - 56}" y="${ym - 13}" width="112" height="26" rx="13" fill="#fff" stroke="${GAPC}" stroke-width="1.6"/>
    <text x="${xm}" y="${ym + 4}" font-size="12" font-weight="800" text-anchor="middle" fill="${GAPC}">gap Δ ${EL.gap_ev} eV</text>`;
  out += `<text x="${W / 2}" y="${H - 8}" font-size="8.5" text-anchor="middle" fill="#9a968a">rank layout; printed energies are the data (GFN2-xTB)</text>`;
  svg.innerHTML = out;
  $('#mokv').innerHTML = `<b>HOMO</b><span>${EL.homo_ev} eV (highest filled orbital)</span><b>LUMO</b><span>${EL.lumo_ev} eV (lowest empty orbital)</span><b>Gap</b><span><b style="color:${GAPC}">${EL.gap_ev} eV</b> <span class="hint">(xTB gaps run 2–3 eV small)</span></span><b>Dipole</b><span>${EL.dipole_debye} D</span><b>Filling</b><span>${EL.n_electrons} e⁻ in ${EL.n_orbitals} orbitals</span><b>View in 3D</b><span class="hint">pick any of these orbitals in the Surface selector to see its isosurface and charge-centroid arrow</span>`;
}
async function applySurface() {
  const mode = $('#surface').value; SURF = mode;
  if (!viewer || !GEOM) return;
  viewer.removeAllSurfaces(); viewer.removeAllShapes(); viewer.clear();
  const isMO = mode.startsWith('mo:');
  const useEL = EL && !EL.error && (mode === 'vdw-q' || mode === 'sas-q');
  let mb = GEOM.molblock, moData = null;
  $('#surfnote').textContent = ''; $('#surflegend').style.display = 'none';
  if (isMO) {
    const [, which, off] = mode.split(':'); const key = which + off;
    $('#surfnote').innerHTML = '<span class="spin"></span> computing orbital…';
    try {
      if (!CUBES[key]) { const r = await fetch(`/api/orbital?smiles=${encodeURIComponent(R.smiles_canonical)}&which=${which}&offset=${off}`); if (!r.ok) throw new Error((await r.json()).detail); CUBES[key] = await r.json(); }
      moData = CUBES[key]; mb = moData.molblock;
    } catch (e) { $('#surfnote').textContent = 'orbital unavailable: ' + e.message; moData = null; }
  } else if (useEL) mb = EL.molblock;
  const model = viewer.addModel(mb, 'mol');
  viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.22 } });
  if (useEL) {
    const qs = EL.charges; const qmax = Math.max(0.08, ...qs.map(Math.abs));
    model.selectedAtoms({}).forEach((a, i) => { a.properties = a.properties || {}; a.properties.q = qs[i] ?? 0; });
    const st = mode === 'vdw-q' ? $3Dmol.SurfaceType.VDW : $3Dmol.SurfaceType.SAS;
    viewer.addSurface(st, { opacity: 0.96, colorfunc: a => divColor(((a.properties && a.properties.q || 0) + qmax) / (2 * qmax)) });
    $('#surfnote').textContent = `GFN2-xTB partial charge (approximate electrostatics); range ±${qmax.toFixed(2)} e`;
    $('#surflegend').style.display = 'inline-flex';
    $('#surfbar').style.background = `linear-gradient(90deg,${[0, .125, .25, .375, .5, .625, .75, .875, 1].map(divColor).join(',')})`;
    $('#surfmin').textContent = `−${qmax.toFixed(2)} e`; $('#surfmax').textContent = `+${qmax.toFixed(2)} e`;
  } else if (moData) {
    viewer.addVolumetricData(moData.cube, 'cube', { isoval: 0.02, color: '#2166ac', opacity: 0.85, smoothness: 8 });
    viewer.addVolumetricData(moData.cube, 'cube', { isoval: -0.02, color: '#b2182b', opacity: 0.85, smoothness: 8 });
    $('#surfnote').textContent = `${moData.label}: ε = ${moData.energy_ev} eV (minimal-basis; nodal structure qualitative)` + (moData.occ ? '' : ' · virtual orbital (unoccupied: no charge, dashed arrow shows its charge centroid)');
    if ($('#showdip').checked) {
      const c = moData.centroid, n = moData.nuc_center;
      const v = [c[0] - n[0], c[1] - n[1], c[2] - n[2]]; const len = Math.hypot(...v);
      if (len > 0.03) {
        const sc = Math.min(4, 1 + 2 * len) / len;
        viewer.addArrow({ start: { x: n[0], y: n[1], z: n[2] }, end: { x: n[0] + v[0] * sc, y: n[1] + v[1] * sc, z: n[2] + v[2] * sc }, radius: 0.09, color: moData.occ ? '#2e2e33' : '#9a968a' });
        const dl = Math.hypot(...moData.orb_dipole_debye);
        $('#surfnote').textContent += moData.occ ? ` · arrow: this orbital's charge centroid (its dipole contribution ${dl.toFixed(1)} D, pointing toward electron density)` : '';
      } else $('#surfnote').textContent += ' · charge centroid coincides with the nuclear centre (symmetric orbital, no dipole contribution)';
    }
  }
  if ($('#showdip').checked && !isMO && EL && !EL.error) {
    const mu = EL.dipole_vec_debye, n = EL.nuc_center; const len = Math.hypot(...mu);
    if (len > 0.05) {
      const sc = Math.min(4, 1 + 0.8 * len) / len;
      viewer.addArrow({ start: { x: n[0], y: n[1], z: n[2] }, end: { x: n[0] + mu[0] * sc, y: n[1] + mu[1] * sc, z: n[2] + mu[2] * sc }, radius: 0.1, color: '#2e2e33' });
      $('#surfnote').textContent += ` · total dipole μ = ${len.toFixed(2)} D (physics convention: arrow points − → +)`;
    } else $('#surfnote').textContent += ' · dipole ≈ 0 by symmetry';
  }
  viewer.zoomTo(); viewer.render();
}
$('#surface').onchange = applySurface;
$('#showdip').onchange = applySurface;
function render3D(g) {
  $('#geomnote').textContent = `${g.forcefield}${g.energy_kcal != null ? ', E = ' + g.energy_kcal + ' kcal/mol' : ''}. ${g.note}`;
  try {
    if (!viewer) viewer = $3Dmol.createViewer($('#viewer3d'), { backgroundColor: 'white' });
    viewer.clear(); viewer.addModel(g.molblock, 'mol'); viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.22 } }); viewer.zoomTo(); viewer.render();
  } catch (e) { $('#viewer3d').innerHTML = '<div class="hint" style="padding:10px">3D viewer unavailable offline (3Dmol.js CDN)</div>'; }
  renderGeomTables();
  loadElectronic();
}
function geomAtomLabel(i) {   // element + locant of the SELECTED name (primed for substituents); falls back to element+index
  const mol = DEP; const sym = GEOM && GEOM.symbols ? GEOM.symbols[i] : null;
  const loc = SEL && SEL.locants && SEL.locants[i]; const f = SEL && SEL.fragments && SEL.fragments[i];
  const base = sym || (mol && mol.symbols ? mol.symbols[i] : 'X');
  if (loc && hasNum(SEL) && f !== 'hetero') return base + loc + (f === 'sub' ? '′' : '');
  if (loc && f === 'hetero') return base + '(' + loc + ')';
  return base + (i + 1);
}
function renderGeomTables() {
  if (!GEOM || GEOM.error) return;
  const lbl = i => geomAtomLabel(i);
  $('#bl').innerHTML = '<tr><th>bond</th><th>order</th><th>length / Å</th></tr>' + GEOM.bond_lengths.map(b => `<tr><td class="mono">${lbl(b.a)}–${lbl(b.b)}</td><td>${b.aromatic ? 'arom.' : b.order}</td><td>${b.length.toFixed(3)}</td></tr>`).join('');
  $('#ba').innerHTML = '<tr><th>angle</th><th>°</th></tr>' + GEOM.bond_angles.map(b => `<tr><td class="mono">${lbl(b.a)}–${lbl(b.b)}–${lbl(b.c)}</td><td>${b.angle.toFixed(1)}</td></tr>`).join('');
}
