// Identifier verification, batch mode, exports (Markdown/JSON/CSV/PNG/SVG), methods modal, structure editor (JSME).
// Load order matters only for the top-level statements; all files share the page's global scope.

// ---------------- verify box
$('#vgo').onclick = async () => {
  const name = $('#vname').value.trim(); if (!name || !R) return;
  $('#vout').innerHTML = '<div class="hint" style="margin-top:8px"><span class="spin"></span> checking…</div>';
  const v = await api('/api/verify', { name, inchikey: R.inchikey });
  const head = { identical: '<span class="badge ok">same compound</span>', same_connectivity: '<span class="badge warn">same skeleton, different stereo / charge</span>', different: '<span class="badge bad">different compound</span>', uninterpretable: '<span class="badge info">not interpreted</span>' }[v.status];
  $('#vout').innerHTML = `<div class="verdict ${v.status}">${head} ${esc(v.note)}${v.kind ? `<div class="hint" style="margin-top:4px">read as ${esc(v.kind)} → <span class="mono">${esc(v.smiles)}</span> (${esc(v.formula)}, ${esc(v.inchikey)})</div>` : (v.detail ? `<div class="hint" style="margin-top:4px">${esc(v.detail)}</div>` : '')}</div>`;
};
$('#vname').addEventListener('keydown', e => { if (e.key === 'Enter') $('#vgo').click(); });


// ---------------- export
function reportMarkdown() {
  const L = [];
  L.push(`# ${R.common_name || R.iupac_name || R.formula}`, '', `_Glass report, ${new Date().toISOString().slice(0, 10)}_`, '');
  L.push('## Identity', '', '| | |', '|---|---|', `| Input | \`${R.input}\` (${R.input_kind}) |`, `| Canonical SMILES | \`${R.smiles_canonical}\` |`, `| IUPAC name | ${R.iupac_name || '·'} |`, `| Common name | ${R.common_name || '·'} |`, `| CAS RN | ${R.cas || '·'} |`, `| Formula | ${R.formula} |`, `| InChI | \`${R.inchi}\` |`, `| InChIKey | \`${R.inchikey}\` |`, R.cid ? `| PubChem CID | ${R.cid} |` : '', '');
  if (R.warnings.length) L.push(...R.warnings.map(w => `> ⚠ ${w}`), '');
  L.push('## Names', '', '| status | name | kind |', '|---|---|---|', ...R.names.map(n => `| ${n.verified ? 'verified' : 'unverified'} | ${n.name} | ${n.kind} |`), '');
  if (PROPS) L.push('## Predicted properties', '', '| property | value |', '|---|---|', ...propRows(PROPS).map(([k, v]) => `| ${k} | ${v} |`), '');
  if (R.experimental && Object.keys(R.experimental).length) L.push('## Experimental (PubChem)', '', ...Object.entries(R.experimental).map(([k, v]) => `- **${k}**: ${v.join('; ')}`), '');
  if (SPEC) {
    L.push(`## ¹H NMR (simulated, ${SPEC.h1.solvent}, ${SPEC.h1.mhz} MHz)`, '', '| δ / ppm | mult. | J / Hz | nH | assignment |', '|---|---|---|---|---|', ...SPEC.h1.peaks.map(p => `| ${p.shift.toFixed(2)} | ${p.multiplicity} | ${(p.J || []).join(', ')} | ${p.nH}H | ${p.type} |`), '');
    L.push(`## ¹³C NMR (simulated, ${SPEC.c13.solvent})`, '', '| δ / ppm | DEPT | type |', '|---|---|---|', ...SPEC.c13.peaks.map(p => `| ${p.shift.toFixed(1)} | ${p.dept} | ${p.type} |`), '');
    L.push('## IR (group frequencies)', '', ...SPEC.ir.bands.map(b => `- ${b.center} cm⁻¹: ${b.label}`), '');
    L.push('## UV–Vis', '', ...SPEC.uv.peaks.map(p => `- ${p.lambda_max ?? '·'} nm (ε ≈ ${p.epsilon ?? '·'}): ${p.assignment}`), '');
    if (PL && 'S' in PL) L.push('## Photoluminescence (estimate)', '', `- class: ${PL.class}`, `- absorption max ${PL.lambda_abs_max} nm, 0-0 ${PL.lambda_00} nm, emission max ${PL.lambda_em_max} nm (${PL.colour}), Stokes shift ${PL.stokes_shift_cm} cm⁻¹`, `- Huang–Rhys S = ${PL.S}, effective mode ${PL.omega} cm⁻¹`, ...PL.reasons.map(r => `- ${r}`), '');
  }
  L.push('---', 'Predictions are estimates (see stated errors). Names marked verified were parsed back to a structure with an identical Standard InChIKey.');
  return L.filter(x => x !== undefined).join('\n');
}


// ---------------- methods & limitations
let VALID = null;
async function showMethods(keys) {
  const ks = keys && keys.length ? keys : Object.keys(METHODS);
  if (!VALID) { try { VALID = await api('/api/validation', {}); } catch { VALID = {}; } }
  const extra = k => {
    if (k === 'nmr' && VALID.nmr_nmrshiftdb2 && !VALID.nmr_nmrshiftdb2.error) { const v = VALID.nmr_nmrshiftdb2; const c = v['13C'] && v['13C'].MAE, h = v['1H'] && v['1H'].MAE; return `<div class="row" style="margin-top:6px"><b>Validation</b><span>Database-level predictions are validated on held-out experimental spectra${c != null ? ` (¹³C MAE ${c} ppm` + (h != null ? `, ¹H ${h} ppm)` : ')') : ''}; each peak also carries its own match quality and uncertainty.</span></div>`; }
    if (k === 'ir_modes' && VALID.ir) { const b = VALID.ir; return `<div class="row" style="margin-top:6px"><b>Validation</b><span>Vibrational frequencies are benchmarked against ${b.rows ? b.rows.length : ''} characteristic experimental bands; the higher level lands within roughly ±30 cm⁻¹.</span></div>`; }
    return '';
  };
  $('#mBody').innerHTML = ks.map(k => { const m = METHODS[k]; if (!m) return ''; return `<div class="method"><h4>${esc(m.title)}</h4><div class="row"><b>Level of theory</b><span>${esc(m.theory)}</span><b>Good for</b><span>${esc(m.good_for)}</span><b>Cannot resolve</b><span>${esc(m.fails_at)}</span><b>Typical error</b><span>${esc(m.error)}</span></div>${extra(k)}<div class="refs">${m.refs.map(esc).join(' · ')}</div></div>`; }).join('');
  $('#methodsModal').classList.add('open');
}
document.addEventListener('click', e => { const b = e.target.closest('.ib'); if (b) { e.preventDefault(); showMethods(b.dataset.m.split(',')); } });
$('#methodsLink').onclick = e => { e.preventDefault(); showMethods(null); };
document.addEventListener('keydown', e => { if (e.key === 'Escape') $$('.modal.open').forEach(m => m.classList.remove('open')); });
$$('.modal').forEach(m => m.addEventListener('click', e => { if (e.target === m) m.classList.remove('open'); }));


// ---------------- downloads
function canvasPNG(cv, name) { const a = document.createElement('a'); const c = document.createElement('canvas'); c.width = cv.width; c.height = cv.height; const g = c.getContext('2d'); g.fillStyle = '#fff'; g.fillRect(0, 0, c.width, c.height); g.drawImage(cv, 0, 0); a.href = c.toDataURL('image/png'); a.download = name; a.click(); }
function csvOf(kind) {
  const q = v => '"' + String(v ?? '').replace(/"/g, '""') + '"'; let rows = [];
  if (kind === 'h1' && SPEC) { rows = [['delta_ppm', 'multiplicity', 'J_Hz', 'nH', 'assignment', 'protons'], ...SPEC.h1.peaks.map(p => [p.shift, p.multiplicity, (p.J || []).join(' '), p.nH, p.type, p.atoms.map(a => hLabel(a) || '#' + a).join(' ')]), [], ['x_ppm', 'y'], ...SPEC.h1.x.map((x, i) => [x, SPEC.h1.y[i]])]; }
  if (kind === 'c13' && SPEC) { rows = [['delta_ppm', 'DEPT', 'type', 'atoms'], ...SPEC.c13.peaks.map(p => [p.shift, p.dept, p.type, p.atoms.map(a => atomLabel(a)).join(' ')]), [], ['x_ppm', 'y'], ...SPEC.c13.x.map((x, i) => [x, SPEC.c13.y[i]])]; }
  if (kind === 'ir' && SPEC) { rows = [['wavenumber_cm-1', 'intensity', 'assignment', 'source'], ...SPEC.ir.bands.map(b => [b.center, b.intensity, b.label, 'group frequency']), ...(MODES && MODES.modes ? MODES.modes.map(m => [m.freq, m.intensity, m.label, 'MMFF94s normal mode']) : []), [], ['x_cm-1', 'percent_T'], ...SPEC.ir.x.map((x, i) => [x, SPEC.ir.y[i]])]; }
  if (kind === 'uv' && SPEC) { rows = [['lambda_max_nm', 'epsilon', 'assignment'], ...SPEC.uv.peaks.map(p => [p.lambda_max, p.epsilon, p.assignment]), [], ['x_nm', 'epsilon'], ...SPEC.uv.x.map((x, i) => [x, SPEC.uv.y[i]])]; }
  if (kind === 'pl' && PL && 'S' in PL) { rows = [['band', 'transition', 'lambda_nm', 'wavenumber_cm-1', 'FC_factor'], ...PL.absorption.map(p => ['absorption', p.label, p.nm, p.nu, p.fc]), ...PL.emission.map(p => ['emission', p.label, p.nm, p.nu, p.fc]), [], ['x_nm', 'absorption', 'emission'], ...PL.x.map((x, i) => [x, PL.y_abs[i], PL.y_em[i]])]; }
  return rows.map(r => r.map(q).join(',')).join('\n');
}
document.addEventListener('click', e => {
  const b = e.target.closest('[data-dl]'); if (!b) return;
  const [kind, id, name] = b.dataset.dl.split(':'); const base = (R && (R.common_name || R.formula) || 'chemxref').replace(/[^\w\-]+/g, '_');
  if (kind === 'png') canvasPNG($('#' + id), `${base}_${name}.png`);
  if (kind === 'csv') download(`${base}_${id}.csv`, csvOf(id), 'text/csv');
});
$('#dlFC').onclick = () => { const svg = $('#fc'); download(((R && R.common_name) || 'compound').replace(/[^\w\-]+/g, '_') + '_franck_condon.svg', '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 520 470" width="1040" height="940" style="font-family:Inter,Arial,sans-serif;background:#fff"><style>.lvl{stroke:#9aa0a6;stroke-width:1}.arrow{stroke-width:1.6;opacity:.3}.arrow.on{opacity:1;stroke-width:3}.arrow.base{opacity:.6;stroke-width:2}.abs{stroke:#4f46e5}.em{stroke:#0f766e}.lvl.on{stroke:#f59e0b;stroke-width:2.5}.wf{fill:none;stroke-width:1.2;opacity:.6}</style>' + svg.innerHTML + '</svg>', 'image/svg+xml'); };
$('#dl3d').onclick = () => { if (!viewer) return; const a = document.createElement('a'); a.href = viewer.pngURI(); a.download = ((R && R.common_name) || 'compound').replace(/[^\w\-]+/g, '_') + '_3D.png'; a.click(); };
$('#dlVib').onclick = () => { if (!vibViewer || MODESEL < 0) { toast('select a mode first'); return; } const a = document.createElement('a'); a.href = vibViewer.pngURI(); a.download = `mode_${MODES.modes[MODESEL].freq.toFixed(0)}cm-1.png`; a.click(); };
$('#dlGeom').onclick = () => { if (!GEOM) return; const q = v => '"' + String(v ?? '') + '"'; const rows = [['bond', 'order', 'length_A'], ...GEOM.bond_lengths.map(b => [b.label, b.aromatic ? 'aromatic' : b.order, b.length]), [], ['angle', 'degrees'], ...GEOM.bond_angles.map(b => [b.label, b.angle])]; download('geometry.csv', rows.map(r => r.map(q).join(',')).join('\n'), 'text/csv'); };
$('#dlMol').onclick = () => { if (!DEP) return; download(((R && R.common_name) || 'compound').replace(/[^\w\-]+/g, '_') + '.mol', DEP.molblock, 'chemical/x-mdl-molfile'); };
$('#dlPng').onclick = () => { if (!DEP) return; const svg = structureSVG(); const img = new Image(); const [x, y, w, h] = bbox(DEP.coords, 40); img.onload = () => { const c = document.createElement('canvas'); c.width = w * 3; c.height = h * 3; const g = c.getContext('2d'); g.fillStyle = '#fff'; g.fillRect(0, 0, c.width, c.height); g.drawImage(img, 0, 0, c.width, c.height); const a = document.createElement('a'); a.href = c.toDataURL('image/png'); a.download = ((R && R.common_name) || 'compound').replace(/[^\w\-]+/g, '_') + '.png'; a.click(); }; img.src = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(svg))); };


// ---------------- structure editor (Ketcher standalone, vendored, Apache-2.0)
function ketcherReady() {
  const fr = $('#ketcher');
  return new Promise((res, rej) => {
    let tries = 0;
    const poll = () => { const k = fr.contentWindow && fr.contentWindow.ketcher; if (k && typeof k.getSmiles === 'function') return res(k); if (++tries > 300) return rej(new Error('Ketcher did not initialise (static/vendor/ketcher missing?)')); setTimeout(poll, 100); };
    poll();
  });
}
async function openEditor(blank) {
  $('#editModal').classList.add('open'); $('#editStatus').innerHTML = '<span class="spin"></span> loading editor…';
  const fr = $('#ketcher'); if (fr.getAttribute('src') === 'about:blank') fr.src = '/static/vendor/ketcher/index.html';
  try {
    const k = await ketcherReady();
    if (blank || !DEP) { try { await k.setMolecule(''); } catch (e) { try { k.editor && k.editor.clear(); } catch (e2) {} } $('#editStatus').textContent = 'draw a molecule, then apply'; }
    else { await k.setMolecule(DEP.molblock); $('#editStatus').textContent = 'edit, then apply'; }
  } catch (e) { $('#editStatus').textContent = e.message; }
}
$('#drawBtn').onclick = () => openEditor(true);
$('#editBtn').onclick = async () => {
  if (!DEP) return; $('#editModal').classList.add('open'); $('#editStatus').innerHTML = '<span class="spin"></span> loading editor…';
  const fr = $('#ketcher'); if (fr.getAttribute('src') === 'about:blank') fr.src = '/static/vendor/ketcher/index.html';
  try { const k = await ketcherReady(); await k.setMolecule(DEP.molblock); $('#editStatus').textContent = 'edit, then apply'; } catch (e) { $('#editStatus').textContent = e.message; }
};
$('#editCancel').onclick = () => $('#editModal').classList.remove('open');
$('#editApply').onclick = async () => {
  try {
    const k = await ketcherReady(); const smi = await k.getSmiles(); if (!smi) { $('#editStatus').textContent = 'empty structure'; return; }
    $('#editModal').classList.remove('open'); $('#q').value = smi; resolve(); toast('Structure applied');
  } catch (e) { $('#editStatus').textContent = e.message; }
};
