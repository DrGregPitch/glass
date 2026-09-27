// NMR / IR / UV panels, solvent + frequency controls, MMFF94s normal-mode animation.
// Load order matters only for the top-level statements; all files share the page's global scope.

// ---------------- spectra
async function loadSpectra() {
  if (!R) return;
  const mhz = +$('#mhz').value || 400; $('#c13mhz').textContent = `(¹³C at ${(mhz / 3.976).toFixed(1)} MHz)`;
  try { SPEC = await api('/api/spectra', { smiles: R.smiles_canonical, solvent: $('#solvent').value, mhz, nmr: $('#nmrlevel').value }); } catch (e) { return; }
  levelNote();
  if (SPEC.giao_note) $('#levelnote').textContent += '  ·  ' + SPEC.giao_note;
  SPEC.h1.peaks.sort((a, b) => b.shift - a.shift);
  renderSpectraTables();
  IRDATA = SPEC.ir; drawIR();
  $('#t_ir').innerHTML = '<tr><th>ν̃ / cm⁻¹</th><th>int.</th><th>assignment</th></tr>' + SPEC.ir.bands.map(b => `<tr><td>${b.center}</td><td>${b.intensity > 0.7 ? 's' : b.intensity > 0.4 ? 'm' : 'w'}${b.width > 100 ? ', br' : ''}</td><td>${esc(b.label)}</td></tr>`).join('');
  loadModes().then(loadPhoto);
  const uv = SPEC.uv;
  chart($('#c_uv'), uv.x, uv.y, { xlabel: 'λ / nm', ylabel: 'ε', marks: uv.peaks.filter(p => p.lambda_max).map(p => [p.lambda_max, `${p.lambda_max} nm`]) });
  $('#t_uv').innerHTML = '<tr><th>λmax / nm</th><th>ε (approx.)</th><th>assignment</th></tr>' + uv.peaks.map(p => `<tr><td>${p.lambda_max ?? '·'}</td><td>${p.epsilon ?? '·'}</td><td>${esc(p.assignment)}</td></tr>`).join('');
  $('#uvnote').textContent = (uv.colour ? uv.colour + '. ' : '') + uv.note;
  $('#uvsolv').textContent = uv.solvent_note || '';
  $('#irnote').textContent = (SPEC.ir.solvent_note || '') + ' · red sticks: normal modes (gas phase), click to animate';
  const sv = SPEC.solvent; $('#solvnote').textContent = `${sv.nmr ? 'NMR in ' + sv.nmr + ' · ' : 'no common deuterated form for NMR · '}ε = ${sv.eps}, E_T(30) = ${sv.et30}, Δf = ${sv.delta_f}${sv.protic ? ' · protic' : ''}; applied to NMR exchangeables, IR C=O / X–H, UV–Vis solvatochromism and the emission Stokes shift`;
}
function renderSpectraTables() {
  const h = SPEC.h1, c = SPEC.c13;
  $('#t_h1').innerHTML = '<tr><th>δ / ppm</th><th>mult.</th><th>J / Hz</th><th>nH</th><th>assignment</th><th>protons</th><th>basis</th></tr>' + h.peaks.map((p, i) => `<tr class="pk" data-k="${i}"><td>${p.shift.toFixed(2)}</td><td>${p.multiplicity}</td><td>${(p.J || []).join(', ')}</td><td>${p.nH}H</td><td>${esc(p.type)}${p.exchangeable ? ' <span class="hint">(exch.)</span>' : ''}</td><td class="muted">${p.atoms.map(a => hLabel(a) || ('H on #' + a)).join(', ')}</td><td>${basisCell(p)}</td></tr>`).join('');
  nmrChart($('#c_h1'), $('#t_h1'), h.x, h.y, h.peaks.map(p => ({ x: p.shift, label: p.shift.toFixed(2), atoms: p.atoms })), { vline: h.residual_solvent_ppm, decimals: 2, kind: 'h' });
  $('#t_c13').innerHTML = '<tr><th>δ / ppm</th><th>DEPT</th><th>type</th><th>atoms</th><th>basis</th></tr>' + c.peaks.map((p, i) => `<tr class="pk" data-k="${i}"><td>${p.shift.toFixed(1)}</td><td>${p.dept}</td><td>${esc(p.type)}</td><td class="muted">${p.atoms.map(a => atomLabel(a)).join(', ')}</td><td>${basisCell(p)}</td></tr>`).join('');
  nmrChart($('#c_c13'), $('#t_c13'), c.x, c.y, c.peaks.map(p => ({ x: p.shift, label: p.shift.toFixed(1), atoms: p.atoms })), { vline: c.solvent_c13_ppm, decimals: 1, kind: 'c' });
  renderExpSection();
}

// journal-style experimental section from the predicted peaks (draft; user must verify)
function expSectionText() {
  const h = SPEC.h1, c = SPEC.c13, mhz = +$('#mhz').value || 400;
  const c13mhz = Math.round(mhz * 0.25145);
  const solv = (SPEC.solvent && SPEC.solvent.nmr) || 'solvent';
  const name = R.common_name || R.iupac_name || R.formula;
  const hs = [...h.peaks].sort((a, b) => b.shift - a.shift).map(p => {
    const parts = [p.multiplicity || 'm'];
    if (p.J && p.J.length) parts.push('J = ' + p.J.map(j => j.toFixed(1)).join(', ') + ' Hz');
    parts.push(p.nH + 'H');
    if (p.exchangeable) parts.push(p.type);
    return p.shift.toFixed(2) + ' (' + parts.join(', ') + ')';
  }).join(', ');
  const cs = [...c.peaks].sort((a, b) => b.shift - a.shift).map(p => p.shift.toFixed(1)).join(', ');
  const mw = (typeof PROPS !== 'undefined' && PROPS && PROPS.molar_mass) ? ' (' + PROPS.formula + ', M = ' + PROPS.molar_mass.value.toFixed(2) + ' g/mol)' : '';
  let t = name + mw + '.';
  if (h.peaks.length) t += ' \u00b9H NMR (' + mhz + ' MHz, ' + solv + ') \u03b4 ' + hs + '.';
  if (c.peaks.length) t += ' \u00b9\u00b3C{\u00b9H} NMR (' + c13mhz + ' MHz, ' + solv + ') \u03b4 ' + cs + '.';
  return t;
}
function renderExpSection() {
  const el = $('#expsec'); if (!el || !SPEC || !R) return;
  el.textContent = expSectionText();
}
$('#copyExpSec').onclick = async () => { try { await navigator.clipboard.writeText(expSectionText()); toast('experimental section copied'); } catch (e) { toast('copy failed'); } };
function basisCell(p) { if (p.source === 'giao') return `<span class="badge acc plain" title="quantum shielding, empirically scaled; estimated ± ${p.sigma} ppm">quantum · ±${p.sigma}</span>`; return p.source === 'nmrshiftdb2' ? `<span class="badge ok plain" title="experimental-database match: ${p.n_ref} reference atoms, spread ± ${p.sigma} ppm">exp. db · n=${p.n_ref} · ±${p.sigma}</span>` : '<span class="badge info plain" title="no database match for this environment; additivity increments">increments</span>'; }
function atomLabel(i) { const loc = SEL && SEL.locants && SEL.locants[i]; const f = SEL && SEL.fragments && SEL.fragments[i]; if (!loc || !hasNum(SEL)) return `#${i}`; return f === 'sub' ? loc + '′' : loc; }
$('#respec').onclick = loadSpectra; $('#solvent').onchange = loadSpectra; $('#mhz').onchange = loadSpectra; $('#nmrlevel').onchange = loadSpectra;
$('#viblevel').onchange = () => { loadModes().then(loadPhoto); levelNote(); };
const LEVEL_HINTS = {
  giao: 'NMR: quantum calculation, empirically scaled. Real physics, no database: works for novel environments; per-peak uncertainty shown; ≤16 heavy atoms, no Br/I; takes seconds to ~a minute.',
  hose: 'NMR: nearest experimental environment from a reference database. Best whenever the environment exists in the database; each peak shows its match quality and spread.',
  increments: 'NMR: additivity increments only (¹³C ±3–5 ppm, ¹H ±0.3 ppm). Useful for novel environments and to see the textbook estimate; no database dependence.',
  auto: 'Vibrations: GFN2-xTB when the molecule has ≤45 heavy atoms, otherwise MMFF94s.',
  gfn2: 'Vibrations: GFN2-xTB. Real electronic structure (charges respond to motion → realistic intensities); best for carbonyl, ring and CH bending modes (±25–35 cm⁻¹); systematically low for O–H and off for C≡N and C–Cl.',
  mmff: 'Vibrations: MMFF94s force field. Fast; parametrised X–H and nitrile stretches can be closer to experiment than GFN2, but intensities are only qualitative and conjugated C=C/C=O can be 30–80 cm⁻¹ off.',
};
function levelNote() { $('#levelnote').textContent = LEVEL_HINTS[$('#nmrlevel').value]; $('#levelnote').title = LEVEL_HINTS[$('#viblevel').value]; }


// ---------------- vibrational modes
let IRDATA = null, MODES = null, MODESEL = -1, vibViewer = null; const VIB = { raf: 0 };
function showPeriod() { const real = 33.36 * 1000 / Math.max(MODES.modes[MODESEL].freq, 20); $('#vibperiod').textContent = `animation fixed at 1.05 s per cycle for every mode; the real period of this mode is ${real.toFixed(1)} fs (≈${(1.05e15 / real).toExponential(1)}× slowed)`; }
function stopVib() { if (VIB.raf) cancelAnimationFrame(VIB.raf); VIB.raf = 0; if (vibViewer) { try { vibViewer.stopAnimate(); } catch (e) {} } }
function drawIR() {
  if (!IRDATA) return;
  const peaks = MODES ? MODES.modes.map(m => ({ x: m.freq, h: m.intensity, label: `${m.freq.toFixed(0)}`, atoms: [] })) : [];
  nmrChart($('#c_ir'), $('#t_modes'), IRDATA.x, IRDATA.y, peaks, { xlabel: 'wavenumber / cm⁻¹', ylabel: '%T', ymax: 100, sticks: true, rowSelector: 'tr.mode', decimals: 0, onHover: k => {}, onSelect: k => selectMode(k) });
  if (MODESEL >= 0 && $('#c_ir')._setSel) $('#c_ir')._setSel(MODESEL);
}
async function loadModes() {
  MODES = null; MODESEL = -1; stopVib(); $('#vib3d').style.display = 'none'; $('#vibnote').textContent = ''; $('#vibperiod').textContent = '';
  $('#t_modes').innerHTML = '<tr><td class="hint"><span class="spin"></span> computing normal modes (GFN2-xTB, a few seconds)…</td></tr>'; $('#modelevel').innerHTML = '';
  try { MODES = await api('/api/modes', { smiles: R.smiles_canonical, level: $('#viblevel').value }); } catch (e) { $('#t_modes').innerHTML = `<tr><td class="hint">${esc(e.message)}</td></tr>`; return; }
  if (MODES.error) { $('#t_modes').innerHTML = `<tr><td class="hint">${esc(MODES.error)}</td></tr>`; return; }
  $('#modelevel').innerHTML = `<span class="badge ${MODES.level === 'GFN2-xTB' ? 'ok' : 'warn'} plain">${esc(MODES.level || 'MMFF94s')}</span>${MODES.fallback_reason ? ` <span class="hint">${esc(MODES.fallback_reason)}</span>` : ''}`;
  $('#t_modes').innerHTML = `<tr><th>ν̃ / cm⁻¹</th><th>IR int.</th><th>dominant motion</th></tr>` + MODES.modes.map((m, k) => `<tr class="mode" data-k="${k}"><td>${m.freq.toFixed(0)}</td><td>${m.intensity.toFixed(2)}</td><td>${esc(m.label)}</td></tr>`).join('') + `<tr><td colspan="3" class="hint">${esc(MODES.method)}. ${esc(MODES.note)}</td></tr>`;
  drawIR();
}
function selectMode(k) {
  MODESEL = MODESEL === k ? -1 : k;
  $$('#t_modes tr.mode').forEach(r => r.classList.toggle('sel', +r.dataset.k === MODESEL));
  if ($('#c_ir')._setSel) $('#c_ir')._setSel(MODESEL);
  const box = $('#vib3d');
  if (MODESEL < 0) { box.style.display = 'none'; $('#vibperiod').textContent = ''; $('#vibnote').textContent = ''; stopVib(); return; }
  box.style.display = 'block';
  const m = MODES.modes[MODESEL]; const n = MODES.symbols.length; const frames = 32;
  const maxd = Math.max(...m.vec.map(v => Math.hypot(v[0], v[1], v[2]))) || 1;
  const amp = 0.35 / maxd;
  // multi-frame SDF with the real bond table: XYZ frames made 3Dmol re-perceive bonds
  // from distances each frame, so large modes (e.g. ring stretches) looked like bonds breaking
  const pad = (v, w) => String(v).padStart(w);
  const bondBlock = (MODES.bonds || []).map(b => `${pad(b[0] + 1, 3)}${pad(b[1] + 1, 3)}  1  0`).join('\n');
  const nb = (MODES.bonds || []).length;
  let sdf = '';
  for (let f = 0; f < frames; f++) {
    const s = amp * Math.sin(2 * Math.PI * f / frames);
    const atoms = MODES.symbols.map((el, i) => `${(MODES.xyz[i][0] + s * m.vec[i][0]).toFixed(4).padStart(10)}${(MODES.xyz[i][1] + s * m.vec[i][1]).toFixed(4).padStart(10)}${(MODES.xyz[i][2] + s * m.vec[i][2]).toFixed(4).padStart(10)} ${el.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`).join('\n');
    sdf += `mode ${m.freq}\n  Glass\n\n${pad(n, 3)}${pad(nb, 3)}  0  0  0  0  0  0  0  0999 V2000\n${atoms}\n${bondBlock}\nM  END\n$$$$\n`;
  }
  try {
    if (!vibViewer) vibViewer = $3Dmol.createViewer(box, { backgroundColor: 'white' });
    stopVib(); vibViewer.clear();
    vibViewer.addModelsAsFrames(sdf, 'sdf');
    vibViewer.setStyle({}, { stick: { radius: 0.14 }, sphere: { scale: 0.22 } });
    MODES.symbols.forEach((el, i) => { const v = m.vec[i]; const L = Math.hypot(v[0], v[1], v[2]) / maxd; if (L < 0.12) return; const p = MODES.xyz[i]; const kk = 1.4 / maxd; vibViewer.addArrow({ start: { x: p[0], y: p[1], z: p[2] }, end: { x: p[0] + kk * v[0], y: p[1] + kk * v[1], z: p[2] + kk * v[2] }, radius: 0.05, color: '#c0392b' }); });
    vibViewer.zoomTo(); vibViewer.setFrame(0); vibViewer.render();
    VIB.base = 1050; VIB.phase = 0; VIB.t = performance.now(); let last = -1; showPeriod();
    const tick = now => { VIB.phase += (now - VIB.t) / VIB.base; VIB.t = now; const f = Math.floor((VIB.phase % 1) * frames); if (f !== last) { last = f; vibViewer.setFrame(f); vibViewer.render(); } VIB.raf = requestAnimationFrame(tick); };
    VIB.raf = requestAnimationFrame(tick);
    $('#vibnote').textContent = `mode ${MODESEL + 1}: ${m.freq.toFixed(0)} cm⁻¹ (unscaled ${m.freq_unscaled.toFixed(0)}), ${m.label}; arrows = displacement vectors`;
  } catch (e) { box.innerHTML = '<div class="hint" style="padding:10px">3D viewer unavailable (3Dmol.js CDN)</div>'; }
}

// redraw the charts of a pane after it becomes visible (hidden canvases draw at fallback size)
function redrawPane(n) {
  if (!SPEC) return;
  try {
    if (n === 'nmr') renderSpectraTables();
    else if (n === 'ir') { drawIR(); if (typeof vibViewer !== 'undefined' && vibViewer) { try { vibViewer.resize(); vibViewer.render(); } catch (e) {} } }
    else if (n === 'uv') { const uv = SPEC.uv; chart($('#c_uv'), uv.x, uv.y, { xlabel: 'λ / nm', ylabel: 'ε', marks: uv.peaks.filter(p => p.lambda_max).map(p => [p.lambda_max, `${p.lambda_max} nm`]) }); }
    else if (n === 'pl' && typeof redrawPL === 'function') redrawPL();
  } catch (e) { console.warn('redrawPane', e); }
}
