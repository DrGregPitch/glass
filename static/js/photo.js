// Photoluminescence panel: vibronic spectra, excitation/emission selection, interactive Franck–Condon diagram.
// Load order matters only for the top-level statements; all files share the page's global scope.

// ---------------- photoluminescence + Franck–Condon
let PL = null;
async function loadPhoto() {
  PL = null; $('#plsummary').innerHTML = '<span class="hint"><span class="spin"></span> estimating…</span>'; $('#t_pl').innerHTML = ''; $('#fc').innerHTML = ''; $('#plreasons').innerHTML = ''; $('#plclass').textContent = '';
  try { PL = await api('/api/photo', { smiles: R.smiles_canonical, solvent: $('#solvent').value, level: $('#viblevel').value }); } catch (e) { $('#plsummary').innerHTML = `<span class="hint">${esc(e.message)}</span>`; return; }
  $('#plclass').textContent = PL.class;
  $('#plreasons').innerHTML = PL.reasons.map(r => `<div class="reason">${esc(r)}</div>`).join('');
  $('#plnote').textContent = PL.note;
  if (!('S' in PL)) { $('#plsummary').innerHTML = `<span class="hint" style="padding:8px">${esc(PL.reasons[0])}</span>`; chart($('#c_pl'), [200, 900], [0, 0], { xlabel: 'λ / nm' }); return; }
  const sw = nm => `<span class="sw" style="background:${wlColor(nm)}"></span>`;
  $('#plsummary').innerHTML = [
    ['Absorption max', `${PL.lambda_abs_max} <span class="u">nm</span>`, sw(PL.lambda_abs_max) + PL.colour_abs],
    ['0-0 origin', `${PL.lambda_00} <span class="u">nm</span>`, `${PL.nu_00} cm⁻¹`],
    ['Emission max', `${PL.lambda_em_max} <span class="u">nm</span>`, sw(PL.lambda_em_max) + PL.colour],
    ['Stokes shift', `${PL.stokes_shift_cm} <span class="u">cm⁻¹</span>`, `${PL.stokes_shift_nm} nm`],
    ['Huang–Rhys S', `${PL.S}`, `ω = ${PL.omega} cm⁻¹`],
    ['Phosphorescence', `~${PL.lambda_phos} <span class="u">nm</span>`, 'T₁ estimate, rigid media'],
  ].map(([k, v, e]) => `<div class="tile"><div class="k">${k}</div><div class="v" style="font-size:17px">${v}</div><div class="e">${e}</div></div>`).join('');
  // vibronic table: absorption then emission
  const rows = PL.absorption.map((p, i) => ({ ...p, kind: 'abs', k: i })).concat(PL.emission.map((p, i) => ({ ...p, kind: 'em', k: PL.absorption.length + i })));
  $('#t_pl').innerHTML = '<tr><th>band</th><th>transition</th><th>λ / nm</th><th>ν̃ / cm⁻¹</th><th>rel. FC intensity</th></tr>' + rows.map(r => `<tr class="pk" data-k="${r.k}"><td><span class="badge plain ${r.kind === 'abs' ? 'acc' : 'ok'}">${r.kind === 'abs' ? 'absorption' : 'emission'}</span></td><td class="mono">${r.label}</td><td>${r.nm.toFixed(1)}</td><td>${r.nu}</td><td>${(r.fc / Math.max(...PL.absorption.map(p => p.fc))).toFixed(2)}</td></tr>`).join('');
  const peaks = rows.map(r => ({ x: r.nm, label: r.label, h: r.fc, kind: r.kind, atoms: [] }));
  drawFC();
  PLPEAKS = peaks;
  redrawPL();
}
let PLPEAKS = [];
function redrawPL() {
  if (!PL || !('S' in PL)) return;
  const lex = +$('#lex').value || null;
  const lem = +$('#lem').value || null;
  const vlines = [];
  if (lex) vlines.push({ x: lex, color: '#b45309', label: 'λex' });
  if (lem) vlines.push({ x: lem, color: '#1d4ed8', label: 'λem' });
  nmrChart($('#c_pl'), $('#t_pl'), PL.x, PL.y_abs, PLPEAKS, { xlabel: 'λ / nm', ylabel: 'normalised intensity', ymax: 1.0, reverseX: false, extra: PL.y_em, extraColor: '#8a8a93', fillWavelength: true, decimals: 1, vlines,
    onHover: k => fcHighlight(k), onSelect: k => fcHighlight(k, true) });
  // feedback text
  const interp = (xs, ys, x) => { if (x < xs[0] || x > xs[xs.length - 1]) return 0; let i = 0; while (xs[i + 1] < x) i++; const t = (x - xs[i]) / (xs[i + 1] - xs[i]); return ys[i] + t * (ys[i + 1] - ys[i]); };
  const parts = [];
  const TOL = 15;   // nm: a wavelength only "hits" a vibronic band within this window
  let fcIdx = -1;   // combined index for fcHighlight: [0..n) absorption, [n..) emission
  if (lex) {
    const A = interp(PL.x, PL.y_abs, lex); const nearest = PL.absorption.reduce((b, p) => Math.abs(p.nm - lex) < Math.abs(b.nm - lex) ? p : b);
    parts.push(`Excite at ${lex} nm → relative absorbance ${A.toFixed(2)} (nearest band ${nearest.label} at ${nearest.nm} nm). Kasha's rule: the emission spectrum is unchanged; its brightness scales with this absorbance${A < 0.05 ? '; essentially no excitation here' : ''}.`);
    if (Math.abs(nearest.nm - lex) <= TOL) fcIdx = PL.absorption.indexOf(nearest);
  }
  if (lem && PL.y_em) {
    const E = interp(PL.x, PL.y_em, lem);
    const nearestEm = (PL.emission || []).length ? PL.emission.reduce((b, p) => Math.abs(p.nm - lem) < Math.abs(b.nm - lem) ? p : b) : null;
    parts.push(`Detect at ${lem} nm → relative emission intensity ${E.toFixed(2)}${nearestEm ? ` (nearest vibronic band ${nearestEm.label} at ${nearestEm.nm} nm)` : ''}${E < 0.05 ? '; essentially no emission here' : ''}.`);
    if (fcIdx < 0 && nearestEm && Math.abs(nearestEm.nm - lem) <= TOL) fcIdx = PL.absorption.length + PL.emission.indexOf(nearestEm);
  }
  $('#plfb').textContent = parts.join(' ');
  fcHighlight(fcIdx, true);
}
$('#lex').oninput = redrawPL; $('#lem').oninput = redrawPL; $('#plreset').onclick = () => { $('#lex').value = ''; $('#lem').value = ''; redrawPL(); };
function wlColor(nm) {   // approximate visible-spectrum colour, grey outside
  if (nm < 380 || nm > 750) return '#c9c7bf';
  let r = 0, g = 0, b = 0;
  if (nm < 440) { r = -(nm - 440) / 60; b = 1; } else if (nm < 490) { g = (nm - 440) / 50; b = 1; } else if (nm < 510) { g = 1; b = -(nm - 510) / 20; } else if (nm < 580) { r = (nm - 510) / 70; g = 1; } else if (nm < 645) { r = 1; g = -(nm - 645) / 65; } else { r = 1; }
  const f = nm < 420 ? 0.3 + 0.7 * (nm - 380) / 40 : nm > 700 ? 0.3 + 0.7 * (750 - nm) / 50 : 1;
  return `rgb(${Math.round(255 * r * f)},${Math.round(255 * g * f)},${Math.round(255 * b * f)})`;
}
function drawFC() {
  const svg = $('#fc'); if (!PL || !('S' in PL)) { svg.innerHTML = ''; return; }
  const W = 520, H = 470, L = 60, Rr = 20, T = 30, B = 40;
  const S = PL.S, d = Math.sqrt(2 * S);              // dimensionless displacement of S1 well
  const n = Math.min(5, PL.absorption.length);      // vibrational levels drawn (0…4)
  const qmin = -4.4, qmax = d + 4.4;
  const X = q => L + (q - qmin) / (qmax - qmin) * (W - L - Rr);
  const wS0 = 1.0, wS1 = 0.95;                        // curvature (S1 slightly softer)
  const E00 = 6.5;                                    // separation between wells in units of ħω (drawing scale only)
  const Emax = E00 + wS1 * (n + 1.2);
  const Y = e => T + (1 - e / (Emax + 0.5)) * (H - T - B);
  let out = '';
  // axes
  out += `<line x1="${L}" y1="${T}" x2="${L}" y2="${H - B}" stroke="#aab3bd"/><line x1="${L}" y1="${H - B}" x2="${W - Rr}" y2="${H - B}" stroke="#aab3bd"/>`;
  out += `<text x="${(L + W - Rr) / 2}" y="${H - 10}" font-size="11" text-anchor="middle" fill="#667380">nuclear coordinate q (effective ${PL.omega} cm⁻¹ mode)</text>`;
  out += `<text transform="translate(16 ${(T + H - B) / 2}) rotate(-90)" font-size="11" text-anchor="middle" fill="#667380">energy</text>`;
  // potential curves
  const well = (q0, w, e0, cls) => { let p = ''; for (let i = 0; i <= 120; i++) { const q = qmin + (qmax - qmin) * i / 120; const e = e0 + 0.5 * w * (q - q0) ** 2; if (e > Emax + 0.5) { p += ' '; continue; } p += (p.endsWith(' ') || !p ? 'M' : 'L') + X(q).toFixed(1) + ',' + Y(e).toFixed(1); } return `<path d="${p.replace(/ +M/g, 'M')}" fill="none" stroke="${cls === 'S0' ? '#1c1c1a' : '#6b6b70'}" stroke-width="2"/>`; };
  out += well(0, wS0, 0, 'S0') + well(d, wS1, E00, 'S1');
  out += `<text x="${X(0) + 8}" y="${Y(-0.35)}" font-size="13" font-weight="700" fill="#1c1c1a">S₀</text>`;
  out += `<text x="${X(d) + 8}" y="${Y(E00 - 0.35)}" font-size="13" font-weight="700" fill="#6b6b70">S₁</text>`;
  // vibrational levels + wavefunctions (harmonic oscillator, first n levels)
  const herm = (v, x) => { let h0 = 1, h1 = 2 * x; if (v === 0) return h0; for (let k = 1; k < v; k++) { const h2 = 2 * x * h1 - 2 * k * h0; h0 = h1; h1 = h2; } return h1; };
  const psi = (v, x) => herm(v, x) * Math.exp(-x * x / 2) / Math.sqrt(Math.pow(2, v) * fact(v) * Math.sqrt(Math.PI));
  function fact(k) { let f = 1; for (let i = 2; i <= k; i++) f *= i; return f; }
  const level = (q0, w, e0, v, id, cls) => {
    const e = e0 + w * (v + 0.5); const half = Math.sqrt(2 * e - 2 * e0) / Math.sqrt(w) ; // classical turning points
    let path = '';
    for (let i = 0; i <= 60; i++) { const x = -half - 0.6 + (2 * half + 1.2) * i / 60; const yy = psi(v, x) * 0.30 * w; path += (i ? 'L' : 'M') + X(q0 + x).toFixed(1) + ',' + Y(e + yy).toFixed(1); }
    return `<line class="lvl" id="lvl-${id}" x1="${X(q0 - half)}" y1="${Y(e)}" x2="${X(q0 + half)}" y2="${Y(e)}"/><path class="wf" id="wf-${id}" d="${path}" stroke="${cls === 'S0' ? '#1c1c1a' : '#6b6b70'}"/><text x="${X(q0 + half) + 4}" y="${Y(e) + 4}" font-size="9" fill="#667380">${v}${cls === 'S1' ? '′' : ''}</text>`;
  };
  for (let v = 0; v < n; v++) out += level(0, wS0, 0, v, 'S0-' + v, 'S0') + level(d, wS1, E00, v, 'S1-' + v, 'S1');
  // arrows: absorption from S0 v=0 at q=0 up to S1 v'=k; emission from S1 v'=0 at q=d down to S0 v=k
  out += `<defs><marker id="ah-abs" markerWidth="8" markerHeight="8" refX="4" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#17171a"/></marker><marker id="ah-em" markerWidth="8" markerHeight="8" refX="4" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#8a8a93"/></marker></defs>`;
  for (let k = 0; k < n; k++) {
    const x = X(0) - 6 + k * 0; out += `<line class="arrow abs${k === 0 ? ' base' : ''}" id="arr-${k}" x1="${X(0) - 10 - k * 3}" y1="${Y(0.5 * wS0)}" x2="${X(0) - 10 - k * 3}" y2="${Y(E00 + wS1 * (k + 0.5)) + 6}" marker-end="url(#ah-abs)"/>`;
    out += `<line class="arrow em${k === 0 ? ' base' : ''}" id="arr-${n + k}" x1="${X(d) + 10 + k * 3}" y1="${Y(E00 + 0.5 * wS1)}" x2="${X(d) + 10 + k * 3}" y2="${Y(wS0 * (k + 0.5)) - 6}" marker-end="url(#ah-em)"/>`;
  }
  out += `<text x="${X(0) - 14}" y="${T + 12}" font-size="11" fill="#17171a" text-anchor="end">absorption 0→n′</text><text x="${X(d) + 14}" y="${T + 12}" font-size="11" fill="#8a8a93">emission 0′→n</text>`;
  out += `<text x="${W - Rr}" y="${H - B - 8}" font-size="10.5" fill="#667380" text-anchor="end">Δq = √(2S) = ${d.toFixed(2)}   S = ${S}</text>`;
  svg.innerHTML = out;
  svg.querySelectorAll('.arrow').forEach(a => { const kd = +a.id.split('-')[1]; const k = kd < n ? kd : PL.absorption.length + (kd - n); a.style.cursor = 'pointer'; a.onmouseenter = () => { fcHighlight(k); if ($('#c_pl')._setHover) $('#c_pl')._setHover(k); }; a.onmouseleave = () => { fcHighlight(-1); if ($('#c_pl')._setHover) $('#c_pl')._setHover(-1); }; });
}
function fcHighlight(k, sticky) {
  const n = PL ? PL.absorption.length : 0; const nd = Math.min(5, n);
  $$('#fc .arrow, #fc .lvl, #fc .wf').forEach(e => e.classList.remove('on'));
  if (k < 0) return;
  const isAbs = k < n; const v = isAbs ? k : k - n; if (v >= nd) return;
  const a = $('#arr-' + (isAbs ? v : nd + v)); if (a) a.classList.add('on');
  const ids = isAbs ? ['S0-0', 'S1-' + v] : ['S1-0', 'S0-' + v];
  ids.forEach(id => { const l = $('#lvl-' + id), w = $('#wf-' + id); if (l) l.classList.add('on'); if (w) w.classList.add('on'); });
}
