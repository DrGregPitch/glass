// Interactive periodic table (PubChem public-domain data): category / property colour maps, state-at-temperature, detail card.
let ELEMENTS = null, PT_PIN = null;
const PT_CAT = {
  'Alkali metal': '#ffd6a5', 'Alkaline earth metal': '#ffe8b3', 'Transition metal': '#dbe7ff', 'Post-transition metal': '#d7f2e3', 'Metalloid': '#e9e3ff',
  'Nonmetal': '#d1f0f6', 'Halogen': '#ffe0e0', 'Noble gas': '#e6d9ff', 'Lanthanide': '#ffd9ec', 'Actinide': '#f6d6e9',
};
const PT_UNITS = { en: '', mass: ' u', radius: ' pm', ie: ' eV', ea: ' eV', mp: ' K', bp: ' K', density: ' g/cm³', year: '' };
const PT_NAMES = { en: 'Electronegativity', mass: 'Atomic mass', radius: 'Atomic radius', ie: 'Ionisation energy', ea: 'Electron affinity', mp: 'Melting point', bp: 'Boiling point', density: 'Density', year: 'Discovered' };

async function initElements() {
  if (ELEMENTS) return;
  const d = await api('/api/elements', {}); ELEMENTS = d.elements;
  const g = $('#ptgrid'); let html = '';
  const byPos = {}; ELEMENTS.forEach(e => byPos[e.row + ',' + e.col] = e);
  for (let r = 1; r <= 10; r++) {
    if (r === 8) { html += `<div style="grid-column:1/19;height:8px"></div>`; continue; }
    for (let c = 1; c <= 18; c++) {
      const e = byPos[r + ',' + c];
      if (!e) { if ((r === 9 || r === 10) && c === 2) html += `<div class="ptseries" style="grid-column:1/3;grid-row:${r}">${r === 9 ? 'lanthanides' : 'actinides'}</div>`; else if (!((r === 9 || r === 10) && c === 1)) html += `<div class="el gap" style="grid-row:${r};grid-column:${c}"></div>`; continue; }
      html += `<div class="el" data-z="${e.z}" style="grid-row:${r};grid-column:${c}" title="${e.name} · ${e.config}"><div class="z">${e.z}</div><div class="sy">${e.symbol}</div><div class="nm">${e.name}</div><div class="cf">${esc(e.config)}</div><div class="mv"></div></div>`;
    }
  }
  g.innerHTML = html;
  g.querySelectorAll('.el[data-z]').forEach(el => {
    const e = ELEMENTS[+el.dataset.z - 1];
    el.onmouseenter = () => { if (!PT_PIN) ptCard(e); };
    el.onclick = () => { PT_PIN = e; g.querySelectorAll('.el.pin').forEach(x => x.classList.remove('pin')); el.classList.add('pin'); ptCard(e); openAtom(e); };
  });
  g.onmouseleave = () => { if (PT_PIN) ptCard(PT_PIN); };
  $('#ptprop').onchange = ptColor;
  let raf = 0; $('#ptT').oninput = () => { $('#ptTval').textContent = $('#ptT').value + ' K'; if (!raf) raf = requestAnimationFrame(() => { raf = 0; ptColor(); }); };
  $$('#ptTlab .presets button').forEach(b => b.onclick = () => { $('#ptT').value = b.dataset.t; $('#ptTval').textContent = b.dataset.t + ' K'; ptColor(); });
  $('#ptsearch').oninput = () => { const q = $('#ptsearch').value.trim().toLowerCase(); g.querySelectorAll('.el[data-z]').forEach(el => { const e = ELEMENTS[+el.dataset.z - 1]; const hit = !q || e.symbol.toLowerCase() === q || e.name.toLowerCase().startsWith(q) || String(e.z) === q; el.classList.toggle('dim', !hit); }); };
  ptColor();
  ptCard(ELEMENTS[5]);
}
function ptScale(prop) {
  const vals = ELEMENTS.map(e => prop === 'year' ? (isNaN(+e.year) ? null : +e.year) : e[prop]).filter(v => v != null);
  let lo = Math.min(...vals), hi = Math.max(...vals);
  const log = prop === 'density' || prop === 'mass';
  return { lo, hi, log, f: v => { if (v == null) return null; let t = log ? (Math.log(v + 0.01) - Math.log(lo + 0.01)) / (Math.log(hi + 0.01) - Math.log(lo + 0.01)) : (v - lo) / (hi - lo); return Math.max(0, Math.min(1, t)); } };
}
// viridis (perceptually uniform, colour-blind safe): low = deep purple, high = yellow
const VIRIDIS = [[68, 1, 84], [72, 40, 120], [62, 74, 137], [49, 104, 142], [38, 130, 142], [31, 158, 137], [53, 183, 121], [109, 205, 89], [180, 222, 44], [253, 231, 37]];
function heat(t) {
  const x = t * (VIRIDIS.length - 1), i = Math.min(VIRIDIS.length - 2, Math.floor(x)), f = x - i;
  const c = VIRIDIS[i].map((v, k) => Math.round(v + (VIRIDIS[i + 1][k] - v) * f));
  return `rgb(${c.join(',')})`;
}
const SBLOCK = { s: '#ffd6a5', p: '#d1f0f6', d: '#dbe7ff', f: '#ffd9ec' };
function sblock(e) { if ([2].includes(e.z)) return 's'; if (e.block === 'Lanthanide' || e.block === 'Actinide') return 'f'; if (e.col <= 2) return 's'; if (e.col >= 13) return 'p'; return 'd'; }
function sig(v) { if (v == null) return '·'; if (v === 0) return '0'; if (Math.abs(v) >= 100) return String(Math.round(v)); if (Math.abs(v) >= 1) return String(+v.toFixed(2)); return String(+v.toPrecision(3)); }
function ptColor() {
  const prop = $('#ptprop').value; $('#ptTlab').style.display = prop === 'state' ? '' : 'none';
  const T = +$('#ptT').value; const leg = $('#ptlegend');
  const sc = (prop !== 'block' && prop !== 'sblock' && prop !== 'state') ? ptScale(prop) : null;
  $$('#ptgrid .el[data-z]').forEach(el => {
    const e = ELEMENTS[+el.dataset.z - 1]; const mv = el.querySelector('.mv'); let bg = '#f0efe9', txt = '';
    el.classList.remove('liq-q', 'dark', 'na');
    if (prop === 'block') { bg = PT_CAT[e.block] || '#eee'; txt = e.mass != null ? e.mass.toFixed(e.mass > 100 ? 1 : 2) : ''; }
    else if (prop === 'sblock') { bg = SBLOCK[sblock(e)]; txt = sblock(e) + '-block'; }
    else if (prop === 'state') { if (e.mp == null) { el.classList.add('na'); txt = 'unknown'; } else if (T < e.mp) { bg = '#dbe7ff'; txt = 'solid'; } else if (e.bp != null && T >= e.bp) { bg = '#ffe0e0'; txt = 'gas'; } else if (e.bp == null) { bg = '#d7f2e3'; txt = 'liquid?'; el.classList.add('liq-q'); } else { bg = '#d7f2e3'; txt = 'liquid'; } }
    else {
      const v = prop === 'year' ? (isNaN(+e.year) ? null : +e.year) : e[prop]; const t = sc.f(v);
      if (t == null) { el.classList.add('na'); txt = prop === 'year' && e.year ? e.year.toLowerCase() : 'no data'; }
      else { bg = heat(t); txt = prop === 'year' ? v : sig(v) + PT_UNITS[prop]; if (t < 0.6) el.classList.add('dark'); }
    }
    el.style.background = bg; mv.textContent = txt;
  });
  const fmt = sig;
  if (prop === 'block') leg.innerHTML = Object.entries(PT_CAT).map(([k, c]) => `<span class="sw"><i style="background:${c}"></i>${k}</span>`).join('') + '<span class="hint">· cell value: atomic mass (u)</span>';
  else if (prop === 'sblock') leg.innerHTML = Object.entries(SBLOCK).map(([k, c]) => `<span class="sw"><i style="background:${c}"></i>${k}-block</span>`).join('') + '<span class="hint">· highest-energy subshell being filled (He counted as s)</span>';
  else if (prop === 'state') leg.innerHTML = `<span class="sw"><i style="background:#dbe7ff"></i>solid</span><span class="sw"><i style="background:#d7f2e3"></i>liquid</span><span class="sw"><i style="background:#d7f2e3;background-image:repeating-linear-gradient(45deg,transparent 0 3px,rgba(0,0,0,.15) 3px 4px)"></i>liquid? (bp unknown)</span><span class="sw"><i style="background:#ffe0e0"></i>gas</span><span class="sw"><i style="background:repeating-linear-gradient(45deg,#f3f2ec 0 3px,#e6e4dc 3px 4px)"></i>unknown mp</span> <span class="hint">at 1 atm, from tabulated mp/bp</span>`;
  else {
    const mid = sc.log ? Math.exp((Math.log(sc.lo + 0.01) + Math.log(sc.hi + 0.01)) / 2) : (sc.lo + sc.hi) / 2;
    leg.innerHTML = `<span><b>${PT_NAMES[prop]}</b>${PT_UNITS[prop] ? ' (' + PT_UNITS[prop].trim() + ')' : ''}${sc.log ? ', log scale' : ''}</span><span class="col"><span class="bar" style="background:linear-gradient(90deg,${[0, .2, .4, .6, .8, 1].map(heat).join(',')})"></span><span class="ticks"><span>${fmt(sc.lo)}</span><span>${fmt(mid)}</span><span>${fmt(sc.hi)}</span></span></span><span class="sw"><i style="background:repeating-linear-gradient(45deg,#f3f2ec 0 3px,#e6e4dc 3px 4px)"></i>no data</span>`;
  }
}
function ptCard(e) {
  const c = $('#ptcard'); const cat = PT_CAT[e.block] || '#eee';
  const f = (v, u, d = 3) => v == null ? '<span class="hint">·</span>' : `${+v.toFixed(d)}${u}`;
  const K = v => v == null ? '<span class="hint">·</span>' : `${v} K <span class="hint">(${(v - 273.15).toFixed(0)} °C)</span>`;
  c.innerHTML = `<div class="big"><span class="sym" style="color:${e.cpk ? '#' + e.cpk : 'var(--ink)'};-webkit-text-stroke:.5px rgba(0,0,0,.25)">${e.symbol}</span><span class="zn">Z = ${e.z}</span></div><h3>${e.name}</h3><span class="cat" style="background:${cat}">${e.block}</span>
    <div class="kv"><b>Atomic mass</b><span>${f(e.mass, ' u', 4)}</span><b>Configuration</b><span class="mono" style="font-size:12px">${esc(e.config)}</span><b>Electronegativity</b><span>${f(e.en, '', 2)}</span><b>Atomic radius</b><span>${f(e.radius, ' pm', 0)}</span><b>1st ionisation</b><span>${f(e.ie, ' eV')}</span><b>Electron affinity</b><span>${f(e.ea, ' eV')}</span><b>Oxidation states</b><span>${esc(e.ox || '·')}</span><b>Standard state</b><span>${esc(e.state || '·')}</span><b>Melting point</b><span>${K(e.mp)}</span><b>Boiling point</b><span>${K(e.bp)}</span><b>Density</b><span>${e.density == null ? '<span class="hint">·</span>' : sig(e.density) + ' g/cm³' + (e.density < 0.01 ? ' <span class="hint">(gas, ' + sig(e.density * 1000) + ' g/L)</span>' : '')}</span><b>Discovered</b><span>${esc(e.year || '·')}</span></div>
    <div class="ctl" style="margin-top:12px"><button class="small" data-el="${e.symbol}">open in Glass</button><span class="hint">resolves the elemental species</span></div>`;
  c.querySelector('[data-el]').onclick = () => { $('#q').value = `[${e.symbol}]`; resolve(); };
}

// ---------------- single-atom viewer (shells + hydrogenic orbitals)
let atomViewer = null, ATOM = null;
async function openAtom(e) {
  $('#atomModal').classList.add('open');
  $('#atomTitle').textContent = `${e.name}  ·  ${e.symbol}  ·  Z = ${e.z}`;
  $('#atomSub').textContent = ''; $('#aoSeg').innerHTML = ''; $('#atomShells').innerHTML = '';
  try { ATOM = await api('/api/atom', { z: e.z }); } catch (err) { $('#atomSub').textContent = err.message; return; }
  $('#atomSub').innerHTML = `<span class="mono">${esc(ATOM.config)}</span> · ${ATOM.note}`;
  drawShells(ATOM);
  $('#atomKV').innerHTML = `<b>Configuration</b><span class="mono">${esc(ATOM.config)}</span><b>Valence</b><span>${ATOM.valence.n}${ATOM.valence.l} <span class="hint">(${ATOM.valence.electrons} e⁻)</span></span><b>Z<sub>eff</sub> (valence)</b><span>${ATOM.z_eff}</span>`;
  // subshell buttons: unique (n,l) present, valence last & selected
  const seen = new Set(); const subs = [];
  ATOM.config.replace(/\[[A-Za-z]+\]/g, '').match(/\d[spdf]/g)?.forEach(t => { if (!seen.has(t)) { seen.add(t); subs.push(t); } });
  // include noble-core representative shells too via shells list
  $('#aoSeg').innerHTML = subs.map((t, i) => `<button data-n="${t[0]}" data-l="${t[1]}" class="${i === subs.length - 1 ? 'on' : ''}">${t}</button>`).join('');
  $('#aoSeg').querySelectorAll('button').forEach(b => b.onclick = () => { $('#aoSeg').querySelectorAll('button').forEach(x => x.classList.remove('on')); b.classList.add('on'); showAO(e.z, +b.dataset.n, b.dataset.l); });
  const last = $('#aoSeg').querySelector('button.on');
  if (last) showAO(e.z, +last.dataset.n, last.dataset.l);
}
function drawShells(a) {
  const cx = 120, cy = 120, R0 = 26, dR = (108 - R0) / Math.max(1, a.shells.length);
  const nucCol = '#334155';
  let out = `<circle cx="${cx}" cy="${cy}" r="18" fill="${nucCol}"/><text x="${cx}" y="${cy + 4}" text-anchor="middle" font-size="12" font-weight="700" fill="#fff">${a.z}</text>`;
  a.shells.forEach((sh, i) => {
    const r = R0 + (i + 1) * dR;
    out += `<circle cx="${cx}" cy="${cy}" r="${r}" fill="none" stroke="var(--muted)" stroke-opacity=".4" stroke-width="1"/>`;
    for (let k = 0; k < sh.electrons; k++) {
      const ang = (k / sh.electrons) * 2 * Math.PI - Math.PI / 2;
      out += `<circle cx="${(cx + r * Math.cos(ang)).toFixed(1)}" cy="${(cy + r * Math.sin(ang)).toFixed(1)}" r="2.6" fill="var(--acc)"/>`;
    }
    out += `<text x="${cx + r}" y="${cy - 3}" font-size="8" fill="var(--muted)">${sh.electrons}</text>`;
  });
  $('#atomShells').innerHTML = out;
}
async function showAO(z, n, l) {
  $('#atomONote').innerHTML = '<span class="spin"></span> building orbital…';
  try {
    const j = await api('/api/atom_orbital', { z, n, l });
    if (!atomViewer) atomViewer = $3Dmol.createViewer($('#atomViewer'), { backgroundColor: 'white' });
    try { atomViewer.resize(); } catch (e) {}
    atomViewer.clear();
    atomViewer.addVolumetricData(j.cube, 'cube', { isoval: j.iso, color: '#c2410c', opacity: 0.82, smoothness: 6 });
    atomViewer.addVolumetricData(j.cube, 'cube', { isoval: -j.iso, color: '#1d4ed8', opacity: 0.82, smoothness: 6 });
    atomViewer.zoomTo();
    atomViewer.rotate(72, {x: 1, y: 0.35, z: 0});   // 3/4 view so z-axis orbitals (p_z, d_z2) show their lobes
    atomViewer.render();
    $('#atomONote').textContent = `${j.orbital} orbital: ${j.n_lobes > 1 ? j.n_lobes + ' orientations exist; showing one' : 'spherically symmetric'}. Red/blue = wavefunction sign (phase).`;
  } catch (e) { $('#atomONote').textContent = 'orbital unavailable: ' + e.message; }
}
