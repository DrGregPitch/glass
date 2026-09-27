// Glass front-end: shared state, API client, landing page, resolution and identity/names rendering.
// Load order matters only for the top-level statements; all files share the page's global scope.

const $ = s => document.querySelector(s), $$ = s => [...document.querySelectorAll(s)];
const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
let SOLVS = [], METHODS = {}, R = null, DEP = null, DEPH = null, SEL = null, hoverLoc = null, hoverAtom = -1, PROPS = null, SPEC = null, GEOM = null;
const EXAMPLES = ['anthracene', 'azulene', 'trans-stilbene', 'indigo', 'pyrene', 'coumarin', 'pentacene', 'thiophene'];

async function api(path, params, body) {
  const u = new URL(path, location.origin);
  for (const [k, v] of Object.entries(params || {})) u.searchParams.set(k, v);
  const r = await fetch(u, body ? { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) } : undefined);
  if (!r.ok) { let m = r.statusText; try { m = (await r.json()).detail || m; } catch {} throw new Error(m); }
  return r.json();
}
function toast(msg) { const t = $('#toast'); t.textContent = msg; t.classList.add('show'); clearTimeout(t._h); t._h = setTimeout(() => t.classList.remove('show'), 1400); }
function copyBtn(t) { return `<span class="copy" data-copy="${esc(t)}">copy</span>`; }
document.addEventListener('click', e => {
  const c = e.target.closest('[data-copy]');
  if (c) { navigator.clipboard.writeText(c.dataset.copy); toast('Copied'); }
});
function setStatus(txt, busy) { $('#status').innerHTML = (busy ? '<span class="spin"></span>' : '') + esc(txt); }
function download(name, text, type) { const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name; a.click(); }


// ---------------- modes / landing
function showMode(m) {
  const isEl = m === 'elements';
  $('#tab_mol').classList.toggle('on', !isEl); $('#tab_el').classList.toggle('on', isEl);
  $('.search').style.display = isEl ? 'none' : ''; $('#go').style.display = isEl ? 'none' : ''; $('#usepc').parentElement.style.display = isEl ? 'none' : '';
  $('#elements').style.display = isEl ? 'block' : 'none';
  if (isEl && typeof initElements === 'function') initElements();
  const single = !isEl;
  $('#landing').style.display = single && !R ? 'block' : 'none';
  $$('#sec-identity,#sec-names').forEach(s => s.style.display = single && R ? '' : 'none');
  $('#paneWrap').style.display = single && R ? 'flex' : 'none';
  if (single && R) showPane(PANE);
  $('#subnav').style.display = single && R ? 'block' : 'none';
  }
$('#tab_mol').onclick = () => showMode('single'); $('#tab_el').onclick = () => showMode('elements');
$('#brandHome').onclick = () => {
  R = null; SEL = null; DEP = null; DEPH = null; GEOM = null; PROPS = null; SPEC = null;
  $('#q').value = ''; $('#error').style.display = 'none';
  history.replaceState(null, '', '/');
  document.title = 'Glass';
  renderLanding(); showMode('single');
  window.scrollTo({ top: 0 });
};
function recents() { try { return JSON.parse(localStorage.getItem('cxr_recent') || '[]'); } catch { return []; } }
function pushRecent(q) { try { const r = [q, ...recents().filter(x => x !== q)].slice(0, 8); localStorage.setItem('cxr_recent', JSON.stringify(r)); } catch {} }
function renderLanding() {
  $('#examples').innerHTML = '<span class="hint" style="align-self:center">Try:</span>' + EXAMPLES.map(e => `<span class="chip" data-q="${esc(e)}">${esc(e)}</span>`).join('');
  const r = recents();
  $('#recentbox').style.display = r.length ? 'block' : 'none';
  $('#recent').innerHTML = r.map(e => `<span class="chip" data-q="${esc(e)}">${esc(e)}</span>`).join('');
  $$('.chip[data-q]').forEach(c => c.onclick = () => { $('#q').value = c.dataset.q; resolve(); });
}


// ---------------- file drop (CDXML / MOL / SDF) → parse → resolve
document.addEventListener('dragover', e => { e.preventDefault(); });
document.addEventListener('drop', async e => {
  e.preventDefault();
  const f = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0]; if (!f) return;
  if (!/\.(cdxml|mol|sdf|mdl)$/i.test(f.name)) { toast('Drop a .cdxml or .mol file'); return; }
  toast('Reading ' + f.name + '…');
  try {
    const r = await fetch('/api/read?fmt=' + (/\.cdxml$/i.test(f.name) ? 'cdxml' : 'mol'), { method: 'POST', headers: { 'Content-Type': 'text/plain' }, body: await f.text() });
    if (!r.ok) throw new Error((await r.json()).detail);
    const j = await r.json();
    $('#q').value = j.smiles; resolve();
  } catch (err) { toast('Could not read file: ' + err.message); }
});

// ---------------- resolve
async function resolve() {
  const q = $('#q').value.trim(); if (!q) return;
  showMode('single'); setStatus('resolving…', true); $('#error').style.display = 'none';
  try {
    R = await api('/api/resolve', { q, pubchem: $('#usepc').checked });
  } catch (err) {
    setStatus('', false); showError(q, err.message); return;
  }
  setStatus('', false); pushRecent(q);
  PROPS = SPEC = GEOM = null;
  if (typeof CUBES !== 'undefined') { CUBES = {}; EL = null; }
  showMode('single');
  history.replaceState(null, '', '?q=' + encodeURIComponent(q));
  document.title = `${R.common_name || R.iupac_name || R.formula} · Glass`;
  renderIdentity(); renderNames();
  const smi = R.smiles_canonical;
  $('#nmrdb').href = 'https://www.nmrdb.org/new_predictor/index.shtml?v=v2.157.0&smiles=' + encodeURIComponent(smi);
  $('#struct').innerHTML = '<div class="hint" style="padding:30px"><span class="spin"></span> drawing…</div>';
  $('#props').innerHTML = '<span class="hint"><span class="spin"></span> computing…</span>';
  api('/api/structure', { smiles: smi }).then(renderStructure).catch(e => $('#struct').textContent = e.message);
  api('/api/properties', { smiles: smi }).then(p => { PROPS = p; renderProps(p); if (typeof renderExpSection === 'function') renderExpSection(); }).catch(e => $('#props').textContent = e.message);
  loadSpectra();
  $('#exp').innerHTML = R.cid ? '<span class="spin"></span> loading…' : '<span class="hint">not in PubChem / offline</span>';
  if (R.cid) api('/api/experimental', { cid: R.cid }).then(renderExp).catch(() => $('#exp').textContent = 'unavailable');
  window.scrollTo({ top: 0 });
}
async function showError(q, msg) {
  const e = $('#error'); e.innerHTML = `<b>Could not resolve “${esc(q)}”.</b><div class="hint" style="margin-top:4px">${esc(msg)}</div><div class="sugg"></div>`; e.style.display = 'block';
  if (!$('#usepc').checked) return;
  try {
    const s = await api('/api/suggest', { q });
    if (s.suggestions.length) e.querySelector('.sugg').innerHTML = '<div class="hint" style="margin-top:8px">Did you mean:</div>' + s.suggestions.map(x => `<button class="small" data-q="${esc(x)}">${esc(x)}</button>`).join('');
    e.querySelectorAll('button[data-q]').forEach(b => b.onclick = () => { $('#q').value = b.dataset.q; resolve(); });
  } catch {}
}

function renderIdentity() {
  const k = $('#identity');
  const nv = R.names.filter(n => n.verified).length;
  const title = R.common_name || R.iupac_name || R.formula;
  k.innerHTML = `<div class="hint">${esc(R.input_kind)} · <span class="mono">${esc(R.input)}</span></div>
    <div class="name">${esc(title)}</div>
    <div class="iupac">${R.iupac_name && R.iupac_name !== title ? esc(R.iupac_name) + copyBtn(R.iupac_name) : (R.iupac_name ? '' : '<span class="hint">no verified systematic name yet; this structure is not in PubChem' + (R.generated_note ? '; ' + esc(R.generated_note) : '') + '. Type a candidate name in the checker below and it will be verified against the structure.</span>')}</div>
    <div class="chips-row"><span class="pill"><b>formula</b>${esc(R.formula)}</span>${R.cas ? `<span class="pill"><b>CAS</b>${esc(R.cas)}${copyBtn(R.cas)}</span>` : ''}${R.cid ? `<span class="pill"><b>PubChem</b><a href="https://pubchem.ncbi.nlm.nih.gov/compound/${R.cid}" target="_blank" rel="noopener">CID ${R.cid} ↗</a></span>` : ''}<span class="pill"><b>names</b>${R.names.length} · ${nv} verified</span></div>
    <div class="idrow">
      <b>Canonical SMILES</b><span class="mono">${esc(R.smiles_canonical)}${copyBtn(R.smiles_canonical)}</span>
      <b>IUPAC name</b><span>${R.iupac_name ? esc(R.iupac_name) + copyBtn(R.iupac_name) + ' <span class="badge ok">verified</span> <span class="hint">' + (R.generated_note ? esc(R.generated_note) : 'systematic name; not guaranteed to be the preferred IUPAC name') + '</span>' : '<span class="hint">·</span>'}</span>
      <b>Common name</b><span>${R.common_name ? esc(R.common_name) + copyBtn(R.common_name) : '<span class="hint">·</span>'}</span>
      <b>InChIKey</b><span class="mono">${esc(R.inchikey)}${copyBtn(R.inchikey)}</span>
    </div>`;
  $('#warnings').innerHTML = R.warnings.map(w => `<div class="warnbox">⚠ ${esc(w)}</div>`).join('');
  const rows = [
    ['Canonical SMILES (isomeric)', R.smiles_canonical], ['Canonical SMILES (no stereo)', R.smiles_canonical_nostereo], ['Kekulé SMILES', R.smiles_kekule],
    ['SMILES as entered', R.smiles_input], ['SMILES from OPSIN', R.smiles_opsin], ['SMILES from PubChem', R.smiles_pubchem], ['InChI', R.inchi], ['Standard InChIKey', R.inchikey], ['CAS RN', R.cas],
  ].filter(r => r[1]);
  $('#smiles').innerHTML = rows.map(([l, v]) => `<tr><th style="width:220px;position:static;text-transform:none;letter-spacing:0;font-size:12.5px">${l}</th><td class="mono">${esc(v)}${copyBtn(v)}</td></tr>`).join('');
  $('#thumb').innerHTML = '<span class="spin"></span>';
}

function renderNames() {
  const t = $('#names'); const f = $('#namefilter').value.trim().toLowerCase();
  const showUnv = $('#showunv').checked;
  const hidden = R.names.filter(n => !n.verified).length;
  const kindLabel = { input: 'your input', iupac: 'IUPAC (PubChem)', title: 'common name', synonym: 'synonym', generated: 'generated (SMILES→IUPAC model) · OPSIN-verified', 'opsin-generated': 'OPSIN' };
  t.innerHTML = `<tr><th style="width:90px">status</th><th>name</th><th>kind</th><th>numbering</th><th>note</th></tr>` + R.names.map((n, i) => {
    if (f && !n.name.toLowerCase().includes(f)) return '';
    if (!showUnv && !n.verified) return '';
    const b = n.verified ? '<span class="badge ok">verified</span>' : (n.note.includes('DIFFERENT') ? '<span class="badge bad">mismatch</span>' : '<span class="badge warn">unverified</span>');
    const num = hasNum(n) ? '<span class="badge acc">available</span>' : '';
    return `<tr data-i="${i}"><td>${b}</td><td>${esc(n.name)}${copyBtn(n.name)}</td><td class="muted">${esc(kindLabel[n.kind] || n.kind)}</td><td>${num}</td><td class="hint">${esc(n.note)}</td></tr>`;
  }).join('') + (!showUnv && hidden ? `<tr><td colspan="5" class="hint">${hidden} unverified database synonym(s) hidden; tick “show unverified names” to display them</td></tr>` : '');
  t.querySelectorAll('tr[data-i]').forEach(tr => tr.onclick = () => selectName(+tr.dataset.i));
  if (!SEL || !R.names.includes(SEL)) { const first = R.names.findIndex(hasNum); selectName(first >= 0 ? first : (R.names.length ? 0 : -1)); }
  else $$('#names tr').forEach(tr => tr.classList.toggle('sel', R.names[+tr.dataset.i] === SEL));
}
$('#namefilter').oninput = () => R && renderNames();
$('#showunv').onchange = () => R && renderNames();

const hasNum = n => !!(n && n.fragments && n.fragments.includes('parent'));
function selectName(i) {
  SEL = i >= 0 ? R.names[i] : null;
  $$('#names tr').forEach(tr => tr.classList.toggle('sel', +tr.dataset.i === i));
  renderNameBox(); drawOverlay();
  if (SPEC) renderSpectraTables();
  if (typeof renderGeomTables === 'function') renderGeomTables();
}
function renderNameBox() {
  const box = $('#namebox');
  if (!SEL) { box.innerHTML = '<span class="hint">no name selected</span>'; return; }
  const name = SEL.name; const spans = (SEL.spans && hoverAtom >= 0 && SEL.spans[String(hoverAtom)]) || [];
  let out = '', pos = 0;
  for (const [s, e] of spans.sort((a, b) => a[0] - b[0])) { out += esc(name.slice(pos, s)) + '<mark>' + esc(name.slice(s, e)) + '</mark>'; pos = e; }
  out += esc(name.slice(pos));
  box.innerHTML = out + (hasNum(SEL) ? '' : ' <span class="hint">(no locants: not a systematic name)</span>');
}

// ---------------- pane rail: one working pane at a time (structure default)
let PANE = 'structure';
function showPane(n) {
  PANE = n;
  $$('#paneRail button').forEach(b => b.classList.toggle('on', b.dataset.pane === n));
  $('#sec-structure').style.display = n === 'structure' ? '' : 'none';
  $('#sec-props').style.display = n === 'props' ? '' : 'none';
  $('#sec-exp').style.display = n === 'props' ? '' : 'none';
  const sp = ['nmr', 'ir', 'uv', 'pl'].includes(n);
  $('#sec-spectra').style.display = sp ? '' : 'none';
  const ctl = $('#specctl'); if (ctl) ctl.style.display = sp ? '' : 'none';
  $$('#specctl .nmronly').forEach(el => el.style.display = n === 'nmr' ? '' : 'none');
  if (sp) {
    const ic = { nmr: 'δ', ir: 'ν̃', uv: 'λ', pl: '✦' };
    const el = $('#specIc'); if (el) el.textContent = ic[n] || '∿';
    const map = { nmr: ['c_h1', 'c_c13', 'expsec'], ir: ['c_ir'], uv: ['c_uv'], pl: ['c_pl'] };
    $$('#sec-spectra .spec').forEach(d => d.style.display = 'none');
    (map[n] || []).forEach(id => { const el = document.getElementById(id); if (el) { const sd = el.closest('.spec'); if (sd) sd.style.display = ''; } });
    if (typeof redrawPane === 'function') requestAnimationFrame(() => redrawPane(n));
  }
}
$$('#paneRail button').forEach(b => b.onclick = () => showPane(b.dataset.pane));

// landing feature words: hover a word, its description appears
const FEATS = {
  verify: `A name is marked <span class="badge ok">verified</span> only if it parses back to this exact structure (identical Standard InChIKey).<br>Anything that cannot be verified stays clearly flagged.`,
  numbering: `Atom numbers come from the selected IUPAC name itself.<br>Hover an atom to see its locant; pick a different name to renumber.`,
  theory: `Identity is exact. Spectra and properties are computed estimates, each shown with its own error bar.<br>For assignment support and sanity checks, not measured data.`,
  docs: `Glass is free software under the AGPL-3.0-or-later; <a href="/api/about" target="_blank">credits, licences and source</a>. Built on open-source chemistry software and public data; the NMR model is derived from nmrshiftdb2 and shared under its database licence.<br>Landing art: the HOMO of pentacene, computed in Glass.`,
};
$$('.featwords button').forEach(b => {
  const show = () => {
    $$('.featwords button').forEach(x => x.classList.toggle('on', x === b));
    const bl = $('#featblurb'); bl.innerHTML = FEATS[b.dataset.key] || ''; bl.style.visibility = 'visible';
  };
  b.onmouseenter = show; b.onfocus = show; b.onclick = show;
});
$('#featwrap').onmouseleave = () => {
  $$('.featwords button').forEach(x => x.classList.remove('on'));
  $('#featblurb').style.visibility = 'hidden';
};
