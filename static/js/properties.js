// Predicted-property tiles and grouped list, experimental data.
// Load order matters only for the top-level statements; all files share the page's global scope.

// ---------------- properties
function propRows(p) {
  const f = o => o && o.value != null ? `${o.value}${o.unit ? ' ' + o.unit : ''}${o.error ? ` (${o.error})` : ''}` : '·';
  const rows = [
    ['Molecular formula', p.formula], ['Molar mass', f(p.molar_mass)], ['Exact mass', f(p.exact_mass)],
    ['Melting point', p.melting_point ? `${p.melting_point.C} °C (${p.melting_point.K} K) · ${p.melting_point.method}, ${p.melting_point.error}` : '·'],
    ['Boiling point', p.boiling_point ? `${p.boiling_point.C} °C (${p.boiling_point.K} K) · ${p.boiling_point.error}` : '·'],
    ['Density', p.density.value != null ? `${p.density.value} g/cm³ · ${p.density.method}; ${p.density.error}` : '·'],
    ['Molar volume', f(p.molar_volume)], ['log P', f(p.logP)], ['TPSA', f(p.TPSA)], ['Molar refractivity', f(p.molar_refractivity)],
    ['H-bond donors / acceptors', `${p.h_bond_donors} / ${p.h_bond_acceptors}`], ['Rotatable bonds', p.rotatable_bonds],
    ['Rings (aromatic)', `${p.rings} (${p.aromatic_rings})`], ['Fraction sp³ C', p.fraction_sp3], ['Stereocentres', p.stereocentres], ['QED drug-likeness', p.QED_druglikeness ?? '·'], ['Formal charge', p.charge],
  ];
  const jb = p.joback;
  if (jb && jb.range_warning) rows.push(['Note', jb.range_warning]);
  if (jb && jb.error) rows.push(['Note', `${jb.error}; Joback not applicable (salts, ions, unusual groups)`]);
  if (jb && !jb.error) for (const k of ['Tc', 'Pc', 'Hvap', 'Hfus', 'Hf']) if (jb[k]) rows.push([{ Tc: 'Critical T', Pc: 'Critical P', Hvap: 'ΔH vap (Tb)', Hfus: 'ΔH fus', Hf: 'ΔH f° (gas)' }[k], f(jb[k])]);
  return rows;
}
function renderProps(p) {
  const t = [];
  t.push(['Molar mass', p.molar_mass.value.toFixed(2), 'g/mol', '']);
  if (p.melting_point) t.push(['Melting pt', p.melting_point.C.toFixed(0), '°C', 'Joback ±25–50 K']);
  if (p.boiling_point) t.push(['Boiling pt', p.boiling_point.C.toFixed(0), '°C', 'Joback ±15–25 K']);
  if (p.density.value != null) t.push(['Density', p.density.value.toFixed(2), 'g/cm³', 'Girolami ±0.1']);
  t.push(['log P', p.logP.value.toFixed(2), '', 'Crippen ±0.7']);
  t.push(['TPSA', p.TPSA.value.toFixed(0), 'Å²', '']);
  $('#tiles').innerHTML = t.map(([k, v, u, e]) => `<div class="tile"><div class="k">${k}</div><div class="v">${v} <span class="u">${u}</span></div>${e ? `<div class="e">${e}</div>` : ''}</div>`).join('');
  const groups = { 'Composition': ['Molecular formula', 'Molar mass', 'Exact mass', 'Formal charge', 'Stereocentres'], 'Phase & thermodynamics': ['Melting point', 'Boiling point', 'Density', 'Molar volume', 'Critical T', 'Critical P', 'ΔH vap (Tb)', 'ΔH fus', 'ΔH f° (gas)'], 'Polarity & shape': ['log P', 'TPSA', 'Molar refractivity', 'H-bond donors / acceptors', 'Rotatable bonds', 'Rings (aromatic)', 'Fraction sp³ C', 'QED drug-likeness'], 'Notes': ['Note'] };
  const rows = propRows(p); let html = '';
  for (const [g, keys] of Object.entries(groups)) {
    const rs = rows.filter(r => keys.includes(r[0])); if (!rs.length) continue;
    html += `<div class="pgroup">${g}</div>` + rs.map(([k, v]) => { const m = String(v).match(/^(.*?)( · | \()(.*)$/); return `<b>${k}</b><span>${m ? esc(m[1]) + ` <span class="hint">${esc(m[2].trim() === '(' ? '(' + m[3] : m[3])}</span>` : esc(v)}</span>`; }).join('');
  }
  $('#props').innerHTML = html;
}
function renderExp(e) {
  const keys = { mp: 'Melting point', bp: 'Boiling point', density: 'Density', solubility: 'Solubility', logp: 'log P' };
  const rows = Object.entries(e).map(([k, v]) => `<b>${keys[k] || k}</b><span>${v.map(esc).join('<br>')}</span>`);
  $('#exp').innerHTML = rows.length ? `<div class="kv">${rows.join('')}</div>` : '<span class="hint">no experimental values on PubChem</span>';
  R.experimental = e;
}
