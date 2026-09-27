// Section navigation and start-up.
// Load order matters only for the top-level statements; all files share the page's global scope.

// ---------------- section nav highlight
const secs = ['sec-identity', 'sec-names'];
window.addEventListener('scroll', () => {
  if (!R) return; let cur = secs[0];
  for (const id of secs) { const el = document.getElementById(id); if (el && el.offsetParent && el.getBoundingClientRect().top < 140) cur = id; }
  $$('#subnav a').forEach(a => a.classList.toggle('on', a.getAttribute('href') === '#' + cur));
}, { passive: true });


// ---------------- init
$('#themebtn').onclick = () => {
  const cur = document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark';
  document.documentElement.dataset.theme = cur;
  try { localStorage.setItem('glass_theme', cur); } catch (e) {}
  if (typeof SPEC !== 'undefined' && SPEC) renderSpectraTables();     // charts re-read theme colors
  if (typeof drawIR === 'function') drawIR();
  if (typeof redrawPL === 'function' && typeof PL !== 'undefined' && PL) redrawPL();
};

(async () => {
  renderLanding(); showMode('single');
  const sols = await api('/api/solvents', {}); SOLVS = sols; $('#solvent').innerHTML = sols.map(s => `<option value="${s.name}">${s.name}${s.nmr ? ' (' + s.nmr + ')' : ''}</option>`).join('');
  METHODS = await api('/api/methods', {});
  $('#go').onclick = resolve; $('#q').addEventListener('keydown', e => { if (e.key === 'Enter') resolve(); });
  document.addEventListener('keydown', e => { if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === 'k') { e.preventDefault(); $('#q').focus(); $('#q').select(); } });
  const q = new URLSearchParams(location.search).get('q'); if (q) { $('#q').value = q; resolve(); }
})();
