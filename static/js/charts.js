function chTheme() { const cs = getComputedStyle(document.documentElement); const v = n => cs.getPropertyValue(n).trim(); return { axis: v('--ch-axis') || '#aab3bd', text: v('--ch-text') || '#667380', faint: v('--ch-faint') || '#c9c5b9', curve: v('--ch-curve') || '#0d9488' }; }
// Canvas chart engine: zoom/pan, peak hover & selection linked to tables and structures; static chart helper.
// Load order matters only for the top-level statements; all files share the page's global scope.

// ---------------- interactive chart engine
function nmrChart(cv, table, xs, ys, peaks, o) {
  const st = { xs, ys, peaks, o, lo: Math.min(...xs), hi: Math.max(...xs), hover: -1, sel: -1, drag: null };
  const rowSel = o.rowSelector || 'tr.pk'; st.full = [st.lo, st.hi];
  const L = 44, Rr = 10, T = 18, B = 30;
  const rev = o.reverseX !== false;
  const geom = () => { const W = cv.clientWidth || 600, H = cv.clientHeight || 300; return { W, H, X: x => rev ? L + (st.hi - x) / (st.hi - st.lo) * (W - L - Rr) : L + (x - st.lo) / (st.hi - st.lo) * (W - L - Rr), invX: px => rev ? st.hi - (px - L) / (W - L - Rr) * (st.hi - st.lo) : st.lo + (px - L) / (W - L - Rr) * (st.hi - st.lo) }; };
  const rows = [...table.querySelectorAll(rowSel)];
  const setHover = k => { st.hover = k; rows.forEach(r => r.classList.toggle('hl', +r.dataset.k === k)); draw(); const p = peaks[k]; if (o.onHover) o.onHover(k); else if (p) highlightAtoms(p.atoms, o.kind); else highlightAtoms([], o.kind); };
  const setSel = k => { st.sel = st.sel === k ? -1 : k; if (o.onSelect) { o.onSelect(k); } else rows.forEach(r => r.classList.toggle('sel', +r.dataset.k === st.sel)); draw(); };
  cv._setSel = k => { st.sel = k; draw(); };
  cv._setHover = k => { st.hover = k; rows.forEach(r => r.classList.toggle('hl', +r.dataset.k === k)); draw(); };
  rows.forEach(r => { r.onmouseenter = () => setHover(+r.dataset.k); r.onmouseleave = () => setHover(-1); r.onclick = () => { setSel(+r.dataset.k); const p = peaks[+r.dataset.k]; const w = (st.hi - st.lo); if (p.x < st.lo || p.x > st.hi) { st.lo = p.x - w / 2; st.hi = p.x + w / 2; draw(); } }; });
  function draw() {
    const dpr = devicePixelRatio || 1; const { W, H, X } = geom();
    cv.width = W * dpr; cv.height = H * dpr; const g = cv.getContext('2d'); g.scale(dpr, dpr); g.clearRect(0, 0, W, H);
    let i0 = 0, i1 = xs.length - 1; while (i0 < i1 && xs[i0 + 1] < st.lo) i0++; while (i1 > i0 && xs[i1 - 1] > st.hi) i1--;
    let ymax = 1e-9; for (let i = i0; i <= i1; i++) if (ys[i] > ymax) ymax = ys[i];
    if (o.ymax != null) ymax = o.ymax / 1.08;
    const Y = y => T + (1 - y / (ymax * 1.08)) * (H - T - B);
    g.strokeStyle = chTheme().axis; g.lineWidth = 1; g.beginPath(); g.moveTo(L, T); g.lineTo(L, H - B); g.lineTo(W - Rr, H - B); g.stroke();
    g.fillStyle = chTheme().text; g.font = '10px sans-serif'; g.textAlign = 'center';
    const span = st.hi - st.lo; const step = niceStep(span / 8);
    for (let x = Math.ceil(st.lo / step) * step; x <= st.hi + 1e-9; x += step) { g.fillText(+x.toFixed(3), X(x), H - B + 12); g.beginPath(); g.moveTo(X(x), H - B); g.lineTo(X(x), H - B + 3); g.stroke(); }
    g.fillText(o.xlabel || 'δ / ppm', (L + W - Rr) / 2, H - 4);
    if (o.ylabel) { g.save(); g.translate(10, (T + H - B) / 2); g.rotate(-Math.PI / 2); g.fillText(o.ylabel, 0, 0); g.restore(); }
    if (o.vlines) for (const v of o.vlines) { if (v.x < st.lo || v.x > st.hi) continue; g.strokeStyle = v.color; g.setLineDash([5, 4]); g.lineWidth = 1.5; g.beginPath(); g.moveTo(X(v.x), T); g.lineTo(X(v.x), H - B); g.stroke(); g.setLineDash([]); g.fillStyle = v.color; g.font = 'bold 10px sans-serif'; g.fillText(v.label, X(v.x), T + 46); }
    if (o.vline != null && o.vline >= st.lo && o.vline <= st.hi) { g.strokeStyle = chTheme().faint; g.setLineDash([3, 3]); g.beginPath(); g.moveTo(X(o.vline), T); g.lineTo(X(o.vline), H - B); g.stroke(); g.setLineDash([]); g.fillStyle = chTheme().faint; g.fillText('solvent', X(o.vline), T + 8); }
    peaks.forEach((p, k) => { if (k !== st.hover && k !== st.sel) return; g.fillStyle = k === st.sel ? 'rgba(255,200,0,.28)' : 'rgba(255,224,138,.35)'; const hw = Math.max(6, (W - L - Rr) * 0.004); g.fillRect(X(p.x) - hw, T, 2 * hw, H - T - B); });
    if (o.fillWavelength) {   // visible-spectrum band under the plot
      for (let i = i0; i < i1; i++) { if (xs[i] < 380 || xs[i] > 750) continue; g.fillStyle = wlColor(xs[i]); g.globalAlpha = 0.18; g.fillRect(X(xs[i]), H - B - 8, Math.max(1, X(xs[i + 1]) - X(xs[i]) + 0.5), 8); }
      g.globalAlpha = 1;
    }
    const series = [[ys, chTheme().curve]]; if (o.extra) series.push([o.extra, o.extraColor || '#0f766e']);
    for (const [yy, col] of series) {
      if (o.fillWavelength) { g.beginPath(); g.moveTo(X(xs[Math.max(0, i0 - 1)]), Y(0)); for (let i = Math.max(0, i0 - 1); i <= Math.min(xs.length - 1, i1 + 1); i++) g.lineTo(X(xs[i]), Y(yy[i])); g.lineTo(X(xs[Math.min(xs.length - 1, i1 + 1)]), Y(0)); g.closePath(); g.fillStyle = col; g.globalAlpha = 0.10; g.fill(); g.globalAlpha = 1; }
      g.strokeStyle = col; g.lineWidth = 1.4; g.beginPath();
      for (let i = Math.max(0, i0 - 1); i <= Math.min(xs.length - 1, i1 + 1); i++) { const px = X(xs[i]), py = Y(yy[i]); i === Math.max(0, i0 - 1) ? g.moveTo(px, py) : g.lineTo(px, py); } g.stroke();
    }
    if (o.extra) { g.font = '10px sans-serif'; g.textAlign = 'right'; g.fillStyle = '#4f46e5'; g.fillText('▬ absorption', W - Rr - 4, T + 24); g.fillStyle = '#0f766e'; g.fillText('▬ emission', W - Rr - 4, T + 36); g.textAlign = 'center'; }
    if (o.extra) peaks.forEach((p, k) => { if (p.x < st.lo || p.x > st.hi) return; const on = k === st.hover || k === st.sel; g.strokeStyle = on ? '#f59e0b' : (p.kind === 'abs' ? 'rgba(79,70,229,.45)' : 'rgba(15,118,110,.45)'); g.lineWidth = on ? 3 : 1; const px = X(p.x); g.beginPath(); g.moveTo(px, H - B); g.lineTo(px, Y(p.h / Math.max(...peaks.filter(q => q.kind === p.kind).map(q => q.h)))); g.stroke(); });
    if (o.sticks) peaks.forEach((p, k) => { if (p.x < st.lo || p.x > st.hi) return; const on = k === st.hover || k === st.sel; g.strokeStyle = on ? '#ff7f00' : 'rgba(196,55,43,.7)'; g.lineWidth = on ? 3 : 1.2; const px = X(p.x); g.beginPath(); g.moveTo(px, H - B); g.lineTo(px, H - B - p.h * (H - T - B) * 0.9); g.stroke(); });
    g.font = '9px sans-serif'; let lastx = -1e9;
    peaks.slice().sort((a, b) => b.x - a.x).forEach(p => { if (p.x < st.lo || p.x > st.hi) return; const px = X(p.x); const k = peaks.indexOf(p); const on = k === st.hover || k === st.sel; if (o.sticks && !on && (p.h || 0) < 0.15) return; if (!on && Math.abs(px - lastx) < 24) return; lastx = px; g.fillStyle = on ? '#b7791f' : '#8a6d1f'; g.font = on ? 'bold 10px sans-serif' : '9px sans-serif'; g.fillText(p.label, px, T + 8); });
    if (st.drag) { g.fillStyle = 'rgba(31,94,255,.12)'; g.fillRect(Math.min(st.drag[0], st.drag[1]), T, Math.abs(st.drag[1] - st.drag[0]), H - T - B); }
  }
  function nearest(px) { const { X } = geom(); let best = -1, bd = 14; peaks.forEach((p, k) => { const d = Math.abs(X(p.x) - px); if (d < bd) { bd = d; best = k; } }); return best; }
  cv.onmousemove = e => { const px = e.offsetX; if (st.drag) { st.drag[1] = px; draw(); return; } const k = nearest(px); if (k !== st.hover) setHover(k); };
  cv.onmouseleave = () => { if (st.hover !== -1) setHover(-1); };
  cv.onmousedown = e => { st.drag = [e.offsetX, e.offsetX]; };
  if (cv._mu) window.removeEventListener('mouseup', cv._mu);
  cv._mu = () => { if (!st.drag) return; const [a, b] = st.drag; st.drag = null; const { invX } = geom(); if (Math.abs(a - b) < 4) { const k = nearest(a); if (k >= 0) setSel(k); } else { const x1 = invX(Math.min(a, b)), x2 = invX(Math.max(a, b)); st.lo = Math.min(x1, x2); st.hi = Math.max(x1, x2); } draw(); };
  window.addEventListener('mouseup', cv._mu);
  cv.ondblclick = () => { [st.lo, st.hi] = st.full; draw(); };
  cv.onwheel = e => { if (!(e.ctrlKey || e.metaKey || e.shiftKey)) return; e.preventDefault(); const { invX } = geom(); const cx = invX(e.offsetX); const f = e.deltaY > 0 ? 1.25 : 0.8; let lo = cx - (cx - st.lo) * f, hi = cx + (st.hi - cx) * f; lo = Math.max(st.full[0], lo); hi = Math.min(st.full[1], hi); if (hi - lo > 0.02) { st.lo = lo; st.hi = hi; draw(); } };
  cv.onselectstart = () => false;
  draw();
  if (!cv._ro) { cv._ro = new ResizeObserver(() => draw()); cv._ro.observe(cv); }
}
function niceStep(raw) { const p = Math.pow(10, Math.floor(Math.log10(raw))); const m = raw / p; return (m < 1.5 ? 1 : m < 3 ? 2 : m < 7 ? 5 : 10) * p; }
function highlightAtoms(atoms, kind) {
  if (!DEP) return;
  const heavy = kind === 'h' ? [] : atoms;
  $$('#overlay .hit, #struct2 .hit').forEach(c => c.classList.toggle('on', heavy.includes(+c.dataset.i)));
  if (DEPH) { const hs = kind === 'h' ? DEPH.h_parent.map((p, i) => atoms.includes(p) && DEPH.is_h[i] ? i : -1).filter(i => i >= 0) : []; $$('#structH .hit').forEach(c => c.classList.toggle('on', hs.includes(+c.dataset.i))); }
}
function chart(cv, xs, ys, o) {
  const dpr = devicePixelRatio || 1; const W = cv.clientWidth || 500, H = cv.clientHeight || 230;
  cv.width = W * dpr; cv.height = H * dpr; const g = cv.getContext('2d'); g.scale(dpr, dpr);
  const L = 44, Rr = 10, T = 14, B = 30; g.clearRect(0, 0, W, H);
  const xmin = Math.min(...xs), xmax = Math.max(...xs); const ymax = o.ymax ?? Math.max(...ys, 1e-9) * 1.08, ymin = o.ymin ?? 0;
  const X = x => o.reverse ? L + (xmax - x) / (xmax - xmin) * (W - L - Rr) : L + (x - xmin) / (xmax - xmin) * (W - L - Rr);
  const Y = y => T + (1 - (y - ymin) / (ymax - ymin)) * (H - T - B);
  g.strokeStyle = chTheme().axis; g.lineWidth = 1; g.beginPath(); g.moveTo(L, T); g.lineTo(L, H - B); g.lineTo(W - Rr, H - B); g.stroke();
  g.fillStyle = chTheme().text; g.font = '10px sans-serif'; g.textAlign = 'center';
  const nt = 8; for (let i = 0; i <= nt; i++) { const x = xmin + (xmax - xmin) * i / nt; g.fillText(+x.toFixed(x < 20 ? 1 : 0), X(x), H - B + 12); }
  g.fillText(o.xlabel || '', (L + W - Rr) / 2, H - 4);
  if (o.ylabel) { g.save(); g.translate(10, (T + H - B) / 2); g.rotate(-Math.PI / 2); g.fillText(o.ylabel, 0, 0); g.restore(); }
  g.strokeStyle = chTheme().curve; g.lineWidth = 1.3; g.beginPath();
  xs.forEach((x, i) => { const px = X(x), py = Y(ys[i]); i ? g.lineTo(px, py) : g.moveTo(px, py); }); g.stroke();
  if (o.marks) { g.fillStyle = '#b7791f'; g.font = '9px sans-serif'; let lastx = -1e9; for (const [x, lab] of o.marks.sort((a, b) => a[0] - b[0])) { const px = X(x); if (Math.abs(px - lastx) < 22) continue; lastx = px; g.fillText(lab, px, T + 8); } }
}
