// Page wiring: starts the solver worker, draws frames, handles controls.
import { Renderer, FIELDS, lutGradientCSS } from './render.js';

const $ = (id) => document.getElementById(id);
const canvas = $('flow');
const svg = $('walls');
const plot = $('plot');
const probe = $('probe');
const renderer = new Renderer(canvas);

const state = {
  field: 'mach',
  mach: 6,
  phi: 0.6,
  running: true,
  grid: null,
  frame: null,
  diag: null,
  dirty: false,
  hover: null,
};

// ---------- worker ----------
const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
worker.onerror = (e) => {
  setStatus('bad', 'Simulation failed to start', 'Your browser could not start the solver. Try a recent Chrome, Edge, Firefox or Safari, and open the page through a web server (not file://).');
  console.error(e);
};
worker.onmessage = (e) => {
  const m = e.data;
  if (m.type === 'grid') {
    state.grid = m;
    renderer.setGrid(m);
    drawWalls(m.geom);
    drawStations(m.geom);
    state.frame = null;
  } else if (m.type === 'frame') {
    state.frame = m.fields;
    state.diag = m.diag;
    state.dirty = true;
  }
};
worker.postMessage({ type: 'set', mach: state.mach, phi: state.phi });
worker.postMessage({ type: 'init', nx: +$('grid').value });

// ---------- theme colors for the canvas ----------
function cssColor(name) {
  const probeEl = document.createElement('span');
  probeEl.style.color = `var(${name})`;
  document.body.appendChild(probeEl);
  const c = getComputedStyle(probeEl).color;
  probeEl.remove();
  const m = c.match(/[\d.]+/g) || [128, 128, 128];
  return [+m[0], +m[1], +m[2]];
}
function refreshTheme() {
  renderer.setWallColor(cssColor('--wall'));
  renderer.setSchlierenColors(cssColor('--plot-bg'), cssColor('--ink'));
  state.dirty = true;
}
refreshTheme();
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', refreshTheme);
new MutationObserver(refreshTheme).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

// ---------- static overlays ----------
const SVGNS = 'http://www.w3.org/2000/svg';
function el(name, attrs, parent) {
  const n = document.createElementNS(SVGNS, name);
  for (const k in attrs) n.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(n);
  return n;
}

function drawWalls(geom) {
  const S = 1000;
  const W = geom.Lx * S, H = geom.Ly * S;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.innerHTML = '';
  const X = (x) => (x * S).toFixed(1);
  const Y = (y) => ((geom.Ly - (y - geom.y0)) * S).toFixed(1);
  const defs = el('defs', {}, svg);
  const pat = el('pattern', { id: 'hatch', width: 26, height: 26, patternUnits: 'userSpaceOnUse', patternTransform: 'rotate(45)' }, defs);
  el('rect', { width: 26, height: 26, class: 'solid' }, pat);
  el('line', { x1: 0, y1: 0, x2: 0, y2: 26, class: 'hatch-line' }, pat);
  for (const s of geom.solids) {
    const pts = s.pts.map(([x, y]) => `${X(x)},${Y(y)}`).join(' ');
    el('polygon', { points: pts, fill: 'url(#hatch)' }, svg);
    el('polygon', { points: pts, class: 'edge' }, svg);
  }
  // fuel injector band
  const body = geom.solids[0].pts;
  const yAt = (x) => {
    for (let k = 1; k < body.length - 1; k++) {
      const [xa, ya] = body[k], [xb, yb] = body[k + 1];
      if (x >= xa && x <= xb) return ya + (yb - ya) * (x - xa) / (xb - xa);
    }
    return body[1][1];
  };
  const inj = geom.injector;
  const yb = yAt(inj.x0);
  el('rect', { x: X(inj.x0), y: Y(geom.yCowl), width: ((inj.x1 - inj.x0) * S).toFixed(1), height: ((geom.yCowl - yb) * S).toFixed(1), class: 'inj' }, svg);
}

function drawStations(geom) {
  const box = $('stations');
  box.innerHTML = '';
  for (const s of geom.sections) {
    const d = document.createElement('div');
    d.className = 'sec';
    d.style.left = `${(s.x0 / geom.Lx) * 100}%`;
    d.style.width = `${((s.x1 - s.x0) / geom.Lx) * 100}%`;
    d.textContent = s.name;
    box.appendChild(d);
  }
  const axis = $('axis');
  axis.innerHTML = '';
  for (let x = 0; x <= geom.Lx; x += 0.5) {
    const t = document.createElement('span');
    t.style.left = `${(x / geom.Lx) * 100}%`;
    t.textContent = x === geom.Lx ? `${x} m` : String(x);
    axis.appendChild(t);
  }
}

// ---------- readouts ----------
function setStatus(kind, label, detail) {
  $('status').dataset.state = kind;
  $('status-label').textContent = label;
  $('status-detail').textContent = detail;
}

const fmt1 = (v) => (Number.isFinite(v) ? v.toFixed(1) : '—');
const fmt2 = (v) => (Number.isFinite(v) ? v.toFixed(2) : '—');
const fmtK = (v) => (Number.isFinite(v) ? Math.round(v).toLocaleString('en-US') : '—');

function updateReadouts(d) {
  const ign = state.grid?.fuel.TignK ?? 850;
  const F = d.forceKN;
  const tEl = $('thrust');
  tEl.textContent = (F > 0 ? '+' : F < 0 ? '−' : '') + Math.abs(F).toFixed(2);
  tEl.parentElement.classList.toggle('pos', F > 0.02);
  tEl.parentElement.classList.toggle('neg', F < -0.02);
  const fill = $('thrust-fill');
  const w = Math.min(Math.abs(F) / 10, F >= 0 ? 0.7 : 0.3) * 100;
  fill.style.width = `${w}%`;
  fill.style.left = F >= 0 ? '30%' : `${30 - w}%`;
  fill.classList.toggle('neg', F < 0);

  $('iso').textContent = fmt2(d.isoM);
  $('tiso').textContent = fmtK(d.TisoK);
  $('tcomb').textContent = d.phi > 0 && d.eta > 0.02 ? fmtK(d.TcombK) : '—';
  $('eta').textContent = d.phi > 0 ? Math.round(d.eta * 100) : '—';
  $('exitm').textContent = fmt2(d.exitM);
  $('time').textContent = fmt1(d.tPhysMs);
  $('sps').textContent = d.running ? `${Math.round(d.stepsPerSec)} steps/s` : 'paused';

  if (d.settling) {
    setStatus('settling', 'Starting up', 'Flow is still sweeping through the engine. Readouts settle in a few seconds.');
  } else if (d.throatM < 1.1 || d.isoMin < 0.9) {
    setStatus('bad', 'Inlet unstarted', 'Pressure from the combustor pushed a shock out of the duct, so air spills around the engine. Lower the fuel, then restart the flow.');
  } else if (d.phi === 0) {
    setStatus('warn', 'Fuel off', 'The inlet works, but with no heat added the engine only makes drag.');
  } else if (d.eta < 0.15) {
    setStatus('warn', 'No ignition', `Air before the injector peaks at ${fmtK(d.TisoK)} K. Hydrogen needs about ${fmtK(ign)} K here. Fly faster to heat the air more.`);
  } else if (d.exitMin < 1) {
    setStatus('warn', 'Close to thermal choking', 'So much heat is added that parts of the combustor flow have dropped below Mach 1. More fuel could unstart the inlet.');
  } else if (F <= 0) {
    setStatus('warn', 'Burning, but net drag', 'Combustion is supersonic, but at this speed and fuel setting the heat added does not overcome the drag.');
  } else {
    setStatus('good', 'Supersonic combustion', `Burning ${Math.round(d.eta * 100)}% of the fuel with the flow supersonic throughout. Net thrust +${F.toFixed(2)} kN/m.`);
  }
}

function updateLegend(range) {
  const def = FIELDS[state.field];
  $('legend-label').textContent = def.label + (def.unit ? ` (${def.unit})` : '');
  const bar = $('legend-bar');
  const lo = $('legend-lo'), hi = $('legend-hi');
  if (def.lut) {
    bar.style.background = lutGradientCSS(def.lut);
    if (range) { lo.textContent = def.fmt(range[0]); hi.textContent = def.fmt(range[1]); }
  } else if (def.gray) {
    bar.style.background = 'linear-gradient(90deg, var(--plot-bg), var(--ink))';
    lo.textContent = 'uniform';
    hi.textContent = 'steep density change';
  } else {
    bar.style.background = 'linear-gradient(90deg, var(--plot-bg) 0%, rgb(20,150,170) 45%, rgb(20,150,170) 55%, rgb(255,110,20) 100%)';
    lo.textContent = 'unburnt H₂';
    hi.textContent = 'burning';
  }
}

function updateProbe() {
  const h = state.hover;
  if (!h || !state.frame) { probe.hidden = true; return; }
  const s = renderer.sample(state.frame, h.x, h.y);
  if (!s) { probe.hidden = true; return; }
  probe.hidden = false;
  if (s.wall) {
    probe.innerHTML = `<b>Wall</b>`;
  } else {
    probe.innerHTML =
      `<b>Mach</b> ${s.M.toFixed(2)}<br>` +
      `<b>T</b> ${fmtK(s.TK)} K<br>` +
      `<b>p</b> ${s.pKPa < 10 ? s.pKPa.toFixed(2) : s.pKPa.toFixed(1)} kPa<br>` +
      `<b>speed</b> ${fmtK(s.speed)} m/s` +
      (s.phiLocal > 0.01 ? `<br><b>fuel φ</b> ${s.phiLocal.toFixed(2)}${s.burning ? ' · burning' : ''}` : '');
  }
  const r = plot.getBoundingClientRect();
  const flip = h.px > r.width * 0.7;
  probe.style.left = `${h.px}px`;
  probe.style.top = `${h.py}px`;
  probe.style.transform = flip ? 'translate(calc(-100% - 12px), -50%)' : 'translate(12px, -50%)';
}

// ---------- draw loop ----------
let lastFrame = null;
function tick() {
  if (state.dirty && state.frame) {
    state.dirty = false;
    const range = renderer.draw(state.frame, state.field, state.diag);
    updateLegend(range);
    if (state.diag) updateReadouts(state.diag);
    updateProbe();
    if (lastFrame !== state.frame) {
      lastFrame = state.frame;
      worker.postMessage({ type: 'ack' });
    }
  }
  requestAnimationFrame(tick);
}
requestAnimationFrame(tick);

// ---------- controls ----------
document.querySelectorAll('.field-tabs button').forEach((b) => {
  b.addEventListener('click', () => {
    state.field = b.dataset.field;
    document.querySelectorAll('.field-tabs button').forEach((o) => o.setAttribute('aria-selected', String(o === b)));
    state.dirty = true;
  });
});

$('mach').addEventListener('input', (e) => {
  state.mach = +e.target.value;
  $('mach-out').textContent = state.mach.toFixed(1);
  $('lede-mach').textContent = `Mach ${state.mach.toFixed(1).replace(/\.0$/, '')}`;
  worker.postMessage({ type: 'set', mach: state.mach });
});
$('phi').addEventListener('input', (e) => {
  state.phi = +e.target.value;
  $('phi-out').textContent = state.phi.toFixed(2);
  worker.postMessage({ type: 'set', phi: state.phi });
});
$('run').addEventListener('click', () => {
  state.running = !state.running;
  $('run').textContent = state.running ? 'Pause' : 'Resume';
  worker.postMessage({ type: 'run', on: state.running });
});
$('reset').addEventListener('click', () => worker.postMessage({ type: 'reset' }));
$('grid').addEventListener('change', (e) => worker.postMessage({ type: 'init', nx: +e.target.value }));

plot.addEventListener('pointermove', (e) => {
  if (!state.grid) return;
  const r = plot.getBoundingClientRect();
  const fx = (e.clientX - r.left) / r.width, fy = (e.clientY - r.top) / r.height;
  const g = state.grid.geom;
  state.hover = { x: fx * g.Lx, y: g.y0 + (1 - fy) * g.Ly, px: e.clientX - r.left, py: e.clientY - r.top };
  updateProbe();
});
plot.addEventListener('pointerleave', () => { state.hover = null; probe.hidden = true; });

// restore controls if the browser kept form values
for (const id of ['mach', 'phi']) $(id).dispatchEvent(new Event('input'));
