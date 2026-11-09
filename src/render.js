// Dark-tunnel renderer: neon colormaps, flame bloom, glowing walls, particles.
import { toKelvin, toKPa, toMS } from './solver.js';

// ---------- colormaps (hand-tuned stops, interpolated into 256-entry LUTs) ----------
const STOPS = {
  // fast = bright: indigo -> electric blue -> cyan
  ice: ['#020308', '#140a33', '#30188f', '#3446e8', '#2380f0', '#1fb2f2', '#5fd8ff'],
  // cold = black, hot = white: wine -> crimson -> orange -> white-hot
  fire: ['#020203', '#2a0614', '#6e0b2a', '#c0182f', '#f2501b', '#ffa531', '#ffe08a', '#fffaf0'],
  // low = black: violet -> magenta -> amber
  neon: ['#030208', '#1e0b3d', '#4e1688', '#9b1fa8', '#e2337f', '#ff7a45', '#ffd38a'],
  // eddy viscosity: deep teal -> aqua -> lime
  aurora: ['#020507', '#062124', '#0a4f55', '#0e8f8a', '#2fd3a0', '#a6f56b', '#f2ffb0'],
};

function hex(h) { const n = parseInt(h.slice(1), 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
function buildLUT(stops) {
  const cs = stops.map(hex);
  const lut = new Uint8ClampedArray(256 * 3);
  for (let k = 0; k < 256; k++) {
    const t = (k / 255) * (cs.length - 1);
    const i = Math.min(cs.length - 2, Math.floor(t)), f = t - i;
    for (let ch = 0; ch < 3; ch++) lut[k * 3 + ch] = cs[i][ch] + (cs[i + 1][ch] - cs[i][ch]) * f;
  }
  return lut;
}
export const LUTS = Object.fromEntries(Object.entries(STOPS).map(([k, v]) => [k, buildLUT(v)]));
export function lutGradientCSS(name) {
  const s = STOPS[name];
  return `linear-gradient(90deg, ${s.map((c, i) => `${c} ${((i / (s.length - 1)) * 100).toFixed(0)}%`).join(', ')})`;
}

const TUNNEL_BG = [3, 4, 9];

// value below which a fraction q of the samples fall (histogram estimate)
const HIST = new Uint32Array(1024);
function percentile(vals, n, lo, hi, q) {
  if (!(hi > lo) || n === 0) return hi;
  HIST.fill(0);
  const sc = 1023 / (hi - lo);
  for (let k = 0; k < n; k++) HIST[Math.min(1023, Math.max(0, ((vals[k] - lo) * sc) | 0))]++;
  let acc = 0;
  const target = q * n;
  for (let b = 0; b < 1024; b++) { acc += HIST[b]; if (acc >= target) return lo + (b + 1) / sc; }
  return hi;
}
const WALL_FILL = [14, 18, 24];

// ---------- fields ----------
const mach = (f, c) => Math.hypot(f.u[c], f.v[c]) / Math.sqrt(f.ga[c] * f.p[c] / f.r[c]);
export const FIELDS = {
  temperature: {
    label: 'Static temperature', unit: 'K', lut: 'fire', glow: 0.45, curve: 0.75,
    value: (f, c) => toKelvin(f.th[c]),
    range: (d, lo, hi) => [200, Math.max(1200, hi)],
    fmt: (v) => Math.round(v).toLocaleString('en-US'),
  },
  mach: {
    label: 'Mach number', unit: '', lut: 'ice',
    value: mach,
    range: (d) => [0, Math.max(1, d.mach || 6)],
    fmt: (v) => v.toFixed(1),
  },
  pressure: {
    label: 'Pressure (log scale)', unit: 'kPa', lut: 'neon', log: true,
    value: (f, c) => toKPa(f.p[c]),
    range: (d, lo, hi) => [toKPa(1 / 1.4) * 0.5, Math.max(toKPa(1 / 1.4) * 10, hi)],
    fmt: (v) => (v < 10 ? v.toFixed(1) : Math.round(v).toLocaleString('en-US')),
  },
  turbulence: {
    label: 'Eddy viscosity / molecular viscosity (log)', unit: '', lut: 'aurora', log: true,
    value: (f, c) => 1 + f.mt[c],
    range: () => [1, 1000],
    fmt: (v) => (v < 2 ? '0' : Math.round(v).toLocaleString('en-US')),
  },
  schlieren: { label: 'Schlieren (density gradient)', unit: '', schlieren: true },
  fuel: { label: 'Fuel and flame', unit: '', fuel: true, glow: 1 },
};

// ---------- renderer ----------
export class Renderer {
  constructor(flowCanvas, fxCanvas, trailCanvas) {
    this.flow = flowCanvas;
    this.ctx = flowCanvas.getContext('2d', { alpha: false });
    this.fx = fxCanvas;
    this.fxCtx = fxCanvas.getContext('2d');
    this.trail = trailCanvas;
    this.trailCtx = trailCanvas.getContext('2d');
    this.hiSmooth = {};
    this.heat = document.createElement('canvas');
    this.heatCtx = this.heat.getContext('2d');
    this.bloomSmall = document.createElement('canvas');
    this.bloomCtx = this.bloomSmall.getContext('2d');
    this.particles = null;
    this.particlesOn = true;
  }

  setGrid(grid) {
    this.grid = grid;
    this.flow.width = grid.nx;
    this.flow.height = grid.ny;
    this.img = this.ctx.createImageData(grid.nx, grid.ny);
    this.val = new Float32Array(grid.NX * grid.NY);
    this.heat.width = grid.nx;
    this.heat.height = grid.ny;
    this.heatImg = this.heatCtx.createImageData(grid.nx, grid.ny);
    this.bloomSmall.width = Math.max(8, Math.round(grid.nx / 5));
    this.bloomSmall.height = Math.max(4, Math.round(grid.ny / 5));
    this.hiSmooth = {};
    this.particles = null;
    this.resize();
  }

  // match the overlay canvases to their on-screen size
  resize() {
    const r = this.fx.getBoundingClientRect();
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    const w = Math.max(1, Math.round(r.width * dpr)), h = Math.max(1, Math.round(r.height * dpr));
    for (const cv of [this.fx, this.trail]) if (cv.width !== w || cv.height !== h) { cv.width = w; cv.height = h; }
    this.dpr = dpr;
  }

  idxAt(x, y) {
    const g = this.grid;
    const i = 2 + Math.min(g.nx - 1, Math.max(0, Math.floor((x / g.geom.Lx) * g.nx)));
    const j = 2 + Math.min(g.ny - 1, Math.max(0, Math.floor(((y - g.geom.y0) / g.geom.Ly) * g.ny)));
    return j * g.NX + i;
  }

  // ---- main field ----
  draw(frame, fieldKey, diag) {
    const g = this.grid;
    if (!g || !frame) return null;
    const { nx, ny, NX, cellType } = g;
    const f = frame, px = this.img.data, val = this.val;
    const def = FIELDS[fieldKey];
    const PAD = 2;
    let range = null;

    if (def.schlieren || def.fuel) {
      let gmax = 1e-9;
      for (let j = PAD; j < ny + PAD; j++)
        for (let i = PAD; i < nx + PAD; i++) {
          const c = j * NX + i;
          if (cellType[c] !== 0) { val[c] = 0; continue; }
          const gx = f.r[c + 1] - f.r[c - 1], gy = f.r[c + NX] - f.r[c - NX];
          const v = Math.hypot(gx, gy) / f.r[c];
          val[c] = v;
          if (v > gmax) gmax = v;
        }
      const k = 14 / gmax;
      const zmax = g.fuel.fst * Math.max(0.3, diag?.phi ?? 1);
      for (let j = 0; j < ny; j++) {
        const row = (ny - 1 - j) * nx;
        for (let i = 0; i < nx; i++) {
          const c = (j + PAD) * NX + i + PAD, o = (row + i) * 4;
          if (cellType[c] !== 0) { px[o] = WALL_FILL[0]; px[o + 1] = WALL_FILL[1]; px[o + 2] = WALL_FILL[2]; px[o + 3] = 255; continue; }
          const s = 1 - Math.exp(-k * val[c]);   // 0 = uniform, 1 = sharp gradient
          let R, G, B;
          if (def.schlieren) {
            // cyan shocks glowing on black, the strongest ones burn white
            R = TUNNEL_BG[0] + s * s * 230; G = TUNNEL_BG[1] + s * 225; B = TUNNEL_BG[2] + s * 246;
          } else {
            R = TUNNEL_BG[0] + s * 34; G = TUNNEL_BG[1] + s * 46; B = TUNNEL_BG[2] + s * 70;
            // unburnt hydrogen: teal plumes (local fuel/air relative to the target)
            const zf = f.y[c] / (1 - f.z[c] + 1e-9);
            const a = Math.min(1, zf / zmax);
            R += (20 - R) * a * 0.9; G += (220 - G) * a * 0.9; B += (200 - B) * a * 0.9;
            // burnt products: faint warm haze
            const pb = Math.min(1, (f.z[c] - f.y[c]) / (zmax * 0.6));
            R += (120 - R) * pb * 0.35; G += (40 - G) * pb * 0.35; B += (60 - B) * pb * 0.35;
          }
          px[o] = R; px[o + 1] = G; px[o + 2] = B; px[o + 3] = 255;
        }
      }
    } else {
      let lo = Infinity, hi = -Infinity, n = 0;
      const samp = this.samp || (this.samp = new Float32Array(this.val.length));
      for (let j = PAD; j < ny + PAD; j++)
        for (let i = PAD; i < nx + PAD; i++) {
          const c = j * NX + i;
          if (cellType[c] !== 0) continue;
          const v = def.value(f, c);
          val[c] = v;
          samp[n++] = v;
          if (v < lo) lo = v;
          if (v > hi) hi = v;
        }
      // scale to the 99.5th percentile so one stagnant hot spot can't wash out the rest
      hi = percentile(samp, n, lo, hi, 0.995);
      const prev = this.hiSmooth[fieldKey];
      const hs = prev === undefined ? hi : prev + (hi - prev) * 0.15;
      this.hiSmooth[fieldKey] = hs;
      let [a, b] = def.range(diag || {}, lo, hs);
      range = [a, b];
      if (def.log) { a = Math.log(a); b = Math.log(b); }
      const lut = LUTS[def.lut];
      const span = 255 / Math.max(1e-9, b - a);
      const curve = def.curve || 1;
      for (let j = 0; j < ny; j++) {
        const row = (ny - 1 - j) * nx;
        for (let i = 0; i < nx; i++) {
          const c = (j + PAD) * NX + i + PAD, o = (row + i) * 4;
          if (cellType[c] !== 0) { px[o] = WALL_FILL[0]; px[o + 1] = WALL_FILL[1]; px[o + 2] = WALL_FILL[2]; px[o + 3] = 255; continue; }
          let v = val[c];
          if (def.log) v = Math.log(Math.max(1e-12, v));
          let t = (v - a) * span;
          if (curve !== 1 && t > 0) t = 255 * Math.pow(t / 255, curve);
          const kk = Math.max(0, Math.min(255, t | 0)) * 3;
          px[o] = lut[kk]; px[o + 1] = lut[kk + 1]; px[o + 2] = lut[kk + 2]; px[o + 3] = 255;
        }
      }
    }
    this.ctx.putImageData(this.img, 0, 0);
    this.drawFx(frame, def);
    return range;
  }

  // ---- glow layer: flame bloom + hot walls ----
  drawFx(frame, def) {
    const g = this.grid;
    this.resize();
    const ctx = this.fxCtx, W = this.fx.width, H = this.fx.height;
    ctx.clearRect(0, 0, W, H);
    ctx.globalCompositeOperation = 'lighter';

    if (def.glow) {
      const { nx, ny, NX } = g;
      const hr = frame.hr, d = this.heatImg.data;
      const hs = this.hsamp || (this.hsamp = new Float32Array(hr.length));
      let n = 0, hmax = 1e-9;
      for (let c = 0; c < hr.length; c++) if (hr[c] > 0) { hs[n++] = hr[c]; if (hr[c] > hmax) hmax = hr[c]; }
      hmax = Math.max(1e-9, percentile(hs, n, 0, hmax, 0.97));
      for (let j = 0; j < ny; j++) {
        const row = (ny - 1 - j) * nx;
        for (let i = 0; i < nx; i++) {
          const c = (j + 2) * NX + i + 2, o = (row + i) * 4;
          const h = Math.min(1, Math.sqrt(Math.max(0, hr[c]) / hmax));
          // orange at the edges of the flame, white-hot at its core
          d[o] = 255; d[o + 1] = 90 + 165 * h * h; d[o + 2] = 20 + 200 * h * h * h; d[o + 3] = 255 * Math.min(1, h * 1.4);
        }
      }
      this.heatCtx.putImageData(this.heatImg, 0, 0);
      const b = this.bloomCtx, bw = this.bloomSmall.width, bh = this.bloomSmall.height;
      b.clearRect(0, 0, bw, bh);
      b.imageSmoothingEnabled = true;
      b.drawImage(this.heat, 0, 0, bw, bh);
      ctx.imageSmoothingEnabled = true;
      const gs = def.glow;
      ctx.globalAlpha = gs;
      ctx.drawImage(this.bloomSmall, 0, 0, W, H);      // wide soft bloom
      ctx.drawImage(this.bloomSmall, 0, 0, W, H);
      ctx.globalAlpha = 0.75 * gs;
      ctx.drawImage(this.heat, 0, 0, W, H);            // sharper core
      ctx.globalAlpha = 1;
    }

    // hot walls: heat flux into the wall, dark red -> orange -> white
    const wq = frame.wq, wl = g.walls;
    if (wq && wq.length) {
      let qmax = 1e-9;
      for (let k = 0; k < wq.length; k++) if (wq[k] > qmax) qmax = wq[k];
      qmax = Math.max(1e-9, percentile(wq, wq.length, 0, qmax, 0.95));
      const sx = W / g.geom.Lx, sy = H / g.geom.Ly;
      const half = 0.75 * (g.geom.Lx / g.nx) * sx;
      ctx.lineCap = 'round';
      for (let pass = 0; pass < 2; pass++) {
        ctx.lineWidth = (pass === 0 ? 7 : 2.2) * this.dpr;
        for (let k = 0; k < wq.length; k++) {
          const q = Math.min(1, Math.max(0, wq[k]) / qmax);
          if (q < 0.04) continue;
          const x = wl.x[k] * sx, y = H - (wl.y[k] - g.geom.y0) * sy;
          const tx = -wl.ny[k], ty = wl.nx[k];
          const a = pass === 0 ? 0.18 * q : 0.85 * q;
          const gch = Math.round(60 + 190 * q * q), bch = Math.round(20 + 200 * q * q * q);
          ctx.strokeStyle = `rgba(255,${gch},${bch},${a.toFixed(3)})`;
          ctx.beginPath();
          ctx.moveTo(x - tx * half, y + ty * half);
          ctx.lineTo(x + tx * half, y - ty * half);
          ctx.stroke();
        }
      }
    }
    ctx.globalCompositeOperation = 'source-over';
  }

  // ---- particles: tracers advected through the latest velocity field ----
  setParticles(on) {
    this.particlesOn = on;
    if (!on) this.trailCtx.clearRect(0, 0, this.trail.width, this.trail.height);
  }

  spawn(p, k, anywhere) {
    const g = this.grid, L = g.geom;
    for (let tries = 0; tries < 20; tries++) {
      const x = anywhere ? Math.random() * L.Lx : Math.random() * 0.08 * L.Lx;
      const y = L.y0 + Math.random() * L.Ly;
      if (g.cellType[this.idxAt(x, y)] === 0) { p.x[k] = x; p.y[k] = y; p.age[k] = Math.random() * 40; return; }
    }
    p.x[k] = -1;
  }

  stepParticles(frame, dtSec, mach) {
    const g = this.grid;
    if (!g || !frame || !this.particlesOn) return;
    const n = Math.round(Math.min(1300, 800 * (this.fx.width / 1600) + 250));
    if (!this.particles || this.particles.x.length !== n) {
      this.particles = { x: new Float32Array(n), y: new Float32Array(n), age: new Float32Array(n) };
      for (let k = 0; k < n; k++) this.spawn(this.particles, k, true);
    }
    const p = this.particles, L = g.geom, W = this.trail.width, H = this.trail.height;
    const ctx = this.trailCtx;
    // fade old trails
    ctx.globalCompositeOperation = 'destination-out';
    ctx.fillStyle = 'rgba(0,0,0,0.22)';
    ctx.fillRect(0, 0, W, H);
    ctx.globalCompositeOperation = 'source-over';
    // freestream particles cross the 4 m engine in about 2.4 s
    const k = Math.min(0.05, dtSec) * L.Lx / (2.4 * Math.max(1, mach));
    const sx = W / L.Lx, sy = H / L.Ly;
    const f = frame;
    ctx.lineWidth = 1 * this.dpr;
    ctx.strokeStyle = 'rgba(190,232,255,0.38)';
    ctx.beginPath();
    for (let q = 0; q < p.x.length; q++) {
      if (p.x[q] < 0) { this.spawn(p, q, false); continue; }
      const c = this.idxAt(p.x[q], p.y[q]);
      const u = f.u[c], v = f.v[c];
      const x1 = p.x[q] + u * k, y1 = p.y[q] + v * k;
      p.age[q] += 1;
      const out = x1 < 0 || x1 >= L.Lx || y1 < L.y0 || y1 >= L.y0 + L.Ly;
      if (out || g.cellType[this.idxAt(x1, y1)] !== 0 || p.age[q] > 420 || !(Number.isFinite(x1))) {
        this.spawn(p, q, Math.random() < 0.35);
        continue;
      }
      ctx.moveTo(p.x[q] * sx, H - (p.y[q] - L.y0) * sy);
      ctx.lineTo(x1 * sx, H - (y1 - L.y0) * sy);
      p.x[q] = x1; p.y[q] = y1;
    }
    ctx.stroke();
  }

  // ---- probe ----
  sample(frame, x, y) {
    const g = this.grid;
    if (!g || !frame) return null;
    const c = this.idxAt(x, y);
    if (g.cellType[c] !== 0) return { wall: true };
    const f = frame;
    return {
      M: mach(f, c),
      TK: toKelvin(f.th[c]),
      pKPa: toKPa(f.p[c]),
      speed: toMS(Math.hypot(f.u[c], f.v[c])),
      phiLocal: (f.z[c] / Math.max(1e-9, 1 - f.z[c])) / g.fuel.fst,
      burnt: f.z[c] > 1e-6 ? 1 - f.y[c] / f.z[c] : 0,
      mt: f.mt[c],
    };
  }
}
