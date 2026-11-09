// Turns a flow snapshot into pixels.
import { GAMMA, toKelvin, toKPa, toMS } from './solver.js';

const G = GAMMA;

// Polynomial fits of the matplotlib / Google colormaps (Matt Zucker, CC0)
const POLY = {
  turbo: [
    [0.1140890109226559, 0.06288340699912215, 0.2248337216805064],
    [6.716419496985708, 3.182286745507602, 7.571581586103393],
    [-66.09402360453038, -4.9279827041226, -10.09439367561635],
    [228.7660791526501, 25.04986699771073, -91.54105330182436],
    [-334.8351565777451, -69.31749712757485, 288.5858850615712],
    [218.7637218434795, 67.52150567819112, -305.2045772184957],
    [-52.88903478218835, -21.54527364654712, 110.5174647748972],
  ],
  inferno: [
    [0.0002189403691192265, 0.001651004631001012, -0.01948089843709184],
    [0.1065134194856116, 0.5639564367884091, 3.932712388889277],
    [11.60249308247187, -3.972853965665698, -15.9423941062914],
    [-41.70399613139459, 17.43639888205313, 44.35414519872813],
    [77.162935699427, -33.40235894210092, -81.80730925738993],
    [-71.31942824499214, 32.62606426397723, 73.20951985803202],
    [25.13112622477341, -12.24266895238567, -23.07032500287172],
  ],
  plasma: [
    [0.05873234392399702, 0.02333670892565664, 0.5433401826748754],
    [2.176514634195958, 0.2383834171260182, 0.7539604599784036],
    [-2.689460476458034, -7.455851135738909, 3.110799939717086],
    [6.130348345893603, 42.3461881477227, -28.51885465332158],
    [-11.10743619062271, -82.66631109428045, 60.13984767418263],
    [10.02306557647065, 71.41361770095349, -54.07218655560067],
    [-3.658713842777788, -22.93153465461149, 18.19190778539828],
  ],
};

function buildLUT(name) {
  const c = POLY[name];
  const lut = new Uint8ClampedArray(256 * 3);
  for (let k = 0; k < 256; k++) {
    const t = k / 255;
    for (let ch = 0; ch < 3; ch++) {
      let v = 0;
      for (let d = 6; d >= 0; d--) v = v * t + c[d][ch];
      lut[k * 3 + ch] = Math.round(255 * Math.max(0, Math.min(1, v)));
    }
  }
  return lut;
}

export const LUTS = { turbo: buildLUT('turbo'), inferno: buildLUT('inferno'), plasma: buildLUT('plasma') };

export function lutGradientCSS(name) {
  const lut = LUTS[name];
  const stops = [];
  for (let k = 0; k <= 10; k++) {
    const i = Math.round(k * 25.5);
    stops.push(`rgb(${lut[i * 3]},${lut[i * 3 + 1]},${lut[i * 3 + 2]}) ${k * 10}%`);
  }
  return `linear-gradient(90deg, ${stops.join(',')})`;
}

// Field definitions. Each returns the scalar to color, plus how to label it.
export const FIELDS = {
  mach: {
    label: 'Mach number', unit: '', lut: 'turbo',
    value: (f, c) => Math.hypot(f.u[c], f.v[c]) / Math.sqrt(G * f.p[c] / f.r[c]),
    range: (d) => [0, Math.max(1, d.mach)],
    fmt: (v) => v.toFixed(1),
  },
  temperature: {
    label: 'Static temperature', unit: 'K', lut: 'inferno',
    value: (f, c) => toKelvin(f.p[c] / f.r[c]),
    range: (d, lo, hi) => [200, Math.max(900, hi)],
    fmt: (v) => Math.round(v).toLocaleString('en-US'),
  },
  pressure: {
    label: 'Pressure (log scale)', unit: 'kPa', lut: 'plasma', log: true,
    value: (f, c) => toKPa(f.p[c]),
    range: (d, lo, hi) => [toKPa(1 / G) * 0.5, Math.max(toKPa(1 / G) * 10, hi)],
    fmt: (v) => (v < 10 ? v.toFixed(1) : Math.round(v).toLocaleString('en-US')),
  },
  schlieren: {
    label: 'Schlieren (density gradient)', unit: '', gray: true,
  },
  fuel: {
    label: 'Fuel and flame', unit: '', fuel: true,
  },
};

export class Renderer {
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d', { alpha: false });
    this.hiSmooth = {};
  }

  setGrid(grid) {
    this.grid = grid;
    this.canvas.width = grid.nx;
    this.canvas.height = grid.ny;
    this.img = this.ctx.createImageData(grid.nx, grid.ny);
    this.val = new Float32Array(grid.NX * grid.NY);
    this.hiSmooth = {};
  }

  // wall color from the page theme, as [r,g,b]
  setWallColor(rgb) { this.wall = rgb; }
  setSchlierenColors(bg, ink) { this.sBg = bg; this.sInk = ink; }

  draw(frame, fieldKey, diag) {
    const g = this.grid;
    if (!g || !frame) return null;
    const { nx, ny, NX, cellType } = g;
    const f = frame;
    const px = this.img.data;
    const wall = this.wall || [128, 128, 128];
    const PAD = 2;
    const def = FIELDS[fieldKey];

    if (def.gray || def.fuel) {
      // numerical schlieren: |grad rho| / rho, shown on an exponential scale
      const val = this.val;
      let gmax = 1e-9;
      for (let j = PAD; j < ny + PAD; j++)
        for (let i = PAD; i < nx + PAD; i++) {
          const c = j * NX + i;
          if (cellType[c] !== 0) { val[c] = 0; continue; }
          const gx = (f.r[c + 1] - f.r[c - 1]) * 0.5, gy = (f.r[c + NX] - f.r[c - NX]) * 0.5;
          const v = Math.hypot(gx, gy) / f.r[c];
          val[c] = v;
          if (v > gmax) gmax = v;
        }
      const k = 12 / gmax;
      const bg = this.sBg || [236, 238, 240], ink = this.sInk || [20, 24, 30];
      const Ymax = g.fuel.Ystoich * Math.max(0.2, diag?.phi ?? 1);
      let hmax = 1e-9;
      if (def.fuel) for (let c = 0; c < f.hr.length; c++) if (f.hr[c] > hmax) hmax = f.hr[c];
      for (let j = 0; j < ny; j++) {
        const row = (ny - 1 - j) * nx;
        for (let i = 0; i < nx; i++) {
          const c = (j + PAD) * NX + i + PAD;
          const o = (row + i) * 4;
          if (cellType[c] !== 0) { px[o] = wall[0]; px[o + 1] = wall[1]; px[o + 2] = wall[2]; px[o + 3] = 255; continue; }
          const s = Math.exp(-k * val[c]);
          let R = ink[0] + (bg[0] - ink[0]) * s, Gc = ink[1] + (bg[1] - ink[1]) * s, B = ink[2] + (bg[2] - ink[2]) * s;
          if (def.fuel) {
            // soften shocks, tint unburnt fuel teal, flame orange
            R = bg[0] + (R - bg[0]) * 0.45; Gc = bg[1] + (Gc - bg[1]) * 0.45; B = bg[2] + (B - bg[2]) * 0.45;
            const a = Math.min(1, f.y[c] / Ymax) * 0.85;
            R += (20 - R) * a; Gc += (150 - Gc) * a; B += (170 - B) * a;
            const h = Math.min(1, Math.pow(f.hr[c] / hmax, 0.5));
            R += (255 - R) * h; Gc += (110 - Gc) * h; B += (20 - B) * h;
          }
          px[o] = R; px[o + 1] = Gc; px[o + 2] = B; px[o + 3] = 255;
        }
      }
      this.ctx.putImageData(this.img, 0, 0);
      return null;
    }

    // colormapped scalar
    const val = this.val;
    let lo = Infinity, hi = -Infinity;
    for (let j = PAD; j < ny + PAD; j++)
      for (let i = PAD; i < nx + PAD; i++) {
        const c = j * NX + i;
        if (cellType[c] !== 0) continue;
        const v = def.value(f, c);
        val[c] = v;
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
    // smooth the auto range so the colors don't flicker
    const prev = this.hiSmooth[fieldKey];
    const hs = prev === undefined ? hi : prev + (hi - prev) * 0.15;
    this.hiSmooth[fieldKey] = hs;
    let [a, b] = def.range(diag || {}, lo, hs);
    if (def.log) { a = Math.log(a); b = Math.log(b); }
    const lut = LUTS[def.lut];
    const span = 255 / Math.max(1e-9, b - a);
    for (let j = 0; j < ny; j++) {
      const row = (ny - 1 - j) * nx;
      for (let i = 0; i < nx; i++) {
        const c = (j + PAD) * NX + i + PAD;
        const o = (row + i) * 4;
        if (cellType[c] !== 0) { px[o] = wall[0]; px[o + 1] = wall[1]; px[o + 2] = wall[2]; px[o + 3] = 255; continue; }
        let v = val[c];
        if (def.log) v = Math.log(Math.max(1e-12, v));
        const k = Math.max(0, Math.min(255, ((v - a) * span) | 0)) * 3;
        px[o] = lut[k]; px[o + 1] = lut[k + 1]; px[o + 2] = lut[k + 2]; px[o + 3] = 255;
      }
    }
    this.ctx.putImageData(this.img, 0, 0);
    return def.log ? [Math.exp(a), Math.exp(b)] : [a, b];
  }

  sample(frame, x, y) {
    const g = this.grid;
    if (!g || !frame) return null;
    const PAD = 2;
    const i = PAD + Math.min(g.nx - 1, Math.max(0, Math.floor(x / g.geom.Lx * g.nx)));
    const j = PAD + Math.min(g.ny - 1, Math.max(0, Math.floor((y - g.geom.y0) / g.geom.Ly * g.ny)));
    const c = j * g.NX + i;
    if (g.cellType[c] !== 0) return { wall: true };
    const f = frame;
    const a = Math.sqrt(G * f.p[c] / f.r[c]);
    return {
      M: Math.hypot(f.u[c], f.v[c]) / a,
      TK: toKelvin(f.p[c] / f.r[c]),
      pKPa: toKPa(f.p[c]),
      speed: toMS(Math.hypot(f.u[c], f.v[c])),
      phiLocal: f.y[c] / g.fuel.Ystoich,
      burning: f.hr[c] > 1e-3,
    };
  }
}
