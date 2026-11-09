// 2D compressible Euler solver with a fuel scalar and simple combustion.
//
// Equations (conservative form):  dU/dt + dF/dx + dG/dy = S
//   U = [rho, rho*u, rho*v, E, rho*Y]
//   Y = unburnt hydrogen mass fraction (passive scalar that can burn)
//
// Method
//   - finite volume on a uniform Cartesian grid
//   - MUSCL reconstruction of primitive variables (van Leer limiter)
//   - HLLC approximate Riemann solver
//   - 2-stage SSP Runge-Kutta (Heun) time stepping, CFL-limited dt
//   - walls via a ghost-cell immersed boundary: each solid cell next to the
//     flow mirrors the flow at its image point across the true wall line,
//     so the ramps are smooth slip walls instead of staircases
//
// Units: everything is non-dimensional. Freestream density = 1 and
// freestream speed of sound = 1, so freestream pressure = 1/gamma and the
// freestream x-velocity equals the flight Mach number.

export const GAMMA = 1.4;
const G = GAMMA;
const GM1 = G - 1;
const PAD = 2;

// Freestream at about 30 km altitude, for converting to physical units.
export const ATMOS = { T: 226.5, p: 1197, a: 301.7, rho: 0.0184, altitude_km: 30 };

// Hydrogen-air chemistry, deliberately simple
export const FUEL = {
  Ystoich: 0.0283,                         // H2 mass fraction at stoichiometric (phi = 1)
  // Effective heat release per kg of H2. The real value is 120 MJ/kg, but this
  // gas has constant specific heat and no dissociation, so the full value would
  // overshoot flame temperatures (~7000 K). 50 MJ/kg gives a realistic ~2900 K.
  Q: 50e6 / (ATMOS.a * ATMOS.a),
  Tign: 850 / (G * ATMOS.T),               // ignition temperature 850 K, non-dim (T = p/rho)
  dTign: 50 / (G * ATMOS.T),               // width of the ignition switch
  Tcap: 3000 / (G * ATMOS.T),              // above ~3000 K dissociation stops net heat release
  rate: 40,                                // burn rate once ignited (per unit time)
  inject: 120,                             // how fast the injector brings Y up to target
};

export const toKelvin = (Tnd) => Tnd * G * ATMOS.T;
export const toKPa = (pnd) => pnd * G * ATMOS.p / 1000;
export const toMS = (vnd) => vnd * ATMOS.a;

function pointInPoly(x, y, pts) {
  let inside = false;
  for (let a = 0, b = pts.length - 1; a < pts.length; b = a++) {
    const [xa, ya] = pts[a], [xb, yb] = pts[b];
    if ((ya > y) !== (yb > y) && x < (xb - xa) * (y - ya) / (yb - ya) + xa) inside = !inside;
  }
  return inside;
}

function nearestOnPoly(x, y, pts) {
  let best = Infinity, bx = 0, by = 0;
  for (let a = 0, b = pts.length - 1; a < pts.length; b = a++) {
    const [xa, ya] = pts[b], [xb, yb] = pts[a];
    const ex = xb - xa, ey = yb - ya;
    const L2 = ex * ex + ey * ey;
    let t = L2 > 0 ? ((x - xa) * ex + (y - ya) * ey) / L2 : 0;
    t = Math.max(0, Math.min(1, t));
    const px = xa + t * ex, py = ya + t * ey;
    const d2 = (x - px) ** 2 + (y - py) ** 2;
    if (d2 < best) { best = d2; bx = px; by = py; }
  }
  return { x: bx, y: by, d2: best };
}

// van Leer limiter slope
function vl(a, b) {
  const ab = a * b;
  return ab > 0 ? (2 * ab) / (a + b) : 0;
}

export class Solver {
  constructor(geom, nx, ny) {
    this.geom = geom;
    this.nx = nx;
    this.ny = ny;
    const NX = nx + 2 * PAD, NY = ny + 2 * PAD;
    this.NX = NX;
    this.NY = NY;
    const N = NX * NY;
    this.N = N;
    this.dx = geom.Lx / nx;
    this.dy = geom.Ly / ny;

    const mk = () => [0, 1, 2, 3, 4].map(() => new Float64Array(N));
    this.U = mk();   // conserved, current
    this.U0 = mk();  // conserved, start of step
    this.W = mk();   // primitives rho,u,v,p,Y
    this.Fx = mk();
    this.Fy = mk();
    this.f = new Float64Array(5);  // scratch flux
    this.hr = new Float32Array(N); // heat release rate (for display)
    this.type = new Uint8Array(N); // 0 fluid, 1 ghost wall, 2 deep wall, 3 domain pad

    this.mach = 6;
    this.phi = 0.4;          // fuel equivalence ratio
    this.cfl = 0.4;
    this.t = 0;
    this.steps = 0;
    this.burnRate = 0;       // smoothed fuel burned per unit time
    this.injRate = 0;        // smoothed fuel injected per unit time

    this.buildMask();
    this.buildGhosts();
    this.buildZones();
    this.reset();
  }

  idx(i, j) { return j * this.NX + i; }
  xc(i) { return (i - PAD + 0.5) * this.dx; }
  yc(j) { return this.geom.y0 + (j - PAD + 0.5) * this.dy; }

  isSolidPoint(x, y) {
    for (const s of this.geom.solids) if (pointInPoly(x, y, s.pts)) return true;
    return false;
  }

  buildMask() {
    const { NX, NY, type } = this;
    for (let j = 0; j < NY; j++) {
      for (let i = 0; i < NX; i++) {
        const c = j * NX + i;
        const pad = i < PAD || i >= NX - PAD || j < PAD || j >= NY - PAD;
        if (pad) { type[c] = 3; continue; }
        type[c] = this.isSolidPoint(this.xc(i), this.yc(j)) ? 2 : 0;
      }
    }
    // solid cells within reach of a fluid stencil become ghost cells
    const mark = new Uint8Array(this.N);
    for (let j = PAD; j < NY - PAD; j++) {
      for (let i = PAD; i < NX - PAD; i++) {
        const c = j * NX + i;
        if (type[c] !== 2) continue;
        let near = false;
        for (let dj = -2; dj <= 2 && !near; dj++)
          for (let di = -2; di <= 2; di++) {
            const ii = i + di, jj = j + dj;
            if (ii < 0 || jj < 0 || ii >= NX || jj >= NY) continue;
            if (type[jj * NX + ii] === 0) { near = true; break; }
          }
        if (near) mark[c] = 1;
      }
    }
    for (let c = 0; c < this.N; c++) if (mark[c]) type[c] = 1;
  }

  buildGhosts() {
    const { NX, NY, type, dx, dy } = this;
    const list = [];
    for (let c = 0; c < this.N; c++) if (type[c] === 1) list.push(c);
    const n = list.length;
    this.gCell = Int32Array.from(list);
    this.gNbr = new Int32Array(n * 4);
    this.gW = new Float64Array(n * 4);
    this.gNx = new Float64Array(n);
    this.gNy = new Float64Array(n);

    for (let g = 0; g < n; g++) {
      const c = list[g];
      const i = c % NX, j = (c / NX) | 0;
      const x = this.xc(i), y = this.yc(j);
      // nearest wall point over all solids
      let best = null;
      for (const s of this.geom.solids) {
        if (!pointInPoly(x, y, s.pts)) continue;
        const q = nearestOnPoly(x, y, s.pts);
        if (!best || q.d2 < best.d2) best = q;
      }
      if (!best) best = { x, y: y + dy, d2: dy * dy };
      let nxv = best.x - x, nyv = best.y - y;
      const len = Math.hypot(nxv, nyv) || 1e-12;
      nxv /= len; nyv /= len;
      this.gNx[g] = nxv;
      this.gNy[g] = nyv;
      // image point: mirror of the ghost center across the wall
      const xi = 2 * best.x - x, yi = 2 * best.y - y;
      const fi = xi / dx - 0.5 + PAD;
      const fj = (yi - this.geom.y0) / dy - 0.5 + PAD;
      const i0 = Math.floor(fi), j0 = Math.floor(fj);
      const tx = fi - i0, ty = fj - j0;
      const cand = [
        [i0, j0, (1 - tx) * (1 - ty)],
        [i0 + 1, j0, tx * (1 - ty)],
        [i0, j0 + 1, (1 - tx) * ty],
        [i0 + 1, j0 + 1, tx * ty],
      ];
      let wsum = 0;
      for (let k = 0; k < 4; k++) {
        let [ii, jj, w] = cand[k];
        ii = Math.max(0, Math.min(NX - 1, ii));
        jj = Math.max(0, Math.min(NY - 1, jj));
        const cc = jj * NX + ii;
        if (type[cc] !== 0) w = 0;
        this.gNbr[g * 4 + k] = cc;
        this.gW[g * 4 + k] = w;
        wsum += w;
      }
      if (wsum < 1e-6) {
        // fall back to the closest fluid cell
        let bd = Infinity, bc = c;
        for (let dj = -4; dj <= 4; dj++)
          for (let di = -4; di <= 4; di++) {
            const ii = i + di, jj = j + dj;
            if (ii < 0 || jj < 0 || ii >= NX || jj >= NY) continue;
            const cc = jj * NX + ii;
            if (type[cc] !== 0) continue;
            const d = di * di + dj * dj;
            if (d < bd) { bd = d; bc = cc; }
          }
        for (let k = 0; k < 4; k++) { this.gNbr[g * 4 + k] = bc; this.gW[g * 4 + k] = k === 0 ? 1 : 0; }
      } else {
        for (let k = 0; k < 4; k++) this.gW[g * 4 + k] /= wsum;
      }
    }
  }

  buildZones() {
    const { NX, NY, type, geom } = this;
    const inj = [];
    for (let j = PAD; j < NY - PAD; j++)
      for (let i = PAD; i < NX - PAD; i++) {
        const c = j * NX + i;
        if (type[c] !== 0) continue;
        const x = this.xc(i), y = this.yc(j);
        if (x >= geom.injector.x0 && x <= geom.injector.x1 && y < geom.yCowl) inj.push(c);
      }
    this.injCells = Int32Array.from(inj);
    this.fluidCells = Int32Array.from([...Array(this.N).keys()].filter((c) => type[c] === 0));
  }

  freestream() {
    const p = 1 / G;
    return { r: 1, u: this.mach, v: 0, p, E: p / GM1 + 0.5 * this.mach * this.mach };
  }

  reset() {
    const fs = this.freestream();
    const [r, mu, mv, E, rY] = this.U;
    r.fill(fs.r); mu.fill(fs.r * fs.u); mv.fill(0); E.fill(fs.E); rY.fill(0);
    const W = this.W;
    W[0].fill(1); W[1].fill(fs.u); W[2].fill(0); W[3].fill(fs.p); W[4].fill(0);
    this.hr.fill(0);
    this.t = 0;
    this.steps = 0;
    this.burnRate = 0;
    this.injRate = 0;
  }

  // conserved -> primitive (fluid cells), then fill walls and domain edges
  prims(U) {
    const r = U[0], mu = U[1], mv = U[2], E = U[3], rY = U[4];
    const W = this.W, Wr = W[0], Wu = W[1], Wv = W[2], Wp = W[3], Wy = W[4];
    const fc = this.fluidCells;
    for (let k = 0; k < fc.length; k++) {
      const c = fc[k];
      const rr = r[c], u = mu[c] / rr, v = mv[c] / rr;
      Wr[c] = rr; Wu[c] = u; Wv[c] = v;
      Wp[c] = GM1 * (E[c] - 0.5 * rr * (u * u + v * v));
      Wy[c] = rY[c] / rr;
    }
    // immersed-boundary ghost cells: mirror normal velocity
    const { gCell, gNbr, gW, gNx, gNy } = this;
    for (let g = 0; g < gCell.length; g++) {
      let ar = 0, au = 0, av = 0, ap = 0, ay = 0;
      for (let k = 0; k < 4; k++) {
        const w = gW[g * 4 + k];
        if (w === 0) continue;
        const cc = gNbr[g * 4 + k];
        ar += w * Wr[cc]; au += w * Wu[cc]; av += w * Wv[cc]; ap += w * Wp[cc]; ay += w * Wy[cc];
      }
      const nx = gNx[g], ny = gNy[g];
      const vn = au * nx + av * ny;
      const c = gCell[g];
      Wr[c] = ar; Wu[c] = au - 2 * vn * nx; Wv[c] = av - 2 * vn * ny; Wp[c] = ap; Wy[c] = ay;
    }
    // domain edges
    const { NX, NY } = this;
    for (let i = PAD; i < NX - PAD; i++) {
      for (let q = 0; q < 5; q++) {
        const A = W[q];
        A[i] = A[NX + i] = A[2 * NX + i];                                   // bottom: zero-gradient
        const top = (NY - PAD - 1) * NX + i;
        A[top + NX] = A[top + 2 * NX] = A[top];                             // top: zero-gradient
      }
    }
    const fs = this.freestream();
    for (let j = 0; j < NY; j++) {
      const row = j * NX;
      // left: supersonic inflow
      for (let i = 0; i < PAD; i++) {
        const c = row + i;
        Wr[c] = 1; Wu[c] = fs.u; Wv[c] = 0; Wp[c] = fs.p; Wy[c] = 0;
      }
      // right: extrapolation outflow
      const last = row + NX - PAD - 1;
      for (let q = 0; q < 5; q++) { const A = W[q]; A[last + 1] = A[last + 2] = A[last]; }
    }
  }

  // HLLC flux for one face. Normal velocity u, tangential v. Writes the
  // flux into the scratch array f = [mass, normal mom, tangential mom, energy, scalar].
  hllc(rL, uL, vL, pL, yL, rR, uR, vR, pR, yR, f) {
    const cL = Math.sqrt(G * pL / rL), cR = Math.sqrt(G * pR / rR);
    const SL = Math.min(uL - cL, uR - cR);
    const SR = Math.max(uL + cL, uR + cR);
    const EL = pL / GM1 + 0.5 * rL * (uL * uL + vL * vL);
    const ER = pR / GM1 + 0.5 * rR * (uR * uR + vR * vR);
    if (SL >= 0) {
      const m = rL * uL;
      f[0] = m; f[1] = m * uL + pL; f[2] = m * vL; f[3] = (EL + pL) * uL; f[4] = m * yL;
      return;
    }
    if (SR <= 0) {
      const m = rR * uR;
      f[0] = m; f[1] = m * uR + pR; f[2] = m * vR; f[3] = (ER + pR) * uR; f[4] = m * yR;
      return;
    }
    const dL = rL * (SL - uL), dR = rR * (SR - uR);
    const SM = (pR - pL + uL * dL - uR * dR) / (dL - dR);
    if (SM >= 0) {
      const k = dL / (SL - SM);
      const m = rL * uL;
      const Es = k * (EL / rL + (SM - uL) * (SM + pL / dL));
      f[0] = m + SL * (k - rL);
      f[1] = m * uL + pL + SL * (k * SM - m);
      f[2] = m * vL + SL * (k - rL) * vL;
      f[3] = (EL + pL) * uL + SL * (Es - EL);
      f[4] = m * yL + SL * (k - rL) * yL;
    } else {
      const k = dR / (SR - SM);
      const m = rR * uR;
      const Es = k * (ER / rR + (SM - uR) * (SM + pR / dR));
      f[0] = m + SR * (k - rR);
      f[1] = m * uR + pR + SR * (k * SM - m);
      f[2] = m * vR + SR * (k - rR) * vR;
      f[3] = (ER + pR) * uR + SR * (Es - ER);
      f[4] = m * yR + SR * (k - rR) * yR;
    }
  }

  // Flux through every face between cells a (left/bottom) and b = a + s,
  // where s = 1 for x-faces and s = NX for y-faces. Wn/Wt are the normal and
  // tangential velocity arrays. Results go to F0..F4 at index a, with the
  // momentum components put back in (x, y) order.
  sweep(s, isX) {
    const { NX, NY, type } = this;
    const Wr = this.W[0], Wp = this.W[3], Wy = this.W[4];
    const Wn = isX ? this.W[1] : this.W[2];
    const Wt = isX ? this.W[2] : this.W[1];
    const F = isX ? this.Fx : this.Fy;
    const F0 = F[0], Fm = isX ? F[1] : F[2], Ft = isX ? F[2] : F[1], F3 = F[3], F4 = F[4];
    const f = this.f;
    const jStart = isX ? PAD : PAD - 1, iStart = isX ? PAD - 1 : PAD;
    for (let j = jStart; j < NY - PAD; j++) {
      const row = j * NX;
      for (let i = iStart; i < NX - PAD; i++) {
        const a = row + i, b = a + s;
        if (type[a] !== 0 && type[b] !== 0) continue;
        const a0 = a - s, b1 = b + s;
        let rL = Wr[a] + 0.5 * vl(Wr[a] - Wr[a0], Wr[b] - Wr[a]);
        let uL = Wn[a] + 0.5 * vl(Wn[a] - Wn[a0], Wn[b] - Wn[a]);
        let vL = Wt[a] + 0.5 * vl(Wt[a] - Wt[a0], Wt[b] - Wt[a]);
        let pL = Wp[a] + 0.5 * vl(Wp[a] - Wp[a0], Wp[b] - Wp[a]);
        let yL = Wy[a] + 0.5 * vl(Wy[a] - Wy[a0], Wy[b] - Wy[a]);
        let rR = Wr[b] - 0.5 * vl(Wr[b] - Wr[a], Wr[b1] - Wr[b]);
        let uR = Wn[b] - 0.5 * vl(Wn[b] - Wn[a], Wn[b1] - Wn[b]);
        let vR = Wt[b] - 0.5 * vl(Wt[b] - Wt[a], Wt[b1] - Wt[b]);
        let pR = Wp[b] - 0.5 * vl(Wp[b] - Wp[a], Wp[b1] - Wp[b]);
        let yR = Wy[b] - 0.5 * vl(Wy[b] - Wy[a], Wy[b1] - Wy[b]);
        if (!(rL > 1e-6 && pL > 1e-7)) { rL = Wr[a]; uL = Wn[a]; vL = Wt[a]; pL = Wp[a]; yL = Wy[a]; }
        if (!(rR > 1e-6 && pR > 1e-7)) { rR = Wr[b]; uR = Wn[b]; vR = Wt[b]; pR = Wp[b]; yR = Wy[b]; }
        this.hllc(rL, uL, vL, pL, yL, rR, uR, vR, pR, yR, f);
        F0[a] = f[0]; Fm[a] = f[1]; Ft[a] = f[2]; F3[a] = f[3]; F4[a] = f[4];
      }
    }
  }

  maxDt() {
    const [Wr, Wu, Wv, Wp] = this.W;
    const fc = this.fluidCells;
    let s = 1e-12;
    const idx = 1 / this.dx, idy = 1 / this.dy;
    for (let k = 0; k < fc.length; k++) {
      const c = fc[k];
      const a = Math.sqrt(G * Wp[c] / Wr[c]);
      const v = (Math.abs(Wu[c]) + a) * idx + (Math.abs(Wv[c]) + a) * idy;
      if (v > s) s = v;
    }
    return this.cfl / s;
  }

  // dU = L(U) * dt added to target: target = base + dt * L(U)  (U already in this.W)
  rhs(dt, base, target, record) {
    const { NX, Fx, Fy, dx, dy } = this;
    this.sweep(1, true);
    this.sweep(NX, false);
    const fc = this.fluidCells;
    const kx = dt / dx, ky = dt / dy;
    for (let q = 0; q < 5; q++) {
      const B = base[q], T = target[q], FX = Fx[q], FY = Fy[q];
      for (let k = 0; k < fc.length; k++) {
        const c = fc[k];
        T[c] = B[c] - kx * (FX[c] - FX[c - 1]) - ky * (FY[c] - FY[c - NX]);
      }
    }
    // sources: fuel injection + combustion
    const [Wr, , , Wp, Wy] = this.W;
    const Tg = target[3], Ty = target[4];
    const Ytarget = this.phi * FUEL.Ystoich;
    let inj = 0, burn = 0;
    if (Ytarget > 0) {
      const ic = this.injCells;
      for (let k = 0; k < ic.length; k++) {
        const c = ic[k];
        const dY = Ytarget - Wy[c];
        if (dY > 0) {
          const s = FUEL.inject * Wr[c] * dY;
          Ty[c] += dt * s;
          inj += s;
        }
      }
    }
    const Tign = FUEL.Tign, dT = FUEL.dTign, A = FUEL.rate, Q = FUEL.Q, Tcap = FUEL.Tcap;
    const hr = this.hr;
    for (let k = 0; k < fc.length; k++) {
      const c = fc[k];
      const Y = Wy[c];
      if (Y < 1e-7) { if (record) hr[c] = 0; continue; }
      const T = Wp[c] / Wr[c];
      const sw = 1 / (1 + Math.exp(-(T - Tign) / dT));
      const cap = Math.min(1, Math.max(0, (Tcap - T) / (0.1 * Tcap)));
      const w = A * Wr[c] * Y * sw * cap;
      Ty[c] -= dt * w;
      Tg[c] += dt * Q * w;
      burn += w;
      if (record) hr[c] = Q * w;
    }
    if (record) {
      const dV = dx * dy;
      const a = Math.min(1, dt / 0.3);
      this.injRate += a * (inj * dV - this.injRate);
      this.burnRate += a * (burn * dV - this.burnRate);
    }
  }

  fixup(U) {
    const [r, mu, mv, E, rY] = U;
    const fc = this.fluidCells;
    for (let k = 0; k < fc.length; k++) {
      const c = fc[k];
      if (!(r[c] > 1e-4)) { r[c] = 1e-4; mu[c] = 0; mv[c] = 0; }
      const ke = 0.5 * (mu[c] * mu[c] + mv[c] * mv[c]) / r[c];
      const pmin = 1e-5;
      if (!(E[c] - ke > pmin / GM1)) E[c] = ke + pmin / GM1;
      if (!(rY[c] >= 0)) rY[c] = 0;
    }
  }

  step() {
    const { U, U0 } = this;
    for (let q = 0; q < 5; q++) U0[q].set(U[q]);
    this.prims(U);
    const dt = this.maxDt();
    // stage 1: U1 = U0 + dt L(U0)
    this.rhs(dt, U0, U, false);
    this.fixup(U);
    // stage 2: U = 0.5 U0 + 0.5 (U1 + dt L(U1))
    this.prims(U);
    this.rhs(dt, U, U, true);
    const fc = this.fluidCells;
    for (let q = 0; q < 5; q++) {
      const A = U[q], B = U0[q];
      for (let k = 0; k < fc.length; k++) { const c = fc[k]; A[c] = 0.5 * (A[c] + B[c]); }
    }
    this.fixup(U);
    this.t += dt;
    this.steps++;
    return dt;
  }

  // --- diagnostics (use the latest primitives) ---
  refreshPrims() { this.prims(this.U); }

  columnAt(x) { return PAD + Math.min(this.nx - 1, Math.max(0, Math.floor(x / this.dx))); }
  rowAt(y) { return PAD + Math.min(this.ny - 1, Math.max(0, Math.floor((y - this.geom.y0) / this.dy))); }

  // mass-flow-weighted Mach number across the duct at station x
  ductMach(x) {
    const i = this.columnAt(x);
    const [Wr, Wu, Wv, Wp] = this.W;
    let num = 0, den = 0, minM = Infinity;
    for (let j = PAD; j < this.NY - PAD; j++) {
      const y = this.yc(j);
      if (y > this.geom.yCowl) break;
      const c = j * this.NX + i;
      if (this.type[c] !== 0) continue;
      const a = Math.sqrt(G * Wp[c] / Wr[c]);
      const M = Math.hypot(Wu[c], Wv[c]) / a;
      const m = Math.abs(Wr[c] * Wu[c]);
      num += m * M; den += m;
      if (M < minM) minM = M;
    }
    return { mean: den > 0 ? num / den : 0, min: minM };
  }

  // maximum static temperature in the duct downstream of the injector
  maxT(x0, x1) {
    const [Wr, , , Wp] = this.W;
    let T = 0;
    const i0 = this.columnAt(x0), i1 = this.columnAt(x1);
    for (let j = PAD; j < this.NY - PAD; j++) {
      if (this.yc(j) > this.geom.yCowl) break;
      for (let i = i0; i <= i1; i++) {
        const c = j * this.NX + i;
        if (this.type[c] === 0) T = Math.max(T, Wp[c] / Wr[c]);
      }
    }
    return T;
  }

  // Net axial force on the walls from pressure (positive = thrust, i.e. pushing
  // the vehicle upstream, toward -x). Uses the wall-adjacent fluid pressure on
  // every face between a fluid cell and a wall cell. Inviscid, so no skin friction.
  axialForce() {
    const { NX, NY, type, dy } = this;
    const Wp = this.W[3];
    const pinf = 1 / G;
    let F = 0;
    for (let j = PAD; j < NY - PAD; j++) {
      const row = j * NX;
      for (let i = PAD; i < NX - PAD - 1; i++) {
        const a = row + i, b = a + 1;
        const fa = type[a] === 0, fb = type[b] === 0;
        if (fa && !fb && type[b] !== 3) F -= (Wp[a] - pinf) * dy; // wall on the right: pushed +x (drag)
        else if (!fa && fb && type[a] !== 3) F += (Wp[b] - pinf) * dy; // wall on the left: pushed -x (thrust)
      }
    }
    return F; // non-dim force per unit span; physical N/m = F * G * ATMOS.p * 1 m
  }

  sample(x, y) {
    const i = this.columnAt(x), j = this.rowAt(y);
    const c = j * this.NX + i;
    if (this.type[c] !== 0) return null;
    const [Wr, Wu, Wv, Wp, Wy] = this.W;
    const a = Math.sqrt(G * Wp[c] / Wr[c]);
    return {
      M: Math.hypot(Wu[c], Wv[c]) / a,
      T: Wp[c] / Wr[c],
      p: Wp[c],
      rho: Wr[c],
      u: Wu[c],
      v: Wv[c],
      Y: Wy[c],
    };
  }
}
