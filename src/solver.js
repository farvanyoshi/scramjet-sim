// 2D compressible Reynolds-averaged Navier-Stokes solver for the scramjet.
//
// Conserved variables (per unit volume):
//   U = [rho, rho*u, rho*v, E, rho*nu~, rho*Z, rho*Y]
//     nu~ : Spalart-Allmaras turbulence variable (eddy viscosity ~ rho*nu~*fv1)
//     Z   : mixture fraction (all hydrogen injected, burnt or not)
//     Y   : unburnt hydrogen mass fraction
//
// Physics
//   - inviscid fluxes: MUSCL (van Leer limiter) + HLLC Riemann solver
//   - viscous stresses and heat conduction (Sutherland viscosity, Pr = 0.72)
//   - Spalart-Allmaras one-equation turbulence model
//   - wall functions: Spalding's law of the wall for skin friction, Reynolds
//     analogy for heat transfer, walls held at 1000 K (fuel-cooled)
//   - thermally perfect air (vibrational energy of N2 and O2)
//   - hydrogen injected from wall slots on the combustor floor and cowl,
//     mixing by turbulent diffusion, one-step H2 + O2 reaction limited by
//     available oxygen and by dissociation above ~2700 K
//   - 2-stage SSP Runge-Kutta in time
//   - walls: ghost-cell immersed boundary (mirror state at the image point)
//
// Units: see gas.js. Freestream rho = 1, a = 1, so u_inf = Mach number.

import {
  ATMOS, TS, P_SCALE, Q_SCALE, toKelvin, fromKelvin,
  thetaFromE, eOf, gammaOf, muOf, E_MIN, TH_MIN,
  PR, PR_T, SC, SC_T, SA, CHEM, WALL,
} from './gas.js';

export { toKelvin, toKPa, toMS, ATMOS, CHEM } from './gas.js';

const PAD = 2;
const NV = 7;
const VISC_ON = 0.5;   // evaluate viscous fluxes where mu_t > 2 mu
const G0 = 1.4;
export const MODELS = ['turbulent', 'inviscid'];

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

export class Solver {
  constructor(geom, nx, ny, model = 'turbulent') {
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
    this.model = model;

    const mk = (n) => Array.from({ length: n }, () => new Float64Array(N));
    this.U = mk(NV);
    this.U0 = mk(NV);
    this.W = mk(NV);            // rho, u, v, p, nu~, Z, Y
    this.Th = new Float64Array(N);  // theta = p/rho
    this.Ge = new Float64Array(N);  // effective gamma 1 + p/(rho e)
    this.Ga = new Float64Array(N);  // true gamma (sound speed)
    this.Mu = new Float64Array(N);  // laminar viscosity
    this.Mt = new Float64Array(N);  // eddy viscosity
    const mk32 = (n) => Array.from({ length: n }, () => new Float32Array(N));
    this.S = mk32(NV);          // limited slopes (reused per direction)
    this.Fx = mk32(NV);
    this.Fy = mk32(NV);
    this.f = new Float64Array(5);
    this.hr = new Float32Array(N);  // heat release rate (display)
    this.type = new Uint8Array(N);  // 0 fluid, 1 ghost wall, 2 deep wall, 3 domain pad
    this.dist = new Float64Array(N); // distance to nearest wall

    this.mach = 6;
    this.phi = 0.6;
    this.cfl = 0.5;
    this.t = 0;
    this.steps = 0;
    this.mdotAir = 0;
    this.burnRate = 0;
    this.injRate = 0;

    this.buildMask();
    this.buildGhosts();
    this.buildWallCells();
    this.buildInjectors();
    this.fluidCells = Int32Array.from([...Array(N).keys()].filter((c) => this.type[c] === 0));
    this.reset();
  }

  xc(i) { return (i - PAD + 0.5) * this.dx; }
  yc(j) { return this.geom.y0 + (j - PAD + 0.5) * this.dy; }

  // ---------------------------------------------------------------- setup
  buildMask() {
    const { NX, NY, type } = this;
    for (let j = 0; j < NY; j++) {
      for (let i = 0; i < NX; i++) {
        const c = j * NX + i;
        if (i < PAD || i >= NX - PAD || j < PAD || j >= NY - PAD) { type[c] = 3; continue; }
        const x = this.xc(i), y = this.yc(j);
        let solid = false;
        for (const s of this.geom.solids) if (pointInPoly(x, y, s.pts)) solid = true;
        type[c] = solid ? 2 : 0;
      }
    }
    // solid cells within two cells of the flow become ghost cells
    const mark = new Uint8Array(this.N);
    for (let j = PAD; j < NY - PAD; j++)
      for (let i = PAD; i < NX - PAD; i++) {
        const c = j * NX + i;
        if (type[c] !== 2) continue;
        outer: for (let dj = -2; dj <= 2; dj++)
          for (let di = -2; di <= 2; di++) {
            const ii = i + di, jj = j + dj;
            if (ii >= 0 && jj >= 0 && ii < NX && jj < NY && type[jj * NX + ii] === 0) { mark[c] = 1; break outer; }
          }
      }
    for (let c = 0; c < this.N; c++) if (mark[c]) type[c] = 1;

    // wall distance for every non-wall cell
    for (let j = 0; j < NY; j++)
      for (let i = 0; i < NX; i++) {
        const c = j * NX + i;
        if (type[c] === 1 || type[c] === 2) { this.dist[c] = 0; continue; }
        const x = this.xc(i), y = this.yc(j);
        let d2 = Infinity;
        for (const s of this.geom.solids) d2 = Math.min(d2, nearestOnPoly(x, y, s.pts).d2);
        this.dist[c] = Math.min(1, Math.sqrt(d2));
      }
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
      const xi = 2 * best.x - x, yi = 2 * best.y - y;
      const fi = xi / dx - 0.5 + PAD;
      const fj = (yi - this.geom.y0) / dy - 0.5 + PAD;
      const i0 = Math.floor(fi), j0 = Math.floor(fj);
      const tx = fi - i0, ty = fj - j0;
      const cand = [[i0, j0, (1 - tx) * (1 - ty)], [i0 + 1, j0, tx * (1 - ty)], [i0, j0 + 1, (1 - tx) * ty], [i0 + 1, j0 + 1, tx * ty]];
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

  // Fluid cells that touch a wall get the wall-function treatment.
  buildWallCells() {
    const { NX, type, dx, dy } = this;
    const cells = [], nxs = [], nys = [], ds = [], aV = [], px = [], py = [];
    for (let c = 0; c < this.N; c++) {
      if (type[c] !== 0) continue;
      const gl = type[c - 1] === 1, gr = type[c + 1] === 1, gb = type[c - NX] === 1, gt = type[c + NX] === 1;
      if (!(gl || gr || gb || gt)) continue;
      const i = c % NX, j = (c / NX) | 0;
      const x = this.xc(i), y = this.yc(j);
      let best = null;
      for (const s of this.geom.solids) {
        const q = nearestOnPoly(x, y, s.pts);
        if (!best || q.d2 < best.d2) best = q;
      }
      const d = Math.max(0.15 * Math.min(dx, dy), Math.sqrt(best.d2));
      let nx = (x - best.x) / d, ny = (y - best.y) / d;
      const len = Math.hypot(nx, ny) || 1;
      nx /= len; ny /= len;
      // staircase faces weighted by |n| add up to the true wall length
      const A = ((gl ? 1 : 0) + (gr ? 1 : 0)) * dy * Math.abs(nx) + ((gb ? 1 : 0) + (gt ? 1 : 0)) * dx * Math.abs(ny);
      if (A <= 0) continue;
      cells.push(c); nxs.push(nx); nys.push(ny); ds.push(d); aV.push(A / (dx * dy)); px.push(best.x); py.push(best.y);
    }
    const n = cells.length;
    this.wCell = Int32Array.from(cells);
    this.wNx = Float64Array.from(nxs);
    this.wNy = Float64Array.from(nys);
    this.wD = Float64Array.from(ds);
    this.wAV = Float64Array.from(aV);      // wall area / cell volume
    this.wPx = Float64Array.from(px);      // wall point (for drawing)
    this.wPy = Float64Array.from(py);
    this.wUp = new Float64Array(n).fill(20); // u+ from last solve (warm start)
    this.wTau = new Float64Array(n);       // wall shear stress
    this.wQ = new Float32Array(n);         // wall heat flux (into wall)
    this.wTx = new Float64Array(n);        // x-component of flow direction at the wall
  }

  // Wall-slot injectors: cells within the jet penetration height of the
  // combustor floor and of the cowl underside, over the injector x-range.
  buildInjectors() {
    const { NX, NY, type, geom } = this;
    const inj = geom.injector;
    const body = geom.solids[0].pts, cowl = geom.solids[1].pts;
    const lists = [[], []], normals = [[0, 0], [0, 0]];
    for (let j = PAD; j < NY - PAD; j++)
      for (let i = PAD; i < NX - PAD; i++) {
        const c = j * NX + i;
        if (type[c] !== 0) continue;
        const x = this.xc(i), y = this.yc(j);
        if (x < inj.x0 || x > inj.x1 || y > geom.yCowl) continue;
        const qb = nearestOnPoly(x, y, body), qc = nearestOnPoly(x, y, cowl);
        const db = Math.sqrt(qb.d2), dc = Math.sqrt(qc.d2);
        const side = db <= dc ? 0 : 1;
        if (Math.min(db, dc) > inj.penetration) continue;
        const q = side === 0 ? qb : qc, dd = side === 0 ? db : dc;
        lists[side].push(c);
        normals[side][0] += (x - q.x) / (dd || 1);
        normals[side][1] += (y - q.y) / (dd || 1);
      }
    this.injSides = lists.map((l, k) => {
      const [a, b] = normals[k];
      const L = Math.hypot(a, b) || 1;
      return { cells: Int32Array.from(l), nx: a / L, ny: b / L };
    });
    this.jet = {
      V: 1200 / ATMOS.a,                 // sonic H2 jet, ~1200 m/s
      th: fromKelvin(250),               // static temperature of the jet
    };
  }

  freestream() {
    const th = 1 / G0;                   // 226.5 K
    const e = eOf(th);
    const nu = muOf(th);                 // rho = 1, so nu = mu
    return { r: 1, u: this.mach, v: 0, p: th, th, e, E: e + 0.5 * this.mach * this.mach, nut: 3 * nu, mu: nu };
  }

  setModel(model) {
    this.model = model;
  }

  reset() {
    const fs = this.freestream();
    const U = this.U, W = this.W;
    U[0].fill(1); U[1].fill(fs.u); U[2].fill(0); U[3].fill(fs.E); U[4].fill(fs.nut); U[5].fill(0); U[6].fill(0);
    W[0].fill(1); W[1].fill(fs.u); W[2].fill(0); W[3].fill(fs.p); W[4].fill(fs.nut); W[5].fill(0); W[6].fill(0);
    this.Th.fill(fs.th); this.Ge.fill(1.4); this.Ga.fill(1.4); this.Mu.fill(fs.mu); this.Mt.fill(0);
    this.hr.fill(0);
    this.wUp.fill(20); this.wTau.fill(0); this.wQ.fill(0);
    this.t = 0;
    this.steps = 0;
    this.mdotAir = 0;
    this.burnRate = 0;
    this.injRate = 0;
  }

  // ---------------------------------------------------------------- primitives
  prims(U) {
    const r = U[0], mu = U[1], mv = U[2], E = U[3], rn = U[4], rz = U[5], ry = U[6];
    const W = this.W, Wr = W[0], Wu = W[1], Wv = W[2], Wp = W[3], Wn = W[4], Wz = W[5], Wy = W[6];
    const { Th, Ge, Ga, Mu, Mt } = this;
    const turb = this.model === 'turbulent';
    const cv13 = SA.cv1 ** 3;
    const fc = this.fluidCells;
    for (let k = 0; k < fc.length; k++) {
      const c = fc[k];
      const rr = r[c], u = mu[c] / rr, v = mv[c] / rr;
      let e = E[c] / rr - 0.5 * (u * u + v * v);
      if (e < E_MIN) e = E_MIN;
      const th = thetaFromE(e);
      Wr[c] = rr; Wu[c] = u; Wv[c] = v; Wp[c] = rr * th;
      Th[c] = th; Ge[c] = 1 + th / e; Ga[c] = gammaOf(th);
      const m = muOf(th);
      Mu[c] = m;
      const nt = rn[c] / rr;
      Wn[c] = nt; Wz[c] = rz[c] / rr; Wy[c] = ry[c] / rr;
      if (turb) {
        const chi = (rr * nt) / m, c3 = chi * chi * chi;
        Mt[c] = rr * nt * c3 / (c3 + cv13);
      } else Mt[c] = 0;
    }
    // immersed-boundary ghost cells: mirror normal velocity at the image point
    const { gCell, gNbr, gW, gNx, gNy } = this;
    const ARR = [Wr, Wp, Wn, Wz, Wy, Ge, Ga, Mu, Mt];
    for (let g = 0; g < gCell.length; g++) {
      const c = gCell[g];
      let au = 0, av = 0;
      for (let q = 0; q < ARR.length; q++) ARR[q][c] = 0;
      for (let k = 0; k < 4; k++) {
        const w = gW[g * 4 + k];
        if (w === 0) continue;
        const cc = gNbr[g * 4 + k];
        au += w * Wu[cc]; av += w * Wv[cc];
        for (let q = 0; q < ARR.length; q++) ARR[q][c] += w * ARR[q][cc];
      }
      const vn = au * gNx[g] + av * gNy[g];
      Wu[c] = au - 2 * vn * gNx[g];
      Wv[c] = av - 2 * vn * gNy[g];
      Th[c] = Wp[c] / Wr[c];
    }
    // domain edges
    const { NX, NY } = this;
    const ALL = [Wr, Wu, Wv, Wp, Wn, Wz, Wy, Th, Ge, Ga, Mu, Mt];
    for (let q = 0; q < ALL.length; q++) {
      const A = ALL[q];
      for (let i = PAD; i < NX - PAD; i++) {
        A[i] = A[NX + i] = A[2 * NX + i];
        const top = (NY - PAD - 1) * NX + i;
        A[top + NX] = A[top + 2 * NX] = A[top];
      }
      for (let j = 0; j < NY; j++) {
        const last = j * NX + NX - PAD - 1;
        A[last + 1] = A[last + 2] = A[last];
      }
    }
    const fs = this.freestream();
    const fsMt = turb ? (() => { const chi = fs.nut / fs.mu, c3 = chi ** 3; return fs.nut * c3 / (c3 + cv13); })() : 0;
    for (let j = 0; j < NY; j++)
      for (let i = 0; i < PAD; i++) {
        const c = j * NX + i;
        Wr[c] = 1; Wu[c] = fs.u; Wv[c] = 0; Wp[c] = fs.p; Wn[c] = fs.nut; Wz[c] = 0; Wy[c] = 0;
        Th[c] = fs.th; Ge[c] = 1 + fs.th / fs.e; Ga[c] = 1.4; Mu[c] = fs.mu; Mt[c] = fsMt;
      }
  }

  // ---------------------------------------------------------------- fluxes
  hllc(rL, uL, vL, pL, gL, aL, rR, uR, vR, pR, gR, aR, f) {
    const cL = Math.sqrt(aL * pL / rL), cR = Math.sqrt(aR * pR / rR);
    const SL = Math.min(uL - cL, uR - cR);
    const SR = Math.max(uL + cL, uR + cR);
    const EL = pL / (gL - 1) + 0.5 * rL * (uL * uL + vL * vL);
    const ER = pR / (gR - 1) + 0.5 * rR * (uR * uR + vR * vR);
    if (SL >= 0) {
      const m = rL * uL;
      f[0] = m; f[1] = m * uL + pL; f[2] = m * vL; f[3] = (EL + pL) * uL;
      return true;
    }
    if (SR <= 0) {
      const m = rR * uR;
      f[0] = m; f[1] = m * uR + pR; f[2] = m * vR; f[3] = (ER + pR) * uR;
      return false;
    }
    const dL = rL * (SL - uL), dR = rR * (SR - uR);
    const SM = (pR - pL + uL * dL - uR * dR) / (dL - dR);
    if (SM >= 0) {
      const k = dL / (SL - SM), m = rL * uL;
      const Es = k * (EL / rL + (SM - uL) * (SM + pL / dL));
      f[0] = m + SL * (k - rL);
      f[1] = m * uL + pL + SL * (k * SM - m);
      f[2] = m * vL + SL * (k - rL) * vL;
      f[3] = (EL + pL) * uL + SL * (Es - EL);
      return true;
    }
    const k = dR / (SR - SM), m = rR * uR;
    const Es = k * (ER / rR + (SM - uR) * (SM + pR / dR));
    f[0] = m + SR * (k - rR);
    f[1] = m * uR + pR + SR * (k * SM - m);
    f[2] = m * vR + SR * (k - rR) * vR;
    f[3] = (ER + pR) * uR + SR * (Es - ER);
    return false;
  }

  // Limited slopes of all primitives along stride s (1 = x, NX = y)
  slopes(s) {
    const N = this.N;
    for (let q = 0; q < NV; q++) {
      const A = this.W[q], Sq = this.S[q];
      for (let c = s; c < N - s; c++) {
        const a = A[c] - A[c - s], b = A[c + s] - A[c];
        const ab = a * b;
        Sq[c] = ab > 0 ? (2 * ab) / (a + b) : 0;
      }
    }
  }

  // Inviscid + viscous flux through every face between a and b = a + s.
  sweep(isX) {
    const { NX, NY, type, Th, Ge, Ga, Mu, Mt } = this;
    const s = isX ? 1 : NX, t = isX ? NX : 1;
    const hn = isX ? this.dx : this.dy, ht = isX ? this.dy : this.dx;
    this.slopes(s);
    const W = this.W, Sl = this.S;
    const Wr = W[0], Wp = W[3], Wnu = W[4], Wz = W[5], Wy = W[6];
    const Wn = isX ? W[1] : W[2], Wt = isX ? W[2] : W[1];
    const Sr = Sl[0], Sp = Sl[3], Snu = Sl[4], Sz = Sl[5], Sy = Sl[6];
    const Sn = isX ? Sl[1] : Sl[2], St = isX ? Sl[2] : Sl[1];
    const F = isX ? this.Fx : this.Fy;
    const F0 = F[0], Fm = isX ? F[1] : F[2], Ft = isX ? F[2] : F[1], F3 = F[3], F4 = F[4], F5 = F[5], F6 = F[6];
    const f = this.f;
    const visc = this.model !== 'inviscid';
    const turb = this.model === 'turbulent';
    const ihn = 1 / hn, iht4 = 1 / (4 * ht);
    const iPr = 1 / PR, iPrt = 1 / PR_T, iSc = 1 / SC, iSct = 1 / SC_T, iSig = 1 / SA.sigma;
    const jStart = isX ? PAD : PAD - 1, iStart = isX ? PAD - 1 : PAD;
    for (let j = jStart; j < NY - PAD; j++) {
      const row = j * NX;
      for (let i = iStart; i < NX - PAD; i++) {
        const a = row + i, b = a + s;
        if (type[a] !== 0 && type[b] !== 0) continue;
        let rL = Wr[a] + 0.5 * Sr[a], uL = Wn[a] + 0.5 * Sn[a], vL = Wt[a] + 0.5 * St[a], pL = Wp[a] + 0.5 * Sp[a];
        let rR = Wr[b] - 0.5 * Sr[b], uR = Wn[b] - 0.5 * Sn[b], vR = Wt[b] - 0.5 * St[b], pR = Wp[b] - 0.5 * Sp[b];
        let hoL = 0.5, hoR = 0.5;
        if (!(rL > 1e-6 && pL > 1e-7)) { rL = Wr[a]; uL = Wn[a]; vL = Wt[a]; pL = Wp[a]; hoL = 0; }
        if (!(rR > 1e-6 && pR > 1e-7)) { rR = Wr[b]; uR = Wn[b]; vR = Wt[b]; pR = Wp[b]; hoR = 0; }
        const left = this.hllc(rL, uL, vL, pL, Ge[a], Ga[a], rR, uR, vR, pR, Ge[b], Ga[b], f);
        const m = f[0];
        let fm = f[1], ft = f[2], fe = f[3], fn, fz, fy;
        // passive scalars ride on the mass flux, taken from the upwind side
        if (left) { fn = m * (Wnu[a] + hoL * Snu[a]); fz = m * (Wz[a] + hoL * Sz[a]); fy = m * (Wy[a] + hoL * Sy[a]); }
        else { fn = m * (Wnu[b] - hoR * Snu[b]); fz = m * (Wz[b] - hoR * Sz[b]); fy = m * (Wy[b] - hoR * Sy[b]); }

        // Viscous terms. Outside turbulent regions the molecular viscosity is
        // far below what the grid resolves (cell Reynolds number ~30,000), so
        // they are only evaluated where the eddy viscosity is active.
        if (visc && Mt[a] + Mt[b] > 4 * (Mu[a] + Mu[b]) * VISC_ON) {
          const mu = 0.5 * (Mu[a] + Mu[b]), mut = 0.5 * (Mt[a] + Mt[b]);
          const me = mu + mut;
          const dun = (Wn[b] - Wn[a]) * ihn, dut = (Wt[b] - Wt[a]) * ihn;
          const dunT = (Wn[a + t] - Wn[a - t] + Wn[b + t] - Wn[b - t]) * iht4;
          const dutT = (Wt[a + t] - Wt[a - t] + Wt[b + t] - Wt[b - t]) * iht4;
          const tnn = me * ((4 / 3) * dun - (2 / 3) * dutT);
          const tnt = me * (dut + dunT);
          const gam = 0.5 * (Ga[a] + Ga[b]);
          const cp = gam / (gam - 1);
          const k = cp * (mu * iPr + mut * iPrt);
          const un = 0.5 * (Wn[a] + Wn[b]), ut = 0.5 * (Wt[a] + Wt[b]);
          fm -= tnn;
          ft -= tnt;
          fe -= un * tnn + ut * tnt + k * (Th[b] - Th[a]) * ihn;
          const D = mu * iSc + mut * iSct;
          fz -= D * (Wz[b] - Wz[a]) * ihn;
          fy -= D * (Wy[b] - Wy[a]) * ihn;
          if (turb) {
            const rf = 0.5 * (Wr[a] + Wr[b]), nf = 0.5 * (Wnu[a] + Wnu[b]);
            fn -= iSig * (mu + rf * nf) * (Wnu[b] - Wnu[a]) * ihn;
          }
        }
        F0[a] = m; Fm[a] = fm; Ft[a] = ft; F3[a] = fe; F4[a] = fn; F5[a] = fz; F6[a] = fy;
      }
    }
  }

  maxDt() {
    const W = this.W, Wr = W[0], Wu = W[1], Wv = W[2], Wp = W[3];
    const { Ga, Mu, Mt } = this;
    const fc = this.fluidCells;
    const idx = 1 / this.dx, idy = 1 / this.dy;
    const vfac = 2 * (idx * idx + idy * idy) * (1.4 / PR);
    const visc = this.model !== 'inviscid';
    let s = 1e-12;
    for (let k = 0; k < fc.length; k++) {
      const c = fc[k];
      const a = Math.sqrt(Ga[c] * Wp[c] / Wr[c]);
      let v = (Math.abs(Wu[c]) + a) * idx + (Math.abs(Wv[c]) + a) * idy;
      if (visc) v += vfac * (Mu[c] + Mt[c]) / Wr[c];
      if (v > s) s = v;
    }
    return this.cfl / s;
  }

  // ---------------------------------------------------------------- RHS
  // target = base + dt * L(U), where the primitives of U are already in this.W
  rhs(dt, base, target, record) {
    const { NX, Fx, Fy, dx, dy } = this;
    this.sweep(true);
    this.sweep(false);
    const fc = this.fluidCells;
    const kx = dt / dx, ky = dt / dy;
    for (let q = 0; q < NV; q++) {
      const B = base[q], T = target[q], FX = Fx[q], FY = Fy[q];
      for (let k = 0; k < fc.length; k++) {
        const c = fc[k];
        T[c] = B[c] - kx * (FX[c] - FX[c - 1]) - ky * (FY[c] - FY[c - NX]);
      }
    }
    if (this.model === 'turbulent') this.saSource(dt, target);
    if (this.model !== 'inviscid') this.wallFunctions(dt, target, record);
    this.inject(dt, target, record);
    this.chemistry(dt, target, record);
  }

  // Spalart-Allmaras production / destruction / cross-diffusion
  saSource(dt, T) {
    const { NX, dist, Mu } = this;
    const W = this.W, Wr = W[0], Wu = W[1], Wv = W[2], Wn = W[4];
    const i2x = 1 / (2 * this.dx), i2y = 1 / (2 * this.dy);
    const { cb1, cb2, sigma, kappa, cw1, cw2, cw3, cv1 } = SA;
    const k2 = kappa * kappa, cv13 = cv1 ** 3, cw36 = cw3 ** 6;
    const fwSmall = Math.cbrt(Math.sqrt((1 + cw36) / cw36));
    const Tn = T[4];
    const fc = this.fluidCells;
    for (let k = 0; k < fc.length; k++) {
      const c = fc[k];
      const r = Wr[c], nt = Wn[c] > 0 ? Wn[c] : 0;
      const d = dist[c] > 1e-4 ? dist[c] : 1e-4;
      const om = Math.abs((Wv[c + 1] - Wv[c - 1]) * i2x - (Wu[c + NX] - Wu[c - NX]) * i2y);
      const nu = Mu[c] / r;
      const chi = nt / nu, c3 = chi * chi * chi;
      const fv1 = c3 / (c3 + cv13);
      const fv2 = 1 - chi / (1 + chi * fv1);
      const kd2 = k2 * d * d;
      let St = om + (nt / kd2) * fv2;
      if (St < 0.3 * om) St = 0.3 * om;
      if (St < 1e-10) St = 1e-10;
      let rr = nt / (St * kd2);
      if (rr > 10) rr = 10;
      let fw;
      if (rr < 1e-3) {
        fw = (1 - cw2) * rr * fwSmall;            // limit for small r (far from walls)
      } else {
        const r2 = rr * rr, r6 = r2 * r2 * r2;
        const g = rr + cw2 * (r6 - rr);
        const g2 = g * g, g6 = g2 * g2 * g2;
        fw = g * Math.cbrt(Math.sqrt((1 + cw36) / (g6 + cw36)));
      }
      const gx = (Wn[c + 1] - Wn[c - 1]) * i2x, gy = (Wn[c + NX] - Wn[c - NX]) * i2y;
      const prod = r * cb1 * St * nt + r * (cb2 / sigma) * (gx * gx + gy * gy);
      const Drate = cw1 * fw * nt / (d * d);       // destruction, treated implicitly
      Tn[c] = (Tn[c] + dt * prod) / (1 + dt * Drate);
    }
  }

  // Wall functions: friction and heat flux at cells next to a wall
  wallFunctions(dt, T, record) {
    const { wCell, wNx, wNy, wD, wAV, wUp, wTau, wQ, wTx, Th, Ga } = this;
    const W = this.W, Wr = W[0], Wu = W[1], Wv = W[2], Wp = W[3];
    const thwMax = fromKelvin(WALL.T);
    const kap = SA.kappa, eKB = Math.exp(-kap * WALL.B);
    const rec = Math.cbrt(PR), pr23 = Math.pow(PR, 2 / 3);
    const turb = this.model === 'turbulent';
    const Tmx = T[1], Tmy = T[2], TE = T[3], Tn = T[4];
    for (let w = 0; w < wCell.length; w++) {
      const c = wCell[w];
      const nx = wNx[w], ny = wNy[w], d = wD[w];
      const u = Wu[c], v = Wv[c];
      const vn = u * nx + v * ny;
      const tx = u - vn * nx, ty = v - vn * ny;
      const U = Math.hypot(tx, ty);
      const th = Th[c];
      const gam = Ga[c], cp = gam / (gam - 1);
      const thaw = th + rec * U * U / (2 * cp);
      // a fuel-cooled wall never runs hotter than the gas recovery temperature
      const thw = thaw < thwMax ? thaw : thwMax;
      const ths = Math.max(TH_MIN, 0.5 * (thw + th) + 0.22 * (thaw - th));   // Eckert reference temperature
      const rs = Wp[c] / ths;
      const nus = muOf(ths) / rs;
      let tau = 0, q = 0;
      if (U > 1e-6) {
        let up;
        if (turb) {
          // Spalding's law of the wall, solved for u+ by Newton iteration
          const Red = U * d / nus;
          up = wUp[w];
          for (let it = 0; it < 6; it++) {
            const ku = kap * up, ek = Math.exp(ku);
            const sp = up + eKB * (ek - 1 - ku - 0.5 * ku * ku - ku * ku * ku / 6);
            const h = Red / up - sp;
            const dh = -Red / (up * up) - 1 - eKB * kap * (ek - 1 - ku - 0.5 * ku * ku);
            let nu2 = up - h / dh;
            if (!(nu2 > 0.5)) nu2 = 0.5;
            if (nu2 > 80) nu2 = 80;
            if (Math.abs(nu2 - up) < 1e-6 * up) { up = nu2; break; }
            up = nu2;
          }
          wUp[w] = up;
          const ut = U / up;
          tau = rs * ut * ut;
        } else {
          tau = rs * nus * U / d;   // linear (laminar) profile
        }
        q = (tau / U) * cp * (thaw - thw) / pr23;
        const f = dt * tau * wAV[w] / U;
        Tmx[c] -= f * tx;
        Tmy[c] -= f * ty;
        TE[c] -= dt * q * wAV[w];
        if (turb) {
          // log-layer eddy viscosity in the wall cell: mu_t = kappa d sqrt(rho tau_w)
          const nt = kap * d * Math.sqrt(tau / Wr[c]);
          Tn[c] = T[0][c] * nt;
        }
      }
      if (record) { wTau[w] = tau; wQ[w] = q; wTx[w] = U > 1e-6 ? tx / U : 0; }
    }
  }

  // Hydrogen jets from the wall slots
  inject(dt, T, record) {
    const mdot = this.phi * CHEM.fst * this.mdotAir;   // kg/s per m span, non-dim
    if (record) this.injRateNow = mdot;
    if (!(mdot > 0)) return;
    const V = this.jet.V, th = this.jet.th;
    const h0 = eOf(th) + th + 0.5 * V * V;             // jet stagnation enthalpy (air-like)
    const vol = this.dx * this.dy;
    for (const side of this.injSides) {
      const n = side.cells.length;
      if (!n) continue;
      const sm = (0.5 * mdot) / (n * vol) * dt;
      for (let k = 0; k < n; k++) {
        const c = side.cells[k];
        T[0][c] += sm;
        T[1][c] += sm * V * side.nx;
        T[2][c] += sm * V * side.ny;
        T[3][c] += sm * h0;
        T[5][c] += sm;
        T[6][c] += sm;
      }
    }
  }

  // One-step H2 + O2 -> H2O with an Arrhenius rate, limited by available
  // oxygen and by dissociation at very high temperature. Integrated exactly
  // over the step for a frozen rate constant, so it is stable when fast.
  chemistry(dt, T, record) {
    const W = this.W, Wr = W[0], Wz = W[5], Wy = W[6];
    const { Th, hr } = this;
    const { A, Ta, Q, YO2air, sO2, Tdiss, dTdiss } = CHEM;
    const TaS = Ta / TS, thD = fromKelvin(Tdiss), dthD = fromKelvin(dTdiss);
    const eTop = eOf(thD + dthD);
    const TE = T[3], TY = T[6];
    const fc = this.fluidCells;
    let burnSum = 0;
    for (let k = 0; k < fc.length; k++) {
      const c = fc[k];
      const Y = Wy[c];
      if (Y < 1e-7) { if (record) hr[c] = 0; continue; }
      const Z = Wz[c], th = Th[c], r = Wr[c];
      const YO = YO2air * (1 - Z) - sO2 * (Z - Y);
      if (YO <= 1e-7) { if (record) hr[c] = 0; continue; }
      let taper = (thD + dthD - th) / (2 * dthD);
      if (taper > 1) taper = 1;
      if (taper <= 0) { if (record) hr[c] = 0; continue; }
      const kr = A * r * YO * Math.exp(-TaS / th) * taper;
      let burned = r * Y * (1 - Math.exp(-kr * dt));
      const maxO = (r * YO) / sO2;
      if (burned > maxO) burned = maxO;
      // don't release more heat in one step than it takes to reach the
      // dissociation limit (keeps fast burning from overshooting)
      const maxQ = r * (eTop - eOf(th)) / Q;
      if (burned > maxQ) burned = maxQ > 0 ? maxQ : 0;
      TY[c] -= burned;
      TE[c] += Q * burned;
      burnSum += burned;
      if (record) hr[c] = (Q * burned) / dt;
    }
    if (record) {
      const vol = this.dx * this.dy;
      const a = Math.min(1, dt / 0.3);
      this.burnRate += a * ((burnSum * vol) / dt - this.burnRate);
      this.injRate += a * ((this.injRateNow || 0) - this.injRate);
    }
  }

  fixup(U) {
    const r = U[0], mu = U[1], mv = U[2], E = U[3], rn = U[4], rz = U[5], ry = U[6];
    const fc = this.fluidCells;
    for (let k = 0; k < fc.length; k++) {
      const c = fc[k];
      if (!(r[c] > 1e-4)) { r[c] = 1e-4; mu[c] = 0; mv[c] = 0; }
      const ke = 0.5 * (mu[c] * mu[c] + mv[c] * mv[c]) / r[c];
      if (!(E[c] - ke > r[c] * E_MIN)) E[c] = ke + r[c] * E_MIN;
      if (!(rn[c] >= 0)) rn[c] = 0;
      if (!(rz[c] >= 0)) rz[c] = 0;
      if (rz[c] > r[c]) rz[c] = r[c];
      if (!(ry[c] >= 0)) ry[c] = 0;
      if (ry[c] > rz[c]) ry[c] = rz[c];
    }
  }

  // captured air mass flow through the throat (for fuel metering)
  measureAirflow() {
    const i = this.columnAt(this.geom.stations.throat);
    const Wr = this.W[0], Wu = this.W[1];
    let m = 0;
    for (let j = PAD; j < this.NY - PAD; j++) {
      if (this.yc(j) > this.geom.yCowl) break;
      const c = j * this.NX + i;
      if (this.type[c] === 0) m += Wr[c] * Math.max(0, Wu[c]) * this.dy;
    }
    return m;
  }

  step() {
    const { U, U0 } = this;
    for (let q = 0; q < NV; q++) U0[q].set(U[q]);
    this.prims(U);
    if (this.steps % 10 === 0) {
      const m = this.measureAirflow();
      this.mdotAir += (m - this.mdotAir) * (this.mdotAir === 0 ? 1 : 0.2);
    }
    const dt = this.maxDt();
    this.rhs(dt, U0, U, false);
    this.fixup(U);
    this.prims(U);
    this.rhs(dt, U, U, true);
    const fc = this.fluidCells;
    for (let q = 0; q < NV; q++) {
      const A = U[q], B = U0[q];
      for (let k = 0; k < fc.length; k++) { const c = fc[k]; A[c] = 0.5 * (A[c] + B[c]); }
    }
    this.fixup(U);
    this.t += dt;
    this.steps++;
    return dt;
  }

  // ---------------------------------------------------------------- diagnostics
  refreshPrims() { this.prims(this.U); }
  columnAt(x) { return PAD + Math.min(this.nx - 1, Math.max(0, Math.floor(x / this.dx))); }
  rowAt(y) { return PAD + Math.min(this.ny - 1, Math.max(0, Math.floor((y - this.geom.y0) / this.dy))); }

  ductMach(x) {
    const i = this.columnAt(x);
    const Wr = this.W[0], Wu = this.W[1], Wv = this.W[2], Wp = this.W[3];
    // the core excludes the boundary layers (cells within 2.5 cells of a wall),
    // which are always partly subsonic once friction is on
    const core = 2.5 * this.dy;
    let num = 0, den = 0, minM = Infinity;
    for (let j = PAD; j < this.NY - PAD; j++) {
      if (this.yc(j) > this.geom.yCowl) break;
      const c = j * this.NX + i;
      if (this.type[c] !== 0) continue;
      const a = Math.sqrt(this.Ga[c] * Wp[c] / Wr[c]);
      const M = Math.hypot(Wu[c], Wv[c]) / a;
      const m = Math.abs(Wr[c] * Wu[c]);
      num += m * M; den += m;
      if (M < minM && this.dist[c] > core) minM = M;
    }
    return { mean: den > 0 ? num / den : 0, min: minM };
  }

  // mass-flux-weighted static temperature across the duct at station x
  ductT(x) {
    const i = this.columnAt(x);
    const Wr = this.W[0], Wu = this.W[1];
    let num = 0, den = 0;
    for (let j = PAD; j < this.NY - PAD; j++) {
      if (this.yc(j) > this.geom.yCowl) break;
      const c = j * this.NX + i;
      if (this.type[c] !== 0) continue;
      const m = Math.abs(Wr[c] * Wu[c]);
      num += m * this.Th[c]; den += m;
    }
    return den > 0 ? num / den : 0;
  }

  maxT(x0, x1) {
    let th = 0;
    const i0 = this.columnAt(x0), i1 = this.columnAt(x1);
    for (let j = PAD; j < this.NY - PAD; j++) {
      if (this.yc(j) > this.geom.yCowl) break;
      for (let i = i0; i <= i1; i++) {
        const c = j * this.NX + i;
        if (this.type[c] === 0) th = Math.max(th, this.Th[c]);
      }
    }
    return th;
  }

  // Pressure force on all walls, x-component (positive = thrust, toward -x)
  pressureThrust() {
    const { NX, NY, type, dy } = this;
    const Wp = this.W[3];
    const pinf = 1 / G0;
    let F = 0;
    for (let j = PAD; j < NY - PAD; j++) {
      const row = j * NX;
      for (let i = PAD; i < NX - PAD - 1; i++) {
        const a = row + i, b = a + 1;
        const fa = type[a] === 0, fb = type[b] === 0;
        if (fa && !fb && type[b] !== 3) F -= (Wp[a] - pinf) * dy;
        else if (!fa && fb && type[a] !== 3) F += (Wp[b] - pinf) * dy;
      }
    }
    return F;
  }

  // Skin-friction drag (x-component, positive = drag) and heat into the walls
  wallLoads() {
    const { wTau, wQ, wTx, wAV, wCell } = this;
    const vol = this.dx * this.dy;
    let drag = 0, heat = 0;
    for (let w = 0; w < wCell.length; w++) {
      const A = wAV[w] * vol;
      drag += wTau[w] * A * wTx[w];
      heat += wQ[w] * A;
    }
    return { drag, heat };
  }
}

export const SCALES = { force: P_SCALE, heat: Q_SCALE };
