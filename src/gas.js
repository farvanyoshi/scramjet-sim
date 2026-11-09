// Gas model: thermally perfect air + transport properties + H2-air chemistry.
//
// Non-dimensionalization (shared with solver.js)
//   density   / rho_inf          (freestream = 1)
//   velocity  / a_inf            (freestream speed of sound = 1, so u_inf = Mach)
//   pressure  / (rho_inf a_inf^2)
//   length    / 1 m
//   time      / (1 m / a_inf)    (= 3.31 ms)
//   "theta" = p / rho = R T / a_inf^2, so T [K] = theta * TS
//
// Thermally perfect air: translation + rotation (cv = 5/2 R) plus the
// vibrational energy of N2 and O2 as harmonic oscillators. That makes the
// specific heat rise and gamma fall from 1.40 at room temperature to about
// 1.29 at 3000 K, like real air (dissociation is handled in the chemistry).

export const R_AIR = 287.05;                       // J/(kg K)
export const ATMOS = {                              // US Standard Atmosphere, 30 km
  altitude_km: 30,
  T: 226.5,                                         // K
  p: 1197,                                          // Pa
  rho: 0.01841,                                     // kg/m^3
};
ATMOS.a = Math.sqrt(1.4 * R_AIR * ATMOS.T);        // 301.7 m/s
export const TS = (ATMOS.a * ATMOS.a) / R_AIR;     // K per unit theta (317.1 K)
export const P_SCALE = ATMOS.rho * ATMOS.a * ATMOS.a;          // Pa per unit pressure
export const Q_SCALE = P_SCALE * ATMOS.a;                       // W/m^2 per unit heat flux
export const MU_SCALE = ATMOS.rho * ATMOS.a * 1.0;             // Pa s per unit viscosity

export const toKelvin = (theta) => theta * TS;
export const fromKelvin = (TK) => TK / TS;
export const toKPa = (p) => (p * P_SCALE) / 1000;
export const toMS = (v) => v * ATMOS.a;

// --- vibrational modes of air ---
const MODES = [
  { Y: 0.7671, R: 296.8, thv: 3395 },   // N2
  { Y: 0.2329, R: 259.8, thv: 2239 },   // O2
].map((m) => ({ w: (m.Y * m.R) / R_AIR, X: m.thv / TS }));

function eExact(th) {
  let e = 2.5 * th;
  for (const m of MODES) { const q = Math.exp(-m.X / th); e += m.w * m.X * q / (1 - q); }
  return e;
}
function cvExact(th) {
  let cv = 2.5;
  for (const m of MODES) {
    const x = m.X / th, q = Math.exp(-x);   // overflow-safe form
    cv += m.w * (x * x * q) / ((1 - q) * (1 - q));
  }
  return cv;
}

// Sutherland's law for air
function muKelvin(TK) {
  return 1.716e-5 * Math.pow(TK / 273.15, 1.5) * (273.15 + 110.4) / (TK + 110.4);
}

// --- lookup tables ---
export const TH_MIN = 0.15;   // 48 K
export const TH_MAX = 30;     // 9500 K
const DTH = 0.005;
const NTH = Math.ceil((TH_MAX) / DTH) + 2;
export const tabGamma = new Float64Array(NTH);   // gamma(theta)
export const tabMu = new Float64Array(NTH);      // mu(theta), non-dim
export const tabCp = new Float64Array(NTH);      // cp/R
export const tabE = new Float64Array(NTH);       // e(theta)
for (let k = 0; k < NTH; k++) {
  const th = Math.max(1e-3, k * DTH);
  const cv = cvExact(th);
  tabGamma[k] = (cv + 1) / cv;
  tabCp[k] = cv + 1;
  tabMu[k] = muKelvin(th * TS) / MU_SCALE;
  tabE[k] = eExact(th);
}
export const INV_DTH = 1 / DTH;

const E_MAX = eExact(TH_MAX);
const DE = 0.002;
const NE = Math.ceil(E_MAX / DE) + 2;
const tabTh = new Float64Array(NE);               // theta(e)
{
  let th = 1e-3;
  for (let k = 0; k < NE; k++) {
    const e = k * DE;
    for (let it = 0; it < 30; it++) {
      const f = eExact(th) - e;
      th -= f / cvExact(th);
      if (th < 1e-3) th = 1e-3;
      if (Math.abs(f) < 1e-12) break;
    }
    tabTh[k] = th;
  }
}
export const INV_DE = 1 / DE;
export const E_MIN = eExact(TH_MIN);

// theta from specific internal energy (table + linear interpolation)
export function thetaFromE(e) {
  let f = e * INV_DE;
  if (f < 0) f = 0;
  const k = f | 0;
  if (k >= NE - 1) return TH_MAX;
  const t = f - k;
  return tabTh[k] + t * (tabTh[k + 1] - tabTh[k]);
}

// table lookup on theta (linear)
export function lookup(tab, th) {
  let f = th * INV_DTH;
  if (f < 0) f = 0;
  const k = f | 0;
  if (k >= NTH - 1) return tab[NTH - 1];
  const t = f - k;
  return tab[k] + t * (tab[k + 1] - tab[k]);
}

export const eOf = (th) => lookup(tabE, th);
export const gammaOf = (th) => lookup(tabGamma, th);
export const muOf = (th) => lookup(tabMu, th);
export const cpOf = (th) => lookup(tabCp, th);

// --- transport constants ---
export const PR = 0.72;       // Prandtl number
export const PR_T = 0.9;      // turbulent Prandtl number
export const SC = 0.22;       // H2-air Schmidt number (hydrogen diffuses fast)
export const SC_T = 0.7;      // turbulent Schmidt number

// --- Spalart-Allmaras constants ---
export const SA = {
  cb1: 0.1355, cb2: 0.622, sigma: 2 / 3, kappa: 0.41,
  cw2: 0.3, cw3: 2, cv1: 7.1,
};
SA.cw1 = SA.cb1 / (SA.kappa * SA.kappa) + (1 + SA.cb2) / SA.sigma;

// --- hydrogen-air chemistry: one global step  2 H2 + O2 -> 2 H2O ---
export const CHEM = {
  fst: 0.02916,                     // stoichiometric fuel/air mass ratio (1/34.3)
  YO2air: 0.2329,                   // O2 mass fraction in air
  sO2: 7.937,                       // kg O2 burned per kg H2
  Q: 120e6 / (ATMOS.a * ATMOS.a),   // lower heating value of H2, 120 MJ/kg, non-dim
  A: 0,                             // pre-exponential factor (set below)
  Ta: 0,                            // activation temperature, K
  Tdiss: 2600,                      // K: above this, products stay dissociated
  dTdiss: 300,                      // K: width of that taper
};
CHEM.Ystoich = CHEM.fst / (1 + CHEM.fst);
// Calibrated with test/run.js so autoignition sets in around 950-1000 K
// within the ~0.5 ms the air spends in the combustor.
CHEM.Ta = 20000;
CHEM.A = 7.0e9;

// --- engine walls ---
export const WALL = {
  T: 1000,                          // K, fuel-cooled wall temperature
  B: 5.0,                           // log-law constant
};
