// Headless check: run the solver and print the readouts as the flow develops.
// Usage: node test/run.js [mach] [phi] [steps] [nx] [model]
//   model: turbulent (default) | inviscid
import { buildEngine } from '../src/geometry.js';
import { Solver, toKelvin, SCALES } from '../src/solver.js';

const mach = +(process.argv[2] ?? 6);
const phi = +(process.argv[3] ?? 0.6);
const steps = +(process.argv[4] ?? 4000);
const nx = +(process.argv[5] ?? 320);
const model = process.argv[6] ?? 'turbulent';
const geom = buildEngine();
const ny = Math.round(nx * geom.Ly / geom.Lx);
const s = new Solver(geom, nx, ny, model);
s.mach = mach; s.phi = 0; s.reset();
const t0 = performance.now();
let bad = false;
const every = Math.max(1, Math.floor(steps / 6));
for (let n = 1; n <= steps; n++) {
  if (n === Math.floor(steps / 3)) s.phi = phi;
  s.step();
  if (n % every === 0 || n === steps) {
    s.refreshPrims();
    let nan = 0;
    for (const c of s.fluidCells) if (!Number.isFinite(s.U[3][c]) || !Number.isFinite(s.U[4][c])) nan++;
    if (nan) bad = true;
    const st = geom.stations;
    const thr = s.ductMach(st.throat), iso = s.ductMach(st.isolator), ce = s.ductMach(st.combustorEnd);
    const Fp = s.pressureThrust() * SCALES.force / 1000;
    const wl = s.wallLoads();
    const Df = wl.drag * SCALES.force / 1000;
    const Q = wl.heat * SCALES.heat / 1e6;
    let mtMax = 0;
    for (const c of s.fluidCells) mtMax = Math.max(mtMax, s.Mt[c] / s.Mu[c]);
    console.log(`step ${n} t=${s.t.toFixed(3)} thrM=${thr.mean.toFixed(2)}(min ${thr.min.toFixed(2)}) isoM=${iso.mean.toFixed(2)}(min ${iso.min.toFixed(2)}) exitM=${ce.mean.toFixed(2)}(min ${ce.min.toFixed(2)}) ` +
      `Tiso=${toKelvin(s.ductT(2.25)).toFixed(0)}K(max ${toKelvin(s.maxT(1.7, 2.28)).toFixed(0)}) Tcomb=${toKelvin(s.maxT(2.3, 3.0)).toFixed(0)}K ` +
      `Fp=${Fp.toFixed(2)} Df=${Df.toFixed(2)} net=${(Fp - Df).toFixed(2)}kN/m heat=${Q.toFixed(2)}MW/m ` +
      `eta=${(s.burnRate / Math.max(1e-12, s.injRate)).toFixed(2)} mair=${s.mdotAir.toFixed(3)} mut/mu=${mtMax.toFixed(0)} nan=${nan}`);
  }
}
const ms = (performance.now() - t0) / steps;
console.log(`grid ${nx}x${ny} model=${model}: ${ms.toFixed(2)} ms/step, wall cells=${s.wCell.length}, injector cells=${s.injSides.map((x) => x.cells.length).join('+')}`);
if (bad) { console.error('NaN detected'); process.exit(1); }
