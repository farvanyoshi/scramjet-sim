// Headless check: run the solver for a while and print the readouts.
// Usage: node test/run.js [mach] [phi] [steps] [nx]
import { buildEngine } from '../src/geometry.js';
import { Solver, toKelvin, toKPa, ATMOS, GAMMA } from '../src/solver.js';

const mach = +(process.argv[2] ?? 6);
const phi = +(process.argv[3] ?? 0.4);
const steps = +(process.argv[4] ?? 6000);
const nx = +(process.argv[5] ?? 320);
const geom = buildEngine();
const ny = Math.round(nx * geom.Ly / geom.Lx);
const s = new Solver(geom, nx, ny);
s.mach = mach; s.phi = 0; s.reset();
const t0 = performance.now();
let bad = false;
for (let n = 1; n <= steps; n++) {
  if (n === Math.floor(steps / 3)) s.phi = phi;
  s.step();
  if (n % Math.floor(steps / 6) === 0 || n === steps) {
    s.refreshPrims();
    let nan = 0;
    for (const c of s.fluidCells) if (!Number.isFinite(s.U[3][c])) nan++;
    if (nan) bad = true;
    const st = geom.stations;
    const thr = s.ductMach(st.throat), iso = s.ductMach(st.isolator), ce = s.ductMach(st.combustorEnd);
    const F = s.axialForce() * GAMMA * ATMOS.p / 1000;
    console.log(`step ${n} t=${s.t.toFixed(3)} throatM=${thr.mean.toFixed(2)}(min ${thr.min.toFixed(2)}) isoM=${iso.mean.toFixed(2)} combExitM=${ce.mean.toFixed(2)}(min ${ce.min.toFixed(2)}) ` +
      `Tiso=${toKelvin(s.maxT(1.7, 2.25)).toFixed(0)}K Tmax=${toKelvin(s.maxT(2.3, 3.0)).toFixed(0)}K force=${F.toFixed(2)}kN/m eta=${(s.burnRate / Math.max(1e-9, s.injRate)).toFixed(2)} nan=${nan}`);
  }
}
const ms = (performance.now() - t0) / steps;
console.log(`grid ${nx}x${ny}, ${ms.toFixed(2)} ms/step, ghosts=${s.gCell.length}`);
if (bad) { console.error('NaN detected'); process.exit(1); }
