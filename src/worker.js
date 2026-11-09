// Runs the solver off the main thread so the page stays responsive.
// The main thread sends commands; the worker steps the solver as fast as it
// can and sends back a snapshot of the flow whenever the page asks for one.
import { buildEngine } from './geometry.js';
import { Solver, GAMMA, ATMOS, FUEL, toKelvin } from './solver.js';

let solver = null;
let running = true;
let wantFrame = true;
let mach = 6, phi = 0.6;
let stepTimes = [];

function init(nx) {
  const geom = buildEngine();
  const ny = Math.round(nx * geom.Ly / geom.Lx);
  solver = new Solver(geom, nx, ny);
  solver.mach = mach;
  solver.phi = phi;
  solver.reset();
  postMessage({
    type: 'grid',
    nx, ny, NX: solver.NX, NY: solver.NY,
    cellType: solver.type.slice(),
    geom: { Lx: geom.Lx, Ly: geom.Ly, y0: geom.y0, solids: geom.solids, sections: geom.sections, stations: geom.stations, injector: geom.injector, yCowl: geom.yCowl },
    fuel: { Ystoich: FUEL.Ystoich, Q: FUEL.Q, TignK: toKelvin(FUEL.Tign) },
  });
  wantFrame = true;
}

function diagnostics() {
  const s = solver;
  const st = s.geom.stations;
  const throat = s.ductMach(st.throat);
  const iso = s.ductMach(st.isolator);
  const exit = s.ductMach(st.combustorEnd);
  const TisoK = toKelvin(s.maxT(st.throat, s.geom.injector.x0 - 0.02));
  const TcombK = toKelvin(s.maxT(s.geom.injector.x0, st.combustorEnd));
  const forceKN = s.axialForce() * GAMMA * ATMOS.p / 1000;
  const eta = s.injRate > 1e-9 ? Math.min(1, s.burnRate / s.injRate) : 0;
  const flowThrough = s.geom.Lx / Math.max(1, s.mach);
  return {
    t: s.t,
    tPhysMs: s.t * (1 / ATMOS.a) * 1000,   // 1 length unit = 1 m
    steps: s.steps,
    settling: s.t < 1.5 * flowThrough,
    throatM: throat.mean,
    isoM: iso.mean,
    isoMin: iso.min,
    exitM: exit.mean,
    exitMin: exit.min,
    TisoK,
    TcombK,
    forceKN,
    eta,
    mach: s.mach,
    phi: s.phi,
  };
}

function sendFrame() {
  const s = solver;
  s.refreshPrims();
  const N = s.N;
  const out = {};
  const names = ['r', 'u', 'v', 'p', 'y'];
  const buffers = [];
  for (let q = 0; q < 5; q++) {
    const a = new Float32Array(N);
    a.set(s.W[q]);
    out[names[q]] = a;
    buffers.push(a.buffer);
  }
  const hr = s.hr.slice();
  buffers.push(hr.buffer);
  out.hr = hr;
  const now = performance.now();
  stepTimes = stepTimes.filter((x) => now - x.t < 1000);
  const sps = stepTimes.length > 1 ? (stepTimes[stepTimes.length - 1].n - stepTimes[0].n) / ((now - stepTimes[0].t) / 1000) : 0;
  postMessage({ type: 'frame', fields: out, diag: { ...diagnostics(), stepsPerSec: sps, running } }, buffers);
  wantFrame = false;
}

function loop() {
  if (solver && running) {
    const t0 = performance.now();
    while (performance.now() - t0 < 30) solver.step();
    stepTimes.push({ t: performance.now(), n: solver.steps });
  }
  if (solver && wantFrame) sendFrame();
  setTimeout(loop, running ? 0 : 50);
}

onmessage = (e) => {
  const m = e.data;
  switch (m.type) {
    case 'init': init(m.nx); break;
    case 'set':
      if (m.mach !== undefined) { mach = m.mach; if (solver) solver.mach = mach; }
      if (m.phi !== undefined) { phi = m.phi; if (solver) solver.phi = phi; }
      break;
    case 'reset': if (solver) { solver.reset(); stepTimes = []; } wantFrame = true; break;
    case 'run': running = m.on; wantFrame = true; break;
    case 'step': if (solver) { for (let k = 0; k < (m.n || 1); k++) solver.step(); } wantFrame = true; break;
    case 'ack': wantFrame = true; break;
  }
};

loop();
