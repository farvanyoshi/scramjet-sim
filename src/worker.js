// Runs the solver off the main thread so the page stays responsive.
// The main thread sends commands; the worker steps the solver as fast as it
// can and sends back a snapshot of the flow whenever the page asks for one.
import { buildEngine } from './geometry.js';
import { Solver, SCALES, toKelvin, CHEM, ATMOS } from './solver.js';

let solver = null;
let running = true;
let wantFrame = true;
let mach = 6, phi = 0.6, model = 'turbulent';
let stepTimes = [];

function init(nx) {
  const geom = buildEngine();
  const ny = Math.round(nx * geom.Ly / geom.Lx);
  solver = new Solver(geom, nx, ny, model);
  solver.mach = mach;
  solver.phi = phi;
  solver.reset();
  stepTimes = [];
  postMessage({
    type: 'grid',
    nx, ny, NX: solver.NX, NY: solver.NY,
    cellType: solver.type.slice(),
    walls: {
      x: Float32Array.from(solver.wPx), y: Float32Array.from(solver.wPy),
      nx: Float32Array.from(solver.wNx), ny: Float32Array.from(solver.wNy),
    },
    geom: {
      Lx: geom.Lx, Ly: geom.Ly, y0: geom.y0, solids: geom.solids, sections: geom.sections,
      stations: geom.stations, injector: geom.injector, yCowl: geom.yCowl,
    },
    fuel: { fst: CHEM.fst },
  });
  wantFrame = true;
}

function diagnostics() {
  const s = solver;
  const st = s.geom.stations;
  const throat = s.ductMach(st.throat);
  const iso = s.ductMach(st.isolator);
  const exit = s.ductMach(st.combustorEnd);
  const Fp = (s.pressureThrust() * SCALES.force) / 1000;           // kN/m
  const wl = s.wallLoads();
  const Df = s.model === 'inviscid' ? 0 : (wl.drag * SCALES.force) / 1000;
  const heat = s.model === 'inviscid' ? 0 : (wl.heat * SCALES.heat) / 1e6; // MW/m
  const eta = s.injRate > 1e-9 ? Math.min(1, s.burnRate / s.injRate) : 0;
  const flowThrough = s.geom.Lx / Math.max(1, s.mach);
  return {
    t: s.t,
    tPhysMs: (s.t / ATMOS.a) * 1000,
    steps: s.steps,
    settling: s.t < 1.2 * flowThrough,
    throatM: throat.mean,
    isoM: iso.mean,
    isoMin: iso.min,
    exitM: exit.mean,
    exitMin: exit.min,
    TcoreK: toKelvin(s.ductT(s.geom.injector.x0 - 0.05)),
    TcombK: toKelvin(s.maxT(s.geom.injector.x0, st.combustorEnd)),
    pressureKN: Fp,
    frictionKN: Df,
    forceKN: Fp - Df,
    heatMW: heat,
    eta,
    mach: s.mach,
    phi: s.phi,
    model: s.model,
  };
}

function sendFrame() {
  const s = solver;
  s.refreshPrims();
  const N = s.N;
  const pack = (src) => { const a = new Float32Array(N); a.set(src); return a; };
  const fields = {
    r: pack(s.W[0]), u: pack(s.W[1]), v: pack(s.W[2]), p: pack(s.W[3]),
    z: pack(s.W[5]), y: pack(s.W[6]), th: pack(s.Th), ga: pack(s.Ga),
    mt: new Float32Array(N),
    hr: s.hr.slice(),
    wq: s.wQ.slice(),
  };
  for (let c = 0; c < N; c++) fields.mt[c] = s.Mu[c] > 0 ? s.Mt[c] / s.Mu[c] : 0;
  const buffers = Object.values(fields).map((a) => a.buffer);
  const now = performance.now();
  stepTimes = stepTimes.filter((x) => now - x.t < 1500);
  const sps = stepTimes.length > 1
    ? (stepTimes[stepTimes.length - 1].n - stepTimes[0].n) / ((stepTimes[stepTimes.length - 1].t - stepTimes[0].t) / 1000)
    : 0;
  postMessage({ type: 'frame', fields, diag: { ...diagnostics(), stepsPerSec: sps, running } }, buffers);
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
      if (m.model !== undefined) { model = m.model; if (solver) solver.setModel(model); }
      wantFrame = true;
      break;
    case 'reset': if (solver) { solver.reset(); stepTimes = []; } wantFrame = true; break;
    case 'run': running = m.on; wantFrame = true; break;
    case 'ack': wantFrame = true; break;
  }
};

loop();
