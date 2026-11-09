// Scramjet engine shape, in simulation length units (1 unit = 1 m for the
// physical readouts). Air flows left to right. The engine is 2D: think of a
// slice through the middle of a wide, flat engine.
//
//   body (bottom): forebody ramps -> isolator floor -> combustor -> nozzle ramp
//   cowl (top):    a thin wedge that closes the duct from above
//
// Everything below the body line and inside the cowl polygon is solid wall.

const deg = Math.PI / 180;

export function buildEngine() {
  const Lx = 4.0;   // domain length
  const y0 = 0.0;   // domain bottom
  const Ly = 0.9;   // domain height

  // --- body (lower wall) ---
  const xLE = 0.25;               // start of first ramp
  const yLE = 0.12;               // flat plate height before the ramps
  const xR2 = 1.0;                // start of second ramp
  const xThroat = 1.6;            // end of compression = start of isolator
  const xComb = 2.3;              // isolator -> combustor
  const xNoz = 3.0;               // combustor -> nozzle
  const ramp1 = 9 * deg;          // first ramp angle
  const ramp2 = 18 * deg;         // second ramp (total turning)
  const combDiv = 3 * deg;        // combustor floor drops slightly (avoids choking)
  const nozDiv = 16 * deg;        // nozzle expansion ramp

  const yR2 = yLE + (xR2 - xLE) * Math.tan(ramp1);
  const yThroat = yR2 + (xThroat - xR2) * Math.tan(ramp2);
  const yNoz = yThroat - (xNoz - xComb) * Math.tan(combDiv);
  const xEnd = Lx + 0.2;
  const yEnd = yNoz - (xEnd - xNoz) * Math.tan(nozDiv);

  const body = [
    [-0.2, -0.2],
    [-0.2, yLE],
    [xLE, yLE],
    [xR2, yR2],
    [xThroat, yThroat],
    [xComb, yThroat],
    [xNoz, yNoz],
    [xEnd, yEnd],
    [xEnd, -0.2],
  ];

  // --- cowl (upper wall of the duct) ---
  const xLip = 1.55;
  const yCowl = 0.60;             // underside of the cowl (flat)
  const cowlAngle = 7 * deg;      // outer wedge angle at the lip
  const cowlThick = 0.06;
  const xCowlTop = xLip + cowlThick / Math.tan(cowlAngle);
  const xCowlEnd = 3.05;
  const cowl = [
    [xLip, yCowl],
    [xCowlTop, yCowl + cowlThick],
    [2.85, yCowl + cowlThick],
    [xCowlEnd, yCowl + 0.015],
    [xCowlEnd, yCowl],
  ];

  // Lower-wall height at x (for zone tests and labels)
  function yBody(x) {
    for (let k = 1; k < body.length - 1; k++) {
      const [xa, ya] = body[k], [xb, yb] = body[k + 1];
      if (x >= xa && x <= xb) return ya + (yb - ya) * (x - xa) / (xb - xa);
    }
    return yLE;
  }

  return {
    Lx, Ly, y0,
    solids: [
      { name: 'body', pts: body },
      { name: 'cowl', pts: cowl },
    ],
    yBody,
    yCowl,
    xLip,
    // Fuel (hydrogen) is injected across the duct in this x-range.
    injector: { x0: 2.30, x1: 2.36 },
    // Stations used for the readouts
    stations: { throat: xThroat + 0.05, isolator: 1.95, combustorEnd: xNoz - 0.02 },
    sections: [
      { name: 'Inlet', x0: xLE, x1: xThroat },
      { name: 'Isolator', x0: xThroat, x1: xComb },
      { name: 'Combustor', x0: xComb, x1: xNoz },
      { name: 'Nozzle', x0: xNoz, x1: Lx },
    ],
    angles: { ramp1: 9, ramp2: 18, cowl: 7, nozzle: 16 },
  };
}
