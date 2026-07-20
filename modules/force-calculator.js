/**
 * Force Calculator Module
 *
 * Van der Waals (Lennard-Jones), electrostatics (Coulomb), hydrogen bonding,
 * and harmonic bonds.
 *
 * STABILITY FIXES (so the MD integrator doesn't explode):
 *  - 1-2 (bonded) and 1-3 (angle) neighbor pairs are excluded from the
 *    non-bonded / H-bond terms. Previously covalently bonded atoms ~1.5 Å
 *    apart were hit with enormous LJ repulsion, blowing up the simulation.
 *  - A spatial cell list replaces the O(n^2) pair loop so large structures
 *    scale roughly linearly instead of freezing.
 *  - Per-atom force magnitude is capped so any residual close contact can't
 *    produce an infinite kick.
 */

export class ForceCalculator {
  constructor(parameters = {}) {
    this.lennardJones = parameters.lennardJones || {
      sigma: {
        'C-C': 3.5, 'C-O': 3.18, 'C-N': 3.25, 'C-S': 3.7,
        'O-O': 3.0, 'N-O': 3.1, 'O-S': 3.54,
        'N-N': 3.25, 'N-S': 3.58,
        'S-S': 3.8, 'C-H': 3.08, 'H-O': 2.96,
        'H-N': 3.1, 'H-S': 3.6,
        'H-H': 2.6,
      },
      epsilon: {
        'C-C': 0.066, 'C-O': 0.101, 'C-N': 0.077, 'C-S': 0.25,
        'O-O': 0.152, 'N-O': 0.134, 'O-S': 0.274,
        'N-N': 0.11, 'N-S': 0.297,
        'S-S': 0.557, 'C-H': 0.03, 'H-O': 0.041,
        'H-N': 0.032, 'H-S': 0.125,
        'H-H': 0.015,
      },
      cutoff: 12.0,
    };

    this.electrostatic = parameters.electrostatic || {
      coulombConstant: 332.06,
      dielectric: parameters.dielectric || 4.0,
      cutoff: 12.0,
    };

    this.hydrogenBond = parameters.hydrogenBond || {
      maxDistance: 3.5, minAngle: 120, energy: -5.0, energyRange: 0.5,
    };

    this.harmonicBond = parameters.harmonicBond || {
      forceConstant: 700,
      defaultLengths: {
        'C-C': 1.54, 'C-O': 1.43, 'C-N': 1.47, 'C-H': 1.09,
        'H-O': 0.96, 'H-N': 1.01,
      },
    };

    // Largest short-range cutoff used for the neighbor search.
    this.pairCutoff = Math.max(this.lennardJones.cutoff, this.electrostatic.cutoff);
    // Cap on per-atom force magnitude (kcal/mol/Å) to keep the integrator sane.
    this.maxForce = parameters.maxForce || 500;
  }

  /* ---------------------------------------------------------------------- */

  /** Build a set of excluded (1-2 and 1-3) index pairs from the bond list. */
  buildExclusions(atoms, bonds) {
    const n = atoms.length;
    const set = new Set();
    const key = (i, j) => (i < j ? i * n + j : j * n + i);
    const adj = Array.from({ length: n }, () => []);

    if (bonds) {
      for (const b of bonds) {
        const i = b.atom1Index, j = b.atom2Index;
        if (i == null || j == null || i < 0 || j < 0 || i >= n || j >= n) continue;
        adj[i].push(j);
        adj[j].push(i);
        set.add(key(i, j)); // 1-2 (directly bonded)
      }
      // 1-3 (share a common neighbor → angle pair)
      for (let i = 0; i < n; i++) {
        for (const j of adj[i]) {
          for (const k of adj[j]) {
            if (k !== i) set.add(key(i, k));
          }
        }
      }
    }
    return { has: (i, j) => set.has(key(i, j)) };
  }

  /**
   * Iterate unique atom pairs within `cutoff` using a uniform spatial grid.
   * Calls cb(i, j, r, dx, dy, dz) with i < j and dx = xj - xi.
   */
  _forEachNeighborPair(atoms, cutoff, cb) {
    const n = atoms.length;
    if (n === 0) return;

    let minx = Infinity, miny = Infinity, minz = Infinity;
    for (const a of atoms) {
      const p = a.position;
      if (p.x < minx) minx = p.x;
      if (p.y < miny) miny = p.y;
      if (p.z < minz) minz = p.z;
    }
    const cs = Math.max(cutoff, 1e-3);
    const cells = new Map();
    const coords = new Array(n);

    for (let i = 0; i < n; i++) {
      const p = atoms[i].position;
      const cx = Math.floor((p.x - minx) / cs);
      const cy = Math.floor((p.y - miny) / cs);
      const cz = Math.floor((p.z - minz) / cs);
      coords[i] = [cx, cy, cz];
      const k = cx + ',' + cy + ',' + cz;
      let arr = cells.get(k);
      if (!arr) { arr = []; cells.set(k, arr); }
      arr.push(i);
    }

    const cutoff2 = cutoff * cutoff;
    for (let i = 0; i < n; i++) {
      const [cx, cy, cz] = coords[i];
      const pi = atoms[i].position;
      for (let dx = -1; dx <= 1; dx++) {
        for (let dy = -1; dy <= 1; dy++) {
          for (let dz = -1; dz <= 1; dz++) {
            const arr = cells.get((cx + dx) + ',' + (cy + dy) + ',' + (cz + dz));
            if (!arr) continue;
            for (const j of arr) {
              if (j <= i) continue;
              const pj = atoms[j].position;
              const ddx = pj.x - pi.x, ddy = pj.y - pi.y, ddz = pj.z - pi.z;
              const r2 = ddx * ddx + ddy * ddy + ddz * ddz;
              if (r2 <= cutoff2) cb(i, j, Math.sqrt(r2), ddx, ddy, ddz);
            }
          }
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */

  calculateAllForces(atoms, bonds) {
    const forces = atoms.map(() => ({ x: 0, y: 0, z: 0 }));
    const excluded = this.buildExclusions(atoms, bonds);

    // Single short-range pass: LJ + electrostatic + H-bond.
    this._forEachNeighborPair(atoms, this.pairCutoff, (i, j, r, dx, dy, dz) => {
      if (r < 0.1 || excluded.has(i, j)) return;
      const ai = atoms[i], aj = atoms[j];
      let fx = 0, fy = 0, fz = 0;

      if (r < this.lennardJones.cutoff) {
        const c = this.calculateLJForce(ai, aj, r) / r;
        fx += c * dx; fy += c * dy; fz += c * dz;
      }
      if (ai.charge && aj.charge && r < this.electrostatic.cutoff) {
        const c = this.calculateElectrostaticForce(ai, aj, r) / r;
        fx += c * dx; fy += c * dy; fz += c * dz;
      }
      if (this.isHydrogenBond(ai, aj) && r < this.hydrogenBond.maxDistance + 1) {
        const r_opt = 2.8;
        const fm = -50 * (r - r_opt) * Math.exp(-Math.pow(r - r_opt, 2) / 0.5);
        if (Math.abs(fm) > 0.01) {
          const c = fm / r;
          fx += c * dx; fy += c * dy; fz += c * dz;
        }
      }

      forces[i].x += fx; forces[i].y += fy; forces[i].z += fz;
      forces[j].x -= fx; forces[j].y -= fy; forces[j].z -= fz;
    });

    if (bonds && bonds.length > 0) this.calculateBondForces(atoms, bonds, forces);
    this.capForces(forces);
    return forces;
  }

  /** Clamp each atom's force magnitude to maxForce. */
  capForces(forces) {
    const max = this.maxForce;
    for (const f of forces) {
      const m = Math.sqrt(f.x * f.x + f.y * f.y + f.z * f.z);
      if (m > max) {
        const s = max / m;
        f.x *= s; f.y *= s; f.z *= s;
      }
    }
  }

  calculateLJForce(atom_i, atom_j, r) {
    if (r < 0.1) return 0;
    const key = this.getLJKey(atom_i.element, atom_j.element);
    const sigma = this.lennardJones.sigma[key];
    const epsilon = this.lennardJones.epsilon[key];
    if (!sigma || !epsilon) return 0;
    const sr6 = Math.pow(sigma / r, 6);
    const sr12 = sr6 * sr6;
    return (24 * epsilon * (2 * sr12 - sr6)) / r;
  }

  getLJKey(elem1, elem2) {
    const pair = [elem1, elem2].sort().join('-');
    if (this.lennardJones.sigma[pair]) return pair;
    const reverse = [elem2, elem1].join('-');
    if (this.lennardJones.sigma[reverse]) return reverse;
    return pair;
  }

  calculateElectrostaticForce(atom_i, atom_j, r) {
    if (r < 0.1) return 0;
    const q1 = atom_i.charge || 0;
    const q2 = atom_j.charge || 0;
    const k = this.electrostatic.coulombConstant;
    const epsilon = this.electrostatic.dielectric;
    return (k * q1 * q2) / (epsilon * r * r);
  }

  calculateBondForces(atoms, bonds, forces) {
    bonds.forEach((bond) => {
      const a_i = atoms[bond.atom1Index];
      const a_j = atoms[bond.atom2Index];
      if (!a_i || !a_j) return;

      const dx = a_j.position.x - a_i.position.x;
      const dy = a_j.position.y - a_i.position.y;
      const dz = a_j.position.z - a_i.position.z;
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (r < 0.01) return;

      const r0 = bond.length || this.getDefaultBondLength(a_i.element, a_j.element) || 1.5;
      const k = bond.forceConstant || this.harmonicBond.forceConstant;
      const forceMag = -k * (r - r0);

      const fx = (forceMag / r) * dx;
      const fy = (forceMag / r) * dy;
      const fz = (forceMag / r) * dz;

      forces[bond.atom1Index].x -= fx;
      forces[bond.atom1Index].y -= fy;
      forces[bond.atom1Index].z -= fz;
      forces[bond.atom2Index].x += fx;
      forces[bond.atom2Index].y += fy;
      forces[bond.atom2Index].z += fz;
    });
  }

  getDefaultBondLength(elem1, elem2) {
    const key = [elem1, elem2].sort().join('-');
    return this.harmonicBond.defaultLengths[key] || null;
  }

  isHydrogenBond(atom1, atom2) {
    const set = ['N', 'O'];
    return set.includes(atom1.element) && set.includes(atom2.element);
  }

  /** Kept for API compatibility; exclusions now handle bonded-pair skipping. */
  areBonded() {
    return false;
  }

  calculatePotentialEnergy(atoms, bonds) {
    let PE = 0;
    const excluded = this.buildExclusions(atoms, bonds);

    this._forEachNeighborPair(atoms, this.pairCutoff, (i, j, r) => {
      if (r < 0.1 || excluded.has(i, j)) return;
      const ai = atoms[i], aj = atoms[j];
      if (r < this.lennardJones.cutoff) PE += this.calculateLJPotential(ai, aj, r);
      if (ai.charge && aj.charge && r < this.electrostatic.cutoff) {
        PE += this.calculateElectrostaticPotential(ai, aj, r);
      }
    });

    if (bonds) {
      bonds.forEach((bond) => {
        const a_i = atoms[bond.atom1Index];
        const a_j = atoms[bond.atom2Index];
        if (!a_i || !a_j) return;
        const dx = a_j.position.x - a_i.position.x;
        const dy = a_j.position.y - a_i.position.y;
        const dz = a_j.position.z - a_i.position.z;
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        const r0 = bond.length || this.getDefaultBondLength(a_i.element, a_j.element) || 1.5;
        const k = bond.forceConstant || this.harmonicBond.forceConstant;
        PE += 0.5 * k * (r - r0) ** 2;
      });
    }

    return PE;
  }

  calculateLJPotential(atom_i, atom_j, r) {
    const key = this.getLJKey(atom_i.element, atom_j.element);
    const sigma = this.lennardJones.sigma[key];
    const epsilon = this.lennardJones.epsilon[key];
    if (!sigma || !epsilon) return 0;
    const sr6 = Math.pow(sigma / r, 6);
    const sr12 = sr6 * sr6;
    return 4 * epsilon * (sr12 - sr6);
  }

  calculateElectrostaticPotential(atom_i, atom_j, r) {
    const q1 = atom_i.charge || 0;
    const q2 = atom_j.charge || 0;
    const k = this.electrostatic.coulombConstant;
    const epsilon = this.electrostatic.dielectric;
    return (k * q1 * q2) / (epsilon * r);
  }
}

export default ForceCalculator;
