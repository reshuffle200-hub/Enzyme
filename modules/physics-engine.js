/**
 * Molecular Dynamics Simulation Engine
 *
 * Velocity Verlet integration with Lennard-Jones, electrostatics, H-bonds,
 * harmonic bonds, and a Langevin thermostat. Optional trajectory recording.
 *
 * FIX: the constructor previously called structure.getAllAtoms()/getAllBonds(),
 * which don't exist on the plain object returned by PDBFetcher. It now reads
 * structure.atoms / structure.bonds directly.
 */

import { ForceCalculator } from './force-calculator.js';

export class PhysicsEngine {
  constructor(structure, options = {}) {
    this.structure = structure;
    // FIX: read the arrays directly rather than calling non-existent methods.
    this.atoms = structure.atoms || [];
    this.bonds = structure.bonds || [];

    this.timestep = options.timestep || 0.001;
    this.temperature = options.temperature || 300;
    this.frictionCoeff = options.frictionCoeff || 0.1;
    this.constrainHeavyAtoms = options.constrainHeavyAtoms || false;
    this.constrainedResidues = options.constrainedResidues || [];
    this.maxVelocity = options.maxVelocity || 50; // Å/ps clamp, prevents runaway

    this.recordTrajectory = options.recordTrajectory !== false;
    this.trajectoryInterval = options.trajectoryInterval || 100;
    this.trajectory = [];
    this.stepCount = 0;

    this.forceCalc = new ForceCalculator(options.forceParameters || {});

    this.isRunning = false;
    this.isPaused = false;
    this.energyStats = { kinetic: 0, potential: 0, total: 0, history: [] };

    this.initializeVelocities();

    this.onStepComplete = null;
    this.onEnergyUpdate = null;
  }

  initializeVelocities() {
    const k_B = 8.314e-3;
    this.atoms.forEach((atom) => {
      if (!atom.velocity) {
        atom.velocity = { x: 0, y: 0, z: 0 };
        atom.force = { x: 0, y: 0, z: 0 };
        atom.mass = atom.mass || this.getAtomicMass(atom.element);
      }
      const sigma = Math.sqrt((k_B * this.temperature) / atom.mass);
      atom.velocity.x = this.randomGaussian(0, sigma);
      atom.velocity.y = this.randomGaussian(0, sigma);
      atom.velocity.z = this.randomGaussian(0, sigma);
    });
  }

  randomGaussian(mean, sigma) {
    let u1, u2;
    do {
      u1 = Math.random();
      u2 = Math.random();
    } while (u1 <= 1e-6);
    const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    return mean + sigma * z;
  }

  getAtomicMass(element) {
    const masses = {
      H: 1.008, C: 12.011, N: 14.007, O: 15.999,
      S: 32.065, P: 30.974, Cl: 35.45, F: 18.998,
      Br: 79.904, I: 126.9,
    };
    return masses[element] || 12.0;
  }

  step() {
    if (this.isPaused) return;

    this.atoms.forEach((atom) => {
      atom.force.x = 0;
      atom.force.y = 0;
      atom.force.z = 0;
    });

    const forces = this.forceCalc.calculateAllForces(this.atoms, this.bonds);
    this.atoms.forEach((atom, idx) => {
      const f = forces[idx];
      atom.force.x = f.x;
      atom.force.y = f.y;
      atom.force.z = f.z;
    });

    this.velocityVerletStep();
    this.updateEnergies();
    this.applyTemperatureControl();

    if (this.recordTrajectory && this.stepCount % this.trajectoryInterval === 0) {
      this.saveSnapshot();
    }

    this.stepCount++;

    if (this.onStepComplete) {
      this.onStepComplete({
        step: this.stepCount,
        energy: this.energyStats.total,
        kinetic: this.energyStats.kinetic,
        potential: this.energyStats.potential,
        temperature: this.calculateInstantaneousTemperature(),
      });
    }
  }

  velocityVerletStep() {
    const dt = this.timestep;
    const dt2half = 0.5 * dt * dt;

    this.atoms.forEach((atom) => {
      if (this.isConstrained(atom)) return;

      const ax = atom.force.x / atom.mass;
      const ay = atom.force.y / atom.mass;
      const az = atom.force.z / atom.mass;

      atom.position.x += atom.velocity.x * dt + ax * dt2half;
      atom.position.y += atom.velocity.y * dt + ay * dt2half;
      atom.position.z += atom.velocity.z * dt + az * dt2half;

      if (!atom.velocityOld) atom.velocityOld = { x: 0, y: 0, z: 0 };
      atom.velocityOld.x = atom.velocity.x;
      atom.velocityOld.y = atom.velocity.y;
      atom.velocityOld.z = atom.velocity.z;

      atom.velocity.x += ax * dt * 0.5;
      atom.velocity.y += ay * dt * 0.5;
      atom.velocity.z += az * dt * 0.5;
    });

    const forces = this.forceCalc.calculateAllForces(this.atoms, this.bonds);

    this.atoms.forEach((atom, idx) => {
      if (this.isConstrained(atom)) return;
      const f = forces[idx];
      const ax = f.x / atom.mass;
      const ay = f.y / atom.mass;
      const az = f.z / atom.mass;
      atom.velocity.x += ax * dt * 0.5;
      atom.velocity.y += ay * dt * 0.5;
      atom.velocity.z += az * dt * 0.5;
    });

    this.applyLangevinThermostat();
  }

  applyLangevinThermostat() {
    const gamma = this.frictionCoeff;
    const dt = this.timestep;
    const k_B = 8.314e-3;
    const expTerm = Math.exp(-gamma * dt);
    const sqrtTerm = Math.sqrt(k_B * this.temperature * (1 - expTerm * expTerm));

    this.atoms.forEach((atom) => {
      if (this.isConstrained(atom)) return;
      const sigma = sqrtTerm / Math.sqrt(atom.mass);
      atom.velocity.x = atom.velocity.x * expTerm + this.randomGaussian(0, sigma);
      atom.velocity.y = atom.velocity.y * expTerm + this.randomGaussian(0, sigma);
      atom.velocity.z = atom.velocity.z * expTerm + this.randomGaussian(0, sigma);

      // Safety clamp so no atom can run away.
      const v = Math.sqrt(
        atom.velocity.x ** 2 + atom.velocity.y ** 2 + atom.velocity.z ** 2
      );
      if (v > this.maxVelocity) {
        const s = this.maxVelocity / v;
        atom.velocity.x *= s;
        atom.velocity.y *= s;
        atom.velocity.z *= s;
      }
    });
  }

  /**
   * Berendsen velocity rescaling: gently pulls the instantaneous temperature
   * toward the target each step. The per-step scale is clamped so the motion
   * stays smooth rather than jerky.
   */
  applyTemperatureControl() {
    const Tcur = this.calculateInstantaneousTemperature();
    if (!(Tcur > 0)) return;
    const tau = 0.02; // coupling time (ps); smaller = tighter tracking
    let lambda = Math.sqrt(1 + (this.timestep / tau) * (this.temperature / Tcur - 1));
    lambda = Math.max(0.7, Math.min(1.4, lambda));
    this.atoms.forEach((atom) => {
      if (this.isConstrained(atom)) return;
      atom.velocity.x *= lambda;
      atom.velocity.y *= lambda;
      atom.velocity.z *= lambda;
    });
  }

  isConstrained(atom) {
    if (this.constrainHeavyAtoms && atom.element === 'H') return false;
    return this.constrainedResidues.includes(atom.residueId);
  }

  updateEnergies() {
    const k_B = 8.314e-3;
    let KE = 0;
    this.atoms.forEach((atom) => {
      const v2 = atom.velocity.x ** 2 + atom.velocity.y ** 2 + atom.velocity.z ** 2;
      KE += 0.5 * atom.mass * v2;
    });

    const PE = this.forceCalc.calculatePotentialEnergy(this.atoms, this.bonds);

    this.energyStats.kinetic = KE;
    this.energyStats.potential = PE;
    this.energyStats.total = KE + PE;
    this.energyStats.history.push({
      step: this.stepCount,
      KE,
      PE,
      total: KE + PE,
      temp: this.calculateInstantaneousTemperature(),
    });

    if (this.onEnergyUpdate) this.onEnergyUpdate(this.energyStats);
  }

  calculateInstantaneousTemperature() {
    const k_B = 8.314e-3;
    const N_atoms = this.atoms.length || 1;
    const temp = (2 * this.energyStats.kinetic) / (3 * N_atoms * k_B);
    return Math.max(0, temp);
  }

  saveSnapshot() {
    const snapshot = {
      step: this.stepCount,
      time: this.stepCount * this.timestep,
      positions: this.atoms.map((a) => ({ ...a.position })),
      velocities: this.atoms.map((a) => ({ ...a.velocity })),
      energy: { ...this.energyStats },
      temperature: this.calculateInstantaneousTemperature(),
    };
    this.trajectory.push(snapshot);
  }

  start() {
    this.isRunning = true;
    this.isPaused = false;
    this.simulationLoop();
  }

  pause() {
    this.isPaused = !this.isPaused;
  }

  stop() {
    this.isRunning = false;
    this.isPaused = false;
  }

  simulationLoop = () => {
    if (!this.isRunning) return;
    const stepsPerFrame = 10;
    for (let i = 0; i < stepsPerFrame; i++) this.step();
    requestAnimationFrame(this.simulationLoop);
  };

  runSteps(nSteps) {
    for (let i = 0; i < nSteps; i++) this.step();
  }

  minimizeEnergy(maxSteps = 1000, tolerance = 1e-4) {
    console.log('Starting energy minimization...');
    const learningRate = 0.001;
    let prevEnergy = Infinity;

    for (let step = 0; step < maxSteps; step++) {
      this.atoms.forEach((atom) => {
        atom.force.x = 0;
        atom.force.y = 0;
        atom.force.z = 0;
      });

      const forces = this.forceCalc.calculateAllForces(this.atoms, this.bonds);
      this.atoms.forEach((atom, idx) => {
        if (this.isConstrained(atom)) return;
        const f = forces[idx];
        atom.position.x -= f.x * learningRate;
        atom.position.y -= f.y * learningRate;
        atom.position.z -= f.z * learningRate;
      });

      const PE = this.forceCalc.calculatePotentialEnergy(this.atoms, this.bonds);
      const energyDelta = Math.abs(PE - prevEnergy);

      if (step % 100 === 0) {
        console.log(`  Step ${step}: E = ${PE.toFixed(2)} kcal/mol, ΔE = ${energyDelta.toFixed(4)}`);
      }
      if (energyDelta < tolerance) {
        console.log(`Converged at step ${step}, E = ${PE.toFixed(2)} kcal/mol`);
        break;
      }
      prevEnergy = PE;
    }
  }

  getStats() {
    return {
      stepCount: this.stepCount,
      simTime: this.stepCount * this.timestep,
      temperature: this.calculateInstantaneousTemperature(),
      energy: this.energyStats.total,
      isRunning: this.isRunning,
      isPaused: this.isPaused,
      trajectoryLength: this.trajectory.length,
    };
  }

  exportTrajectory() {
    return {
      metadata: {
        timestep: this.timestep,
        nAtoms: this.atoms.length,
        nFrames: this.trajectory.length,
        atomLabels: this.atoms.map((a) => `${a.element}${a.id ?? a.serial ?? ''}`),
      },
      trajectory: this.trajectory,
    };
  }

  reset() {
    this.stepCount = 0;
    this.trajectory = [];
    this.energyStats = { kinetic: 0, potential: 0, total: 0, history: [] };
    this.initializeVelocities();
  }
}

export default PhysicsEngine;
