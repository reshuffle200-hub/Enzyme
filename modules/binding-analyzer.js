/**
 * Binding Analyzer Module
 *
 * Analyzes protein–ligand interactions: hydrogen bonds, electrostatics,
 * hydrophobic contacts, shape complementarity, steric clashes, and an
 * empirical binding score / ΔG / Kd estimate.
 *
 * FIX: hydrogen-bond and electrostatic entries now carry the actual atom
 * positions (atom1/atom2 objects) so the visualization engine can draw the
 * interaction lines. Previously they stored only atom *names*, which couldn't
 * be mapped back to 3D coordinates.
 */

export class BindingAnalyzer {
  constructor(structure, inhibitor, options = {}) {
    this.structure = structure;
    this.inhibitor = inhibitor;
    this.trajectory = null;

    this.hbondParams = options.hbondParams || {
      maxDistance: 3.5, minAngle: 120, optDistance: 2.8, weight: 5.0,
    };
    this.electrostaticParams = options.electrostaticParams || { weight: 1.0, dielectric: 4.0 };
    this.hydrophobicParams = options.hydrophobicParams || { weight: 0.5, contactDistance: 4.0 };
    this.shapeParams = options.shapeParams || { weight: 2.0, overlapPenalty: 1.0 };

    this.interactionHistory = { hbonds: [], electrostatic: [], hydrophobic: [], steric: [] };
  }

  analyzeBinding() {
    console.log('Starting binding analysis...');

    const metrics = {
      timestamp: new Date().toISOString(),
      metrics: {
        hydrogenBonds: this.analyzeHydrogenBonds(),
        electrostaticInteractions: this.analyzeElectrostatic(),
        hydrophobicInteractions: this.analyzeHydrophobic(),
        shapeComplementarity: this.analyzeShapeComplementarity(),
        stericClashes: this.analyzeStericClashes(),
        interfaceArea: this.calculateInterfaceArea(),
        contactResidues: this.identifyContactResidues(),
      },
    };

    metrics.bindingScore = this.calculateBindingScore(metrics.metrics);
    metrics.estimatedKd = this.estimateBindingAffinity(metrics.bindingScore);
    metrics.estimatedDeltaG = this.calculateDeltaG(metrics.metrics);

    if (this.trajectory) metrics.stability = this.analyzeStability(this.trajectory);

    metrics.prediction = this.makePrediction(metrics);
    return metrics;
  }

  analyzeHydrogenBonds() {
    const hbonds = [];
    const inhibitorAtoms = this.inhibitor.atoms;
    const proteinAtoms = this.structure.atoms;

    const donors = ['N', 'O'];
    const acceptors = ['O', 'N'];

    for (const lig of inhibitorAtoms) {
      for (const prot of proteinAtoms) {
        let isDonor1 = donors.includes(lig.element);
        let isAcceptor2 = acceptors.includes(prot.element);

        const dx = prot.position.x - lig.position.x;
        const dy = prot.position.y - lig.position.y;
        const dz = prot.position.z - lig.position.z;
        const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

        if (distance > this.hbondParams.maxDistance) continue;

        if ((isDonor1 && isAcceptor2) || (!isDonor1 && !isAcceptor2)) {
          if (distance < this.hbondParams.optDistance * 1.2) {
            hbonds.push({
              ligAtom: lig.name,
              protAtom: prot.name,
              protResidue: prot.resName,
              distance: distance.toFixed(2),
              energy: this.estimateHbondEnergy(distance),
              geometry: 'putative',
              // Positions so the viewer can draw the dashed line:
              atom1: { position: { x: lig.position.x, y: lig.position.y, z: lig.position.z } },
              atom2: { position: { x: prot.position.x, y: prot.position.y, z: prot.position.z } },
            });
          }
        }
      }
    }

    return {
      count: hbonds.length,
      bonds: hbonds,
      energy: hbonds.reduce((sum, h) => sum + h.energy, 0),
    };
  }

  estimateHbondEnergy(distance) {
    const r_opt = this.hbondParams.optDistance;
    const sigma = 0.3;
    const maxEnergy = this.hbondParams.weight;
    const energy = -maxEnergy * Math.exp(-Math.pow(distance - r_opt, 2) / (2 * sigma * sigma));
    return parseFloat(energy.toFixed(2));
  }

  analyzeElectrostatic() {
    let totalEnergy = 0;
    const interactions = [];
    const inhibitorAtoms = this.inhibitor.atoms;
    const proteinAtoms = this.structure.atoms;

    const k = 332.06;
    const epsilon = this.electrostaticParams.dielectric;
    const cutoff = 15.0;

    for (const lig of inhibitorAtoms) {
      if (!lig.charge) continue;
      for (const prot of proteinAtoms) {
        if (!prot.charge) continue;

        const dx = prot.position.x - lig.position.x;
        const dy = prot.position.y - lig.position.y;
        const dz = prot.position.z - lig.position.z;
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);

        if (r > cutoff || r < 0.1) continue;

        const energy = (k * lig.charge * prot.charge) / (epsilon * r);

        if (Math.abs(energy) > 0.1) {
          interactions.push({
            ligAtom: lig.name,
            ligCharge: lig.charge.toFixed(2),
            protAtom: prot.name,
            protCharge: prot.charge.toFixed(2),
            protResidue: prot.resName,
            distance: r.toFixed(2),
            energy: energy.toFixed(2),
            type: lig.charge * prot.charge > 0 ? 'repulsion' : 'attraction',
            atom1: { position: { x: lig.position.x, y: lig.position.y, z: lig.position.z } },
            atom2: { position: { x: prot.position.x, y: prot.position.y, z: prot.position.z } },
          });
          totalEnergy += energy;
        }
      }
    }

    interactions.sort((a, b) => Math.abs(parseFloat(b.energy)) - Math.abs(parseFloat(a.energy)));

    return {
      energy: parseFloat(totalEnergy.toFixed(2)),
      count: interactions.length,
      interactions: interactions.slice(0, 10),
    };
  }

  analyzeHydrophobic() {
    const hydrophobicElements = ['C', 'S'];
    const inhibitorAtoms = this.inhibitor.atoms;
    const proteinAtoms = this.structure.atoms;

    let contactCount = 0;
    let clusterScore = 0;
    const contacts = [];

    for (const lig of inhibitorAtoms) {
      if (!hydrophobicElements.includes(lig.element)) continue;
      for (const prot of proteinAtoms) {
        if (!hydrophobicElements.includes(prot.element)) continue;
        if (prot.element === 'C' && prot.resName === 'GLY') continue;

        const dx = prot.position.x - lig.position.x;
        const dy = prot.position.y - lig.position.y;
        const dz = prot.position.z - lig.position.z;
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);

        if (r < this.hydrophobicParams.contactDistance) {
          contactCount++;
          clusterScore += Math.exp(-r / 2.0);
          if (contacts.length < 20) {
            contacts.push({
              atom1: { position: { x: lig.position.x, y: lig.position.y, z: lig.position.z } },
              atom2: { position: { x: prot.position.x, y: prot.position.y, z: prot.position.z } },
            });
          }
        }
      }
    }

    return {
      contactCount,
      contacts,
      clusterScore: parseFloat(clusterScore.toFixed(2)),
      strength: contactCount > 5 ? 'strong' : contactCount > 2 ? 'moderate' : 'weak',
    };
  }

  analyzeShapeComplementarity() {
    const inhibitorAtoms = this.inhibitor.atoms;
    const proteinAtoms = this.structure.atoms;
    const vdwRadii = { H: 1.2, C: 1.7, N: 1.55, O: 1.52, S: 1.8, P: 1.8, F: 1.47, Cl: 1.77 };

    let overlapVolume = 0;
    let complementaryContact = 0;

    for (const lig of inhibitorAtoms) {
      const r_lig = vdwRadii[lig.element] || 1.7;
      for (const prot of proteinAtoms) {
        const r_prot = vdwRadii[prot.element] || 1.7;

        const dx = prot.position.x - lig.position.x;
        const dy = prot.position.y - lig.position.y;
        const dz = prot.position.z - lig.position.z;
        const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

        const contactDistance = r_lig + r_prot;
        if (distance < contactDistance * 0.8) {
          overlapVolume += Math.pow(contactDistance - distance, 2);
        } else if (distance < contactDistance * 1.1) {
          complementaryContact++;
        }
      }
    }

    const complementarity =
      Math.max(0, complementaryContact - overlapVolume * 10) / Math.max(1, complementaryContact);

    return {
      complementarity: parseFloat(Math.min(1.0, complementarity).toFixed(3)),
      contactCount: complementaryContact,
      overlapPenalty: parseFloat(overlapVolume.toFixed(2)),
    };
  }

  analyzeStericClashes() {
    const vdwRadii = { H: 1.2, C: 1.7, N: 1.55, O: 1.52, S: 1.8, P: 1.8, F: 1.47, Cl: 1.77 };
    const inhibitorAtoms = this.inhibitor.atoms;
    const proteinAtoms = this.structure.atoms;

    let clashCount = 0;
    let clashPenalty = 0;

    for (const lig of inhibitorAtoms) {
      const r_lig = vdwRadii[lig.element] || 1.7;
      for (const prot of proteinAtoms) {
        const r_prot = vdwRadii[prot.element] || 1.7;

        const dx = prot.position.x - lig.position.x;
        const dy = prot.position.y - lig.position.y;
        const dz = prot.position.z - lig.position.z;
        const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

        const minDistance = r_lig + r_prot;
        const clashThreshold = minDistance * 0.8;

        if (distance < clashThreshold) {
          clashCount++;
          clashPenalty += Math.pow(minDistance - distance, 2);
        }
      }
    }

    return {
      count: clashCount,
      penalty: parseFloat(clashPenalty.toFixed(2)),
      severity:
        clashCount > 5 ? 'severe' : clashCount > 2 ? 'moderate' : clashCount > 0 ? 'minor' : 'none',
    };
  }

  calculateInterfaceArea() {
    const contactDistance = 4.5;
    const contactAtoms = new Set();
    for (const lig of this.inhibitor.atoms) {
      for (const prot of this.structure.atoms) {
        const dx = prot.position.x - lig.position.x;
        const dy = prot.position.y - lig.position.y;
        const dz = prot.position.z - lig.position.z;
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (r < contactDistance) contactAtoms.add(prot.serial);
      }
    }
    return parseFloat((contactAtoms.size * 4).toFixed(1));
  }

  identifyContactResidues() {
    const contactDistance = 4.0;
    const residueSet = new Set();
    for (const lig of this.inhibitor.atoms) {
      for (const prot of this.structure.atoms) {
        const dx = prot.position.x - lig.position.x;
        const dy = prot.position.y - lig.position.y;
        const dz = prot.position.z - lig.position.z;
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (r < contactDistance) residueSet.add(`${prot.resName}${prot.resSeq}`);
      }
    }
    return Array.from(residueSet);
  }

  calculateBindingScore(metrics) {
    let score = 0;
    score += Math.min(7.5, metrics.hydrogenBonds.count * 1.5);
    score += metrics.shapeComplementarity.complementarity * 2;

    const esEnergy = metrics.electrostaticInteractions.energy;
    if (esEnergy < -10) score += 0.5;
    else if (esEnergy > 10) score -= 1;

    score += Math.min(1.5, metrics.hydrophobicInteractions.clusterScore * 0.3);
    score -= Math.min(2, metrics.stericClashes.penalty / 10);

    return parseFloat(Math.max(0, Math.min(10, score)).toFixed(2));
  }

  estimateBindingAffinity(bindingScore) {
    const R = 1.987e-3;
    const T = 298;
    const deltaG = -bindingScore * 1.2;
    const Kd = Math.exp(-deltaG / (R * T));
    return Kd; // dissociation constant (M)
  }

  calculateDeltaG(metrics) {
    let deltaG = 0;
    deltaG -= metrics.hydrogenBonds.count * 5;
    deltaG += metrics.electrostaticInteractions.energy * 0.5;
    deltaG -= metrics.hydrophobicInteractions.contactCount * 0.5;
    deltaG -= metrics.shapeComplementarity.complementarity * 2;
    deltaG += metrics.stericClashes.penalty;
    deltaG += 2.0;
    return parseFloat(deltaG.toFixed(2));
  }

  analyzeStability(trajectory) {
    if (!trajectory || trajectory.length < 2) return { status: 'insufficient_data' };
    let stableFrames = 0;
    for (const frame of trajectory) {
      let avgDistance = 0;
      let count = 1;
      avgDistance /= count;
      if (avgDistance < 5.0) stableFrames++;
    }
    const stability = stableFrames / trajectory.length;
    return {
      stability: parseFloat((stability * 100).toFixed(1)),
      status: stability > 0.8 ? 'stable' : stability > 0.5 ? 'moderate' : 'unstable',
    };
  }

  makePrediction(metrics) {
    const score = metrics.bindingScore;
    const Kd = metrics.estimatedKd;

    let prediction = '';
    if (score > 7) prediction = 'Excellent binding - Strong inhibitor candidate';
    else if (score > 5) prediction = 'Good binding - Moderate inhibitor potential';
    else if (score > 3) prediction = 'Weak binding - Optimization needed';
    else prediction = 'Very weak or no binding - Requires redesign';

    if (Kd < 1e-9) prediction += ' (Kd: sub-nanomolar)';
    else if (Kd < 1e-6) prediction += ' (Kd: nanomolar)';
    else if (Kd < 1e-3) prediction += ' (Kd: micromolar)';
    else prediction += ' (Kd: weak)';

    return prediction;
  }

  exportJSON() {
    return JSON.stringify(this.analyzeBinding(), null, 2);
  }

  exportCSV() {
    const analysis = this.analyzeBinding();
    const m = analysis.metrics;

    let csv = 'Metric,Value,Unit\n';
    csv += `Hydrogen Bonds,${m.hydrogenBonds.count},count\n`;
    csv += `H-Bond Energy,${m.hydrogenBonds.energy},kcal/mol\n`;
    csv += `Electrostatic Energy,${m.electrostaticInteractions.energy},kcal/mol\n`;
    csv += `Shape Complementarity,${m.shapeComplementarity.complementarity},score (0-1)\n`;
    csv += `Hydrophobic Contacts,${m.hydrophobicInteractions.contactCount},count\n`;
    csv += `Steric Clashes,${m.stericClashes.count},count\n`;
    csv += `Interface Area,${m.interfaceArea},A^2\n`;
    csv += `Binding Score,${analysis.bindingScore},/10\n`;
    csv += `Estimated DeltaG,${analysis.estimatedDeltaG},kcal/mol\n`;
    csv += `Estimated Kd,${analysis.estimatedKd},M\n`;

    return csv;
  }
}

export default BindingAnalyzer;
