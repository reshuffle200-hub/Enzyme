/**
 * Inhibitor Builder Module
 *
 * Constructs and customizes small-molecule inhibitors from prebuilt scaffolds.
 */

export class InhibitorBuilder {
  constructor() {
    this.atoms = [];
    this.bonds = [];
    this.scaffolds = this.initializeScaffolds();
  }

  loadPreset(scaffoldName) {
    const scaffold = this.scaffolds[scaffoldName];
    if (!scaffold) {
      throw new Error(`Unknown scaffold: ${scaffoldName}`);
    }
    this.atoms = JSON.parse(JSON.stringify(scaffold.atoms));
    this.bonds = JSON.parse(JSON.stringify(scaffold.bonds));
    console.log(`Loaded scaffold: ${scaffoldName} (${this.atoms.length} atoms)`);
    return this;
  }

  initializeScaffolds() {
    return {
      'nucleoside-analog': {
        name: 'Nucleoside Analog',
        description: 'Mimics natural nucleotide, inhibits RNA polymerase',
        atoms: this.buildNucleosideAnalog(),
        bonds: this.buildNucleosideBonds(),
      },
      'protease-inhibitor': {
        name: 'Protease Inhibitor',
        description: 'Peptide-like inhibitor targeting active site',
        atoms: this.buildProteaseInhibitor(),
        bonds: this.buildProteaseInhibitorBonds(),
      },
      'allosteric-inhibitor': {
        name: 'Allosteric Inhibitor',
        description: 'Non-competitive inhibitor binding to regulatory site',
        atoms: this.buildAllosericInhibitor(),
        bonds: this.buildAllosericInhibitorBonds(),
      },
      'generic-small-molecule': {
        name: 'Generic Small Molecule',
        description: 'Simple drug-like molecule',
        atoms: this.buildGenericMolecule(),
        bonds: this.buildGenericMoleculeBonds(),
      },
    };
  }

  buildNucleosideAnalog() {
    return [
      { id: 1, element: 'C', name: 'C1', charge: 0, position: { x: 0, y: 0, z: 0 } },
      { id: 2, element: 'C', name: 'C2', charge: 0, position: { x: 1.5, y: 0, z: 0 } },
      { id: 3, element: 'C', name: 'C3', charge: 0, position: { x: 2.0, y: 1.3, z: 0 } },
      { id: 4, element: 'C', name: 'C4', charge: 0, position: { x: 0.8, y: 2.0, z: 0 } },
      { id: 5, element: 'O', name: 'O5', charge: -0.3, position: { x: -0.4, y: 1.2, z: 0 } },
      { id: 6, element: 'N', name: 'N1', charge: -0.2, position: { x: 2.2, y: -1.3, z: 0 } },
      { id: 7, element: 'C', name: 'C1N', charge: 0.1, position: { x: 3.5, y: -1.8, z: 0 } },
      { id: 8, element: 'N', name: 'N2', charge: -0.2, position: { x: 4.5, y: -1.0, z: 0 } },
      { id: 9, element: 'C', name: 'C2N', charge: 0.1, position: { x: 5.7, y: -1.4, z: 0 } },
      { id: 10, element: 'O', name: 'O2N', charge: -0.5, position: { x: 6.7, y: -0.7, z: 0 } },
      { id: 11, element: 'N', name: 'N3', charge: -0.2, position: { x: 5.5, y: -2.6, z: 0 } },
      { id: 12, element: 'P', name: 'P', charge: 0.8, position: { x: -1.8, y: 1.7, z: 0.8 } },
      { id: 13, element: 'O', name: 'OP1', charge: -0.6, position: { x: -2.8, y: 0.8, z: 1.2 } },
      { id: 14, element: 'O', name: 'OP2', charge: -0.6, position: { x: -2.2, y: 3.0, z: 1.2 } },
    ];
  }

  buildNucleosideBonds() {
    return [
      { atom1: 1, atom2: 2, type: 'single' },
      { atom1: 2, atom2: 3, type: 'single' },
      { atom1: 3, atom2: 4, type: 'single' },
      { atom1: 4, atom2: 5, type: 'single' },
      { atom1: 5, atom2: 1, type: 'single' },
      { atom1: 1, atom2: 6, type: 'single' },
      { atom1: 6, atom2: 7, type: 'double' },
      { atom1: 7, atom2: 8, type: 'single' },
      { atom1: 8, atom2: 9, type: 'single' },
      { atom1: 9, atom2: 10, type: 'double' },
      { atom1: 9, atom2: 11, type: 'single' },
      { atom1: 5, atom2: 12, type: 'single' },
      { atom1: 12, atom2: 13, type: 'double' },
      { atom1: 12, atom2: 14, type: 'double' },
    ];
  }

  buildProteaseInhibitor() {
    return [
      { id: 1, element: 'C', name: 'C1', charge: 0.1, position: { x: 0, y: 0, z: 0 } },
      { id: 2, element: 'O', name: 'O1', charge: -0.5, position: { x: -1.2, y: 0.3, z: 0 } },
      { id: 3, element: 'N', name: 'N1', charge: -0.3, position: { x: 0.2, y: -1.3, z: 0 } },
      { id: 4, element: 'C', name: 'CA1', charge: 0.1, position: { x: 1.5, y: 0.7, z: 0 } },
      { id: 5, element: 'C', name: 'CB1', charge: 0, position: { x: 2.5, y: 0.1, z: 1.0 } },
      { id: 6, element: 'C', name: 'CG1', charge: 0, position: { x: 3.9, y: 0.7, z: 1.2 } },
      { id: 7, element: 'C', name: 'CD1', charge: 0, position: { x: 4.9, y: 0.0, z: 2.0 } },
      { id: 8, element: 'C', name: 'C2', charge: 0.1, position: { x: 1.3, y: 2.2, z: 0 } },
      { id: 9, element: 'O', name: 'O2', charge: -0.5, position: { x: 0.2, y: 2.8, z: 0 } },
      { id: 10, element: 'N', name: 'N2', charge: -0.3, position: { x: 2.5, y: 2.8, z: 0 } },
      { id: 11, element: 'C', name: 'CB2', charge: 0, position: { x: 3.5, y: 2.3, z: 0.7 } },
      { id: 12, element: 'O', name: 'OG2', charge: -0.4, position: { x: 4.7, y: 2.8, z: 0.5 } },
      { id: 13, element: 'H', name: 'HG2', charge: 0.3, position: { x: 5.2, y: 2.4, z: 1.1 } },
      { id: 14, element: 'C', name: 'C3', charge: 0.3, position: { x: 2.4, y: 3.0, z: -1.2 } },
      { id: 15, element: 'O', name: 'OXT', charge: -0.6, position: { x: 3.6, y: 3.5, z: -1.2 } },
    ];
  }

  buildProteaseInhibitorBonds() {
    return [
      { atom1: 1, atom2: 2, type: 'double' },
      { atom1: 1, atom2: 3, type: 'single' },
      { atom1: 1, atom2: 4, type: 'single' },
      { atom1: 4, atom2: 5, type: 'single' },
      { atom1: 5, atom2: 6, type: 'single' },
      { atom1: 6, atom2: 7, type: 'single' },
      { atom1: 4, atom2: 8, type: 'single' },
      { atom1: 8, atom2: 9, type: 'double' },
      { atom1: 8, atom2: 10, type: 'single' },
      { atom1: 10, atom2: 11, type: 'single' },
      { atom1: 11, atom2: 12, type: 'single' },
      { atom1: 12, atom2: 13, type: 'single' },
      { atom1: 10, atom2: 14, type: 'single' },
      { atom1: 14, atom2: 15, type: 'double' },
    ];
  }

  buildAllosericInhibitor() {
    return [
      { id: 1, element: 'C', name: 'C1', charge: -0.1, position: { x: 0, y: 0, z: 0 } },
      { id: 2, element: 'C', name: 'C2', charge: -0.1, position: { x: 1.2, y: 0.7, z: 0 } },
      { id: 3, element: 'C', name: 'C3', charge: -0.1, position: { x: 2.4, y: 0, z: 0 } },
      { id: 4, element: 'C', name: 'C4', charge: -0.1, position: { x: 2.4, y: -1.4, z: 0 } },
      { id: 5, element: 'C', name: 'C5', charge: -0.1, position: { x: 1.2, y: -2.1, z: 0 } },
      { id: 6, element: 'C', name: 'C6', charge: -0.1, position: { x: 0, y: -1.4, z: 0 } },
      { id: 7, element: 'C', name: 'LINK', charge: 0.1, position: { x: 3.7, y: 0.7, z: 0 } },
      { id: 8, element: 'N', name: 'NL', charge: -0.2, position: { x: 4.8, y: 0.1, z: 0 } },
      { id: 9, element: 'C', name: 'C7', charge: -0.1, position: { x: 6.0, y: 0.8, z: 0 } },
      { id: 10, element: 'C', name: 'C8', charge: -0.1, position: { x: 7.2, y: 0.1, z: 0 } },
      { id: 11, element: 'C', name: 'C9', charge: -0.1, position: { x: 8.4, y: 0.8, z: 0 } },
      { id: 12, element: 'C', name: 'C10', charge: -0.1, position: { x: 8.4, y: 2.2, z: 0 } },
      { id: 13, element: 'C', name: 'C11', charge: -0.1, position: { x: 7.2, y: 2.9, z: 0 } },
      { id: 14, element: 'C', name: 'C12', charge: -0.1, position: { x: 6.0, y: 2.2, z: 0 } },
      { id: 15, element: 'O', name: 'O1', charge: -0.4, position: { x: 9.6, y: 0.1, z: 0 } },
      { id: 16, element: 'H', name: 'HO1', charge: 0.3, position: { x: 10.3, y: 0.6, z: 0 } },
    ];
  }

  buildAllosericInhibitorBonds() {
    return [
      { atom1: 1, atom2: 2, type: 'aromatic' },
      { atom1: 2, atom2: 3, type: 'aromatic' },
      { atom1: 3, atom2: 4, type: 'aromatic' },
      { atom1: 4, atom2: 5, type: 'aromatic' },
      { atom1: 5, atom2: 6, type: 'aromatic' },
      { atom1: 6, atom2: 1, type: 'aromatic' },
      { atom1: 3, atom2: 7, type: 'single' },
      { atom1: 7, atom2: 8, type: 'double' },
      { atom1: 8, atom2: 9, type: 'single' },
      { atom1: 9, atom2: 10, type: 'aromatic' },
      { atom1: 10, atom2: 11, type: 'aromatic' },
      { atom1: 11, atom2: 12, type: 'aromatic' },
      { atom1: 12, atom2: 13, type: 'aromatic' },
      { atom1: 13, atom2: 14, type: 'aromatic' },
      { atom1: 14, atom2: 9, type: 'aromatic' },
      { atom1: 11, atom2: 15, type: 'single' },
      { atom1: 15, atom2: 16, type: 'single' },
    ];
  }

  buildGenericMolecule() {
    return [
      { id: 1, element: 'C', name: 'C1', charge: 0, position: { x: 0, y: 0, z: 0 } },
      { id: 2, element: 'C', name: 'C2', charge: 0, position: { x: 1.5, y: 0, z: 0 } },
      { id: 3, element: 'O', name: 'O1', charge: -0.4, position: { x: 2.3, y: -0.8, z: 0 } },
      { id: 4, element: 'N', name: 'N1', charge: -0.2, position: { x: 2.2, y: 1.2, z: 0 } },
      { id: 5, element: 'C', name: 'C3', charge: 0, position: { x: 3.6, y: 1.2, z: 0 } },
      { id: 6, element: 'C', name: 'C4', charge: 0, position: { x: 4.4, y: 0, z: 0 } },
      { id: 7, element: 'H', name: 'H1', charge: 0.2, position: { x: 5.8, y: 0, z: 0 } },
    ];
  }

  buildGenericMoleculeBonds() {
    return [
      { atom1: 1, atom2: 2, type: 'single' },
      { atom1: 2, atom2: 3, type: 'double' },
      { atom1: 2, atom2: 4, type: 'single' },
      { atom1: 4, atom2: 5, type: 'single' },
      { atom1: 5, atom2: 6, type: 'single' },
      { atom1: 6, atom2: 7, type: 'single' },
    ];
  }

  addAtom(element, position, charge = 0) {
    const id = Math.max(0, ...this.atoms.map((a) => a.id)) + 1;
    this.atoms.push({ id, element, name: `${element}${id}`, charge, position });
    return this;
  }

  removeAtom(atomId) {
    this.atoms = this.atoms.filter((a) => a.id !== atomId);
    this.bonds = this.bonds.filter((b) => b.atom1 !== atomId && b.atom2 !== atomId);
    return this;
  }

  addBond(atom1Id, atom2Id, type = 'single') {
    this.bonds.push({ atom1: atom1Id, atom2: atom2Id, type });
    return this;
  }

  setCharge(atomId, charge) {
    const atom = this.atoms.find((a) => a.id === atomId);
    if (atom) atom.charge = charge;
    return this;
  }

  moveAtom(atomId, dx, dy, dz) {
    const atom = this.atoms.find((a) => a.id === atomId);
    if (atom) {
      atom.position.x += dx;
      atom.position.y += dy;
      atom.position.z += dz;
    }
    return this;
  }

  translate(dx, dy, dz) {
    this.atoms.forEach((atom) => {
      atom.position.x += dx;
      atom.position.y += dy;
      atom.position.z += dz;
    });
    return this;
  }

  rotate(angle) {
    const cos = Math.cos(angle);
    const sin = Math.sin(angle);
    this.atoms.forEach((atom) => {
      const x = atom.position.x;
      const y = atom.position.y;
      atom.position.x = x * cos - y * sin;
      atom.position.y = x * sin + y * cos;
    });
    return this;
  }

  build() {
    return {
      atoms: this.atoms,
      bonds: this.bonds,
      properties: {
        atomCount: this.atoms.length,
        bondCount: this.bonds.length,
        totalCharge: this.atoms.reduce((sum, a) => sum + a.charge, 0),
        centerOfMass: this.calculateCenterOfMass(),
      },
    };
  }

  calculateCenterOfMass() {
    const com = { x: 0, y: 0, z: 0 };
    let totalMass = 0;
    const masses = { H: 1.008, C: 12.011, N: 14.007, O: 15.999, S: 32.065, P: 30.974 };
    this.atoms.forEach((atom) => {
      const mass = masses[atom.element] || 12;
      com.x += atom.position.x * mass;
      com.y += atom.position.y * mass;
      com.z += atom.position.z * mass;
      totalMass += mass;
    });
    if (totalMass > 0) {
      com.x /= totalMass;
      com.y /= totalMass;
      com.z /= totalMass;
    }
    return com;
  }

  exportJSON() {
    return JSON.stringify(this.build(), null, 2);
  }
}

export default InhibitorBuilder;
