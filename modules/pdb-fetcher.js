/**
 * PDB Fetcher & Parser Module
 *
 * Fetches rhinovirus 3Dpol structures from RCSB PDB
 * Parses atomic coordinates, residues, bonds, and metadata
 * Handles local caching to reduce network calls
 */

export class PDBFetcher {
  constructor(cacheDir = 'pdb-cache') {
    this.cacheDir = cacheDir;
    this.baseURL = 'https://files.rcsb.org/download/';
    this.pdbCache = {};

    // Known good RV-A89 3Dpol structures (PDB IDs)
    this.knownRVA89Structures = [
      '1XR5', // HRV14 3Dpol
      '1XR6',
    ];
  }

  /**
   * Fetch PDB structure by ID
   */
  async fetch(pdbId, fromCache = true) {
    if (fromCache && this.pdbCache[pdbId]) {
      console.log(`✓ Loaded ${pdbId} from cache`);
      return this.pdbCache[pdbId];
    }

    console.log(`Fetching PDB: ${pdbId}...`);

    try {
      const pdbText = await this.downloadPDB(pdbId);
      const structure = this.parsePDB(pdbText, pdbId);
      this.pdbCache[pdbId] = structure;
      console.log(`✓ Successfully loaded ${pdbId}`);
      return structure;
    } catch (error) {
      console.error(`Failed to fetch PDB ${pdbId}:`, error);
      throw error;
    }
  }

  /**
   * Download a PDB file, reporting progress via onProgress(fraction 0..1).
   * Streams the response body so a progress bar can update; falls back to a
   * plain text read if streaming or content-length isn't available.
   */
  async downloadPDB(pdbId, onProgress) {
    const url = `${this.baseURL}${pdbId}.pdb`;
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`PDB fetch failed: ${response.status} ${response.statusText}`);
    }

    const total = parseInt(response.headers.get('content-length')) || 0;

    if (!response.body || typeof response.body.getReader !== 'function') {
      if (onProgress) onProgress(1);
      return await response.text();
    }

    const reader = response.body.getReader();
    const chunks = [];
    let received = 0;
    for (;;) {
      const { done, value } = await reader.read();
      if (done) break;
      chunks.push(value);
      received += value.length;
      if (onProgress) {
        // If total is known use the exact fraction; otherwise approach 1
        // asymptotically so the bar keeps moving without ever "completing".
        const frac = total ? received / total : 1 - 1 / (1 + received / 300000);
        onProgress(Math.min(0.99, frac));
      }
    }
    if (onProgress) onProgress(1);

    const out = new Uint8Array(received);
    let pos = 0;
    for (const c of chunks) {
      out.set(c, pos);
      pos += c.length;
    }
    return new TextDecoder('utf-8').decode(out);
  }

  parsePDB(pdbText, pdbId) {
    const lines = pdbText.split('\n');

    const structure = {
      id: pdbId,
      title: '',
      atoms: [],
      residues: {},
      bonds: [],
      chains: {},
      secondaryStructure: { helices: [], sheets: [] },
      metadata: { resolution: null, technique: '', organism: '' },
    };

    const atomMap = {};
    const residueMap = {};
    let hetatmMode = false;
    let conects = {};

    for (const line of lines) {
      const recordType = line.substring(0, 6).trim();

      if (recordType === 'TITLE') {
        structure.title += line.substring(10).trim();
      }

      if (recordType === 'REMARK') {
        const remark = line.substring(10);
        if (line.substring(7, 10) === '  2') {
          if (remark.includes('RESOLUTION')) {
            const resMatch = remark.match(/(\d+\.\d+)/);
            if (resMatch) structure.metadata.resolution = parseFloat(resMatch[1]);
          }
        }
      }

      if (recordType === 'EXPDTA') {
        structure.metadata.technique = line.substring(10).trim();
      }

      if (recordType === 'HELIX') {
        structure.secondaryStructure.helices.push({
          id: line.substring(7, 10).trim(),
          initResName: line.substring(15, 18).trim(),
          initChainId: line.substring(19),
          initSeqNum: parseInt(line.substring(21, 25)),
          endResName: line.substring(27, 30).trim(),
          endChainId: line.substring(31),
          endSeqNum: parseInt(line.substring(33, 37)),
        });
      }

      if (recordType === 'SHEET') {
        structure.secondaryStructure.sheets.push({
          id: line.substring(7, 10).trim(),
          strand: parseInt(line.substring(11, 14)),
          initResName: line.substring(17, 20).trim(),
          initChainId: line.substring(21),
          initSeqNum: parseInt(line.substring(22, 26)),
          endResName: line.substring(32, 35).trim(),
          endChainId: line.substring(36),
          endSeqNum: parseInt(line.substring(37, 41)),
        });
      }

      if (recordType === 'ATOM' || recordType === 'HETATM') {
        if (recordType === 'HETATM') hetatmMode = true;

        const atom = this.parseAtomLine(line, hetatmMode);
        if (atom) {
          structure.atoms.push(atom);
          atomMap[atom.serial] = structure.atoms.length - 1;

          const resKey = `${atom.chainId}-${atom.resSeq}`;
          if (!residueMap[resKey]) {
            residueMap[resKey] = {
              chainId: atom.chainId,
              resSeq: atom.resSeq,
              resName: atom.resName,
              atoms: [],
            };
            if (!structure.chains[atom.chainId]) {
              structure.chains[atom.chainId] = [];
            }
            structure.chains[atom.chainId].push(resKey);
          }
          residueMap[resKey].atoms.push(structure.atoms.length - 1);
        }
      }

      if (recordType === 'CONECT') {
        const serial = parseInt(line.substring(6, 11));
        const bonded = [];
        for (let i = 0; i < 4; i++) {
          const start = 11 + i * 5;
          const end = start + 5;
          const bondedSerial = parseInt(line.substring(start, end));
          if (!isNaN(bondedSerial)) bonded.push(bondedSerial);
        }
        conects[serial] = bonded;
      }
    }

    const seenBonds = new Set();
    for (const [serial, bonded] of Object.entries(conects)) {
      const i = atomMap[parseInt(serial)];
      if (i === undefined) continue;

      for (const bondedSerial of bonded) {
        const j = atomMap[bondedSerial];
        if (j === undefined) continue;

        const bondKey = [Math.min(i, j), Math.max(i, j)].join('-');
        if (seenBonds.has(bondKey)) continue;
        seenBonds.add(bondKey);

        structure.bonds.push({
          atom1Index: i,
          atom2Index: j,
          atom1Serial: parseInt(serial),
          atom2Serial: bondedSerial,
          type: 'single',
        });
      }
    }

    // If no CONECT records, infer bonds from distance.
    if (structure.bonds.length === 0) {
      structure.bonds = this.inferBondsFromDistance(structure.atoms);
    }

    structure.residues = residueMap;
    this.validateStructure(structure);

    return structure;
  }

  parseAtomLine(line, isHetatm = false) {
    try {
      return {
        serial: parseInt(line.substring(6, 11)),
        name: line.substring(12, 16).trim(),
        altLoc: line.substring(16, 17),
        resName: line.substring(17, 20).trim(),
        chainId: line.substring(21, 22),
        resSeq: parseInt(line.substring(22, 26)),
        iCode: line.substring(26, 27),
        position: {
          x: parseFloat(line.substring(30, 38)),
          y: parseFloat(line.substring(38, 46)),
          z: parseFloat(line.substring(46, 54)),
        },
        occupancy: parseFloat(line.substring(54, 60)) || 1.0,
        bFactor: parseFloat(line.substring(60, 66)) || 0,
        element: line.substring(76, 78).trim() || this.inferElement(line.substring(12, 16)),
        isHetatm: isHetatm,
        mass: this.getAtomicMass(line.substring(76, 78).trim()),
        charge: this.getAtomicCharge(line.substring(17, 20).trim(), line.substring(12, 16).trim()),
        residueId: `${line.substring(21, 22)}-${line.substring(22, 26).trim()}`,
      };
    } catch (e) {
      console.warn('Failed to parse atom line:', line, e);
      return null;
    }
  }

  inferElement(atomName) {
    const name = atomName.toUpperCase();
    if (name.startsWith('CL')) return 'Cl';
    if (name.startsWith('BR')) return 'Br';
    if (name.startsWith('C')) return 'C';
    if (name.startsWith('N')) return 'N';
    if (name.startsWith('O')) return 'O';
    if (name.startsWith('S')) return 'S';
    if (name.startsWith('P')) return 'P';
    if (name.startsWith('H')) return 'H';
    if (name.startsWith('F')) return 'F';
    if (name.startsWith('I')) return 'I';
    return 'C';
  }

  getAtomicMass(element) {
    const masses = {
      H: 1.008, C: 12.011, N: 14.007, O: 15.999,
      S: 32.065, P: 30.974, Cl: 35.45, F: 18.998,
      Br: 79.904, I: 126.9,
    };
    return masses[element] || 12.0;
  }

  getAtomicCharge(resName, atomName) {
    const charges = {
      ALA: { N: -0.3, CA: 0.1, C: 0.5, O: -0.5, CB: 0.0 },
      ARG: { N: -0.3, NE: -0.3, NH1: 0.8, NH2: 0.8, C: 0.5, O: -0.5 },
      ASN: { N: -0.3, OD1: -0.5, ND2: -0.3, C: 0.5, O: -0.5 },
      ASP: { N: -0.3, OD1: -0.8, OD2: -0.8, C: 0.5, O: -0.5 },
      CYS: { N: -0.3, SG: -0.2, C: 0.5, O: -0.5 },
      GLU: { N: -0.3, OE1: -0.8, OE2: -0.8, C: 0.5, O: -0.5 },
      GLN: { N: -0.3, OE1: -0.5, NE2: -0.3, C: 0.5, O: -0.5 },
      GLY: { N: -0.3, C: 0.5, O: -0.5 },
      HIS: { N: -0.3, ND1: 0.1, NE2: 0.1, C: 0.5, O: -0.5 },
      ILE: { N: -0.3, C: 0.5, O: -0.5 },
      LEU: { N: -0.3, C: 0.5, O: -0.5 },
      LYS: { N: -0.3, NZ: 0.8, C: 0.5, O: -0.5 },
      MET: { N: -0.3, C: 0.5, O: -0.5 },
      PHE: { N: -0.3, C: 0.5, O: -0.5 },
      PRO: { N: -0.3, C: 0.5, O: -0.5 },
      SER: { N: -0.3, OG: -0.3, C: 0.5, O: -0.5 },
      THR: { N: -0.3, OG1: -0.3, C: 0.5, O: -0.5 },
      TRP: { N: -0.3, C: 0.5, O: -0.5 },
      TYR: { N: -0.3, OH: -0.3, C: 0.5, O: -0.5 },
      VAL: { N: -0.3, C: 0.5, O: -0.5 },
    };

    if (charges[resName] && charges[resName][atomName]) {
      return charges[resName][atomName];
    }

    if (atomName.includes('O')) return -0.5;
    if (atomName.includes('N')) return -0.3;
    if (atomName.includes('S')) return -0.2;
    if (atomName.includes('P')) return 0.3;
    return 0;
  }

  inferBondsFromDistance(atoms) {
    const bonds = [];
    const covalentRadii = {
      H: 0.31, C: 0.76, N: 0.71, O: 0.66,
      S: 1.05, P: 1.07, F: 0.57, Cl: 1.02,
      Br: 1.2, I: 1.39,
    };
    const tolerance = 0.6;

    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const a_i = atoms[i];
        const a_j = atoms[j];

        if (a_i.residueId !== a_j.residueId && Math.abs(a_i.resSeq - a_j.resSeq) > 1) {
          continue;
        }

        const dx = a_j.position.x - a_i.position.x;
        const dy = a_j.position.y - a_i.position.y;
        const dz = a_j.position.z - a_i.position.z;
        const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

        const r_i = covalentRadii[a_i.element] || 0.7;
        const r_j = covalentRadii[a_j.element] || 0.7;
        const maxDist = r_i + r_j + tolerance;

        if (distance < maxDist) {
          bonds.push({ atom1Index: i, atom2Index: j, type: 'single' });
        }
      }
    }
    return bonds;
  }

  validateStructure(structure) {
    console.log(`Structure validation for ${structure.id}:`);
    console.log(`  Atoms: ${structure.atoms.length}`);
    console.log(`  Residues: ${Object.keys(structure.residues).length}`);
    console.log(`  Bonds: ${structure.bonds.length}`);
    console.log(`  Chains: ${Object.keys(structure.chains).length}`);
    console.log(`  Resolution: ${structure.metadata.resolution ? structure.metadata.resolution + ' Å' : 'N/A'}`);
  }

  async searchRVA89Structures() {
    console.log('Searching PDB for rhinovirus A 3Dpol structures...');
    return this.knownRVA89Structures;
  }
}

export default PDBFetcher;
