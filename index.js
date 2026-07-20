/**
 * Rhinovirus 3Dpol MD Simulation - Main Application Controller
 *
 * Connects all modules together:
 * - PDB fetching
 * - Physics engine
 * - Binding analysis
 * - Visualization
 *
 * CHANGES vs. the original controller:
 * - The 3D view now plays back the simulation as a recorded trajectory:
 *   each completed physics step is streamed into the visualization engine's
 *   trajectory up front, then plays it back smoothly (see captureTrajectory
 *   + Start/Pause/Reset handlers).
 * - Hydrophobic contacts are drawn in addition to H-bonds/electrostatics.
 * - All visualization calls use the method names the engine now exposes.
 */

import { PDBFetcher } from './modules/pdb-fetcher.js';
import { VisualizationEngine } from './modules/visualization-engine.js';
import { PhysicsEngine } from './modules/physics-engine.js';
import { InhibitorBuilder } from './modules/inhibitor-builder.js';
import { BindingAnalyzer } from './modules/binding-analyzer.js';

// Application state
let appState = {
  structure: null,
  inhibitor: null,
  physics: null,
  visualization: null,
  bindingAnalyzer: null,
  pdbFetcher: new PDBFetcher(),
  isRunning: false,
  trajCount: 0, // how many recorded frames have been streamed to the viewer
};

// DOM elements
const canvas = document.getElementById('canvas');
const pdbIdInput = document.getElementById('pdbId');
const loadBtn = document.getElementById('loadBtn');
const addInhibitorBtn = document.getElementById('addInhibitorBtn');
const inhibitorType = document.getElementById('inhibitorType');
const startBtn = document.getElementById('startBtn');
const pauseBtn = document.getElementById('pauseBtn');
const resetBtn = document.getElementById('resetBtn');
const viewMode = document.getElementById('viewMode');
const resetCameraBtn = document.getElementById('resetCameraBtn');
const zoomFitBtn = document.getElementById('zoomFitBtn');
const analyzeBtn = document.getElementById('analyzeBtn');
const exportJsonBtn = document.getElementById('exportJsonBtn');
const exportCsvBtn = document.getElementById('exportCsvBtn');
const exportPdbBtn = document.getElementById('exportPdbBtn');

const simTime = document.getElementById('simTime');
const stepCount = document.getElementById('stepCount');
const tempValue = document.getElementById('tempValue');
const energyValue = document.getElementById('energyValue');
const loadStatus = document.getElementById('loadStatus');
const bindingResult = document.getElementById('bindingResult');

// Loading overlay elements
const loader = document.getElementById('loader');
const loaderLabel = document.getElementById('loaderLabel');
const loaderFill = document.getElementById('loaderFill');
const loaderPct = document.getElementById('loaderPct');

// ============================================================================
// EVENT LISTENERS
// ============================================================================

loadBtn.addEventListener('click', async () => {
  const pdbId = pdbIdInput.value.trim().toUpperCase();

  if (!pdbId || pdbId.length !== 4) {
    showStatus('Invalid PDB ID (must be 4 characters)', 'error');
    return;
  }

  try {
    loadBtn.disabled = true;
    loadBtn.textContent = '⏳ Loading...';
    showLoader('Connecting…', 2);
    await nextFrame();

    // 1) Download (0–60%) — cached structures skip straight to building.
    let structure;
    if (appState.pdbFetcher.pdbCache[pdbId]) {
      structure = appState.pdbFetcher.pdbCache[pdbId];
      setProgress(70, 'Loaded from cache…');
      await nextFrame();
    } else {
      const pdbText = await appState.pdbFetcher.downloadPDB(pdbId, (frac) => {
        setProgress(2 + frac * 58, 'Downloading structure…');
      });

      // 2) Parse (60–70%)
      setProgress(62, 'Parsing atoms…');
      await nextFrame(); // let the bar paint before the blocking parse
      structure = appState.pdbFetcher.parsePDB(pdbText, pdbId);
      appState.pdbFetcher.pdbCache[pdbId] = structure;
    }
    appState.structure = structure;
    // Remember the crystal pose so every capture starts fresh from it.
    appState.initialPositions = appState.structure.atoms.map((a) => ({
      x: a.position.x, y: a.position.y, z: a.position.z,
    }));

    // 3) Build the 3D scene (70–90%)
    setProgress(72, 'Building 3D scene…');
    await nextFrame();
    appState.visualization = new VisualizationEngine(canvas, appState.structure);
    appState.visualization.renderStructure(appState.structure);
    appState.visualization.zoomToFit();

    // 4) Initialize the simulation (90–100%)
    setProgress(92, 'Initializing simulation…');
    await nextFrame();
    appState.physics = new PhysicsEngine(appState.structure, {
      temperature: 300,
      timestep: 0.001,
      recordTrajectory: false, // we capture manually in captureTrajectory()
    });

    setProgress(100, 'Ready');
    await nextFrame();
    hideLoader();

    showStatus(`✓ Loaded ${pdbId} (${appState.structure.atoms.length} atoms)`, 'success');
    updateMetrics();
  } catch (error) {
    hideLoader();
    showStatus(`✗ Error: ${error.message}`, 'error');
    console.error('PDB loading error:', error);
  } finally {
    loadBtn.disabled = false;
    loadBtn.textContent = 'Load PDB';
  }
});

addInhibitorBtn.addEventListener('click', () => {
  if (!appState.structure) {
    showStatus('Load a protein structure first', 'error');
    return;
  }

  try {
    const builder = new InhibitorBuilder();
    const scaffoldMap = {
      nucleoside: 'nucleoside-analog',
      protease: 'protease-inhibitor',
      allosteric: 'allosteric-inhibitor',
      generic: 'generic-small-molecule',
    };

    appState.inhibitor = builder
      .loadPreset(scaffoldMap[inhibitorType.value])
      .build();

    // Visualize inhibitor (orange carbons, CPK heteroatoms).
    if (appState.visualization) {
      appState.visualization.renderLigand(appState.inhibitor, 0xffaa00);
    }

    // Initialize binding analyzer
    appState.bindingAnalyzer = new BindingAnalyzer(appState.structure, appState.inhibitor);

    showStatus(`✓ Added inhibitor (${appState.inhibitor.atoms.length} atoms)`, 'success');
  } catch (error) {
    showStatus(`✗ Error: ${error.message}`, 'error');
    console.error('Inhibitor error:', error);
  }
});

startBtn.addEventListener('click', async () => {
  if (!appState.physics) {
    showStatus('Load a structure first', 'error');
    return;
  }
  if (appState.capturing) return;

  appState.capturing = true;
  startBtn.disabled = true;
  pauseBtn.disabled = true;

  try {
    // CAPTURE PHASE: run the whole simulation up front, recording each frame.
    // This is the compute-heavy part, shown behind the progress bar.
    const frames = await captureTrajectory();

    // PLAYBACK PHASE: hand the recorded frames to the viewer and play them
    // back smoothly. No physics runs now, so rendering stays fluid.
    appState.visualization.loadTrajectory(frames);
    appState.visualization.play();

    pauseBtn.disabled = false;
    pauseBtn.textContent = '⏸ Pause';
    showStatus('Playing recorded trajectory', 'success');
  } catch (error) {
    hideLoader();
    startBtn.disabled = false;
    showStatus(`✗ Simulation error: ${error.message}`, 'error');
    console.error('Capture error:', error);
  } finally {
    appState.capturing = false;
  }
});

pauseBtn.addEventListener('click', () => {
  const viz = appState.visualization;
  if (!viz) return;

  // Pause/resume the recorded PLAYBACK (physics already finished).
  if (viz.isPlaying) {
    viz.pause();
    pauseBtn.textContent = '▶ Resume';
    showStatus('Paused', 'info');
  } else {
    viz.play();
    pauseBtn.textContent = '⏸ Pause';
    showStatus('Playing', 'info');
  }
});

resetBtn.addEventListener('click', () => {
  // Stop playback and restore the structure to its crystal pose.
  if (appState.visualization) appState.visualization.resetTrajectory?.();
  if (appState.physics) {
    appState.physics.stop();
    appState.physics.reset();
  }
  appState.isRunning = false;

  startBtn.disabled = false;
  pauseBtn.disabled = true;
  pauseBtn.textContent = '⏸ Pause';
  updateMetrics();
  showStatus('Simulation reset', 'info');
});

viewMode.addEventListener('change', (e) => {
  if (appState.visualization) {
    appState.visualization.setViewMode(e.target.value);
  }
});

resetCameraBtn.addEventListener('click', () => {
  if (appState.visualization) {
    appState.visualization.resetCamera();
  }
});

zoomFitBtn.addEventListener('click', () => {
  if (appState.visualization) {
    appState.visualization.zoomToFit();
  }
});

analyzeBtn.addEventListener('click', () => {
  if (!appState.bindingAnalyzer) {
    showStatus('Add an inhibitor first', 'error');
    return;
  }

  try {
    const results = appState.bindingAnalyzer.analyzeBinding();

    let html = `
      <strong>Binding Analysis</strong><br>
      Score: ${results.bindingScore}/10<br>
      H-Bonds: ${results.metrics.hydrogenBonds.count}<br>
      ΔG: ${results.estimatedDeltaG.toFixed(2)} kcal/mol<br>
      Kd: ${results.estimatedKd.toExponential(1)} M<br>
      <em>${results.prediction}</em>
    `;

    bindingResult.innerHTML = html;
    bindingResult.style.display = 'block';

    // Visualize interactions. These setters are additive on the engine, so
    // drawing electrostatics does not wipe the hydrogen bonds.
    if (appState.visualization) {
      const m = results.metrics;
      appState.visualization.drawHydrogenBonds(m.hydrogenBonds.bonds);
      appState.visualization.drawElectrostaticInteractions(
        m.electrostaticInteractions.interactions
      );
      // Hydrophobic contacts if the analyzer provides them.
      const hydro =
        m.hydrophobicInteractions?.contacts ||
        m.hydrophobicContacts?.contacts ||
        m.hydrophobic?.contacts;
      if (hydro) appState.visualization.drawHydrophobicInteractions(hydro);
    }

    console.log('Binding results:', results);
  } catch (error) {
    showStatus(`✗ Analysis error: ${error.message}`, 'error');
    console.error('Analysis error:', error);
  }
});

exportJsonBtn.addEventListener('click', () => {
  if (!appState.physics) {
    showStatus('Run simulation first', 'error');
    return;
  }

  try {
    const trajectory = appState.physics.exportTrajectory();
    const data = JSON.stringify(trajectory, null, 2);
    downloadFile(data, 'trajectory.json', 'application/json');
    showStatus('✓ Exported JSON', 'success');
  } catch (error) {
    showStatus(`✗ Export error: ${error.message}`, 'error');
  }
});

exportCsvBtn.addEventListener('click', () => {
  if (!appState.bindingAnalyzer) {
    showStatus('Run analysis first', 'error');
    return;
  }

  try {
    const csv = appState.bindingAnalyzer.exportCSV();
    downloadFile(csv, 'binding-results.csv', 'text/csv');
    showStatus('✓ Exported CSV', 'success');
  } catch (error) {
    showStatus(`✗ Export error: ${error.message}`, 'error');
  }
});

exportPdbBtn.addEventListener('click', () => {
  if (!appState.structure) {
    showStatus('Load structure first', 'error');
    return;
  }

  try {
    // Simple PDB export (structure only, no trajectory)
    const pdbLines = [];
    appState.structure.atoms.forEach((atom, idx) => {
      const line = `ATOM  ${String(atom.serial).padStart(5)}  ${String(atom.name).padEnd(4)}${String(atom.resName).padEnd(3)} ${atom.chainId}${String(atom.resSeq).padStart(4)}    ${Number(atom.position.x).toFixed(3).padStart(8)}${Number(atom.position.y).toFixed(3).padStart(8)}${Number(atom.position.z).toFixed(3).padStart(8)}  1.00  0.00           ${atom.element.padEnd(2)}\n`;
      pdbLines.push(line);
    });
    pdbLines.push('END\n');

    downloadFile(pdbLines.join(''), 'structure.pdb', 'text/plain');
    showStatus('✓ Exported PDB', 'success');
  } catch (error) {
    showStatus(`✗ Export error: ${error.message}`, 'error');
  }
});

// ============================================================================
// TRAJECTORY CAPTURE (compute up front, then play back smoothly)
// ============================================================================

// Tunables: total recorded frames and physics steps advanced per frame.
// More frames = longer/smoother playback but a longer capture.
const CAPTURE_FRAMES = 200;
const STEPS_PER_FRAME = 2;

/**
 * Run the whole simulation synchronously in chunks, recording one snapshot per
 * frame, and return the trajectory. The chunked `await nextFrame()` keeps the
 * progress bar animating and the tab responsive during the compute-heavy part.
 * This is what removes the choppiness: all the physics happens here, before any
 * playback, instead of competing with rendering frame-by-frame.
 */
async function captureTrajectory() {
  const physics = appState.physics;

  // Start every capture from the crystal pose with fresh velocities, so runs
  // are deterministic and don't drift after repeated Start presses.
  if (appState.initialPositions) {
    physics.atoms.forEach((a, i) => {
      const p = appState.initialPositions[i];
      if (p) {
        a.position.x = p.x;
        a.position.y = p.y;
        a.position.z = p.z;
      }
    });
  }
  physics.reset();

  const snapshot = () => ({
    positions: physics.atoms.map((a) => ({
      x: a.position.x, y: a.position.y, z: a.position.z,
    })),
  });

  const frames = [snapshot()]; // frame 0 = starting pose
  showLoader('Simulating…', 0);
  await nextFrame();

  for (let f = 1; f <= CAPTURE_FRAMES; f++) {
    for (let s = 0; s < STEPS_PER_FRAME; s++) physics.step();
    frames.push(snapshot());

    // Yield periodically so the bar paints and the page doesn't freeze.
    if (f % 4 === 0 || f === CAPTURE_FRAMES) {
      setProgress((f / CAPTURE_FRAMES) * 100, `Simulating… ${f}/${CAPTURE_FRAMES}`);
      updateMetrics();
      await nextFrame();
    }
  }

  hideLoader();
  return frames;
}

// ============================================================================
// LOADING OVERLAY
// ============================================================================

/** Resolve after the browser has painted (double rAF), so the bar updates. */
function nextFrame() {
  return new Promise((resolve) =>
    requestAnimationFrame(() => requestAnimationFrame(resolve))
  );
}

function showLoader(label, pct) {
  loader.style.display = 'flex';
  setProgress(pct || 0, label || 'Loading…');
}

function hideLoader() {
  loader.style.display = 'none';
}

function setProgress(pct, label) {
  const p = Math.max(0, Math.min(100, pct));
  loaderFill.style.width = p + '%';
  loaderPct.textContent = Math.round(p) + '%';
  if (label) loaderLabel.textContent = label;
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

function updateMetrics() {
  if (!appState.physics) return;

  const stats = appState.physics.getStats();
  simTime.textContent = `${stats.simTime.toFixed(2)} ps`;
  stepCount.textContent = stats.stepCount;
  tempValue.textContent = `${stats.temperature.toFixed(0)} K`;
  energyValue.textContent = `${stats.energy.toFixed(2)} kcal/mol`;
}

function showStatus(message, type = 'info') {
  loadStatus.textContent = message;
  loadStatus.className = `status ${type}`;
  loadStatus.style.display = 'block';

  if (type === 'success' || type === 'error') {
    setTimeout(() => {
      loadStatus.style.display = 'none';
    }, 5000);
  }
}

function downloadFile(content, filename, type) {
  const blob = new Blob([content], { type });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

// ============================================================================
// INITIALIZATION
// ============================================================================

console.log('✓ Rhinovirus 3Dpol MD Simulation Ready');
console.log('  1. Enter PDB ID and load structure (default 1XR5 = HRV14 3Dpol)');
console.log('  2. Add inhibitor');
console.log('  3. Start simulation (view plays back the recorded trajectory)');
console.log('  4. Analyze binding');
