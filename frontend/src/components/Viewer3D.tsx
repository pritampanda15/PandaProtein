import React, { useEffect, useRef, useState } from 'react';
import * as $3Dmol from '3dmol';
import {
  Box,
  Camera,
  Eye,
  Flame,
  Focus,
  HelpCircle,
  Maximize2,
  Minimize2,
  Ruler,
  Trash2,
  X,
} from 'lucide-react';
import type {
  DockingLiveOverlay,
  HoleAnalysisData,
  MembraneInfo,
  MoleChannel,
  MoleColorBy,
  SelectedResidue,
  StructureSessionState,
} from '../types/pdbprep';

export interface ChainHighlight {
  ids: string[];
  tone: 'select' | 'delete';
}

interface DistanceMeasurement {
  id: number;
  label1: string;
  label2: string;
  p1: [number, number, number];
  p2: [number, number, number];
  dist: number;
}

interface PendingMeasureAtom {
  label: string;
  pos: [number, number, number];
}

interface Viewer3DProps {
  session: StructureSessionState;
  darkMode: boolean;
  selectedResidue: SelectedResidue | null;
  onSelectResidue: (res: SelectedResidue) => void;
  focusTrigger: number;
  highlightedChains: ChainHighlight | null;
  chainFocusTrigger: number;
  onClearChainHighlight: () => void;
  moleChannels?: MoleChannel[];
  moleColorBy?: MoleColorBy;
  holeData?: HoleAnalysisData | null;
  dockingOverlay?: DockingLiveOverlay | null;
  dockingFocusTrigger?: number;
  diffOverlayPdb?: string | null;
  membraneOverlay?: MembraneInfo | null;
  expanded: boolean;
  onToggleExpand: () => void;
}

type RepStyle = 'cartoon' | 'sticks' | 'surface' | 'ballstick';
type ColorMode = 'chain' | 'ss' | 'bfactor' | 'issues';

export const Viewer3D: React.FC<Viewer3DProps> = ({
  session,
  darkMode,
  selectedResidue,
  onSelectResidue,
  focusTrigger,
  highlightedChains,
  chainFocusTrigger,
  onClearChainHighlight,
  moleChannels = [],
  moleColorBy = 'id',
  holeData = null,
  dockingOverlay = null,
  dockingFocusTrigger = 0,
  diffOverlayPdb = null,
  membraneOverlay = null,
  expanded,
  onToggleExpand,
}) => {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const viewerRef = useRef<any>(null);
  const lastZoomSessionRef = useRef<string | null>(null);
  const defaultViewRef = useRef<number[] | null>(null);
  const dockingModelRef = useRef<any>(null);
  const dockingShapesRef = useRef<any[]>([]);
  const dockingLabelsRef = useRef<any[]>([]);
  const [sceneVersion, setSceneVersion] = useState(0);

  const [rep, setRep] = useState<RepStyle>('cartoon');
  const [colorMode, setColorMode] = useState<ColorMode>('chain');
  const [showWaters, setShowWaters] = useState(true);
  const [showHydrogens, setShowHydrogens] = useState(false);
  const [showClashes, setShowClashes] = useState(true);
  const [showDockingBox, setShowDockingBox] = useState(false);
  const [showPockets, setShowPockets] = useState(false);
  const [showInteractions, setShowInteractions] = useState(true);
  const [showLegend, setShowLegend] = useState(false);
  const [measureMode, setMeasureMode] = useState(false);
  const [pendingAtom, setPendingAtom] = useState<PendingMeasureAtom | null>(null);
  const [measurements, setMeasurements] = useState<DistanceMeasurement[]>([]);

  const measureModeRef = useRef(measureMode);
  const pendingAtomRef = useRef(pendingAtom);
  measureModeRef.current = measureMode;
  pendingAtomRef.current = pendingAtom;

  // Clear distance measurements when switching to a brand-new session
  useEffect(() => {
    setMeasurements([]);
    setPendingAtom(null);
  }, [session.session_id]);

  // Initialize 3Dmol viewer once
  useEffect(() => {
    if (!containerRef.current) return;
    const bg = darkMode ? '#111116' : '#ECEEF2';
    const viewer = $3Dmol.createViewer(containerRef.current, {
      backgroundColor: bg,
      antialias: true,
    });
    viewerRef.current = viewer;

    const ro = new ResizeObserver(() => {
      try {
        viewer.resize();
        viewer.render();
      } catch {
        // ignore resize errors
      }
    });
    ro.observe(containerRef.current);

    return () => {
      ro.disconnect();
      viewerRef.current = null;
    };
  }, []);

  // Update background on theme toggle
  useEffect(() => {
    if (viewerRef.current) {
      viewerRef.current.setBackgroundColor(darkMode ? '#111116' : '#ECEEF2');
      viewerRef.current.render();
    }
  }, [darkMode]);

  // Render structure, styles, overlays, clashes & docking box
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer || !session?.pdb_text) return;

    viewer.removeAllModels();
    viewer.removeAllShapes();
    viewer.removeAllLabels();
    viewer.removeAllSurfaces();
    dockingModelRef.current = null;
    dockingShapesRef.current = [];
    dockingLabelsRef.current = [];

    // CRITICAL: 3Dmol.js strips hydrogens by default unless { keepH: true } is passed!
    const model = viewer.addModel(session.pdb_text, 'pdb', { keepH: true });

    // Issue color lookup
    const issueMap = new Map<string, string>();
    for (const r of session.validation.residues) {
      if (r.issues.includes('clash') || r.issues.includes('rama_outlier') || r.issues.includes('cis_nonpro') || r.issues.includes('chiral_outlier')) {
        issueMap.set(`${r.chain}:${r.resseq}`, '#FF3B30');
      } else if (r.issues.includes('rotamer_outlier') || r.issues.includes('missing_atoms') || r.issues.includes('altloc') || r.issues.includes('nonstandard')) {
        issueMap.set(`${r.chain}:${r.resseq}`, '#FF9500');
      }
    }

    const proteinColorSpec = (): Record<string, unknown> => {
      if (colorMode === 'chain') return { colorscheme: 'chainHetatm' };
      if (colorMode === 'ss') return { colorscheme: 'ssPyMOL' };
      if (colorMode === 'bfactor') {
        return {
          colorfunc: (atom: any) => {
            const b = atom.b || 20;
            if (b < 20) return '#007AFF';
            if (b < 35) return '#34C759';
            if (b < 55) return '#FF9500';
            return '#FF3B30';
          },
        };
      }
      // Issues overlay
      return {
        colorfunc: (atom: any) => {
          const key = `${atom.chain}:${atom.resi}`;
          return issueMap.get(key) || (darkMode ? '#6C7380' : '#A0A8B4');
        },
      };
    };

    const cSpec = proteinColorSpec();
    const proteinSel = showHydrogens ? { hetflag: false } : { hetflag: false, elem: 'H', invert: false };

    // Base protein representation (ATOM records)
    if (rep === 'cartoon') {
      viewer.setStyle({ hetflag: false }, { cartoon: { ...cSpec, thickness: 0.38, opacity: 0.96 } });
    } else if (rep === 'sticks') {
      viewer.setStyle({ hetflag: false }, { stick: { ...cSpec, radius: 0.16 } });
    } else if (rep === 'ballstick') {
      viewer.setStyle(
        { hetflag: false },
        { stick: { ...cSpec, radius: 0.12 }, sphere: { ...cSpec, scale: 0.24 } }
      );
    } else if (rep === 'surface') {
      viewer.setStyle({ hetflag: false }, { cartoon: { ...cSpec, opacity: 0.72 } });
      viewer.addSurface($3Dmol.SurfaceType.VDW, {
        opacity: 0.78,
        color: darkMode ? '#2C3E50' : '#DCE4F0',
        ...cSpec,
      }, { hetflag: false });
    }

    const hStickColor = darkMode ? '#F8FAFC' : '#CBD5E1';
    const hSphereColor = darkMode ? '#FFFFFF' : '#E2E8F0';
    if (!showHydrogens) {
      viewer.setStyle({ elem: ['H', 'D'] }, {});
    } else if (rep === 'ballstick') {
      viewer.setStyle(
        { elem: ['H', 'D'] },
        { stick: { radius: 0.10, color: hStickColor }, sphere: { scale: 0.21, color: hSphereColor } }
      );
    } else if (rep === 'sticks') {
      viewer.setStyle(
        { elem: ['H', 'D'] },
        { stick: { radius: 0.12, color: hStickColor }, sphere: { scale: 0.14, color: hSphereColor } }
      );
    } else {
      // Cartoon / Surface with H toggle explicitly ON: show hydrogens clearly
      viewer.addStyle(
        { elem: ['H', 'D'] },
        { stick: { radius: 0.09, color: hStickColor }, sphere: { scale: 0.16, color: hSphereColor } }
      );
    }

    // Always render Ligands / Cofactors crisply in licorice sticks (use `not`, NEVER `invert: true`)
    const waterAndIonNames = [
      'HOH', 'WAT', 'SOL', 'DOD',
      'ZN', 'MG', 'CA', 'FE', 'FE2', 'MN', 'CU', 'CO', 'NI', 'NA', 'K', 'CL', 'BR', 'I',
    ];
    viewer.setStyle(
      { hetflag: true, not: { resn: waterAndIonNames } },
      { stick: { radius: 0.24, colorscheme: 'cyanCarbon' }, sphere: { scale: 0.22, colorscheme: 'cyanCarbon' } }
    );

    // Metal ions as metallic spheres
    viewer.setStyle(
      { hetflag: true, resn: ['ZN', 'MG', 'CA', 'FE', 'FE2', 'MN', 'CU', 'CO', 'NI', 'NA', 'K', 'CL'] },
      { sphere: { scale: 0.65, color: '#AF52DE' } }
    );

    // Crystallographic waters
    if (showWaters) {
      viewer.setStyle({ hetflag: true, resn: ['HOH', 'WAT', 'SOL', 'DOD'] }, { sphere: { scale: 0.22, color: '#5AC8FA', opacity: 0.75 } });
    } else {
      viewer.setStyle({ hetflag: true, resn: ['HOH', 'WAT', 'SOL', 'DOD'] }, {});
    }

    // Chain Spotlight Highlight (dims other chains and illuminates selected / delete-marked chains!)
    const activeChainIds = highlightedChains?.ids?.length ? highlightedChains.ids : [];
    if (activeChainIds.length > 0) {
      const isDelete = highlightedChains?.tone === 'delete';
      const spotColor = isDelete ? '#FF3B30' : '#007AFF';
      const ghostColor = darkMode ? '#52525B' : '#B8BEC8';

      // Dim non-selected chains using `not: { chain }` (NEVER `invert: true`)
      viewer.setStyle(
        { hetflag: false, not: { chain: activeChainIds } },
        { cartoon: { color: ghostColor, opacity: 0.22, thickness: 0.28 } }
      );
      viewer.setStyle(
        { hetflag: true, not: { chain: activeChainIds } },
        { stick: { color: ghostColor, radius: 0.14, opacity: 0.22 }, sphere: { color: ghostColor, scale: 0.18, opacity: 0.22 } }
      );

      // Spotlight selected chain(s) in vivid high-contrast color + bold representation
      const hlProtSel = { chain: activeChainIds, hetflag: false };
      if (rep === 'sticks' || rep === 'ballstick') {
        viewer.setStyle(hlProtSel, { stick: { color: spotColor, radius: 0.22 } });
      } else {
        viewer.setStyle(hlProtSel, {
          cartoon: { color: spotColor, thickness: 0.5, opacity: 1.0 },
        });
      }
      // Also spotlight heterogens/waters/ions belonging to the selected chain(s)
      viewer.setStyle(
        { chain: activeChainIds, hetflag: true },
        {
          stick: { color: spotColor, radius: 0.26 },
          sphere: { color: spotColor, scale: 0.36, opacity: 0.95 },
        }
      );

      // Floating 3D label centered on each highlighted chain
      for (const cid of activeChainIds) {
        viewer.addLabel(
          isDelete ? `Delete: Chain ${cid}` : `Chain ${cid}`,
          {
            font: '-apple-system, sans-serif',
            fontSize: 13,
            fontColor: '#FFFFFF',
            backgroundColor: isDelete ? 'rgba(255, 59, 48, 0.92)' : 'rgba(0, 122, 255, 0.92)',
            borderRadius: 8,
            backgroundOpacity: 0.92,
            inFront: true,
          },
          { chain: cid }
        );
      }
    }

    // Highlight selected residue with high-contrast sticks + halo label
    if (selectedResidue) {
      const sel = { chain: selectedResidue.chain, resi: selectedResidue.resseq };
      viewer.addStyle(sel, {
        stick: { radius: 0.24, colorscheme: 'orangeCarbon' },
        sphere: { scale: 0.22, colorscheme: 'orangeCarbon' },
      });
      // Also show neighbors within 4.5 Å in thin translucent sticks (COOT local neighborhood!)
      viewer.addStyle(
        { within: { distance: 4.5, sel }, hetflag: false },
        { line: { linewidth: 1.5, colorscheme: 'default' } }
      );
      viewer.addLabel(
        `${selectedResidue.resname || ''} ${selectedResidue.chain}:${selectedResidue.resseq}`,
        {
          font: '-apple-system, sans-serif',
          fontSize: 12,
          fontColor: '#FFFFFF',
          backgroundColor: 'rgba(0, 122, 255, 0.88)',
          borderRadius: 6,
          backgroundOpacity: 0.9,
          inFront: true,
        },
        { chain: selectedResidue.chain, resi: selectedResidue.resseq, atom: 'CA' }
      );
    }

    // Render 3D Steric Clash Markers (red glowing spheres at clash midpoints)
    if (showClashes && session.validation.clashes.clashes.length > 0) {
      const topClashes = session.validation.clashes.clashes.slice(0, 25);
      for (const cl of topClashes) {
        viewer.addSphere({
          center: { x: cl.midpoint[0], y: cl.midpoint[1], z: cl.midpoint[2] },
          radius: cl.severity === 'severe' ? 0.55 : 0.38,
          color: '#FF3B30',
          opacity: 0.65,
        });
      }
    }

    // Render 3D AutoDock Vina Grid Box wireframe
    if (showDockingBox && session.ligands_info.docking_boxes.length > 0) {
      const box = session.ligands_info.docking_boxes[0];
      viewer.addBox({
        center: { x: box.center[0], y: box.center[1], z: box.center[2] },
        dimensions: { w: box.size[0], h: box.size[1], d: box.size[2] },
        color: '#007AFF',
        wireframe: true,
        opacity: 0.85,
      });
    }

    // Render 3D Detected Druggable Binding Pockets (translucent alpha-sphere clouds)
    if (showPockets && session.pockets && session.pockets.length > 0) {
      const palette = ['#5AC8FA', '#AF52DE', '#34C759', '#FF9500', '#007AFF'];
      session.pockets.slice(0, 3).forEach((pk, pi) => {
        const col = palette[pi % palette.length];
        for (const sp of pk.spheres) {
          viewer.addSphere({
            center: { x: sp[0], y: sp[1], z: sp[2] },
            radius: pk.probe_radius || 1.15,
            color: col,
            opacity: 0.34,
          });
        }
        viewer.addLabel(
          `Pocket #${pk.rank} (${pk.volume_A3} Å³)`,
          {
            font: '-apple-system, sans-serif',
            fontSize: 11,
            fontColor: '#FFFFFF',
            backgroundColor: 'rgba(88, 86, 214, 0.88)',
            borderRadius: 6,
            inFront: true,
          },
          { x: pk.center[0], y: pk.center[1], z: pk.center[2] }
        );
      });
    }

    // Render 3D Disulfide Bridges (gold), Metal Coordination (purple), and Ligand H-Bonds (cyan dashed)
    if (showInteractions && session.overlays) {
      for (const ss of session.overlays.disulfides) {
        viewer.addCylinder({
          start: { x: ss.p1[0], y: ss.p1[1], z: ss.p1[2] },
          end: { x: ss.p2[0], y: ss.p2[1], z: ss.p2[2] },
          radius: 0.16,
          color: '#FFCC00',
          fromCap: 1,
          toCap: 1,
        });
      }
      for (const mb of session.overlays.metal_bonds) {
        viewer.addCylinder({
          start: { x: mb.p1[0], y: mb.p1[1], z: mb.p1[2] },
          end: { x: mb.p2[0], y: mb.p2[1], z: mb.p2[2] },
          radius: 0.08,
          dashed: true,
          color: '#AF52DE',
        });
      }
      for (const hb of session.overlays.ligand_hbonds) {
        viewer.addCylinder({
          start: { x: hb.p1[0], y: hb.p1[1], z: hb.p1[2] },
          end: { x: hb.p2[0], y: hb.p2[1], z: hb.p2[2] },
          radius: 0.07,
          dashed: true,
          color: '#007AFF',
        });
      }
    }

    // Render MOLE 3D Tunnels & Pores (overlapping translucent spheres along centerline + bottleneck highlight!)
    if (moleChannels && moleChannels.length > 0) {
      const sliceColor = (ch: MoleChannel, pt: { hydropathy: number; polarity: number; charge: number }) => {
        if (moleColorBy === 'hydropathy') {
          if (pt.hydropathy >= 1.5) return '#FF9500';
          if (pt.hydropathy >= 0.4) return '#FFCC00';
          if (pt.hydropathy <= -1.5) return '#007AFF';
          if (pt.hydropathy <= -0.4) return '#5AC8FA';
          return '#A1A1A6';
        }
        if (moleColorBy === 'polarity') {
          if (pt.polarity >= 25) return '#007AFF';
          if (pt.polarity >= 8) return '#5AC8FA';
          return '#FF9500';
        }
        if (moleColorBy === 'charge') {
          if (pt.charge > 0) return '#007AFF';
          if (pt.charge < 0) return '#FF3B30';
          return '#D1D5DB';
        }
        return ch.color;
      };

      for (const ch of moleChannels) {
        for (const pt of ch.profile) {
          viewer.addSphere({
            center: { x: pt.x, y: pt.y, z: pt.z },
            radius: Math.max(0.45, Math.min(4.0, pt.radius)),
            color: sliceColor(ch, pt),
            opacity: 0.62,
          });
        }
        // Highlight narrowest bottleneck in bright red
        const [bx, by, bz] = ch.bottleneck_point;
        viewer.addSphere({
          center: { x: bx, y: by, z: bz },
          radius: Math.max(0.55, ch.bottleneck_radius_A + 0.12),
          color: '#FF3B30',
          opacity: 0.85,
        });
        viewer.addLabel(
          `${ch.id} (${ch.bottleneck_radius_A} Å)`,
          {
            font: '-apple-system, sans-serif',
            fontSize: 11,
            fontColor: '#FFFFFF',
            backgroundColor: ch.kind === 'pore' ? 'rgba(175, 82, 222, 0.90)' : 'rgba(0, 122, 255, 0.90)',
            borderRadius: 6,
            inFront: true,
          },
          { x: bx, y: by, z: bz }
        );
      }
    }

    // Render Classic HOLE 3D Pore Surface (Red < 1.15 Å | Green 1.15–2.30 Å | Blue > 2.30 Å)
    if (holeData && holeData.profile && holeData.profile.length > 0) {
      for (const pt of holeData.profile) {
        viewer.addSphere({
          center: { x: pt.x, y: pt.y, z: pt.z },
          radius: Math.max(0.32, Math.min(6.0, pt.radius_A)),
          color: pt.color,
          opacity: 0.58,
        });
      }
      const [cxH, cyH, czH] = holeData.summary.constriction_point;
      viewer.addLabel(
        `HOLE Gate: ${holeData.summary.min_radius_A} Å`,
        {
          font: '-apple-system, sans-serif',
          fontSize: 11,
          fontColor: '#FFFFFF',
          backgroundColor:
            holeData.summary.min_radius_A < 1.15
              ? 'rgba(255, 59, 48, 0.92)'
              : holeData.summary.min_radius_A <= 2.3
              ? 'rgba(52, 199, 89, 0.92)'
              : 'rgba(0, 122, 255, 0.92)',
          borderRadius: 6,
          inFront: true,
        },
        { x: cxH, y: cyH, z: czH }
      );
    }

    // Render Before vs. After / Reference Superposition Ghost Overlay (translucent silver cartoon)
    if (diffOverlayPdb) {
      const ghostModel = viewer.addModel(diffOverlayPdb, 'pdb');
      ghostModel.setStyle(
        { hetflag: false },
        { cartoon: { color: darkMode ? '#9CA3AF' : '#64748B', opacity: 0.36, thickness: 0.28 } }
      );
      ghostModel.setStyle({ hetflag: true }, {});
    }

    // Render Realistic 3D Phospholipid Bilayer (Outer/Inner Headgroup Spheres + Dual Acyl Tails + Rim Planes)
    if (membraneOverlay) {
      try {
        const [nx0, ny0, nz0] = membraneOverlay.normal;
        const nLen = Math.hypot(nx0, ny0, nz0) || 1;
        const n: [number, number, number] = [nx0 / nLen, ny0 / nLen, nz0 / nLen];
        const [cx, cy, cz] = membraneOverlay.center;
        const [tx, ty, tz] = membraneOverlay.top_leaflet_center;
        const [bx, by, bz] = membraneOverlay.bottom_leaflet_center;
        const halfT = (membraneOverlay.thickness_A || 30) / 2;
        const rDisc = Math.max(20, Math.min(48, membraneOverlay.disc_radius_A || 28));

        // Build orthonormal basis (u, v) spanning the membrane plane perpendicular to n
        const refAxis: [number, number, number] = Math.abs(n[0]) < 0.85 ? [1, 0, 0] : [0, 1, 0];
        const ux0 = refAxis[1] * n[2] - refAxis[2] * n[1];
        const uy0 = refAxis[2] * n[0] - refAxis[0] * n[2];
        const uz0 = refAxis[0] * n[1] - refAxis[1] * n[0];
        const uLen = Math.hypot(ux0, uy0, uz0) || 1;
        const u: [number, number, number] = [ux0 / uLen, uy0 / uLen, uz0 / uLen];
        const v: [number, number, number] = [
          n[1] * u[2] - n[2] * u[1],
          n[2] * u[0] - n[0] * u[2],
          n[0] * u[1] - n[1] * u[0],
        ];

        // Collect protein heavy-atom (u, v) footprints inside the membrane slab so lipids wrap AROUND the protein!
        const protAtoms: any[] =
          typeof model.selectedAtoms === 'function'
            ? model.selectedAtoms({ hetflag: false, not: { elem: ['H', 'D'] } }) || []
            : [];
        const cavityPts: Array<[number, number]> = [];
        for (let i = 0; i < protAtoms.length; i += 2) {
          const a = protAtoms[i];
          const dx = a.x - cx;
          const dy = a.y - cy;
          const dz = a.z - cz;
          const projN = dx * n[0] + dy * n[1] + dz * n[2];
          if (Math.abs(projN) <= halfT + 2.5) {
            cavityPts.push([dx * u[0] + dy * u[1] + dz * u[2], dx * v[0] + dy * v[1] + dz * v[2]]);
          }
        }

        // Synthesize a single lightweight pseudo-PDB of coarse-grained phospholipids (headgroup P + dual acyl tails)
        const fmtCoord = (val: number) =>
          (Math.abs(val) >= 999 ? val.toFixed(2) : val.toFixed(3)).padStart(8, ' ');
        const pdbLines: string[] = [];
        const conectLines: string[] = [];
        let atomSerial = 1;
        let resSeq = 1;

        const addAtomLine = (
          serial: number,
          name: string,
          resn: string,
          chain: string,
          rseq: number,
          px: number,
          py: number,
          pz: number,
          elem: string
        ) => {
          pdbLines.push(
            `HETATM${String(serial).padStart(5, ' ')} ${name.padEnd(4, ' ')} ${resn.padEnd(3, ' ')} ${chain}${String(
              rseq % 9999
            ).padStart(4, ' ')}    ${fmtCoord(px)}${fmtCoord(py)}${fmtCoord(pz)}  1.00 20.00          ${elem.padStart(2, ' ')}`
          );
        };

        const spacing = 5.4;
        const rowStep = spacing * 0.866;
        const clearSq = 4.6 * 4.6;
        const maxLipidsPerLeaflet = 150;
        let lipidCount = 0;

        const gridSites: Array<[number, number]> = [];
        for (let row = -Math.ceil(rDisc / rowStep); row <= Math.ceil(rDisc / rowStep); row++) {
          const gv = row * rowStep;
          const xOffset = (row & 1) ? spacing * 0.5 : 0;
          for (let col = -Math.ceil(rDisc / spacing); col <= Math.ceil(rDisc / spacing); col++) {
            const gu = col * spacing + xOffset;
            if (gu * gu + gv * gv > (rDisc - 1.5) * (rDisc - 1.5)) continue;
            // Carve out transmembrane protein footprint so lipids hug the protein cleanly
            let blocked = false;
            for (let k = 0; k < cavityPts.length; k++) {
              const du = gu - cavityPts[k][0];
              const dv = gv - cavityPts[k][1];
              if (du * du + dv * dv < clearSq) {
                blocked = true;
                break;
              }
            }
            if (!blocked) gridSites.push([gu, gv]);
          }
        }

        // Uniformly subsample if patch is very large
        const stride = Math.max(1, Math.ceil(gridSites.length / maxLipidsPerLeaflet));
        for (let idx = 0; idx < gridSites.length; idx += stride) {
          const [gu, gv] = gridSites[idx];
          lipidCount++;
          const phase = (gu * 0.73 + gv * 1.19) % 6.28;
          const wx = Math.cos(phase) * 0.75;
          const wy = Math.sin(phase) * 0.75;

          for (const side of [1, -1] as const) {
            const resn = side === 1 ? 'LOU' : 'LIN'; // Outer vs Inner Leaflet
            const chId = side === 1 ? 'U' : 'L';

            const toWorld = (lu: number, lv: number, ln: number): [number, number, number] => [
              cx + u[0] * lu + v[0] * lv + n[0] * ln,
              cy + u[1] * lu + v[1] * lv + n[1] * ln,
              cz + u[2] * lu + v[2] * lv + n[2] * ln,
            ];

            const pHead = toWorld(gu, gv, side * halfT);
            const pA1 = toWorld(gu + 0.9, gv + 0.3, side * (halfT - 3.6));
            const pA2 = toWorld(gu + 0.9 + wx, gv + 0.3 + wy, side * (halfT - 7.8));
            const pA3 = toWorld(gu + 0.7 - wx * 0.5, gv + 0.4 + wy * 0.5, side * (halfT - 12.2));

            const pB1 = toWorld(gu - 0.9, gv - 0.3, side * (halfT - 3.6));
            const pB2 = toWorld(gu - 0.9 - wx, gv - 0.3 - wy, side * (halfT - 7.8));
            const pB3 = toWorld(gu - 0.7 + wx * 0.5, gv - 0.4 - wy * 0.5, side * (halfT - 12.2));

            const s0 = atomSerial++;
            const s1 = atomSerial++;
            const s2 = atomSerial++;
            const s3 = atomSerial++;
            const s4 = atomSerial++;
            const s5 = atomSerial++;
            const s6 = atomSerial++;

            addAtomLine(s0, 'P', resn, chId, resSeq, pHead[0], pHead[1], pHead[2], 'P');
            addAtomLine(s1, 'C1', resn, chId, resSeq, pA1[0], pA1[1], pA1[2], 'C');
            addAtomLine(s2, 'C2', resn, chId, resSeq, pA2[0], pA2[1], pA2[2], 'C');
            addAtomLine(s3, 'C3', resn, chId, resSeq, pA3[0], pA3[1], pA3[2], 'C');
            addAtomLine(s4, 'C4', resn, chId, resSeq, pB1[0], pB1[1], pB1[2], 'C');
            addAtomLine(s5, 'C5', resn, chId, resSeq, pB2[0], pB2[1], pB2[2], 'C');
            addAtomLine(s6, 'C6', resn, chId, resSeq, pB3[0], pB3[1], pB3[2], 'C');

            conectLines.push(`CONECT${String(s0).padStart(5, ' ')}${String(s1).padStart(5, ' ')}${String(s4).padStart(5, ' ')}`);
            conectLines.push(`CONECT${String(s1).padStart(5, ' ')}${String(s2).padStart(5, ' ')}`);
            conectLines.push(`CONECT${String(s2).padStart(5, ' ')}${String(s3).padStart(5, ' ')}`);
            conectLines.push(`CONECT${String(s4).padStart(5, ' ')}${String(s5).padStart(5, ' ')}`);
            conectLines.push(`CONECT${String(s5).padStart(5, ' ')}${String(s6).padStart(5, ' ')}`);
            resSeq++;
          }
        }

        if (lipidCount > 0) {
          const lipidModel = viewer.addModel([...pdbLines, ...conectLines, 'END'].join('\n'), 'pdb', {
            assignBonds: false,
          });
          // Combine sphere + stick inside single setStyle calls (GLModel does NOT have addStyle!)
          lipidModel.setStyle(
            { resn: 'LOU', atom: 'P' },
            {
              sphere: { radius: 1.45, color: '#FF9500', opacity: 0.85 },
              stick: { radius: 0.24, color: '#F59E0B', opacity: 0.42 },
            }
          );
          lipidModel.setStyle(
            { resn: 'LIN', atom: 'P' },
            {
              sphere: { radius: 1.45, color: '#38BDF8', opacity: 0.85 },
              stick: { radius: 0.24, color: '#38BDF8', opacity: 0.42 },
            }
          );
          lipidModel.setStyle(
            { resn: ['LOU', 'LIN'], elem: 'C' },
            { stick: { radius: 0.22, color: darkMode ? '#CBD5E1' : '#94A3B8', opacity: 0.36 } }
          );
        }

        // Draw crisp outer & inner leaflet perimeter rings (thin tubular hoops)
        const ringSegs = 40;
        for (const [leafCenter, ringCol] of [
          [[tx, ty, tz] as [number, number, number], '#FF9500'],
          [[bx, by, bz] as [number, number, number], '#38BDF8'],
        ] as const) {
          for (let k = 0; k < ringSegs; k++) {
            const a1 = (2 * Math.PI * k) / ringSegs;
            const a2 = (2 * Math.PI * (k + 1)) / ringSegs;
            const p1 = {
              x: leafCenter[0] + (u[0] * Math.cos(a1) + v[0] * Math.sin(a1)) * rDisc,
              y: leafCenter[1] + (u[1] * Math.cos(a1) + v[1] * Math.sin(a1)) * rDisc,
              z: leafCenter[2] + (u[2] * Math.cos(a1) + v[2] * Math.sin(a1)) * rDisc,
            };
            const p2 = {
              x: leafCenter[0] + (u[0] * Math.cos(a2) + v[0] * Math.sin(a2)) * rDisc,
              y: leafCenter[1] + (u[1] * Math.cos(a2) + v[1] * Math.sin(a2)) * rDisc,
              z: leafCenter[2] + (u[2] * Math.cos(a2) + v[2] * Math.sin(a2)) * rDisc,
            };
            viewer.addCylinder({
              start: p1,
              end: p2,
              radius: 0.28,
              color: ringCol,
              opacity: 0.72,
              fromCap: 0,
              toCap: 0,
            });
          }
        }

        // Place clean leaflet labels on the outer rim so they never overlap the protein
        viewer.addLabel(
          `Outer Leaflet (+${halfT.toFixed(0)} Å)`,
          {
            font: '-apple-system, sans-serif',
            fontSize: 11,
            fontColor: '#FFFFFF',
            backgroundColor: 'rgba(255, 149, 0, 0.90)',
            borderRadius: 6,
            inFront: true,
          },
          { x: tx + u[0] * rDisc, y: ty + u[1] * rDisc, z: tz + u[2] * rDisc }
        );
        viewer.addLabel(
          `Inner Leaflet (−${halfT.toFixed(0)} Å)`,
          {
            font: '-apple-system, sans-serif',
            fontSize: 11,
            fontColor: '#FFFFFF',
            backgroundColor: 'rgba(0, 122, 255, 0.90)',
            borderRadius: 6,
            inFront: true,
          },
          { x: bx + u[0] * rDisc, y: by + u[1] * rDisc, z: bz + u[2] * rDisc }
        );
      } catch (err) {
        console.error('3D Lipid Bilayer render error:', err);
      }
    }

    // Render 2-Click Distance Ruler Measurements & Pending Anchor
    if (pendingAtom) {
      viewer.addSphere({
        center: { x: pendingAtom.pos[0], y: pendingAtom.pos[1], z: pendingAtom.pos[2] },
        radius: 0.48,
        color: '#FFCC00',
        opacity: 0.9,
      });
      viewer.addLabel(
        `Start: ${pendingAtom.label}`,
        {
          font: '-apple-system, sans-serif',
          fontSize: 11,
          fontColor: '#1D1D1F',
          backgroundColor: 'rgba(255, 204, 0, 0.95)',
          borderRadius: 6,
          inFront: true,
        },
        { x: pendingAtom.pos[0], y: pendingAtom.pos[1], z: pendingAtom.pos[2] }
      );
    }
    for (const m of measurements) {
      viewer.addCylinder({
        start: { x: m.p1[0], y: m.p1[1], z: m.p1[2] },
        end: { x: m.p2[0], y: m.p2[1], z: m.p2[2] },
        radius: 0.10,
        dashed: true,
        color: '#FF9500',
      });
      viewer.addSphere({
        center: { x: m.p1[0], y: m.p1[1], z: m.p1[2] },
        radius: 0.30,
        color: '#FF9500',
        opacity: 0.9,
      });
      viewer.addSphere({
        center: { x: m.p2[0], y: m.p2[1], z: m.p2[2] },
        radius: 0.30,
        color: '#FF9500',
        opacity: 0.9,
      });
      viewer.addLabel(
        `${m.dist.toFixed(2)} Å`,
        {
          font: '-apple-system, monospace',
          fontSize: 12,
          fontColor: '#FFFFFF',
          backgroundColor: 'rgba(255, 149, 0, 0.94)',
          borderRadius: 6,
          inFront: true,
        },
        {
          x: (m.p1[0] + m.p2[0]) / 2,
          y: (m.p1[1] + m.p2[1]) / 2,
          z: (m.p1[2] + m.p2[2]) / 2,
        }
      );
    }

    // Click handler: either 2-Click Distance Ruler OR COOT Residue Selection
    model.setClickable({}, true, (atom: any) => {
      if (!atom) return;
      if (measureModeRef.current) {
        const pos: [number, number, number] = [Number(atom.x), Number(atom.y), Number(atom.z)];
        const label = `${atom.resn} ${atom.chain}:${atom.resi}.${atom.atom}`;
        const prev = pendingAtomRef.current;
        if (!prev) {
          setPendingAtom({ label, pos });
        } else {
          const dx = pos[0] - prev.pos[0];
          const dy = pos[1] - prev.pos[1];
          const dz = pos[2] - prev.pos[2];
          const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
          setMeasurements((cur) => [
            ...cur,
            { id: Date.now(), label1: prev.label, label2: label, p1: prev.pos, p2: pos, dist },
          ]);
          setPendingAtom(null);
        }
        return;
      }
      onSelectResidue({
        chain: atom.chain,
        resseq: atom.resi,
        icode: '',
        resname: atom.resn,
        key: `${atom.chain}:${atom.resi}`,
      });
    });

    if (lastZoomSessionRef.current !== session.session_id) {
      lastZoomSessionRef.current = session.session_id;
      viewer.setView([0, 0, 0, 0, 0, 0, 0, 1]);
      viewer.zoomTo({ model: 0 });
      try {
        defaultViewRef.current = viewer.getView();
      } catch {
        defaultViewRef.current = null;
      }
    }
    viewer.render();
    setSceneVersion((v) => v + 1);
  }, [
    session.pdb_text,
    session.session_id,
    rep,
    colorMode,
    showWaters,
    showHydrogens,
    showClashes,
    showDockingBox,
    showPockets,
    showInteractions,
    selectedResidue?.chain,
    selectedResidue?.resseq,
    highlightedChains?.tone,
    highlightedChains?.ids?.join(','),
    moleChannels.map((c) => c.id).join(','),
    moleColorBy,
    holeData ? `${holeData.axis_label}_${holeData.profile.length}_${holeData.summary.min_radius_A}` : '',
    diffOverlayPdb,
    membraneOverlay ? `${membraneOverlay.tilt_from_z_deg}_${membraneOverlay.top_leaflet_center.join(',')}` : '',
    pendingAtom ? pendingAtom.label : '',
    measurements.length,
    darkMode,
  ]);

  // Dedicated zero-flicker 60-FPS WebGL layer for Live 3D Molecular Docking (ligand + grid box + live H-bonds/clashes)
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;

    // Remove ONLY previous docking layer artifacts without touching the receptor protein!
    try {
      if (dockingModelRef.current) {
        viewer.removeModel(dockingModelRef.current);
        dockingModelRef.current = null;
      }
      for (const sh of dockingShapesRef.current) {
        viewer.removeShape(sh);
      }
      dockingShapesRef.current = [];
      for (const lb of dockingLabelsRef.current) {
        viewer.removeLabel(lb);
      }
      dockingLabelsRef.current = [];

      if (!dockingOverlay) {
        viewer.render();
        return;
      }

      const [bx, by, bz] = dockingOverlay.boxCenter;
      const [sx, sy, sz] = dockingOverlay.boxSize;
      const boxCol = dockingOverlay.isNewBest
        ? '#34C759'
        : dockingOverlay.running
        ? '#0A84FF'
        : '#007AFF';

      if (dockingOverlay.showBox) {
        const boxShape = viewer.addBox({
          center: { x: bx, y: by, z: bz },
          dimensions: { w: sx, h: sy, d: sz },
          color: boxCol,
          wireframe: true,
          opacity: 0.9,
        });
        if (boxShape) dockingShapesRef.current.push(boxShape);

        // Crisp glowing corner nodes on the docking grid box
        const hx = sx / 2;
        const hy = sy / 2;
        const hz = sz / 2;
        for (const dx of [-hx, hx]) {
          for (const dy of [-hy, hy]) {
            for (const dz of [-hz, hz]) {
              const sp = viewer.addSphere({
                center: { x: bx + dx, y: by + dy, z: bz + dz },
                radius: 0.35,
                color: boxCol,
                opacity: 0.85,
              });
              if (sp) dockingShapesRef.current.push(sp);
            }
          }
        }
      }

      if (dockingOverlay.posePdb) {
        const ligModel = viewer.addModel(dockingOverlay.posePdb, 'pdb', { keepH: true });
        dockingModelRef.current = ligModel;
        const ligColorScheme = dockingOverlay.isNewBest ? 'greenCarbon' : 'magentaCarbon';
        ligModel.setStyle(
          {},
          {
            stick: { radius: 0.24, colorscheme: ligColorScheme },
            sphere: { scale: 0.26, colorscheme: ligColorScheme },
          }
        );
        // Highlight polar hydrogens cleanly
        ligModel.setStyle(
          { elem: 'H' },
          {
            stick: { radius: 0.11, color: '#FFFFFF' },
            sphere: { scale: 0.16, color: '#FFFFFF' },
          }
        );

        // Live Hydrogen Bonds (glowing cyan dashed cylinders)
        for (const hb of dockingOverlay.hbonds || []) {
          const cyl = viewer.addCylinder({
            start: { x: hb.p1[0], y: hb.p1[1], z: hb.p1[2] },
            end: { x: hb.p2[0], y: hb.p2[1], z: hb.p2[2] },
            radius: 0.10,
            dashed: true,
            color: '#00D4FF',
          });
          if (cyl) dockingShapesRef.current.push(cyl);
        }

        // Live Steric Clash Flashes (red dashed cylinders when trial pose clips a wall)
        for (const cl of dockingOverlay.clashes || []) {
          const cyl = viewer.addCylinder({
            start: { x: cl.p1[0], y: cl.p1[1], z: cl.p1[2] },
            end: { x: cl.p2[0], y: cl.p2[1], z: cl.p2[2] },
            radius: 0.12,
            dashed: true,
            color: '#FF3B30',
          });
          if (cyl) dockingShapesRef.current.push(cyl);
        }

        if (dockingOverlay.score !== null && dockingOverlay.score !== undefined) {
          const lblText =
            dockingOverlay.label ||
            `${dockingOverlay.isNewBest ? '★ NEW BEST ' : ''}ΔG = ${dockingOverlay.score.toFixed(2)} kcal/mol`;
          const lb = viewer.addLabel(
            lblText,
            {
              font: '-apple-system, monospace',
              fontSize: 12,
              fontColor: '#FFFFFF',
              backgroundColor: dockingOverlay.isNewBest
                ? 'rgba(52, 199, 89, 0.94)'
                : dockingOverlay.running
                ? 'rgba(0, 122, 255, 0.92)'
                : 'rgba(88, 86, 214, 0.92)',
              borderRadius: 8,
              inFront: true,
            },
            { x: bx, y: by + sy / 2 + 1.4, z: bz }
          );
          if (lb) dockingLabelsRef.current.push(lb);
        }
      }

      viewer.render();
    } catch (err) {
      console.error('Live docking 3D overlay error:', err);
    }
  }, [dockingOverlay, sceneVersion, darkMode]);

  // Smooth camera zoom directly into the active 3D Docking Grid Box (using exact predicate so protein never flies off-screen!)
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer || !dockingOverlay || dockingFocusTrigger === 0) return;
    try {
      const [bx, by, bz] = dockingOverlay.boxCenter;
      const rad = Math.max(...dockingOverlay.boxSize) * 0.65 + 4.0;
      const radSq = rad * rad;
      const pocketSel = {
        model: 0,
        predicate: (atom: any) => {
          const dx = atom.x - bx;
          const dy = atom.y - by;
          const dz = atom.z - bz;
          return dx * dx + dy * dy + dz * dz <= radSq;
        },
      };
      const matched = typeof viewer.selectedAtoms === 'function' ? viewer.selectedAtoms(pocketSel) : [];
      if (matched && matched.length > 0) {
        viewer.zoomTo(pocketSel, 420);
      } else {
        viewer.zoomTo({ model: 0 }, 380);
      }
      viewer.render();
    } catch {
      try {
        viewer.zoomTo({ model: 0 }, 380);
        viewer.render();
      } catch {
        // ignore
      }
    }
  }, [dockingFocusTrigger]);

  // Smoothly zoom to highlighted chain(s) when triggered
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer || chainFocusTrigger === 0) return;
    try {
      if (highlightedChains && highlightedChains.ids.length > 0) {
        viewer.zoomTo({ chain: highlightedChains.ids }, 420);
      } else {
        viewer.zoomTo({}, 380);
      }
      viewer.render();
    } catch {
      // ignore
    }
  }, [chainFocusTrigger, highlightedChains?.ids?.join(',')]);

  // Smoothly zoom to selected residue when explicitly triggered
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer || !selectedResidue || focusTrigger === 0) return;
    try {
      viewer.zoomTo({ chain: selectedResidue.chain, resi: selectedResidue.resseq }, 450);
      viewer.render();
    } catch {
      // ignore
    }
  }, [focusTrigger, selectedResidue?.chain, selectedResidue?.resseq]);

  const handleReset = () => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    try {
      if (defaultViewRef.current) {
        viewer.setView(defaultViewRef.current);
      } else {
        viewer.setView([0, 0, 0, 0, 0, 0, 0, 1]);
      }
      viewer.zoomTo({ model: 0 }, 320);
      viewer.render();
    } catch {
      viewer.zoomTo({}, 320);
      viewer.render();
    }
  };

  // Press "R" anywhere outside text inputs to Reframe back to the default view
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      const tag = (e.target as HTMLElement)?.tagName;
      if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return;
      if (e.key.toLowerCase() === 'r') {
        handleReset();
      }
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, []);

  const handleSnapshot = () => {
    if (!viewerRef.current) return;
    const png = viewerRef.current.pngURI();
    const a = document.createElement('a');
    a.href = png;
    a.download = `${session.name.replace(/\.\w+$/, '')}_3d.png`;
    a.click();
  };

  return (
    <div className="relative flex flex-col h-full w-full overflow-hidden rounded-2xl border border-black/[0.07] dark:border-white/[0.09] shadow-apple-sm bg-[#ECEEF2] dark:bg-[#111116]">
      {/* Top Floating Apple Glass HUD */}
      <div className="absolute top-3 left-3 right-3 z-10 flex flex-wrap items-center justify-between gap-2 pointer-events-none">
        <div className="flex flex-wrap items-center gap-1.5 pointer-events-auto">
          {/* Representation Segmented Control */}
          <div className="apple-glass rounded-xl p-0.5 flex items-center shadow-apple-sm">
            {(['cartoon', 'sticks', 'ballstick', 'surface'] as RepStyle[]).map((r) => (
              <button
                key={r}
                onClick={() => {
                  setRep(r);
                  if (r === 'sticks' || r === 'ballstick') {
                    setShowHydrogens(true);
                  }
                }}
                className={`px-2.5 py-1 rounded-[9px] text-[11.5px] font-medium capitalize transition ${
                  rep === r
                    ? 'bg-[#007AFF] text-white shadow-sm'
                    : 'text-[#1D1D1F] dark:text-[#E5E5EA] hover:bg-black/5 dark:hover:bg-white/10'
                }`}
              >
                {r === 'ballstick' ? 'Ball & Stick' : r}
              </button>
            ))}
          </div>

          {/* Color Scheme Selector */}
          <select
            value={colorMode}
            onChange={(e) => setColorMode(e.target.value as ColorMode)}
            className="apple-glass rounded-xl px-2.5 py-1.5 text-[11.5px] font-medium shadow-apple-sm focus:outline-none"
          >
            <option value="chain">Color: Chain</option>
            <option value="ss">Color: Secondary Structure</option>
            <option value="bfactor">Color: B-Factor / pLDDT</option>
            <option value="issues">Color: QA Issues Heatmap</option>
          </select>
        </div>

        {/* Right Utility Controls */}
        <div className="flex items-center gap-1 pointer-events-auto apple-glass rounded-xl p-1 shadow-apple-sm">
          <button
            onClick={() => setShowClashes((v) => !v)}
            title="Toggle 3D Steric Clash Spheres"
            className={`px-2 py-1 rounded-lg text-[11px] font-medium flex items-center gap-1 transition ${
              showClashes ? 'bg-[#FF3B30]/15 text-[#D70015] dark:text-[#FF453A]' : 'text-[#86868B]'
            }`}
          >
            <Flame className="w-3.5 h-3.5" />
            <span>Clashes</span>
          </button>
          <button
            onClick={() => setShowDockingBox((v) => !v)}
            title="Toggle AutoDock Vina Grid Box"
            className={`px-2 py-1 rounded-lg text-[11px] font-medium flex items-center gap-1 transition ${
              showDockingBox ? 'bg-[#007AFF]/15 text-[#007AFF] dark:text-[#0A84FF]' : 'text-[#86868B]'
            }`}
          >
            <Box className="w-3.5 h-3.5" />
            <span>Grid</span>
          </button>
          <button
            onClick={() => setShowPockets((v) => !v)}
            title="Toggle 3D Detected Druggable Binding Pockets"
            className={`px-2 py-1 rounded-lg text-[11px] font-medium transition ${
              showPockets ? 'bg-[#5856D6]/20 text-[#5856D6] dark:text-[#7D7AFF]' : 'text-[#86868B]'
            }`}
          >
            Pockets
          </button>
          <button
            onClick={() => setShowInteractions((v) => !v)}
            title="Toggle Disulfide, Metal & H-Bond Cylinders"
            className={`px-2 py-1 rounded-lg text-[11px] font-medium transition ${
              showInteractions ? 'bg-[#FFCC00]/25 text-[#9E7700] dark:text-[#FFD60A]' : 'text-[#86868B]'
            }`}
          >
            S–S/HB
          </button>
          <button
            onClick={() => setShowHydrogens((v) => !v)}
            title={
              session.metadata.counts.hydrogens === 0
                ? '0 hydrogens in structure — run PROPKA or Deep Prep to build hydrogens'
                : `Toggle ${session.metadata.counts.hydrogens} Hydrogens`
            }
            className={`px-2 py-1 rounded-lg text-[11px] font-medium transition ${
              showHydrogens && session.metadata.counts.hydrogens > 0
                ? 'bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158]'
                : showHydrogens
                ? 'bg-[#FF9500]/15 text-[#B26200] dark:text-[#FF9F0A]'
                : 'text-[#86868B]'
            }`}
          >
            H ({session.metadata.counts.hydrogens})
          </button>
          <button
            onClick={() => setShowWaters((v) => !v)}
            title="Toggle Waters"
            className={`px-2 py-1 rounded-lg text-[11px] font-medium transition ${
              showWaters ? 'bg-[#5AC8FA]/20 text-[#0071A4] dark:text-[#5AC8FA]' : 'text-[#86868B]'
            }`}
          >
            H₂O
          </button>
          <button
            onClick={() => {
              setMeasureMode((v) => {
                if (v) setPendingAtom(null);
                return !v;
              });
            }}
            title="2-Click 3D Distance Ruler (click any two atoms to measure Å)"
            className={`px-2 py-1 rounded-lg text-[11px] font-medium flex items-center gap-1 transition ${
              measureMode
                ? 'bg-[#FF9500] text-white shadow-sm'
                : measurements.length > 0
                ? 'bg-[#FF9500]/15 text-[#B26200] dark:text-[#FF9F0A]'
                : 'text-[#86868B]'
            }`}
          >
            <Ruler className="w-3.5 h-3.5" />
            <span>{measurements.length > 0 ? `Ruler (${measurements.length})` : 'Ruler'}</span>
          </button>
          <button
            onClick={() => setShowLegend((v) => !v)}
            title="Show 3D Visual Overlay Legend"
            className={`p-1.5 rounded-lg transition ${
              showLegend ? 'bg-[#007AFF] text-white' : 'hover:bg-black/5 dark:hover:bg-white/10 text-[#6E6E73] dark:text-[#A1A1A6]'
            }`}
          >
            <HelpCircle className="w-3.5 h-3.5" />
          </button>
          <div className="w-px h-4 bg-black/10 dark:bg-white/10 mx-0.5" />
          <button
            onClick={handleReset}
            title="Reframe protein to default view (Shortcut: R)"
            className="px-2 py-1 rounded-lg text-[11px] font-medium flex items-center gap-1 bg-black/[0.05] dark:bg-white/[0.08] hover:bg-[#007AFF] hover:text-white transition"
          >
            <Focus className="w-3.5 h-3.5" />
            <span>Reframe</span>
          </button>
          <button onClick={handleSnapshot} title="Save PNG Snapshot" className="p-1.5 rounded-lg hover:bg-black/5 dark:hover:bg-white/10">
            <Camera className="w-3.5 h-3.5" />
          </button>
          <button onClick={onToggleExpand} title={expanded ? 'Split View' : 'Expand 3D Viewport'} className="p-1.5 rounded-lg hover:bg-black/5 dark:hover:bg-white/10">
            {expanded ? <Minimize2 className="w-3.5 h-3.5" /> : <Maximize2 className="w-3.5 h-3.5" />}
          </button>
        </div>
      </div>

      {/* Collapsible 3D Overlay Legend Popover */}
      {showLegend && (
        <div className="absolute top-16 right-3 z-20 w-64 apple-glass rounded-2xl p-3 shadow-apple-lg border border-black/10 dark:border-white/10 space-y-2 text-[11px]">
          <div className="flex items-center justify-between font-semibold text-[11.5px]">
            <span>3D Viewport Legend</span>
            <button onClick={() => setShowLegend(false)} className="p-0.5 rounded hover:bg-black/10 dark:hover:bg-white/10">
              <X className="w-3.5 h-3.5" />
            </button>
          </div>
          <div className="space-y-1.5 text-[#1D1D1F] dark:text-[#E5E5EA]">
            {[
              ['#FF3B30', 'rounded-full', 'Red/Pink Spheres — Steric Clashes'],
              ['#FFCC00', 'rounded-sm', 'Gold Rods — Disulfide Bridges (S–S)'],
              ['#AF52DE', 'rounded-sm', 'Purple Dashes — Metal Coordination'],
              ['#007AFF', 'rounded-sm', 'Blue Dashes — Ligand H-Bonds'],
              ['#5AC8FA', 'rounded-full', 'Cyan Spheres — Waters & Pockets'],
              ['#FF3B30', 'rounded-full ring-2 ring-[#007AFF]', 'Red Tunnel Sphere — MOLE Bottleneck'],
              ['#34C759', 'rounded-full ring-2 ring-[#FF3B30]', 'R/G/B Pore — HOLE (<1.15 / 1.15–2.3 / >2.3 Å)'],
              ['#9CA3AF', 'rounded-sm', 'Silver Ghost — Reference Superposition'],
              ['#FF9500', 'rounded-full', 'Orange/Cyan Lipids — 30 Å Lipid Bilayer'],
              ['#FF9500', 'rounded-full', 'Amber Dashed Line — Distance Ruler'],
            ].map(([col, shape, text]) => (
              <div key={text} className="flex items-center gap-2">
                <span className={`w-2.5 h-2.5 shrink-0 ${shape}`} style={{ backgroundColor: col }} />
                <span className="leading-tight">{text}</span>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* 3Dmol WebGL Canvas */}
      <div ref={containerRef} className="w-full h-full min-h-[300px]" />

      {/* Bottom Status Pill */}
      <div className="absolute bottom-3 left-3 right-3 z-10 flex items-center justify-between gap-2 pointer-events-none">
        <div className="apple-glass px-3 py-1 rounded-full text-[11px] text-[#6E6E73] dark:text-[#A1A1A6] shadow-apple-sm flex items-center gap-2">
          {measureMode ? (
            <>
              <Ruler className="w-3 h-3 text-[#FF9500]" />
              <span className="text-[#B26200] dark:text-[#FF9F0A] font-medium">
                {pendingAtom
                  ? `Ruler: Click 2nd atom to measure from ${pendingAtom.label}`
                  : 'Ruler Active: Click any two atoms in 3D to measure distance (Å)'}
              </span>
            </>
          ) : (
            <>
              <Eye className="w-3 h-3 text-[#007AFF]" />
              <span>Click any atom in 3D or select a chain to spotlight it</span>
            </>
          )}
        </div>
        <div className="flex items-center gap-1.5 pointer-events-auto">
          {measurements.length > 0 && (
            <button
              onClick={() => {
                setMeasurements([]);
                setPendingAtom(null);
              }}
              className="apple-glass px-2.5 py-1 rounded-full text-[11px] font-semibold text-[#B26200] dark:text-[#FF9F0A] shadow-apple-sm flex items-center gap-1 hover:bg-[#FF9500]/10"
              title="Clear all 3D distance measurements"
            >
              <Trash2 className="w-3 h-3" />
              <span>Clear {measurements.length} Dist</span>
            </button>
          )}
          {diffOverlayPdb && (
            <span className="apple-glass px-2.5 py-1 rounded-full text-[11px] font-medium text-[#6E6E73] dark:text-[#E5E5EA] shadow-apple-sm">
              Ghost Overlay ON
            </span>
          )}
          {membraneOverlay && (
            <span className="apple-glass px-2.5 py-1 rounded-full text-[11px] font-medium text-[#B26200] dark:text-[#FF9F0A] shadow-apple-sm">
              Bilayer 30 Å ON
            </span>
          )}
          {highlightedChains && highlightedChains.ids.length > 0 && (
            <button
              onClick={onClearChainHighlight}
              className={`apple-glass px-3 py-1 rounded-full text-[11.5px] font-semibold shadow-apple-sm flex items-center gap-1.5 transition hover:scale-[1.02] ${
                highlightedChains.tone === 'delete'
                  ? 'text-[#D70015] dark:text-[#FF453A] border-[#FF3B30]/40'
                  : 'text-[#007AFF] dark:text-[#0A84FF] border-[#007AFF]/40'
              }`}
              title="Clear chain spotlight and show all chains"
            >
              <span>
                {highlightedChains.tone === 'delete' ? 'Marked to Delete:' : 'Spotlighting:'}{' '}
                Chain {highlightedChains.ids.join(', ')}
              </span>
              <X className="w-3.5 h-3.5" />
            </button>
          )}
          {selectedResidue && (
            <div className="apple-glass px-3 py-1 rounded-full text-[11.5px] font-semibold text-[#007AFF] dark:text-[#0A84FF] shadow-apple-sm font-mono">
              {selectedResidue.resname || ''} {selectedResidue.chain}:{selectedResidue.resseq}
            </div>
          )}
        </div>
      </div>
    </div>
  );
};
