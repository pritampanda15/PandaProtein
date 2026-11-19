import React, { useEffect, useRef, useState } from 'react';
import * as $3Dmol from '3dmol';
import {
  Box,
  Camera,
  Eye,
  Flame,
  Maximize2,
  Minimize2,
  RotateCcw,
  X,
} from 'lucide-react';
import type { SelectedResidue, StructureSessionState } from '../types/pdbprep';

export interface ChainHighlight {
  ids: string[];
  tone: 'select' | 'delete';
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
  expanded,
  onToggleExpand,
}) => {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const viewerRef = useRef<any>(null);

  const [rep, setRep] = useState<RepStyle>('cartoon');
  const [colorMode, setColorMode] = useState<ColorMode>('chain');
  const [showWaters, setShowWaters] = useState(true);
  const [showHydrogens, setShowHydrogens] = useState(false);
  const [showClashes, setShowClashes] = useState(true);
  const [showDockingBox, setShowDockingBox] = useState(false);
  const [showPockets, setShowPockets] = useState(false);
  const [showInteractions, setShowInteractions] = useState(true);

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

    // Click handler to pick any residue directly in 3D (COOT click-to-center!)
    model.setClickable({}, true, (atom: any) => {
      if (!atom) return;
      onSelectResidue({
        chain: atom.chain,
        resseq: atom.resi,
        icode: '',
        resname: atom.resn,
        key: `${atom.chain}:${atom.resi}`,
      });
    });

    if (focusTrigger === 0 && chainFocusTrigger === 0 && !selectedResidue && activeChainIds.length === 0) {
      viewer.zoomTo();
    }
    viewer.render();
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
    darkMode,
  ]);

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
    if (viewerRef.current) {
      viewerRef.current.zoomTo({}, 350);
      viewerRef.current.render();
    }
  };

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
          <div className="w-px h-4 bg-black/10 dark:bg-white/10 mx-0.5" />
          <button onClick={handleReset} title="Reset Camera" className="p-1.5 rounded-lg hover:bg-black/5 dark:hover:bg-white/10">
            <RotateCcw className="w-3.5 h-3.5" />
          </button>
          <button onClick={handleSnapshot} title="Save PNG Snapshot" className="p-1.5 rounded-lg hover:bg-black/5 dark:hover:bg-white/10">
            <Camera className="w-3.5 h-3.5" />
          </button>
          <button onClick={onToggleExpand} title={expanded ? 'Split View' : 'Expand 3D Viewport'} className="p-1.5 rounded-lg hover:bg-black/5 dark:hover:bg-white/10">
            {expanded ? <Minimize2 className="w-3.5 h-3.5" /> : <Maximize2 className="w-3.5 h-3.5" />}
          </button>
        </div>
      </div>

      {/* 3Dmol WebGL Canvas */}
      <div ref={containerRef} className="w-full h-full min-h-[300px]" />

      {/* Bottom Status Pill */}
      <div className="absolute bottom-3 left-3 right-3 z-10 flex items-center justify-between gap-2 pointer-events-none">
        <div className="apple-glass px-3 py-1 rounded-full text-[11px] text-[#6E6E73] dark:text-[#A1A1A6] shadow-apple-sm flex items-center gap-2">
          <Eye className="w-3 h-3 text-[#007AFF]" />
          <span>Click any atom in 3D or select a chain to spotlight it</span>
        </div>
        <div className="flex items-center gap-1.5 pointer-events-auto">
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
