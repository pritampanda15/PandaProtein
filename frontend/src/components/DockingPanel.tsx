import React, { useEffect, useMemo, useRef, useState } from 'react';
import {
  Beaker,
  Box,
  CheckCircle2,
  Crosshair,
  Download,
  Flame,
  Pause,
  Play,
  RotateCcw,
  ShieldCheck,
  Sparkles,
  Square,
  Upload,
  Zap,
} from 'lucide-react';
import { api } from '../api/client';
import type {
  DockingFrameEvent,
  DockingLiveOverlay,
  DockingMode,
  DockingPresetLigand,
  PreparedDockingLigand,
  SelectedResidue,
  StructureSessionState,
} from '../types/pdbprep';

interface DockingPanelProps {
  session: StructureSessionState;
  busy: boolean;
  onStateUpdate: (state: StructureSessionState, toast: string) => void;
  onUpdateDockingOverlay: (overlay: DockingLiveOverlay | null) => void;
  onFocusDockingBox: () => void;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
  onSuccess: (msg: string) => void;
}

export const DockingPanel: React.FC<DockingPanelProps> = ({
  session,
  busy,
  onStateUpdate,
  onUpdateDockingOverlay,
  onFocusDockingBox,
  onSelectResidue,
  onSetBusy,
  onError,
  onSuccess,
}) => {
  // Step 2: Ligand State
  const [ligSourceTab, setLigSourceTab] = useState<'smiles' | 'upload' | 'receptor'>('smiles');
  const [presets, setPresets] = useState<DockingPresetLigand[]>([]);
  const [vinaAvailable, setVinaAvailable] = useState(false);
  const [smilesInput, setSmilesInput] = useState(
    'CC1(C2C1C(N(C2)C(=O)C(C(C)(C)C)NC(=O)C(F)(F)F)C(=O)NC(CC3CCNC3=O)C#N)C'
  );
  const [ligNameInput, setLigNameInput] = useState('Nirmatrelvir (Paxlovid)');
  const [preparedLigand, setPreparedLigand] = useState<PreparedDockingLigand | null>(null);
  const fileRef = useRef<HTMLInputElement | null>(null);

  // Step 3: Interactive 3D Grid Box State
  const defaultBox = useMemo(() => {
    if (session.ligands_info.docking_boxes.length > 0) {
      const b = session.ligands_info.docking_boxes[0];
      return { center: b.center as [number, number, number], size: [20, 20, 20] as [number, number, number] };
    }
    if (session.pockets && session.pockets.length > 0) {
      const p = session.pockets[0];
      return { center: p.center as [number, number, number], size: [22, 22, 22] as [number, number, number] };
    }
    return { center: [0, 0, 0] as [number, number, number], size: [22, 22, 22] as [number, number, number] };
  }, [session.session_id]);

  const [boxCenter, setBoxCenter] = useState<[number, number, number]>(defaultBox.center);
  const [boxSize, setBoxSize] = useState<[number, number, number]>(defaultBox.size);
  const [showBox3D, setShowBox3D] = useState(true);

  // Step 4: Live Action Docking Stream State
  const [speed, setSpeed] = useState<'cinematic' | 'normal' | 'fast'>('normal');
  const [exhaustiveness, setExhaustiveness] = useState<number>(8);
  const [useVinaVerify, setUseVinaVerify] = useState(false);
  const [isDocking, setIsDocking] = useState(false);
  const [statusMsg, setStatusMsg] = useState<string>('');
  const [frames, setFrames] = useState<DockingFrameEvent[]>([]);
  const [activeFrameIdx, setActiveFrameIdx] = useState<number>(-1);
  const [modes, setModes] = useState<DockingMode[]>([]);
  const [selectedMode, setSelectedMode] = useState<number | null>(null);
  const [resultId, setResultId] = useState<string | null>(null);
  const [vinaCliResult, setVinaCliResult] = useState<{ binary: string; best_affinity_kcal_mol: number | null; stdout: string } | null>(null);
  const [isReplaying, setIsReplaying] = useState(false);

  const abortRef = useRef<AbortController | null>(null);

  // Load curated ligand presets & detect optional vina CLI once
  useEffect(() => {
    api
      .getDockingPresets()
      .then((res) => {
        setPresets(res.presets || []);
        setVinaAvailable(!!res.vina_available);
      })
      .catch(() => {
        // ignore
      });
  }, []);

  // Reset grid center when switching structure session
  useEffect(() => {
    setBoxCenter(defaultBox.center);
    setBoxSize(defaultBox.size);
    setFrames([]);
    setModes([]);
    setSelectedMode(null);
    setResultId(null);
    setActiveFrameIdx(-1);
  }, [session.session_id]);

  // Keep live 3D Grid Box & Active Pose synced to Viewer3D in real time
  const activeFrame = activeFrameIdx >= 0 && frames[activeFrameIdx] ? frames[activeFrameIdx] : null;
  const activeModeObj = selectedMode !== null ? modes.find((m) => m.mode === selectedMode) || null : null;

  useEffect(() => {
    const posePdb = activeModeObj
      ? activeModeObj.pose_pdb
      : activeFrame
      ? activeFrame.pose_pdb
      : null;
    const score = activeModeObj
      ? activeModeObj.affinity_kcal_mol
      : activeFrame
      ? activeFrame.score
      : null;
    const hbonds = activeModeObj
      ? activeModeObj.hbonds
      : activeFrame
      ? activeFrame.hbonds
      : [];
    const clashes = activeModeObj
      ? activeModeObj.clashes
      : activeFrame
      ? activeFrame.clashes
      : [];
    const label = activeModeObj
      ? `Mode #${activeModeObj.mode} · ΔG = ${activeModeObj.affinity_kcal_mol.toFixed(2)} kcal/mol (${activeModeObj.ki})`
      : activeFrame
      ? `${activeFrame.is_new_best ? '★ BEST ' : ''}Step ${activeFrame.step}/${activeFrame.total_steps} · ΔG = ${activeFrame.score.toFixed(2)} kcal/mol`
      : undefined;

    onUpdateDockingOverlay({
      boxCenter,
      boxSize,
      showBox: showBox3D,
      posePdb,
      score,
      isNewBest: activeModeObj ? activeModeObj.mode === 1 : activeFrame?.is_new_best,
      running: isDocking,
      hbonds,
      clashes,
      label,
    });
  }, [
    boxCenter[0],
    boxCenter[1],
    boxCenter[2],
    boxSize[0],
    boxSize[1],
    boxSize[2],
    showBox3D,
    activeFrameIdx,
    selectedMode,
    isDocking,
    modes.length,
  ]);

  // Cleanup SSE stream on unmount
  useEffect(() => {
    return () => {
      abortRef.current?.abort();
    };
  }, []);

  // Trajectory Movie Replay interval
  useEffect(() => {
    if (!isReplaying || frames.length === 0) return;
    const timer = window.setInterval(() => {
      setActiveFrameIdx((cur) => {
        if (cur >= frames.length - 1) {
          setIsReplaying(false);
          return cur;
        }
        return cur + 1;
      });
    }, 95);
    return () => window.clearInterval(timer);
  }, [isReplaying, frames.length]);

  // Receptor Readiness Diagnostics
  const hasHydrogens = session.metadata.counts.hydrogens > 0;
  const zeroClashes = session.validation.clashes.total_clashes === 0;
  const noMissingAtoms = session.pdbfixer.missing_atoms.length === 0;
  const receptorReady = hasHydrogens && zeroClashes && noMissingAtoms;

  const handlePrepareReceptor = async () => {
    onSetBusy('Preparing receptor (PDBFixer + COOT + PROPKA Polar H)…');
    try {
      const r = await api.deepPrepPreset(session.session_id, 'docking_receptor', session.propka.ph);
      onStateUpdate(r.state, `Receptor prepared for docking (${r.actions.length} steps applied)`);
    } catch (e: any) {
      onError(e?.message || 'Receptor preparation failed');
    } finally {
      onSetBusy(null);
    }
  };

  const handleStripPocketLigands = async () => {
    onSetBusy('Clearing co-crystallized ligands & bulk waters to open binding pocket…');
    try {
      const r = await api.stripDockingLigands(session.session_id);
      onStateUpdate(r.state, 'Pocket cleared — ready for docking!');
    } catch (e: any) {
      onError(e?.message || 'Failed to clear pocket');
    } finally {
      onSetBusy(null);
    }
  };

  const handlePrepareSmiles = async (smi?: string, customName?: string) => {
    const targetSmi = (smi ?? smilesInput).trim();
    if (!targetSmi) return;
    onSetBusy('Preparing 3D Ligand (RDKit Protonation + ETKDGv3 Conformers + MMFF94)…');
    try {
      const res = await api.prepareDockingLigand(session.session_id, {
        smiles: targetSmi,
        name: customName || ligNameInput || 'Custom Ligand',
        ph: session.propka.ph,
        num_conformers: 8,
      });
      setPreparedLigand(res);
      setModes([]);
      setFrames([]);
      setSelectedMode(null);
      setActiveFrameIdx(-1);
      onSuccess(`Prepared ${res.name} (${res.conformers_generated} 3D conformers, ${res.rotatable_bonds} rot bonds)`);
    } catch (e: any) {
      onError(e?.message || 'Ligand preparation failed');
    } finally {
      onSetBusy(null);
    }
  };

  const handleUploadLigand = async (file: File) => {
    onSetBusy(`Preparing uploaded ligand ${file.name}…`);
    try {
      const res = await api.prepareDockingLigand(session.session_id, {
        file,
        ph: session.propka.ph,
        num_conformers: 8,
      });
      setPreparedLigand(res);
      setModes([]);
      setFrames([]);
      setSelectedMode(null);
      setActiveFrameIdx(-1);
      onSuccess(`Prepared ${res.name} (${res.conformers_generated} 3D conformers)`);
    } catch (e: any) {
      onError(e?.message || 'Ligand upload failed');
    } finally {
      onSetBusy(null);
    }
  };

  const handleExtractSessionLigand = async (ligId: string, box?: [number, number, number]) => {
    onSetBusy(`Extracting & preparing ${ligId} for redocking…`);
    try {
      const res = await api.prepareDockingLigand(session.session_id, {
        session_ligand_id: ligId,
        ph: session.propka.ph,
        num_conformers: 8,
      });
      setPreparedLigand(res);
      if (box) {
        setBoxCenter(box);
      }
      setModes([]);
      setFrames([]);
      setSelectedMode(null);
      setActiveFrameIdx(-1);
      onSuccess(`Prepared ${res.name} for redocking — tip: click "Open Pocket" above to remove the original copy before docking!`);
    } catch (e: any) {
      onError(e?.message || 'Extraction failed');
    } finally {
      onSetBusy(null);
    }
  };

  const handleStartLiveDocking = async () => {
    if (!preparedLigand) {
      onError('Please prepare a ligand in Step 2 first.');
      return;
    }
    abortRef.current?.abort();
    const controller = new AbortController();
    abortRef.current = controller;

    setIsDocking(true);
    setIsReplaying(false);
    setFrames([]);
    setActiveFrameIdx(-1);
    setModes([]);
    setSelectedMode(null);
    setResultId(null);
    setVinaCliResult(null);
    setShowBox3D(true);
    setStatusMsg('Initializing receptor grid & live docking stream…');
    onFocusDockingBox();

    try {
      const res = await fetch(`/api/sessions/${session.session_id}/docking/run-stream`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          ligand_key: preparedLigand.ligand_key,
          box_center: boxCenter,
          box_size: boxSize,
          exhaustiveness,
          speed,
          use_vina_verify: useVinaVerify,
          include_hetatms: false,
        }),
        signal: controller.signal,
      });

      if (!res.ok || !res.body) {
        let msg = `Docking stream failed (${res.status})`;
        try {
          const j = await res.json();
          if (j?.detail) msg = j.detail;
        } catch {
          // ignore
        }
        throw new Error(msg);
      }

      const reader = res.body.getReader();
      const decoder = new TextDecoder();
      let buffer = '';
      const accFrames: DockingFrameEvent[] = [];

      while (true) {
        const { done, value } = await reader.read();
        if (done) break;
        buffer += decoder.decode(value, { stream: true });
        const chunks = buffer.split('\n\n');
        buffer = chunks.pop() || '';

        for (const chunk of chunks) {
          const line = chunk
            .split('\n')
            .find((l) => l.startsWith('data: '));
          if (!line) continue;
          const evt = JSON.parse(line.slice(6));

          if (evt.type === 'error') {
            throw new Error(evt.message || 'Docking error');
          } else if (evt.type === 'init') {
            setStatusMsg(`Exploring ${evt.pocket_atoms} pocket atoms across ${evt.total_steps} live steps…`);
          } else if (evt.type === 'status') {
            setStatusMsg(evt.message || '');
          } else if (evt.type === 'frame') {
            accFrames.push(evt as DockingFrameEvent);
            const nextLen = accFrames.length;
            setFrames([...accFrames]);
            setActiveFrameIdx(nextLen - 1);
            setStatusMsg(evt.phase);
          } else if (evt.type === 'complete') {
            setModes(evt.modes || []);
            setResultId(evt.result_id || null);
            setVinaCliResult(evt.vina_cli || null);
            if (evt.modes && evt.modes.length > 0) {
              setSelectedMode(1);
            }
            setStatusMsg('Converged — Top binding modes ranked!');
            onSuccess(
              `Docking converged! Best Affinity: ${evt.best_affinity_kcal_mol} kcal/mol (Est. Ki ${evt.best_ki})`
            );
          }
        }
      }
    } catch (e: any) {
      if (e?.name !== 'AbortError') {
        onError(e?.message || 'Live docking failed');
        setStatusMsg('');
      } else {
        setStatusMsg('Docking paused/stopped by user');
      }
    } finally {
      setIsDocking(false);
    }
  };

  const handleStopDocking = () => {
    abortRef.current?.abort();
    setIsDocking(false);
  };

  // Current displayed telemetry (either selected ranked mode or current live/replayed frame)
  const displayTerms = activeModeObj?.terms || activeFrame?.terms || null;
  const displayScore = activeModeObj
    ? activeModeObj.affinity_kcal_mol
    : activeFrame
    ? activeFrame.score
    : null;
  const displayBest = activeModeObj
    ? modes[0]?.affinity_kcal_mol ?? activeModeObj.affinity_kcal_mol
    : activeFrame
    ? activeFrame.best_score
    : null;
  const displayKi = activeModeObj
    ? activeModeObj.ki
    : activeFrame
    ? activeFrame.ki
    : '—';
  const displayHbonds = activeModeObj?.hbonds || activeFrame?.hbonds || [];
  const displayContacts = activeModeObj?.contacts || activeFrame?.contacts || [];

  // Live SVG Energy Descent Chart coordinates
  const svgW = 500;
  const svgH = 148;
  const padL = 42;
  const padR = 14;
  const padT = 14;
  const padB = 24;

  const chartScores = frames.map((f) => Math.max(-14, Math.min(4, f.score)));
  const chartBests = frames.map((f) => Math.max(-14, Math.min(4, f.best_score)));
  const minY = Math.min(-10, ...(chartBests.length ? chartBests : [-8])) - 0.6;
  const maxY = Math.max(1.0, ...(chartScores.length ? chartScores.map((s) => Math.min(3.5, s)) : [1])) + 0.5;
  const maxStep = Math.max(10, frames[0]?.total_steps || frames.length || 40);
  const cxPlot = (st: number) => padL + (st / maxStep) * (svgW - padL - padR);
  const cyPlot = (val: number) => padT + (1 - (val - minY) / Math.max(1, maxY - minY)) * (svgH - padT - padB);

  const rawPath = frames
    .map((f, i) => `${i === 0 ? 'M' : 'L'} ${cxPlot(f.step).toFixed(1)} ${cyPlot(Math.max(minY, Math.min(maxY, f.score))).toFixed(1)}`)
    .join(' ');
  const bestPath = frames
    .map((f, i) => `${i === 0 ? 'M' : 'L'} ${cxPlot(f.step).toFixed(1)} ${cyPlot(Math.max(minY, Math.min(maxY, f.best_score))).toFixed(1)}`)
    .join(' ');

  return (
    <div className="space-y-4 pb-8">
      {/* Apple Hero Header + 4-Step Live Pipeline Ribbon */}
      <div className="apple-card p-4 space-y-3.5">
        <div className="flex flex-wrap items-center justify-between gap-3">
          <div>
            <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#007AFF] dark:text-[#0A84FF] flex items-center gap-1.5">
              <Zap className="w-3.5 h-3.5" />
              <span>Real-Time Interactive 3D Docking Studio</span>
            </div>
            <h2 className="text-[17px] font-semibold tracking-[-0.02em]">
              Watch Your Ligand Fit Live Inside the 3D Binding Pocket
            </h2>
          </div>

          <div className="flex items-center gap-1.5">
            {isDocking ? (
              <button
                onClick={handleStopDocking}
                className="apple-btn-primary !bg-[#FF3B30] hover:!bg-[#D70015]"
              >
                <Square className="w-3.5 h-3.5 fill-current" />
                <span>Stop Stream</span>
              </button>
            ) : (
              <button
                disabled={busy || !preparedLigand}
                onClick={handleStartLiveDocking}
                className="apple-btn-primary shadow-sm"
              >
                <Play className="w-3.5 h-3.5 fill-current" />
                <span>Start Live Docking</span>
              </button>
            )}
          </div>
        </div>

        {/* 4-Stage Apple Progress Stepper */}
        <div className="grid grid-cols-4 gap-1.5 text-[11px]">
          {[
            { num: 1, title: 'Receptor', ok: receptorReady, sub: hasHydrogens ? `${session.metadata.counts.hydrogens} H` : 'Add Polar H' },
            { num: 2, title: 'Ligand', ok: !!preparedLigand, sub: preparedLigand ? preparedLigand.formula : 'Upload / SMILES' },
            { num: 3, title: 'Grid Box', ok: showBox3D, sub: `${boxSize[0].toFixed(0)}×${boxSize[1].toFixed(0)}×${boxSize[2].toFixed(0)} Å` },
            { num: 4, title: 'Live Dock', ok: modes.length > 0, sub: isDocking ? 'Streaming…' : modes.length ? `${modes[0].affinity_kcal_mol} kcal` : 'Ready' },
          ].map((st) => (
            <div
              key={st.num}
              className={`px-2.5 py-2 rounded-xl border transition ${
                st.ok
                  ? 'bg-[#34C759]/10 border-[#34C759]/25'
                  : 'bg-black/[0.025] dark:bg-white/[0.035] border-transparent'
              }`}
            >
              <div className="flex items-center justify-between">
                <span className="font-semibold">
                  {st.num}. {st.title}
                </span>
                {st.ok ? (
                  <CheckCircle2 className="w-3.5 h-3.5 text-[#34C759]" />
                ) : (
                  <span className="w-2 h-2 rounded-full bg-[#86868B]/40" />
                )}
              </div>
              <div className="text-[10px] font-mono text-[#86868B] truncate mt-0.5">{st.sub}</div>
            </div>
          ))}
        </div>
      </div>

      {/* ==================== STEP 1: RECEPTOR PREPARATION ==================== */}
      <div className="apple-card p-4 space-y-3">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div className="flex items-center gap-2 text-[13.5px] font-semibold">
            <ShieldCheck className={`w-4 h-4 ${receptorReady ? 'text-[#34C759]' : 'text-[#FF9500]'}`} />
            <span>1. Receptor Preparation ({session.name})</span>
          </div>
          <div className="flex flex-wrap items-center gap-1.5">
            {session.ligands_info.ligands.length > 0 && (
              <button
                disabled={busy || isDocking}
                onClick={handleStripPocketLigands}
                className="apple-btn-secondary text-[11.5px]"
                title="Remove co-crystallized organic ligand(s) & bulk waters so the active site is open for docking"
              >
                Open Pocket (Strip {session.ligands_info.ligands.length} Co-Ligand)
              </button>
            )}
            <button
              disabled={busy || isDocking}
              onClick={handlePrepareReceptor}
              className={`text-[11.5px] ${receptorReady ? 'apple-btn-secondary' : 'apple-btn-primary'}`}
            >
              <Sparkles className="w-3.5 h-3.5" />
              <span>{receptorReady ? 'Re-Prepare Receptor' : 'Prepare Receptor (Polar H & Fix)'}</span>
            </button>
          </div>
        </div>

        <div className="flex flex-wrap gap-1.5 text-[11px] font-mono">
          <span
            className={`apple-pill ${
              hasHydrogens
                ? 'bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158]'
                : 'bg-[#FF9500]/15 text-[#B26200] dark:text-[#FF9F0A]'
            }`}
          >
            Hydrogens: {session.metadata.counts.hydrogens} {hasHydrogens ? '✓' : '(Missing)'}
          </span>
          <span
            className={`apple-pill ${
              zeroClashes
                ? 'bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158]'
                : 'bg-[#FF3B30]/15 text-[#D70015] dark:text-[#FF453A]'
            }`}
          >
            Clashes: {session.validation.clashes.total_clashes}
          </span>
          <span className="apple-pill bg-black/[0.05] dark:bg-white/[0.08] text-[#6E6E73] dark:text-[#A1A1A6]">
            PROPKA pH {session.propka.ph}
          </span>
          <span className="apple-pill bg-black/[0.05] dark:bg-white/[0.08] text-[#6E6E73] dark:text-[#A1A1A6]">
            Waters: {session.metadata.counts.waters}
          </span>
        </div>
      </div>

      {/* ==================== STEP 2: LIGAND UPLOAD & PREPARATION ==================== */}
      <div className="apple-card p-4 space-y-3">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div className="flex items-center gap-2 text-[13.5px] font-semibold">
            <Beaker className="w-4 h-4 text-[#AF52DE]" />
            <span>2. Ligand Upload &amp; 3D Preparation</span>
          </div>

          <div className="apple-segmented">
            {(
              [
                ['smiles', 'SMILES & Presets'],
                ['upload', 'Upload File'],
                ['receptor', `From PDB (${session.ligands_info.ligands.length})`],
              ] as const
            ).map(([k, lbl]) => (
              <button
                key={k}
                onClick={() => setLigSourceTab(k)}
                className={`apple-segmented-item ${ligSourceTab === k ? 'apple-segmented-active' : 'text-[#6E6E73]'}`}
              >
                {lbl}
              </button>
            ))}
          </div>
        </div>

        {ligSourceTab === 'smiles' && (
          <div className="space-y-2.5">
            {presets.length > 0 && (
              <div className="flex flex-wrap items-center gap-1.5">
                <span className="text-[11px] text-[#86868B]">1-Click Drugs:</span>
                {presets.map((p) => (
                  <button
                    key={p.id}
                    disabled={busy || isDocking}
                    onClick={() => {
                      setSmilesInput(p.smiles);
                      setLigNameInput(p.name);
                      void handlePrepareSmiles(p.smiles, p.name);
                    }}
                    className="px-2.5 py-1 rounded-full text-[11px] font-medium bg-[#AF52DE]/10 text-[#8944AB] dark:text-[#BF5AF2] hover:bg-[#AF52DE]/20 transition"
                  >
                    {p.name.split(' (')[0]}
                  </button>
                ))}
              </div>
            )}
            <div className="flex flex-wrap gap-2">
              <input
                value={smilesInput}
                onChange={(e) => setSmilesInput(e.target.value)}
                placeholder="Paste ligand SMILES (e.g. CC(=O)Oc1ccccc1C(=O)O)…"
                className="apple-input flex-1 min-w-[200px] font-mono text-[12px]"
              />
              <button
                disabled={busy || isDocking || !smilesInput.trim()}
                onClick={() => handlePrepareSmiles()}
                className="apple-btn-primary text-[12px] shrink-0"
              >
                <Sparkles className="w-3.5 h-3.5" />
                <span>Prepare 3D Ligand</span>
              </button>
            </div>
          </div>
        )}

        {ligSourceTab === 'upload' && (
          <div
            onClick={() => fileRef.current?.click()}
            className="cursor-pointer rounded-2xl border border-dashed border-black/15 dark:border-white/15 hover:border-[#007AFF] p-5 text-center transition bg-black/[0.015] dark:bg-white/[0.02]"
          >
            <input
              ref={fileRef}
              type="file"
              accept=".sdf,.mol,.mol2,.pdb,.pdbqt,.smi"
              className="hidden"
              onChange={(e) => {
                const f = e.target.files?.[0];
                if (f) void handleUploadLigand(f);
                e.target.value = '';
              }}
            />
            <Upload className="w-5 h-5 mx-auto text-[#007AFF] mb-1.5" />
            <div className="text-[12.5px] font-semibold">
              Click to Upload Ligand (.sdf, .mol2, .pdb, .pdbqt, .smi)
            </div>
            <div className="text-[11px] text-[#86868B] mt-0.5">
              Auto-protonates at pH {session.propka.ph}, builds 3D ETKDGv3 conformers &amp; MMFF94-minimizes
            </div>
          </div>
        )}

        {ligSourceTab === 'receptor' && (
          <div className="space-y-1.5">
            {session.ligands_info.ligands.length === 0 ? (
              <div className="text-[12px] text-[#86868B] py-2">
                No co-crystallized ligands in the current receptor. Use SMILES or Upload File above!
              </div>
            ) : (
              session.ligands_info.ligands.map((l) => (
                <div
                  key={l.id}
                  className="flex items-center justify-between gap-2 p-2.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.04]"
                >
                  <div className="text-[12px]">
                    <span className="font-mono font-bold">{l.resname}</span> ({l.chain}:{l.resseq}) · {l.heavy_atoms} heavy atoms
                  </div>
                  <button
                    disabled={busy || isDocking}
                    onClick={() => handleExtractSessionLigand(l.id, l.center as [number, number, number])}
                    className="apple-btn-primary !py-1 !px-2.5 text-[11.5px]"
                  >
                    Extract &amp; Prepare
                  </button>
                </div>
              ))
            )}
          </div>
        )}

        {/* Prepared Ligand Summary Card */}
        {preparedLigand && (
          <div className="p-3 rounded-2xl bg-[#AF52DE]/[0.06] border border-[#AF52DE]/20 grid grid-cols-1 md:grid-cols-12 gap-3 items-center">
            {preparedLigand.svg && (
              <div
                className="md:col-span-4 rounded-xl bg-white p-1.5 border border-black/[0.06] flex items-center justify-center overflow-hidden"
                dangerouslySetInnerHTML={{ __html: preparedLigand.svg }}
              />
            )}
            <div className={`${preparedLigand.svg ? 'md:col-span-8' : 'md:col-span-12'} space-y-2`}>
              <div className="flex flex-wrap items-center justify-between gap-1">
                <div className="text-[13px] font-bold">{preparedLigand.name}</div>
                <span className="apple-pill bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158] font-mono">
                  {preparedLigand.conformers_generated} 3D Conformers Ready
                </span>
              </div>
              <div className="grid grid-cols-3 sm:grid-cols-6 gap-1.5 text-center text-[11px]">
                {[
                  ['Formula', preparedLigand.formula],
                  ['MW', `${preparedLigand.mw}`],
                  ['cLogP', `${preparedLigand.clogp}`],
                  ['TPSA', `${preparedLigand.tpsa}`],
                  ['HBD/A', `${preparedLigand.hbd}/${preparedLigand.hba}`],
                  ['Rot Bonds', `${preparedLigand.rotatable_bonds}`],
                ].map(([k, v]) => (
                  <div key={k} className="p-1.5 rounded-lg bg-white/75 dark:bg-black/25">
                    <div className="text-[9.5px] text-[#86868B]">{k}</div>
                    <div className="font-mono font-semibold">{v}</div>
                  </div>
                ))}
              </div>
            </div>
          </div>
        )}
      </div>

      {/* ==================== STEP 3: INTERACTIVE 3D GRID BOX ==================== */}
      <div className="apple-card p-4 space-y-3">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div className="flex items-center gap-2 text-[13.5px] font-semibold">
            <Box className="w-4 h-4 text-[#007AFF]" />
            <span>3. 3D Docking Grid Box</span>
          </div>
          <div className="flex flex-wrap items-center gap-1.5">
            <button
              onClick={() => setShowBox3D((v) => !v)}
              className="apple-btn-secondary !py-1 !px-2.5 text-[11.5px]"
            >
              {showBox3D ? 'Hide Box' : 'Show Box'}
            </button>
            <button
              onClick={() => {
                setShowBox3D(true);
                onFocusDockingBox();
              }}
              className="apple-btn-secondary !py-1 !px-2.5 text-[11.5px]"
            >
              <Crosshair className="w-3.5 h-3.5 text-[#007AFF]" />
              Focus Grid in 3D
            </button>
          </div>
        </div>

        {/* Auto-Center Quick Presets */}
        <div className="flex flex-wrap items-center gap-1.5 text-[11px]">
          <span className="text-[#86868B]">Auto-Center:</span>
          {session.ligands_info.docking_boxes.map((b, idx) => (
            <button
              key={`lb-${idx}`}
              onClick={() => {
                setBoxCenter(b.center as [number, number, number]);
                setBoxSize([
                  Math.max(16, Math.min(30, b.size[0])),
                  Math.max(16, Math.min(30, b.size[1])),
                  Math.max(16, Math.min(30, b.size[2])),
                ]);
                setShowBox3D(true);
                onFocusDockingBox();
              }}
              className="px-2.5 py-1 rounded-full bg-[#007AFF]/10 text-[#007AFF] dark:text-[#0A84FF] font-medium hover:bg-[#007AFF]/20 transition"
            >
              {b.label}
            </button>
          ))}
          {(session.pockets || []).slice(0, 3).map((pk) => (
            <button
              key={pk.id}
              onClick={() => {
                setBoxCenter(pk.center as [number, number, number]);
                setBoxSize([20, 20, 20]);
                setShowBox3D(true);
                onFocusDockingBox();
              }}
              className="px-2.5 py-1 rounded-full bg-[#5856D6]/10 text-[#5856D6] dark:text-[#7D7AFF] font-medium hover:bg-[#5856D6]/20 transition"
            >
              Pocket #{pk.rank} ({pk.volume_A3} Å³)
            </button>
          ))}
        </div>

        {/* Live Sliders for Center X,Y,Z and Size X,Y,Z */}
        <div className="grid grid-cols-2 sm:grid-cols-6 gap-2.5 text-[11px]">
          {(['X', 'Y', 'Z'] as const).map((axis, idx) => (
            <div key={`c-${axis}`}>
              <div className="flex justify-between mb-1">
                <span className="text-[#86868B]">Center {axis}</span>
                <span className="font-mono font-semibold">{boxCenter[idx].toFixed(1)}</span>
              </div>
              <input
                type="number"
                step={0.5}
                value={Number(boxCenter[idx].toFixed(1))}
                onChange={(e) => {
                  const next: [number, number, number] = [...boxCenter];
                  next[idx] = parseFloat(e.target.value) || 0;
                  setBoxCenter(next);
                }}
                className="apple-input !py-1 font-mono text-[11.5px]"
              />
            </div>
          ))}
          {(['X', 'Y', 'Z'] as const).map((axis, idx) => (
            <div key={`s-${axis}`}>
              <div className="flex justify-between mb-1">
                <span className="text-[#86868B]">Size {axis}</span>
                <span className="font-mono font-semibold">{boxSize[idx].toFixed(0)} Å</span>
              </div>
              <input
                type="range"
                min={10}
                max={34}
                step={1}
                value={boxSize[idx]}
                onChange={(e) => {
                  const next: [number, number, number] = [...boxSize];
                  next[idx] = parseFloat(e.target.value);
                  setBoxSize(next);
                }}
                className="w-full apple-slider mt-1.5"
              />
            </div>
          ))}
        </div>
      </div>

      {/* ==================== STEP 4: LIVE DOCKING ENGINE & TELEMETRY ==================== */}
      <div className="apple-card p-4 space-y-4">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div>
            <div className="text-[13.5px] font-semibold flex items-center gap-2">
              <Flame className={`w-4 h-4 ${isDocking ? 'text-[#FF9500] animate-pulse' : 'text-[#007AFF]'}`} />
              <span>4. Live 3D Docking Simulation (Vina 6-Term Physics)</span>
            </div>
            {statusMsg && (
              <div className="text-[11.5px] font-mono text-[#007AFF] dark:text-[#0A84FF] mt-0.5">
                {statusMsg}
              </div>
            )}
          </div>

          <div className="flex flex-wrap items-center gap-2 text-[11.5px]">
            <select
              disabled={isDocking}
              value={speed}
              onChange={(e) => setSpeed(e.target.value as any)}
              className="apple-select !py-1 text-[11.5px]"
            >
              <option value="cinematic">Speed: Cinematic (Watch Every Move)</option>
              <option value="normal">Speed: Normal (Live)</option>
              <option value="fast">Speed: Fast</option>
            </select>
            <select
              disabled={isDocking}
              value={exhaustiveness}
              onChange={(e) => setExhaustiveness(parseInt(e.target.value, 10))}
              className="apple-select !py-1 text-[11.5px]"
            >
              <option value={4}>Quick (Exh 4)</option>
              <option value={8}>Standard (Exh 8)</option>
              <option value={12}>Deep (Exh 12)</option>
            </select>
            {vinaAvailable && (
              <label className="flex items-center gap-1.5 px-2 py-1 rounded-xl bg-black/[0.03] dark:bg-white/[0.04] cursor-pointer">
                <input
                  type="checkbox"
                  disabled={isDocking}
                  checked={useVinaVerify}
                  onChange={(e) => setUseVinaVerify(e.target.checked)}
                  className="accent-[#007AFF]"
                />
                <span>+ Vina CLI Check</span>
              </label>
            )}
          </div>
        </div>

        {/* Live Action Hero Telemetry Strip */}
        {(frames.length > 0 || modes.length > 0) && (
          <>
            <div className="grid grid-cols-2 sm:grid-cols-4 gap-2.5">
              <div className="p-3 rounded-2xl bg-black/[0.03] dark:bg-white/[0.04]">
                <div className="text-[10.5px] text-[#86868B]">
                  {activeModeObj ? `Mode #${activeModeObj.mode} Affinity` : 'Live Trial ΔG'}
                </div>
                <div className="text-[20px] font-bold font-mono tabular-nums text-[#007AFF] dark:text-[#0A84FF]">
                  {displayScore !== null ? `${displayScore.toFixed(2)}` : '—'}
                  <span className="text-[11px] font-normal text-[#86868B] ml-1">kcal/mol</span>
                </div>
              </div>
              <div className="p-3 rounded-2xl bg-[#34C759]/10 border border-[#34C759]/20">
                <div className="text-[10.5px] text-[#1B8A38] dark:text-[#30D158] font-medium">
                  Best Binding Affinity
                </div>
                <div className="text-[20px] font-bold font-mono tabular-nums text-[#1B8A38] dark:text-[#30D158]">
                  {displayBest !== null ? `${displayBest.toFixed(2)}` : '—'}
                  <span className="text-[11px] font-normal ml-1">kcal/mol</span>
                </div>
              </div>
              <div className="p-3 rounded-2xl bg-black/[0.03] dark:bg-white/[0.04]">
                <div className="text-[10.5px] text-[#86868B]">Estimated Ki (298 K)</div>
                <div className="text-[18px] font-bold font-mono tabular-nums mt-0.5">{displayKi}</div>
              </div>
              <div className="p-3 rounded-2xl bg-black/[0.03] dark:bg-white/[0.04]">
                <div className="text-[10.5px] text-[#86868B]">Live Pocket H-Bonds</div>
                <div className="text-[18px] font-bold font-mono tabular-nums text-[#00A3FF] mt-0.5">
                  {displayHbonds.length}
                  {activeFrame && isDocking && (
                    <span className="text-[10.5px] font-normal text-[#86868B] ml-1.5">
                      T={activeFrame.temperature}
                    </span>
                  )}
                </div>
              </div>
            </div>

            {/* Live Streaming SVG Affinity Trajectory Chart */}
            {frames.length > 0 && (
              <div className="rounded-2xl bg-black/[0.025] dark:bg-white/[0.03] p-2.5 space-y-2">
                <div className="flex items-center justify-between text-[11px]">
                  <span className="font-semibold">Live Binding Affinity Optimization Trajectory</span>
                  <div className="flex items-center gap-3 font-mono text-[10.5px]">
                    <span className="text-[#007AFF]">─ Trial Pose ΔG</span>
                    <span className="text-[#34C759]">━ Best-So-Far ΔG</span>
                  </div>
                </div>
                <svg viewBox={`0 0 ${svgW} ${svgH}`} className="w-full h-36 overflow-visible">
                  {/* Horizontal kcal/mol reference grid */}
                  {[-10, -8, -6, -4, -2, 0].map((val) =>
                    val >= minY && val <= maxY ? (
                      <g key={val}>
                        <line
                          x1={padL}
                          y1={cyPlot(val)}
                          x2={svgW - padR}
                          y2={cyPlot(val)}
                          stroke="currentColor"
                          className="text-black/10 dark:text-white/10"
                          strokeDasharray="3 3"
                        />
                        <text
                          x={padL - 5}
                          y={cyPlot(val) + 3}
                          textAnchor="end"
                          fontSize="8.5"
                          className="fill-[#86868B]"
                          fontFamily="monospace"
                        >
                          {val}
                        </text>
                      </g>
                    ) : null
                  )}

                  {/* Live trial curve */}
                  {rawPath && (
                    <path
                      d={rawPath}
                      fill="none"
                      stroke="#007AFF"
                      strokeWidth="1.4"
                      opacity="0.45"
                    />
                  )}
                  {/* Best-so-far descent curve */}
                  {bestPath && (
                    <path
                      d={bestPath}
                      fill="none"
                      stroke="#34C759"
                      strokeWidth="2.5"
                      strokeLinecap="round"
                    />
                  )}
                  {/* Current active step vertical marker */}
                  {activeFrame && (
                    <g>
                      <line
                        x1={cxPlot(activeFrame.step)}
                        y1={padT}
                        x2={cxPlot(activeFrame.step)}
                        y2={svgH - padB}
                        stroke="#007AFF"
                        strokeWidth="1"
                        strokeDasharray="3 2"
                        opacity="0.5"
                      />
                      <circle
                        cx={cxPlot(activeFrame.step)}
                        cy={cyPlot(Math.max(minY, Math.min(maxY, activeFrame.best_score)))}
                        r="4.5"
                        fill="#34C759"
                        stroke="#FFFFFF"
                        strokeWidth="1.5"
                      />
                    </g>
                  )}
                  <text
                    x={svgW / 2}
                    y={svgH - 4}
                    textAnchor="middle"
                    fontSize="8.5"
                    className="fill-[#86868B]"
                    fontFamily="monospace"
                  >
                    Live Docking Step (Global Exploration → Simulated Annealing → Local 6-DOF Minimization)
                  </text>
                </svg>

                {/* Trajectory Replay Movie Scrubber (available once docking stops/completes) */}
                {!isDocking && frames.length > 1 && (
                  <div className="flex items-center gap-2.5 pt-1 border-t border-black/[0.06] dark:border-white/[0.06]">
                    <button
                      onClick={() => {
                        setSelectedMode(null);
                        if (activeFrameIdx >= frames.length - 1) setActiveFrameIdx(0);
                        setIsReplaying((v) => !v);
                      }}
                      className="apple-btn-secondary !py-1 !px-2.5 text-[11px] shrink-0"
                    >
                      {isReplaying ? <Pause className="w-3 h-3" /> : <Play className="w-3 h-3 fill-current" />}
                      <span>{isReplaying ? 'Pause' : 'Replay Movie'}</span>
                    </button>
                    <input
                      type="range"
                      min={0}
                      max={frames.length - 1}
                      value={Math.max(0, activeFrameIdx)}
                      onChange={(e) => {
                        setIsReplaying(false);
                        setSelectedMode(null);
                        setActiveFrameIdx(parseInt(e.target.value, 10));
                      }}
                      className="flex-1 apple-slider"
                    />
                    <span className="text-[10.5px] font-mono text-[#86868B] shrink-0">
                      Frame {Math.max(1, activeFrameIdx + 1)}/{frames.length}
                    </span>
                  </div>
                )}
              </div>
            )}

            {/* Live 6-Term AutoDock Vina Score Decomposition */}
            {displayTerms && (
              <div className="space-y-2">
                <div className="text-[12px] font-semibold">
                  Live AutoDock Vina Energy Decomposition (How the Score Is Generated)
                </div>
                <div className="grid grid-cols-1 sm:grid-cols-2 gap-2 text-[11px]">
                  {[
                    { label: 'Gauss1 Steric Fit', val: displayTerms.gauss1, color: '#007AFF', maxAbs: 6.0 },
                    { label: 'Gauss2 Long-Range', val: displayTerms.gauss2, color: '#5856D6', maxAbs: 4.0 },
                    { label: 'Hydrophobic Contacts', val: displayTerms.hydrophobic, color: '#FF9500', maxAbs: 3.5 },
                    { label: 'Hydrogen Bonding', val: displayTerms.hbond, color: '#00C7BE', maxAbs: 4.0 },
                    { label: 'Steric Repulsion (+Penalty)', val: displayTerms.repulsion, color: '#FF3B30', maxAbs: 4.0 },
                    { label: 'Torsional Entropy (+Penalty)', val: displayTerms.torsion, color: '#8E8E93', maxAbs: 3.0 },
                  ].map((item) => {
                    const pct = Math.min(100, (Math.abs(item.val) / item.maxAbs) * 100);
                    return (
                      <div key={item.label} className="p-2 rounded-xl bg-black/[0.025] dark:bg-white/[0.035] space-y-1">
                        <div className="flex justify-between font-mono">
                          <span className="text-[#6E6E73] dark:text-[#A1A1A6] font-sans">{item.label}</span>
                          <span className="font-semibold">
                            {item.val > 0 ? '+' : ''}
                            {item.val.toFixed(2)} kcal
                          </span>
                        </div>
                        <div className="h-1.5 rounded-full bg-black/[0.07] dark:bg-white/[0.08] overflow-hidden">
                          <div
                            className="h-full rounded-full transition-all duration-150"
                            style={{ width: `${pct}%`, backgroundColor: item.color }}
                          />
                        </div>
                      </div>
                    );
                  })}
                </div>
              </div>
            )}

            {/* Live Pocket H-Bonds & Contacts */}
            {displayContacts.length > 0 && (
              <div className="space-y-1.5">
                <div className="text-[12px] font-semibold">Live Active-Site Interactions</div>
                <div className="flex flex-wrap gap-1.5">
                  {displayContacts.map((ct) => (
                    <button
                      key={`${ct.key}-${ct.type}`}
                      onClick={() =>
                        onSelectResidue(
                          { chain: ct.chain, resseq: ct.resseq, icode: '', resname: ct.resname, key: ct.key },
                          true
                        )
                      }
                      className={`px-2.5 py-1 rounded-xl text-[11px] font-mono border transition ${
                        ct.type === 'hbond'
                          ? 'bg-[#00D4FF]/15 border-[#00A3FF]/30 text-[#0071A4] dark:text-[#5AC8FA]'
                          : ct.type === 'hydrophobic'
                          ? 'bg-[#FF9500]/15 border-[#FF9500]/25 text-[#B26200] dark:text-[#FF9F0A]'
                          : 'bg-black/[0.035] dark:bg-white/[0.05] border-transparent'
                      }`}
                    >
                      {ct.residue} ({ct.type === 'hbond' ? 'H-Bond' : ct.type}, {ct.distance_A} Å)
                    </button>
                  ))}
                </div>
              </div>
            )}
          </>
        )}

        {/* Top Ranked Binding Modes & Export Downloads */}
        {modes.length > 0 && (
          <div className="space-y-3 pt-2 border-t border-black/[0.07] dark:border-white/[0.08]">
            <div className="flex flex-wrap items-center justify-between gap-2">
              <div className="text-[13px] font-semibold">
                Top {modes.length} Ranked Binding Modes (Click Any Pose to Inspect in 3D)
              </div>
              {resultId && (
                <div className="flex flex-wrap items-center gap-1.5">
                  <button
                    onClick={() => api.downloadDockedResult(resultId, 'complex_pdb', selectedMode || 1)}
                    className="apple-btn-primary !py-1 !px-2.5 text-[11.5px]"
                  >
                    <Download className="w-3.5 h-3.5" /> Complex .PDB
                  </button>
                  <button
                    onClick={() => api.downloadDockedResult(resultId, 'poses_sdf', selectedMode || 1)}
                    className="apple-btn-secondary !py-1 !px-2.5 text-[11.5px]"
                  >
                    <Download className="w-3.5 h-3.5" /> All Poses .SDF
                  </button>
                  <button
                    onClick={() => api.downloadDockedResult(resultId, 'mode_pdb', selectedMode || 1)}
                    className="apple-btn-secondary !py-1 !px-2 text-[11.5px]"
                  >
                    Mode #{selectedMode || 1} .PDB
                  </button>
                </div>
              )}
            </div>

            {vinaCliResult && vinaCliResult.best_affinity_kcal_mol !== null && (
              <div className="p-2.5 rounded-xl bg-[#34C759]/10 border border-[#34C759]/25 text-[11.5px] font-mono text-[#1B8A38] dark:text-[#30D158]">
                AutoDock Vina CLI Verified Best Affinity: {vinaCliResult.best_affinity_kcal_mol} kcal/mol ({vinaCliResult.binary})
              </div>
            )}

            <div className="grid grid-cols-1 sm:grid-cols-3 gap-2">
              {modes.map((m) => {
                const active = selectedMode === m.mode;
                return (
                  <button
                    key={m.mode}
                    onClick={() => {
                      setIsReplaying(false);
                      setSelectedMode(m.mode);
                    }}
                    className={`text-left p-2.5 rounded-xl border transition ${
                      active
                        ? 'bg-[#007AFF]/10 border-[#007AFF]/45 ring-1 ring-[#007AFF]/20'
                        : 'bg-black/[0.025] dark:bg-white/[0.035] border-transparent hover:border-black/15 dark:hover:border-white/15'
                    }`}
                  >
                    <div className="flex items-center justify-between">
                      <span className="text-[12px] font-semibold">Mode #{m.mode}</span>
                      <span className="font-mono font-bold text-[13px] text-[#34C759]">
                        {m.affinity_kcal_mol.toFixed(2)} kcal
                      </span>
                    </div>
                    <div className="text-[10.5px] font-mono text-[#86868B] mt-1">
                      Ki {m.ki} · {m.hbonds.length} H-B · RMSD {m.rmsd_from_best_A} Å
                    </div>
                  </button>
                );
              })}
            </div>
          </div>
        )}

        {/* Helpful empty state before first run */}
        {frames.length === 0 && modes.length === 0 && !isDocking && (
          <div className="rounded-2xl bg-black/[0.02] dark:bg-white/[0.025] p-5 text-center space-y-1.5">
            <RotateCcw className="w-5 h-5 mx-auto text-[#007AFF] opacity-80" />
            <div className="text-[12.5px] font-semibold">Ready for Live 3D Docking</div>
            <div className="text-[11.5px] text-[#86868B] max-w-md mx-auto">
              Prepare your ligand in Step 2 (try 1-click <strong>Nirmatrelvir</strong> on <code>6LU7</code>!), adjust the glowing blue 3D Grid Box, and press <strong>Start Live Docking</strong> to watch the ligand explore, form H-bonds, and converge in real time.
            </div>
          </div>
        )}
      </div>
    </div>
  );
};
