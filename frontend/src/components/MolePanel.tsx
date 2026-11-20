import React, { useEffect, useMemo, useState } from 'react';
import {
  Dna,
  Download,
  Eye,
  EyeOff,
  Play,
  SlidersHorizontal,
  Sparkles,
  Waypoints,
} from 'lucide-react';
import { api } from '../api/client';
import type {
  MoleAnalysisData,
  MoleChannel,
  MoleColorBy,
  MoleMutationScanResult,
  SelectedResidue,
  StructureSessionState,
} from '../types/pdbprep';

interface MolePanelProps {
  session: StructureSessionState;
  busy: boolean;
  moleData: MoleAnalysisData | null;
  visibleIds: string[];
  colorBy: MoleColorBy;
  onUpdateMole: (data: MoleAnalysisData, visibleIds: string[]) => void;
  onToggleVisible: (id: string) => void;
  onSetAllVisible: (ids: string[]) => void;
  onChangeColorBy: (c: MoleColorBy) => void;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
  onJumpToCootMutate: (res: SelectedResidue) => void;
  onApplyMutation: (chain: string, resseq: number, icode: string, targetAA: string) => Promise<void>;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
}

export const MolePanel: React.FC<MolePanelProps> = ({
  session,
  busy,
  moleData,
  visibleIds,
  colorBy,
  onUpdateMole,
  onToggleVisible,
  onSetAllVisible,
  onChangeColorBy,
  onSelectResidue,
  onJumpToCootMutate,
  onApplyMutation,
  onSetBusy,
  onError,
}) => {
  const [mode, setMode] = useState<'both' | 'tunnels' | 'pores'>('both');
  const [startMode, setStartMode] = useState<'auto' | 'ligand' | 'pocket' | 'residues' | 'coordinates'>('auto');
  const [customResidues, setCustomResidues] = useState('A:41, A:145');
  const [coordX, setCoordX] = useState('0');
  const [coordY, setCoordY] = useState('0');
  const [coordZ, setCoordZ] = useState('0');
  const [bottleneckRadius, setBottleneckRadius] = useState(0.9);
  const [probeRadius, setProbeRadius] = useState(4.2);
  const [minLength, setMinLength] = useState(5.0);
  const [maxTunnels, setMaxTunnels] = useState(6);
  const [includeHetatms, setIncludeHetatms] = useState(false);

  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [hoverIdx, setHoverIdx] = useState<number | null>(null);
  const [scanTargetAA, setScanTargetAA] = useState<'ALA' | 'GLY' | 'VAL'>('ALA');
  const [scanResult, setScanResult] = useState<MoleMutationScanResult | null>(null);

  const buildPayload = () => {
    const start_residues =
      startMode === 'residues'
        ? customResidues
            .split(/[,\s]+/)
            .map((s) => s.trim())
            .filter(Boolean)
        : null;
    const start_point =
      startMode === 'coordinates'
        ? [parseFloat(coordX) || 0, parseFloat(coordY) || 0, parseFloat(coordZ) || 0]
        : null;
    return {
      mode,
      bottleneck_radius: bottleneckRadius,
      probe_radius: probeRadius,
      min_length: minLength,
      max_tunnels: maxTunnels,
      start_mode: startMode,
      start_residues,
      start_point,
      include_hetatms: includeHetatms,
    };
  };

  const runAnalysis = async () => {
    onSetBusy('Computing 3D Voronoi Tunnels & Pores (MOLE)…');
    try {
      const data = await api.analyzeMole(session.session_id, buildPayload());
      const all = [...data.tunnels, ...data.pores];
      // Show top 3 tunnels + first pore cleanly in 3D by default
      const defaultVis = [
        ...data.tunnels.slice(0, 3).map((t) => t.id),
        ...data.pores.slice(0, 1).map((p) => p.id),
      ];
      onUpdateMole(data, defaultVis);
      setSelectedId(all[0]?.id || null);
      setScanResult(null);
    } catch (e: any) {
      onError(e?.message || 'MOLE tunnel calculation failed');
    } finally {
      onSetBusy(null);
    }
  };

  // Auto-run MOLE once when opening the tab for the current structure if not computed yet
  useEffect(() => {
    if (!moleData && !busy) {
      void runAnalysis();
    }
  }, [session.session_id]);

  const allChannels: MoleChannel[] = useMemo(() => {
    if (!moleData) return [];
    return [...moleData.tunnels, ...moleData.pores];
  }, [moleData]);

  const activeChannel = useMemo(() => {
    if (!allChannels.length) return null;
    return allChannels.find((c) => c.id === selectedId) || allChannels[0];
  }, [allChannels, selectedId]);

  const handleScanMutations = async () => {
    if (!activeChannel) return;
    const bnKeys = activeChannel.bottleneck_residues.length
      ? activeChannel.bottleneck_residues
      : activeChannel.lining_residues.slice(0, 6).map((r) => r.key);
    if (!bnKeys.length) return;
    onSetBusy(`Scanning ${bnKeys.length} bottleneck residues → ${scanTargetAA}…`);
    try {
      const res = await api.scanMoleMutations(session.session_id, {
        bottleneck_residues: bnKeys,
        bottleneck_point: activeChannel.bottleneck_point,
        baseline_bottleneck_A: activeChannel.bottleneck_radius_A,
        target_aa: scanTargetAA,
        include_hetatms: includeHetatms,
      });
      setScanResult(res);
    } catch (e: any) {
      onError(e?.message || 'Bottleneck mutation scan failed');
    } finally {
      onSetBusy(null);
    }
  };

  const handleExport = async (fmt: 'csv' | 'residues_csv' | 'json' | 'pdb') => {
    onSetBusy(`Exporting MOLE tunnels (${fmt.toUpperCase()})…`);
    try {
      await api.exportMole(
        session.session_id,
        { ...buildPayload(), format: fmt },
        fmt === 'residues_csv' ? 'mole_residues.csv' : `mole_tunnels.${fmt}`
      );
    } catch (e: any) {
      onError(e?.message || 'Export failed');
    } finally {
      onSetBusy(null);
    }
  };

  // Build SVG MOLE Radius + Hydropathy Profile Plot
  const svgW = 520;
  const svgH = 165;
  const padL = 38;
  const padR = 16;
  const padT = 16;
  const padB = 34;

  const profile = activeChannel?.profile || [];
  const maxDist = Math.max(5, activeChannel?.length_A || 10);
  const maxRad = Math.max(3.5, ...(profile.map((p) => p.radius)));
  const px = (d: number) => padL + (d / maxDist) * (svgW - padL - padR);
  const py = (r: number) => padT + (1 - r / maxRad) * (svgH - padT - padB);

  const linePath = profile
    .map((p, i) => `${i === 0 ? 'M' : 'L'} ${px(p.distance).toFixed(1)} ${py(p.radius).toFixed(1)}`)
    .join(' ');
  const areaPath =
    profile.length > 1
      ? `${linePath} L ${px(profile[profile.length - 1].distance).toFixed(1)} ${py(0).toFixed(1)} L ${px(
          profile[0].distance
        ).toFixed(1)} ${py(0).toFixed(1)} Z`
      : '';

  const hydroColor = (h: number) => {
    if (h >= 1.5) return '#FF9500';
    if (h >= 0.4) return '#FFCC00';
    if (h <= -1.5) return '#007AFF';
    if (h <= -0.4) return '#5AC8FA';
    return '#98989D';
  };

  const hoveredPt = hoverIdx !== null && profile[hoverIdx] ? profile[hoverIdx] : null;

  return (
    <div className="space-y-4 pb-6">
      {/* Header */}
      <div className="apple-card p-4 flex flex-wrap items-center justify-between gap-3">
        <div>
          <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#007AFF] dark:text-[#0A84FF]">
            MOLE 2.5 Voronoi Void-Space Analysis
          </div>
          <h2 className="text-[16px] font-semibold tracking-[-0.015em]">
            Substrate Tunnels, Transmembrane Pores &amp; Bottleneck Engineering
          </h2>
        </div>
        <div className="flex flex-wrap items-center gap-1.5">
          {moleData && allChannels.length > 0 && (
            <>
              <button
                disabled={busy}
                onClick={() => handleExport('csv')}
                className="apple-btn-secondary text-[12px]"
                title="Download profile CSV with per-slice & bottleneck residues"
              >
                <Download className="w-3.5 h-3.5" /> Profile CSV
              </button>
              <button
                disabled={busy}
                onClick={() => handleExport('residues_csv')}
                className="apple-btn-secondary text-[12px]"
                title="Download residue-by-residue lining & bottleneck CSV table"
              >
                <Download className="w-3.5 h-3.5" /> Residues CSV
              </button>
              <button
                disabled={busy}
                onClick={() => handleExport('pdb')}
                className="apple-btn-secondary text-[12px]"
                title="Download PyMOL / ChimeraX pseudo-atom PDB with radii in B-factor"
              >
                <Download className="w-3.5 h-3.5" /> PyMOL PDB
              </button>
              <button disabled={busy} onClick={() => handleExport('json')} className="apple-btn-secondary text-[12px]">
                JSON
              </button>
            </>
          )}
          <button disabled={busy} onClick={runAnalysis} className="apple-btn-primary">
            <Play className="w-3.5 h-3.5 fill-current" />
            <span>Compute Tunnels</span>
          </button>
        </div>
      </div>

      {/* Controls Card */}
      <div className="apple-card p-4 space-y-3.5">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div className="flex items-center gap-2 text-[13px] font-semibold">
            <SlidersHorizontal className="w-4 h-4 text-[#007AFF]" />
            <span>Voronoi Channel Parameters</span>
          </div>

          <div className="apple-segmented">
            {(
              [
                ['both', 'Tunnels & Pores'],
                ['tunnels', 'Tunnels'],
                ['pores', 'Pores'],
              ] as const
            ).map(([k, lbl]) => (
              <button
                key={k}
                onClick={() => setMode(k)}
                className={`apple-segmented-item ${mode === k ? 'apple-segmented-active' : 'text-[#6E6E73]'}`}
              >
                {lbl}
              </button>
            ))}
          </div>
        </div>

        {/* Origin Picker + 3D Coloring */}
        <div className="grid grid-cols-1 sm:grid-cols-12 gap-2.5 text-[12px]">
          <div className="sm:col-span-5">
            <label className="text-[10.5px] text-[#86868B] block mb-1">Tunnel Starting Origin</label>
            <select
              value={startMode}
              onChange={(e) => setStartMode(e.target.value as any)}
              className="apple-select w-full"
            >
              <option value="auto">Auto (Active-Site Ligand / Cavity)</option>
              <option value="ligand">Co-Crystallized Ligands Only</option>
              <option value="pocket">Deepest Detected Cavities</option>
              <option value="residues">Custom Catalytic Residues</option>
              <option value="coordinates">Custom 3D Coordinates (X, Y, Z)</option>
            </select>
          </div>

          <div className="sm:col-span-4">
            <label className="text-[10.5px] text-[#86868B] block mb-1">3D Tunnel Color Scheme</label>
            <select
              value={colorBy}
              onChange={(e) => onChangeColorBy(e.target.value as MoleColorBy)}
              className="apple-select w-full"
            >
              <option value="id">Distinct Tunnel ID</option>
              <option value="hydropathy">Hydropathy (Blue=Polar · Orange=Hydrophobic)</option>
              <option value="polarity">Polarity (Zimmerman)</option>
              <option value="charge">Formal Charge (Red=- · Blue=+)</option>
            </select>
          </div>

          <div className="sm:col-span-3 flex items-end">
            <label className="w-full flex items-center justify-between px-2.5 py-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04] cursor-pointer text-[11.5px]">
              <span>Include Ligands</span>
              <input
                type="checkbox"
                checked={includeHetatms}
                onChange={(e) => setIncludeHetatms(e.target.checked)}
                className="accent-[#007AFF] w-3.5 h-3.5"
              />
            </label>
          </div>
        </div>

        {/* Conditional custom origin inputs */}
        {startMode === 'residues' && (
          <div>
            <label className="text-[10.5px] text-[#86868B] block mb-1">
              Starting Residue(s) — comma-separated (e.g. A:41, A:145 for 6LU7 catalytic dyad)
            </label>
            <input
              value={customResidues}
              onChange={(e) => setCustomResidues(e.target.value)}
              placeholder="A:41, A:145"
              className="apple-input font-mono"
            />
          </div>
        )}
        {startMode === 'coordinates' && (
          <div className="grid grid-cols-3 gap-2">
            <input value={coordX} onChange={(e) => setCoordX(e.target.value)} placeholder="X (Å)" className="apple-input font-mono" />
            <input value={coordY} onChange={(e) => setCoordY(e.target.value)} placeholder="Y (Å)" className="apple-input font-mono" />
            <input value={coordZ} onChange={(e) => setCoordZ(e.target.value)} placeholder="Z (Å)" className="apple-input font-mono" />
          </div>
        )}

        {/* Sliders */}
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 text-[11.5px]">
          <div>
            <div className="flex justify-between mb-1">
              <span className="text-[#86868B]">Bottleneck Min</span>
              <span className="font-mono font-semibold">{bottleneckRadius.toFixed(1)} Å</span>
            </div>
            <input
              type="range"
              min={0.5}
              max={2.5}
              step={0.1}
              value={bottleneckRadius}
              onChange={(e) => setBottleneckRadius(parseFloat(e.target.value))}
              className="w-full apple-slider"
            />
          </div>
          <div>
            <div className="flex justify-between mb-1">
              <span className="text-[#86868B]">Surface Probe</span>
              <span className="font-mono font-semibold">{probeRadius.toFixed(1)} Å</span>
            </div>
            <input
              type="range"
              min={2.5}
              max={6.5}
              step={0.1}
              value={probeRadius}
              onChange={(e) => setProbeRadius(parseFloat(e.target.value))}
              className="w-full apple-slider"
            />
          </div>
          <div>
            <div className="flex justify-between mb-1">
              <span className="text-[#86868B]">Min Length</span>
              <span className="font-mono font-semibold">{minLength.toFixed(0)} Å</span>
            </div>
            <input
              type="range"
              min={3}
              max={20}
              step={1}
              value={minLength}
              onChange={(e) => setMinLength(parseFloat(e.target.value))}
              className="w-full apple-slider"
            />
          </div>
          <div>
            <div className="flex justify-between mb-1">
              <span className="text-[#86868B]">Max Tunnels</span>
              <span className="font-mono font-semibold">{maxTunnels}</span>
            </div>
            <input
              type="range"
              min={1}
              max={10}
              step={1}
              value={maxTunnels}
              onChange={(e) => setMaxTunnels(parseInt(e.target.value, 10))}
              className="w-full apple-slider"
            />
          </div>
        </div>
      </div>

      {/* Results KPIs + Channel List */}
      {moleData && (
        <>
          <div className="grid grid-cols-2 sm:grid-cols-4 gap-2.5">
            <div className="apple-card p-3">
              <div className="text-[11px] text-[#86868B]">Tunnels Found</div>
              <div className="text-[18px] font-bold font-mono text-[#007AFF] dark:text-[#0A84FF]">
                {moleData.summary.total_tunnels}
              </div>
            </div>
            <div className="apple-card p-3">
              <div className="text-[11px] text-[#86868B]">Trans-Pores</div>
              <div className="text-[18px] font-bold font-mono text-[#AF52DE] dark:text-[#BF5AF2]">
                {moleData.summary.total_pores}
              </div>
            </div>
            <div className="apple-card p-3">
              <div className="text-[11px] text-[#86868B]">Widest Bottleneck</div>
              <div className="text-[18px] font-bold font-mono text-[#34C759]">
                {moleData.summary.widest_bottleneck_A || 0} Å
              </div>
            </div>
            <div className="apple-card p-3">
              <div className="text-[11px] text-[#86868B]">Voronoi Vertices</div>
              <div className="text-[18px] font-bold font-mono">
                {moleData.summary.voronoi_nodes.toLocaleString()}
              </div>
            </div>
          </div>

          {allChannels.length === 0 ? (
            <div className="apple-card p-5 text-[13px] text-[#86868B]">
              No tunnels or pores met the minimum bottleneck ({bottleneckRadius} Å) and length ({minLength} Å) criteria. Try lowering Bottleneck Min to 0.7 Å or selecting custom catalytic residues.
            </div>
          ) : (
            <>
              {/* Tunnels & Pores List */}
              <div className="apple-card p-4 space-y-2.5">
                <div className="flex flex-wrap items-center justify-between gap-2">
                  <div className="flex items-center gap-2 text-[13.5px] font-semibold">
                    <Waypoints className="w-4 h-4 text-[#007AFF]" />
                    <span>Detected Channels ({allChannels.length}) — Click to Inspect Profile</span>
                  </div>
                  <div className="flex items-center gap-1.5 text-[11.5px]">
                    <button
                      onClick={() => onSetAllVisible(allChannels.map((c) => c.id))}
                      className="apple-btn-secondary !py-1 !px-2 text-[11px]"
                    >
                      Show All in 3D
                    </button>
                    <button
                      onClick={() => onSetAllVisible([])}
                      className="apple-btn-secondary !py-1 !px-2 text-[11px]"
                    >
                      Hide All
                    </button>
                  </div>
                </div>

                <div className="space-y-1.5 max-h-60 overflow-y-auto">
                  {allChannels.map((ch) => {
                    const isSel = activeChannel?.id === ch.id;
                    const isVis = visibleIds.includes(ch.id);
                    return (
                      <div
                        key={ch.id}
                        onClick={() => setSelectedId(ch.id)}
                        className={`cursor-pointer flex flex-wrap items-center justify-between gap-2 p-2.5 rounded-xl border transition ${
                          isSel
                            ? 'bg-[#007AFF]/[0.08] border-[#007AFF]/40 ring-1 ring-[#007AFF]/20'
                            : 'bg-black/[0.025] dark:bg-white/[0.035] border-transparent hover:border-black/10 dark:hover:border-white/10'
                        }`}
                      >
                        <div className="flex items-center gap-2.5 min-w-0">
                          <span
                            className="w-3 h-3 rounded-full shrink-0 shadow-sm"
                            style={{ backgroundColor: ch.color }}
                          />
                          <div className="min-w-0">
                            <div className="flex items-center gap-1.5 text-[12.5px] font-semibold">
                              <span>{ch.name}</span>
                              <span
                                className={`apple-pill uppercase ${
                                  ch.kind === 'pore'
                                    ? 'bg-[#AF52DE]/15 text-[#8944AB] dark:text-[#BF5AF2]'
                                    : 'bg-[#007AFF]/15 text-[#007AFF] dark:text-[#0A84FF]'
                                }`}
                              >
                                {ch.kind}
                              </span>
                              <span className="apple-pill bg-black/[0.05] dark:bg-white/[0.08] text-[#6E6E73] dark:text-[#A1A1A6]">
                                {ch.properties.character}
                              </span>
                            </div>
                            <div className="text-[11px] font-mono text-[#86868B] mt-0.5 truncate">
                              Len {ch.length_A} Å · Bottleneck{' '}
                              <span className="text-[#FF3B30] font-semibold">{ch.bottleneck_radius_A} Å</span> · Throughput{' '}
                              {ch.throughput} · Charge {ch.properties.net_charge >= 0 ? '+' : ''}
                              {ch.properties.net_charge} · {ch.origin_label}
                            </div>
                          </div>
                        </div>

                        <button
                          onClick={(e) => {
                            e.stopPropagation();
                            onToggleVisible(ch.id);
                          }}
                          title={isVis ? 'Hide in 3D' : 'Show in 3D'}
                          className={`p-1.5 rounded-lg transition ${
                            isVis
                              ? 'bg-[#007AFF] text-white'
                              : 'bg-black/[0.05] dark:bg-white/[0.08] text-[#86868B]'
                          }`}
                        >
                          {isVis ? <Eye className="w-3.5 h-3.5" /> : <EyeOff className="w-3.5 h-3.5" />}
                        </button>
                      </div>
                    );
                  })}
                </div>
              </div>

              {/* Active Channel MOLE 2D Profile Chart + Physicochemical Properties */}
              {activeChannel && (
                <div className="apple-card p-4 space-y-3.5">
                  <div className="flex flex-wrap items-center justify-between gap-2">
                    <div>
                      <div className="text-[13.5px] font-semibold flex items-center gap-2">
                        <span
                          className="w-2.5 h-2.5 rounded-full"
                          style={{ backgroundColor: activeChannel.color }}
                        />
                        <span>{activeChannel.name} — Radius &amp; Hydropathy Profile</span>
                      </div>
                      <div className="text-[11px] font-mono text-[#86868B] mt-0.5">
                        Bottleneck {activeChannel.bottleneck_radius_A} Å at {activeChannel.bottleneck_distance_A} Å · Max{' '}
                        {activeChannel.max_radius_A} Å · Tortuosity {activeChannel.tortuosity} · Mutability{' '}
                        {activeChannel.properties.mutability}
                      </div>
                    </div>
                    {hoveredPt && (
                      <div className="px-2.5 py-1 rounded-xl bg-black/[0.04] dark:bg-white/[0.07] text-[11px] font-mono">
                        d={hoveredPt.distance} Å · r={hoveredPt.radius} Å · H={hoveredPt.hydropathy} · Q=
                        {hoveredPt.charge >= 0 ? '+' : ''}
                        {hoveredPt.charge}
                      </div>
                    )}
                  </div>

                  {/* Interactive SVG Radius + Hydropathy Bar */}
                  <div className="rounded-2xl bg-black/[0.025] dark:bg-white/[0.03] p-2">
                    <svg
                      viewBox={`0 0 ${svgW} ${svgH}`}
                      className="w-full h-40 overflow-visible"
                      onMouseLeave={() => setHoverIdx(null)}
                    >
                      {/* Horizontal grid lines */}
                      {[0, 1, 2, 3].map((rVal) =>
                        rVal <= maxRad ? (
                          <g key={rVal}>
                            <line
                              x1={padL}
                              y1={py(rVal)}
                              x2={svgW - padR}
                              y2={py(rVal)}
                              stroke="currentColor"
                              className="text-black/10 dark:text-white/10"
                              strokeDasharray="3 3"
                            />
                            <text
                              x={padL - 6}
                              y={py(rVal) + 3}
                              textAnchor="end"
                              fontSize="9"
                              className="fill-[#86868B]"
                              fontFamily="monospace"
                            >
                              {rVal}Å
                            </text>
                          </g>
                        ) : null
                      )}

                      {/* Shaded radius area */}
                      {areaPath && <path d={areaPath} fill={activeChannel.color} opacity="0.16" />}
                      {linePath && (
                        <path
                          d={linePath}
                          fill="none"
                          stroke={activeChannel.color}
                          strokeWidth="2.4"
                          strokeLinecap="round"
                        />
                      )}

                      {/* Bottleneck vertical marker */}
                      <line
                        x1={px(activeChannel.bottleneck_distance_A)}
                        y1={padT}
                        x2={px(activeChannel.bottleneck_distance_A)}
                        y2={svgH - padB}
                        stroke="#FF3B30"
                        strokeWidth="1.5"
                        strokeDasharray="4 3"
                      />
                      <circle
                        cx={px(activeChannel.bottleneck_distance_A)}
                        cy={py(activeChannel.bottleneck_radius_A)}
                        r="4.5"
                        fill="#FF3B30"
                        stroke="#FFFFFF"
                        strokeWidth="1.5"
                      />
                      <text
                        x={Math.min(svgW - 95, px(activeChannel.bottleneck_distance_A) + 6)}
                        y={padT + 10}
                        fontSize="9.5"
                        fill="#FF3B30"
                        fontWeight="600"
                        fontFamily="monospace"
                      >
                        Bottleneck {activeChannel.bottleneck_radius_A} Å
                      </text>

                      {/* Hydropathy color ribbon below X-axis */}
                      {profile.map((pt, i) => {
                        const segW = (svgW - padL - padR) / Math.max(profile.length, 1);
                        return (
                          <rect
                            key={i}
                            x={px(pt.distance) - segW / 2}
                            y={svgH - padB + 6}
                            width={Math.max(2, segW + 0.5)}
                            height={7}
                            rx={1}
                            fill={hydroColor(pt.hydropathy)}
                            opacity="0.88"
                            onMouseEnter={() => setHoverIdx(i)}
                            className="cursor-crosshair"
                          />
                        );
                      })}

                      {/* Interactive hover hit targets */}
                      {profile.map((pt, i) => (
                        <circle
                          key={`hit-${i}`}
                          cx={px(pt.distance)}
                          cy={py(pt.radius)}
                          r={hoverIdx === i ? 4 : 7}
                          fill={hoverIdx === i ? activeChannel.color : 'transparent'}
                          className="cursor-crosshair"
                          onMouseEnter={() => setHoverIdx(i)}
                        />
                      ))}

                      <text
                        x={svgW / 2}
                        y={svgH - 4}
                        textAnchor="middle"
                        fontSize="9"
                        className="fill-[#86868B]"
                        fontFamily="monospace"
                      >
                        Distance along {activeChannel.kind} (0 → {activeChannel.length_A} Å) · Hydropathy bar (Blue=Polar · Orange=Hydrophobic)
                      </text>
                    </svg>
                  </div>

                  {/* Physicochemical Summary Strip */}
                  <div className="grid grid-cols-2 sm:grid-cols-5 gap-2 text-center">
                    <div className="p-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
                      <div className="text-[10px] text-[#86868B]">Hydropathy</div>
                      <div className="text-[13px] font-bold font-mono">{activeChannel.properties.hydropathy}</div>
                    </div>
                    <div className="p-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
                      <div className="text-[10px] text-[#86868B]">Polarity</div>
                      <div className="text-[13px] font-bold font-mono">{activeChannel.properties.polarity}</div>
                    </div>
                    <div className="p-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
                      <div className="text-[10px] text-[#86868B]">Net Charge</div>
                      <div className="text-[13px] font-bold font-mono">
                        {activeChannel.properties.net_charge >= 0 ? '+' : ''}
                        {activeChannel.properties.net_charge} ({activeChannel.properties.positive_residues}+/
                        {activeChannel.properties.negative_residues}-)
                      </div>
                    </div>
                    <div className="p-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
                      <div className="text-[10px] text-[#86868B]">Mutability</div>
                      <div className="text-[13px] font-bold font-mono">{activeChannel.properties.mutability}</div>
                    </div>
                    <div className="p-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
                      <div className="text-[10px] text-[#86868B]">Lining Res</div>
                      <div className="text-[13px] font-bold font-mono">{activeChannel.lining_residues.length}</div>
                    </div>
                  </div>

                  {/* Automated Bottleneck Mutation Scanner */}
                  <div className="p-3 rounded-2xl bg-[#AF52DE]/[0.06] border border-[#AF52DE]/20 space-y-2.5">
                    <div className="flex flex-wrap items-center justify-between gap-2">
                      <div>
                        <div className="text-[12.5px] font-semibold flex items-center gap-1.5 text-[#8944AB] dark:text-[#BF5AF2]">
                          <Sparkles className="w-3.5 h-3.5" />
                          <span>In-Silico Bottleneck Mutation Scanner</span>
                        </div>
                        <div className="text-[11px] text-[#86868B]">
                          Tests single-point mutations at bottleneck residues to predict tunnel widening (Δr in Å).
                        </div>
                      </div>
                      <div className="flex items-center gap-1.5">
                        <select
                          value={scanTargetAA}
                          onChange={(e) => setScanTargetAA(e.target.value as 'ALA' | 'GLY' | 'VAL')}
                          className="apple-select text-[11.5px] !py-1"
                        >
                          <option value="ALA">Mutate → ALA</option>
                          <option value="GLY">Mutate → GLY</option>
                          <option value="VAL">Mutate → VAL</option>
                        </select>
                        <button
                          disabled={busy}
                          onClick={handleScanMutations}
                          className="apple-btn-primary !py-1 !px-2.5 text-[11.5px] !bg-[#AF52DE] hover:!bg-[#8944AB]"
                        >
                          <Dna className="w-3.5 h-3.5" />
                          Scan Bottleneck
                        </button>
                      </div>
                    </div>

                    {scanResult && (
                      <div className="space-y-1.5 pt-1">
                        {scanResult.candidates.length === 0 ? (
                          <div className="text-[11.5px] text-[#86868B]">
                            No mutatable non-target residues found at this bottleneck.
                          </div>
                        ) : (
                          scanResult.candidates.map((cand) => (
                            <div
                              key={cand.mutation}
                              className="flex flex-wrap items-center justify-between gap-2 p-2 rounded-xl bg-white/80 dark:bg-black/25 text-[11.5px]"
                            >
                              <div className="flex items-center gap-2 min-w-0">
                                <span className="font-mono font-bold text-[#1D1D1F] dark:text-white">
                                  {cand.mutation}
                                </span>
                                <span
                                  className={`apple-pill font-mono ${
                                    cand.delta_bottleneck_A > 0.05
                                      ? 'bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158]'
                                      : 'bg-black/[0.06] dark:bg-white/[0.08] text-[#86868B]'
                                  }`}
                                >
                                  {cand.delta_bottleneck_A >= 0 ? '+' : ''}
                                  {cand.delta_bottleneck_A} Å
                                </span>
                                <span className="font-mono text-[#86868B]">
                                  {cand.baseline_bottleneck_A} → {cand.new_bottleneck_A} Å ({cand.effect})
                                </span>
                              </div>
                              <button
                                disabled={busy}
                                onClick={async () => {
                                  await onApplyMutation(cand.chain, cand.resseq, cand.icode, cand.mut_resname);
                                  setScanResult(null);
                                  await runAnalysis();
                                }}
                                className="apple-btn-primary !py-0.5 !px-2 text-[11px]"
                              >
                                Apply Mutation
                              </button>
                            </div>
                          ))
                        )}
                      </div>
                    )}
                  </div>

                  {/* Lining & Bottleneck Residues Table (Click to Focus 3D or Mutate in COOT!) */}
                  <div className="space-y-1.5">
                    <div className="text-[12.5px] font-semibold">
                      Tunnel-Lining &amp; Bottleneck Residues (Click to Inspect or Mutate in COOT)
                    </div>
                    <div className="max-h-56 overflow-y-auto divide-y divide-black/[0.05] dark:divide-white/[0.06] text-[12px]">
                      {activeChannel.lining_residues.map((lr) => (
                        <div
                          key={lr.key}
                          className="flex items-center justify-between gap-2 py-1.5 px-2 hover:bg-black/[0.025] dark:hover:bg-white/[0.03] rounded-lg transition"
                        >
                          <button
                            onClick={() =>
                              onSelectResidue(
                                {
                                  chain: lr.chain,
                                  resseq: lr.resseq,
                                  icode: lr.icode,
                                  resname: lr.resname,
                                  key: lr.key,
                                },
                                true
                              )
                            }
                            className="flex items-center gap-2 text-left min-w-0"
                          >
                            <span className="font-mono font-semibold">
                              {lr.resname} {lr.key}
                            </span>
                            {lr.is_bottleneck && (
                              <span className="apple-pill bg-[#FF3B30]/15 text-[#D70015] dark:text-[#FF453A]">
                                Bottleneck
                              </span>
                            )}
                            <span className="text-[11px] font-mono text-[#86868B]">
                              {lr.start_dist}–{lr.end_dist} Å · H={lr.hydropathy}
                              {lr.backbone_contact ? ' · BB' : ' · SC'}
                            </span>
                          </button>

                          <button
                            onClick={() =>
                              onJumpToCootMutate({
                                chain: lr.chain,
                                resseq: lr.resseq,
                                icode: lr.icode,
                                resname: lr.resname,
                                key: lr.key,
                              })
                            }
                            title="Open this bottleneck/lining residue in COOT Workbench to mutate & redesign tunnel"
                            className="apple-btn-secondary !py-0.5 !px-2 text-[11px] shrink-0"
                          >
                            <Dna className="w-3 h-3 text-[#AF52DE]" />
                            COOT Mutate
                          </button>
                        </div>
                      ))}
                    </div>
                  </div>
                </div>
              )}
            </>
          )}
        </>
      )}
    </div>
  );
};
