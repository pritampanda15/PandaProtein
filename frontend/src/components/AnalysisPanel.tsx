import React, { useEffect, useRef, useState } from 'react';
import {
  Compass,
  Eye,
  EyeOff,
  GitCompare,
  Layers,
  Magnet,
  RefreshCw,
  Upload,
} from 'lucide-react';
import { api } from '../api/client';
import type {
  InteractionsMembraneData,
  MembraneInfo,
  SelectedResidue,
  StructuralDiffData,
  StructureSessionState,
} from '../types/pdbprep';

interface AnalysisPanelProps {
  session: StructureSessionState;
  busy: boolean;
  diffOverlayActive: boolean;
  membraneOverlayActive: boolean;
  onToggleDiffOverlay: (pdbText: string | null) => void;
  onToggleMembraneOverlay: (mem: MembraneInfo | null) => void;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
  onStateUpdate: (state: StructureSessionState, toast: string) => void;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
}

export const AnalysisPanel: React.FC<AnalysisPanelProps> = ({
  session,
  busy,
  diffOverlayActive,
  membraneOverlayActive,
  onToggleDiffOverlay,
  onToggleMembraneOverlay,
  onSelectResidue,
  onStateUpdate,
  onSetBusy,
  onError,
}) => {
  const [diffData, setDiffData] = useState<StructuralDiffData | null>(null);
  const [interData, setInterData] = useState<InteractionsMembraneData | null>(null);
  const [netTab, setNetTab] = useState<'salt' | 'pi' | 'interface'>('salt');
  const fileRef = useRef<HTMLInputElement | null>(null);

  const loadAll = async () => {
    try {
      const [d, inter] = await Promise.all([
        api.getStructuralDiff(session.session_id),
        api.getInteractionsMembrane(session.session_id),
      ]);
      setDiffData(d);
      setInterData(inter);
    } catch (e: any) {
      onError(e?.message || 'Structural analysis failed');
    }
  };

  useEffect(() => {
    void loadAll();
  }, [session.session_id, session.pdb_text]);

  const handleAlignUpload = async (files: FileList | null) => {
    if (!files || !files.length) return;
    onSetBusy(`Aligning ${files[0].name} via Kabsch Cα superposition…`);
    try {
      const d = await api.alignReferenceFile(session.session_id, files[0]);
      setDiffData(d);
      onToggleDiffOverlay(d.aligned_ref_pdb);
    } catch (e: any) {
      onError(e?.message || 'Superposition failed');
    } finally {
      onSetBusy(null);
    }
  };

  const handleOrientZ = async () => {
    onSetBusy('Aligning membrane / pore axis along +Z…');
    try {
      const r = await api.orientAlongZ(session.session_id);
      onStateUpdate(r.state, r.result.message);
    } catch (e: any) {
      onError(e?.message || 'Z-axis orientation failed');
    } finally {
      onSetBusy(null);
    }
  };

  // Per-residue shift SVG bar chart
  const shifts = diffData?.per_residue || [];
  const maxShift = Math.max(0.5, ...(shifts.map((s) => s.max_atom_shift_A)));

  return (
    <div className="space-y-4 pb-6">
      {/* Header */}
      <div className="apple-card p-4 flex flex-wrap items-center justify-between gap-3">
        <div>
          <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#007AFF] dark:text-[#0A84FF]">
            Superposition · Bilayer Orientation · Non-Covalent Networks
          </div>
          <h2 className="text-[16px] font-semibold tracking-[-0.015em]">
            Before/After RMSD Diff, Transmembrane Belt &amp; Salt-Bridge Analyzer
          </h2>
        </div>
        <button disabled={busy} onClick={loadAll} className="apple-btn-secondary text-[12px]">
          <RefreshCw className="w-3.5 h-3.5" /> Refresh
        </button>
      </div>

      {/* 1. Kabsch 3D Structural Superposition & RMSD Diff */}
      <div className="apple-card p-4 space-y-3.5">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div className="flex items-center gap-2 text-[13.5px] font-semibold">
            <GitCompare className="w-4 h-4 text-[#007AFF]" />
            <span>3D Structural Superposition &amp; RMSD ({diffData?.ref_label || 'Original'})</span>
          </div>
          <div className="flex items-center gap-1.5">
            <input
              ref={fileRef}
              type="file"
              accept=".pdb,.ent,.cif,.mmcif,.gro"
              className="hidden"
              onChange={(e) => handleAlignUpload(e.target.files)}
            />
            <button
              disabled={busy}
              onClick={() => fileRef.current?.click()}
              className="apple-btn-secondary text-[12px]"
              title="Upload a second structure (e.g. AlphaFold or crystal) to align in 3D"
            >
              <Upload className="w-3.5 h-3.5" /> Compare 2nd PDB
            </button>
            {diffData && (
              <button
                onClick={() => onToggleDiffOverlay(diffOverlayActive ? null : diffData.aligned_ref_pdb)}
                className={diffOverlayActive ? 'apple-btn-primary text-[12px]' : 'apple-btn-secondary text-[12px]'}
              >
                {diffOverlayActive ? <EyeOff className="w-3.5 h-3.5" /> : <Eye className="w-3.5 h-3.5" />}
                <span>{diffOverlayActive ? 'Hide 3D Ghost' : 'Overlay in 3D'}</span>
              </button>
            )}
          </div>
        </div>

        {diffData && (
          <>
            <div className="grid grid-cols-2 sm:grid-cols-4 gap-2 text-center">
              <div className="p-2.5 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
                <div className="text-[10px] text-[#86868B]">Cα RMSD</div>
                <div className="text-[16px] font-bold font-mono text-[#007AFF] dark:text-[#0A84FF]">
                  {diffData.ca_rmsd_A} Å
                </div>
              </div>
              <div className="p-2.5 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
                <div className="text-[10px] text-[#86868B]">Backbone RMSD</div>
                <div className="text-[16px] font-bold font-mono">{diffData.backbone_rmsd_A} Å</div>
              </div>
              <div className="p-2.5 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
                <div className="text-[10px] text-[#86868B]">All-Heavy RMSD</div>
                <div className="text-[16px] font-bold font-mono">{diffData.all_heavy_rmsd_A} Å</div>
              </div>
              <div className="p-2.5 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
                <div className="text-[10px] text-[#86868B]">Max Atom Shift</div>
                <div className="text-[16px] font-bold font-mono text-[#FF9500]">{diffData.max_shift_A} Å</div>
              </div>
            </div>

            {/* Per-residue displacement sparkline */}
            {shifts.length > 0 && (
              <div className="space-y-1.5">
                <div className="text-[11.5px] font-medium text-[#86868B]">
                  Per-Residue Max Heavy-Atom Displacement (matched {diffData.matched_ca} Cα)
                </div>
                <svg viewBox="0 0 500 84" className="w-full h-20 rounded-xl bg-black/[0.025] dark:bg-white/[0.03]">
                  {shifts.map((s, idx) => {
                    const x = 14 + (idx / Math.max(shifts.length - 1, 1)) * 472;
                    const h = (s.max_atom_shift_A / maxShift) * 62;
                    return (
                      <line
                        key={s.key}
                        x1={x}
                        y1={74}
                        x2={x}
                        y2={74 - h}
                        stroke={s.max_atom_shift_A > 1.0 ? '#FF9500' : '#007AFF'}
                        strokeWidth={Math.max(1.2, 460 / Math.max(shifts.length, 40))}
                      >
                        <title>{`${s.resname} ${s.key}: max shift ${s.max_atom_shift_A} Å (Cα ${s.ca_shift_A} Å)`}</title>
                      </line>
                    );
                  })}
                </svg>
              </div>
            )}

            {/* Most-shifted residues clickable list */}
            {diffData.top_shifted.some((r) => r.max_atom_shift_A > 0.05) && (
              <div className="flex flex-wrap gap-1.5">
                {diffData.top_shifted
                  .filter((r) => r.max_atom_shift_A > 0.05)
                  .slice(0, 10)
                  .map((r) => (
                    <button
                      key={r.key}
                      onClick={() =>
                        onSelectResidue({ chain: r.chain, resseq: r.resseq, icode: '', resname: r.resname, key: r.key }, true)
                      }
                      className="px-2 py-1 rounded-lg text-[11px] font-mono bg-[#FF9500]/10 text-[#B26200] dark:text-[#FF9F0A] hover:bg-[#FF9500]/20 transition"
                    >
                      {r.resname} {r.key} ({r.max_atom_shift_A} Å)
                    </button>
                  ))}
              </div>
            )}
          </>
        )}
      </div>

      {/* 2. Transmembrane Hydrophobic Belt & Z-Axis Pore Orientation */}
      {interData?.membrane && (
        <div className="apple-card p-4 space-y-3">
          <div className="flex flex-wrap items-center justify-between gap-2">
            <div>
              <div className="flex items-center gap-2 text-[13.5px] font-semibold">
                <Compass className="w-4 h-4 text-[#007AFF]" />
                <span>Membrane Hydrophobic Belt &amp; Z-Axis Orientation</span>
                <span className="apple-pill bg-[#AF52DE]/15 text-[#007AFF] dark:text-[#0A84FF]">
                  {interData.membrane.classification}
                </span>
              </div>
              <p className="text-[11.5px] text-[#86868B] mt-0.5">
                Estimates the 30 Å hydrocarbon bilayer core normal and aligns channel/pore axes cleanly along +Z for GROMACS / CHARMM-GUI membrane MD.
              </p>
            </div>
            <div className="flex items-center gap-1.5">
              <button
                onClick={() =>
                  onToggleMembraneOverlay(membraneOverlayActive ? null : interData.membrane)
                }
                className={membraneOverlayActive ? 'apple-btn-primary text-[12px]' : 'apple-btn-secondary text-[12px]'}
              >
                <Layers className="w-3.5 h-3.5" />
                <span>{membraneOverlayActive ? 'Hide 3D Bilayer' : 'Show 3D Bilayer'}</span>
              </button>
              <button disabled={busy} onClick={handleOrientZ} className="apple-btn-primary text-[12px]">
                Align Axis → +Z
              </button>
            </div>
          </div>

          <div className="grid grid-cols-2 sm:grid-cols-4 gap-2 text-center">
            <div className="p-2.5 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#86868B]">Tilt from +Z</div>
              <div className="text-[15px] font-bold font-mono">{interData.membrane.tilt_from_z_deg}°</div>
            </div>
            <div className="p-2.5 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#86868B]">Belt Hydrophobic</div>
              <div className="text-[15px] font-bold font-mono text-[#FF9500]">
                {interData.membrane.belt_hydrophobic_pct}%
              </div>
            </div>
            <div className="p-2.5 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#86868B]">Flanks Hydrophobic</div>
              <div className="text-[15px] font-bold font-mono">{interData.membrane.flank_hydrophobic_pct}%</div>
            </div>
            <div className="p-2.5 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#86868B]">Slab Thickness</div>
              <div className="text-[15px] font-bold font-mono">{interData.membrane.thickness_A} Å</div>
            </div>
          </div>
        </div>
      )}

      {/* 3. Non-Covalent Interaction Networks (Salt Bridges, Pi-Pi Stacking, Inter-Chain Interface) */}
      {interData && (
        <div className="apple-card p-4 space-y-3">
          <div className="flex flex-wrap items-center justify-between gap-2">
            <div className="flex items-center gap-2 text-[13.5px] font-semibold">
              <Magnet className="w-4 h-4 text-[#34C759]" />
              <span>Non-Covalent Interaction Networks</span>
            </div>
            <div className="apple-segmented">
              <button
                onClick={() => setNetTab('salt')}
                className={`apple-segmented-item ${netTab === 'salt' ? 'apple-segmented-active' : 'text-[#6E6E73]'}`}
              >
                Salt Bridges ({interData.counts.salt_bridges})
              </button>
              <button
                onClick={() => setNetTab('pi')}
                className={`apple-segmented-item ${netTab === 'pi' ? 'apple-segmented-active' : 'text-[#6E6E73]'}`}
              >
                π–π Stacking ({interData.counts.pi_stacking})
              </button>
              <button
                onClick={() => setNetTab('interface')}
                className={`apple-segmented-item ${netTab === 'interface' ? 'apple-segmented-active' : 'text-[#6E6E73]'}`}
              >
                Interfaces ({interData.counts.interface_contacts})
              </button>
            </div>
          </div>

          <div className="max-h-60 overflow-y-auto space-y-1 text-[12px]">
            {netTab === 'salt' &&
              (interData.salt_bridges.length === 0 ? (
                <div className="text-[#86868B] py-2">No salt bridges within 4.0 Å detected.</div>
              ) : (
                interData.salt_bridges.map((sb, i) => (
                  <button
                    key={i}
                    onClick={() =>
                      onSelectResidue({ chain: sb.acid_chain, resseq: sb.acid_resseq, icode: '' }, true)
                    }
                    className="w-full flex items-center justify-between px-3 py-1.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03] hover:bg-[#007AFF]/10 transition text-left font-mono"
                  >
                    <span>
                      <span className="text-[#FF3B30] font-semibold">{sb.acid}</span> ({sb.acid_atom}) ⟷{' '}
                      <span className="text-[#007AFF] font-semibold">{sb.base}</span> ({sb.base_atom})
                    </span>
                    <span className="flex items-center gap-1.5">
                      {sb.inter_chain && (
                        <span className="apple-pill bg-[#AF52DE]/15 text-[#007AFF] dark:text-[#0A84FF]">Inter-chain</span>
                      )}
                      <span className="text-[#86868B]">{sb.distance_A} Å</span>
                    </span>
                  </button>
                ))
              ))}

            {netTab === 'pi' &&
              (interData.pi_stacking.length === 0 ? (
                <div className="text-[#86868B] py-2">No aromatic π–π stacking pairs within 5.6 Å detected.</div>
              ) : (
                interData.pi_stacking.map((pi, i) => (
                  <button
                    key={i}
                    onClick={() =>
                      onSelectResidue({ chain: pi.chain1, resseq: pi.resseq1, icode: '' }, true)
                    }
                    className="w-full flex items-center justify-between px-3 py-1.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03] hover:bg-[#007AFF]/10 transition text-left font-mono"
                  >
                    <span>
                      <span className="font-semibold">{pi.res1}</span> ⟷{' '}
                      <span className="font-semibold">{pi.res2}</span>
                    </span>
                    <span className="text-[#86868B]">
                      {pi.geometry} · {pi.distance_A} Å ({pi.angle_deg}°)
                    </span>
                  </button>
                ))
              ))}

            {netTab === 'interface' &&
              (interData.interface_contacts.length === 0 ? (
                <div className="text-[#86868B] py-2">
                  No inter-chain protein–protein interface contacts (single-chain or separated domains).
                </div>
              ) : (
                interData.interface_contacts.map((ic, i) => (
                  <button
                    key={i}
                    onClick={() =>
                      onSelectResidue({ chain: ic.chain1, resseq: ic.resseq1, icode: '' }, true)
                    }
                    className="w-full flex items-center justify-between px-3 py-1.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03] hover:bg-[#007AFF]/10 transition text-left font-mono"
                  >
                    <span>
                      <span className="font-semibold">{ic.res1}</span> ⟷{' '}
                      <span className="font-semibold">{ic.res2}</span>
                    </span>
                    <span className="text-[#86868B]">{ic.distance_A} Å</span>
                  </button>
                ))
              ))}
          </div>
        </div>
      )}
    </div>
  );
};
