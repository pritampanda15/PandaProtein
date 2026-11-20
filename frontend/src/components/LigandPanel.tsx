import React, { useState } from 'react';
import { Atom, Beaker, Download, Sparkles } from 'lucide-react';
import { api } from '../api/client';
import type { SelectedResidue, StructureSessionState } from '../types/pdbprep';

interface LigandPanelProps {
  session: StructureSessionState;
  busy: boolean;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
}

export const LigandPanel: React.FC<LigandPanelProps> = ({
  session,
  busy,
  onSelectResidue,
  onSetBusy,
  onError,
}) => {
  const { ligands, metals, buffers } = session.ligands_info;
  const [confs, setConfs] = useState<
    Record<string, { conformers: Array<{ rank: number; energy_kcal_mol: number; delta_e: number }>; sdf: string }>
  >({});

  const handleGenConformers = async (ligId: string) => {
    onSetBusy(`Generating 3D ETKDGv3 conformers for ${ligId}…`);
    try {
      const resp = await api.generateConformers(session.session_id, ligId, 5);
      setConfs((prev) => ({ ...prev, [ligId]: { conformers: resp.conformers, sdf: resp.sdf } }));
    } catch (e: any) {
      onError(e?.message || 'Conformer generation failed');
    } finally {
      onSetBusy(null);
    }
  };

  const downloadSdf = (ligId: string, sdf: string) => {
    const blob = new Blob([sdf], { type: 'chemical/x-mdl-sdfile' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${ligId.replace(/:/g, '_')}_conformers.sdf`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const handleExportLigand = async (ligId: string, fmt: 'pdbqt' | 'mol2' | 'sdf' | 'pdb') => {
    onSetBusy(`Preparing ${ligId} as ${fmt.toUpperCase()} (pH ${session.propka.ph})…`);
    try {
      await api.exportLigand(session.session_id, ligId, fmt, session.propka.ph);
    } catch (e: any) {
      onError(e?.message || `Ligand ${fmt.toUpperCase()} export failed`);
    } finally {
      onSetBusy(null);
    }
  };

  return (
    <div className="space-y-4 pb-6">
      <div className="apple-card p-4">
        <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#AF52DE] dark:text-[#BF5AF2]">
          RDKit Cheminformatics & Metal Coordination Analyzer
        </div>
        <h2 className="text-[16px] font-semibold tracking-[-0.015em]">
          Co-Crystallized Ligands, Cofactors, Binding Pocket Contacts & Catalytic Metals
        </h2>
      </div>

      {/* Drug-like Ligands & Cofactors */}
      {ligands.length === 0 ? (
        <div className="apple-card p-5 text-[13px] text-[#86868B]">
          No organic ligands or cofactors detected in this structure (apo protein).
        </div>
      ) : (
        ligands.map((lig) => (
          <div key={lig.id} className="apple-card p-4 space-y-4">
            <div className="flex flex-wrap items-start justify-between gap-3">
              <div>
                <div className="flex items-center gap-2">
                  <Beaker className="w-4 h-4 text-[#AF52DE]" />
                  <span className="text-[16px] font-bold font-mono">
                    {lig.resname} ({lig.chain}:{lig.resseq})
                  </span>
                  <span className="apple-pill bg-[#AF52DE]/15 text-[#8944AB] dark:text-[#BF5AF2] uppercase">
                    {lig.category}
                  </span>
                  {lig.drug_like !== undefined && (
                    <span
                      className={`apple-pill ${
                        lig.drug_like
                          ? 'bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158]'
                          : 'bg-[#FF9500]/15 text-[#B26200] dark:text-[#FF9F0A]'
                      }`}
                    >
                      {lig.drug_like ? 'Lipinski Ro5 Pass' : `${lig.lipinski_violations} Ro5 Violation(s)`}
                    </span>
                  )}
                </div>
                {lig.smiles && (
                  <div className="text-[11px] font-mono text-[#86868B] mt-1 break-all">
                    SMILES: {lig.smiles}
                  </div>
                )}
              </div>

              <div className="flex items-center gap-1.5">
                <button
                  onClick={() =>
                    onSelectResidue({ chain: lig.chain, resseq: lig.resseq, icode: lig.icode, resname: lig.resname }, true)
                  }
                  className="apple-btn-secondary text-[12px]"
                >
                  Focus Pocket
                </button>
                <button
                  disabled={busy}
                  onClick={() => handleGenConformers(lig.id)}
                  className="apple-btn-primary text-[12px]"
                >
                  <Sparkles className="w-3.5 h-3.5" />
                  ETKDGv3 Conformers
                </button>
              </div>
            </div>

            {/* 2D Depiction + RDKit Physicochemical Properties */}
            <div className="grid grid-cols-1 md:grid-cols-12 gap-3 items-center">
              {lig.svg && (
                <div
                  className="md:col-span-5 rounded-2xl bg-white p-2 border border-black/[0.06] flex items-center justify-center overflow-hidden"
                  dangerouslySetInnerHTML={{ __html: lig.svg }}
                />
              )}
              <div className={`${lig.svg ? 'md:col-span-7' : 'md:col-span-12'} grid grid-cols-3 gap-2 text-center`}>
                {[
                  ['Formula', lig.formula || '—'],
                  ['MW (Da)', lig.mw ?? '—'],
                  ['cLogP', lig.clogp ?? '—'],
                  ['TPSA (Å²)', lig.tpsa ?? '—'],
                  ['HBD / HBA', lig.hbd !== undefined ? `${lig.hbd} / ${lig.hba}` : '—'],
                  ['Rot Bonds', lig.rotatable_bonds ?? '—'],
                ].map(([k, v]) => (
                  <div key={String(k)} className="p-2.5 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
                    <div className="text-[10px] text-[#86868B]">{k}</div>
                    <div className="text-[13px] font-semibold font-mono mt-0.5">{v}</div>
                  </div>
                ))}
              </div>
            </div>

            {/* 1-Click Docking & MD Ligand Preparation Exports */}
            <div className="p-3 rounded-2xl bg-black/[0.025] dark:bg-white/[0.035] border border-black/[0.05] dark:border-white/[0.06] flex flex-wrap items-center justify-between gap-2">
              <div>
                <div className="text-[11.5px] font-semibold">Docking &amp; MD Ligand Prep (pH {session.propka.ph})</div>
                <div className="text-[10.5px] text-[#86868B]">
                  3D protonated ligand with Gasteiger charges, Sybyl types &amp; AutoDock torsion tree
                </div>
              </div>
              <div className="flex flex-wrap items-center gap-1.5">
                <button
                  disabled={busy}
                  onClick={() => handleExportLigand(lig.id, 'pdbqt')}
                  className="apple-btn-primary !py-1 !px-2.5 text-[11.5px]"
                  title="AutoDock Vina / GNINA ligand .pdbqt with ROOT/BRANCH torsion tree"
                >
                  <Download className="w-3 h-3" /> .PDBQT (Vina)
                </button>
                <button
                  disabled={busy}
                  onClick={() => handleExportLigand(lig.id, 'mol2')}
                  className="apple-btn-secondary !py-1 !px-2.5 text-[11.5px]"
                  title="Tripos .mol2 with Sybyl atom types & Gasteiger charges (CGenFF / ACPYPE / AnteChamber)"
                >
                  <Download className="w-3 h-3" /> .MOL2 (MD)
                </button>
                <button
                  disabled={busy}
                  onClick={() => handleExportLigand(lig.id, 'sdf')}
                  className="apple-btn-secondary !py-1 !px-2 text-[11.5px]"
                >
                  .SDF
                </button>
                <button
                  disabled={busy}
                  onClick={() => handleExportLigand(lig.id, 'pdb')}
                  className="apple-btn-secondary !py-1 !px-2 text-[11.5px]"
                >
                  .PDB
                </button>
              </div>
            </div>

            {/* Generated Conformers Result */}
            {confs[lig.id] && (
              <div className="p-3 rounded-2xl bg-[#34C759]/10 border border-[#34C759]/25 space-y-2">
                <div className="flex items-center justify-between">
                  <span className="text-[12px] font-semibold text-[#1B8A38] dark:text-[#30D158]">
                    Generated {confs[lig.id].conformers.length} MMFF94-Minimized 3D Conformers
                  </span>
                  <button
                    onClick={() => downloadSdf(lig.id, confs[lig.id].sdf)}
                    className="apple-btn-primary !py-1 !px-2.5 !bg-[#34C759] hover:!bg-[#248A3D] text-[11.5px]"
                  >
                    <Download className="w-3 h-3" />
                    Download .SDF
                  </button>
                </div>
                <div className="flex flex-wrap gap-1.5 text-[11px] font-mono">
                  {confs[lig.id].conformers.map((c) => (
                    <span key={c.rank} className="px-2 py-0.5 rounded-lg bg-white/75 dark:bg-black/25">
                      #{c.rank}: {c.energy_kcal_mol} kcal/mol (Δ{c.delta_e})
                    </span>
                  ))}
                </div>
              </div>
            )}

            {/* Binding Pocket Contacts */}
            {lig.contacts.length > 0 && (
              <div className="space-y-1.5">
                <div className="text-[12px] font-semibold">Binding Pocket Contacts (≤ 4.2 Å)</div>
                <div className="flex flex-wrap gap-1.5">
                  {lig.contacts.slice(0, 18).map((ct, idx) => (
                    <button
                      key={idx}
                      onClick={() =>
                        onSelectResidue({ chain: ct.chain, resseq: ct.resseq, icode: ct.icode, resname: ct.resname }, true)
                      }
                      className={`px-2.5 py-1 rounded-xl text-[11px] font-mono border transition ${
                        ct.type === 'hbond'
                          ? 'bg-[#007AFF]/10 border-[#007AFF]/25 text-[#007AFF] dark:text-[#0A84FF]'
                          : ct.type === 'hydrophobic'
                          ? 'bg-[#FF9500]/10 border-[#FF9500]/25 text-[#B26200] dark:text-[#FF9F0A]'
                          : 'bg-black/[0.035] dark:bg-white/[0.05] border-transparent'
                      }`}
                    >
                      {ct.resname} {ct.key} ({ct.type === 'hbond' ? 'H-bond' : ct.type}, {ct.distance_A} Å)
                    </button>
                  ))}
                </div>
              </div>
            )}
          </div>
        ))
      )}

      {/* Metal Coordination Spheres */}
      {metals.length > 0 && (
        <div className="apple-card p-4 space-y-3">
          <div className="flex items-center gap-2 text-[13.5px] font-semibold">
            <Atom className="w-4 h-4 text-[#AF52DE]" />
            <span>Catalytic & Structural Metal Sites ({metals.length})</span>
          </div>
          <div className="space-y-2.5">
            {metals.map((m) => (
              <div key={m.id} className="p-3 rounded-2xl bg-[#AF52DE]/8 border border-[#AF52DE]/20 space-y-2">
                <div className="flex items-center justify-between">
                  <div className="font-mono font-bold text-[13.5px]">
                    {m.resname} ({m.chain}:{m.resseq}) — {m.geometry}
                  </div>
                  <button
                    onClick={() => onSelectResidue({ chain: m.chain, resseq: m.resseq, icode: '', resname: m.resname }, true)}
                    className="apple-btn-secondary !py-1 text-[11.5px]"
                  >
                    Center Metal
                  </button>
                </div>
                <div className="flex flex-wrap gap-1.5">
                  {m.coordinators.map((co, idx) => (
                    <button
                      key={idx}
                      onClick={() =>
                        onSelectResidue({ chain: co.chain, resseq: co.resseq, icode: '', resname: co.resname }, true)
                      }
                      className="px-2.5 py-1 rounded-xl bg-white/75 dark:bg-white/10 text-[11.5px] font-mono hover:bg-[#AF52DE]/20 transition"
                    >
                      {co.resname} {co.key}.{co.atom} ({co.distance_A} Å)
                    </button>
                  ))}
                </div>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Detected 3D Druggable Binding Pockets (Cavity Finder for Apo & Holo Proteins) */}
      {session.pockets && session.pockets.length > 0 && (
        <div className="apple-card p-4 space-y-3">
          <div className="flex flex-wrap items-center justify-between gap-2">
            <div>
              <div className="text-[13.5px] font-semibold flex items-center gap-2">
                <Sparkles className="w-4 h-4 text-[#5856D6]" />
                <span>Detected 3D Druggable Binding Pockets ({session.pockets.length})</span>
              </div>
              <p className="text-[11.5px] text-[#86868B] mt-0.5">
                Toggle &ldquo;Pockets&rdquo; in the 3D Viewport to inspect cavities live; all pockets are auto-available as Vina docking boxes.
              </p>
            </div>
          </div>

          <div className="space-y-2.5">
            {session.pockets.map((pk) => (
              <div
                key={pk.id}
                className="p-3 rounded-2xl bg-[#5856D6]/[0.06] border border-[#5856D6]/20 space-y-2"
              >
                <div className="flex flex-wrap items-center justify-between gap-2">
                  <div className="flex items-center gap-2">
                    <span className="font-mono font-bold text-[13px]">Pocket #{pk.rank}</span>
                    <span className="apple-pill bg-[#5856D6]/15 text-[#5856D6] dark:text-[#7D7AFF] font-mono">
                      Druggability {pk.druggability_score}%
                    </span>
                    <span className="apple-pill bg-black/[0.05] dark:bg-white/[0.08] font-mono">
                      {pk.volume_A3} Å³
                    </span>
                    <span className="apple-pill bg-[#FF9500]/12 text-[#B26200] dark:text-[#FF9F0A] font-mono">
                      {pk.hydrophobic_pct}% hydrophobic
                    </span>
                  </div>
                  <span className="text-[10.5px] font-mono text-[#86868B]">
                    Center ({pk.center.join(', ')}) · Box {pk.size.join('×')} Å
                  </span>
                </div>
                <div className="flex flex-wrap gap-1">
                  {pk.lining_residues.slice(0, 16).map((lr) => (
                    <button
                      key={lr.key}
                      onClick={() =>
                        onSelectResidue(
                          { chain: lr.chain, resseq: lr.resseq, icode: lr.icode, resname: lr.resname, key: lr.key },
                          true
                        )
                      }
                      className="px-2 py-0.5 rounded-lg text-[10.5px] font-mono bg-white/80 dark:bg-white/10 hover:bg-[#5856D6]/20 transition"
                    >
                      {lr.resname} {lr.key}
                    </button>
                  ))}
                </div>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Crystallization Buffers Alert */}
      {buffers.length > 0 && (
        <div className="apple-card p-4 flex items-center justify-between">
          <div>
            <div className="text-[13px] font-semibold text-[#B26200] dark:text-[#FF9F0A]">
              {buffers.length} Crystallization Buffer Artifact(s) Present
            </div>
            <div className="text-[11.5px] text-[#86868B] font-mono mt-0.5">
              {buffers.map((b) => `${b.resname} (${b.chain}:${b.resseq})`).join(', ')} — Strip in PDBFixer before MD/Docking
            </div>
          </div>
        </div>
      )}
    </div>
  );
};
