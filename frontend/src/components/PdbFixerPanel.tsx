import React, { useEffect, useState } from 'react';
import { CheckCircle2, Play, Trash2, Wrench } from 'lucide-react';
import { api } from '../api/client';
import type { SelectedResidue, StructureSessionState } from '../types/pdbprep';
import type { ChainHighlight } from './Viewer3D';

interface PdbFixerPanelProps {
  session: StructureSessionState;
  busy: boolean;
  onHighlightChains: (hl: ChainHighlight | null, zoom?: boolean) => void;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
  onStateUpdate: (state: StructureSessionState, toast: string) => void;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
}

export const PdbFixerPanel: React.FC<PdbFixerPanelProps> = ({
  session,
  busy,
  onHighlightChains,
  onSelectResidue,
  onStateUpdate,
  onSetBusy,
  onError,
}) => {
  const diag = session.pdbfixer;

  const [replaceNonstandard, setReplaceNonstandard] = useState(true);
  const [resolveAltlocs, setResolveAltlocs] = useState<'highest_occupancy' | 'lowest_bfactor' | 'keep_A' | 'keep_all'>(
    'highest_occupancy'
  );
  const [addMissingAtoms, setAddMissingAtoms] = useState(true);
  const [addMissingResidues, setAddMissingResidues] = useState(true);
  const [maxLoopLength, setMaxLoopLength] = useState(12);
  const [capTermini, setCapTermini] = useState<'none' | 'ACE_NME' | 'charged'>('none');
  const [assignDisulfides, setAssignDisulfides] = useState(true);
  const [disulfideCutoff, setDisulfideCutoff] = useState(2.5);
  const [removeWaters, setRemoveWaters] = useState<'bulk' | 'all' | 'none'>('bulk');
  const [waterKeepDist, setWaterKeepDist] = useState(3.5);
  const [removeBuffers, setRemoveBuffers] = useState(true);
  const [keepLigands, setKeepLigands] = useState(true);
  const [keepCofactors, setKeepCofactors] = useState(true);
  const [keepMetals, setKeepMetals] = useState(true);
  const [renumber, setRenumber] = useState(false);
  const [selectedChains, setSelectedChains] = useState<string[]>(diag.chains.map((c) => c.id));

  useEffect(() => {
    setSelectedChains(diag.chains.map((c) => c.id));
  }, [session.session_id, diag.chains.length]);

  const uncheckedChains = diag.chains.map((c) => c.id).filter((id) => !selectedChains.includes(id));

  const handleDeleteUncheckedChains = async () => {
    if (uncheckedChains.length === 0) return;
    onSetBusy(`Deleting chain(s) ${uncheckedChains.join(', ')}…`);
    try {
      const resp = await api.deleteChains(session.session_id, uncheckedChains);
      onStateUpdate(resp.state, resp.result.message);
    } catch (e: any) {
      onError(e?.message || 'Failed to delete chain(s)');
    } finally {
      onSetBusy(null);
    }
  };

  const toggleChain = (cid: string) => {
    const next = selectedChains.includes(cid)
      ? selectedChains.length > 1
        ? selectedChains.filter((c) => c !== cid)
        : selectedChains
      : [...selectedChains, cid];
    setSelectedChains(next);
    const removed = diag.chains.map((c) => c.id).filter((id) => !next.includes(id));
    if (removed.length > 0) {
      onHighlightChains({ ids: removed, tone: 'delete' }, true);
    } else {
      onHighlightChains({ ids: [cid], tone: 'select' }, true);
    }
  };

  const handleRun = async () => {
    onSetBusy('Running Deep PDBFixer Pipeline…');
    try {
      const resp = await api.runPdbFixer(session.session_id, {
        replace_nonstandard: replaceNonstandard,
        resolve_altlocs: resolveAltlocs,
        add_missing_atoms: addMissingAtoms,
        add_missing_residues: addMissingResidues,
        max_loop_length: maxLoopLength,
        cap_termini: capTermini,
        assign_disulfides: assignDisulfides,
        disulfide_cutoff: disulfideCutoff,
        remove_waters: removeWaters,
        water_keep_distance: waterKeepDist,
        remove_crystallization_artifacts: removeBuffers,
        keep_ligands: keepLigands,
        keep_cofactors: keepCofactors,
        keep_metals: keepMetals,
        keep_chains: selectedChains.length ? selectedChains : null,
        renumber_residues: renumber,
        fix_insertion_codes: true,
        standardize_atom_names: true,
      });
      onStateUpdate(resp.state, `PDBFixer completed (${resp.result.actions.length} repairs applied)`);
    } catch (e: any) {
      onError(e?.message || 'PDBFixer execution failed');
    } finally {
      onSetBusy(null);
    }
  };

  return (
    <div className="space-y-4 pb-6">
      <div className="apple-card p-4 flex flex-wrap items-center justify-between gap-3">
        <div>
          <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#007AFF] dark:text-[#0A84FF]">
            PDBFixer Structural Repair
          </div>
          <h2 className="text-[16px] font-semibold tracking-[-0.015em]">
            Deep Topology, Missing Loops, Atoms, Disulfides & Solvent Curation
          </h2>
        </div>
        <button disabled={busy} onClick={handleRun} className="apple-btn-primary">
          <Play className="w-3.5 h-3.5 fill-current" />
          <span>Run PDBFixer</span>
        </button>
      </div>

      {/* Live Pre-Fix Diagnostics Grid */}
      <div className="grid grid-cols-2 sm:grid-cols-3 gap-2.5">
        {[
          { label: 'Non-Standard Residues', val: diag.nonstandard_residues.length, warn: diag.nonstandard_residues.length > 0 },
          { label: 'Missing Heavy Atoms', val: diag.missing_atoms.length, warn: diag.missing_atoms.length > 0 },
          { label: 'Missing Loops / Breaks', val: diag.missing_loops.length, warn: diag.missing_loops.length > 0 },
          { label: 'AltLoc Conformers', val: diag.altlocs.length, warn: diag.altlocs.length > 0 },
          { label: 'Disulfide Candidates', val: diag.disulfides.length, warn: false },
          { label: 'Bulk / Struct Waters', val: `${diag.waters.bulk} / ${diag.waters.structural}`, warn: diag.waters.bulk > 0 },
        ].map((card) => (
          <div key={card.label} className="apple-card p-3">
            <div className="text-[11px] text-[#86868B]">{card.label}</div>
            <div
              className={`text-[18px] font-bold font-mono mt-0.5 ${
                card.warn ? 'text-[#FF9500]' : 'text-[#34C759]'
              }`}
            >
              {card.val}
            </div>
          </div>
        ))}
      </div>

      {/* Detailed Detected Anomalies (clickable to zoom 3D) */}
      {(diag.nonstandard_residues.length > 0 ||
        diag.missing_atoms.length > 0 ||
        diag.missing_loops.length > 0 ||
        diag.altlocs.length > 0) && (
        <div className="apple-card p-4 space-y-3">
          <h3 className="text-[13px] font-semibold">Detected Structural Defects (Click to Focus in 3D)</h3>
          <div className="space-y-1.5 text-[12px]">
            {diag.nonstandard_residues.map((ns) => (
              <button
                key={`ns-${ns.chain}:${ns.resseq}`}
                onClick={() => onSelectResidue({ chain: ns.chain, resseq: ns.resseq, icode: ns.icode, resname: ns.resname }, true)}
                className="w-full flex items-center justify-between px-3 py-1.5 rounded-xl bg-[#FF9500]/10 hover:bg-[#FF9500]/18 transition text-left"
              >
                <span className="font-mono font-semibold">
                  Non-standard: {ns.resname} {ns.chain}:{ns.resseq}
                </span>
                <span className="text-[#B26200] dark:text-[#FF9F0A] font-medium">→ Convert to {ns.replacement}</span>
              </button>
            ))}
            {diag.missing_atoms.map((ma) => (
              <button
                key={`ma-${ma.chain}:${ma.resseq}-${ma.missing.join('')}`}
                onClick={() => onSelectResidue({ chain: ma.chain, resseq: ma.resseq, icode: ma.icode, resname: ma.resname }, true)}
                className="w-full flex items-center justify-between px-3 py-1.5 rounded-xl bg-[#FF3B30]/10 hover:bg-[#FF3B30]/18 transition text-left"
              >
                <span className="font-mono font-semibold">
                  Incomplete: {ma.resname} {ma.chain}:{ma.resseq}
                </span>
                <span className="font-mono text-[#D70015] dark:text-[#FF453A]">Missing: {ma.missing.join(', ')}</span>
              </button>
            ))}
            {diag.missing_loops.map((ml, i) => (
              <button
                key={`ml-${i}`}
                onClick={() => onSelectResidue({ chain: ml.chain, resseq: ml.start_resseq, icode: '' }, true)}
                className="w-full flex items-center justify-between px-3 py-1.5 rounded-xl bg-[#AF52DE]/10 hover:bg-[#AF52DE]/18 transition text-left"
              >
                <span className="font-mono font-semibold">
                  Loop Break: Chain {ml.chain} ({ml.start_resseq} → {ml.end_resseq})
                </span>
                <span className="text-[#8944AB] dark:text-[#BF5AF2]">
                  {ml.missing_count} missing ({ml.sequence.join('-') || `${ml.gap_distance_A} Å gap`})
                </span>
              </button>
            ))}
            {diag.altlocs.map((al) => (
              <button
                key={`al-${al.chain}:${al.resseq}`}
                onClick={() => onSelectResidue({ chain: al.chain, resseq: al.resseq, icode: al.icode, resname: al.resname }, true)}
                className="w-full flex items-center justify-between px-3 py-1.5 rounded-xl bg-[#007AFF]/10 hover:bg-[#007AFF]/18 transition text-left"
              >
                <span className="font-mono font-semibold">
                  AltLoc: {al.resname} {al.chain}:{al.resseq}
                </span>
                <span className="font-mono text-[#007AFF] dark:text-[#0A84FF]">
                  {al.conformers.map((c) => `${c.id}(occ=${c.occ}, B=${c.b_factor})`).join(' / ')}
                </span>
              </button>
            ))}
          </div>
        </div>
      )}

      {/* Configuration Controls */}
      <div className="apple-card p-4 space-y-4">
        <div className="flex items-center gap-2 text-[13.5px] font-semibold">
          <Wrench className="w-4 h-4 text-[#007AFF]" />
          <span>Repair & Curation Parameters</span>
        </div>

        {/* Chains selector */}
        <div>
          <div className="flex items-center justify-between gap-2">
            <label className="text-[12px] font-medium text-[#6E6E73] dark:text-[#A1A1A6]">
              Retain Polymer Chains (click to unselect chains you want to remove)
            </label>
            {uncheckedChains.length > 0 && (
              <button
                disabled={busy}
                onClick={handleDeleteUncheckedChains}
                className="apple-btn-danger !py-1 !px-2.5 text-[11.5px]"
              >
                <Trash2 className="w-3.5 h-3.5" />
                Delete {uncheckedChains.join(', ')} Now
              </button>
            )}
          </div>
          <div className="flex flex-wrap gap-1.5 mt-1.5">
            {diag.chains.map((ch) => {
              const on = selectedChains.includes(ch.id);
              return (
                <button
                  key={ch.id}
                  onClick={() => toggleChain(ch.id)}
                  className={`px-3 py-1 rounded-xl text-[12px] font-mono font-semibold border transition ${
                    on
                      ? 'bg-[#007AFF] text-white border-[#007AFF]'
                      : 'bg-black/[0.04] dark:bg-white/[0.06] border-transparent text-[#86868B]'
                  }`}
                >
                  Chain {ch.id} ({ch.protein_residues} res)
                </button>
              );
            })}
          </div>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-3 text-[12.5px]">
          <label className="flex items-center justify-between p-2.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03] cursor-pointer">
            <span>Convert Non-Standard Residues (MSE→MET, SEP→SER…)</span>
            <input
              type="checkbox"
              checked={replaceNonstandard}
              onChange={(e) => setReplaceNonstandard(e.target.checked)}
              className="accent-[#007AFF] w-4 h-4"
            />
          </label>

          <label className="flex items-center justify-between p-2.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03] cursor-pointer">
            <span>Rebuild Missing Heavy Atoms & C-term OXT</span>
            <input
              type="checkbox"
              checked={addMissingAtoms}
              onChange={(e) => setAddMissingAtoms(e.target.checked)}
              className="accent-[#007AFF] w-4 h-4"
            />
          </label>

          <div className="flex items-center justify-between p-2.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03]">
            <span>Rebuild Missing Loops (≤ {maxLoopLength} res)</span>
            <div className="flex items-center gap-2">
              <input
                type="number"
                min={1}
                max={30}
                value={maxLoopLength}
                onChange={(e) => setMaxLoopLength(Number(e.target.value))}
                className="w-14 apple-input !py-0.5 text-center font-mono"
              />
              <input
                type="checkbox"
                checked={addMissingResidues}
                onChange={(e) => setAddMissingResidues(e.target.checked)}
                className="accent-[#007AFF] w-4 h-4"
              />
            </div>
          </div>

          <div className="flex items-center justify-between p-2.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03]">
            <span>AltLoc Resolution</span>
            <select
              value={resolveAltlocs}
              onChange={(e) => setResolveAltlocs(e.target.value as any)}
              className="apple-select text-[12px]"
            >
              <option value="highest_occupancy">Highest Occupancy</option>
              <option value="lowest_bfactor">Lowest B-Factor</option>
              <option value="keep_A">Always Conformer A</option>
              <option value="keep_all">Keep All (No Change)</option>
            </select>
          </div>

          <div className="flex items-center justify-between p-2.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03]">
            <span>Terminal Capping</span>
            <select
              value={capTermini}
              onChange={(e) => setCapTermini(e.target.value as any)}
              className="apple-select text-[12px]"
            >
              <option value="none">Standard Charged (NH3+ / COO-)</option>
              <option value="ACE_NME">Neutral Caps (ACE / NME)</option>
            </select>
          </div>

          <div className="flex items-center justify-between p-2.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03]">
            <span>Disulfide Bridges (SG–SG ≤ {disulfideCutoff} Å)</span>
            <div className="flex items-center gap-2">
              <input
                type="number"
                step={0.1}
                min={1.8}
                max={3.2}
                value={disulfideCutoff}
                onChange={(e) => setDisulfideCutoff(Number(e.target.value))}
                className="w-16 apple-input !py-0.5 text-center font-mono"
              />
              <input
                type="checkbox"
                checked={assignDisulfides}
                onChange={(e) => setAssignDisulfides(e.target.checked)}
                className="accent-[#007AFF] w-4 h-4"
              />
            </div>
          </div>

          <div className="flex items-center justify-between p-2.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03]">
            <span>Crystallographic Waters</span>
            <select
              value={removeWaters}
              onChange={(e) => setRemoveWaters(e.target.value as any)}
              className="apple-select text-[12px]"
            >
              <option value="bulk">Smart (Keep Active-Site ≤{waterKeepDist}Å, Strip Bulk)</option>
              <option value="all">Remove All Waters</option>
              <option value="none">Keep All Waters</option>
            </select>
          </div>

          <label className="flex items-center justify-between p-2.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03] cursor-pointer">
            <span>Strip Buffers (SO4, EDO, GOL, PEG…)</span>
            <input
              type="checkbox"
              checked={removeBuffers}
              onChange={(e) => setRemoveBuffers(e.target.checked)}
              className="accent-[#007AFF] w-4 h-4"
            />
          </label>

          <div className="flex items-center justify-between p-2.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03]">
            <span>Retain Biological Heterogens</span>
            <div className="flex items-center gap-3 text-[11.5px]">
              <label className="flex items-center gap-1 cursor-pointer">
                <input type="checkbox" checked={keepLigands} onChange={(e) => setKeepLigands(e.target.checked)} className="accent-[#007AFF]" />
                Ligands
              </label>
              <label className="flex items-center gap-1 cursor-pointer">
                <input type="checkbox" checked={keepCofactors} onChange={(e) => setKeepCofactors(e.target.checked)} className="accent-[#007AFF]" />
                Cofactors
              </label>
              <label className="flex items-center gap-1 cursor-pointer">
                <input type="checkbox" checked={keepMetals} onChange={(e) => setKeepMetals(e.target.checked)} className="accent-[#007AFF]" />
                Metals
              </label>
            </div>
          </div>

          <label className="flex items-center justify-between p-2.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03] cursor-pointer">
            <span>Renumber Residues Sequentially (1..N)</span>
            <input
              type="checkbox"
              checked={renumber}
              onChange={(e) => setRenumber(e.target.checked)}
              className="accent-[#007AFF] w-4 h-4"
            />
          </label>
        </div>

        {/* Disulfides Table */}
        {diag.disulfides.length > 0 && (
          <div className="pt-2 border-t border-black/[0.06] dark:border-white/[0.08] space-y-1.5">
            <div className="text-[12px] font-semibold">Detected Disulfide Pairs</div>
            {diag.disulfides.map((ss, idx) => (
              <div
                key={idx}
                className="flex items-center justify-between px-3 py-1.5 rounded-xl bg-[#34C759]/10 text-[12px] font-mono"
              >
                <span>
                  {ss.resname1} {ss.cys1} ⟷ {ss.resname2} {ss.cys2}
                </span>
                <span className="flex items-center gap-1.5 text-[#1B8A38] dark:text-[#30D158] font-semibold">
                  <CheckCircle2 className="w-3.5 h-3.5" />
                  {ss.distance_A} Å ({ss.bonded ? 'S–S Bond' : 'Reduced'})
                </span>
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
};
