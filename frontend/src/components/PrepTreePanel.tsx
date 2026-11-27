import React, { useMemo, useRef, useState } from 'react';
import {
  AlertTriangle,
  Atom,
  CheckCircle2,
  ChevronDown,
  ChevronRight,
  Database,
  Dna,
  Download,
  FileCheck2,
  FileJson,
  Layers,
  RefreshCw,
  RotateCcw,
  Search,
  ShieldCheck,
  Sparkles,
  Wrench,
} from 'lucide-react';
import { api } from '../api/client';
import type {
  NavTab,
  ProvenanceVerifyResult,
  SelectedResidue,
  StructureSessionState,
} from '../types/pdbprep';

interface PrepTreePanelProps {
  session: StructureSessionState;
  onUpdate: (next: StructureSessionState, message: string) => void;
  onError: (msg: string) => void;
  onBusy: (label: string | null) => void;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
  onJumpTab: (tab: NavTab) => void;
}

type TreeFilter = 'all' | 'issues' | 'heterogens' | 'gaps';

export const PrepTreePanel: React.FC<PrepTreePanelProps> = ({
  session,
  onUpdate,
  onError,
  onBusy,
  onSelectResidue,
  onJumpTab,
}) => {
  const [filter, setFilter] = useState<TreeFilter>('all');
  const [query, setQuery] = useState('');
  const [dedupAxialIons, setDedupAxialIons] = useState(true);
  const [strictGapGate, setStrictGapGate] = useState(true);
  const [expandedChains, setExpandedChains] = useState<Record<string, boolean>>({});
  const [expandedWaters, setExpandedWaters] = useState<Record<string, boolean>>({});
  const [ionChargeDrafts, setIonChargeDrafts] = useState<Record<string, string>>({});
  const [provPreview, setProvPreview] = useState<Record<string, unknown> | null>(null);
  const [verifyResult, setVerifyResult] = useState<ProvenanceVerifyResult | null>(null);
  const verifyFileRef = useRef<HTMLInputElement | null>(null);

  const assemblies = session.assemblies || [];
  const activeAsm = session.active_assembly_id || 'ASU';
  const prepTree = session.prep_tree;
  const dockingPrefs = session.docking_prefs || {
    residue_roles: {},
    chain_roles: {},
    ion_charges: {},
    ion_substitutions: {},
    include_nucleic: true,
  };
  const rebuiltGaps = session.missing_regions?.rebuilt_gaps || [];
  const ad4Subs = prepTree?.ad4_substitutes || [];

  const isChainOpen = (cid: string, idx: number) =>
    expandedChains[cid] !== undefined ? expandedChains[cid] : idx < 2;

  const toggleChain = (cid: string, idx: number) => {
    setExpandedChains((prev) => ({ ...prev, [cid]: !isChainOpen(cid, idx) }));
  };

  const handleSwitchAssembly = async (asmId: string) => {
    if (asmId === activeAsm) return;
    onBusy(`Building ${asmId === 'ASU' ? 'Asymmetric Unit (ASU)' : `Biomol ${asmId}`}…`);
    try {
      const r = await api.switchAssembly(session.session_id, asmId, dedupAxialIons);
      onUpdate(r.state, r.result.message);
    } catch (e: any) {
      onError(e?.message || 'Failed to switch biological assembly');
    } finally {
      onBusy(null);
    }
  };

  const handleUpdatePrefs = async (
    payload: Parameters<typeof api.updateDockingPrefs>[1],
    toastOverride?: string
  ) => {
    try {
      const r = await api.updateDockingPrefs(session.session_id, payload);
      onUpdate(r.state, toastOverride || r.result.message);
    } catch (e: any) {
      onError(e?.message || 'Failed to update docking preferences');
    }
  };

  const handleRebuildGap = async (
    chain: string,
    startResseq: number,
    endResseq: number,
    sequence3: string[],
    strict: boolean
  ) => {
    onBusy(`Rebuilding gap ${chain}:${startResseq}-${endResseq} (testing 3 seeds)…`);
    try {
      const r = await api.rebuildGap(session.session_id, {
        chain,
        start_resseq: startResseq,
        end_resseq: endResseq,
        sequence_3: sequence3,
        strict_geometry: strict,
        max_seeds: 3,
      });
      if (r.result.refused) {
        onUpdate(r.state, `Gap ${r.result.gap_id} refused by geometry gate`);
        onError(r.result.reason || 'Loop failed geometry validation');
      } else {
        onUpdate(r.state, r.result.message || `Rebuilt gap ${r.result.gap_id}`);
      }
    } catch (e: any) {
      onError(e?.message || 'Gap rebuild failed');
    } finally {
      onBusy(null);
    }
  };

  const handlePickAltloc = async (chain: string, resseq: number, icode: string, altloc: string) => {
    onBusy(`Selecting AltLoc ${altloc} on ${chain}:${resseq}…`);
    try {
      const r = await api.cootAltloc(session.session_id, { chain, resseq, icode, altloc });
      onUpdate(r.state, r.result.message);
    } catch (e: any) {
      onError(e?.message || 'AltLoc update failed');
    } finally {
      onBusy(null);
    }
  };

  const handleVerifyRcsb = async () => {
    onBusy('Querying RCSB Data API for official entry ligands…');
    try {
      const r = await api.lookupRcsbLigands(session.session_id, undefined, true);
      if (r.result.source === 'rcsb_data_api') {
        onUpdate(
          r.state,
          `RCSB verified ${r.result.pdb_id}: ${r.result.ligands.length} official ligand component(s)`
        );
      } else {
        onUpdate(r.state, 'Used offline CCD classification (entry not found on RCSB or offline)');
      }
    } catch (e: any) {
      onError(e?.message || 'RCSB lookup failed');
    } finally {
      onBusy(null);
    }
  };

  const handleExport = async (target: 'pdbqt' | 'rec_h_pdb' | 'provenance_json', suffix: string) => {
    const stem = session.name.replace(/\.\w+$/, '');
    onBusy(`Generating ${stem}${suffix}…`);
    try {
      await api.exportFile(session.session_id, { target }, `${stem}${suffix}`);
    } catch (e: any) {
      onError(e?.message || `Failed to export ${target}`);
    } finally {
      onBusy(null);
    }
  };

  const handleExportLigandPdbqt = async (ligKey: string, resname: string) => {
    onBusy(`Preparing ${resname} (${ligKey}) as AutoDock .pdbqt…`);
    try {
      await api.exportLigand(session.session_id, ligKey, 'pdbqt', session.propka.ph);
    } catch (e: any) {
      onError(e?.message || 'Ligand PDBQT export failed');
    } finally {
      onBusy(null);
    }
  };

  const handleInspectProvenance = async () => {
    if (provPreview) {
      setProvPreview(null);
      return;
    }
    try {
      const data = await api.getProvenance(session.session_id);
      setProvPreview(data);
    } catch (e: any) {
      onError(e?.message || 'Failed to load provenance');
    }
  };

  const handleVerifyProvenanceFile = async (file: File) => {
    try {
      const text = await file.text();
      const parsed = JSON.parse(text);
      const res = await api.verifyProvenance(session.session_id, parsed);
      setVerifyResult(res);
    } catch (e: any) {
      onError(e?.message || 'Invalid provenance JSON file');
    }
  };

  const qLower = query.trim().toLowerCase();
  const filteredChains = useMemo(() => {
    if (!prepTree) return [];
    return prepTree.chains.filter((ch) => {
      if (!qLower) return true;
      if (ch.chain.toLowerCase().includes(qLower)) return true;
      const allTexts = [
        ...ch.ligands.map((x) => `${x.resname} ${x.key} ${x.rcsb_name || ''}`),
        ...ch.cofactors.map((x) => `${x.resname} ${x.key}`),
        ...ch.additives.map((x) => `${x.resname} ${x.key}`),
        ...ch.metals_ions.map((x) => `${x.resname} ${x.key} ${x.element}`),
        ...ch.altlocs.map((x) => `${x.resname} ${x.key}`),
        ...ch.missing_gaps.map((x) => `${x.id} ${x.sequence_1}`),
        ...ch.modified.map((x) => `${x.resname} ${x.key}`),
        ...ch.unknown_atom_types.map((x) => `${x.resname} ${x.key} ${x.unsupported_elements.join(' ')}`),
      ];
      return allTexts.some((t) => t.toLowerCase().includes(qLower));
    });
  }, [prepTree, qLower]);

  if (!prepTree) {
    return (
      <div className="p-6 text-center text-[13px] text-[#86868B]">
        Preparation tree data is loading…
      </div>
    );
  }

  const s = prepTree.summary;

  return (
    <div className="space-y-4">
      {/* Header Banner */}
      <div className="apple-card p-4">
        <div className="flex flex-wrap items-start justify-between gap-2">
          <div>
            <div className="flex items-center gap-2">
              <span className="inline-flex items-center justify-center w-7 h-7 rounded-xl bg-[#007AFF]/12 text-[#007AFF] dark:text-[#0A84FF]">
                <Layers className="w-4 h-4" />
              </span>
              <h2 className="text-[15px] font-semibold tracking-tight">
                Receptor &amp; Docking Preparation Tree
              </h2>
            </div>
            <p className="text-[12px] text-[#6E6E73] dark:text-[#98989D] mt-1">
              ADapt-grade preparation: switch biological assemblies, rebuild individual missing gaps with geometry checks, customize ion formal charges &amp; AD4 substitutions, assign per-item Receptor/Ligand roles, and export MD5-verified PDBQT + provenance.
            </p>
          </div>
          <div className="flex items-center gap-1.5 text-[11px] font-mono text-[#6E6E73] dark:text-[#98989D]">
            {session.input_md5 && (
              <span className="px-2 py-0.5 rounded-md bg-black/[0.04] dark:bg-white/[0.06]" title={`Input MD5: ${session.input_md5}`}>
                MD5: {session.input_md5.slice(0, 8)}
              </span>
            )}
          </div>
        </div>

        {/* Summary KPI Chips — Unified Monochrome + Semantic Alert */}
        <div className="grid grid-cols-3 sm:grid-cols-6 gap-2 mt-3">
          {[
            { label: 'Chains', val: s.protein_chains, tone: 'text-[#1D1D1F] dark:text-white' },
            { label: 'DNA/RNA', val: s.nucleic_residues, tone: s.nucleic_residues > 0 ? 'text-[#007AFF] dark:text-[#0A84FF]' : 'text-[#86868B]' },
            { label: 'Gaps/Ends', val: s.internal_gaps + s.missing_termini, tone: s.internal_gaps > 0 ? 'text-[#B25000] dark:text-[#FF9F0A]' : 'text-[#86868B]' },
            { label: 'Metals/Ions', val: s.metals_ions, tone: s.metals_ions > 0 ? 'text-[#1D1D1F] dark:text-white' : 'text-[#86868B]' },
            { label: 'AD4 Alerts', val: s.unsupported_ad4_residues, tone: s.unsupported_ad4_residues > 0 ? 'text-[#D70015] dark:text-[#FF453A]' : 'text-[#1B8A38] dark:text-[#30D158]' },
            { label: 'Ligands', val: s.ligands + s.cofactors, tone: 'text-[#007AFF] dark:text-[#0A84FF]' },
          ].map((k) => (
            <div key={k.label} className="rounded-xl bg-black/[0.03] dark:bg-white/[0.04] px-2.5 py-1.5 text-center">
              <div className={`text-[15px] font-semibold tabular-nums ${k.tone}`}>{k.val}</div>
              <div className="text-[10px] font-medium text-[#86868B] uppercase tracking-wider">{k.label}</div>
            </div>
          ))}
        </div>
      </div>

      {/* 1. Biological Assembly Switcher */}
      <div className="apple-card p-4 space-y-2.5">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div>
            <div className="text-[12.5px] font-semibold">Biological Assembly (REMARK 350 / mmCIF)</div>
            <div className="text-[11px] text-[#86868B]">
              Choose between the deposited Asymmetric Unit (ASU) and crystallographic biological assemblies.
            </div>
          </div>
          <label className="flex items-center gap-1.5 text-[11.5px] text-[#6E6E73] dark:text-[#98989D] cursor-pointer select-none">
            <input
              type="checkbox"
              checked={dedupAxialIons}
              onChange={(e) => setDedupAxialIons(e.target.checked)}
              className="rounded accent-[#007AFF]"
            />
            Deduplicate axial ions
          </label>
        </div>

        <div className="flex flex-wrap gap-1.5">
          {assemblies.map((asm) => (
            <button
              key={asm.id}
              onClick={() => handleSwitchAssembly(asm.id)}
              title={asm.description}
              className={`px-3 py-1.5 rounded-xl text-[12px] font-medium transition flex items-center gap-1.5 border ${
                asm.is_active
                  ? 'bg-[#007AFF] text-white border-[#007AFF] shadow-sm'
                  : 'bg-black/[0.03] dark:bg-white/[0.05] text-[#1D1D1F] dark:text-[#E5E5EA] border-black/[0.07] dark:border-white/[0.08] hover:border-[#007AFF]/50'
              }`}
            >
              <span>{asm.label}</span>
              {asm.has_symmetry_copies && (
                <span
                  className={`text-[10px] px-1.5 py-0.2 rounded-md ${
                    asm.is_active ? 'bg-white/20 text-white' : 'pp-badge-neutral'
                  }`}
                >
                  ×{asm.operator_count} sym
                </span>
              )}
            </button>
          ))}
        </div>

        {(session.symmetry_copy_chains?.length || 0) > 0 && (
          <div className="text-[11.5px] rounded-xl bg-[#007AFF]/[0.08] text-[#007AFF] dark:text-[#0A84FF] px-3 py-2 flex items-center justify-between">
            <span>
              Symmetry-generated chains active: <strong>{session.symmetry_copy_chains?.join(', ')}</strong> (toggle <strong>Sym Tint</strong> in 3D viewer)
            </span>
          </div>
        )}
      </div>

      {/* 2. Global Options & RCSB Ligand Verification */}
      <div className="apple-card p-3.5 space-y-2.5">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div className="flex flex-wrap items-center gap-3">
            <label className="flex items-center gap-1.5 text-[12px] font-medium cursor-pointer select-none">
              <input
                type="checkbox"
                checked={dockingPrefs.include_nucleic}
                onChange={(e) => handleUpdatePrefs({ include_nucleic: e.target.checked })}
                className="rounded accent-[#007AFF]"
              />
              <Dna className="w-3.5 h-3.5 text-[#6E6E73] dark:text-[#98989D]" />
              Include DNA / RNA in Receptor
            </label>

            <label className="flex items-center gap-1.5 text-[12px] font-medium cursor-pointer select-none">
              <input
                type="checkbox"
                checked={strictGapGate}
                onChange={(e) => setStrictGapGate(e.target.checked)}
                className="rounded accent-[#007AFF]"
              />
              <ShieldCheck className="w-3.5 h-3.5 text-[#6E6E73] dark:text-[#98989D]" />
              Strict Loop Geometry Gate
            </label>
          </div>

          <div className="flex items-center gap-1.5">
            <button
              onClick={handleVerifyRcsb}
              className="px-2.5 py-1.5 rounded-xl text-[11.5px] font-medium bg-[#007AFF]/12 text-[#007AFF] dark:text-[#0A84FF] hover:bg-[#007AFF]/20 transition flex items-center gap-1.5"
              title="Query RCSB Data API (cached locally) for official ligand vs. additive classification"
            >
              <Database className="w-3.5 h-3.5" />
              {prepTree.rcsb_source === 'rcsb_data_api' ? `RCSB Verified (${prepTree.rcsb_pdb_id})` : 'Verify Ligands via RCSB'}
            </button>
            <button
              onClick={() => handleUpdatePrefs({ reset_all: true })}
              className="px-2.5 py-1.5 rounded-xl text-[11.5px] font-medium bg-black/[0.04] dark:bg-white/[0.06] text-[#6E6E73] dark:text-[#98989D] hover:text-[#1D1D1F] dark:hover:text-white transition flex items-center gap-1"
              title="Reset all custom roles, ion charges, and substitutions"
            >
              <RotateCcw className="w-3 h-3" />
              Reset
            </button>
          </div>
        </div>

        {/* Search & Filter Row */}
        <div className="flex flex-wrap items-center gap-2 pt-1">
          <div className="relative flex-1 min-w-[160px]">
            <Search className="w-3.5 h-3.5 absolute left-2.5 top-1/2 -translate-y-1/2 text-[#86868B]" />
            <input
              type="text"
              value={query}
              onChange={(e) => setQuery(e.target.value)}
              placeholder="Filter chain, ligand, ion, gap…"
              className="w-full pl-8 pr-3 py-1.5 rounded-xl text-[12px] bg-black/[0.03] dark:bg-white/[0.05] border border-black/[0.06] dark:border-white/[0.08] focus:outline-none focus:border-[#007AFF]"
            />
          </div>
          <div className="flex items-center gap-1 bg-black/[0.04] dark:bg-white/[0.06] p-0.5 rounded-xl">
            {(
              [
                ['all', 'All'],
                ['issues', 'Issues & AD4'],
                ['heterogens', 'Ligands & Ions'],
                ['gaps', 'Gaps & AltLocs'],
              ] as Array<[TreeFilter, string]>
            ).map(([id, label]) => (
              <button
                key={id}
                onClick={() => setFilter(id)}
                className={`px-2.5 py-1 rounded-[9px] text-[11px] font-medium transition ${
                  filter === id
                    ? 'bg-white dark:bg-[#2C2C2E] text-[#1D1D1F] dark:text-white shadow-sm'
                    : 'text-[#6E6E73] dark:text-[#98989D]'
                }`}
              >
                {label}
              </button>
            ))}
          </div>
        </div>
      </div>

      {/* 3. Hierarchical Per-Chain Tree */}
      <div className="space-y-2.5">
        {filteredChains.map((ch, idx) => {
          const open = isChainOpen(ch.chain, idx);
          const showGaps = filter === 'all' || filter === 'gaps' || filter === 'issues';
          const showAltlocs = filter === 'all' || filter === 'gaps' || filter === 'issues';
          const showScMod = filter === 'all' || filter === 'issues';
          const showUnknown = filter === 'all' || filter === 'issues' || filter === 'heterogens';
          const showIons = filter === 'all' || filter === 'heterogens' || filter === 'issues';
          const showHet = filter === 'all' || filter === 'heterogens';

          return (
            <div key={ch.chain} className="apple-card overflow-hidden">
              {/* Chain Row Header */}
              <div className="px-3.5 py-2.5 flex items-center justify-between gap-2 bg-black/[0.02] dark:bg-white/[0.03] border-b border-black/[0.05] dark:border-white/[0.06]">
                <button
                  onClick={() => toggleChain(ch.chain, idx)}
                  className="flex items-center gap-2 text-left flex-1 min-w-0"
                >
                  {open ? <ChevronDown className="w-4 h-4 text-[#86868B]" /> : <ChevronRight className="w-4 h-4 text-[#86868B]" />}
                  <span className="text-[13px] font-semibold">Chain {ch.chain}</span>
                  {ch.is_symmetry_copy && (
                    <span className="px-1.5 py-0.5 rounded-md text-[10px] font-semibold pp-badge-accent">
                      Sym Copy
                    </span>
                  )}
                  <span className="text-[11px] text-[#86868B] truncate">
                    {ch.protein_count} aa
                    {ch.nucleic_count > 0 ? ` · ${ch.nucleic_count} nt` : ''}
                    {ch.ligands.length > 0 ? ` · ${ch.ligands.length} lig` : ''}
                    {ch.metals_ions.length > 0 ? ` · ${ch.metals_ions.length} ion` : ''}
                    {ch.waters.length > 0 ? ` · ${ch.waters.length} H₂O` : ''}
                  </span>
                </button>

                <label className="flex items-center gap-1.5 text-[11.5px] font-medium cursor-pointer shrink-0">
                  <input
                    type="checkbox"
                    checked={ch.role === 'receptor'}
                    onChange={(e) =>
                      handleUpdatePrefs({
                        chain_roles: { [ch.chain]: e.target.checked ? 'receptor' : 'exclude' },
                      })
                    }
                    className="rounded accent-[#007AFF]"
                  />
                  <span>Receptor</span>
                </label>
              </div>

              {open && (
                <div className="p-3 space-y-3 text-[12px]">
                  {/* Unsupported AutoDock 4 Atom Types */}
                  {showUnknown && ch.unknown_atom_types.length > 0 && (
                    <div className="rounded-xl border border-[#FF3B30]/30 bg-[#FF3B30]/[0.06] p-2.5 space-y-2">
                      <div className="flex items-center gap-1.5 text-[11.5px] font-semibold text-[#D70015] dark:text-[#FF453A]">
                        <AlertTriangle className="w-3.5 h-3.5" />
                        Unsupported AutoDock 4 Atom Types ({ch.unknown_atom_types.length})
                      </div>
                      <div className="space-y-1.5">
                        {ch.unknown_atom_types.map((u) => (
                          <div key={u.key} className="flex flex-wrap items-center justify-between gap-2 bg-white/70 dark:bg-black/25 rounded-lg px-2.5 py-1.5">
                            <button
                              onClick={() => onSelectResidue({ chain: u.chain, resseq: u.resseq, icode: u.icode, resname: u.resname, key: u.key }, true)}
                              className="text-left hover:underline font-medium"
                            >
                              {u.resname} <span className="font-mono text-[#6E6E73]">{u.key}</span>{' '}
                              <span className="text-[10.5px] text-[#FF3B30] font-mono">[{u.unsupported_elements.join(', ')}]</span>
                            </button>
                            {u.is_ion ? (
                              <select
                                value={u.substituted_element || ''}
                                onChange={(e) =>
                                  handleUpdatePrefs({
                                    ion_substitutions: { [u.key]: e.target.value || null },
                                    ...(e.target.value ? { residue_roles: { [u.key]: 'receptor' } } : {}),
                                  })
                                }
                                className="rounded-lg px-2 py-1 text-[11px] bg-white dark:bg-[#2C2C2E] border border-black/10 dark:border-white/15"
                              >
                                <option value="">Exclude (untyped)</option>
                                {ad4Subs.map((sub) => (
                                  <option key={sub.element} value={sub.element}>
                                    Substitute → {sub.label}
                                  </option>
                                ))}
                              </select>
                            ) : (
                              <span className="text-[10.5px] text-[#86868B]">Auto-excluded from .pdbqt</span>
                            )}
                          </div>
                        ))}
                      </div>
                    </div>
                  )}

                  {/* Missing Residues (REMARK 465 Internal Gaps & Termini) */}
                  {showGaps && (ch.missing_gaps.length > 0 || ch.missing_termini.length > 0) && (
                    <div className="space-y-1.5">
                      <div className="text-[11px] font-semibold uppercase tracking-wider text-[#FF9500]">
                        Missing Residues (REMARK 465 / SEQRES)
                      </div>
                      {ch.missing_gaps.map((g) => (
                        <div
                          key={g.id}
                          className="rounded-xl border border-[#FF9500]/25 bg-[#FF9500]/[0.05] px-3 py-2 space-y-1.5"
                        >
                          <div className="flex flex-wrap items-center justify-between gap-2">
                            <button
                              onClick={() =>
                                onSelectResidue(
                                  { chain: g.chain, resseq: g.start_resseq, icode: '', resname: g.start_resname },
                                  true
                                )
                              }
                              className="text-left hover:underline"
                            >
                              <span className="font-semibold text-[#B26200] dark:text-[#FF9F0A]">Gap:</span>{' '}
                              <span className="font-mono text-[11.5px]">{g.label}</span>{' '}
                              <span className="text-[11px] text-[#86868B]">
                                ({g.missing_count} res · {g.gap_distance_A} Å · {g.source})
                              </span>
                            </button>
                            <div className="flex items-center gap-1.5">
                              {g.can_rebuild && (
                                <button
                                  onClick={() =>
                                    handleRebuildGap(g.chain, g.start_resseq, g.end_resseq, g.sequence_3, strictGapGate)
                                  }
                                  className="px-2.5 py-1 rounded-lg text-[11px] font-medium bg-[#FF9500] text-white hover:bg-[#E08600] transition flex items-center gap-1"
                                >
                                  <Wrench className="w-3 h-3" />
                                  Rebuild Gap
                                </button>
                              )}
                              {g.status === 'refused' && strictGapGate && (
                                <button
                                  onClick={() =>
                                    handleRebuildGap(g.chain, g.start_resseq, g.end_resseq, g.sequence_3, false)
                                  }
                                  className="px-2 py-1 rounded-lg text-[10.5px] font-medium bg-black/[0.06] dark:bg-white/[0.08] text-[#6E6E73] dark:text-[#AEAEB2]"
                                  title="Override geometry gate and insert best seed anyway"
                                >
                                  Force
                                </button>
                              )}
                            </div>
                          </div>
                          {g.refusal_reason && (
                            <div className="text-[11px] text-[#D70015] dark:text-[#FF453A]">
                              ⚠ {g.refusal_reason}
                            </div>
                          )}
                        </div>
                      ))}

                      {ch.missing_termini.map((t) => (
                        <div
                          key={t.id}
                          className="rounded-xl bg-black/[0.025] dark:bg-white/[0.04] px-3 py-1.5 flex items-center justify-between gap-2"
                        >
                          <button
                            onClick={() =>
                              onSelectResidue(
                                { chain: t.chain, resseq: t.anchor_resseq, icode: '', resname: t.anchor_resname },
                                true
                              )
                            }
                            className="text-left hover:underline"
                          >
                            <span className="font-semibold text-[#86868B]">{t.terminus}-term:</span>{' '}
                            <span className="font-mono text-[11px]">{t.sequence_1.slice(0, 18)}{t.sequence_1.length > 18 ? '…' : ''}</span>{' '}
                            <span className="text-[10.5px] text-[#86868B]">({t.missing_count} res · {t.range_label})</span>
                          </button>
                          <span className="text-[10px] text-[#86868B]">Info only (3D marker)</span>
                        </div>
                      ))}
                    </div>
                  )}

                  {/* Rebuilt Loops Summary for this chain */}
                  {showGaps && rebuiltGaps.filter((rg) => rg.chain === ch.chain).length > 0 && (
                    <div className="space-y-1">
                      <div className="text-[11px] font-semibold uppercase tracking-wider text-[#6E6E73] dark:text-[#98989D]">
                        Rebuilt Loops (Highlighted in 3D)
                      </div>
                      {rebuiltGaps
                        .filter((rg) => rg.chain === ch.chain)
                        .map((rg) => (
                          <div
                            key={rg.gap_id}
                            className="rounded-xl bg-[#007AFF]/[0.08] text-[#007AFF] dark:text-[#0A84FF] px-3 py-1.5 flex flex-wrap items-center justify-between gap-2 text-[11px]"
                          >
                            <span>
                              <strong>{rg.gap_id}</strong> ({rg.sequence_1}) · Seed {rg.seed}
                            </span>
                            {rg.geometry && (
                              <span className="font-mono text-[10.5px]">
                                C–N {rg.geometry.cn_min_A.toFixed(2)}–{rg.geometry.cn_max_A.toFixed(2)} Å · Contact {rg.geometry.min_contact_A.toFixed(2)} Å
                              </span>
                            )}
                          </div>
                        ))}
                    </div>
                  )}

                  {/* Alternate Locations */}
                  {showAltlocs && ch.altlocs.length > 0 && (
                    <div className="space-y-1.5">
                      <div className="text-[11px] font-semibold uppercase tracking-wider text-[#007AFF]">
                        Alternate Locations ({ch.altlocs.length})
                      </div>
                      <div className="space-y-1">
                        {ch.altlocs.map((al) => (
                          <div
                            key={al.key}
                            className="rounded-xl bg-black/[0.025] dark:bg-white/[0.04] px-3 py-1.5 flex flex-wrap items-center justify-between gap-2"
                          >
                            <button
                              onClick={() =>
                                onSelectResidue({ chain: al.chain, resseq: al.resseq, icode: al.icode, resname: al.resname, key: al.key }, true)
                              }
                              className="font-medium hover:underline"
                            >
                              {al.resname} <span className="font-mono text-[#6E6E73]">{al.key}</span>
                            </button>
                            <div className="flex items-center gap-1">
                              {al.conformers.map((cf) => (
                                <button
                                  key={cf.id}
                                  onClick={() => handlePickAltloc(al.chain, al.resseq, al.icode, cf.id)}
                                  className="px-2 py-0.5 rounded-lg text-[11px] font-medium bg-white dark:bg-[#2C2C2E] border border-black/10 dark:border-white/15 hover:border-[#007AFF] hover:text-[#007AFF] transition"
                                  title={`Pick conformer ${cf.id} (occupancy ${Math.round(cf.occ * 100)}%, B=${cf.b_factor})`}
                                >
                                  {cf.id} · {Math.round(cf.occ * 100)}%
                                </button>
                              ))}
                            </div>
                          </div>
                        ))}
                      </div>
                    </div>
                  )}

                  {/* Incomplete Side Chains & Modified Residues */}
                  {showScMod && (ch.incomplete_sidechains.length > 0 || ch.modified.length > 0) && (
                    <div className="space-y-1.5">
                      <div className="flex items-center justify-between">
                        <span className="text-[11px] font-semibold uppercase tracking-wider text-[#6E6E73] dark:text-[#98989D]">
                          Side Chains &amp; Modified ({ch.incomplete_sidechains.length + ch.modified.length})
                        </span>
                        <button
                          onClick={() => onJumpTab('pdbfixer')}
                          className="text-[11px] text-[#007AFF] hover:underline"
                        >
                          Open PDBFixer →
                        </button>
                      </div>
                      <div className="flex flex-wrap gap-1.5">
                        {ch.incomplete_sidechains.map((sc) => (
                          <button
                            key={`sc-${sc.key}`}
                            onClick={() => {
                              onSelectResidue({ chain: sc.chain, resseq: sc.resseq, icode: sc.icode, resname: sc.resname, key: sc.key }, true);
                              onJumpTab('coot');
                            }}
                            className="px-2 py-1 rounded-lg text-[11px] pp-badge-warn hover:opacity-85"
                            title={`Missing atoms: ${sc.missing.join(', ')} (click to edit in COOT)`}
                          >
                            {sc.resname} {sc.key} (-{sc.missing.length} atoms)
                          </button>
                        ))}
                        {ch.modified.map((md) => (
                          <button
                            key={`mod-${md.key}`}
                            onClick={() =>
                              onSelectResidue({ chain: md.chain, resseq: md.resseq, icode: md.icode, resname: md.resname, key: md.key }, true)
                            }
                            className="px-2 py-1 rounded-lg text-[11px] pp-badge-neutral hover:opacity-85"
                          >
                            {md.resname} {md.key} → {md.replacement}
                          </button>
                        ))}
                      </div>
                    </div>
                  )}

                  {/* Metals & Manual Ion Charges */}
                  {showIons && ch.metals_ions.length > 0 && (
                    <div className="space-y-1.5">
                      <div className="flex items-center gap-1.5 text-[11px] font-semibold uppercase tracking-wider text-[#6E6E73] dark:text-[#98989D]">
                        <Atom className="w-3.5 h-3.5 text-[#007AFF]" />
                        Metals &amp; Manual Ion Charges ({ch.metals_ions.length})
                      </div>
                      <div className="space-y-1.5">
                        {ch.metals_ions.map((ion) => {
                          const draft = ionChargeDrafts[ion.key] ?? ion.charge.toFixed(1);
                          return (
                            <div
                              key={ion.key}
                              className="rounded-xl bg-black/[0.025] dark:bg-white/[0.04] px-3 py-2 flex flex-wrap items-center justify-between gap-2"
                            >
                              <button
                                onClick={() =>
                                  onSelectResidue({ chain: ion.chain, resseq: ion.resseq, icode: ion.icode, resname: ion.resname, key: ion.key }, true)
                                }
                                className="text-left hover:underline font-medium flex items-center gap-1.5"
                              >
                                <span>{ion.resname}</span>
                                <span className="font-mono text-[#6E6E73]">{ion.key}</span>
                                {ion.substituted_element && (
                                  <span className="px-1.5 py-0.5 rounded-md text-[10px] bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158]">
                                    → {ion.substituted_element}
                                  </span>
                                )}
                                {!ion.ad4_supported && (
                                  <span className="px-1.5 py-0.5 rounded-md text-[10px] bg-[#FF3B30]/15 text-[#D70015] dark:text-[#FF453A]">
                                    AD4 untyped
                                  </span>
                                )}
                              </button>

                              <div className="flex flex-wrap items-center gap-1.5">
                                {/* Editable Formal Charge */}
                                <div className="flex items-center gap-1 bg-white dark:bg-[#2C2C2E] rounded-lg border border-black/10 dark:border-white/15 px-1.5 py-0.5">
                                  <span className="text-[10px] text-[#86868B]">q=</span>
                                  <input
                                    type="number"
                                    step="0.5"
                                    min="-4"
                                    max="6"
                                    value={draft}
                                    onChange={(e) =>
                                      setIonChargeDrafts((p) => ({ ...p, [ion.key]: e.target.value }))
                                    }
                                    onBlur={() => {
                                      const parsed = parseFloat(draft);
                                      if (!Number.isNaN(parsed) && Math.abs(parsed - ion.charge) > 1e-4) {
                                        handleUpdatePrefs({ ion_charges: { [ion.key]: parsed } });
                                      }
                                    }}
                                    className="w-12 text-[11px] font-mono bg-transparent focus:outline-none text-right"
                                    title={`CCD default formal charge: ${ion.default_charge >= 0 ? '+' : ''}${ion.default_charge}`}
                                  />
                                </div>

                                {/* Ion Substitution Picker */}
                                <select
                                  value={ion.substituted_element || 'ORIGINAL'}
                                  onChange={(e) => {
                                    const val = e.target.value;
                                    handleUpdatePrefs({
                                      ion_substitutions: { [ion.key]: val === 'ORIGINAL' ? null : val },
                                      ...(val !== 'ORIGINAL' ? { residue_roles: { [ion.key]: 'receptor' } } : {}),
                                    });
                                  }}
                                  className="rounded-lg px-1.5 py-1 text-[10.5px] bg-white dark:bg-[#2C2C2E] border border-black/10 dark:border-white/15"
                                  title="Substitute ion with an AutoDock 4-supported metal ion"
                                >
                                  <option value="ORIGINAL">Keep {ion.element}</option>
                                  {ad4Subs.map((sub) => (
                                    <option key={sub.element} value={sub.element}>
                                      Sub → {sub.element}
                                    </option>
                                  ))}
                                </select>

                                {/* Receptor / Exclude Role */}
                                <div className="flex rounded-lg bg-black/[0.05] dark:bg-white/[0.08] p-0.5">
                                  {(['receptor', 'exclude'] as const).map((role) => (
                                    <button
                                      key={role}
                                      onClick={() => handleUpdatePrefs({ residue_roles: { [ion.key]: role } })}
                                      className={`px-2 py-0.5 rounded-md text-[10.5px] font-medium capitalize transition ${
                                        ion.role === role
                                          ? role === 'receptor'
                                            ? 'bg-[#007AFF] text-white'
                                            : 'bg-[#6E6E73] text-white'
                                          : 'text-[#6E6E73] dark:text-[#98989D]'
                                      }`}
                                    >
                                      {role === 'receptor' ? 'Rec' : 'Excl'}
                                    </button>
                                  ))}
                                </div>
                              </div>
                            </div>
                          );
                        })}
                      </div>
                    </div>
                  )}

                  {/* Ligands, Cofactors & Additives (Per-Item Receptor | Ligand | Exclude) */}
                  {showHet && (ch.ligands.length > 0 || ch.cofactors.length > 0 || ch.additives.length > 0) && (
                    <div className="space-y-2">
                      {[
                        { title: 'Ligands', items: ch.ligands, accent: 'text-[#6E6E73] dark:text-[#98989D]' },
                        { title: 'Cofactors', items: ch.cofactors, accent: 'text-[#6E6E73] dark:text-[#98989D]' },
                        { title: 'Additives / Buffers', items: ch.additives, accent: 'text-[#86868B]' },
                      ]
                        .filter((grp) => grp.items.length > 0)
                        .map((grp) => (
                          <div key={grp.title} className="space-y-1">
                            <div className={`text-[11px] font-semibold uppercase tracking-wider ${grp.accent}`}>
                              {grp.title} ({grp.items.length})
                            </div>
                            {grp.items.map((item) => (
                              <div
                                key={item.key}
                                className="rounded-xl bg-black/[0.025] dark:bg-white/[0.04] px-3 py-1.5 flex flex-wrap items-center justify-between gap-2"
                              >
                                <button
                                  onClick={() =>
                                    onSelectResidue(
                                      { chain: item.chain, resseq: item.resseq, icode: item.icode, resname: item.resname, key: item.key },
                                      true
                                    )
                                  }
                                  className="text-left hover:underline min-w-0 flex-1"
                                >
                                  <span className="font-semibold">{item.resname}</span>{' '}
                                  <span className="font-mono text-[#6E6E73]">{item.key}</span>{' '}
                                  <span className="text-[10.5px] text-[#86868B]">({item.heavy_atoms} heavy)</span>
                                  {item.rcsb_name && (
                                    <div className="text-[10.5px] text-[#6E6E73] dark:text-[#98989D] truncate max-w-[230px]">
                                      {item.rcsb_name}
                                    </div>
                                  )}
                                </button>

                                <div className="flex items-center gap-1.5 shrink-0">
                                  {item.role === 'ligand' && (
                                    <button
                                      onClick={() => handleExportLigandPdbqt(item.key, item.resname)}
                                      className="px-2 py-0.5 rounded-lg text-[10.5px] font-medium pp-badge-accent hover:opacity-85"
                                      title="Download this ligand prepared as AutoDock .pdbqt"
                                    >
                                      .pdbqt
                                    </button>
                                  )}
                                  <div className="flex rounded-lg bg-black/[0.05] dark:bg-white/[0.08] p-0.5">
                                    {(['receptor', 'ligand', 'exclude'] as const).map((role) => (
                                      <button
                                        key={role}
                                        onClick={() => handleUpdatePrefs({ residue_roles: { [item.key]: role } })}
                                        className={`px-2 py-0.5 rounded-md text-[10.5px] font-medium transition ${
                                          item.role === role
                                            ? role === 'receptor'
                                              ? 'bg-[#007AFF] text-white'
                                              : role === 'ligand'
                                              ? 'bg-[#1D1D1F] dark:bg-white text-white dark:text-[#1D1D1F]'
                                              : 'bg-[#6E6E73] text-white'
                                            : 'text-[#6E6E73] dark:text-[#98989D]'
                                        }`}
                                      >
                                        {role === 'receptor' ? 'Rec' : role === 'ligand' ? 'Lig' : 'Excl'}
                                      </button>
                                    ))}
                                  </div>
                                </div>
                              </div>
                            ))}
                          </div>
                        ))}
                    </div>
                  )}

                  {/* Waters (Expandable Per-Water Receptor Toggle) */}
                  {showHet && ch.waters.length > 0 && (
                    <div className="space-y-1">
                      <div className="flex items-center justify-between">
                        <button
                          onClick={() => setExpandedWaters((p) => ({ ...p, [ch.chain]: !p[ch.chain] }))}
                          className="text-[11px] font-semibold uppercase tracking-wider text-[#6E6E73] dark:text-[#98989D] flex items-center gap-1"
                        >
                          {expandedWaters[ch.chain] ? <ChevronDown className="w-3.5 h-3.5" /> : <ChevronRight className="w-3.5 h-3.5" />}
                          Waters ({ch.waters.length} · {ch.waters.filter((w) => w.role === 'receptor').length} kept in receptor)
                        </button>
                      </div>
                      {expandedWaters[ch.chain] && (
                        <div className="max-h-36 overflow-y-auto grid grid-cols-2 gap-1 pt-1">
                          {ch.waters.slice(0, 80).map((w) => (
                            <div
                              key={w.key}
                              className="flex items-center justify-between rounded-lg bg-black/[0.025] dark:bg-white/[0.04] px-2 py-1 text-[11px]"
                            >
                              <button
                                onClick={() =>
                                  onSelectResidue({ chain: w.chain, resseq: w.resseq, icode: w.icode, resname: w.resname, key: w.key }, true)
                                }
                                className="font-mono hover:underline"
                              >
                                {w.resname} {w.key}
                              </button>
                              <input
                                type="checkbox"
                                checked={w.role === 'receptor'}
                                onChange={(e) =>
                                  handleUpdatePrefs({
                                    residue_roles: { [w.key]: e.target.checked ? 'receptor' : 'exclude' },
                                  })
                                }
                                className="rounded accent-[#007AFF]"
                                title="Include this structural water in receptor PDBQT"
                              />
                            </div>
                          ))}
                        </div>
                      )}
                    </div>
                  )}

                  {/* Empty state when chain is completely clean */}
                  {ch.unknown_atom_types.length === 0 &&
                    ch.missing_gaps.length === 0 &&
                    ch.missing_termini.length === 0 &&
                    ch.altlocs.length === 0 &&
                    ch.incomplete_sidechains.length === 0 &&
                    ch.modified.length === 0 &&
                    ch.metals_ions.length === 0 &&
                    ch.ligands.length === 0 &&
                    ch.cofactors.length === 0 &&
                    ch.additives.length === 0 &&
                    ch.waters.length === 0 && (
                      <div className="text-[11.5px] text-[#34C759] flex items-center gap-1.5">
                        <CheckCircle2 className="w-3.5 h-3.5" />
                        Clean polymer chain — ready for receptor PDBQT export.
                      </div>
                    )}
                </div>
              )}
            </div>
          );
        })}
      </div>

      {/* 4. Docking Outputs & Cryptographic Provenance (REMARK 850 + MD5 session.json) */}
      <div className="apple-card p-4 space-y-3">
        <div className="flex items-center justify-between">
          <div>
            <div className="text-[13px] font-semibold flex items-center gap-1.5">
              <Sparkles className="w-4 h-4 text-[#007AFF]" />
              Docking Outputs &amp; Cryptographic Provenance
            </div>
            <div className="text-[11px] text-[#86868B]">
              Vina 1.2 column-aligned `.pdbqt`, all-atom `_rec_H.pdb` with PROPKA `USER MOD` log, `REMARK 850` headers, hex serials &gt;99,999, and MD5 `session.json`.
            </div>
          </div>
          <button
            onClick={handleInspectProvenance}
            className="text-[11.5px] font-medium text-[#007AFF] hover:underline flex items-center gap-1"
          >
            <RefreshCw className="w-3 h-3" />
            {provPreview ? 'Hide JSON' : 'Inspect Provenance'}
          </button>
        </div>

        <div className="grid grid-cols-2 gap-2">
          <button
            onClick={() => handleExport('pdbqt', '_receptor.pdbqt')}
            className="apple-btn-primary !py-2.5 text-[12px]"
          >
            <Download className="w-3.5 h-3.5" />
            Receptor .pdbqt
          </button>
          <button
            onClick={() => handleExport('rec_h_pdb', '_rec_H.pdb')}
            className="apple-btn-secondary !py-2.5 text-[12px] border border-black/[0.08] dark:border-white/[0.1]"
          >
            <Download className="w-3.5 h-3.5 text-[#007AFF]" />
            Protonated _rec_H.pdb
          </button>
          <button
            onClick={() => handleExport('provenance_json', '_provenance_session.json')}
            className="apple-btn-secondary !py-2 text-[12px]"
          >
            <FileJson className="w-3.5 h-3.5 text-[#6E6E73] dark:text-[#98989D]" />
            Download session.json
          </button>
          <button
            onClick={() => verifyFileRef.current?.click()}
            className="apple-btn-secondary !py-2 text-[12px]"
          >
            <FileCheck2 className="w-3.5 h-3.5 text-[#6E6E73] dark:text-[#98989D]" />
            Verify session.json
          </button>
          <input
            ref={verifyFileRef}
            type="file"
            accept=".json"
            className="hidden"
            onChange={(e) => {
              const f = e.target.files?.[0];
              if (f) handleVerifyProvenanceFile(f);
              e.target.value = '';
            }}
          />
        </div>

        {verifyResult && (
          <div
            className={`rounded-xl px-3 py-2 text-[11.5px] border ${
              verifyResult.verified
                ? 'bg-[#34C759]/10 border-[#34C759]/30 text-[#1B8A38] dark:text-[#30D158]'
                : 'bg-[#FF3B30]/10 border-[#FF3B30]/30 text-[#D70015] dark:text-[#FF453A]'
            }`}
          >
            <div className="font-semibold">
              {verifyResult.verified ? '✓ Cryptographic MD5 Verification Passed' : '✕ Provenance Checksum Mismatch'}
            </div>
            <div className="text-[10.5px] font-mono mt-0.5 opacity-90">
              Input MD5: {verifyResult.input_md5_match ? 'MATCH' : 'MISMATCH'} · Structure MD5:{' '}
              {verifyResult.current_structure_md5_match ? 'MATCH' : 'MISMATCH'} · Assembly:{' '}
              {verifyResult.assembly_match ? 'MATCH' : 'DIFFERS'}
            </div>
          </div>
        )}

        {provPreview && (
          <pre className="max-h-56 overflow-auto rounded-xl bg-black/[0.04] dark:bg-black/40 p-3 text-[10.5px] font-mono text-[#1D1D1F] dark:text-[#E5E5EA]">
            {JSON.stringify(provPreview, null, 2)}
          </pre>
        )}
      </div>
    </div>
  );
};
