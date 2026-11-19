import React, { useEffect, useMemo, useState } from 'react';
import {
  ArrowDown,
  ArrowUp,
  Check,
  Code2,
  Eye,
  FileEdit,
  GitMerge,
  Hash,
  Layers,
  RefreshCw,
  Tag,
  Trash2,
} from 'lucide-react';
import { api } from '../api/client';
import type { StructureSessionState } from '../types/pdbprep';
import type { ChainHighlight } from './Viewer3D';

interface RenameChainsPanelProps {
  session: StructureSessionState;
  busy: boolean;
  highlightedChains: ChainHighlight | null;
  onHighlightChains: (hl: ChainHighlight | null, zoom?: boolean) => void;
  onStateUpdate: (state: StructureSessionState, toast: string) => void;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
}

export const RenameChainsPanel: React.FC<RenameChainsPanelProps> = ({
  session,
  busy,
  highlightedChains,
  onHighlightChains,
  onStateUpdate,
  onSetBusy,
  onError,
}) => {
  const chains = session.metadata.chains;

  // 1. File & Title
  const [filename, setFilename] = useState(session.name);
  const [title, setTitle] = useState(session.metadata.title);

  // 2. Chain Rename Mapping & Multi-Select Delete
  const [chainMap, setChainMap] = useState<Record<string, string>>({});
  const [selectedForDelete, setSelectedForDelete] = useState<string[]>([]);

  // 3. Merge Chains
  const [mergeSrc, setMergeSrc] = useState(chains[1]?.id || chains[0]?.id || '');
  const [mergeTgt, setMergeTgt] = useState(chains[0]?.id || '');
  const [mergeRenumber, setMergeRenumber] = useState(true);

  // 4. Per-chain Renumbering
  const [renumChain, setRenumChain] = useState(chains[0]?.id || '');
  const [renumMode, setRenumMode] = useState<'start' | 'offset'>('start');
  const [startNum, setStartNum] = useState(1);
  const [offsetVal, setOffsetVal] = useState(0);
  const [includeHet, setIncludeHet] = useState(false);

  // 5. Residue / Ligand Code Rename
  const detectedHets = useMemo(() => {
    const list: Array<{ label: string; resname: string; chain: string; resseq: number }> = [];
    for (const h of session.pdbfixer.heterogens) {
      list.push({
        label: `${h.resname} (${h.chain}:${h.resseq}) — ${h.category}`,
        resname: h.resname,
        chain: h.chain,
        resseq: h.resseq,
      });
    }
    return list;
  }, [session.pdbfixer.heterogens]);

  const [oldResCode, setOldResCode] = useState(detectedHets[0]?.resname || 'BEN');
  const [newResCode, setNewResCode] = useState('LIG');
  const [scopeChain, setScopeChain] = useState<string>('');
  const [moveToChain, setMoveToChain] = useState<string>('');

  // 6. Field Normalization
  const [waterChain, setWaterChain] = useState('W');
  const [uniformB, setUniformB] = useState('20.0');

  // 7. Live PDB preview filter
  const [pdbFilter, setPdbFilter] = useState<'head' | 'hetatm' | 'all_chains'>('head');

  useEffect(() => {
    setFilename(session.name);
    setTitle(session.metadata.title);
    const initMap: Record<string, string> = {};
    for (const c of session.metadata.chains) {
      initMap[c.id] = c.id;
    }
    setChainMap(initMap);
    setSelectedForDelete([]);
    if (session.metadata.chains.length > 0) {
      setRenumChain((prev) =>
        session.metadata.chains.some((c) => c.id === prev) ? prev : session.metadata.chains[0].id
      );
      setMergeTgt((prev) =>
        session.metadata.chains.some((c) => c.id === prev) ? prev : session.metadata.chains[0].id
      );
      setMergeSrc((prev) =>
        session.metadata.chains.some((c) => c.id === prev)
          ? prev
          : session.metadata.chains[1]?.id || session.metadata.chains[0].id
      );
    }
  }, [session.session_id, session.name, session.metadata.title, session.metadata.chains]);

  const toggleChainSpotlight = (cid: string) => {
    const alreadySingle =
      highlightedChains?.tone === 'select' &&
      highlightedChains.ids.length === 1 &&
      highlightedChains.ids[0] === cid;
    if (alreadySingle) {
      onHighlightChains(null, true);
    } else {
      onHighlightChains({ ids: [cid], tone: 'select' }, true);
    }
  };

  const toggleDeleteCheck = (cid: string) => {
    const next = selectedForDelete.includes(cid)
      ? selectedForDelete.filter((x) => x !== cid)
      : [...selectedForDelete, cid];
    setSelectedForDelete(next);
    if (next.length > 0) {
      onHighlightChains({ ids: next, tone: 'delete' }, true);
    } else {
      onHighlightChains(null, true);
    }
  };

  const runOp = async (
    label: string,
    fn: () => Promise<{ result: { message?: string }; state: StructureSessionState }>
  ) => {
    onSetBusy(label);
    try {
      const resp = await fn();
      onHighlightChains(null, false);
      onStateUpdate(resp.state, resp.result?.message || label);
    } catch (e: any) {
      onError(e?.message || 'Rename operation failed');
    } finally {
      onSetBusy(null);
    }
  };

  const hasChainChanges = Object.entries(chainMap).some(([oldId, newId]) => newId.trim() && oldId !== newId.trim());

  const moveChainOrder = (index: number, dir: -1 | 1) => {
    const nextIdx = index + dir;
    if (nextIdx < 0 || nextIdx >= chains.length) return;
    const order = chains.map((c) => c.id);
    [order[index], order[nextIdx]] = [order[nextIdx], order[index]];
    runOp('Reordering chains…', () => api.reorderChains(session.session_id, order));
  };

  const pdbPreviewText = useMemo(() => {
    const lines = session.pdb_text.split('\n');
    if (pdbFilter === 'hetatm') {
      const filtered = lines.filter((l) => l.startsWith('HETATM') || l.startsWith('CONECT') || l.startsWith('SSBOND'));
      return filtered.slice(0, 60).join('\n') || '(No HETATM records)';
    }
    if (pdbFilter === 'all_chains') {
      // Show first 3 lines of each chain + TER lines
      const out: string[] = [];
      const seenPerChain: Record<string, number> = {};
      for (const l of lines) {
        if (l.startsWith('TITLE') || l.startsWith('HEADER') || l.startsWith('SSBOND') || l.startsWith('TER')) {
          out.push(l);
        } else if (l.startsWith('ATOM') || l.startsWith('HETATM')) {
          const ch = l.slice(21, 22).trim() || '_';
          seenPerChain[ch] = (seenPerChain[ch] || 0) + 1;
          if (seenPerChain[ch] <= 4) out.push(l);
        }
      }
      return out.slice(0, 65).join('\n');
    }
    return lines.slice(0, 45).join('\n');
  }, [session.pdb_text, pdbFilter]);

  return (
    <div className="space-y-4 pb-6">
      {/* Header */}
      <div className="apple-card p-4">
        <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#007AFF] dark:text-[#0A84FF]">
          Structure Organizer · Live Coordinate Editor
        </div>
        <h2 className="text-[16px] font-semibold tracking-[-0.015em]">
          PDB File, Chain IDs, Residue Numbering & Ligand Renaming
        </h2>
        <p className="text-[12px] text-[#86868B] mt-0.5">
          Every change applies immediately to the live protein in memory, updates the 3D viewer, and supports COOT Undo/Redo.
        </p>
      </div>

      {/* 1. PDB File & Header Renaming */}
      <div className="apple-card p-4 space-y-3">
        <div className="flex items-center gap-2 text-[13.5px] font-semibold">
          <FileEdit className="w-4 h-4 text-[#007AFF]" />
          <span>PDB File Name & Header Title</span>
        </div>
        <div className="grid grid-cols-1 sm:grid-cols-12 gap-2.5">
          <div className="sm:col-span-5">
            <label className="text-[11px] text-[#86868B] block mb-1">Structure Filename (used for all exports)</label>
            <input
              value={filename}
              onChange={(e) => setFilename(e.target.value)}
              placeholder="e.g. 6LU7_prepared.pdb"
              className="apple-input font-mono"
            />
          </div>
          <div className="sm:col-span-5">
            <label className="text-[11px] text-[#86868B] block mb-1">PDB TITLE Header</label>
            <input
              value={title}
              onChange={(e) => setTitle(e.target.value)}
              placeholder="Structure description"
              className="apple-input"
            />
          </div>
          <div className="sm:col-span-2 flex items-end">
            <button
              disabled={busy || !filename.trim()}
              onClick={() =>
                runOp(`Renaming file to ${filename}…`, () =>
                  api.renameFile(session.session_id, { filename: filename.trim(), title: title.trim() })
                )
              }
              className="apple-btn-primary w-full"
            >
              <Check className="w-3.5 h-3.5" />
              Save
            </button>
          </div>
        </div>
      </div>

      {/* 2. Chain Renaming, Swapping, Reordering & Deletion */}
      <div className="apple-card p-4 space-y-3">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div className="flex items-center gap-2 text-[13.5px] font-semibold">
            <Layers className="w-4 h-4 text-[#5856D6]" />
            <span>Chain IDs & Ordering ({chains.length} chain{chains.length === 1 ? '' : 's'})</span>
          </div>
          <div className="flex flex-wrap items-center gap-1.5">
            {selectedForDelete.length > 0 && (
              <button
                disabled={busy || selectedForDelete.length >= chains.length}
                onClick={() =>
                  runOp(`Deleting chain(s) ${selectedForDelete.join(', ')}…`, () =>
                    api.deleteChains(session.session_id, selectedForDelete)
                  )
                }
                className="apple-btn-danger text-[12px]"
                title={
                  selectedForDelete.length >= chains.length
                    ? 'Keep at least one chain'
                    : 'Delete all checked chains'
                }
              >
                <Trash2 className="w-3.5 h-3.5" />
                Delete Selected ({selectedForDelete.length})
              </button>
            )}
            <button
              disabled={busy}
              onClick={() =>
                runOp('Standardizing chain IDs (A, B, C…)', () => api.standardizeChains(session.session_id))
              }
              className="apple-btn-secondary text-[12px]"
              title="Rename chains sequentially to A, B, C..."
            >
              Auto A, B, C…
            </button>
            <button
              disabled={busy || !hasChainChanges}
              onClick={() =>
                runOp('Applying chain renames…', () => api.renameChains(session.session_id, chainMap))
              }
              className="apple-btn-primary text-[12px]"
            >
              Apply Chain Renames
            </button>
          </div>
        </div>

        <p className="text-[11.5px] text-[#86868B]">
          Click any chain to spotlight &amp; zoom it in 3D (blue), or check chains to highlight them in red before deleting.
        </p>

        <div className="space-y-1.5">
          {chains.map((ch, idx) => {
            const isChecked = selectedForDelete.includes(ch.id);
            const isSpotlit =
              highlightedChains?.tone === 'select' && highlightedChains.ids.includes(ch.id);
            const otherChains = chains.filter((c) => c.id !== ch.id).map((c) => c.id);
            return (
              <div
                key={ch.id}
                onClick={() => toggleChainSpotlight(ch.id)}
                className={`cursor-pointer flex flex-wrap items-center justify-between gap-2 p-2.5 rounded-xl border transition ${
                  isChecked
                    ? 'bg-[#FF3B30]/[0.08] border-[#FF3B30]/40 ring-1 ring-[#FF3B30]/25'
                    : isSpotlit
                    ? 'bg-[#007AFF]/[0.08] border-[#007AFF]/45 ring-1 ring-[#007AFF]/25'
                    : 'bg-black/[0.025] dark:bg-white/[0.035] border-black/[0.04] dark:border-white/[0.05] hover:border-[#007AFF]/30'
                }`}
              >
                <div className="flex items-center gap-2.5">
                  {chains.length > 1 && (
                    <input
                      type="checkbox"
                      checked={isChecked}
                      onClick={(e) => e.stopPropagation()}
                      onChange={() => toggleDeleteCheck(ch.id)}
                      className="accent-[#FF3B30] w-4 h-4 cursor-pointer"
                      title={`Mark chain ${ch.id} for deletion (highlights red in 3D)`}
                    />
                  )}
                  <span
                    className={`w-8 h-8 rounded-xl font-mono font-bold text-[13px] flex items-center justify-center transition ${
                      isChecked
                        ? 'bg-[#FF3B30] text-white shadow-sm'
                        : isSpotlit
                        ? 'bg-[#007AFF] text-white shadow-sm'
                        : 'bg-[#007AFF]/12 text-[#007AFF] dark:text-[#0A84FF]'
                    }`}
                  >
                    {ch.id}
                  </span>
                  <span className="text-[#86868B] text-[12px]">→</span>
                  <input
                    value={chainMap[ch.id] ?? ch.id}
                    maxLength={4}
                    onClick={(e) => e.stopPropagation()}
                    onFocus={() => onHighlightChains({ ids: [ch.id], tone: 'select' }, true)}
                    onChange={(e) =>
                      setChainMap((prev) => ({ ...prev, [ch.id]: e.target.value.toUpperCase() }))
                    }
                    className="w-16 apple-input !py-1 text-center font-mono font-bold uppercase"
                    title="New Chain ID (1 letter recommended for PDB)"
                  />
                  <div className="text-[11.5px] leading-tight">
                    <div className="font-medium flex items-center gap-1.5">
                      <span>
                        {ch.protein_residues} protein res · {ch.waters} H₂O · {ch.atom_count} atoms
                      </span>
                      {isSpotlit && (
                        <span className="apple-pill bg-[#007AFF]/15 text-[#007AFF] dark:text-[#0A84FF]">
                          3D Spotlight
                        </span>
                      )}
                    </div>
                    <div className="text-[10.5px] text-[#86868B] font-mono">
                      {ch.seq_range ? `Range ${ch.seq_range[0]}–${ch.seq_range[1]}` : 'Heterogen-only'}
                      {ch.heterogens.length > 0 ? ` · HET: ${ch.heterogens.join(', ')}` : ''}
                    </div>
                  </div>
                </div>

                <div className="flex items-center gap-1" onClick={(e) => e.stopPropagation()}>
                  <button
                    onClick={() => toggleChainSpotlight(ch.id)}
                    title={isSpotlit ? 'Clear 3D spotlight' : `Spotlight Chain ${ch.id} in 3D`}
                    className={`p-1.5 rounded-lg transition ${
                      isSpotlit
                        ? 'bg-[#007AFF] text-white'
                        : 'hover:bg-black/5 dark:hover:bg-white/10 text-[#007AFF]'
                    }`}
                  >
                    <Eye className="w-3.5 h-3.5" />
                  </button>
                  <button
                    disabled={busy || idx === 0}
                    onClick={() => moveChainOrder(idx, -1)}
                    title="Move chain up"
                    className="p-1.5 rounded-lg hover:bg-black/5 dark:hover:bg-white/10 disabled:opacity-30"
                  >
                    <ArrowUp className="w-3.5 h-3.5" />
                  </button>
                  <button
                    disabled={busy || idx === chains.length - 1}
                    onClick={() => moveChainOrder(idx, 1)}
                    title="Move chain down"
                    className="p-1.5 rounded-lg hover:bg-black/5 dark:hover:bg-white/10 disabled:opacity-30"
                  >
                    <ArrowDown className="w-3.5 h-3.5" />
                  </button>
                  {chains.length > 1 && (
                    <button
                      disabled={busy}
                      onClick={() =>
                        runOp(`Keeping only chain ${ch.id}…`, () =>
                          api.deleteChains(session.session_id, otherChains)
                        )
                      }
                      title={`Delete all other chains and keep only ${ch.id}`}
                      className="apple-btn-secondary !py-1 !px-2 text-[11px]"
                    >
                      Keep Only
                    </button>
                  )}
                  <button
                    disabled={busy || chains.length <= 1}
                    onClick={() =>
                      runOp(`Deleting chain ${ch.id}…`, () => api.deleteChain(session.session_id, ch.id))
                    }
                    title={chains.length <= 1 ? 'At least 1 chain must remain' : `Delete chain ${ch.id}`}
                    className="apple-btn-danger !py-1 !px-2 text-[11px]"
                  >
                    <Trash2 className="w-3.5 h-3.5" />
                    <span>Delete</span>
                  </button>
                </div>
              </div>
            );
          })}
        </div>

        {/* Merge Chains Sub-bar */}
        {chains.length >= 2 && (
          <div className="pt-2.5 border-t border-black/[0.06] dark:border-white/[0.08] flex flex-wrap items-center justify-between gap-2 text-[12px]">
            <div className="flex flex-wrap items-center gap-2">
              <GitMerge className="w-4 h-4 text-[#AF52DE]" />
              <span className="font-semibold">Merge Chain</span>
              <select value={mergeSrc} onChange={(e) => setMergeSrc(e.target.value)} className="apple-select font-mono">
                {chains.map((c) => (
                  <option key={c.id} value={c.id}>
                    {c.id}
                  </option>
                ))}
              </select>
              <span className="text-[#86868B]">into</span>
              <select value={mergeTgt} onChange={(e) => setMergeTgt(e.target.value)} className="apple-select font-mono">
                {chains.map((c) => (
                  <option key={c.id} value={c.id}>
                    {c.id}
                  </option>
                ))}
              </select>
              <label className="flex items-center gap-1.5 text-[11.5px] text-[#6E6E73] dark:text-[#A1A1A6] cursor-pointer ml-1">
                <input
                  type="checkbox"
                  checked={mergeRenumber}
                  onChange={(e) => setMergeRenumber(e.target.checked)}
                  className="accent-[#007AFF]"
                />
                Renumber continuously
              </label>
            </div>
            <button
              disabled={busy || !mergeSrc || !mergeTgt || mergeSrc === mergeTgt}
              onClick={() =>
                runOp(`Merging chain ${mergeSrc} → ${mergeTgt}…`, () =>
                  api.mergeChains(session.session_id, {
                    source_chain: mergeSrc,
                    target_chain: mergeTgt,
                    renumber_merged: mergeRenumber,
                  })
                )
              }
              className="apple-btn-secondary text-[12px]"
            >
              Merge Chains
            </button>
          </div>
        )}
      </div>

      {/* 3. Per-Chain Residue Renumbering & Offset */}
      <div className="apple-card p-4 space-y-3">
        <div className="flex items-center gap-2 text-[13.5px] font-semibold">
          <Hash className="w-4 h-4 text-[#34C759]" />
          <span>Per-Chain Residue Numbering & Offset</span>
        </div>

        <div className="grid grid-cols-1 sm:grid-cols-12 gap-2.5 items-end text-[12px]">
          <div className="sm:col-span-3">
            <label className="text-[11px] text-[#86868B] block mb-1">Chain</label>
            <select value={renumChain} onChange={(e) => setRenumChain(e.target.value)} className="apple-select w-full font-mono">
              {chains.map((c) => (
                <option key={c.id} value={c.id}>
                  Chain {c.id}
                </option>
              ))}
            </select>
          </div>
          <div className="sm:col-span-4">
            <label className="text-[11px] text-[#86868B] block mb-1">Mode</label>
            <select
              value={renumMode}
              onChange={(e) => setRenumMode(e.target.value as 'start' | 'offset')}
              className="apple-select w-full"
            >
              <option value="start">Renumber Sequentially From</option>
              <option value="offset">Shift Numbers by Offset (+/-)</option>
            </select>
          </div>
          <div className="sm:col-span-2">
            <label className="text-[11px] text-[#86868B] block mb-1">
              {renumMode === 'start' ? 'Start #' : 'Offset'}
            </label>
            {renumMode === 'start' ? (
              <input
                type="number"
                min={1}
                value={startNum}
                onChange={(e) => setStartNum(Number(e.target.value))}
                className="apple-input !py-1.5 font-mono"
              />
            ) : (
              <input
                type="number"
                value={offsetVal}
                onChange={(e) => setOffsetVal(Number(e.target.value))}
                className="apple-input !py-1.5 font-mono"
              />
            )}
          </div>
          <div className="sm:col-span-3">
            <button
              disabled={busy || !renumChain}
              onClick={() =>
                runOp(`Renumbering chain ${renumChain}…`, () =>
                  api.renumberChain(session.session_id, {
                    chain_id: renumChain,
                    mode: renumMode,
                    start_num: startNum,
                    offset: offsetVal,
                    include_heterogens: includeHet,
                    clear_icodes: true,
                  })
                )
              }
              className="apple-btn-primary w-full"
            >
              Apply Numbering
            </button>
          </div>
        </div>

        <label className="inline-flex items-center gap-2 text-[11.5px] text-[#6E6E73] dark:text-[#A1A1A6] cursor-pointer">
          <input
            type="checkbox"
            checked={includeHet}
            onChange={(e) => setIncludeHet(e.target.checked)}
            className="accent-[#007AFF]"
          />
          Also renumber heterogens & waters in this chain (clears insertion codes)
        </label>
      </div>

      {/* 4. Ligand / Heterogen & Residue Code Renaming */}
      <div className="apple-card p-4 space-y-3">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div className="flex items-center gap-2 text-[13.5px] font-semibold">
            <Tag className="w-4 h-4 text-[#FF9500]" />
            <span>Rename Ligand / Residue 3-Letter Codes & Reassign Chain</span>
          </div>
          <button
            disabled={busy}
            onClick={() =>
              runOp('Canonicalizing residue names…', () => api.canonicalizeResidueNames(session.session_id))
            }
            className="apple-btn-secondary text-[11.5px]"
            title="Convert HID/HIE/HIP/CYX/CYM/ASH/GLH/LYN back to standard PDB HIS/CYS/ASP/GLU/LYS"
          >
            <RefreshCw className="w-3.5 h-3.5" />
            Reset AMBER/GROMACS Names → Standard PDB
          </button>
        </div>

        {detectedHets.length > 0 && (
          <div className="flex flex-wrap items-center gap-1.5">
            <span className="text-[11px] text-[#86868B] mr-1">Quick pick heterogen:</span>
            {detectedHets.map((h, i) => (
              <button
                key={`${h.chain}-${h.resseq}-${i}`}
                onClick={() => {
                  setOldResCode(h.resname);
                  setScopeChain(h.chain);
                  setNewResCode(h.resname === 'BEN' || h.resname === 'UNL' ? 'LIG' : h.resname);
                }}
                className="px-2 py-0.5 rounded-lg text-[11px] font-mono bg-[#FF9500]/12 text-[#B26200] dark:text-[#FF9F0A] hover:bg-[#FF9500]/22 transition"
              >
                {h.resname} ({h.chain}:{h.resseq})
              </button>
            ))}
          </div>
        )}

        <div className="grid grid-cols-2 sm:grid-cols-12 gap-2.5 items-end text-[12px]">
          <div className="sm:col-span-3">
            <label className="text-[11px] text-[#86868B] block mb-1">Current Code</label>
            <input
              value={oldResCode}
              maxLength={4}
              onChange={(e) => setOldResCode(e.target.value.toUpperCase())}
              placeholder="e.g. BEN"
              className="apple-input font-mono uppercase"
            />
          </div>
          <div className="sm:col-span-3">
            <label className="text-[11px] text-[#86868B] block mb-1">New Code</label>
            <input
              value={newResCode}
              maxLength={4}
              onChange={(e) => setNewResCode(e.target.value.toUpperCase())}
              placeholder="e.g. LIG"
              className="apple-input font-mono uppercase"
            />
          </div>
          <div className="sm:col-span-2">
            <label className="text-[11px] text-[#86868B] block mb-1">In Chain (opt)</label>
            <select value={scopeChain} onChange={(e) => setScopeChain(e.target.value)} className="apple-select w-full font-mono">
              <option value="">All</option>
              {chains.map((c) => (
                <option key={c.id} value={c.id}>
                  {c.id}
                </option>
              ))}
            </select>
          </div>
          <div className="sm:col-span-2">
            <label className="text-[11px] text-[#86868B] block mb-1">Move to Chain</label>
            <input
              value={moveToChain}
              maxLength={4}
              onChange={(e) => setMoveToChain(e.target.value.toUpperCase())}
              placeholder="e.g. L"
              className="apple-input font-mono uppercase"
              title="Leave blank to keep in current chain"
            />
          </div>
          <div className="sm:col-span-2">
            <button
              disabled={busy || !oldResCode.trim() || !newResCode.trim()}
              onClick={() =>
                runOp(`Renaming ${oldResCode} → ${newResCode}…`, () =>
                  api.renameResidueCode(session.session_id, {
                    old_resname: oldResCode.trim(),
                    new_resname: newResCode.trim(),
                    chain_id: scopeChain || null,
                    move_to_chain: moveToChain.trim() || null,
                  })
                )
              }
              className="apple-btn-primary w-full"
            >
              Rename
            </button>
          </div>
        </div>

        {/* Field Cleanup Utilities */}
        <div className="pt-3 border-t border-black/[0.06] dark:border-white/[0.08] flex flex-wrap items-center gap-2">
          <button
            disabled={busy}
            onClick={() =>
              runOp('Resetting occupancies to 1.00…', () =>
                api.normalizeFields(session.session_id, { reset_occupancy: true })
              )
            }
            className="apple-btn-secondary text-[12px]"
          >
            Set All Occ = 1.00
          </button>

          <div className="flex items-center gap-1.5">
            <input
              type="number"
              value={uniformB}
              onChange={(e) => setUniformB(e.target.value)}
              className="w-16 apple-input !py-1 text-center font-mono text-[12px]"
              title="Uniform B-factor value"
            />
            <button
              disabled={busy}
              onClick={() =>
                runOp(`Setting B-factors to ${uniformB}…`, () =>
                  api.normalizeFields(session.session_id, {
                    reset_occupancy: false,
                    reset_bfactor: parseFloat(uniformB) || 20.0,
                  })
                )
              }
              className="apple-btn-secondary text-[12px]"
            >
              Set Uniform B
            </button>
          </div>

          <div className="flex items-center gap-1.5">
            <input
              value={waterChain}
              maxLength={2}
              onChange={(e) => setWaterChain(e.target.value.toUpperCase())}
              className="w-11 apple-input !py-1 text-center font-mono text-[12px] uppercase"
              title="Target water chain ID"
            />
            <button
              disabled={busy || !waterChain.trim()}
              onClick={() =>
                runOp(`Consolidating waters to chain ${waterChain}…`, () =>
                  api.normalizeFields(session.session_id, {
                    reset_occupancy: false,
                    strip_waters_to_chain: waterChain.trim(),
                  })
                )
              }
              className="apple-btn-secondary text-[12px]"
            >
              Move All Waters → Chain {waterChain || 'W'}
            </button>
          </div>
        </div>
      </div>

      {/* 5. Live PDB Coordinate Preview */}
      <div className="apple-card p-4 space-y-2.5">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div className="flex items-center gap-2 text-[13px] font-semibold">
            <Code2 className="w-4 h-4 text-[#007AFF]" />
            <span>Live PDB Coordinate Preview</span>
          </div>
          <div className="apple-segmented">
            {(
              [
                ['head', 'Header + First 45'],
                ['all_chains', 'Chain Boundaries'],
                ['hetatm', 'HETATM / Ligands'],
              ] as const
            ).map(([k, lbl]) => (
              <button
                key={k}
                onClick={() => setPdbFilter(k)}
                className={`apple-segmented-item !text-[11px] !px-2 ${
                  pdbFilter === k ? 'apple-segmented-active' : 'text-[#6E6E73]'
                }`}
              >
                {lbl}
              </button>
            ))}
          </div>
        </div>
        <pre className="p-3 rounded-xl bg-black/[0.035] dark:bg-black/40 text-[11px] font-mono overflow-x-auto max-h-56 leading-relaxed">
          {pdbPreviewText}
        </pre>
      </div>
    </div>
  );
};
