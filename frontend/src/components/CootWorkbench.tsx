import React, { useEffect, useMemo, useState } from 'react';
import {
  CheckCircle2,
  Dna,
  FlipHorizontal2,
  GitCompare,
  RefreshCw,
  Search,
  ShieldAlert,
  Sparkles,
  Trash2,
  Wand2,
} from 'lucide-react';
import { api } from '../api/client';
import type {
  ResidueIndexItem,
  RotamerExplorerData,
  SelectedResidue,
  StructureSessionState,
} from '../types/pdbprep';

interface CootWorkbenchProps {
  session: StructureSessionState;
  selectedResidue: SelectedResidue | null;
  busy: boolean;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
  onStateUpdate: (state: StructureSessionState, toastMsg: string) => void;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
}

const STANDARD_20 = [
  'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
  'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
  'LEU', 'LYS', 'MET', 'PHE', 'PRO',
  'SER', 'THR', 'TRP', 'TYR', 'VAL',
];

export const CootWorkbench: React.FC<CootWorkbenchProps> = ({
  session,
  selectedResidue,
  busy,
  onSelectResidue,
  onStateUpdate,
  onSetBusy,
  onError,
}) => {
  const [filterMode, setFilterMode] = useState<'all' | 'issues' | 'rotamers' | 'clashes'>('issues');
  const [search, setSearch] = useState('');
  const [rotamerData, setRotamerData] = useState<RotamerExplorerData | null>(null);
  const [loadingRotamers, setLoadingRotamers] = useState(false);
  const [showMutateGrid, setShowMutateGrid] = useState(false);
  const [zoneRadius, setZoneRadius] = useState(2);
  const [lastRefineSummary, setLastRefineSummary] = useState<{
    before_issues: number;
    after_issues: number;
    before_score: number;
    after_score: number;
    before_clashscore: number;
    after_clashscore: number;
    actions: string[];
  } | null>(null);

  const activeRes: ResidueIndexItem | null = useMemo(() => {
    if (!selectedResidue) return session.validation.residues[0] || null;
    return (
      session.validation.residues.find(
        (r) => r.chain === selectedResidue.chain && r.resseq === selectedResidue.resseq
      ) ||
      session.validation.residues[0] ||
      null
    );
  }, [session.validation.residues, selectedResidue]);

  const filteredResidues = useMemo(() => {
    return session.validation.residues.filter((r) => {
      if (filterMode === 'issues' && r.issues.length === 0) return false;
      if (filterMode === 'rotamers' && r.rotamer_status !== 'outlier') return false;
      if (filterMode === 'clashes' && r.clash_count === 0) return false;
      if (search.trim()) {
        const q = search.trim().toUpperCase();
        return (
          r.key.toUpperCase().includes(q) ||
          r.resname.toUpperCase().includes(q) ||
          String(r.resseq).includes(q)
        );
      }
      return true;
    });
  }, [session.validation.residues, filterMode, search]);

  // Fetch ranked Dunbrack/Richardson rotamers for active residue
  useEffect(() => {
    if (!activeRes) {
      setRotamerData(null);
      return;
    }
    let cancelled = false;
    setLoadingRotamers(true);
    api
      .cootListRotamers(session.session_id, activeRes.chain, activeRes.resseq, activeRes.icode)
      .then((d) => {
        if (!cancelled) setRotamerData(d);
      })
      .catch(() => {
        if (!cancelled) setRotamerData(null);
      })
      .finally(() => {
        if (!cancelled) setLoadingRotamers(false);
      });
    return () => {
      cancelled = true;
    };
  }, [session.session_id, session.pdb_text, activeRes?.chain, activeRes?.resseq]);

  const runAction = async (label: string, fn: () => Promise<{ result: { message?: string }; state: StructureSessionState }>) => {
    onSetBusy(label);
    try {
      const resp = await fn();
      onStateUpdate(resp.state, resp.result?.message || label);
    } catch (e: any) {
      onError(e?.message || 'COOT operation failed');
    } finally {
      onSetBusy(null);
    }
  };

  return (
    <div className="space-y-3.5 pb-6">
      {/* Header + Auto-Refine All */}
      <div className="apple-card p-4 flex flex-wrap items-center justify-between gap-3">
        <div>
          <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#007AFF] dark:text-[#0A84FF]">
            COOT Interactive Model Builder
          </div>
          <h2 className="text-[16px] font-semibold tracking-[-0.015em]">
            Real-Space Residue Inspector & Editing Suite
          </h2>
        </div>
        <button
          disabled={busy}
          onClick={() =>
            runAction('Running COOT Auto-Refine & Re-Validating…', async () => {
              const r = await api.cootAutoFix(session.session_id, {
                fix_clashes: true,
                fix_rotamers: true,
                fix_nqh_flips: true,
                fix_cis_nonpro: true,
                fix_ramachandran: true,
                regularize_outliers: true,
              });
              setLastRefineSummary({
                before_issues: r.result.before_issues,
                after_issues: r.result.after_issues,
                before_score: r.result.before_score,
                after_score: r.result.after_score,
                before_clashscore: r.result.before_clashscore,
                after_clashscore: r.result.after_clashscore,
                actions: r.result.actions,
              });
              return {
                result: {
                  message: `COOT Refined & Re-Validated: Issues ${r.result.before_issues} → ${r.result.after_issues} · Quality ${r.result.before_score} → ${r.result.after_score}`,
                },
                state: r.state,
              };
            })
          }
          className="apple-btn-primary"
        >
          <Sparkles className="w-3.5 h-3.5" />
          <span>COOT Auto-Refine All</span>
        </button>
      </div>

      {/* Live Re-Validation Summary Banner after COOT Auto-Refine */}
      {lastRefineSummary && (
        <div className="apple-card p-3.5 bg-[#34C759]/[0.08] border-[#34C759]/30 space-y-1.5">
          <div className="flex flex-wrap items-center justify-between gap-2 text-[12.5px]">
            <span className="font-semibold text-[#1B8A38] dark:text-[#30D158]">
              Refined &amp; Re-Validated: Issues {lastRefineSummary.before_issues} → {lastRefineSummary.after_issues} · Quality{' '}
              {lastRefineSummary.before_score} → {lastRefineSummary.after_score} · Clashscore{' '}
              {lastRefineSummary.before_clashscore} → {lastRefineSummary.after_clashscore}
            </span>
            <button
              onClick={() => setLastRefineSummary(null)}
              className="text-[11px] text-[#6E6E73] hover:text-[#1D1D1F] dark:hover:text-white"
            >
              Dismiss
            </button>
          </div>
          {lastRefineSummary.actions.length > 0 && (
            <div className="text-[11px] text-[#6E6E73] dark:text-[#A1A1A6]">
              {lastRefineSummary.actions.join(' · ')}
            </div>
          )}
        </div>
      )}

      {/* Residue Picker Bar */}
      <div className="apple-card p-3.5 space-y-2.5">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div className="apple-segmented">
            {(
              [
                ['issues', `Issues (${session.validation.residues.filter((r) => r.issues.length > 0).length})`],
                ['all', `All (${session.validation.residues.length})`],
                ['clashes', 'Clashes'],
                ['rotamers', 'Rotamers'],
              ] as const
            ).map(([k, lbl]) => (
              <button
                key={k}
                onClick={() => setFilterMode(k)}
                className={`apple-segmented-item ${filterMode === k ? 'apple-segmented-active' : 'text-[#6E6E73]'}`}
              >
                {lbl}
              </button>
            ))}
          </div>
          <div className="relative w-44">
            <Search className="w-3.5 h-3.5 text-[#86868B] absolute left-2.5 top-2.5" />
            <input
              value={search}
              onChange={(e) => setSearch(e.target.value)}
              placeholder="Find e.g. PHE 14"
              className="apple-input !pl-8 !py-1.5 text-[12px]"
            />
          </div>
        </div>

        {/* Scrollable horizontal residue pills */}
        <div className="flex items-center gap-1.5 overflow-x-auto pb-1">
          {filteredResidues.length === 0 ? (
            <div className="flex items-center gap-2 px-3 py-2 rounded-xl bg-[#34C759]/15 border border-[#34C759]/30 text-[12.5px] text-[#1B8A38] dark:text-[#30D158] font-semibold">
              <CheckCircle2 className="w-4 h-4 shrink-0" />
              <span>
                All Clear — 0 {filterMode === 'issues' ? 'issues' : filterMode} across all{' '}
                {session.validation.residues.length} residues (Quality {session.validation.quality_score}%)
              </span>
            </div>
          ) : (
            filteredResidues.slice(0, 120).map((r) => {
              const isSel = activeRes?.key === r.key;
              const hasSevere =
                r.issues.includes('clash') ||
                r.issues.includes('rama_outlier') ||
                r.issues.includes('cis_nonpro') ||
                r.issues.includes('chiral_outlier');
              return (
                <button
                  key={r.key}
                  onClick={() =>
                    onSelectResidue(
                      { chain: r.chain, resseq: r.resseq, icode: r.icode, resname: r.resname, key: r.key },
                      true
                    )
                  }
                  className={`shrink-0 px-2.5 py-1.5 rounded-xl text-[11.5px] font-mono border transition ${
                    isSel
                      ? 'bg-[#007AFF] text-white border-[#007AFF] shadow-sm'
                      : hasSevere
                      ? 'bg-[#FF3B30]/10 border-[#FF3B30]/30 text-[#C92519] dark:text-[#FF453A] hover:bg-[#FF3B30]/18'
                      : r.issues.length > 0
                      ? 'bg-[#FF9500]/10 border-[#FF9500]/30 text-[#B26200] dark:text-[#FF9F0A] hover:bg-[#FF9500]/18'
                      : 'bg-[#34C759]/10 border-[#34C759]/25 text-[#1B8A38] dark:text-[#30D158] hover:bg-[#34C759]/20'
                  }`}
                >
                  <span className="font-semibold">{r.resname}</span> {r.key}
                </button>
              );
            })
          )}
        </div>
      </div>

      {/* Selected Residue Deep COOT Inspector */}
      {activeRes && (
        <div className="apple-card p-4 space-y-4">
          <div className="flex flex-wrap items-start justify-between gap-2 pb-3 border-b border-black/[0.06] dark:border-white/[0.08]">
            <div>
              <div className="flex items-center gap-2">
                <span className="text-[18px] font-bold font-mono tracking-tight">
                  {activeRes.resname} {activeRes.chain}:{activeRes.resseq}
                  {activeRes.icode}
                </span>
                {activeRes.issues.length === 0 ? (
                  <span className="apple-pill bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158]">
                    Geometry Ideal
                  </span>
                ) : (
                  activeRes.issues.map((iss) => (
                    <span
                      key={iss}
                      className="apple-pill bg-[#FF3B30]/12 text-[#D70015] dark:text-[#FF453A]"
                    >
                      {iss.replace('_', ' ')}
                    </span>
                  ))
                )}
              </div>
              <div className="flex flex-wrap items-center gap-3 mt-1.5 text-[11.5px] font-mono text-[#6E6E73] dark:text-[#A1A1A6]">
                <span>φ: {activeRes.phi !== null ? `${activeRes.phi}°` : '—'}</span>
                <span>ψ: {activeRes.psi !== null ? `${activeRes.psi}°` : '—'}</span>
                <span>ω: {activeRes.omega !== null ? `${activeRes.omega}°` : '—'}</span>
                <span>B: {activeRes.mean_b} Å²</span>
                <span>Occ: {activeRes.min_occ}</span>
                <span>Clashes: {activeRes.clash_count}</span>
              </div>
            </div>

            <button
              onClick={() =>
                onSelectResidue(
                  {
                    chain: activeRes.chain,
                    resseq: activeRes.resseq,
                    icode: activeRes.icode,
                    resname: activeRes.resname,
                    key: activeRes.key,
                  },
                  true
                )
              }
              className="apple-btn-secondary text-[12px]"
            >
              Center 3D
            </button>
          </div>

          {/* COOT Tool Action Palette */}
          <div>
            <div className="text-[11px] font-semibold uppercase tracking-[0.05em] text-[#86868B] mb-2">
              COOT Model-Building Actions
            </div>
            <div className="grid grid-cols-2 sm:grid-cols-3 gap-2">
              <button
                disabled={busy}
                onClick={() =>
                  runAction(`Auto-fitting rotamer for ${activeRes.key}…`, () =>
                    api.cootApplyRotamer(session.session_id, {
                      chain: activeRes.chain,
                      resseq: activeRes.resseq,
                      icode: activeRes.icode,
                      rotamer_index: null,
                    })
                  )
                }
                className="apple-btn-secondary justify-start"
              >
                <Wand2 className="w-3.5 h-3.5 text-[#007AFF]" />
                <span>Auto-Fit Rotamer</span>
              </button>

              <button
                disabled={busy}
                onClick={() => setShowMutateGrid((v) => !v)}
                className={`apple-btn-secondary justify-start ${showMutateGrid ? '!bg-[#007AFF]/15 !text-[#007AFF]' : ''}`}
              >
                <Dna className="w-3.5 h-3.5 text-[#AF52DE]" />
                <span>Mutate Residue</span>
              </button>

              <button
                disabled={busy}
                onClick={() =>
                  runAction(`COOT Pepflip ${activeRes.key}…`, () =>
                    api.cootPepflip(session.session_id, {
                      chain: activeRes.chain,
                      resseq: activeRes.resseq,
                      icode: activeRes.icode,
                    })
                  )
                }
                className="apple-btn-secondary justify-start"
              >
                <FlipHorizontal2 className="w-3.5 h-3.5 text-[#5856D6]" />
                <span>180° Pepflip</span>
              </button>

              <button
                disabled={busy}
                onClick={() =>
                  runAction(`Cis ↔ Trans flip on ${activeRes.key}…`, () =>
                    api.cootCisTrans(session.session_id, {
                      chain: activeRes.chain,
                      resseq: activeRes.resseq,
                      icode: activeRes.icode,
                    })
                  )
                }
                className="apple-btn-secondary justify-start"
              >
                <GitCompare className="w-3.5 h-3.5 text-[#FF9500]" />
                <span>Cis ↔ Trans</span>
              </button>

              <button
                disabled={busy || !['ASN', 'GLN', 'HIS'].includes(activeRes.std_resname)}
                onClick={() =>
                  runAction(`180° NQH flip on ${activeRes.key}…`, () =>
                    api.cootNqhFlip(session.session_id, {
                      chain: activeRes.chain,
                      resseq: activeRes.resseq,
                      icode: activeRes.icode,
                    })
                  )
                }
                className="apple-btn-secondary justify-start"
                title="Flip ASN/GLN amide or HIS imidazole 180°"
              >
                <RefreshCw className="w-3.5 h-3.5 text-[#34C759]" />
                <span>NQH 180° Flip</span>
              </button>

              <button
                disabled={busy}
                onClick={() =>
                  runAction(`Fixing Cα chirality on ${activeRes.key}…`, () =>
                    api.cootFixChirality(session.session_id, {
                      chain: activeRes.chain,
                      resseq: activeRes.resseq,
                      icode: activeRes.icode,
                    })
                  )
                }
                className="apple-btn-secondary justify-start"
              >
                <ShieldAlert className="w-3.5 h-3.5 text-[#FF3B30]" />
                <span>Fix L-Chirality</span>
              </button>
            </div>
          </div>

          {/* Point Mutation Grid (collapsible) */}
          {showMutateGrid && (
            <div className="p-3 rounded-2xl bg-[#AF52DE]/8 border border-[#AF52DE]/20 space-y-2">
              <div className="text-[12px] font-semibold text-[#8944AB] dark:text-[#BF5AF2]">
                Select Target Amino Acid (Mutate & Auto-Fit Clash-Free Rotamer)
              </div>
              <div className="grid grid-cols-5 sm:grid-cols-10 gap-1.5">
                {STANDARD_20.map((aa) => (
                  <button
                    key={aa}
                    disabled={busy}
                    onClick={() => {
                      setShowMutateGrid(false);
                      runAction(`Mutating ${activeRes.key} → ${aa}…`, () =>
                        api.cootMutate(session.session_id, {
                          chain: activeRes.chain,
                          resseq: activeRes.resseq,
                          icode: activeRes.icode,
                          target_resname: aa,
                        })
                      );
                    }}
                    className={`py-1.5 rounded-lg text-[11.5px] font-mono font-semibold transition ${
                      aa === activeRes.std_resname
                        ? 'bg-[#AF52DE] text-white'
                        : 'bg-white/80 dark:bg-white/10 hover:bg-[#AF52DE]/20'
                    }`}
                  >
                    {aa}
                  </button>
                ))}
              </div>
            </div>
          )}

          {/* AltLoc Picker if dual conformers exist */}
          {activeRes.altlocs.length > 1 && (
            <div className="p-3 rounded-2xl bg-[#FF9500]/10 border border-[#FF9500]/25 flex flex-wrap items-center justify-between gap-2">
              <div className="text-[12.5px] font-medium text-[#B26200] dark:text-[#FF9F0A]">
                Alternate Conformations Detected ({activeRes.altlocs.join(', ')}) — Pick Single Conformer:
              </div>
              <div className="flex gap-1.5">
                {activeRes.altlocs.map((al) => (
                  <button
                    key={al}
                    disabled={busy}
                    onClick={() =>
                      runAction(`Resolving AltLoc ${activeRes.key} → ${al}…`, () =>
                        api.cootAltloc(session.session_id, {
                          chain: activeRes.chain,
                          resseq: activeRes.resseq,
                          icode: activeRes.icode,
                          altloc: al,
                        })
                      )
                    }
                    className="apple-btn-primary !py-1 !px-3 !bg-[#FF9500] hover:!bg-[#E08600]"
                  >
                    Keep {al}
                  </button>
                ))}
              </div>
            </div>
          )}

          {/* COOT Regularize Zone Control */}
          <div className="p-3 rounded-2xl bg-black/[0.03] dark:bg-white/[0.04] flex flex-wrap items-center justify-between gap-3">
            <div>
              <div className="text-[12.5px] font-semibold">COOT Regularize Zone (Local Idealization)</div>
              <div className="text-[11px] text-[#86868B]">
                Zone {activeRes.chain}:{Math.max(1, activeRes.resseq - zoneRadius)}–{activeRes.resseq + zoneRadius} (harmonic Engh & Huber + peptide relaxation)
              </div>
            </div>
            <div className="flex items-center gap-2">
              <select
                value={zoneRadius}
                onChange={(e) => setZoneRadius(Number(e.target.value))}
                className="apple-select text-[12px]"
              >
                <option value={1}>±1 res (3-mer)</option>
                <option value={2}>±2 res (5-mer)</option>
                <option value={4}>±4 res (9-mer)</option>
              </select>
              <button
                disabled={busy}
                onClick={() =>
                  runAction('Regularizing zone…', () =>
                    api.cootRegularize(session.session_id, {
                      chain: activeRes.chain,
                      start_resseq: Math.max(1, activeRes.resseq - zoneRadius),
                      end_resseq: activeRes.resseq + zoneRadius,
                      iterations: 60,
                    })
                  )
                }
                className="apple-btn-primary !py-1.5"
              >
                Regularize Zone
              </button>
              <button
                disabled={busy}
                onClick={() =>
                  runAction(`Deleting ${activeRes.key}…`, () =>
                    api.cootDeleteResidue(session.session_id, {
                      chain: activeRes.chain,
                      resseq: activeRes.resseq,
                      icode: activeRes.icode,
                    })
                  )
                }
                className="apple-btn-danger !py-1.5"
                title="Delete residue"
              >
                <Trash2 className="w-3.5 h-3.5" />
              </button>
            </div>
          </div>

          {/* Dunbrack / Richardson Rotamer Library Explorer */}
          <div className="space-y-2">
            <div className="flex items-center justify-between">
              <div className="text-[12.5px] font-semibold">
                Dunbrack / Richardson Rotamer Library ({activeRes.std_resname})
              </div>
              {rotamerData?.current && (
                <div className="text-[11.5px] font-mono text-[#86868B]">
                  Current χ: [{rotamerData.current.chis.join(', ')}°] · {rotamerData.current.status}
                </div>
              )}
            </div>

            {loadingRotamers ? (
              <div className="text-[12px] text-[#86868B] py-2">Scoring rotamers against 3D environment…</div>
            ) : !rotamerData || rotamerData.rotamers.length === 0 ? (
              <div className="text-[12px] text-[#86868B] py-1">
                {activeRes.std_resname} has no rotatable side-chain χ angles.
              </div>
            ) : (
              <div className="grid grid-cols-1 sm:grid-cols-2 gap-1.5">
                {rotamerData.rotamers.map((rot) => (
                  <button
                    key={rot.index}
                    disabled={busy}
                    onClick={() =>
                      runAction(`Applying rotamer ${rot.name}…`, () =>
                        api.cootApplyRotamer(session.session_id, {
                          chain: activeRes.chain,
                          resseq: activeRes.resseq,
                          icode: activeRes.icode,
                          rotamer_index: rot.index,
                        })
                      )
                    }
                    className="flex items-center justify-between px-3 py-2 rounded-xl bg-black/[0.03] hover:bg-[#007AFF]/10 dark:bg-white/[0.04] dark:hover:bg-[#007AFF]/15 border border-transparent hover:border-[#007AFF]/30 transition text-left"
                  >
                    <div>
                      <div className="text-[12.5px] font-semibold font-mono">
                        {rot.name}{' '}
                        <span className="text-[11px] font-normal text-[#86868B]">({rot.probability}%)</span>
                      </div>
                      <div className="text-[10.5px] font-mono text-[#86868B]">
                        χ: [{rot.chis.join(', ')}°]
                      </div>
                    </div>
                    <span
                      className={`apple-pill font-mono ${
                        rot.clashes === 0
                          ? 'bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158]'
                          : 'bg-[#FF3B30]/15 text-[#D70015] dark:text-[#FF453A]'
                      }`}
                    >
                      {rot.clashes === 0 ? '0 clashes' : `${rot.clashes} clash`}
                    </span>
                  </button>
                ))}
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
};
