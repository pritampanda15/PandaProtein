import React from 'react';
import { AlertTriangle, Flame, RefreshCw, Sparkles, Wand2 } from 'lucide-react';
import { api } from '../api/client';
import type { SelectedResidue, StructureSessionState } from '../types/pdbprep';

interface ClashesPanelProps {
  session: StructureSessionState;
  busy: boolean;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
  onStateUpdate: (state: StructureSessionState, toast: string) => void;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
}

export const ClashesPanel: React.FC<ClashesPanelProps> = ({
  session,
  busy,
  onSelectResidue,
  onStateUpdate,
  onSetBusy,
  onError,
}) => {
  const { clashes, rotamers, nqh_flips } = session.validation;

  const runAction = async (label: string, fn: () => Promise<{ result: { message?: string }; state: StructureSessionState }>) => {
    onSetBusy(label);
    try {
      const resp = await fn();
      onStateUpdate(resp.state, resp.result?.message || label);
    } catch (e: any) {
      onError(e?.message || 'Action failed');
    } finally {
      onSetBusy(null);
    }
  };

  return (
    <div className="space-y-4 pb-6">
      <div className="apple-card p-4 flex flex-wrap items-center justify-between gap-3">
        <div>
          <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#FF3B30] dark:text-[#FF453A]">
            MolProbity All-Atom Contacts & COOT Side-Chain Repair
          </div>
          <h2 className="text-[16px] font-semibold tracking-[-0.015em]">
            Steric Clashes, Rotamer Outliers & NQH 180° Flips
          </h2>
        </div>
        <button
          disabled={busy}
          onClick={() =>
            runAction('Resolving steric clashes, rotamers & NQH flips…', async () => {
              const r = await api.cootAutoFix(session.session_id, {
                fix_clashes: true,
                fix_rotamers: true,
                fix_nqh_flips: true,
                fix_cis_nonpro: false,
                fix_ramachandran: false,
                regularize_outliers: false,
              });
              return {
                result: { message: `Clashes & rotamers resolved (Quality ${r.result.before_score} → ${r.result.after_score})` },
                state: r.state,
              };
            })
          }
          className="apple-btn-primary"
        >
          <Sparkles className="w-3.5 h-3.5" />
          <span>Auto-Resolve All Clashes</span>
        </button>
      </div>

      {/* KPI Strip */}
      <div className="grid grid-cols-3 gap-2.5">
        <div className="apple-card p-3.5">
          <div className="text-[11px] text-[#86868B]">Clashscore</div>
          <div
            className={`text-[20px] font-bold font-mono ${
              clashes.clashscore <= 5 ? 'text-[#34C759]' : clashes.clashscore <= 15 ? 'text-[#FF9500]' : 'text-[#FF3B30]'
            }`}
          >
            {clashes.clashscore}
          </div>
          <div className="text-[10.5px] text-[#86868B]">{clashes.total_clashes} overlaps &gt;0.4 Å</div>
        </div>
        <div className="apple-card p-3.5">
          <div className="text-[11px] text-[#86868B]">Rotamer Outliers</div>
          <div className={`text-[20px] font-bold font-mono ${rotamers.counts.outlier === 0 ? 'text-[#34C759]' : 'text-[#FF9500]'}`}>
            {rotamers.counts.outlier}
          </div>
          <div className="text-[10.5px] text-[#86868B]">{rotamers.favored_pct}% favored</div>
        </div>
        <div className="apple-card p-3.5">
          <div className="text-[11px] text-[#86868B]">NQH Flip Candidates</div>
          <div className="text-[20px] font-bold font-mono text-[#007AFF] dark:text-[#0A84FF]">
            {nqh_flips.length}
          </div>
          <div className="text-[10.5px] text-[#86868B]">Asn / Gln / His H-bond</div>
        </div>
      </div>

      {/* Rotamer Outliers */}
      <div className="apple-card p-4 space-y-2.5">
        <div className="flex items-center gap-2 text-[13.5px] font-semibold">
          <AlertTriangle className="w-4 h-4 text-[#FF9500]" />
          <span>Rotamer Outliers ({rotamers.outliers.length})</span>
        </div>
        {rotamers.outliers.length === 0 ? (
          <div className="text-[12px] text-[#34C759]">No side-chain rotamer outliers detected.</div>
        ) : (
          <div className="space-y-1.5 max-h-56 overflow-y-auto">
            {rotamers.outliers.map((rot) => (
              <div
                key={rot.key}
                className="flex items-center justify-between p-2.5 rounded-xl bg-[#FF9500]/10 border border-[#FF9500]/20 text-[12px]"
              >
                <button
                  onClick={() =>
                    onSelectResidue(
                      { chain: rot.chain, resseq: rot.resseq, icode: rot.icode, resname: rot.resname, key: rot.key },
                      true
                    )
                  }
                  className="text-left"
                >
                  <div className="font-mono font-semibold">
                    {rot.resname} {rot.key}
                  </div>
                  <div className="text-[11px] font-mono text-[#86868B]">
                    χ: [{rot.chis.join(', ')}°] · RMSD {rot.rmsd_deg}° from {rot.nearest_rotamer}
                  </div>
                </button>
                <button
                  disabled={busy}
                  onClick={() =>
                    runAction(`Auto-fitting rotamer on ${rot.key}…`, () =>
                      api.cootApplyRotamer(session.session_id, {
                        chain: rot.chain,
                        resseq: rot.resseq,
                        icode: rot.icode,
                        rotamer_index: null,
                      })
                    )
                  }
                  className="apple-btn-primary !py-1 !px-2.5 text-[11.5px]"
                >
                  <Wand2 className="w-3 h-3" />
                  Auto-Fit
                </button>
              </div>
            ))}
          </div>
        )}
      </div>

      {/* Asn / Gln / His 180° Flips */}
      <div className="apple-card p-4 space-y-2.5">
        <div className="flex items-center gap-2 text-[13.5px] font-semibold">
          <RefreshCw className="w-4 h-4 text-[#007AFF]" />
          <span>Asn / Gln / His 180° Flip Recommendations ({nqh_flips.length})</span>
        </div>
        {nqh_flips.length === 0 ? (
          <div className="text-[12px] text-[#34C759]">
            All Asn/Gln/His side-chain orientations are already H-bond optimal.
          </div>
        ) : (
          <div className="space-y-1.5">
            {nqh_flips.map((nqh) => (
              <div
                key={nqh.key}
                className="flex items-center justify-between p-2.5 rounded-xl bg-[#007AFF]/8 border border-[#007AFF]/20 text-[12px]"
              >
                <button
                  onClick={() =>
                    onSelectResidue(
                      { chain: nqh.chain, resseq: nqh.resseq, icode: nqh.icode, resname: nqh.resname, key: nqh.key },
                      true
                    )
                  }
                  className="text-left"
                >
                  <div className="font-mono font-semibold">
                    {nqh.resname} {nqh.key}{' '}
                    <span className="text-[#007AFF]">(+{nqh.hbond_gain} H-bond score)</span>
                  </div>
                  <div className="text-[11px] text-[#86868B]">{nqh.reason}</div>
                </button>
                <button
                  disabled={busy}
                  onClick={() =>
                    runAction(`Flipping ${nqh.key} 180°…`, () =>
                      api.cootNqhFlip(session.session_id, {
                        chain: nqh.chain,
                        resseq: nqh.resseq,
                        icode: nqh.icode,
                      })
                    )
                  }
                  className="apple-btn-primary !py-1 !px-2.5 text-[11.5px]"
                >
                  Flip 180°
                </button>
              </div>
            ))}
          </div>
        )}
      </div>

      {/* Worst Atomic Steric Clashes Table */}
      <div className="apple-card p-4 space-y-2.5">
        <div className="flex items-center gap-2 text-[13.5px] font-semibold">
          <Flame className="w-4 h-4 text-[#FF3B30]" />
          <span>Worst Atomic Steric Clashes (Click to Inspect in 3D)</span>
        </div>
        {clashes.clashes.length === 0 ? (
          <div className="text-[12px] text-[#34C759]">Zero steric clashes detected!</div>
        ) : (
          <div className="max-h-72 overflow-y-auto space-y-1">
            {clashes.clashes.slice(0, 40).map((cl, idx) => (
              <div
                key={idx}
                className="flex items-center justify-between px-3 py-2 rounded-xl bg-black/[0.025] dark:bg-white/[0.03] hover:bg-[#FF3B30]/10 transition text-[12px]"
              >
                <button
                  onClick={() =>
                    onSelectResidue(
                      {
                        chain: cl.atom1.chain,
                        resseq: cl.atom1.resseq,
                        icode: cl.atom1.icode,
                        resname: cl.atom1.resname,
                        key: cl.atom1.key,
                      },
                      true
                    )
                  }
                  className="text-left font-mono"
                >
                  <span className="font-semibold">
                    {cl.atom1.resname} {cl.atom1.key}.{cl.atom1.atom}
                  </span>{' '}
                  <span className="text-[#FF3B30]">⚡</span>{' '}
                  <span className="font-semibold">
                    {cl.atom2.resname} {cl.atom2.key}.{cl.atom2.atom}
                  </span>
                </button>
                <div className="flex items-center gap-2">
                  <span className="font-mono text-[11px] text-[#86868B]">{cl.distance} Å</span>
                  <span
                    className={`apple-pill font-mono ${
                      cl.severity === 'severe'
                        ? 'bg-[#FF3B30]/15 text-[#D70015] dark:text-[#FF453A]'
                        : 'bg-[#FF9500]/15 text-[#B26200] dark:text-[#FF9F0A]'
                    }`}
                  >
                    −{cl.overlap} Å
                  </span>
                  <button
                    disabled={busy}
                    onClick={() =>
                      runAction(`Auto-fitting ${cl.atom1.key}…`, () =>
                        api.cootApplyRotamer(session.session_id, {
                          chain: cl.atom1.chain,
                          resseq: cl.atom1.resseq,
                          icode: cl.atom1.icode,
                          rotamer_index: null,
                        })
                      )
                    }
                    className="apple-btn-secondary !py-0.5 !px-2 text-[11px]"
                  >
                    Fix
                  </button>
                </div>
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
};
