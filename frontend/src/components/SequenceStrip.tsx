import React, { useEffect, useState } from 'react';
import type { SelectedResidue, SequenceResidue } from '../types/pdbprep';

interface SequenceStripProps {
  sequenceChains: Record<string, SequenceResidue[]>;
  selectedResidue: SelectedResidue | null;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
}

export const SequenceStrip: React.FC<SequenceStripProps> = ({
  sequenceChains,
  selectedResidue,
  onSelectResidue,
}) => {
  const chainIds = Object.keys(sequenceChains || {});
  const [activeChain, setActiveChain] = useState<string>(chainIds[0] || '');

  useEffect(() => {
    if (selectedResidue?.chain && sequenceChains[selectedResidue.chain]) {
      setActiveChain(selectedResidue.chain);
    } else if (!sequenceChains[activeChain] && chainIds.length > 0) {
      setActiveChain(chainIds[0]);
    }
  }, [selectedResidue?.chain, chainIds.join(',')]);

  if (chainIds.length === 0) return null;
  const residues = sequenceChains[activeChain] || [];

  const toneClass = (sev: SequenceResidue['severity'], isSel: boolean) => {
    if (isSel) return 'bg-[#007AFF] text-white ring-2 ring-[#007AFF]/40 scale-110 z-10';
    if (sev === 'error') return 'bg-[#FF3B30]/20 text-[#C92519] dark:text-[#FF453A] hover:bg-[#FF3B30]/30';
    if (sev === 'warn') return 'bg-[#FF9500]/20 text-[#B26200] dark:text-[#FF9F0A] hover:bg-[#FF9500]/30';
    if (sev === 'metal') return 'bg-[#AF52DE]/20 text-[#8944AB] dark:text-[#BF5AF2] hover:bg-[#AF52DE]/30';
    if (sev === 'pocket') return 'bg-[#5AC8FA]/25 text-[#0071A4] dark:text-[#5AC8FA] hover:bg-[#5AC8FA]/40';
    return 'bg-black/[0.04] dark:bg-white/[0.06] text-[#1D1D1F] dark:text-[#E5E5EA] hover:bg-[#007AFF]/15';
  };

  return (
    <div className="apple-card px-3 py-2 shrink-0 space-y-1.5 select-none">
      <div className="flex items-center justify-between gap-2">
        <div className="flex items-center gap-1 overflow-x-auto">
          <span className="text-[10.5px] font-semibold uppercase tracking-wider text-[#86868B] mr-1">
            Sequence
          </span>
          {chainIds.map((cid) => (
            <button
              key={cid}
              onClick={() => setActiveChain(cid)}
              className={`px-2 py-0.5 rounded-lg text-[11px] font-mono font-semibold transition ${
                activeChain === cid
                  ? 'bg-[#007AFF] text-white'
                  : 'bg-black/[0.045] dark:bg-white/[0.07] text-[#6E6E73] dark:text-[#A1A1A6] hover:bg-black/10'
              }`}
            >
              {cid} ({sequenceChains[cid]?.length || 0})
            </button>
          ))}
        </div>

        <div className="hidden sm:flex items-center gap-2 text-[10px] text-[#86868B]">
          <span className="flex items-center gap-1">
            <span className="w-2 h-2 rounded-full bg-[#FF3B30]" /> Issue
          </span>
          <span className="flex items-center gap-1">
            <span className="w-2 h-2 rounded-full bg-[#FF9500]" /> Warn
          </span>
          <span className="flex items-center gap-1">
            <span className="w-2 h-2 rounded-full bg-[#5AC8FA]" /> Pocket
          </span>
          <span className="flex items-center gap-1">
            <span className="w-2 h-2 rounded-full bg-[#AF52DE]" /> Metal
          </span>
        </div>
      </div>

      {/* Horizontal 1-letter Sequence Ribbon */}
      <div className="flex items-end gap-[2px] overflow-x-auto pb-1 pt-2.5">
        {residues.map((r, idx) => {
          const isSel =
            selectedResidue?.chain === r.chain && selectedResidue?.resseq === r.resseq;
          const showTick = idx === 0 || r.resseq % 10 === 0;
          return (
            <div key={r.key} className="relative flex flex-col items-center shrink-0">
              {showTick && (
                <span className="absolute -top-3 text-[8.5px] font-mono text-[#86868B] leading-none">
                  {r.resseq}
                </span>
              )}
              <button
                onClick={() =>
                  onSelectResidue(
                    { chain: r.chain, resseq: r.resseq, icode: r.icode, resname: r.resname, key: r.key },
                    true
                  )
                }
                title={`${r.resname} ${r.key}${r.issues.length ? ` — ${r.issues.join(', ')}` : ''}`}
                className={`w-[20px] h-[22px] rounded-[5px] text-[11px] font-mono font-bold transition ${toneClass(
                  r.severity,
                  isSel
                )}`}
              >
                {r.aa1}
              </button>
            </div>
          );
        })}
      </div>
    </div>
  );
};
