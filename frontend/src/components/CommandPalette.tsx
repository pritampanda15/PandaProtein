import React, { useEffect, useMemo, useRef, useState } from 'react';
import {
  Beaker,
  Command,
  Compass,
  Layers,
  Microscope,
  Search,
  Sparkles,
  Wrench,
} from 'lucide-react';
import type { NavTab, SelectedResidue, StructureSessionState } from '../types/pdbprep';
import type { ChainHighlight } from './Viewer3D';

interface CommandPaletteProps {
  open: boolean;
  session: StructureSessionState | null;
  onClose: () => void;
  onNavigate: (tab: NavTab) => void;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
  onHighlightChains: (hl: ChainHighlight | null, zoom?: boolean) => void;
  onQuickAction: (action: 'deep_prep' | 'coot_refine' | 'open_import' | 'toggle_theme') => void;
}

interface CmdItem {
  id: string;
  group: string;
  title: string;
  subtitle: string;
  icon: React.ReactNode;
  run: () => void;
}

export const CommandPalette: React.FC<CommandPaletteProps> = ({
  open,
  session,
  onClose,
  onNavigate,
  onSelectResidue,
  onHighlightChains,
  onQuickAction,
}) => {
  const [query, setQuery] = useState('');
  const [cursor, setCursor] = useState(0);
  const inputRef = useRef<HTMLInputElement | null>(null);

  useEffect(() => {
    if (open) {
      setQuery('');
      setCursor(0);
      setTimeout(() => inputRef.current?.focus(), 30);
    }
  }, [open]);

  const items = useMemo<CmdItem[]>(() => {
    const q = query.trim().toLowerCase();
    const list: CmdItem[] = [
      {
        id: 'act-deep',
        group: 'Actions',
        title: 'Run Deep Prep Pipeline (PDBFixer + COOT + PROPKA)',
        subtitle: 'Full automated preparation for MD & Docking',
        icon: <Sparkles className="w-4 h-4 text-[#007AFF]" />,
        run: () => onQuickAction('deep_prep'),
      },
      {
        id: 'act-coot',
        group: 'Actions',
        title: 'Run COOT Auto-Refine All',
        subtitle: 'Fix clashes, rotamers, cis-peptides, Rama & geometry',
        icon: <Microscope className="w-4 h-4 text-[#34C759]" />,
        run: () => onQuickAction('coot_refine'),
      },
      {
        id: 'act-load',
        group: 'Actions',
        title: 'Load PDB / mmCIF / .gro or Fetch RCSB / AlphaFold',
        subtitle: 'Open Structure Import dialog',
        icon: <Compass className="w-4 h-4 text-[#5856D6]" />,
        run: () => onQuickAction('open_import'),
      },
      {
        id: 'act-theme',
        group: 'Actions',
        title: 'Toggle Light / Dark Appearance',
        subtitle: 'Switch Apple Studio theme',
        icon: <Command className="w-4 h-4 text-[#86868B]" />,
        run: () => onQuickAction('toggle_theme'),
      },
    ];

    const tabs: Array<[NavTab, string, string]> = [
      ['overview', 'Go to Overview & Readiness Audit', 'Readiness rings & Before→After diff'],
      ['coot', 'Go to COOT Workbench', 'Mutate, rotamers, pepflip, cis-trans, regularize'],
      ['chains', 'Go to Chains & Renaming', 'Rename PDB file, swap/delete/merge chains, renumber'],
      ['pdbfixer', 'Go to PDBFixer Pipeline', 'Missing atoms/loops, AltLocs, disulfides, waters'],
      ['propka', 'Go to PROPKA 3 & Protonation', 'Interactive pH titration curve & H-bond optimization'],
      ['validation', 'Go to Ramachandran & Geometry', 'Interactive Rama plot, Engh & Huber RMSZ'],
      ['clashes', 'Go to Clashes & Rotamers', 'MolProbity clashscore & NQH 180° flips'],
      ['ligands', 'Go to Ligands, Pockets & Metals', 'RDKit 2D structures, cavities & ETKDGv3 conformers'],
      ['export', 'Go to MD & Docking Export', 'Live GROMACS test, PDBQT, .gro & ZIP bundle'],
    ];
    for (const [id, title, sub] of tabs) {
      list.push({
        id: `tab-${id}`,
        group: 'Navigation',
        title,
        subtitle: sub,
        icon: <Wrench className="w-4 h-4 text-[#6E6E73]" />,
        run: () => onNavigate(id),
      });
    }

    if (session) {
      for (const ch of session.metadata.chains) {
        list.push({
          id: `ch-${ch.id}`,
          group: 'Chains',
          title: `Spotlight Chain ${ch.id} in 3D`,
          subtitle: `${ch.protein_residues} protein res · ${ch.atom_count} atoms`,
          icon: <Layers className="w-4 h-4 text-[#007AFF]" />,
          run: () => {
            onNavigate('chains');
            onHighlightChains({ ids: [ch.id], tone: 'select' }, true);
          },
        });
      }

      for (const lig of session.ligands_info.ligands) {
        list.push({
          id: `lig-${lig.id}`,
          group: 'Ligands & Pockets',
          title: `Focus Ligand ${lig.resname} (${lig.chain}:${lig.resseq})`,
          subtitle: lig.formula ? `${lig.formula} · MW ${lig.mw} Da` : `${lig.heavy_atoms} heavy atoms`,
          icon: <Beaker className="w-4 h-4 text-[#AF52DE]" />,
          run: () => {
            onNavigate('ligands');
            onSelectResidue({ chain: lig.chain, resseq: lig.resseq, icode: lig.icode, resname: lig.resname }, true);
          },
        });
      }

      // Residue search (only include matching residues when user types something)
      if (q.length >= 2) {
        let addedRes = 0;
        for (const r of session.validation.residues) {
          const hay = `${r.resname} ${r.key} ${r.chain}${r.resseq} ${r.issues.join(' ')}`.toLowerCase();
          if (hay.includes(q)) {
            list.push({
              id: `res-${r.key}`,
              group: 'Residues',
              title: `${r.resname} ${r.key}`,
              subtitle: r.issues.length ? `Issues: ${r.issues.join(', ')}` : `B=${r.mean_b} Å² · Ideal geometry`,
              icon: <Microscope className="w-4 h-4 text-[#FF9500]" />,
              run: () => {
                onNavigate('coot');
                onSelectResidue(
                  { chain: r.chain, resseq: r.resseq, icode: r.icode, resname: r.resname, key: r.key },
                  true
                );
              },
            });
            addedRes += 1;
            if (addedRes >= 18) break;
          }
        }
      }
    }

    if (!q) return list.slice(0, 14);
    return list
      .filter((it) => `${it.title} ${it.subtitle} ${it.group}`.toLowerCase().includes(q))
      .slice(0, 18);
  }, [query, session]);

  useEffect(() => {
    setCursor(0);
  }, [query]);

  if (!open) return null;

  const execute = (item?: CmdItem) => {
    if (!item) return;
    onClose();
    item.run();
  };

  return (
    <div
      onClick={onClose}
      className="fixed inset-0 z-50 bg-black/30 backdrop-blur-md flex items-start justify-center pt-[12vh] px-4"
    >
      <div
        onClick={(e) => e.stopPropagation()}
        className="w-full max-w-lg apple-card overflow-hidden shadow-apple-lg"
      >
        <div className="flex items-center gap-2.5 px-4 py-3 border-b border-black/[0.07] dark:border-white/[0.08]">
          <Search className="w-4 h-4 text-[#86868B] shrink-0" />
          <input
            ref={inputRef}
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === 'Escape') onClose();
              else if (e.key === 'ArrowDown') {
                e.preventDefault();
                setCursor((c) => Math.min(items.length - 1, c + 1));
              } else if (e.key === 'ArrowUp') {
                e.preventDefault();
                setCursor((c) => Math.max(0, c - 1));
              } else if (e.key === 'Enter') {
                e.preventDefault();
                execute(items[cursor]);
              }
            }}
            placeholder="Jump to residue (HIS 18, A:41), chain, tool, or action…"
            className="w-full bg-transparent text-[14px] focus:outline-none placeholder:text-[#86868B]"
          />
          <kbd className="px-1.5 py-0.5 rounded-md text-[10px] font-mono bg-black/[0.06] dark:bg-white/[0.08] text-[#86868B]">
            ESC
          </kbd>
        </div>

        <div className="max-h-80 overflow-y-auto p-2 space-y-0.5">
          {items.length === 0 ? (
            <div className="py-8 text-center text-[12.5px] text-[#86868B]">
              No matching residues, chains, or actions
            </div>
          ) : (
            items.map((item, idx) => (
              <button
                key={item.id}
                onMouseEnter={() => setCursor(idx)}
                onClick={() => execute(item)}
                className={`w-full flex items-center justify-between gap-3 px-3 py-2 rounded-xl text-left transition ${
                  idx === cursor
                    ? 'bg-[#007AFF] text-white'
                    : 'hover:bg-black/[0.04] dark:hover:bg-white/[0.06]'
                }`}
              >
                <div className="flex items-center gap-2.5 min-w-0">
                  <span className={idx === cursor ? 'text-white' : ''}>{item.icon}</span>
                  <div className="min-w-0">
                    <div className="text-[13px] font-medium truncate">{item.title}</div>
                    <div
                      className={`text-[11px] truncate ${
                        idx === cursor ? 'text-white/80' : 'text-[#86868B]'
                      }`}
                    >
                      {item.subtitle}
                    </div>
                  </div>
                </div>
                <span
                  className={`text-[10px] uppercase tracking-wider font-semibold shrink-0 ${
                    idx === cursor ? 'text-white/80' : 'text-[#86868B]'
                  }`}
                >
                  {item.group}
                </span>
              </button>
            ))
          )}
        </div>
      </div>
    </div>
  );
};
