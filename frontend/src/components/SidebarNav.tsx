import React from 'react';
import {
  Activity,
  AlertTriangle,
  Beaker,
  Compass,
  Download,
  Droplets,
  FileEdit,
  FolderKanban,
  GitCompare,
  Layers,
  Microscope,
  Waypoints,
  Wrench,
} from 'lucide-react';
import type { NavTab, StructureSessionState } from '../types/pdbprep';

interface SidebarNavProps {
  activeTab: NavTab;
  onSelect: (tab: NavTab) => void;
  session: StructureSessionState | null;
  moleCount?: number;
}

export const SidebarNav: React.FC<SidebarNavProps> = ({ activeTab, onSelect, session, moleCount }) => {
  const fixerIssues = session
    ? session.pdbfixer.nonstandard_residues.length +
      session.pdbfixer.missing_atoms.length +
      session.pdbfixer.missing_loops.length +
      session.pdbfixer.altlocs.length
    : 0;

  const valIssues = session
    ? session.validation.ramachandran.counts.outlier +
      session.validation.geometry.cis_peptides.filter((c) => c.severity === 'outlier').length +
      session.validation.geometry.chiral_outliers.length
    : 0;

  const clashRotIssues = session
    ? session.validation.clashes.total_clashes + session.validation.rotamers.counts.outlier
    : 0;

  const sections: Array<{
    title: string;
    items: Array<{ id: NavTab; label: string; icon: React.ReactNode; badge?: number | string; badgeTone?: string }>;
  }> = [
    {
      title: 'Workspace',
      items: [
        { id: 'overview', label: 'Overview & Audit', icon: <Compass className="w-4 h-4" /> },
        { id: 'coot', label: 'COOT Workbench', icon: <Microscope className="w-4 h-4" /> },
        {
          id: 'chains',
          label: 'Chains & Renaming',
          icon: <FileEdit className="w-4 h-4" />,
          badge: session ? `${session.metadata.chain_count} ch` : undefined,
          badgeTone: 'bg-[#5856D6]/15 text-[#5856D6] dark:text-[#7D7AFF]',
        },
      ],
    },
    {
      title: 'Deep Preparation',
      items: [
        {
          id: 'pdbfixer',
          label: 'PDBFixer Pipeline',
          icon: <Wrench className="w-4 h-4" />,
          badge: fixerIssues > 0 ? fixerIssues : undefined,
          badgeTone: 'bg-[#FF9500]/15 text-[#B26200] dark:text-[#FF9F0A]',
        },
        {
          id: 'propka',
          label: 'PROPKA 3 & pH',
          icon: <Droplets className="w-4 h-4" />,
          badge: session ? `pH ${session.propka.ph}` : undefined,
          badgeTone: 'bg-[#007AFF]/12 text-[#007AFF] dark:text-[#0A84FF]',
        },
      ],
    },
    {
      title: 'Crystallographic QA',
      items: [
        {
          id: 'validation',
          label: 'Ramachandran & Geom',
          icon: <Activity className="w-4 h-4" />,
          badge: valIssues > 0 ? valIssues : undefined,
          badgeTone: 'bg-[#FF3B30]/15 text-[#D70015] dark:text-[#FF453A]',
        },
        {
          id: 'clashes',
          label: 'Clashes & Rotamers',
          icon: <AlertTriangle className="w-4 h-4" />,
          badge: clashRotIssues > 0 ? clashRotIssues : undefined,
          badgeTone: 'bg-[#FF3B30]/15 text-[#D70015] dark:text-[#FF453A]',
        },
        {
          id: 'ligands',
          label: 'Ligands & Metals',
          icon: <Beaker className="w-4 h-4" />,
          badge: session ? session.ligands_info.ligands.length + session.ligands_info.metals.length : undefined,
          badgeTone: 'bg-[#AF52DE]/15 text-[#8944AB] dark:text-[#BF5AF2]',
        },
        {
          id: 'tunnels',
          label: 'Tunnels & Pores',
          icon: <Waypoints className="w-4 h-4" />,
          badge: moleCount && moleCount > 0 ? moleCount : 'MOLE',
          badgeTone: 'bg-[#007AFF]/15 text-[#007AFF] dark:text-[#0A84FF]',
        },
        {
          id: 'analysis',
          label: 'Compare & Nets',
          icon: <GitCompare className="w-4 h-4" />,
          badge: 'RMSD',
          badgeTone: 'bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158]',
        },
      ],
    },
    {
      title: 'Simulation & Docking',
      items: [
        { id: 'export', label: 'MD & Docking Export', icon: <Download className="w-4 h-4" /> },
        {
          id: 'batch',
          label: 'Batch Prep Queue',
          icon: <FolderKanban className="w-4 h-4" />,
          badge: 'ZIP',
          badgeTone: 'bg-[#5856D6]/15 text-[#5856D6] dark:text-[#7D7AFF]',
        },
      ],
    },
  ];

  return (
    <aside className="w-[224px] shrink-0 apple-glass border-y-0 border-l-0 flex flex-col justify-between p-3 select-none overflow-y-auto">
      <nav className="space-y-5">
        {sections.map((sec) => (
          <div key={sec.title}>
            <div className="px-2.5 mb-1.5 text-[10.5px] font-semibold uppercase tracking-[0.06em] text-[#86868B]">
              {sec.title}
            </div>
            <div className="space-y-0.5">
              {sec.items.map((item) => {
                const active = activeTab === item.id;
                return (
                  <button
                    key={item.id}
                    onClick={() => onSelect(item.id)}
                    className={`w-full flex items-center justify-between px-2.5 py-2 rounded-xl text-[13px] font-medium transition-all duration-150 ${
                      active
                        ? 'bg-[#007AFF] text-white shadow-sm'
                        : 'text-[#1D1D1F] dark:text-[#E5E5EA] hover:bg-black/[0.045] dark:hover:bg-white/[0.07]'
                    }`}
                  >
                    <span className="flex items-center gap-2.5 truncate">
                      <span className={active ? 'text-white' : 'text-[#6E6E73] dark:text-[#98989D]'}>
                        {item.icon}
                      </span>
                      <span className="truncate">{item.label}</span>
                    </span>
                    {item.badge !== undefined && (
                      <span
                        className={`px-1.5 py-0.5 rounded-full text-[10.5px] font-semibold leading-none ${
                          active
                            ? 'bg-white/25 text-white'
                            : item.badgeTone || 'bg-black/10 text-[#6E6E73]'
                        }`}
                      >
                        {item.badge}
                      </span>
                    )}
                  </button>
                );
              })}
            </div>
          </div>
        ))}
      </nav>

      {/* Bottom Quick Telemetry Card */}
      {session && (
        <div className="mt-4 p-3 rounded-2xl bg-black/[0.035] dark:bg-white/[0.05] border border-black/[0.05] dark:border-white/[0.07] space-y-2">
          <div className="flex items-center gap-1.5 text-[11px] font-semibold text-[#6E6E73] dark:text-[#98989D]">
            <Layers className="w-3.5 h-3.5" />
            <span>Structure Telemetry</span>
          </div>
          <div className="grid grid-cols-2 gap-1.5 text-[11.5px]">
            <div>
              <div className="text-[#86868B] text-[10px]">Residues</div>
              <div className="font-semibold font-mono">{session.metadata.counts.protein_residues}</div>
            </div>
            <div>
              <div className="text-[#86868B] text-[10px]">Atoms (H)</div>
              <div className="font-semibold font-mono">
                {session.metadata.counts.atoms} ({session.metadata.counts.hydrogens})
              </div>
            </div>
            <div>
              <div className="text-[#86868B] text-[10px]">Clashscore</div>
              <div className="font-semibold font-mono">{session.validation.clashes.clashscore}</div>
            </div>
            <div>
              <div className="text-[#86868B] text-[10px]">Net Charge</div>
              <div className="font-semibold font-mono">
                {session.propka.net_charge_discrete >= 0 ? '+' : ''}
                {session.propka.net_charge_discrete}
              </div>
            </div>
          </div>
        </div>
      )}
    </aside>
  );
};
