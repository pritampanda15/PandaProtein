import React from 'react';
import {
  Download,
  FolderOpen,
  Moon,
  Redo2,
  Search,
  Sparkles,
  Sun,
  Undo2,
  Loader2,
} from 'lucide-react';
import type { StructureSessionState } from '../types/pdbprep';

interface TopBarProps {
  session: StructureSessionState | null;
  darkMode: boolean;
  busy: string | null;
  onToggleTheme: () => void;
  onOpenImport: () => void;
  onOpenCommandPalette: () => void;
  onUndo: () => void;
  onRedo: () => void;
  onQuickDeepPrep: () => void;
  onGoExport: () => void;
}

export const TopBar: React.FC<TopBarProps> = ({
  session,
  darkMode,
  busy,
  onToggleTheme,
  onOpenImport,
  onOpenCommandPalette,
  onUndo,
  onRedo,
  onQuickDeepPrep,
  onGoExport,
}) => {
  const overall = session?.readiness.overall_score ?? 0;
  const scoreColor =
    overall >= 85
      ? 'bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158]'
      : overall >= 60
      ? 'bg-[#FF9500]/15 text-[#B26200] dark:text-[#FF9F0A]'
      : 'bg-[#FF3B30]/15 text-[#C92519] dark:text-[#FF453A]';
  const dotColor =
    overall >= 85 ? 'bg-[#34C759]' : overall >= 60 ? 'bg-[#FF9500]' : 'bg-[#FF3B30]';

  return (
    <header className="h-14 shrink-0 px-4 flex items-center justify-between apple-glass border-x-0 border-t-0 z-30 select-none">
      {/* Left: Brand & Active Structure */}
      <div className="flex items-center gap-3 min-w-0">
        <div className="flex items-center gap-2.5">
          <img
            src="/logo.png"
            alt="PandaProtein"
            className="w-8 h-8 rounded-[9px] shadow-sm border border-black/[0.06] dark:border-white/10"
          />
          <div className="leading-tight">
            <div className="text-[14px] font-semibold tracking-[-0.015em]">PandaProtein</div>
            <div className="text-[10.5px] text-[#86868B] font-medium">
              PDBFixer · PROPKA 3 · COOT · MOLE/HOLE
            </div>
          </div>
        </div>

        {session && (
          <>
            <div className="h-5 w-px bg-black/10 dark:bg-white/10 mx-1" />
            <div className="flex items-center gap-2 min-w-0">
              <span className="text-[13px] font-medium truncate max-w-[200px]">
                {session.name}
              </span>
              <span className="apple-pill bg-black/[0.05] dark:bg-white/[0.08] text-[#6E6E73] dark:text-[#A1A1A6] truncate max-w-[180px]">
                {session.source}
              </span>
              {session.metadata.resolution && (
                <span className="apple-pill bg-[#007AFF]/10 text-[#007AFF] dark:text-[#0A84FF] font-mono">
                  {session.metadata.resolution} Å
                </span>
              )}
            </div>
          </>
        )}
      </div>

      {/* Center / Status */}
      {busy && (
        <div className="hidden md:flex items-center gap-2 px-3 py-1 rounded-full bg-[#007AFF]/10 text-[#007AFF] dark:text-[#0A84FF] text-[12px] font-medium">
          <Loader2 className="w-3.5 h-3.5 animate-spin" />
          <span>{busy}</span>
        </div>
      )}

      {/* Right: Actions */}
      <div className="flex items-center gap-2">
        <button
          onClick={onOpenCommandPalette}
          title="Quick Jump & Actions (⌘K)"
          className="hidden md:flex items-center gap-2 px-2.5 py-1.5 rounded-xl bg-black/[0.045] hover:bg-black/[0.08] dark:bg-white/[0.07] dark:hover:bg-white/[0.12] text-[12px] text-[#6E6E73] dark:text-[#A1A1A6] transition"
        >
          <Search className="w-3.5 h-3.5" />
          <span className="hidden xl:inline">Jump…</span>
          <kbd className="px-1 py-0.5 rounded text-[10px] font-mono bg-white/70 dark:bg-black/30">⌘K</kbd>
        </button>

        {session && (
          <>
            {/* Readiness pill */}
            <div className={`hidden sm:flex items-center gap-1.5 px-2.5 py-1 rounded-full text-[12px] font-semibold ${scoreColor}`}>
              <span className={`w-2 h-2 rounded-full ${dotColor}`} />
              <span>Readiness {overall}%</span>
            </div>

            {/* COOT Undo / Redo */}
            <div className="flex items-center bg-black/[0.045] dark:bg-white/[0.07] rounded-xl p-0.5">
              <button
                onClick={onUndo}
                disabled={!session.can_undo || !!busy}
                title={`Undo (${session.undo_count})`}
                className="p-1.5 rounded-[9px] hover:bg-white dark:hover:bg-white/10 disabled:opacity-35 transition"
              >
                <Undo2 className="w-3.5 h-3.5" />
              </button>
              <button
                onClick={onRedo}
                disabled={!session.can_redo || !!busy}
                title={`Redo (${session.redo_count})`}
                className="p-1.5 rounded-[9px] hover:bg-white dark:hover:bg-white/10 disabled:opacity-35 transition"
              >
                <Redo2 className="w-3.5 h-3.5" />
              </button>
            </div>

            {/* One-Click Deep Prep */}
            <button
              onClick={onQuickDeepPrep}
              disabled={!!busy}
              className="apple-btn-primary"
              title="Run full PDBFixer + COOT Polish + PROPKA Protonation pipeline"
            >
              <Sparkles className="w-3.5 h-3.5" />
              <span className="hidden lg:inline">Deep Prep</span>
            </button>

            <button onClick={onGoExport} className="apple-btn-secondary" title="Export MD & Docking Files">
              <Download className="w-3.5 h-3.5" />
              <span className="hidden xl:inline">Export</span>
            </button>
          </>
        )}

        <button onClick={onOpenImport} className="apple-btn-secondary">
          <FolderOpen className="w-3.5 h-3.5" />
          <span>Load</span>
        </button>

        <button
          onClick={onToggleTheme}
          className="p-2 rounded-xl bg-black/[0.045] hover:bg-black/[0.08] dark:bg-white/[0.07] dark:hover:bg-white/[0.12] transition"
          title={darkMode ? 'Switch to Light Mode' : 'Switch to Dark Mode'}
        >
          {darkMode ? <Sun className="w-4 h-4 text-[#FF9F0A]" /> : <Moon className="w-4 h-4 text-[#5856D6]" />}
        </button>
      </div>
    </header>
  );
};
