import React from 'react';
import {
  AlertCircle,
  CheckCircle2,
  Clock,
  Cpu,
  Crosshair,
  Microscope,
  Sparkles,
  XCircle,
} from 'lucide-react';
import type { NavTab, StructureSessionState } from '../types/pdbprep';

interface OverviewPanelProps {
  session: StructureSessionState;
  busy: boolean;
  onRunPreset: (preset: 'md_production' | 'docking_receptor' | 'crystallographic_polish') => void;
  onNavigate: (tab: NavTab) => void;
}

const ScoreRing: React.FC<{ label: string; score: number; subtitle: string; color: string; onClick?: () => void }> = ({
  label,
  score,
  subtitle,
  color,
  onClick,
}) => {
  const r = 34;
  const circ = 2 * Math.PI * r;
  const offset = circ - (Math.min(100, Math.max(0, score)) / 100) * circ;

  return (
    <button
      onClick={onClick}
      className="apple-card p-4 flex items-center gap-3.5 text-left hover:scale-[1.01] transition-transform"
    >
      <div className="relative w-20 h-20 shrink-0">
        <svg className="w-20 h-20 -rotate-90" viewBox="0 0 80 80">
          <circle
            cx="40"
            cy="40"
            r={r}
            fill="none"
            stroke="currentColor"
            className="text-black/[0.06] dark:text-white/[0.08]"
            strokeWidth="7"
          />
          <circle
            cx="40"
            cy="40"
            r={r}
            fill="none"
            stroke={color}
            strokeWidth="7"
            strokeLinecap="round"
            strokeDasharray={circ}
            strokeDashoffset={offset}
            className="transition-all duration-700 ease-out"
          />
        </svg>
        <div className="absolute inset-0 flex flex-col items-center justify-center">
          <span className="text-[16px] font-bold font-mono tracking-tight">{Math.round(score)}</span>
          <span className="text-[9.5px] text-[#86868B] font-medium">/ 100</span>
        </div>
      </div>
      <div className="min-w-0">
        <div className="text-[13.5px] font-semibold tracking-[-0.01em] truncate">{label}</div>
        <div className="text-[11.5px] text-[#86868B] mt-0.5 leading-snug">{subtitle}</div>
        <div className="mt-1.5">
          <span
            className={`apple-pill ${
              score >= 85
                ? 'bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158]'
                : score >= 60
                ? 'bg-[#FF9500]/15 text-[#B26200] dark:text-[#FF9F0A]'
                : 'bg-[#FF3B30]/15 text-[#D70015] dark:text-[#FF453A]'
            }`}
          >
            {score >= 85 ? 'Production Ready' : score >= 60 ? 'Needs Polish' : 'Preparation Required'}
          </span>
        </div>
      </div>
    </button>
  );
};

export const OverviewPanel: React.FC<OverviewPanelProps> = ({
  session,
  busy,
  onRunPreset,
  onNavigate,
}) => {
  const { readiness, metadata, validation, propka, pdbfixer, baseline } = session;
  const currentIssues = validation.residues.filter((r) => r.issues.length > 0).length;
  const hasProgress =
    baseline &&
    (Math.abs(readiness.overall_score - baseline.overall_score) >= 0.5 ||
      currentIssues !== baseline.issues_count ||
      metadata.counts.hydrogens !== baseline.hydrogens);

  return (
    <div className="space-y-4 pb-6">
      {/* Header Title Card */}
      <div className="apple-card p-5 space-y-3">
        <div className="flex flex-wrap items-start justify-between gap-3">
          <div className="min-w-0 flex-1">
            <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#007AFF] dark:text-[#0A84FF]">
              Executive Structural Audit
            </div>
            <h1 className="text-[18px] font-semibold tracking-[-0.02em] mt-0.5 leading-snug">
              {metadata.title}
            </h1>
            <div className="flex flex-wrap items-center gap-2 mt-2 text-[12px] text-[#6E6E73] dark:text-[#98989D]">
              <span>{metadata.method}</span>
              {metadata.resolution && <span>· {metadata.resolution} Å</span>}
              {metadata.spacegroup && <span>· {metadata.spacegroup}</span>}
              <span>· {metadata.chain_count} chain(s)</span>
              <span>· pI {propka.isoelectric_point}</span>
              {session.pockets && session.pockets.length > 0 && (
                <span>· {session.pockets.length} druggable pocket(s)</span>
              )}
            </div>
          </div>
        </div>

        {/* Before -> After Live Progress Delta Strip */}
        {hasProgress && baseline && (
          <div className="pt-3 border-t border-black/[0.06] dark:border-white/[0.08] grid grid-cols-2 sm:grid-cols-4 gap-2 text-[11.5px]">
            <div className="px-3 py-2 rounded-xl bg-[#34C759]/[0.08] border border-[#34C759]/20">
              <div className="text-[10px] text-[#6E6E73] dark:text-[#A1A1A6]">Overall Readiness</div>
              <div className="font-mono font-semibold text-[#1B8A38] dark:text-[#30D158]">
                {baseline.overall_score}% → {readiness.overall_score}%
              </div>
            </div>
            <div className="px-3 py-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#6E6E73] dark:text-[#A1A1A6]">Open COOT Issues</div>
              <div className="font-mono font-semibold">
                {baseline.issues_count} → {currentIssues}
              </div>
            </div>
            <div className="px-3 py-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#6E6E73] dark:text-[#A1A1A6]">Clashscore</div>
              <div className="font-mono font-semibold">
                {baseline.clashscore} → {validation.clashes.clashscore}
              </div>
            </div>
            <div className="px-3 py-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#6E6E73] dark:text-[#A1A1A6]">Hydrogens Built</div>
              <div className="font-mono font-semibold text-[#007AFF] dark:text-[#0A84FF]">
                {baseline.hydrogens} → {metadata.counts.hydrogens}
              </div>
            </div>
          </div>
        )}
      </div>

      {/* Three Apple Activity-Style Readiness Rings */}
      <div className="grid grid-cols-1 xl:grid-cols-3 gap-3">
        <ScoreRing
          label="MD Simulation"
          score={readiness.md.score}
          subtitle="GROMACS / OpenMM / AMBER"
          color="#007AFF"
          onClick={() => onNavigate('pdbfixer')}
        />
        <ScoreRing
          label="Molecular Docking"
          score={readiness.docking.score}
          subtitle="AutoDock Vina / GNINA / Smina"
          color="#AF52DE"
          onClick={() => onNavigate('export')}
        />
        <ScoreRing
          label="COOT Quality"
          score={readiness.crystallographic.score}
          subtitle="Ramachandran · Rotamers · Clashes"
          color="#34C759"
          onClick={() => onNavigate('validation')}
        />
      </div>

      {/* One-Click Deep Preparation Presets */}
      <div className="apple-card p-4 space-y-3">
        <div className="flex items-center justify-between">
          <div>
            <h2 className="text-[14px] font-semibold">One-Click Deep Preparation Pipelines</h2>
            <p className="text-[11.5px] text-[#86868B]">
              Chains PDBFixer structural repair + COOT crystallographic refinement + PROPKA 3 protonation
            </p>
          </div>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-3 gap-2.5">
          <button
            disabled={busy}
            onClick={() => onRunPreset('md_production')}
            className="p-3.5 rounded-2xl text-left bg-gradient-to-br from-[#007AFF]/10 to-[#5856D6]/5 border border-[#007AFF]/20 hover:border-[#007AFF]/45 transition group"
          >
            <div className="flex items-center gap-2 text-[#007AFF] dark:text-[#0A84FF] font-semibold text-[13px]">
              <Cpu className="w-4 h-4" />
              <span>MD Production Prep</span>
            </div>
            <p className="text-[11.5px] text-[#6E6E73] dark:text-[#A1A1A6] mt-1 leading-relaxed">
              Full PDBFixer (MSE, loops, atoms, S–S, bulk waters, buffers) + COOT polish + PROPKA pH {propka.ph} hydrogens.
            </p>
          </button>

          <button
            disabled={busy}
            onClick={() => onRunPreset('docking_receptor')}
            className="p-3.5 rounded-2xl text-left bg-gradient-to-br from-[#AF52DE]/10 to-[#FF2D55]/5 border border-[#AF52DE]/20 hover:border-[#AF52DE]/45 transition group"
          >
            <div className="flex items-center gap-2 text-[#AF52DE] dark:text-[#BF5AF2] font-semibold text-[13px]">
              <Crosshair className="w-4 h-4" />
              <span>Docking Receptor Prep</span>
            </div>
            <p className="text-[11.5px] text-[#6E6E73] dark:text-[#A1A1A6] mt-1 leading-relaxed">
              Cleans pocket buffers, resolves AltLocs & side-chain clashes, optimizes NQH H-bonds & polar H for PDBQT.
            </p>
          </button>

          <button
            disabled={busy}
            onClick={() => onRunPreset('crystallographic_polish')}
            className="p-3.5 rounded-2xl text-left bg-gradient-to-br from-[#34C759]/10 to-[#5AC8FA]/5 border border-[#34C759]/20 hover:border-[#34C759]/45 transition group"
          >
            <div className="flex items-center gap-2 text-[#1B8A38] dark:text-[#30D158] font-semibold text-[13px]">
              <Microscope className="w-4 h-4" />
              <span>COOT Auto-Polish</span>
            </div>
            <p className="text-[11.5px] text-[#6E6E73] dark:text-[#A1A1A6] mt-1 leading-relaxed">
              Fixes cis non-Pro peptides, chirality, NQH flips, auto-fits rotamers & regularizes outlier zones only.
            </p>
          </button>
        </div>
      </div>

      {/* Diagnostic Checklists */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-3">
        <div className="apple-card p-4 space-y-2.5">
          <div className="flex items-center justify-between">
            <h3 className="text-[13.5px] font-semibold">MD & PDBFixer Pre-Flight Checks</h3>
            <span className="text-[11px] font-mono text-[#86868B]">{readiness.md.score}%</span>
          </div>
          <div className="space-y-1.5">
            {readiness.md.checks.map((chk) => (
              <div
                key={chk.label}
                className="flex items-start gap-2.5 px-2.5 py-1.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03]"
              >
                {chk.status === 'pass' ? (
                  <CheckCircle2 className="w-4 h-4 text-[#34C759] shrink-0 mt-0.5" />
                ) : chk.status === 'warn' ? (
                  <AlertCircle className="w-4 h-4 text-[#FF9500] shrink-0 mt-0.5" />
                ) : (
                  <XCircle className="w-4 h-4 text-[#FF3B30] shrink-0 mt-0.5" />
                )}
                <div className="min-w-0">
                  <div className="text-[12.5px] font-medium leading-tight">{chk.label}</div>
                  <div className="text-[11px] text-[#86868B] mt-0.5">{chk.detail}</div>
                </div>
              </div>
            ))}
          </div>
        </div>

        <div className="apple-card p-4 space-y-2.5">
          <div className="flex items-center justify-between">
            <h3 className="text-[13.5px] font-semibold">COOT & MolProbity Quality</h3>
            <span className="text-[11px] font-mono text-[#86868B]">{readiness.crystallographic.score}%</span>
          </div>
          <div className="space-y-1.5">
            {readiness.crystallographic.checks.map((chk) => (
              <div
                key={chk.label}
                className="flex items-start gap-2.5 px-2.5 py-1.5 rounded-xl bg-black/[0.025] dark:bg-white/[0.03]"
              >
                {chk.status === 'pass' ? (
                  <CheckCircle2 className="w-4 h-4 text-[#34C759] shrink-0 mt-0.5" />
                ) : chk.status === 'warn' ? (
                  <AlertCircle className="w-4 h-4 text-[#FF9500] shrink-0 mt-0.5" />
                ) : (
                  <XCircle className="w-4 h-4 text-[#FF3B30] shrink-0 mt-0.5" />
                )}
                <div className="min-w-0">
                  <div className="text-[12.5px] font-medium leading-tight">{chk.label}</div>
                  <div className="text-[11px] text-[#86868B] mt-0.5">{chk.detail}</div>
                </div>
              </div>
            ))}
          </div>

          {/* Quick Composition Pills */}
          <div className="pt-2 border-t border-black/[0.06] dark:border-white/[0.08] grid grid-cols-3 gap-2 text-center">
            <div className="p-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
              <div className="text-[14px] font-semibold font-mono">{pdbfixer.waters.total}</div>
              <div className="text-[10.5px] text-[#86868B]">Waters ({pdbfixer.waters.structural} struct)</div>
            </div>
            <div className="p-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
              <div className="text-[14px] font-semibold font-mono">{pdbfixer.disulfides.length}</div>
              <div className="text-[10.5px] text-[#86868B]">Disulfide Pairs</div>
            </div>
            <div className="p-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04]">
              <div className="text-[14px] font-semibold font-mono">{validation.nqh_flips.length}</div>
              <div className="text-[10.5px] text-[#86868B]">NQH Flips</div>
            </div>
          </div>
        </div>
      </div>

      {/* Audit History Trail */}
      <div className="apple-card p-4 space-y-2">
        <div className="flex items-center gap-2 text-[13px] font-semibold">
          <Clock className="w-4 h-4 text-[#007AFF]" />
          <span>COOT Checkpoint & Preparation Log</span>
        </div>
        <div className="max-h-36 overflow-y-auto space-y-1 font-mono text-[11.5px] text-[#6E6E73] dark:text-[#A1A1A6]">
          {session.history
            .slice()
            .reverse()
            .map((entry, idx) => (
              <div key={idx} className="px-2 py-1 rounded-lg bg-black/[0.02] dark:bg-white/[0.025]">
                {entry}
              </div>
            ))}
        </div>
      </div>
    </div>
  );
};
