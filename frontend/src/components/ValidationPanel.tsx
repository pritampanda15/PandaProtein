import React, { useState } from 'react';
import { Activity, GitCompare, ShieldAlert, Sparkles } from 'lucide-react';
import { api } from '../api/client';
import type { SelectedResidue, StructureSessionState } from '../types/pdbprep';

interface ValidationPanelProps {
  session: StructureSessionState;
  busy: boolean;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
  onStateUpdate: (state: StructureSessionState, toast: string) => void;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
}

export const ValidationPanel: React.FC<ValidationPanelProps> = ({
  session,
  busy,
  onSelectResidue,
  onStateUpdate,
  onSetBusy,
  onError,
}) => {
  const { validation } = session;
  const [ramaFilter, setRamaFilter] = useState<'all' | 'general' | 'glycine' | 'proline'>('all');
  const [hoveredKey, setHoveredKey] = useState<string | null>(null);

  const points = validation.ramachandran.points.filter((p) =>
    ramaFilter === 'all' ? true : p.category === ramaFilter
  );

  // Ramachandran SVG coordinates (-180..180 -> 36..296)
  const size = 320;
  const pad = 34;
  const plot = size - pad * 2;
  const rx = (phi: number) => pad + ((phi + 180) / 360) * plot;
  const ry = (psi: number) => pad + ((180 - psi) / 360) * plot;

  // B-factor profile SVG
  const bprof = validation.bfactor_profile;
  const maxB = Math.max(50, ...bprof.map((b) => b.mean_b));
  const bPath = bprof
    .map((b, i) => {
      const x = 30 + (i / Math.max(bprof.length - 1, 1)) * 450;
      const y = 105 - (b.mean_b / maxB) * 90;
      return `${i === 0 ? 'M' : 'L'} ${x.toFixed(1)} ${y.toFixed(1)}`;
    })
    .join(' ');

  const runQuickFix = async (label: string, fn: () => Promise<{ result: { message?: string }; state: StructureSessionState }>) => {
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
      {/* Header */}
      <div className="apple-card p-4 flex flex-wrap items-center justify-between gap-3">
        <div>
          <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#007AFF] dark:text-[#0A84FF]">
            COOT & MolProbity Crystallographic Validation
          </div>
          <h2 className="text-[16px] font-semibold tracking-[-0.015em]">
            Ramachandran Plot, Engh & Huber Geometry, Cis-Peptides & Chirality
          </h2>
        </div>
        <button
          disabled={busy}
          onClick={() =>
            runQuickFix('Refining geometry & Ramachandran outliers…', async () => {
              const r = await api.cootAutoFix(session.session_id, {
                fix_clashes: false,
                fix_rotamers: false,
                fix_nqh_flips: false,
                fix_cis_nonpro: true,
                fix_ramachandran: true,
                regularize_outliers: true,
              });
              return { result: { message: `Geometry polished (Score ${r.result.before_score} → ${r.result.after_score})` }, state: r.state };
            })
          }
          className="apple-btn-primary"
        >
          <Sparkles className="w-3.5 h-3.5" />
          <span>Auto-Refine Geometry</span>
        </button>
      </div>

      {/* KPI Row */}
      <div className="grid grid-cols-2 sm:grid-cols-4 gap-2.5">
        <div className="apple-card p-3">
          <div className="text-[11px] text-[#86868B]">Rama Favored</div>
          <div className="text-[18px] font-bold font-mono text-[#34C759]">
            {validation.ramachandran.favored_pct}%
          </div>
          <div className="text-[10.5px] text-[#86868B]">
            {validation.ramachandran.counts.outlier} outlier(s)
          </div>
        </div>
        <div className="apple-card p-3">
          <div className="text-[11px] text-[#86868B]">Bond RMSZ</div>
          <div className={`text-[18px] font-bold font-mono ${validation.geometry.bond_rmsz > 1.5 ? 'text-[#FF9500]' : ''}`}>
            {validation.geometry.bond_rmsz}
          </div>
          <div className="text-[10.5px] text-[#86868B]">
            {validation.geometry.bond_outliers.length} &gt;4σ outlier(s)
          </div>
        </div>
        <div className="apple-card p-3">
          <div className="text-[11px] text-[#86868B]">Angle RMSZ</div>
          <div className={`text-[18px] font-bold font-mono ${validation.geometry.angle_rmsz > 1.5 ? 'text-[#FF9500]' : ''}`}>
            {validation.geometry.angle_rmsz}
          </div>
          <div className="text-[10.5px] text-[#86868B]">
            {validation.geometry.angle_outliers.length} &gt;4σ outlier(s)
          </div>
        </div>
        <div className="apple-card p-3">
          <div className="text-[11px] text-[#86868B]">Cis / Chiral Issues</div>
          <div
            className={`text-[18px] font-bold font-mono ${
              validation.geometry.cis_peptides.filter((c) => c.severity === 'outlier').length +
                validation.geometry.chiral_outliers.length >
              0
                ? 'text-[#FF3B30]'
                : 'text-[#34C759]'
            }`}
          >
            {validation.geometry.cis_peptides.filter((c) => c.severity === 'outlier').length +
              validation.geometry.chiral_outliers.length}
          </div>
          <div className="text-[10.5px] text-[#86868B]">
            {validation.geometry.cbeta_outliers.length} Cβ dev
          </div>
        </div>
      </div>

      {/* Interactive SVG Ramachandran Plot */}
      <div className="apple-card p-4 grid grid-cols-1 lg:grid-cols-12 gap-4 items-center">
        <div className="lg:col-span-7 flex flex-col items-center">
          <div className="w-full flex items-center justify-between mb-2">
            <span className="text-[13px] font-semibold">Interactive Ramachandran Plot</span>
            <div className="apple-segmented">
              {(['all', 'general', 'glycine', 'proline'] as const).map((k) => (
                <button
                  key={k}
                  onClick={() => setRamaFilter(k)}
                  className={`apple-segmented-item capitalize !text-[11px] !px-2 ${
                    ramaFilter === k ? 'apple-segmented-active' : 'text-[#6E6E73]'
                  }`}
                >
                  {k}
                </button>
              ))}
            </div>
          </div>

          <svg
            viewBox={`0 0 ${size} ${size}`}
            className="w-full max-w-[330px] rounded-2xl bg-black/[0.025] dark:bg-white/[0.03] border border-black/[0.06] dark:border-white/[0.08]"
          >
            {/* Reference favored regions (alpha-R, beta, alpha-L) */}
            <ellipse cx={rx(-64)} cy={ry(-41)} rx={26} ry={22} fill="#007AFF" opacity="0.14" />
            <ellipse cx={rx(-118)} cy={ry(132)} rx={42} ry={32} fill="#007AFF" opacity="0.14" />
            <ellipse cx={rx(58)} cy={ry(38)} rx={18} ry={20} fill="#007AFF" opacity="0.11" />
            <ellipse cx={rx(-75)} cy={ry(-32)} rx={44} ry={40} fill="#007AFF" opacity="0.06" />
            <ellipse cx={rx(-115)} cy={ry(132)} rx={56} ry={48} fill="#007AFF" opacity="0.06" />

            {/* Crosshairs at 0,0 */}
            <line x1={rx(0)} y1={pad} x2={rx(0)} y2={size - pad} stroke="currentColor" className="text-black/15 dark:text-white/15" strokeDasharray="3 3" />
            <line x1={pad} y1={ry(0)} x2={size - pad} y2={ry(0)} stroke="currentColor" className="text-black/15 dark:text-white/15" strokeDasharray="3 3" />

            {/* Frame */}
            <rect x={pad} y={pad} width={plot} height={plot} fill="none" stroke="currentColor" className="text-black/15 dark:text-white/15" />

            {/* Residue points */}
            {points.map((pt) => {
              const cx = rx(pt.phi);
              const cy = ry(pt.psi);
              const fill =
                pt.status === 'outlier' ? '#FF3B30' : pt.status === 'allowed' ? '#FF9500' : '#34C759';
              const isHov = hoveredKey === pt.key;
              return (
                <circle
                  key={pt.key}
                  cx={cx}
                  cy={cy}
                  r={pt.status === 'outlier' ? 4.8 : isHov ? 4.5 : 3.1}
                  fill={fill}
                  stroke={pt.status === 'outlier' || isHov ? '#FFFFFF' : 'none'}
                  strokeWidth="1.2"
                  className="cursor-pointer transition-all"
                  onMouseEnter={() => setHoveredKey(pt.key)}
                  onMouseLeave={() => setHoveredKey(null)}
                  onClick={() =>
                    onSelectResidue(
                      { chain: pt.chain, resseq: pt.resseq, icode: pt.icode, resname: pt.resname, key: pt.key },
                      true
                    )
                  }
                >
                  <title>{`${pt.resname} ${pt.key} (φ=${pt.phi}°, ψ=${pt.psi}°, ${pt.status})`}</title>
                </circle>
              );
            })}

            {/* Axis Labels */}
            <text x={size / 2} y={size - 8} textAnchor="middle" fontSize="10" className="fill-[#86868B]" fontFamily="monospace">
              φ (deg): -180° to +180°
            </text>
            <text
              x={12}
              y={size / 2}
              textAnchor="middle"
              fontSize="10"
              className="fill-[#86868B]"
              fontFamily="monospace"
              transform={`rotate(-90 12 ${size / 2})`}
            >
              ψ (deg)
            </text>
          </svg>
        </div>

        {/* Outlier List beside Plot */}
        <div className="lg:col-span-5 space-y-2">
          <div className="text-[12.5px] font-semibold">Ramachandran Outliers</div>
          {validation.ramachandran.points.filter((p) => p.status === 'outlier').length === 0 ? (
            <div className="p-3 rounded-xl bg-[#34C759]/10 text-[#1B8A38] dark:text-[#30D158] text-[12px] font-medium">
              Zero Ramachandran outliers! All backbone torsions are in favored/allowed regions.
            </div>
          ) : (
            <div className="space-y-1.5 max-h-64 overflow-y-auto">
              {validation.ramachandran.points
                .filter((p) => p.status === 'outlier')
                .map((p) => (
                  <div
                    key={p.key}
                    className="flex items-center justify-between p-2.5 rounded-xl bg-[#FF3B30]/10 border border-[#FF3B30]/20 text-[12px]"
                  >
                    <button
                      onClick={() =>
                        onSelectResidue(
                          { chain: p.chain, resseq: p.resseq, icode: p.icode, resname: p.resname, key: p.key },
                          true
                        )
                      }
                      className="text-left"
                    >
                      <div className="font-mono font-semibold text-[#C92519] dark:text-[#FF453A]">
                        {p.resname} {p.key}
                      </div>
                      <div className="text-[11px] font-mono text-[#86868B]">
                        φ={p.phi}° · ψ={p.psi}° ({p.category})
                      </div>
                    </button>
                    <button
                      disabled={busy}
                      onClick={() =>
                        runQuickFix(`Regularizing ${p.key}…`, () =>
                          api.cootRegularize(session.session_id, {
                            chain: p.chain,
                            start_resseq: Math.max(1, p.resseq - 2),
                            end_resseq: p.resseq + 2,
                          })
                        )
                      }
                      className="apple-btn-secondary !py-1 !px-2 text-[11px]"
                    >
                      Regularize
                    </button>
                  </div>
                ))}
            </div>
          )}
        </div>
      </div>

      {/* Cis / Twisted Peptides & Chirality */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
        <div className="apple-card p-4 space-y-2.5">
          <div className="flex items-center gap-2 text-[13px] font-semibold">
            <GitCompare className="w-4 h-4 text-[#FF9500]" />
            <span>Cis & Twisted Peptides ({validation.geometry.cis_peptides.length})</span>
          </div>
          {validation.geometry.cis_peptides.length === 0 ? (
            <div className="text-[12px] text-[#34C759]">All peptide bonds are planar trans (ω ≈ 180°).</div>
          ) : (
            <div className="space-y-1.5">
              {validation.geometry.cis_peptides.map((cp) => (
                <div
                  key={cp.key}
                  className="flex items-center justify-between p-2.5 rounded-xl bg-black/[0.03] dark:bg-white/[0.04] text-[12px]"
                >
                  <button
                    onClick={() =>
                      onSelectResidue({ chain: cp.chain, resseq: cp.resseq, icode: cp.icode, resname: cp.resname }, true)
                    }
                    className="text-left"
                  >
                    <div className="font-mono font-semibold">
                      {cp.prev_resname}{cp.prev_resseq} → {cp.resname}{cp.resseq} ({cp.chain})
                    </div>
                    <div className="text-[11px] text-[#86868B]">
                      {cp.type} · ω = {cp.omega}°
                    </div>
                  </button>
                  <button
                    disabled={busy}
                    onClick={() =>
                      runQuickFix(`Converting ${cp.key} to trans…`, () =>
                        api.cootCisTrans(session.session_id, {
                          chain: cp.chain,
                          resseq: cp.resseq,
                          icode: cp.icode,
                        })
                      )
                    }
                    className="apple-btn-primary !py-1 !px-2.5 text-[11px]"
                  >
                    Flip to Trans
                  </button>
                </div>
              ))}
            </div>
          )}
        </div>

        <div className="apple-card p-4 space-y-2.5">
          <div className="flex items-center gap-2 text-[13px] font-semibold">
            <ShieldAlert className="w-4 h-4 text-[#FF3B30]" />
            <span>Chirality & Bond Outliers</span>
          </div>
          {validation.geometry.chiral_outliers.length === 0 && validation.geometry.bond_outliers.length === 0 ? (
            <div className="text-[12px] text-[#34C759]">
              Zero chirality inversions and zero &gt;4σ Engh & Huber bond outliers.
            </div>
          ) : (
            <div className="space-y-1.5 max-h-48 overflow-y-auto">
              {validation.geometry.chiral_outliers.map((ch) => (
                <div key={ch.key} className="flex items-center justify-between p-2 rounded-xl bg-[#FF3B30]/10 text-[12px]">
                  <div>
                    <div className="font-mono font-semibold">{ch.resname} {ch.key}</div>
                    <div className="text-[11px] text-[#86868B]">{ch.diagnosis} (vol {ch.chiral_volume})</div>
                  </div>
                  <button
                    disabled={busy}
                    onClick={() =>
                      runQuickFix(`Fixing chirality on ${ch.key}…`, () =>
                        api.cootFixChirality(session.session_id, {
                          chain: ch.chain,
                          resseq: ch.resseq,
                          icode: ch.icode,
                        })
                      )
                    }
                    className="apple-btn-primary !py-1 !px-2.5 text-[11px]"
                  >
                    Fix L-AA
                  </button>
                </div>
              ))}
              {validation.geometry.bond_outliers.slice(0, 10).map((bo, idx) => (
                <button
                  key={`${bo.key}-${idx}`}
                  onClick={() =>
                    onSelectResidue({ chain: bo.chain, resseq: bo.resseq, icode: bo.icode, resname: bo.resname }, true)
                  }
                  className="w-full flex items-center justify-between px-2.5 py-1.5 rounded-xl bg-[#FF9500]/10 hover:bg-[#FF9500]/16 text-[11.5px] font-mono text-left"
                >
                  <span>
                    {bo.resname} {bo.key} ({bo.atoms})
                  </span>
                  <span className="text-[#B26200] dark:text-[#FF9F0A] font-semibold">
                    {bo.observed} Å (ideal {bo.ideal}, Z={bo.z_score})
                  </span>
                </button>
              ))}
            </div>
          )}
        </div>
      </div>

      {/* B-Factor / pLDDT Profile */}
      <div className="apple-card p-4 space-y-2">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2 text-[13px] font-semibold">
            <Activity className="w-4 h-4 text-[#007AFF]" />
            <span>Per-Residue B-Factor / Flexibility Profile</span>
          </div>
          <span className="text-[11px] font-mono text-[#86868B]">Max {Math.round(maxB)} Å²</span>
        </div>
        <svg viewBox="0 0 490 120" className="w-full h-28">
          <path d={bPath} fill="none" stroke="#007AFF" strokeWidth="2" strokeLinecap="round" />
          <text x="24" y="20" textAnchor="end" fontSize="9" className="fill-[#86868B]" fontFamily="monospace">
            {Math.round(maxB)}
          </text>
          <text x="24" y="106" textAnchor="end" fontSize="9" className="fill-[#86868B]" fontFamily="monospace">
            0
          </text>
        </svg>
      </div>
    </div>
  );
};
