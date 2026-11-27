import React, { useEffect, useRef, useState } from 'react';
import { Droplets, ShieldCheck, Zap } from 'lucide-react';
import { api } from '../api/client';
import type { PropkaData, SelectedResidue, StructureSessionState } from '../types/pdbprep';

interface PropkaPanelProps {
  session: StructureSessionState;
  busy: boolean;
  onSelectResidue: (res: SelectedResidue, focus3D?: boolean) => void;
  onStateUpdate: (state: StructureSessionState, toast: string) => void;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
}

export const PropkaPanel: React.FC<PropkaPanelProps> = ({
  session,
  busy,
  onSelectResidue,
  onStateUpdate,
  onSetBusy,
  onError,
}) => {
  const [ph, setPh] = useState<number>(session.propka.ph);
  const [respectMetal, setRespectMetal] = useState(true);
  const [respectSS, setRespectSS] = useState(true);
  const [addHydrogens, setAddHydrogens] = useState(true);
  const [optimizeHbonds, setOptimizeHbonds] = useState(true);
  const [hisDefault, setHisDefault] = useState<'auto' | 'HID' | 'HIE' | 'HIP'>('auto');
  const [overrides, setOverrides] = useState<Record<string, string>>({});
  const [liveData, setLiveData] = useState<PropkaData>(session.propka);
  const debounceRef = useRef<number | null>(null);

  useEffect(() => {
    setLiveData(session.propka);
    setPh(session.propka.ph);
  }, [session.propka]);

  useEffect(() => {
    return () => {
      if (debounceRef.current) window.clearTimeout(debounceRef.current);
    };
  }, []);

  const refreshPreview = (
    nextPh: number,
    nextOverrides: Record<string, string>,
    nextHis: 'auto' | 'HID' | 'HIE' | 'HIP',
    nextMetal: boolean,
    nextSS: boolean,
    immediate = false
  ) => {
    if (debounceRef.current) window.clearTimeout(debounceRef.current);
    const run = async () => {
      try {
        const resp = await api.previewPropka(session.session_id, {
          ph: nextPh,
          respect_metal_coordination: nextMetal,
          respect_disulfides: nextSS,
          his_default: nextHis,
          overrides: nextOverrides,
        });
        setLiveData(resp.propka);
      } catch (e: any) {
        onError(e?.message || 'PROPKA preview failed');
      }
    };
    if (immediate) {
      void run();
    } else {
      debounceRef.current = window.setTimeout(() => {
        void run();
      }, 160);
    }
  };

  const handleApply = async () => {
    onSetBusy(`Applying PROPKA 3 protonation at pH ${ph.toFixed(1)}…`);
    try {
      const resp = await api.applyPropka(session.session_id, {
        ph,
        respect_metal_coordination: respectMetal,
        respect_disulfides: respectSS,
        add_hydrogens: addHydrogens,
        optimize_hbond_network: optimizeHbonds,
        his_default: hisDefault,
        overrides,
      });
      onStateUpdate(
        resp.state,
        `Applied PROPKA protonation at pH ${ph.toFixed(1)} (+${resp.result.hydrogens_added} hydrogens)`
      );
    } catch (e: any) {
      onError(e?.message || 'Failed to apply protonation');
    } finally {
      onSetBusy(null);
    }
  };

  // SVG Titration Curve coordinates
  const curve = liveData.titration_curve || [];
  const maxQ = Math.max(6, ...curve.map((c) => Math.abs(c.net_charge)));
  const svgW = 500;
  const svgH = 140;
  const padL = 36;
  const padR = 14;
  const padT = 14;
  const padB = 24;
  const toX = (p: number) => padL + (p / 14) * (svgW - padL - padR);
  const toY = (q: number) => padT + ((maxQ - q) / (2 * maxQ)) * (svgH - padT - padB);
  const pathD = curve
    .map((pt, i) => `${i === 0 ? 'M' : 'L'} ${toX(pt.ph).toFixed(1)} ${toY(pt.net_charge).toFixed(1)}`)
    .join(' ');

  return (
    <div className="space-y-4 pb-6">
      {/* Header */}
      <div className="apple-card p-4 flex flex-wrap items-center justify-between gap-3">
        <div>
          <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#007AFF] dark:text-[#0A84FF]">
            {liveData.engine} · Electrostatics & H-Bond Network
          </div>
          <h2 className="text-[16px] font-semibold tracking-[-0.015em]">
            pH-Dependent Protonation & Histidine Tautomer Assignment
          </h2>
        </div>
        <button disabled={busy} onClick={handleApply} className="apple-btn-primary">
          <Zap className="w-3.5 h-3.5" />
          <span>Apply Protonation & Build H</span>
        </button>
      </div>

      {/* Interactive pH Slider + Live Titration Curve */}
      <div className="apple-card p-4 space-y-4">
        <div className="flex flex-wrap items-center justify-between gap-3">
          <div className="flex items-center gap-3 flex-1 min-w-[240px]">
            <Droplets className="w-5 h-5 text-[#007AFF]" />
            <div className="flex-1">
              <div className="flex items-center justify-between text-[12.5px] font-semibold mb-1.5">
                <span>Target Physiological / Assay pH</span>
                <span className="font-mono text-[#007AFF] dark:text-[#0A84FF] text-[15px]">pH {ph.toFixed(1)}</span>
              </div>
              <input
                type="range"
                min={0}
                max={14}
                step={0.1}
                value={ph}
                onChange={(e) => {
                  const v = parseFloat(e.target.value);
                  setPh(v);
                  refreshPreview(v, overrides, hisDefault, respectMetal, respectSS);
                }}
                className="w-full apple-slider"
              />
              <div className="flex justify-between text-[10px] font-mono text-[#86868B] mt-1">
                <span>0.0 (Acidic)</span>
                <button onClick={() => { setPh(5.5); refreshPreview(5.5, overrides, hisDefault, respectMetal, respectSS); }} className="hover:text-[#007AFF]">5.5 Endosome</button>
                <button onClick={() => { setPh(7.0); refreshPreview(7.0, overrides, hisDefault, respectMetal, respectSS); }} className="hover:text-[#007AFF]">7.0 Neutral</button>
                <button onClick={() => { setPh(7.4); refreshPreview(7.4, overrides, hisDefault, respectMetal, respectSS); }} className="hover:text-[#007AFF] font-semibold">7.4 Physio</button>
                <span>14.0 (Basic)</span>
              </div>
            </div>
          </div>

          <div className="grid grid-cols-3 gap-2 text-center">
            <div className="px-3 py-2 rounded-xl bg-black/[0.035] dark:bg-white/[0.05]">
              <div className="text-[10px] text-[#86868B]">Isoelectric pI</div>
              <div className="text-[15px] font-bold font-mono">{liveData.isoelectric_point}</div>
            </div>
            <div className="px-3 py-2 rounded-xl bg-black/[0.035] dark:bg-white/[0.05]">
              <div className="text-[10px] text-[#86868B]">Net Charge</div>
              <div className="text-[15px] font-bold font-mono text-[#007AFF] dark:text-[#0A84FF]">
                {liveData.net_charge_discrete >= 0 ? '+' : ''}
                {liveData.net_charge_discrete}
              </div>
            </div>
            <div className="px-3 py-2 rounded-xl bg-black/[0.035] dark:bg-white/[0.05]">
              <div className="text-[10px] text-[#86868B]">Shifted pKa</div>
              <div className="text-[15px] font-bold font-mono text-[#FF9500]">{liveData.shifted_count}</div>
            </div>
          </div>
        </div>

        {/* SVG Titration Curve */}
        <div className="rounded-2xl bg-black/[0.025] dark:bg-white/[0.03] p-2">
          <svg viewBox={`0 0 ${svgW} ${svgH}`} className="w-full h-36 overflow-visible">
            {/* Zero charge line */}
            <line
              x1={padL}
              y1={toY(0)}
              x2={svgW - padR}
              y2={toY(0)}
              stroke="currentColor"
              className="text-black/15 dark:text-white/15"
              strokeDasharray="3 3"
            />
            {/* pI marker */}
            <line
              x1={toX(liveData.isoelectric_point)}
              y1={padT}
              x2={toX(liveData.isoelectric_point)}
              y2={svgH - padB}
              stroke="#007AFF"
              strokeWidth="1.2"
              strokeDasharray="3 3"
            />
            <text
              x={toX(liveData.isoelectric_point) + 4}
              y={padT + 10}
              fontSize="9.5"
              fill="#007AFF"
              fontFamily="monospace"
            >
              pI={liveData.isoelectric_point}
            </text>

            {/* Titration Path */}
            <path d={pathD} fill="none" stroke="#007AFF" strokeWidth="2.4" strokeLinecap="round" />

            {/* Current pH indicator */}
            <line
              x1={toX(ph)}
              y1={padT}
              x2={toX(ph)}
              y2={svgH - padB}
              stroke="#34C759"
              strokeWidth="1.5"
            />
            <circle
              cx={toX(ph)}
              cy={toY(liveData.net_charge_continuous)}
              r="4.5"
              fill="#34C759"
              stroke="#FFFFFF"
              strokeWidth="1.5"
            />

            {/* Axis labels */}
            {[0, 3.5, 7, 10.5, 14].map((tick) => (
              <text
                key={tick}
                x={toX(tick)}
                y={svgH - 6}
                textAnchor="middle"
                fontSize="9.5"
                className="fill-[#86868B]"
                fontFamily="monospace"
              >
                {tick}
              </text>
            ))}
            <text x={padL - 6} y={padT + 4} textAnchor="end" fontSize="9.5" className="fill-[#86868B]" fontFamily="monospace">
              +{Math.round(maxQ)}
            </text>
            <text x={padL - 6} y={toY(0) + 3} textAnchor="end" fontSize="9.5" className="fill-[#86868B]" fontFamily="monospace">
              0
            </text>
            <text x={padL - 6} y={svgH - padB} textAnchor="end" fontSize="9.5" className="fill-[#86868B]" fontFamily="monospace">
              -{Math.round(maxQ)}
            </text>
          </svg>
        </div>

        {/* Protonation Policy Toggles */}
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-2 text-[12px]">
          <label className="flex items-center justify-between px-3 py-2 rounded-xl bg-black/[0.025] dark:bg-white/[0.03] cursor-pointer">
            <span>Lock Metal-Coordinating Residues (CYM / HID / HIE)</span>
            <input
              type="checkbox"
              checked={respectMetal}
              onChange={(e) => {
                setRespectMetal(e.target.checked);
                refreshPreview(ph, overrides, hisDefault, e.target.checked, respectSS);
              }}
              className="accent-[#007AFF] w-4 h-4"
            />
          </label>
          <label className="flex items-center justify-between px-3 py-2 rounded-xl bg-black/[0.025] dark:bg-white/[0.03] cursor-pointer">
            <span>Protect Disulfide Bridges (CYX)</span>
            <input
              type="checkbox"
              checked={respectSS}
              onChange={(e) => {
                setRespectSS(e.target.checked);
                refreshPreview(ph, overrides, hisDefault, respectMetal, e.target.checked);
              }}
              className="accent-[#007AFF] w-4 h-4"
            />
          </label>
          <label className="flex items-center justify-between px-3 py-2 rounded-xl bg-black/[0.025] dark:bg-white/[0.03] cursor-pointer">
            <span>Optimize H-Bond Network (Auto NQH Flips)</span>
            <input
              type="checkbox"
              checked={optimizeHbonds}
              onChange={(e) => setOptimizeHbonds(e.target.checked)}
              className="accent-[#007AFF] w-4 h-4"
            />
          </label>
          <div className="flex items-center justify-between px-3 py-2 rounded-xl bg-black/[0.025] dark:bg-white/[0.03]">
            <span>Histidine Tautomer Strategy</span>
            <select
              value={hisDefault}
              onChange={(e) => {
                const v = e.target.value as any;
                setHisDefault(v);
                refreshPreview(ph, overrides, v, respectMetal, respectSS);
              }}
              className="apple-select text-[11.5px]"
            >
              <option value="auto">Auto (PROPKA + H-bond)</option>
              <option value="HID">Force HID (δ-H)</option>
              <option value="HIE">Force HIE (ε-H)</option>
              <option value="HIP">Force HIP (+1)</option>
            </select>
          </div>
        </div>
      </div>

      {/* Titratable Residues Table with Per-Residue Override */}
      <div className="apple-card p-4 space-y-2.5">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2 text-[13.5px] font-semibold">
            <ShieldCheck className="w-4 h-4 text-[#34C759]" />
            <span>Titratable Residues ({liveData.residues.length})</span>
          </div>
          {Object.keys(overrides).length > 0 && (
            <button
              onClick={() => {
                setOverrides({});
                refreshPreview(ph, {}, hisDefault, respectMetal, respectSS);
              }}
              className="text-[11.5px] text-[#007AFF] hover:underline"
            >
              Reset {Object.keys(overrides).length} override(s)
            </button>
          )}
        </div>

        <div className="overflow-x-auto max-h-[340px] overflow-y-auto">
          <table className="w-full text-left border-collapse text-[12px]">
            <thead className="sticky top-0 bg-white/95 dark:bg-[#1C1C22]/95 text-[10.5px] uppercase tracking-wider text-[#86868B] border-b border-black/[0.07] dark:border-white/[0.08]">
              <tr>
                <th className="py-2 pr-2">Residue</th>
                <th className="py-2 px-2">pKa (Model)</th>
                <th className="py-2 px-2">ΔpKa</th>
                <th className="py-2 px-2">Buried</th>
                <th className="py-2 px-2">Assigned State</th>
                <th className="py-2 pl-2">Microenvironment</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-black/[0.045] dark:divide-white/[0.05]">
              {liveData.residues.map((r) => (
                <tr
                  key={r.key}
                  onClick={() => onSelectResidue({ chain: r.chain, resseq: r.resseq, icode: r.icode, resname: r.resname }, true)}
                  className="hover:bg-black/[0.025] dark:hover:bg-white/[0.03] cursor-pointer transition"
                >
                  <td className="py-2 pr-2 font-mono font-semibold whitespace-nowrap">
                    {r.resname} {r.chain}:{r.resseq}
                  </td>
                  <td className="py-2 px-2 font-mono whitespace-nowrap">
                    <span className="font-semibold">{r.pka.toFixed(2)}</span>{' '}
                    <span className="text-[#86868B]">({r.model_pka.toFixed(1)})</span>
                  </td>
                  <td className="py-2 px-2 font-mono whitespace-nowrap">
                    <span
                      className={
                        r.shifted_significantly
                          ? 'text-[#FF9500] font-semibold'
                          : 'text-[#86868B]'
                      }
                    >
                      {r.delta_pka >= 0 ? '+' : ''}
                      {r.delta_pka.toFixed(2)}
                    </span>
                  </td>
                  <td className="py-2 px-2 font-mono text-[#86868B]">{r.buried_pct.toFixed(0)}%</td>
                  <td className="py-2 px-2" onClick={(e) => e.stopPropagation()}>
                    {r.available_states.length > 1 ? (
                      <select
                        value={r.assigned_state}
                        onChange={(e) => {
                          const next = { ...overrides, [r.base_key]: e.target.value };
                          setOverrides(next);
                          refreshPreview(ph, next, hisDefault, respectMetal, respectSS);
                        }}
                        className={`apple-select !py-0.5 !px-2 text-[11.5px] font-mono font-semibold ${
                          r.is_override ? 'ring-1 ring-[#007AFF]' : ''
                        }`}
                      >
                        {r.available_states.map((st) => (
                          <option key={st} value={st}>
                            {st}
                          </option>
                        ))}
                      </select>
                    ) : (
                      <span className="apple-pill bg-black/[0.05] dark:bg-white/[0.08] font-mono">
                        {r.assigned_state}
                      </span>
                    )}
                  </td>
                  <td className="py-2 pl-2 text-[11px] text-[#6E6E73] dark:text-[#A1A1A6]">
                    {r.metal_coordination.length > 0 && (
                      <span className="apple-pill bg-[#AF52DE]/15 text-[#007AFF] dark:text-[#0A84FF] mr-1.5">
                        Metal
                      </span>
                    )}
                    {r.is_disulfide && (
                      <span className="apple-pill bg-[#34C759]/15 text-[#1B8A38] dark:text-[#30D158] mr-1.5">
                        S–S
                      </span>
                    )}
                    <span>{r.rationale}</span>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
};
