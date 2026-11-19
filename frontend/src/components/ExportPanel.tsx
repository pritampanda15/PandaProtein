import React, { useEffect, useState } from 'react';
import {
  AlertCircle,
  CheckCircle2,
  Copy,
  Download,
  FileArchive,
  FileCode2,
  FileText,
  FlaskConical,
  Terminal,
} from 'lucide-react';
import { api } from '../api/client';
import type { GromacsTestResult, SolvationEstimate, StructureSessionState } from '../types/pdbprep';

interface ExportPanelProps {
  session: StructureSessionState;
  busy: boolean;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
}

export const ExportPanel: React.FC<ExportPanelProps> = ({ session, busy, onSetBusy, onError }) => {
  const [forceField, setForceField] = useState('charmm36-jul2022');
  const [waterModel, setWaterModel] = useState('tip3p');
  const [boxType, setBoxType] = useState('dodecahedron');
  const [boxDist, setBoxDist] = useState(1.0);
  const [saltConc, setSaltConc] = useState(0.15);
  const [padding, setPadding] = useState(6.0);
  const [selectedLig, setSelectedLig] = useState(session.ligands_info.docking_boxes[0]?.ligand_id || '');
  const [vinaConfig, setVinaConfig] = useState('');
  const [copied, setCopied] = useState(false);
  const [gmxResult, setGmxResult] = useState<GromacsTestResult | null>(null);
  const [solvEst, setSolvEst] = useState<SolvationEstimate | null>(session.solvation_estimate || null);
  const [showGmxLog, setShowGmxLog] = useState(false);

  useEffect(() => {
    setSolvEst(session.solvation_estimate || null);
    setGmxResult(null);
  }, [session.session_id, session.solvation_estimate]);

  const handleTestGromacs = async () => {
    onSetBusy(`Running live gmx pdb2gmx (${forceField} / ${waterModel})…`);
    try {
      const resp = await api.testGromacs(session.session_id, {
        force_field: forceField,
        water_model: waterModel,
        box_type: boxType,
        box_distance_nm: boxDist,
        salt_conc_m: saltConc,
      });
      setGmxResult(resp.gromacs);
      setSolvEst(resp.solvation);
    } catch (e: any) {
      onError(e?.message || 'GROMACS pre-flight test failed');
    } finally {
      onSetBusy(null);
    }
  };

  useEffect(() => {
    api
      .getVinaConfig(session.session_id, selectedLig, padding)
      .then((d) => setVinaConfig(d.config))
      .catch(() => {});
  }, [session.session_id, selectedLig, padding, session.pdb_text]);

  const doExport = async (target: string, fallbackName: string, label: string) => {
    onSetBusy(`Building ${label}…`);
    try {
      await api.exportFile(
        session.session_id,
        {
          target,
          force_field: forceField,
          water_model: waterModel,
          box_type: boxType,
          box_distance_nm: boxDist,
          salt_conc_m: saltConc,
          ph: session.propka.ph,
          docking_padding: padding,
          selected_ligand_id: selectedLig || null,
        },
        fallbackName
      );
    } catch (e: any) {
      onError(e?.message || 'Export failed');
    } finally {
      onSetBusy(null);
    }
  };

  return (
    <div className="space-y-4 pb-6">
      <div className="apple-card p-4">
        <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#007AFF] dark:text-[#0A84FF]">
          Production Export Center
        </div>
        <h2 className="text-[16px] font-semibold tracking-[-0.015em]">
          GROMACS MD Simulation Bundle, AutoDock Vina / GNINA PDBQT & Reports
        </h2>
      </div>

      {/* Primary Hero Export: Complete Production Bundle */}
      <div className="apple-card p-5 bg-gradient-to-br from-[#007AFF]/10 via-[#5856D6]/5 to-transparent border-[#007AFF]/25 space-y-4">
        <div className="flex flex-wrap items-center justify-between gap-3">
          <div>
            <div className="flex items-center gap-2 text-[15px] font-semibold">
              <FileArchive className="w-5 h-5 text-[#007AFF]" />
              <span>Complete GROMACS MD + Docking Production Bundle (.ZIP)</span>
            </div>
            <p className="text-[12px] text-[#6E6E73] dark:text-[#A1A1A6] mt-0.5">
              Includes prepared PDB & mmCIF, live <code className="font-mono">gmx pdb2gmx</code> topology (if GROMACS available), em/nvt/npt/md MDP files, executable runner script, receptor PDBQT, vina_config.txt, PQR, and full Markdown/JSON audit reports.
            </p>
          </div>
          <div className="flex items-center gap-2">
            <button
              disabled={busy}
              onClick={handleTestGromacs}
              className="apple-btn-secondary !px-3.5 !py-2.5"
              title="Execute live `gmx pdb2gmx` topology dry-run on your prepared protein"
            >
              <FlaskConical className="w-4 h-4 text-[#34C759]" />
              <span>Test in GROMACS</span>
            </button>
            <button
              disabled={busy}
              onClick={() => doExport('gromacs_bundle', 'pdbprep_bundle.zip', 'GROMACS + Docking Bundle')}
              className="apple-btn-primary !px-5 !py-2.5"
            >
              <Download className="w-4 h-4" />
              <span>Download Full Bundle (.zip)</span>
            </button>
          </div>
        </div>

        {/* GROMACS Parameters */}
        <div className="grid grid-cols-2 sm:grid-cols-5 gap-2.5 pt-2 border-t border-black/[0.06] dark:border-white/[0.08] text-[12px]">
          <div>
            <label className="text-[10.5px] text-[#86868B] block mb-1">Force Field</label>
            <select value={forceField} onChange={(e) => setForceField(e.target.value)} className="apple-select w-full">
              <option value="charmm36-jul2022">CHARMM36m (jul2022 · Recommended)</option>
              <option value="amber99sb-ildn">AMBER99SB-ILDN</option>
              <option value="charmm27">CHARMM27 (Legacy)</option>
              <option value="oplsaa">OPLS-AA/L</option>
              <option value="gromos54a7">GROMOS96 54a7</option>
              <option value="amber03">AMBER03</option>
            </select>
          </div>
          <div>
            <label className="text-[10.5px] text-[#86868B] block mb-1">Water Model</label>
            <select value={waterModel} onChange={(e) => setWaterModel(e.target.value)} className="apple-select w-full">
              <option value="tip3p">TIP3P</option>
              <option value="tip4p">TIP4P</option>
              <option value="spce">SPC/E</option>
              <option value="spc">SPC</option>
            </select>
          </div>
          <div>
            <label className="text-[10.5px] text-[#86868B] block mb-1">Box Geometry</label>
            <select value={boxType} onChange={(e) => setBoxType(e.target.value)} className="apple-select w-full">
              <option value="dodecahedron">Dodecahedron</option>
              <option value="cubic">Cubic</option>
              <option value="octahedron">Octahedron</option>
            </select>
          </div>
          <div>
            <label className="text-[10.5px] text-[#86868B] block mb-1">Margin (nm)</label>
            <input
              type="number"
              step={0.1}
              min={0.5}
              max={2.5}
              value={boxDist}
              onChange={(e) => setBoxDist(Number(e.target.value))}
              className="apple-input !py-1.5 font-mono"
            />
          </div>
          <div>
            <label className="text-[10.5px] text-[#86868B] block mb-1">NaCl (M)</label>
            <input
              type="number"
              step={0.05}
              min={0}
              max={1.0}
              value={saltConc}
              onChange={(e) => setSaltConc(Number(e.target.value))}
              className="apple-input !py-1.5 font-mono"
            />
          </div>
        </div>

        {/* Solvation & Ion System Estimator */}
        {solvEst && solvEst.box_volume_nm3 > 0 && (
          <div className="grid grid-cols-2 sm:grid-cols-5 gap-2 pt-2 border-t border-black/[0.06] dark:border-white/[0.08] text-center">
            <div className="p-2 rounded-xl bg-white/70 dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#86868B]">Box Volume</div>
              <div className="text-[12.5px] font-bold font-mono">{solvEst.box_volume_nm3} nm³</div>
            </div>
            <div className="p-2 rounded-xl bg-white/70 dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#86868B]">Est. Waters</div>
              <div className="text-[12.5px] font-bold font-mono text-[#007AFF] dark:text-[#0A84FF]">
                {solvEst.estimated_waters.toLocaleString()}
              </div>
            </div>
            <div className="p-2 rounded-xl bg-white/70 dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#86868B]">Na⁺ Ions</div>
              <div className="text-[12.5px] font-bold font-mono">{solvEst.na_ions}</div>
            </div>
            <div className="p-2 rounded-xl bg-white/70 dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#86868B]">Cl⁻ Ions</div>
              <div className="text-[12.5px] font-bold font-mono">{solvEst.cl_ions}</div>
            </div>
            <div className="p-2 rounded-xl bg-white/70 dark:bg-white/[0.04]">
              <div className="text-[10px] text-[#86868B]">Total MD Atoms</div>
              <div className="text-[12.5px] font-bold font-mono">{solvEst.total_atoms_est.toLocaleString()}</div>
            </div>
          </div>
        )}

        {/* Live GROMACS pdb2gmx Test Result */}
        {gmxResult && (
          <div
            className={`p-3 rounded-2xl border space-y-2 text-[12px] ${
              gmxResult.passed
                ? 'bg-[#34C759]/10 border-[#34C759]/30'
                : 'bg-[#FF3B30]/10 border-[#FF3B30]/30'
            }`}
          >
            <div className="flex flex-wrap items-center justify-between gap-2">
              <div className="flex items-center gap-2 font-semibold">
                {gmxResult.passed ? (
                  <CheckCircle2 className="w-4 h-4 text-[#34C759]" />
                ) : (
                  <AlertCircle className="w-4 h-4 text-[#FF3B30]" />
                )}
                <span className={gmxResult.passed ? 'text-[#1B8A38] dark:text-[#30D158]' : 'text-[#C92519] dark:text-[#FF453A]'}>
                  {gmxResult.message}
                </span>
              </div>
              {gmxResult.log && (
                <button
                  onClick={() => setShowGmxLog((v) => !v)}
                  className="apple-btn-secondary !py-0.5 !px-2 text-[11px]"
                >
                  <Terminal className="w-3 h-3" />
                  {showGmxLog ? 'Hide Log' : 'View gmx Log'}
                </button>
              )}
            </div>
            {showGmxLog && gmxResult.log && (
              <pre className="p-2.5 rounded-xl bg-black/80 text-[#E5E5EA] text-[10.5px] font-mono max-h-48 overflow-auto">
                {gmxResult.log}
              </pre>
            )}
          </div>
        )}
      </div>

      {/* Individual Format Downloads */}
      <div className="grid grid-cols-2 sm:grid-cols-3 gap-2.5">
        {[
          { target: 'pdb', file: 'prepared.pdb', title: 'Prepared PDB', sub: 'Cleaned coordinates + H', icon: <FileCode2 className="w-4 h-4 text-[#007AFF]" /> },
          { target: 'gro', file: 'prepared.gro', title: 'GROMACS .gro', sub: 'nm coordinates + box', icon: <FileCode2 className="w-4 h-4 text-[#34C759]" /> },
          { target: 'mmcif', file: 'prepared.cif', title: 'Prepared mmCIF', sub: 'PDBx / Crystallographic', icon: <FileCode2 className="w-4 h-4 text-[#5856D6]" /> },
          { target: 'pdbqt', file: 'receptor.pdbqt', title: 'Receptor PDBQT', sub: 'AutoDock Vina / GNINA', icon: <FileCode2 className="w-4 h-4 text-[#AF52DE]" /> },
          { target: 'pqr', file: 'prepared.pqr', title: 'APBS PQR', sub: 'Charges + Bondi Radii', icon: <FileCode2 className="w-4 h-4 text-[#5AC8FA]" /> },
          { target: 'report_md', file: 'report.md', title: 'QA Report (.md)', sub: 'Publication Summary', icon: <FileText className="w-4 h-4 text-[#FF9500]" /> },
          { target: 'report_json', file: 'audit.json', title: 'Full Audit (.json)', sub: 'Machine-readable metrics', icon: <FileText className="w-4 h-4 text-[#86868B]" /> },
        ].map((item) => (
          <button
            key={item.target}
            disabled={busy}
            onClick={() => doExport(item.target, item.file, item.title)}
            className="apple-card p-3.5 text-left hover:border-[#007AFF]/45 transition flex items-start gap-2.5"
          >
            <div className="p-2 rounded-xl bg-black/[0.04] dark:bg-white/[0.06]">{item.icon}</div>
            <div className="min-w-0">
              <div className="text-[13px] font-semibold truncate">{item.title}</div>
              <div className="text-[11px] text-[#86868B] truncate">{item.sub}</div>
            </div>
          </button>
        ))}
      </div>

      {/* AutoDock Vina / GNINA Docking Box Preview */}
      <div className="apple-card p-4 space-y-3">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div>
            <h3 className="text-[13.5px] font-semibold">AutoDock Vina / GNINA Grid Box Configuration</h3>
            <p className="text-[11.5px] text-[#86868B]">
              Enable the "Grid" overlay in the 3D Viewport to inspect this docking box live on the receptor
            </p>
          </div>
          <div className="flex items-center gap-2">
            {session.ligands_info.docking_boxes.length > 0 && (
              <select
                value={selectedLig}
                onChange={(e) => setSelectedLig(e.target.value)}
                className="apple-select text-[12px]"
              >
                {session.ligands_info.docking_boxes.map((b) => (
                  <option key={b.ligand_id} value={b.ligand_id}>
                    {b.label}
                  </option>
                ))}
              </select>
            )}
            <button
              onClick={() => {
                navigator.clipboard.writeText(vinaConfig);
                setCopied(true);
                setTimeout(() => setCopied(false), 1800);
              }}
              className="apple-btn-secondary text-[12px]"
            >
              {copied ? <CheckCircle2 className="w-3.5 h-3.5 text-[#34C759]" /> : <Copy className="w-3.5 h-3.5" />}
              <span>{copied ? 'Copied' : 'Copy Config'}</span>
            </button>
          </div>
        </div>

        <pre className="p-3 rounded-xl bg-black/[0.035] dark:bg-black/35 text-[11.5px] font-mono overflow-x-auto leading-relaxed">
          {vinaConfig || 'Computing docking box…'}
        </pre>
      </div>
    </div>
  );
};
