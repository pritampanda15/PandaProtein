import React, { useRef, useState } from 'react';
import { CheckCircle2, Download, FolderKanban, Play, Trash2, Upload, XCircle } from 'lucide-react';
import { api } from '../api/client';
import type { BatchPrepResult } from '../types/pdbprep';

interface BatchPanelProps {
  busy: boolean;
  onSetBusy: (msg: string | null) => void;
  onError: (msg: string) => void;
}

export const BatchPanel: React.FC<BatchPanelProps> = ({ busy, onSetBusy, onError }) => {
  const [files, setFiles] = useState<File[]>([]);
  const [preset, setPreset] = useState<'md_production' | 'docking_receptor' | 'crystallographic_polish'>('md_production');
  const [ph, setPh] = useState(7.4);
  const [includeMole, setIncludeMole] = useState(false);
  const [result, setResult] = useState<BatchPrepResult | null>(null);
  const fileRef = useRef<HTMLInputElement | null>(null);

  const addFiles = (list: FileList | null) => {
    if (!list) return;
    const incoming = Array.from(list).filter((f) => /\.(pdb|ent|cif|mmcif|gro)$/i.test(f.name));
    setFiles((prev) => {
      const names = new Set(prev.map((p) => p.name));
      return [...prev, ...incoming.filter((f) => !names.has(f.name))];
    });
  };

  const handleRun = async () => {
    if (files.length === 0) return;
    onSetBusy(`Batch preparing ${files.length} structure(s)…`);
    try {
      const res = await api.runBatchPrep(files, preset, ph, includeMole);
      setResult(res);
    } catch (e: any) {
      onError(e?.message || 'Batch preparation failed');
    } finally {
      onSetBusy(null);
    }
  };

  return (
    <div className="space-y-4 pb-6">
      <div className="apple-card p-4 flex flex-wrap items-center justify-between gap-3">
        <div>
          <div className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#007AFF] dark:text-[#0A84FF]">
            High-Throughput Virtual Screening &amp; MD Queue
          </div>
          <h2 className="text-[16px] font-semibold tracking-[-0.015em]">
            Batch Multi-Structure Preparation Pipeline
          </h2>
        </div>
        <button disabled={busy || files.length === 0} onClick={handleRun} className="apple-btn-primary">
          <Play className="w-3.5 h-3.5 fill-current" />
          <span>Run Batch ({files.length})</span>
        </button>
      </div>

      {/* Multi-file Dropzone + Preset Config */}
      <div className="apple-card p-4 space-y-4">
        <div
          onDragOver={(e) => e.preventDefault()}
          onDrop={(e) => {
            e.preventDefault();
            addFiles(e.dataTransfer.files);
          }}
          onClick={() => fileRef.current?.click()}
          className="cursor-pointer rounded-2xl border-2 border-dashed border-black/[0.12] dark:border-white/[0.15] hover:border-[#007AFF]/60 p-6 text-center transition bg-black/[0.02] dark:bg-white/[0.02]"
        >
          <input
            ref={fileRef}
            type="file"
            multiple
            accept=".pdb,.ent,.cif,.mmcif,.gro"
            className="hidden"
            onChange={(e) => addFiles(e.target.files)}
          />
          <Upload className="w-6 h-6 mx-auto text-[#007AFF] mb-2" />
          <div className="text-[13.5px] font-medium">
            Drop multiple .pdb, .cif, or .gro files here (up to 30 structures)
          </div>
          <div className="text-[11.5px] text-[#86868B] mt-0.5">
            Produces prepared .pdb, .gro, receptor .pdbqt, vina_config.txt &amp; unified batch_summary.csv in one ZIP
          </div>
        </div>

        {/* Settings */}
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-2.5 text-[12px]">
          <div>
            <label className="text-[10.5px] text-[#86868B] block mb-1">Preparation Preset</label>
            <select value={preset} onChange={(e) => setPreset(e.target.value as any)} className="apple-select w-full">
              <option value="md_production">Full MD Production (PDBFixer + COOT + PROPKA)</option>
              <option value="docking_receptor">Docking Receptor (Clean + Polar H + PDBQT)</option>
              <option value="crystallographic_polish">COOT Auto-Refine Only</option>
            </select>
          </div>
          <div>
            <label className="text-[10.5px] text-[#86868B] block mb-1">Target pH (PROPKA 3)</label>
            <input
              type="number"
              step={0.1}
              min={0}
              max={14}
              value={ph}
              onChange={(e) => setPh(parseFloat(e.target.value) || 7.4)}
              className="apple-input !py-1.5 font-mono"
            />
          </div>
          <div className="flex items-end">
            <label className="w-full flex items-center justify-between px-3 py-2 rounded-xl bg-black/[0.03] dark:bg-white/[0.04] cursor-pointer">
              <span>Run MOLE Tunnels per file</span>
              <input
                type="checkbox"
                checked={includeMole}
                onChange={(e) => setIncludeMole(e.target.checked)}
                className="accent-[#007AFF] w-4 h-4"
              />
            </label>
          </div>
        </div>

        {/* Queued Files List */}
        {files.length > 0 && (
          <div className="space-y-1.5">
            <div className="flex items-center justify-between text-[12px]">
              <span className="font-semibold flex items-center gap-1.5">
                <FolderKanban className="w-3.5 h-3.5 text-[#007AFF]" />
                Queued Structures ({files.length})
              </span>
              <button onClick={() => setFiles([])} className="text-[#FF3B30] hover:underline text-[11.5px]">
                Clear Queue
              </button>
            </div>
            <div className="flex flex-wrap gap-1.5">
              {files.map((f) => (
                <span
                  key={f.name}
                  className="inline-flex items-center gap-1.5 px-2.5 py-1 rounded-xl bg-black/[0.04] dark:bg-white/[0.06] text-[11.5px] font-mono"
                >
                  <span>{f.name}</span>
                  <button
                    onClick={() => setFiles((prev) => prev.filter((x) => x.name !== f.name))}
                    className="text-[#86868B] hover:text-[#FF3B30]"
                  >
                    <Trash2 className="w-3 h-3" />
                  </button>
                </span>
              ))}
            </div>
          </div>
        )}
      </div>

      {/* Batch Audit Results Table + ZIP Download */}
      {result && (
        <div className="apple-card p-4 space-y-3">
          <div className="flex flex-wrap items-center justify-between gap-2">
            <div>
              <div className="text-[13.5px] font-semibold text-[#1B8A38] dark:text-[#30D158]">
                Batch Complete — {result.succeeded}/{result.total} Structures Prepared
              </div>
              <div className="text-[11px] text-[#86868B]">
                Preset: {result.preset} · pH {result.ph} · Bundle size {(result.zip_size_bytes / 1024).toFixed(1)} KB
              </div>
            </div>
            <button
              onClick={() => api.downloadBatchZip(result.batch_id)}
              className="apple-btn-primary"
            >
              <Download className="w-4 h-4" />
              <span>Download Batch ZIP + CSV</span>
            </button>
          </div>

          <div className="overflow-x-auto">
            <table className="w-full text-left border-collapse text-[12px]">
              <thead className="text-[10.5px] uppercase tracking-wider text-[#86868B] border-b border-black/[0.07] dark:border-white/[0.08]">
                <tr>
                  <th className="py-2 pr-2">File</th>
                  <th className="py-2 px-2">Status</th>
                  <th className="py-2 px-2">Readiness</th>
                  <th className="py-2 px-2">Clashscore</th>
                  <th className="py-2 px-2">Issues</th>
                  <th className="py-2 px-2">+H</th>
                  <th className="py-2 pl-2">Charge / pI</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-black/[0.05] dark:divide-white/[0.06] font-mono">
                {result.rows.map((r) => (
                  <tr key={r.filename}>
                    <td className="py-2 pr-2 font-semibold">{r.filename}</td>
                    <td className="py-2 px-2">
                      {r.status === 'ok' ? (
                        <CheckCircle2 className="w-4 h-4 text-[#34C759]" />
                      ) : (
                        <span className="flex items-center gap-1 text-[#FF3B30]" title={r.error || ''}>
                          <XCircle className="w-4 h-4" /> Error
                        </span>
                      )}
                    </td>
                    <td className="py-2 px-2">
                      {r.status === 'ok' ? `${r.before_readiness}% → ${r.after_readiness}%` : '—'}
                    </td>
                    <td className="py-2 px-2">
                      {r.status === 'ok' ? `${r.before_clashscore} → ${r.after_clashscore}` : '—'}
                    </td>
                    <td className="py-2 px-2">
                      {r.status === 'ok' ? `${r.before_issues} → ${r.after_issues}` : '—'}
                    </td>
                    <td className="py-2 px-2">{r.status === 'ok' ? `+${r.hydrogens_added}` : '—'}</td>
                    <td className="py-2 pl-2">
                      {r.status === 'ok' ? `${(r.net_charge || 0) >= 0 ? '+' : ''}${r.net_charge} / pI ${r.pI}` : r.error}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </div>
  );
};
