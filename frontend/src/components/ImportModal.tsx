import React, { useRef, useState } from 'react';
import { CloudDownload, Dna, Sparkles, Upload, X } from 'lucide-react';

interface ImportModalProps {
  open: boolean;
  busy: boolean;
  onClose: () => void;
  onUpload: (file: File) => Promise<void>;
  onFetchRCSB: (pdbId: string) => Promise<void>;
  onFetchAlphaFold: (uniprotId: string) => Promise<void>;
  onLoadSample: () => Promise<void>;
}

const FAMOUS_PDBS = [
  { id: '6LU7', desc: 'SARS-CoV-2 Mpro + N3' },
  { id: '3HTB', desc: 'T4 Lysozyme + Ligand' },
  { id: '1CRN', desc: 'Crambin (0.54 Å)' },
  { id: '4HHB', desc: 'Deoxyhemoglobin' },
];

export const ImportModal: React.FC<ImportModalProps> = ({
  open,
  busy,
  onClose,
  onUpload,
  onFetchRCSB,
  onFetchAlphaFold,
  onLoadSample,
}) => {
  const [pdbId, setPdbId] = useState('6LU7');
  const [uniprotId, setUniprotId] = useState('P00533');
  const [dragOver, setDragOver] = useState(false);
  const fileRef = useRef<HTMLInputElement | null>(null);

  if (!open) return null;

  const handleFiles = async (files: FileList | null) => {
    if (!files || files.length === 0) return;
    await onUpload(files[0]);
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/35 backdrop-blur-md p-4">
      <div className="w-full max-w-xl apple-card p-6 shadow-apple-lg space-y-5 animate-in fade-in zoom-in-95 duration-150">
        <div className="flex items-center justify-between">
          <div>
            <h2 className="text-[18px] font-semibold tracking-[-0.02em]">Load Protein Structure</h2>
            <p className="text-[12.5px] text-[#86868B] mt-0.5">
              Import local PDB / mmCIF / GROMACS .gro coordinates or fetch directly from RCSB PDB & AlphaFold DB
            </p>
          </div>
          <button
            onClick={onClose}
            className="p-1.5 rounded-full bg-black/[0.05] hover:bg-black/[0.1] dark:bg-white/[0.08] dark:hover:bg-white/[0.14]"
          >
            <X className="w-4 h-4" />
          </button>
        </div>

        {/* Drag & Drop Upload Zone */}
        <div
          onDragOver={(e) => {
            e.preventDefault();
            setDragOver(true);
          }}
          onDragLeave={() => setDragOver(false)}
          onDrop={(e) => {
            e.preventDefault();
            setDragOver(false);
            handleFiles(e.dataTransfer.files);
          }}
          onClick={() => fileRef.current?.click()}
          className={`cursor-pointer rounded-2xl border-2 border-dashed p-6 text-center transition-all ${
            dragOver
              ? 'border-[#007AFF] bg-[#007AFF]/8'
              : 'border-black/[0.12] dark:border-white/[0.14] hover:border-[#007AFF]/60 bg-black/[0.02] dark:bg-white/[0.02]'
          }`}
        >
          <input
            ref={fileRef}
            type="file"
            accept=".pdb,.ent,.cif,.mmcif,.gro"
            className="hidden"
            onChange={(e) => handleFiles(e.target.files)}
          />
          <div className="w-10 h-10 mx-auto rounded-2xl bg-[#007AFF]/12 text-[#007AFF] flex items-center justify-center mb-2.5">
            <Upload className="w-5 h-5" />
          </div>
          <div className="text-[13.5px] font-medium">Drop a .pdb, .cif, or .gro file here, or click to browse</div>
          <div className="text-[11.5px] text-[#86868B] mt-1">
            Supports X-ray, Cryo-EM, NMR ensembles, GROMACS .gro frames, AlphaFold, and Boltz/Chai complexes
          </div>
        </div>

        {/* RCSB PDB & AlphaFold DB Fetchers */}
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
          <div className="p-3.5 rounded-2xl bg-black/[0.03] dark:bg-white/[0.04] border border-black/[0.05] dark:border-white/[0.07] space-y-2.5">
            <div className="flex items-center gap-1.5 text-[12px] font-semibold">
              <CloudDownload className="w-3.5 h-3.5 text-[#007AFF]" />
              <span>Fetch RCSB PDB ID</span>
            </div>
            <div className="flex gap-1.5">
              <input
                value={pdbId}
                onChange={(e) => setPdbId(e.target.value.toUpperCase())}
                placeholder="e.g. 6LU7"
                className="apple-input font-mono uppercase"
              />
              <button
                disabled={busy || pdbId.trim().length < 4}
                onClick={() => onFetchRCSB(pdbId.trim())}
                className="apple-btn-primary shrink-0"
              >
                Fetch
              </button>
            </div>
            <div className="flex flex-wrap gap-1">
              {FAMOUS_PDBS.map((p) => (
                <button
                  key={p.id}
                  disabled={busy}
                  onClick={() => onFetchRCSB(p.id)}
                  title={p.desc}
                  className="px-2 py-0.5 rounded-lg text-[11px] font-mono bg-black/[0.05] hover:bg-[#007AFF]/15 hover:text-[#007AFF] dark:bg-white/[0.07] transition"
                >
                  {p.id}
                </button>
              ))}
            </div>
          </div>

          <div className="p-3.5 rounded-2xl bg-black/[0.03] dark:bg-white/[0.04] border border-black/[0.05] dark:border-white/[0.07] space-y-2.5">
            <div className="flex items-center gap-1.5 text-[12px] font-semibold">
              <Dna className="w-3.5 h-3.5 text-[#AF52DE]" />
              <span>Fetch AlphaFold DB</span>
            </div>
            <div className="flex gap-1.5">
              <input
                value={uniprotId}
                onChange={(e) => setUniprotId(e.target.value.toUpperCase())}
                placeholder="UniProt e.g. P00533"
                className="apple-input font-mono uppercase"
              />
              <button
                disabled={busy || uniprotId.trim().length < 4}
                onClick={() => onFetchAlphaFold(uniprotId.trim())}
                className="apple-btn-primary shrink-0 !bg-[#AF52DE] hover:!bg-[#933CC2]"
              >
                Fetch
              </button>
            </div>
            <div className="text-[11px] text-[#86868B]">
              Downloads latest AlphaFold v4 model with pLDDT confidence in B-factor column.
            </div>
          </div>
        </div>

        {/* Built-in Showcase */}
        <div className="flex items-center justify-between p-3.5 rounded-2xl bg-gradient-to-r from-[#007AFF]/10 via-[#5856D6]/10 to-[#AF52DE]/10 border border-[#007AFF]/20">
          <div>
            <div className="text-[13px] font-semibold flex items-center gap-1.5">
              <Sparkles className="w-3.5 h-3.5 text-[#007AFF]" />
              <span>Metalloenzyme Diagnostic Showcase (Offline Ready)</span>
            </div>
            <div className="text-[11.5px] text-[#6E6E73] dark:text-[#A1A1A6] mt-0.5">
              Includes Zn²⁺ site, Benzamidine ligand, S–S bridge, MSE, missing atoms/loop, AltLoc, cis-peptide & rotamer outlier
            </div>
          </div>
          <button disabled={busy} onClick={onLoadSample} className="apple-btn-primary shrink-0 ml-3">
            Load Demo
          </button>
        </div>
      </div>
    </div>
  );
};
