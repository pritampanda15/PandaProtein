import type {
  BatchPrepResult,
  GromacsTestResult,
  HoleAnalysisData,
  InteractionsMembraneData,
  MoleAnalysisData,
  MoleMutationScanResult,
  PropkaData,
  RotamerExplorerData,
  SolvationEstimate,
  StructuralDiffData,
  StructureSessionState,
} from '../types/pdbprep';

async function request<T>(url: string, options?: RequestInit): Promise<T> {
  const res = await fetch(url, {
    headers: options?.body instanceof FormData ? undefined : { 'Content-Type': 'application/json' },
    ...options,
  });
  if (!res.ok) {
    let msg = `Request failed (${res.status})`;
    try {
      const err = await res.json();
      if (err?.detail) msg = err.detail;
    } catch {
      // ignore
    }
    throw new Error(msg);
  }
  return res.json() as Promise<T>;
}

export const api = {
  health: () => request<{ status: string; version: string; engines: Record<string, unknown> }>('/api/health'),

  loadSample: (sample_id = 'catalytic_showcase') =>
    request<StructureSessionState>('/api/import/sample', {
      method: 'POST',
      body: JSON.stringify({ sample_id }),
    }),

  fetchRCSB: (pdb_id: string, format: 'pdb' | 'cif' = 'pdb') =>
    request<StructureSessionState>('/api/import/rcsb', {
      method: 'POST',
      body: JSON.stringify({ pdb_id, format }),
    }),

  fetchAlphaFold: (uniprot_id: string) =>
    request<StructureSessionState>('/api/import/alphafold', {
      method: 'POST',
      body: JSON.stringify({ uniprot_id }),
    }),

  uploadFile: (file: File) => {
    const fd = new FormData();
    fd.append('file', file);
    return request<StructureSessionState>('/api/import/upload', {
      method: 'POST',
      body: fd,
    });
  },

  undo: (sid: string) =>
    request<{ undone: string | null; state: StructureSessionState }>(`/api/sessions/${sid}/undo`, {
      method: 'POST',
    }),

  redo: (sid: string) =>
    request<{ redone: string | null; state: StructureSessionState }>(`/api/sessions/${sid}/redo`, {
      method: 'POST',
    }),

  deepPrepPreset: (sid: string, preset: 'md_production' | 'docking_receptor' | 'crystallographic_polish', ph = 7.4) =>
    request<{ actions: string[]; state: StructureSessionState }>(`/api/sessions/${sid}/deep-prep`, {
      method: 'POST',
      body: JSON.stringify({ preset, ph }),
    }),

  runPdbFixer: (sid: string, options: Record<string, unknown>) =>
    request<{ result: { actions: string[] }; state: StructureSessionState }>(`/api/sessions/${sid}/pdbfixer/run`, {
      method: 'POST',
      body: JSON.stringify(options),
    }),

  previewPropka: (sid: string, options: Record<string, unknown>) =>
    request<{ propka: PropkaData }>(`/api/sessions/${sid}/propka/preview`, {
      method: 'POST',
      body: JSON.stringify(options),
    }),

  applyPropka: (sid: string, options: Record<string, unknown>) =>
    request<{ result: { actions: string[]; hydrogens_added: number }; state: StructureSessionState }>(
      `/api/sessions/${sid}/propka/apply`,
      {
        method: 'POST',
        body: JSON.stringify(options),
      }
    ),

  cootListRotamers: (sid: string, chain: string, resseq: number, icode = '') =>
    request<RotamerExplorerData>(
      `/api/sessions/${sid}/coot/rotamers?chain=${encodeURIComponent(chain)}&resseq=${resseq}&icode=${encodeURIComponent(icode)}`
    ),

  cootApplyRotamer: (
    sid: string,
    payload: { chain: string; resseq: number; icode?: string; rotamer_index?: number | null; chi_angles?: number[] | null }
  ) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/coot/rotamer`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  cootMutate: (sid: string, payload: { chain: string; resseq: number; icode?: string; target_resname: string }) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/coot/mutate`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  cootPepflip: (sid: string, payload: { chain: string; resseq: number; icode?: string }) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/coot/pepflip`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  cootCisTrans: (sid: string, payload: { chain: string; resseq: number; icode?: string }) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/coot/cis-trans`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  cootNqhFlip: (sid: string, payload: { chain: string; resseq: number; icode?: string }) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/coot/nqh-flip`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  cootAltloc: (sid: string, payload: { chain: string; resseq: number; icode?: string; altloc: string }) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/coot/altloc`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  cootFixChirality: (sid: string, payload: { chain: string; resseq: number; icode?: string }) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/coot/fix-chirality`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  cootRegularize: (sid: string, payload: { chain: string; start_resseq: number; end_resseq: number; iterations?: number }) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/coot/regularize`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  cootDeleteResidue: (sid: string, payload: { chain: string; resseq: number; icode?: string }) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/coot/delete-residue`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  cootAutoFix: (sid: string, payload: Record<string, boolean>) =>
    request<{
      result: {
        actions: string[];
        before_score: number;
        after_score: number;
        before_issues: number;
        after_issues: number;
        before_clashscore: number;
        after_clashscore: number;
        before_rotamers: number;
        after_rotamers: number;
        before_rama: number;
        after_rama: number;
      };
      state: StructureSessionState;
    }>(
      `/api/sessions/${sid}/coot/auto-fix`,
      {
        method: 'POST',
        body: JSON.stringify(payload),
      }
    ),

  renameFile: (sid: string, payload: { filename: string; title?: string }) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/rename/file`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  renameChains: (sid: string, mapping: Record<string, string>) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/rename/chains`, {
      method: 'POST',
      body: JSON.stringify({ mapping }),
    }),

  standardizeChains: (sid: string) =>
    request<{ result: { message: string }; state: StructureSessionState }>(
      `/api/sessions/${sid}/rename/chains/standardize`,
      { method: 'POST' }
    ),

  reorderChains: (sid: string, chain_order: string[]) =>
    request<{ result: { message: string }; state: StructureSessionState }>(
      `/api/sessions/${sid}/rename/chains/reorder`,
      { method: 'POST', body: JSON.stringify({ chain_order }) }
    ),

  mergeChains: (sid: string, payload: { source_chain: string; target_chain: string; renumber_merged: boolean }) =>
    request<{ result: { message: string }; state: StructureSessionState }>(
      `/api/sessions/${sid}/rename/chains/merge`,
      { method: 'POST', body: JSON.stringify(payload) }
    ),

  deleteChain: (sid: string, chain_id: string) =>
    request<{ result: { message: string }; state: StructureSessionState }>(
      `/api/sessions/${sid}/rename/chains/delete`,
      { method: 'POST', body: JSON.stringify({ chain_id }) }
    ),

  deleteChains: (sid: string, chain_ids: string[]) =>
    request<{ result: { message: string }; state: StructureSessionState }>(
      `/api/sessions/${sid}/rename/chains/delete`,
      { method: 'POST', body: JSON.stringify({ chain_ids }) }
    ),

  renumberChain: (
    sid: string,
    payload: {
      chain_id: string;
      mode: 'start' | 'offset';
      start_num: number;
      offset: number;
      include_heterogens: boolean;
      clear_icodes: boolean;
    }
  ) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/rename/renumber`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  renameResidueCode: (
    sid: string,
    payload: {
      old_resname: string;
      new_resname: string;
      chain_id?: string | null;
      resseq?: number | null;
      move_to_chain?: string | null;
    }
  ) =>
    request<{ result: { message: string }; state: StructureSessionState }>(
      `/api/sessions/${sid}/rename/residue-code`,
      { method: 'POST', body: JSON.stringify(payload) }
    ),

  canonicalizeResidueNames: (sid: string) =>
    request<{ result: { message: string }; state: StructureSessionState }>(
      `/api/sessions/${sid}/rename/canonicalize`,
      { method: 'POST' }
    ),

  normalizeFields: (
    sid: string,
    payload: { reset_occupancy: boolean; reset_bfactor?: number | null; strip_waters_to_chain?: string | null }
  ) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/rename/normalize`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  generateConformers: (sid: string, ligand_id: string, num_conformers = 5) =>
    request<{
      ligand_id: string;
      conformers: Array<{ rank: number; conformer_id: number; energy_kcal_mol: number; delta_e: number }>;
      sdf: string;
    }>(`/api/sessions/${sid}/ligands/conformers`, {
      method: 'POST',
      body: JSON.stringify({ ligand_id, num_conformers }),
    }),

  getVinaConfig: (sid: string, ligand_id = '', padding = 6.0) =>
    request<{ config: string }>(
      `/api/sessions/${sid}/vina-config?ligand_id=${encodeURIComponent(ligand_id)}&padding=${padding}`
    ),

  testGromacs: (
    sid: string,
    payload: {
      force_field: string;
      water_model: string;
      box_type: string;
      box_distance_nm: number;
      salt_conc_m: number;
    }
  ) =>
    request<{ gromacs: GromacsTestResult; solvation: SolvationEstimate }>(`/api/sessions/${sid}/gromacs/test`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  analyzeMole: (sid: string, payload: Record<string, unknown>) =>
    request<MoleAnalysisData>(`/api/sessions/${sid}/mole/analyze`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  exportMole: async (sid: string, payload: Record<string, unknown>, fallbackFilename: string) => {
    const res = await fetch(`/api/sessions/${sid}/mole/export`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) throw new Error(`MOLE export failed (${res.status})`);
    const blob = await res.blob();
    const cd = res.headers.get('Content-Disposition') || '';
    const match = cd.match(/filename="([^"]+)"/);
    const filename = match ? match[1] : fallbackFilename;
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  },

  scanMoleMutations: (
    sid: string,
    payload: {
      bottleneck_residues: string[];
      bottleneck_point: [number, number, number];
      baseline_bottleneck_A: number;
      target_aa: string;
      include_hetatms?: boolean;
    }
  ) =>
    request<MoleMutationScanResult>(`/api/sessions/${sid}/mole/scan-mutations`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  analyzeHole: (
    sid: string,
    payload: {
      axis_mode?: 'auto' | 'z' | 'y' | 'x';
      step?: number;
      end_radius?: number;
      custom_point?: number[] | null;
      include_hetatms?: boolean;
    }
  ) =>
    request<HoleAnalysisData>(`/api/sessions/${sid}/hole/analyze`, {
      method: 'POST',
      body: JSON.stringify(payload),
    }),

  exportHole: async (
    sid: string,
    payload: Record<string, unknown>,
    fallbackFilename: string
  ) => {
    const res = await fetch(`/api/sessions/${sid}/hole/export`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) throw new Error(`HOLE export failed (${res.status})`);
    const blob = await res.blob();
    const cd = res.headers.get('Content-Disposition') || '';
    const match = cd.match(/filename="([^"]+)"/);
    const filename = match ? match[1] : fallbackFilename;
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  },

  exportLigand: async (sid: string, ligand_id: string, format: 'pdbqt' | 'mol2' | 'sdf' | 'pdb', ph = 7.4) => {
    const res = await fetch(`/api/sessions/${sid}/ligands/export`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ ligand_id, format, ph }),
    });
    if (!res.ok) {
      let msg = `Ligand export failed (${res.status})`;
      try {
        const j = await res.json();
        if (j?.detail) msg = j.detail;
      } catch {
        // ignore
      }
      throw new Error(msg);
    }
    const blob = await res.blob();
    const cd = res.headers.get('Content-Disposition') || '';
    const match = cd.match(/filename="([^"]+)"/);
    const filename = match ? match[1] : `ligand.${format}`;
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  },

  getStructuralDiff: (sid: string) => request<StructuralDiffData>(`/api/sessions/${sid}/analysis/diff`),

  alignReferenceFile: (sid: string, file: File) => {
    const fd = new FormData();
    fd.append('file', file);
    return request<StructuralDiffData>(`/api/sessions/${sid}/analysis/align-upload`, {
      method: 'POST',
      body: fd,
    });
  },

  getInteractionsMembrane: (sid: string) =>
    request<InteractionsMembraneData>(`/api/sessions/${sid}/analysis/interactions`),

  orientAlongZ: (sid: string, axis?: [number, number, number]) =>
    request<{ result: { message: string }; state: StructureSessionState }>(`/api/sessions/${sid}/analysis/orient-z`, {
      method: 'POST',
      body: JSON.stringify(axis ? { axis } : {}),
    }),

  runBatchPrep: (files: File[], preset: string, ph: number, includeMole: boolean) => {
    const fd = new FormData();
    files.forEach((f) => fd.append('files', f));
    fd.append('preset', preset);
    fd.append('ph', String(ph));
    fd.append('include_mole', String(includeMole));
    return request<BatchPrepResult>('/api/batch/run', { method: 'POST', body: fd });
  },

  downloadBatchZip: async (batchId: string) => {
    const res = await fetch(`/api/batch/${batchId}/download`);
    if (!res.ok) throw new Error(`Batch download failed (${res.status})`);
    const blob = await res.blob();
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `pandaprotein_batch_${batchId}.zip`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  },

  exportFile: async (sid: string, payload: Record<string, unknown>, fallbackFilename: string) => {
    const res = await fetch(`/api/sessions/${sid}/export`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      throw new Error(`Export failed (${res.status})`);
    }
    const blob = await res.blob();
    const cd = res.headers.get('Content-Disposition') || '';
    const match = cd.match(/filename="([^"]+)"/);
    const filename = match ? match[1] : fallbackFilename;
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  },
};
