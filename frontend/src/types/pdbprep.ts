export type NavTab =
  | 'overview'
  | 'coot'
  | 'chains'
  | 'pdbfixer'
  | 'propka'
  | 'validation'
  | 'clashes'
  | 'ligands'
  | 'tunnels'
  | 'docking'
  | 'analysis'
  | 'batch'
  | 'export';

export interface SelectedResidue {
  chain: string;
  resseq: number;
  icode: string;
  resname?: string;
  key?: string;
}

export interface ReadinessCheck {
  label: string;
  status: 'pass' | 'warn' | 'fail';
  detail: string;
  penalty: number;
}

export interface ReadinessBlock {
  score: number;
  ready: boolean;
  checks: ReadinessCheck[];
}

export interface StructureMetadata {
  title: string;
  source: string;
  method: string;
  resolution: number | null;
  spacegroup: string | null;
  unit_cell: {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
  } | null;
  model_count: number;
  chain_count: number;
  chains: Array<{
    id: string;
    protein_residues: number;
    waters: number;
    heterogens: string[];
    atom_count: number;
    seq_range: [number, number] | null;
  }>;
  counts: {
    atoms: number;
    hydrogens: number;
    protein_residues: number;
    waters: number;
    ligands: number;
    metals_ions: number;
    buffers: number;
  };
  nonstandard_residues: Array<{
    chain: string;
    resseq: number;
    icode: string;
    resname: string;
    suggested: string;
  }>;
  altloc_residues: Array<{
    chain: string;
    resseq: number;
    icode: string;
    resname: string;
    altlocs: string[];
  }>;
}

export interface PdbFixerDiagnostics {
  nonstandard_residues: Array<{
    chain: string;
    resseq: number;
    icode: string;
    resname: string;
    replacement: string;
  }>;
  missing_atoms: Array<{
    chain: string;
    resseq: number;
    icode: string;
    resname: string;
    missing: string[];
  }>;
  missing_loops: Array<{
    chain: string;
    start_resseq: number;
    end_resseq: number;
    missing_count: number;
    gap_distance_A: number;
    sequence: string[];
    type: string;
  }>;
  altlocs: Array<{
    chain: string;
    resseq: number;
    icode: string;
    resname: string;
    conformers: Array<{ id: string; occ: number; b_factor: number }>;
  }>;
  disulfides: Array<{
    cys1: string;
    cys2: string;
    resname1: string;
    resname2: string;
    distance_A: number;
    bonded: boolean;
  }>;
  heterogens: Array<{
    chain: string;
    resseq: number;
    icode: string;
    resname: string;
    category: 'ligand' | 'cofactor' | 'metal' | 'buffer';
    heavy_atoms: number;
    keep_recommended: boolean;
  }>;
  waters: {
    total: number;
    structural: number;
    bulk: number;
  };
  insertion_codes: Array<{ chain: string; resseq: number; icode: string; resname: string }>;
  chains: Array<{ id: string; protein_residues: number; first_res: number | null; last_res: number | null }>;
  has_hydrogens: boolean;
  model_count: number;
}

export interface RamaPoint {
  chain: string;
  resseq: number;
  icode: string;
  resname: string;
  key: string;
  phi: number;
  psi: number;
  omega: number | null;
  category: 'general' | 'glycine' | 'proline' | 'pre_proline';
  status: 'favored' | 'allowed' | 'outlier';
}

export interface RotamerEntry {
  chain: string;
  resseq: number;
  icode: string;
  resname: string;
  key: string;
  chis: number[];
  nearest_rotamer: string;
  library_prob: number;
  rmsd_deg: number;
  status: 'favored' | 'allowed' | 'outlier';
}

export interface ClashEntry {
  atom1: { chain: string; resseq: number; icode: string; resname: string; atom: string; key: string };
  atom2: { chain: string; resseq: number; icode: string; resname: string; atom: string; key: string };
  distance: number;
  overlap: number;
  severity: 'severe' | 'moderate';
  midpoint: [number, number, number];
}

export interface ResidueIndexItem {
  chain: string;
  resseq: number;
  icode: string;
  resname: string;
  std_resname: string;
  key: string;
  phi: number | null;
  psi: number | null;
  omega: number | null;
  rama_status: string | null;
  rotamer_status: string | null;
  chis: number[];
  clash_count: number;
  missing_atoms: string[];
  altlocs: string[];
  mean_b: number;
  min_occ: number;
  issues: string[];
}

export interface ValidationData {
  quality_score: number;
  ramachandran: {
    favored_pct: number;
    allowed_pct: number;
    outlier_pct: number;
    counts: { favored: number; allowed: number; outlier: number };
    points: RamaPoint[];
  };
  rotamers: {
    favored_pct: number;
    allowed_pct: number;
    outlier_pct: number;
    counts: { favored: number; allowed: number; outlier: number };
    outliers: RotamerEntry[];
    all: RotamerEntry[];
  };
  geometry: {
    bond_rmsz: number;
    angle_rmsz: number;
    bond_outliers: Array<{
      chain: string;
      resseq: number;
      icode: string;
      resname: string;
      key: string;
      atoms: string;
      observed: number;
      ideal: number;
      z_score: number;
    }>;
    angle_outliers: Array<{
      chain: string;
      resseq: number;
      icode: string;
      resname: string;
      key: string;
      atoms: string;
      observed: number;
      ideal: number;
      z_score: number;
    }>;
    cbeta_outliers: Array<{ chain: string; resseq: number; icode: string; resname: string; key: string; deviation: number }>;
    cis_peptides: Array<{
      chain: string;
      resseq: number;
      icode: string;
      prev_resseq: number;
      prev_resname: string;
      resname: string;
      key: string;
      omega: number;
      type: string;
      severity: 'info' | 'outlier';
    }>;
    chiral_outliers: Array<{
      chain: string;
      resseq: number;
      icode: string;
      resname: string;
      key: string;
      chiral_volume: number;
      expected: string;
      diagnosis: string;
    }>;
  };
  clashes: {
    clashscore: number;
    total_clashes: number;
    total_atoms: number;
    clashes: ClashEntry[];
    residue_clash_counts: Record<string, number>;
  };
  nqh_flips: Array<{
    chain: string;
    resseq: number;
    icode: string;
    resname: string;
    key: string;
    current_score: number;
    flipped_score: number;
    hbond_gain: number;
    reason: string;
  }>;
  bfactor_profile: Array<{
    chain: string;
    resseq: number;
    icode: string;
    resname: string;
    key: string;
    mean_b: number;
    min_occ: number;
  }>;
  residues: ResidueIndexItem[];
}

export interface PropkaResidue {
  chain: string;
  resseq: number;
  icode: string;
  resname: string;
  key: string;
  base_key: string;
  model_pka: number;
  pka: number;
  delta_pka: number;
  buried_pct: number;
  assigned_state: string;
  charge_at_ph: number;
  rationale: string;
  is_disulfide: boolean;
  metal_coordination: string[];
  is_override: boolean;
  shifted_significantly: boolean;
  available_states: string[];
}

export interface PropkaData {
  engine: string;
  ph: number;
  isoelectric_point: number;
  net_charge_discrete: number;
  net_charge_continuous: number;
  titratable_count: number;
  shifted_count: number;
  metal_coordinated_count: number;
  residues: PropkaResidue[];
  titration_curve: Array<{ ph: number; net_charge: number }>;
}

export interface DockingBox {
  ligand_id: string;
  label: string;
  center: [number, number, number];
  size: [number, number, number];
}

export interface LigandItem {
  id: string;
  chain: string;
  resseq: number;
  icode: string;
  resname: string;
  category: 'ligand' | 'cofactor';
  heavy_atoms: number;
  mean_b: number;
  center: [number, number, number];
  docking_box: DockingBox;
  contacts: Array<{
    chain: string;
    resseq: number;
    icode: string;
    resname: string;
    key: string;
    protein_atom: string;
    ligand_atom: string;
    distance_A: number;
    type: 'hbond' | 'hydrophobic' | 'vdw';
  }>;
  formula?: string;
  smiles?: string;
  mw?: number;
  clogp?: number;
  tpsa?: number;
  hbd?: number;
  hba?: number;
  rotatable_bonds?: number;
  rings?: number;
  lipinski_violations?: number;
  drug_like?: boolean;
  svg?: string;
}

export interface MetalItem {
  id: string;
  chain: string;
  resseq: number;
  resname: string;
  element: string;
  center: [number, number, number];
  coordination_number: number;
  geometry: string;
  coordinators: Array<{
    chain: string;
    resseq: number;
    resname: string;
    atom: string;
    key: string;
    distance_A: number;
  }>;
}

export interface LigandsInfo {
  ligands: LigandItem[];
  metals: MetalItem[];
  buffers: Array<{ id: string; chain: string; resseq: number; resname: string; heavy_atoms: number }>;
  docking_boxes: DockingBox[];
}

export interface PocketItem {
  id: string;
  rank: number;
  label: string;
  volume_A3: number;
  druggability_score: number;
  hydrophobic_pct: number;
  center: [number, number, number];
  size: [number, number, number];
  lining_residues: Array<{
    chain: string;
    resseq: number;
    icode: string;
    resname: string;
    key: string;
    hydrophobic: boolean;
  }>;
  spheres: Array<[number, number, number]>;
  probe_radius: number;
}

export interface SequenceResidue {
  chain: string;
  resseq: number;
  icode: string;
  resname: string;
  aa1: string;
  key: string;
  severity: 'clean' | 'warn' | 'error' | 'pocket' | 'metal';
  issues: string[];
}

export interface OverlayBond {
  label: string;
  p1: [number, number, number];
  p2: [number, number, number];
}

export interface SolvationEstimate {
  protein_span_nm: [number, number, number];
  box_diameter_nm: number;
  box_volume_nm3: number;
  estimated_waters: number;
  na_ions: number;
  cl_ions: number;
  net_charge_neutralized?: number;
  total_atoms_est: number;
}

export interface BaselineMetrics {
  overall_score: number;
  md_score: number;
  docking_score: number;
  quality_score: number;
  clashscore: number;
  issues_count: number;
  atoms: number;
  hydrogens: number;
  rama_favored_pct: number;
  rotamer_outliers: number;
}

export interface GromacsTestResult {
  available: boolean;
  passed: boolean;
  return_code?: number;
  force_field: string;
  water_model: string;
  total_mass_amu?: number | null;
  total_charge_e?: number | null;
  topology_atoms?: number | null;
  fatal_error?: string | null;
  message: string;
  log?: string;
}

export interface StructureSessionState {
  session_id: string;
  name: string;
  source: string;
  pdb_text: string;
  can_undo: boolean;
  can_redo: boolean;
  undo_count: number;
  redo_count: number;
  history: string[];
  baseline?: BaselineMetrics;
  metadata: StructureMetadata;
  pdbfixer: PdbFixerDiagnostics;
  validation: ValidationData;
  propka: PropkaData;
  ligands_info: LigandsInfo;
  pockets?: PocketItem[];
  sequence_chains?: Record<string, SequenceResidue[]>;
  overlays?: {
    disulfides: OverlayBond[];
    metal_bonds: OverlayBond[];
    ligand_hbonds: OverlayBond[];
  };
  solvation_estimate?: SolvationEstimate;
  readiness: {
    overall_score: number;
    md: ReadinessBlock;
    docking: ReadinessBlock;
    crystallographic: ReadinessBlock;
  };
}

export interface RotamerOption {
  index: number;
  name: string;
  probability: number;
  chis: number[];
  clashes: number;
  clash_energy: number;
  score: number;
}

export interface RotamerExplorerData {
  key: string;
  resname: string;
  std_resname: string;
  current: {
    chis: number[];
    nearest_rotamer: string;
    library_prob: number;
    rmsd_deg: number;
    status: string;
  } | null;
  rotamers: RotamerOption[];
}

export type MoleColorBy = 'id' | 'hydropathy' | 'polarity' | 'charge';

export interface MoleProfilePoint {
  distance: number;
  x: number;
  y: number;
  z: number;
  radius: number;
  hydropathy: number;
  polarity: number;
  charge: number;
  mutability: number;
}

export interface MoleLiningResidue {
  chain: string;
  resseq: number;
  icode: string;
  resname: string;
  std_resname: string;
  key: string;
  min_dist_A: number;
  start_dist: number;
  end_dist: number;
  backbone_contact: boolean;
  hydropathy: number;
  charge: number;
  is_bottleneck: boolean;
}

export interface MoleChannel {
  id: string;
  kind: 'tunnel' | 'pore';
  name: string;
  origin_label: string;
  color: string;
  length_A: number;
  bottleneck_radius_A: number;
  bottleneck_distance_A: number;
  bottleneck_point: [number, number, number];
  bottleneck_residues: string[];
  max_radius_A: number;
  tortuosity: number;
  throughput: number;
  properties: {
    hydropathy: number;
    polarity: number;
    mutability: number;
    net_charge: number;
    positive_residues: number;
    negative_residues: number;
    character: string;
  };
  profile: MoleProfilePoint[];
  lining_residues: MoleLiningResidue[];
}

export interface MoleAnalysisData {
  tunnels: MoleChannel[];
  pores: MoleChannel[];
  origins: Array<{ label: string; point: [number, number, number] }>;
  params: {
    mode: 'tunnels' | 'pores' | 'both';
    bottleneck_radius: number;
    probe_radius: number;
    min_length: number;
    max_tunnels?: number;
    start_mode?: string;
    include_hetatms?: boolean;
  };
  summary: {
    total_tunnels: number;
    total_pores: number;
    voronoi_nodes: number;
    widest_bottleneck_A?: number;
  };
}

export interface MoleMutationCandidate {
  chain: string;
  resseq: number;
  icode: string;
  key: string;
  wt_resname: string;
  mut_resname: string;
  mutation: string;
  baseline_bottleneck_A: number;
  new_bottleneck_A: number;
  delta_bottleneck_A: number;
  effect: string;
}

export interface MoleMutationScanResult {
  target_aa: string;
  baseline_bottleneck_A: number;
  candidates: MoleMutationCandidate[];
}

export interface StructuralDiffData {
  ref_label: string;
  matched_ca: number;
  matched_heavy: number;
  ca_rmsd_A: number;
  backbone_rmsd_A: number;
  all_heavy_rmsd_A: number;
  max_shift_A: number;
  top_shifted: Array<{
    chain: string;
    resseq: number;
    resname: string;
    key: string;
    ca_shift_A: number;
    max_atom_shift_A: number;
    mean_heavy_shift_A: number;
  }>;
  per_residue: Array<{
    chain: string;
    resseq: number;
    resname: string;
    key: string;
    ca_shift_A: number;
    max_atom_shift_A: number;
    mean_heavy_shift_A: number;
  }>;
  aligned_ref_pdb: string;
}

export interface MembraneInfo {
  normal: [number, number, number];
  center: [number, number, number];
  top_leaflet_center: [number, number, number];
  bottom_leaflet_center: [number, number, number];
  thickness_A: number;
  disc_radius_A: number;
  tilt_from_z_deg: number;
  belt_hydrophobic_pct: number;
  flank_hydrophobic_pct: number;
  classification: string;
}

export interface InteractionsMembraneData {
  salt_bridges: Array<{
    acid: string;
    acid_chain: string;
    acid_resseq: number;
    acid_atom: string;
    base: string;
    base_chain: string;
    base_resseq: number;
    base_atom: string;
    distance_A: number;
    inter_chain: boolean;
  }>;
  pi_stacking: Array<{
    res1: string;
    chain1: string;
    resseq1: number;
    res2: string;
    chain2: string;
    resseq2: number;
    distance_A: number;
    angle_deg: number;
    geometry: string;
    inter_chain: boolean;
  }>;
  interface_contacts: Array<{
    res1: string;
    chain1: string;
    resseq1: number;
    res2: string;
    chain2: string;
    resseq2: number;
    distance_A: number;
  }>;
  counts: {
    salt_bridges: number;
    inter_chain_salt_bridges: number;
    pi_stacking: number;
    interface_contacts: number;
  };
  membrane: MembraneInfo | null;
}

export interface BatchPrepResult {
  batch_id: string;
  total: number;
  succeeded: number;
  failed: number;
  preset: string;
  ph: number;
  zip_size_bytes: number;
  rows: Array<{
    filename: string;
    status: 'ok' | 'error';
    residues?: number;
    before_readiness?: number;
    after_readiness?: number;
    before_clashscore?: number;
    after_clashscore?: number;
    before_issues?: number;
    after_issues?: number;
    hydrogens_added?: number;
    net_charge?: number;
    pI?: number;
    mole_channels?: number;
    widest_bottleneck_A?: number;
    error?: string | null;
  }>;
}

export interface HoleProfilePoint {
  axis_coord_A: number;
  radius_A: number;
  x: number;
  y: number;
  z: number;
  color: string;
  state: 'closed' | 'single_file' | 'hydrated';
  closest_residue: string;
  closest_atom: string;
  lining_residues: string[];
}

export interface HoleLiningResidue {
  chain: string;
  resseq: number;
  icode: string;
  resname: string;
  key: string;
  min_pore_radius_A: number;
  axis_min_A: number;
  axis_max_A: number;
  closest_atom: string;
  backbone_contact: boolean;
  is_constriction: boolean;
}

export interface HoleAnalysisData {
  axis_label: string;
  cvect: [number, number, number];
  cpoint: [number, number, number];
  params?: {
    axis_mode: 'auto' | 'z' | 'y' | 'x';
    step: number;
    end_radius: number;
    include_hetatms: boolean;
  };
  profile: HoleProfilePoint[];
  lining_residues: HoleLiningResidue[];
  summary: {
    length_A: number;
    min_radius_A: number;
    max_radius_A: number;
    mean_radius_A: number;
    constriction_coord_A: number;
    constriction_point: [number, number, number];
    constriction_residues: string[];
    gate_status: string;
    conductance_pS: number;
    pore_resistance_GOhm: number;
    red_fraction_pct: number;
    green_fraction_pct: number;
    blue_fraction_pct: number;
  };
}

export interface DockingPresetLigand {
  id: string;
  name: string;
  smiles: string;
}

export interface PreparedDockingLigand {
  ligand_key: string;
  name: string;
  source: string;
  smiles: string;
  formula: string;
  mw: number;
  clogp: number;
  tpsa: number;
  hbd: number;
  hba: number;
  rotatable_bonds: number;
  heavy_atoms: number;
  total_atoms: number;
  conformers_generated: number;
  svg: string;
  ph: number;
  vina_available: boolean;
}

export interface DockingScoreTerms {
  gauss1: number;
  gauss2: number;
  repulsion: number;
  hydrophobic: number;
  hbond: number;
  torsion: number;
}

export interface DockingHBond {
  residue: string;
  chain: string;
  resseq: number;
  rec_atom: string;
  lig_atom: string;
  distance_A: number;
  p1: [number, number, number];
  p2: [number, number, number];
}

export interface DockingClash {
  residue: string;
  distance_A: number;
  p1: [number, number, number];
  p2: [number, number, number];
}

export interface DockingContact {
  residue: string;
  chain: string;
  resseq: number;
  resname: string;
  key: string;
  type: 'hbond' | 'hydrophobic' | 'vdw';
  distance_A: number;
}

export interface DockingFrameEvent {
  step: number;
  total_steps: number;
  phase: string;
  temperature: number;
  accepted: boolean;
  is_new_best: boolean;
  score: number;
  best_score: number;
  ki: string;
  terms: DockingScoreTerms;
  hbonds: DockingHBond[];
  clashes: DockingClash[];
  contacts: DockingContact[];
  pose_pdb: string;
}

export interface DockingMode {
  mode: number;
  affinity_kcal_mol: number;
  ki: string;
  rmsd_from_best_A: number;
  terms: DockingScoreTerms;
  hbonds: DockingHBond[];
  clashes: DockingClash[];
  contacts: DockingContact[];
  pose_pdb: string;
}

export interface DockingLiveOverlay {
  boxCenter: [number, number, number];
  boxSize: [number, number, number];
  showBox: boolean;
  posePdb: string | null;
  score: number | null;
  isNewBest?: boolean;
  running?: boolean;
  hbonds: DockingHBond[];
  clashes: DockingClash[];
  label?: string;
}
