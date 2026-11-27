from __future__ import annotations

from typing import Any, Dict, List, Optional
from pydantic import BaseModel, Field


class FetchRCSBRequest(BaseModel):
    pdb_id: str = Field(..., min_length=4, max_length=12)
    format: str = Field(default="pdb", pattern="^(pdb|cif)$")


class FetchAlphaFoldRequest(BaseModel):
    uniprot_id: str = Field(..., min_length=4, max_length=20)


class LoadSampleRequest(BaseModel):
    sample_id: str = Field(default="6LU7")


class PdbFixerOptions(BaseModel):
    replace_nonstandard: bool = True
    nonstandard_overrides: Dict[str, str] = Field(default_factory=dict)
    resolve_altlocs: str = Field(
        default="highest_occupancy",
        pattern="^(highest_occupancy|lowest_bfactor|keep_A|keep_all)$",
    )
    add_missing_atoms: bool = True
    add_missing_residues: bool = True
    max_loop_length: int = Field(default=12, ge=1, le=50)
    cap_termini: str = Field(default="none", pattern="^(none|ACE_NME|charged)$")
    assign_disulfides: bool = True
    disulfide_cutoff: float = Field(default=2.5, ge=1.5, le=3.5)
    remove_waters: str = Field(default="bulk", pattern="^(none|all|bulk)$")
    water_keep_distance: float = Field(default=3.5, ge=2.0, le=6.0)
    remove_crystallization_artifacts: bool = True
    keep_ligands: bool = True
    keep_cofactors: bool = True
    keep_metals: bool = True
    keep_chains: Optional[List[str]] = None
    renumber_residues: bool = False
    fix_insertion_codes: bool = True
    standardize_atom_names: bool = True


class PropkaOptions(BaseModel):
    ph: float = Field(default=7.4, ge=0.0, le=14.0)
    respect_metal_coordination: bool = True
    respect_disulfides: bool = True
    add_hydrogens: bool = True
    optimize_hbond_network: bool = True
    his_default: str = Field(default="auto", pattern="^(auto|HID|HIE|HIP)$")
    overrides: Dict[str, str] = Field(default_factory=dict)


class ResidueTarget(BaseModel):
    chain: str
    resseq: int
    icode: str = ""


class MutateRequest(ResidueTarget):
    target_resname: str = Field(..., min_length=3, max_length=3)


class RotamerApplyRequest(ResidueTarget):
    rotamer_index: Optional[int] = None
    chi_angles: Optional[List[float]] = None


class AltlocPickRequest(ResidueTarget):
    altloc: str = Field(..., min_length=1, max_length=1)


class RegularizeZoneRequest(BaseModel):
    chain: str
    start_resseq: int
    end_resseq: int
    iterations: int = Field(default=60, ge=10, le=300)


class AutoFixCootRequest(BaseModel):
    fix_clashes: bool = True
    fix_rotamers: bool = True
    fix_nqh_flips: bool = True
    fix_cis_nonpro: bool = True
    fix_ramachandran: bool = True
    regularize_outliers: bool = True


class ExportRequest(BaseModel):
    target: str = Field(
        default="pdb",
        pattern="^(pdb|mmcif|gro|pdbqt|rec_h_pdb|pqr|gromacs_bundle|ligands_sdf|provenance_json|report_json|report_md)$",
    )
    force_field: str = "amber99sb-ildn"
    water_model: str = "tip3p"
    box_type: str = "dodecahedron"
    box_distance_nm: float = Field(default=1.0, ge=0.5, le=3.0)
    salt_conc_m: float = Field(default=0.15, ge=0.0, le=2.0)
    ph: float = Field(default=7.4, ge=0.0, le=14.0)
    docking_padding: float = Field(default=6.0, ge=2.0, le=15.0)
    selected_ligand_id: Optional[str] = None


class DeepPrepPresetRequest(BaseModel):
    preset: str = Field(
        default="md_production",
        pattern="^(md_production|docking_receptor|crystallographic_polish)$",
    )
    ph: float = Field(default=7.4, ge=0.0, le=14.0)


class RenameFileRequest(BaseModel):
    filename: str = Field(..., min_length=1, max_length=120)
    title: Optional[str] = None


class RenameChainsRequest(BaseModel):
    mapping: Dict[str, str]


class ReorderChainsRequest(BaseModel):
    chain_order: List[str]


class MergeChainsRequest(BaseModel):
    source_chain: str
    target_chain: str
    renumber_merged: bool = True


class DeleteChainRequest(BaseModel):
    chain_id: Optional[str] = None
    chain_ids: Optional[List[str]] = None


class RenumberChainRequest(BaseModel):
    chain_id: str
    mode: str = Field(default="start", pattern="^(start|offset)$")
    start_num: int = Field(default=1, ge=1, le=99999)
    offset: int = Field(default=0, ge=-9999, le=9999)
    include_heterogens: bool = False
    clear_icodes: bool = True


class RenameResidueCodeRequest(BaseModel):
    old_resname: str = Field(..., min_length=1, max_length=4)
    new_resname: str = Field(..., min_length=1, max_length=4)
    chain_id: Optional[str] = None
    resseq: Optional[int] = None
    move_to_chain: Optional[str] = None


class NormalizeFieldsRequest(BaseModel):
    reset_occupancy: bool = False
    reset_bfactor: Optional[float] = None
    strip_waters_to_chain: Optional[str] = None


class MoleAnalyzeRequest(BaseModel):
    mode: str = Field(default="both", pattern="^(tunnels|pores|both)$")
    bottleneck_radius: float = Field(default=0.9, ge=0.4, le=3.5)
    probe_radius: float = Field(default=4.2, ge=2.0, le=8.0)
    min_length: float = Field(default=5.0, ge=2.0, le=30.0)
    max_tunnels: int = Field(default=6, ge=1, le=12)
    start_mode: str = Field(default="auto", pattern="^(auto|ligand|pocket|residues|coordinates)$")
    start_residues: Optional[List[str]] = None
    start_point: Optional[List[float]] = None
    include_hetatms: bool = False


class MoleExportRequest(MoleAnalyzeRequest):
    format: str = Field(default="csv", pattern="^(csv|residues_csv|json|pdb)$")


class HoleAnalyzeRequest(BaseModel):
    axis_mode: str = Field(default="auto", pattern="^(auto|z|y|x)$")
    step: float = Field(default=0.5, ge=0.25, le=2.0)
    end_radius: float = Field(default=10.0, ge=4.5, le=18.0)
    custom_point: Optional[List[float]] = None
    include_hetatms: bool = False


class HoleExportRequest(HoleAnalyzeRequest):
    format: str = Field(default="csv", pattern="^(csv|pdb|json)$")


class SwitchAssemblyRequest(BaseModel):
    assembly_id: str = Field(default="ASU")
    deduplicate_axial_ions: bool = True


class RebuildGapRequest(BaseModel):
    chain: str
    start_resseq: int
    end_resseq: int
    sequence_3: Optional[List[str]] = None
    strict_geometry: bool = True
    max_seeds: int = Field(default=3, ge=1, le=5)


class UpdateDockingPrefsRequest(BaseModel):
    residue_roles: Optional[Dict[str, str]] = None
    chain_roles: Optional[Dict[str, str]] = None
    ion_charges: Optional[Dict[str, float]] = None
    ion_substitutions: Optional[Dict[str, Optional[str]]] = None
    include_nucleic: Optional[bool] = None
    reset_all: bool = False


class RcsbLigandLookupRequest(BaseModel):
    pdb_id: Optional[str] = None
    force_refresh: bool = False


class VerifyProvenanceRequest(BaseModel):
    provenance: Dict[str, Any]
