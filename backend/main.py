"""
PDBPrep Studio — FastAPI Backend Server
Provides full REST API for:
- Structure Import (Upload PDB/mmCIF, RCSB PDB Fetch, AlphaFold DB Fetch, Curated Showcase)
- PDBFixer Deep Diagnostics & Repair
- PROPKA 3 pKa Analysis, Titration Curves & Protonation / H-Bond Network Optimization
- COOT Crystallographic Validation & Interactive Residue Editing (Mutate, Rotamer Picker, Auto-Fit, Pepflip, Cis↔Trans, NQH Flip, AltLoc Picker, Chirality Fix, Zone Regularization, Auto-Polish)
- Ligand, Cofactor & Metal Coordination Analysis + RDKit 2D SVG & 3D ETKDGv3 Conformers
- Multi-Target Readiness Audit & Export (PDB, mmCIF, PDBQT + Vina Config, PQR, GROMACS MD Bundle ZIP, Reports)
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, List

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, PlainTextResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles

from .engines.analysis_engine import (
    analyze_membrane_and_interactions,
    compute_structural_diff,
    orient_structure_along_z,
)
from .engines.batch_engine import get_batch_zip, run_batch_preparation
from .engines.coot_engine import (
    apply_rotamer,
    auto_fix_coot,
    delete_residue,
    fix_chirality,
    flip_cis_trans,
    flip_nqh_residue,
    get_residue_rotamers,
    mutate_residue,
    pepflip_residue,
    pick_altloc,
    regularize_zone,
)
from .engines.docking_engine import (
    PRESET_LIGANDS,
    export_docked_result,
    find_vina_binary,
    prepare_docking_ligand,
    stream_live_docking,
)
from .engines.export_engine import (
    build_gromacs_zip_bundle,
    check_external_tools,
    estimate_solvation_system,
    generate_markdown_report,
    generate_pqr,
    generate_receptor_pdbqt,
    generate_vina_config,
    run_gromacs_preflight_test,
)
from .engines.hole_engine import (
    analyze_hole_pore,
    export_hole_csv,
    export_hole_pdb,
)
from .engines.ligand_engine import export_ligand_format, generate_ligand_conformers_sdf
from .engines.mole_engine import (
    analyze_tunnels_and_pores,
    export_mole_csv,
    export_mole_pdb,
    export_mole_residues_csv,
    scan_bottleneck_mutations,
)
from .engines.pdbfixer_engine import run_pdbfixer
from .engines.propka_engine import analyze_propka, apply_protonation
from .engines.rename_engine import (
    auto_standardize_chain_ids,
    canonicalize_residue_names,
    delete_chain,
    delete_chains_batch,
    merge_chains,
    normalize_structure_fields,
    rename_chains_batch,
    rename_residue_code,
    renumber_chain,
    reorder_chains,
)
from .engines.structure_io import (
    build_showcase_structure,
    fetch_alphafold_structure,
    fetch_rcsb_structure,
    parse_structure_text,
    structure_to_gro,
    structure_to_mmcif,
    structure_to_pdb,
)
from .models import (
    AltlocPickRequest,
    AutoFixCootRequest,
    DeepPrepPresetRequest,
    DeleteChainRequest,
    ExportRequest,
    FetchAlphaFoldRequest,
    FetchRCSBRequest,
    HoleAnalyzeRequest,
    HoleExportRequest,
    LoadSampleRequest,
    MergeChainsRequest,
    MoleAnalyzeRequest,
    MoleExportRequest,
    MutateRequest,
    NormalizeFieldsRequest,
    PdbFixerOptions,
    PropkaOptions,
    RegularizeZoneRequest,
    RenameChainsRequest,
    RenameFileRequest,
    RenameResidueCodeRequest,
    RenumberChainRequest,
    ReorderChainsRequest,
    ResidueTarget,
    RotamerApplyRequest,
)
from .session_manager import STORE

app = FastAPI(
    title="PDBPrep Studio API",
    version="2.0.0",
    description="Deep Protein Preparation, PDBFixer, PROPKA 3 & COOT Crystallographic Workbench",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


def _get_session(sid: str):
    try:
        return STORE.get(sid)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Structure session '{sid}' not found")


# ---------------------------------------------------------------------------
# System & Sessions
# ---------------------------------------------------------------------------

@app.get("/api/health")
def health() -> Dict[str, Any]:
    return {
        "status": "ok",
        "version": "2.0.0",
        "engines": check_external_tools(),
        "sessions": STORE.list_sessions(),
    }


@app.get("/api/sessions")
def list_sessions() -> Dict[str, Any]:
    return {"sessions": STORE.list_sessions()}


@app.get("/api/sessions/{sid}")
def get_session_state(sid: str) -> Dict[str, Any]:
    sess = _get_session(sid)
    return sess.get_full_state()


@app.post("/api/sessions/{sid}/undo")
def undo_session(sid: str) -> Dict[str, Any]:
    sess = _get_session(sid)
    action = sess.undo()
    return {"undone": action, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/redo")
def redo_session(sid: str) -> Dict[str, Any]:
    sess = _get_session(sid)
    action = sess.redo()
    return {"redone": action, "state": sess.get_full_state()}


# ---------------------------------------------------------------------------
# Structure Ingestion (Upload, RCSB, AlphaFold, Curated Showcase)
# ---------------------------------------------------------------------------

@app.post("/api/import/upload")
async def upload_structure(file: UploadFile = File(...)) -> Dict[str, Any]:
    try:
        raw = await file.read()
        text = raw.decode("utf-8", errors="replace")
        filename = file.filename or "uploaded.pdb"
        st = parse_structure_text(text, filename)
        sess = STORE.create(name=filename, source="Local Upload", st=st)
        return sess.get_full_state()
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Failed to parse structure: {exc}")


@app.post("/api/import/rcsb")
async def import_rcsb(req: FetchRCSBRequest) -> Dict[str, Any]:
    try:
        text, filename = await fetch_rcsb_structure(req.pdb_id, req.format)
        st = parse_structure_text(text, filename)
        sess = STORE.create(name=filename, source=f"RCSB PDB ({req.pdb_id.upper()})", st=st)
        return sess.get_full_state()
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@app.post("/api/import/alphafold")
async def import_alphafold(req: FetchAlphaFoldRequest) -> Dict[str, Any]:
    try:
        text, filename = await fetch_alphafold_structure(req.uniprot_id)
        st = parse_structure_text(text, filename)
        sess = STORE.create(name=filename, source=f"AlphaFold DB ({req.uniprot_id.upper()})", st=st)
        return sess.get_full_state()
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@app.post("/api/import/sample")
async def import_sample(req: LoadSampleRequest) -> Dict[str, Any]:
    if req.sample_id.strip().upper() == "6LU7":
        text, filename = await fetch_rcsb_structure("6LU7", "pdb")
        st = parse_structure_text(text, filename)
        sess = STORE.create(name=filename, source="RCSB PDB (6LU7 · Default)", st=st)
        return sess.get_full_state()
    text, filename = build_showcase_structure(req.sample_id)
    st = parse_structure_text(text, filename)
    sess = STORE.create(name=filename, source="Built-in Metalloenzyme Showcase", st=st)
    return sess.get_full_state()


# ---------------------------------------------------------------------------
# One-Click Deep Preparation Presets (PDBFixer + COOT + PROPKA)
# ---------------------------------------------------------------------------

@app.post("/api/sessions/{sid}/deep-prep")
def run_deep_prep_preset(sid: str, req: DeepPrepPresetRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"Deep Prep Preset ({req.preset})")
    actions = []

    if req.preset == "md_production":
        fixer_opts = PdbFixerOptions(
            replace_nonstandard=True,
            resolve_altlocs="highest_occupancy",
            add_missing_atoms=True,
            add_missing_residues=True,
            cap_termini="none",
            assign_disulfides=True,
            remove_waters="bulk",
            remove_crystallization_artifacts=True,
            keep_ligands=True,
            keep_cofactors=True,
            keep_metals=True,
        )
        f_res = run_pdbfixer(sess.structure, fixer_opts.model_dump())
        actions.extend(f_res["actions"])
        c_res = auto_fix_coot(sess.structure)
        actions.extend(c_res["actions"])
        sess.ph = req.ph
        p_res = apply_protonation(sess.structure, PropkaOptions(ph=req.ph).model_dump())
        actions.extend(p_res["actions"])
        c_post = auto_fix_coot(sess.structure)
        actions.extend(c_post["actions"])

    elif req.preset == "docking_receptor":
        fixer_opts = PdbFixerOptions(
            replace_nonstandard=True,
            resolve_altlocs="highest_occupancy",
            add_missing_atoms=True,
            add_missing_residues=True,
            assign_disulfides=True,
            remove_waters="bulk",
            remove_crystallization_artifacts=True,
            keep_ligands=True,
            keep_cofactors=True,
            keep_metals=True,
        )
        f_res = run_pdbfixer(sess.structure, fixer_opts.model_dump())
        actions.extend(f_res["actions"])
        c_res = auto_fix_coot(sess.structure)
        actions.extend(c_res["actions"])
        sess.ph = req.ph
        p_res = apply_protonation(sess.structure, PropkaOptions(ph=req.ph).model_dump())
        actions.extend(p_res["actions"])
        c_post = auto_fix_coot(sess.structure)
        actions.extend(c_post["actions"])

    else:  # crystallographic_polish
        c_res = auto_fix_coot(sess.structure)
        actions.extend(c_res["actions"])

    sess.invalidate_cache()
    sess.record_history(f"Executed Deep Prep Preset '{req.preset}' ({len(actions)} steps)")
    return {"actions": actions, "state": sess.get_full_state()}


# ---------------------------------------------------------------------------
# PDBFixer Deep Pipeline
# ---------------------------------------------------------------------------

@app.post("/api/sessions/{sid}/pdbfixer/run")
def run_pdbfixer_endpoint(sid: str, opts: PdbFixerOptions) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo("PDBFixer Pipeline")
    result = run_pdbfixer(sess.structure, opts.model_dump())
    sess.invalidate_cache()
    sess.record_history(f"PDBFixer completed ({len(result['actions'])} operations)")
    return {"result": result, "state": sess.get_full_state()}


# ---------------------------------------------------------------------------
# PROPKA 3 & Protonation Engine
# ---------------------------------------------------------------------------

@app.post("/api/sessions/{sid}/propka/preview")
def preview_propka_endpoint(sid: str, opts: PropkaOptions) -> Dict[str, Any]:
    """Re-evaluate PROPKA predictions live at a new pH or with overrides without mutating coordinates."""
    sess = _get_session(sid)
    sess.ph = opts.ph
    sess.propka_overrides = dict(opts.overrides)
    sess.his_default = opts.his_default
    sess.respect_metal = opts.respect_metal_coordination
    sess.respect_disulfides = opts.respect_disulfides
    sess.invalidate_cache()
    analysis = analyze_propka(
        sess.structure,
        ph=opts.ph,
        respect_metal_coordination=opts.respect_metal_coordination,
        respect_disulfides=opts.respect_disulfides,
        his_default=opts.his_default,
        overrides=opts.overrides,
    )
    return {"propka": analysis}


@app.post("/api/sessions/{sid}/propka/apply")
def apply_propka_endpoint(sid: str, opts: PropkaOptions) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"PROPKA Protonation (pH {opts.ph})")
    sess.ph = opts.ph
    sess.propka_overrides = dict(opts.overrides)
    sess.his_default = opts.his_default
    sess.respect_metal = opts.respect_metal_coordination
    sess.respect_disulfides = opts.respect_disulfides
    res = apply_protonation(sess.structure, opts.model_dump())
    sess.invalidate_cache()
    sess.record_history(f"Applied PROPKA protonation at pH {opts.ph} (+{res['hydrogens_added']} H)")
    return {"result": res, "state": sess.get_full_state()}


# ---------------------------------------------------------------------------
# COOT Interactive Model-Building & Validation Actions
# ---------------------------------------------------------------------------

@app.get("/api/sessions/{sid}/coot/rotamers")
def coot_list_rotamers(sid: str, chain: str, resseq: int, icode: str = "") -> Dict[str, Any]:
    sess = _get_session(sid)
    try:
        return get_residue_rotamers(sess.structure, chain, resseq, icode)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@app.post("/api/sessions/{sid}/coot/rotamer")
def coot_apply_rotamer(sid: str, req: RotamerApplyRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"COOT Rotamer {req.chain}:{req.resseq}")
    try:
        res = apply_rotamer(sess.structure, req.chain, req.resseq, req.icode, req.rotamer_index, req.chi_angles)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/coot/mutate")
def coot_mutate(sid: str, req: MutateRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"COOT Mutate {req.chain}:{req.resseq} -> {req.target_resname}")
    try:
        res = mutate_residue(sess.structure, req.chain, req.resseq, req.target_resname, req.icode)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/coot/pepflip")
def coot_pepflip(sid: str, req: ResidueTarget) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"COOT Pepflip {req.chain}:{req.resseq}")
    try:
        res = pepflip_residue(sess.structure, req.chain, req.resseq, req.icode)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/coot/cis-trans")
def coot_cis_trans(sid: str, req: ResidueTarget) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"COOT Cis-Trans {req.chain}:{req.resseq}")
    try:
        res = flip_cis_trans(sess.structure, req.chain, req.resseq, req.icode)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/coot/nqh-flip")
def coot_nqh_flip(sid: str, req: ResidueTarget) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"COOT NQH Flip {req.chain}:{req.resseq}")
    try:
        res = flip_nqh_residue(sess.structure, req.chain, req.resseq, req.icode)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/coot/altloc")
def coot_altloc(sid: str, req: AltlocPickRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"COOT AltLoc {req.chain}:{req.resseq} -> {req.altloc}")
    try:
        res = pick_altloc(sess.structure, req.chain, req.resseq, req.altloc, req.icode)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/coot/fix-chirality")
def coot_fix_chirality(sid: str, req: ResidueTarget) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"COOT Fix Chirality {req.chain}:{req.resseq}")
    try:
        res = fix_chirality(sess.structure, req.chain, req.resseq, req.icode)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/coot/regularize")
def coot_regularize(sid: str, req: RegularizeZoneRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"COOT Regularize {req.chain}:{req.start_resseq}-{req.end_resseq}")
    try:
        res = regularize_zone(sess.structure, req.chain, req.start_resseq, req.end_resseq, req.iterations)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/coot/delete-residue")
def coot_delete_residue(sid: str, req: ResidueTarget) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"Delete {req.chain}:{req.resseq}")
    try:
        res = delete_residue(sess.structure, req.chain, req.resseq, req.icode)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/coot/auto-fix")
def coot_auto_fix(sid: str, req: AutoFixCootRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo("COOT Auto-Refine")
    res = auto_fix_coot(
        sess.structure,
        fix_clashes=req.fix_clashes,
        fix_rotamers=req.fix_rotamers,
        fix_nqh_flips=req.fix_nqh_flips,
        fix_cis_nonpro=req.fix_cis_nonpro,
        fix_ramachandran=req.fix_ramachandran,
        regularize_outliers=req.regularize_outliers,
    )
    sess.invalidate_cache()
    sess.record_history(
        f"COOT Auto-Refine: Quality {res['before_score']} → {res['after_score']} (Clashscore {res['before_clashscore']} → {res['after_clashscore']})"
    )
    return {"result": res, "state": sess.get_full_state()}


# ---------------------------------------------------------------------------
# Structure, File, Chain & Residue Renaming Studio
# ---------------------------------------------------------------------------

@app.post("/api/sessions/{sid}/rename/file")
def rename_structure_file(sid: str, req: RenameFileRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    new_name = req.filename.strip()
    if not new_name.lower().endswith((".pdb", ".cif", ".mmcif", ".gro", ".ent")):
        new_name += ".pdb"
    old_name = sess.name
    sess.name = new_name
    sess.structure.name = Path(new_name).stem[:16]
    if req.title is not None and req.title.strip():
        sess.structure.info["_struct.title"] = req.title.strip()
    sess.invalidate_cache()
    msg = f"Renamed file {old_name} → {new_name}"
    sess.record_history(msg)
    return {"result": {"message": msg}, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/rename/chains")
def rename_chains_endpoint(sid: str, req: RenameChainsRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo("Rename Chains")
    try:
        res = rename_chains_batch(sess.structure, req.mapping)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    if not res.get("changed"):
        sess.undo_stack.pop()
    else:
        sess.invalidate_cache()
        sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/rename/chains/standardize")
def standardize_chains_endpoint(sid: str) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo("Standardize Chain IDs (A, B, C…)")
    try:
        res = auto_standardize_chain_ids(sess.structure)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    if not res.get("changed"):
        sess.undo_stack.pop()
    else:
        sess.invalidate_cache()
        sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/rename/chains/reorder")
def reorder_chains_endpoint(sid: str, req: ReorderChainsRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo("Reorder Chains")
    try:
        res = reorder_chains(sess.structure, req.chain_order)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/rename/chains/merge")
def merge_chains_endpoint(sid: str, req: MergeChainsRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"Merge Chain {req.source_chain} → {req.target_chain}")
    try:
        res = merge_chains(sess.structure, req.source_chain, req.target_chain, req.renumber_merged)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/rename/chains/delete")
def delete_chain_endpoint(sid: str, req: DeleteChainRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    ids = list(req.chain_ids or [])
    if req.chain_id:
        ids.append(req.chain_id)
    if not ids:
        raise HTTPException(status_code=400, detail="Provide chain_id or chain_ids to delete")
    sess.push_undo(f"Delete Chain(s) {', '.join(ids)}")
    try:
        res = delete_chains_batch(sess.structure, ids)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/rename/renumber")
def renumber_chain_endpoint(sid: str, req: RenumberChainRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"Renumber Chain {req.chain_id}")
    try:
        res = renumber_chain(
            sess.structure,
            req.chain_id,
            req.mode,
            req.start_num,
            req.offset,
            req.include_heterogens,
            req.clear_icodes,
        )
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/rename/residue-code")
def rename_residue_code_endpoint(sid: str, req: RenameResidueCodeRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo(f"Rename Residue {req.old_resname} → {req.new_resname}")
    try:
        res = rename_residue_code(
            sess.structure,
            req.old_resname,
            req.new_resname,
            req.chain_id,
            req.resseq,
            req.move_to_chain,
        )
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/rename/canonicalize")
def canonicalize_names_endpoint(sid: str) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo("Canonicalize Residue Names")
    res = canonicalize_residue_names(sess.structure)
    if not res.get("changed"):
        sess.undo_stack.pop()
    else:
        sess.invalidate_cache()
        sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/rename/normalize")
def normalize_fields_endpoint(sid: str, req: NormalizeFieldsRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo("Normalize PDB Fields")
    res = normalize_structure_fields(
        sess.structure,
        reset_occupancy=req.reset_occupancy,
        reset_bfactor=req.reset_bfactor,
        strip_waters_to_chain=req.strip_waters_to_chain,
    )
    if not res.get("changed"):
        sess.undo_stack.pop()
    else:
        sess.invalidate_cache()
        sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


# ---------------------------------------------------------------------------
# Ligand Conformers
# ---------------------------------------------------------------------------

@app.post("/api/sessions/{sid}/ligands/conformers")
def ligand_conformers(sid: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    sess = _get_session(sid)
    lig_id = payload.get("ligand_id", "")
    num_confs = int(payload.get("num_conformers", 5))
    try:
        sdf, summary = generate_ligand_conformers_sdf(sess.structure, lig_id, num_confs)
        return {"ligand_id": lig_id, "conformers": summary, "sdf": sdf}
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@app.post("/api/sessions/{sid}/ligands/export")
def ligand_export_endpoint(sid: str, payload: Dict[str, Any]) -> Response:
    """Export a prepared ligand as .pdbqt (Vina/GNINA), .mol2 (CGenFF/ACPYPE), .sdf, or .pdb."""
    sess = _get_session(sid)
    lig_id = str(payload.get("ligand_id", ""))
    fmt = str(payload.get("format", "pdbqt")).lower()
    ph = float(payload.get("ph", sess.ph))
    if fmt not in ("pdbqt", "mol2", "sdf", "pdb"):
        raise HTTPException(status_code=400, detail=f"Unsupported ligand format '{fmt}'")
    try:
        text, filename = export_ligand_format(sess.structure, lig_id, fmt=fmt, ph=ph)
        return PlainTextResponse(
            text,
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))


# ---------------------------------------------------------------------------
# MOLE-Style Tunnels, Channels & Pores Analysis + Bottleneck Mutation Scanner
# ---------------------------------------------------------------------------

@app.post("/api/sessions/{sid}/mole/scan-mutations")
def mole_scan_mutations_endpoint(sid: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    sess = _get_session(sid)
    bn_res = payload.get("bottleneck_residues") or []
    bn_pt = payload.get("bottleneck_point") or [0.0, 0.0, 0.0]
    base_r = float(payload.get("baseline_bottleneck_A", 1.0))
    target_aa = str(payload.get("target_aa", "ALA"))
    include_het = bool(payload.get("include_hetatms", False))
    try:
        return scan_bottleneck_mutations(
            sess.structure,
            bottleneck_residues=bn_res,
            bottleneck_point=bn_pt,
            baseline_bottleneck_A=base_r,
            target_aa=target_aa,
            include_hetatms=include_het,
        )
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

@app.post("/api/sessions/{sid}/mole/analyze")
def mole_analyze_endpoint(sid: str, req: MoleAnalyzeRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    try:
        res = analyze_tunnels_and_pores(
            sess.structure,
            mode=req.mode,
            bottleneck_radius=req.bottleneck_radius,
            probe_radius=req.probe_radius,
            min_length=req.min_length,
            max_tunnels=req.max_tunnels,
            start_mode=req.start_mode,
            start_residues=req.start_residues,
            start_point=req.start_point,
            include_hetatms=req.include_hetatms,
        )
        return res
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"MOLE tunnel analysis failed: {exc}")


@app.post("/api/sessions/{sid}/mole/export")
def mole_export_endpoint(sid: str, req: MoleExportRequest) -> Response:
    sess = _get_session(sid)
    stem = Path(sess.name).stem
    data = analyze_tunnels_and_pores(
        sess.structure,
        mode=req.mode,
        bottleneck_radius=req.bottleneck_radius,
        probe_radius=req.probe_radius,
        min_length=req.min_length,
        max_tunnels=req.max_tunnels,
        start_mode=req.start_mode,
        start_residues=req.start_residues,
        start_point=req.start_point,
        include_hetatms=req.include_hetatms,
    )
    if req.format == "csv":
        return PlainTextResponse(
            export_mole_csv(data),
            media_type="text/csv",
            headers={"Content-Disposition": f'attachment; filename="{stem}_mole_tunnels.csv"'},
        )
    if req.format == "residues_csv":
        return PlainTextResponse(
            export_mole_residues_csv(data),
            media_type="text/csv",
            headers={"Content-Disposition": f'attachment; filename="{stem}_mole_residues.csv"'},
        )
    if req.format == "pdb":
        return PlainTextResponse(
            export_mole_pdb(data),
            headers={"Content-Disposition": f'attachment; filename="{stem}_mole_tunnels.pdb"'},
        )
    return Response(
        content=json.dumps(data, indent=2),
        media_type="application/json",
        headers={"Content-Disposition": f'attachment; filename="{stem}_mole_tunnels.json"'},
    )


# ---------------------------------------------------------------------------
# HOLE Ion-Channel Pore Radius & Conductance Engine (Smart et al. 1996)
# ---------------------------------------------------------------------------

@app.post("/api/sessions/{sid}/hole/analyze")
def hole_analyze_endpoint(sid: str, req: HoleAnalyzeRequest) -> Dict[str, Any]:
    sess = _get_session(sid)
    try:
        return analyze_hole_pore(
            sess.structure,
            axis_mode=req.axis_mode,
            step=req.step,
            end_radius=req.end_radius,
            custom_point=req.custom_point,
            include_hetatms=req.include_hetatms,
        )
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"HOLE pore analysis failed: {exc}")


@app.post("/api/sessions/{sid}/hole/export")
def hole_export_endpoint(sid: str, req: HoleExportRequest) -> Response:
    sess = _get_session(sid)
    stem = Path(sess.name).stem
    data = analyze_hole_pore(
        sess.structure,
        axis_mode=req.axis_mode,
        step=req.step,
        end_radius=req.end_radius,
        custom_point=req.custom_point,
        include_hetatms=req.include_hetatms,
    )
    if req.format == "csv":
        return PlainTextResponse(
            export_hole_csv(data),
            media_type="text/csv",
            headers={"Content-Disposition": f'attachment; filename="{stem}_hole_pore.csv"'},
        )
    if req.format == "pdb":
        return PlainTextResponse(
            export_hole_pdb(data),
            headers={"Content-Disposition": f'attachment; filename="{stem}_hole_pore.pdb"'},
        )
    return Response(
        content=json.dumps(data, indent=2),
        media_type="application/json",
        headers={"Content-Disposition": f'attachment; filename="{stem}_hole_pore.json"'},
    )


# ---------------------------------------------------------------------------
# Live Real-Time 3D Molecular Docking Studio (SSE Streaming Vina Engine)
# ---------------------------------------------------------------------------

@app.get("/api/docking/presets")
def docking_presets_endpoint() -> Dict[str, Any]:
    return {
        "presets": PRESET_LIGANDS,
        "vina_available": find_vina_binary() is not None,
    }


@app.post("/api/sessions/{sid}/docking/prepare-ligand")
async def docking_prepare_ligand_endpoint(
    sid: str,
    file: UploadFile = File(None),
    smiles: str = Form(""),
    session_ligand_id: str = Form(""),
    name: str = Form(""),
    ph: float = Form(7.4),
    num_conformers: int = Form(8),
) -> Dict[str, Any]:
    sess = _get_session(sid)
    try:
        file_text = None
        filename = name or None
        if file is not None:
            raw = await file.read()
            file_text = raw.decode("utf-8", errors="replace")
            filename = file.filename or "ligand.sdf"
        return prepare_docking_ligand(
            st=sess.structure,
            smiles=smiles.strip() or None,
            file_text=file_text,
            filename=filename,
            session_ligand_id=session_ligand_id.strip() or None,
            ph=ph,
            num_conformers=num_conformers,
        )
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Ligand preparation failed: {exc}")


@app.post("/api/sessions/{sid}/docking/strip-ligands")
def docking_strip_ligands_endpoint(sid: str) -> Dict[str, Any]:
    """Strip co-crystallized organic ligands (while keeping protein + metals) to open the binding pocket for docking."""
    sess = _get_session(sid)
    sess.push_undo("Strip Co-Crystallized Ligands for Docking")
    result = run_pdbfixer(
        sess.structure,
        {
            "replace_nonstandard": False,
            "resolve_altlocs": "keep_all",
            "add_missing_atoms": False,
            "add_missing_residues": False,
            "cap_termini": "none",
            "assign_disulfides": False,
            "remove_waters": "all",
            "remove_crystallization_artifacts": True,
            "keep_ligands": False,
            "keep_cofactors": True,
            "keep_metals": True,
        },
    )
    sess.invalidate_cache()
    sess.record_history("Stripped co-crystallized ligands & bulk waters to open pocket for live docking")
    return {"result": result, "state": sess.get_full_state()}


@app.post("/api/sessions/{sid}/docking/run-stream")
def docking_run_stream_endpoint(sid: str, payload: Dict[str, Any]) -> StreamingResponse:
    sess = _get_session(sid)
    ligand_key = str(payload.get("ligand_key", "")).strip()
    box_center = payload.get("box_center")
    box_size = payload.get("box_size") or [20.0, 20.0, 20.0]
    exhaustiveness = int(payload.get("exhaustiveness", 8))
    speed = str(payload.get("speed", "normal"))
    use_vina_verify = bool(payload.get("use_vina_verify", False))
    include_hetatms = bool(payload.get("include_hetatms", False))

    if not ligand_key or not box_center or len(box_center) != 3:
        raise HTTPException(status_code=400, detail="Missing ligand_key or 3D box_center [x, y, z].")

    gen = stream_live_docking(
        st=sess.structure,
        ligand_key=ligand_key,
        box_center=[float(x) for x in box_center],
        box_size=[float(x) for x in box_size],
        exhaustiveness=max(4, min(16, exhaustiveness)),
        speed=speed,
        use_vina_verify=use_vina_verify,
        include_hetatms=include_hetatms,
    )
    return StreamingResponse(
        gen,
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no",
        },
    )


@app.get("/api/docking/{result_id}/download")
def docking_download_endpoint(result_id: str, fmt: str = "complex_pdb", mode: int = 1) -> Response:
    try:
        content, filename, media_type = export_docked_result(result_id, fmt=fmt, mode=mode)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc))
    return PlainTextResponse(
        content,
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


# ---------------------------------------------------------------------------
# Structural Superposition (Before vs After RMSD), Interactions & Membrane Z-Align
# ---------------------------------------------------------------------------

@app.get("/api/sessions/{sid}/analysis/diff")
def analysis_before_after_diff(sid: str) -> Dict[str, Any]:
    sess = _get_session(sid)
    orig_pdb = sess.original_pdb or structure_to_pdb(sess.structure)
    try:
        st_orig = parse_structure_text(orig_pdb, sess.name)
        return compute_structural_diff(st_orig, sess.structure, ref_label=f"Original ({sess.name})")
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@app.post("/api/sessions/{sid}/analysis/align-upload")
async def analysis_align_upload(sid: str, file: UploadFile = File(...)) -> Dict[str, Any]:
    sess = _get_session(sid)
    try:
        raw = await file.read()
        text = raw.decode("utf-8", errors="replace")
        fname = file.filename or "reference.pdb"
        st_ref = parse_structure_text(text, fname)
        return compute_structural_diff(st_ref, sess.structure, ref_label=fname)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Alignment failed: {exc}")


@app.get("/api/sessions/{sid}/analysis/interactions")
def analysis_interactions_endpoint(sid: str) -> Dict[str, Any]:
    sess = _get_session(sid)
    return analyze_membrane_and_interactions(sess.structure)


@app.post("/api/sessions/{sid}/analysis/orient-z")
def analysis_orient_z_endpoint(sid: str, payload: Dict[str, Any] = None) -> Dict[str, Any]:
    sess = _get_session(sid)
    sess.push_undo("Orient Membrane / Pore Axis along +Z")
    custom_axis = (payload or {}).get("axis")
    try:
        res = orient_structure_along_z(sess.structure, custom_axis=custom_axis)
    except Exception as exc:
        sess.undo_stack.pop()
        raise HTTPException(status_code=400, detail=str(exc))
    sess.invalidate_cache()
    sess.record_history(res["message"])
    return {"result": res, "state": sess.get_full_state()}


# ---------------------------------------------------------------------------
# Batch Multi-Structure Preparation Queue
# ---------------------------------------------------------------------------

@app.post("/api/batch/run")
async def batch_run_endpoint(
    files: List[UploadFile] = File(...),
    preset: str = Form("md_production"),
    ph: float = Form(7.4),
    include_mole: bool = Form(False),
) -> Dict[str, Any]:
    if not files:
        raise HTTPException(status_code=400, detail="Upload at least one .pdb, .cif, or .gro file")
    loaded = []
    for f in files[:30]:
        raw = await f.read()
        loaded.append((f.filename or "structure.pdb", raw.decode("utf-8", errors="replace")))
    return run_batch_preparation(loaded, preset=preset, ph=ph, include_mole=include_mole)


@app.get("/api/batch/{batch_id}/download")
def batch_download_endpoint(batch_id: str) -> Response:
    try:
        zbytes, fname = get_batch_zip(batch_id)
        return Response(
            content=zbytes,
            media_type="application/zip",
            headers={"Content-Disposition": f'attachment; filename="{fname}"'},
        )
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc))


# ---------------------------------------------------------------------------
# Multi-Target Exports (PDB, mmCIF, PDBQT, PQR, GROMACS Bundle ZIP, Reports)
# ---------------------------------------------------------------------------

@app.post("/api/sessions/{sid}/export")
def export_structure(sid: str, req: ExportRequest) -> Response:
    sess = _get_session(sid)
    state = sess.get_full_state()
    stem = Path(sess.name).stem

    if req.target == "pdb":
        content = structure_to_pdb(sess.structure)
        return PlainTextResponse(
            content,
            headers={"Content-Disposition": f'attachment; filename="{stem}_prepared.pdb"'},
        )
    if req.target == "mmcif":
        content = structure_to_mmcif(sess.structure)
        return PlainTextResponse(
            content,
            headers={"Content-Disposition": f'attachment; filename="{stem}_prepared.cif"'},
        )
    if req.target == "gro":
        content = structure_to_gro(sess.structure)
        return PlainTextResponse(
            content,
            headers={"Content-Disposition": f'attachment; filename="{stem}_prepared.gro"'},
        )
    if req.target == "pdbqt":
        content = generate_receptor_pdbqt(sess.structure)
        return PlainTextResponse(
            content,
            headers={"Content-Disposition": f'attachment; filename="{stem}_receptor.pdbqt"'},
        )
    if req.target == "pqr":
        content = generate_pqr(sess.structure)
        return PlainTextResponse(
            content,
            headers={"Content-Disposition": f'attachment; filename="{stem}_prepared.pqr"'},
        )
    if req.target == "report_md":
        content = generate_markdown_report(
            sess.name,
            state["metadata"],
            state["readiness"],
            state["validation"],
            state["propka"],
            state["ligands_info"],
            state["history"],
        )
        return PlainTextResponse(
            content,
            headers={"Content-Disposition": f'attachment; filename="{stem}_report.md"'},
        )
    if req.target == "report_json":
        report_obj = {k: v for k, v in state.items() if k != "pdb_text"}
        return Response(
            content=json.dumps(report_obj, indent=2),
            media_type="application/json",
            headers={"Content-Disposition": f'attachment; filename="{stem}_audit.json"'},
        )
    if req.target == "gromacs_bundle":
        zip_bytes = build_gromacs_zip_bundle(
            sess.structure,
            sess.name,
            req.model_dump(),
            state["metadata"],
            state["readiness"],
            state["validation"],
            state["propka"],
            state["ligands_info"],
            state["history"],
        )
        return Response(
            content=zip_bytes,
            media_type="application/zip",
            headers={"Content-Disposition": f'attachment; filename="{stem}_pdbprep_bundle.zip"'},
        )

    raise HTTPException(status_code=400, detail=f"Unsupported export target: {req.target}")


@app.get("/api/sessions/{sid}/vina-config")
def get_vina_config(sid: str, ligand_id: str = "", padding: float = 6.0) -> Dict[str, Any]:
    sess = _get_session(sid)
    state = sess.get_full_state()
    cfg = generate_vina_config(
        state["ligands_info"].get("docking_boxes", []),
        ligand_id or None,
        padding,
    )
    return {"config": cfg}


@app.post("/api/sessions/{sid}/gromacs/test")
def test_gromacs_endpoint(sid: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    """Run live `gmx pdb2gmx` pre-flight verification + solvation/ion calculator."""
    sess = _get_session(sid)
    ff = str(payload.get("force_field", "amber99sb-ildn"))
    water = str(payload.get("water_model", "tip3p"))
    box_type = str(payload.get("box_type", "dodecahedron"))
    box_dist = float(payload.get("box_distance_nm", 1.0))
    salt = float(payload.get("salt_conc_m", 0.15))

    state = sess.get_full_state()
    solv = estimate_solvation_system(
        sess.structure,
        box_type=box_type,
        box_distance_nm=box_dist,
        salt_conc_m=salt,
        net_charge=state["propka"]["net_charge_discrete"],
    )
    gmx_res = run_gromacs_preflight_test(sess.structure, force_field=ff, water_model=water)
    return {"gromacs": gmx_res, "solvation": solv}


# Mount built frontend static assets (prefers live dev frontend/dist, falls back to packaged backend/static)
_DEV_FRONTEND_DIST = Path(__file__).resolve().parent.parent / "frontend" / "dist"
_PKG_FRONTEND_STATIC = Path(__file__).resolve().parent / "static"
_FRONTEND_DIR = _DEV_FRONTEND_DIST if (_DEV_FRONTEND_DIST / "index.html").exists() else _PKG_FRONTEND_STATIC
if (_FRONTEND_DIR / "index.html").exists():
    app.mount("/", StaticFiles(directory=str(_FRONTEND_DIR), html=True), name="frontend")
