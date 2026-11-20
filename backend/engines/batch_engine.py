"""
Batch Multi-Structure Preparation Engine:
Processes multiple uploaded .pdb / .cif / .gro structures sequentially through
PDBFixer + COOT Auto-Refine + PROPKA 3 + MOLE Tunnels, producing a live audit table
and a downloadable unified Batch ZIP bundle (prepared PDBs, .gro, .pdbqt, vina_config.txt, batch_summary.csv).
"""

from __future__ import annotations

import io
import uuid
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Tuple

from .coot_engine import analyze_structure_quality, auto_fix_coot
from .export_engine import (
    compute_readiness_audit,
    generate_receptor_pdbqt,
    generate_vina_config,
)
from .ligand_engine import analyze_ligands_and_metals
from .mole_engine import analyze_tunnels_and_pores
from .pdbfixer_engine import diagnose_pdbfixer, run_pdbfixer
from .propka_engine import analyze_propka, apply_protonation
from .structure_io import parse_structure_text, structure_to_gro, structure_to_pdb

_BATCH_STORE: Dict[str, Tuple[bytes, str]] = {}


def run_batch_preparation(
    files: List[Tuple[str, str]],
    preset: str = "md_production",
    ph: float = 7.4,
    include_mole: bool = False,
) -> Dict[str, Any]:
    """Process a batch of (filename, content) tuples and store downloadable ZIP bundle."""
    batch_id = uuid.uuid4().hex[:10]
    rows: List[Dict[str, Any]] = []
    csv_lines = [
        "filename,status,residues,before_readiness,after_readiness,before_clashscore,after_clashscore,"
        "before_issues,after_issues,hydrogens_added,net_charge_at_ph,isoelectric_point,mole_tunnels,widest_bottleneck_A"
    ]

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for fname, content in files:
            stem = Path(fname).stem
            try:
                st = parse_structure_text(content, fname)
                d0 = diagnose_pdbfixer(st)
                q0 = analyze_structure_quality(st)
                p0 = analyze_propka(st, ph=ph)
                l0 = analyze_ligands_and_metals(st)
                r0 = compute_readiness_audit(d0, q0, p0, l0)
                before_issues = sum(1 for r in q0["residues"] if r["issues"])

                h_added = 0
                if preset in ("md_production", "docking_receptor"):
                    run_pdbfixer(
                        st,
                        {
                            "replace_nonstandard": True,
                            "resolve_altlocs": "highest_occupancy",
                            "add_missing_atoms": True,
                            "add_missing_residues": True,
                            "assign_disulfides": True,
                            "remove_waters": "bulk",
                            "remove_crystallization_artifacts": True,
                            "keep_ligands": True,
                            "keep_cofactors": True,
                            "keep_metals": True,
                        },
                    )
                    auto_fix_coot(st)
                    pres = apply_protonation(
                        st,
                        {"ph": ph, "add_hydrogens": True, "optimize_hbond_network": True},
                    )
                    h_added = pres["hydrogens_added"]
                    auto_fix_coot(st)
                else:
                    auto_fix_coot(st)

                d1 = diagnose_pdbfixer(st)
                q1 = analyze_structure_quality(st)
                p1 = analyze_propka(st, ph=ph)
                l1 = analyze_ligands_and_metals(st)
                r1 = compute_readiness_audit(d1, q1, p1, l1)
                after_issues = sum(1 for r in q1["residues"] if r["issues"])

                mole_count = 0
                widest_bn = 0.0
                if include_mole:
                    mres = analyze_tunnels_and_pores(st, mode="both", max_tunnels=4)
                    mole_count = mres["summary"]["total_tunnels"] + mres["summary"]["total_pores"]
                    widest_bn = mres["summary"].get("widest_bottleneck_A", 0.0)

                # Write prepared files into unified ZIP
                zf.writestr(f"{stem}/{stem}_prepared.pdb", structure_to_pdb(st))
                zf.writestr(f"{stem}/{stem}_prepared.gro", structure_to_gro(st))
                zf.writestr(f"{stem}/{stem}_receptor.pdbqt", generate_receptor_pdbqt(st))
                zf.writestr(f"{stem}/{stem}_vina_config.txt", generate_vina_config(l1.get("docking_boxes", [])))

                row = {
                    "filename": fname,
                    "status": "ok",
                    "residues": len(q1["residues"]),
                    "before_readiness": r0["overall_score"],
                    "after_readiness": r1["overall_score"],
                    "before_clashscore": q0["clashes"]["clashscore"],
                    "after_clashscore": q1["clashes"]["clashscore"],
                    "before_issues": before_issues,
                    "after_issues": after_issues,
                    "hydrogens_added": h_added,
                    "net_charge": p1["net_charge_discrete"],
                    "pI": p1["isoelectric_point"],
                    "mole_channels": mole_count,
                    "widest_bottleneck_A": widest_bn,
                    "error": None,
                }
                rows.append(row)
                csv_lines.append(
                    f"{fname},ok,{row['residues']},{row['before_readiness']},{row['after_readiness']},"
                    f"{row['before_clashscore']},{row['after_clashscore']},{before_issues},{after_issues},"
                    f"{h_added},{row['net_charge']},{row['pI']},{mole_count},{widest_bn}"
                )
            except Exception as exc:
                rows.append({
                    "filename": fname,
                    "status": "error",
                    "error": str(exc),
                })
                csv_lines.append(f"{fname},error,0,0,0,0,0,0,0,0,0,0,0,0")

        zf.writestr("batch_summary.csv", "\n".join(csv_lines) + "\n")

    zip_bytes = buf.getvalue()
    _BATCH_STORE[batch_id] = (zip_bytes, f"pandaprotein_batch_{batch_id}.zip")
    # Keep store bounded to last 10 batch runs
    if len(_BATCH_STORE) > 10:
        oldest = next(iter(_BATCH_STORE))
        del _BATCH_STORE[oldest]

    return {
        "batch_id": batch_id,
        "total": len(files),
        "succeeded": sum(1 for r in rows if r["status"] == "ok"),
        "failed": sum(1 for r in rows if r["status"] != "ok"),
        "preset": preset,
        "ph": ph,
        "rows": rows,
        "zip_size_bytes": len(zip_bytes),
    }


def get_batch_zip(batch_id: str) -> Tuple[bytes, str]:
    if batch_id not in _BATCH_STORE:
        raise KeyError(f"Batch '{batch_id}' not found or expired")
    return _BATCH_STORE[batch_id]
