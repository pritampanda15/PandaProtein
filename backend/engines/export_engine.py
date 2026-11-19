"""
Pre-Flight Readiness Audit & Multi-Target Export Engine:
- Independent Readiness Audits: Molecular Dynamics, Molecular Docking, Crystallographic Quality
- Clean PDB & mmCIF export
- AutoDock Vina / GNINA PDBQT receptor export + auto-computed vina_config.txt
- PQR (APBS electrostatics) export
- Production GROMACS MD Bundle (.zip) with live `gmx pdb2gmx` support + MDP suite + shell runner
- Publication-grade Markdown & JSON reports
"""

from __future__ import annotations

import io
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import gemmi
import numpy as np

from .chemistry_data import (
    ALL_PROTEIN_RESNAMES,
    BACKBONE_CHARGES,
    CRYSTALLIZATION_ARTIFACTS,
    METAL_IONS,
    VDW_RADII,
    WATER_RESNAMES,
)
from .geometry_utils import canonical_resname
from .structure_io import structure_to_gro, structure_to_mmcif, structure_to_pdb


def check_external_tools() -> Dict[str, Any]:
    gmx_path = shutil.which("gmx") or ("/opt/homebrew/bin/gmx" if os.path.exists("/opt/homebrew/bin/gmx") else None)
    obabel_path = shutil.which("obabel") or ("/opt/homebrew/bin/obabel" if os.path.exists("/opt/homebrew/bin/obabel") else None)
    try:
        import propka  # noqa: F401
        has_propka = True
    except Exception:
        has_propka = False
    try:
        import rdkit  # noqa: F401
        has_rdkit = True
    except Exception:
        has_rdkit = False
    return {
        "gemmi": True,
        "propka": has_propka,
        "rdkit": has_rdkit,
        "gromacs": bool(gmx_path),
        "gromacs_path": gmx_path,
        "openbabel": bool(obabel_path),
    }


def compute_readiness_audit(
    fixer_diag: Dict[str, Any],
    quality: Dict[str, Any],
    propka: Dict[str, Any],
    lig_info: Dict[str, Any],
) -> Dict[str, Any]:
    """Compute independent MD, Docking, and Crystallographic readiness checklists & 0-100 scores."""

    n_nonstd = len(fixer_diag.get("nonstandard_residues", []))
    n_missing_atoms = len(fixer_diag.get("missing_atoms", []))
    n_missing_loops = len(fixer_diag.get("missing_loops", []))
    n_altlocs = len(fixer_diag.get("altlocs", []))
    n_buffers = sum(1 for h in fixer_diag.get("heterogens", []) if h["category"] == "buffer")
    has_h = bool(fixer_diag.get("has_hydrogens", False))
    unassigned_ss = sum(
        1 for d in fixer_diag.get("disulfides", [])
        if d["bonded"] and (d["resname1"] != "CYX" or d["resname2"] != "CYX")
    )
    clashscore = float(quality["clashes"]["clashscore"])
    rama_fav = float(quality["ramachandran"]["favored_pct"])
    rama_out = int(quality["ramachandran"]["counts"]["outlier"])
    rot_out = int(quality["rotamers"]["counts"]["outlier"])
    cis_nonpro = sum(1 for c in quality["geometry"]["cis_peptides"] if c["severity"] == "outlier")
    chiral_out = len(quality["geometry"]["chiral_outliers"])
    bond_rmsz = float(quality["geometry"]["bond_rmsz"])

    def _item(label: str, status: str, detail: str, penalty: float) -> Dict[str, Any]:
        return {"label": label, "status": status, "detail": detail, "penalty": penalty}

    # 1. MD Simulation Readiness Checklist
    md_checks = [
        _item(
            "Complete Heavy Atoms",
            "pass" if n_missing_atoms == 0 else "fail",
            "All side-chain & terminal atoms present" if n_missing_atoms == 0 else f"{n_missing_atoms} residue(s) missing heavy atoms",
            0 if n_missing_atoms == 0 else min(22, n_missing_atoms * 7),
        ),
        _item(
            "Continuous Backbone / Loops",
            "pass" if n_missing_loops == 0 else "warn",
            "No internal chain breaks detected" if n_missing_loops == 0 else f"{n_missing_loops} missing loop/chain break(s) — rebuild or cap termini",
            0 if n_missing_loops == 0 else min(16, n_missing_loops * 10),
        ),
        _item(
            "Canonical Residue Types",
            "pass" if n_nonstd == 0 else "fail",
            "All residues recognized by AMBER/CHARMM/GROMACS" if n_nonstd == 0 else f"{n_nonstd} non-standard residue(s) need conversion",
            0 if n_nonstd == 0 else min(15, n_nonstd * 8),
        ),
        _item(
            "Single Conformation (No AltLocs)",
            "pass" if n_altlocs == 0 else "fail",
            "All occupancies resolved to 1.0" if n_altlocs == 0 else f"{n_altlocs} residue(s) have multiple AltLoc conformers",
            0 if n_altlocs == 0 else min(12, n_altlocs * 6),
        ),
        _item(
            "Crystallization Buffers Removed",
            "pass" if n_buffers == 0 else "warn",
            "No buffer artifacts (SO4/EDO/PEG)" if n_buffers == 0 else f"{n_buffers} crystallization buffer artifact(s) present",
            0 if n_buffers == 0 else min(10, n_buffers * 5),
        ),
        _item(
            "Protonation & Hydrogens (PROPKA)",
            "pass" if has_h else "fail",
            f"Hydrogens built (net charge {propka['net_charge_discrete']:+.1f} at pH {propka['ph']})" if has_h else "Hydrogens not yet added — run PROPKA Protonation",
            0 if has_h else 16,
        ),
        _item(
            "Disulfide Bridges (CYX)",
            "pass" if unassigned_ss == 0 else "warn",
            "All S–S bridges assigned" if unassigned_ss == 0 else f"{unassigned_ss} S–S pair(s) within 2.5 Å not yet marked CYX",
            0 if unassigned_ss == 0 else 8,
        ),
        _item(
            "Steric Clashscore",
            "pass" if clashscore <= 5.0 else ("warn" if clashscore <= 15.0 else "fail"),
            f"Clashscore = {clashscore} ({quality['clashes']['total_clashes']} clashes)",
            0 if clashscore <= 5.0 else min(18, (clashscore - 5.0) * 0.9),
        ),
    ]
    md_score = round(max(0.0, 100.0 - sum(c["penalty"] for c in md_checks)), 1)

    # 2. Molecular Docking Readiness Checklist
    has_box = len(lig_info.get("docking_boxes", [])) > 0
    dock_checks = [
        _item(
            "Clean Receptor Heavy Atoms",
            "pass" if (n_missing_atoms == 0 and n_altlocs == 0) else "fail",
            "Complete, single-conformer receptor" if (n_missing_atoms == 0 and n_altlocs == 0) else f"{n_missing_atoms} missing atom(s), {n_altlocs} AltLoc(s)",
            0 if (n_missing_atoms == 0 and n_altlocs == 0) else 24,
        ),
        _item(
            "Polar Hydrogens & Tautomers",
            "pass" if has_h else "fail",
            "Polar hydrogens & HIS tautomers assigned for H-bond scoring" if has_h else "Add polar hydrogens via PROPKA before PDBQT export",
            0 if has_h else 22,
        ),
        _item(
            "Crystallization Buffers Cleared",
            "pass" if n_buffers == 0 else "warn",
            "Binding site free of buffer artifacts" if n_buffers == 0 else f"{n_buffers} crystallization artifact(s) may occlude pocket",
            0 if n_buffers == 0 else 12,
        ),
        _item(
            "Binding Pocket Grid Box",
            "pass" if has_box else "warn",
            f"Auto-computed box ready ({lig_info['docking_boxes'][0]['label']})" if has_box else "No pocket detected yet",
            0 if has_box else 10,
        ),
        _item(
            "Side-Chain Rotamers & Clashes",
            "pass" if (rot_out == 0 and clashscore <= 8.0) else "warn",
            f"{rot_out} rotamer outlier(s), clashscore {clashscore}",
            0 if (rot_out == 0 and clashscore <= 8.0) else min(20, rot_out * 4 + max(0, clashscore - 5)),
        ),
    ]
    dock_score = round(max(0.0, 100.0 - sum(c["penalty"] for c in dock_checks)), 1)

    # 3. Crystallographic Quality (COOT / MolProbity)
    cryst_checks = [
        _item(
            "Ramachandran Favored",
            "pass" if (rama_fav >= 95.0 and rama_out == 0) else ("warn" if rama_out <= 1 else "fail"),
            f"{rama_fav}% favored, {rama_out} outlier(s)",
            0 if rama_out == 0 else min(25, rama_out * 8),
        ),
        _item(
            "Rotamer Outliers",
            "pass" if rot_out == 0 else ("warn" if rot_out <= 2 else "fail"),
            f"{quality['rotamers']['favored_pct']}% favored, {rot_out} outlier(s)",
            0 if rot_out == 0 else min(20, rot_out * 5),
        ),
        _item(
            "Cis / Twisted Peptides",
            "pass" if cis_nonpro == 0 else "fail",
            "No non-Pro cis or twisted peptides" if cis_nonpro == 0 else f"{cis_nonpro} cis-non-Pro / twisted peptide(s)",
            0 if cis_nonpro == 0 else min(20, cis_nonpro * 10),
        ),
        _item(
            "Stereochemistry & Chirality",
            "pass" if chiral_out == 0 else "fail",
            "All Cα centers valid L-amino acids" if chiral_out == 0 else f"{chiral_out} chirality outlier(s)",
            0 if chiral_out == 0 else min(20, chiral_out * 12),
        ),
        _item(
            "Engh & Huber Bond RMSZ",
            "pass" if bond_rmsz <= 1.2 else ("warn" if bond_rmsz <= 2.0 else "fail"),
            f"Bond RMSZ = {bond_rmsz}, Angle RMSZ = {quality['geometry']['angle_rmsz']}",
            0 if bond_rmsz <= 1.2 else min(15, (bond_rmsz - 1.2) * 8),
        ),
    ]
    cryst_score = quality["quality_score"]

    overall = round((md_score + dock_score + cryst_score) / 3.0, 1)
    return {
        "overall_score": overall,
        "md": {"score": md_score, "ready": md_score >= 85.0, "checks": md_checks},
        "docking": {"score": dock_score, "ready": dock_score >= 85.0, "checks": dock_checks},
        "crystallographic": {"score": cryst_score, "ready": cryst_score >= 85.0, "checks": cryst_checks},
    }


def _estimate_atom_charge(resname: str, atom_name: str, element: str) -> float:
    if atom_name in BACKBONE_CHARGES:
        return BACKBONE_CHARGES[atom_name]
    r = resname.upper()
    if element == "H":
        return 0.25 if atom_name.startswith(("HN", "H1", "H2", "H3", "HD1", "HE2", "HG", "HH", "HZ", "HE")) else 0.08
    if r == "ASP" and atom_name in ("OD1", "OD2"):
        return -0.78
    if r == "GLU" and atom_name in ("OE1", "OE2"):
        return -0.78
    if r == "LYS" and atom_name == "NZ":
        return -0.35
    if r == "ARG" and atom_name in ("NH1", "NH2"):
        return -0.60
    if r == "ARG" and atom_name == "CZ":
        return 0.80
    if r == "ZN":
        return 2.00
    if r in ("MG", "CA", "FE2", "MN"):
        return 2.00
    if element == "O":
        return -0.45
    if element == "N":
        return -0.35
    if element == "S":
        return -0.20
    return 0.05


def _autodock_atom_type(resname: str, atom_name: str, element: str) -> str:
    el = element.upper()
    if el == "H":
        return "HD"
    if el == "O":
        return "OA"
    if el == "N":
        return "NA" if atom_name in ("ND1", "NE2") else "N"
    if el == "S":
        return "SA"
    if el == "C":
        if canonical_resname(resname) in ("PHE", "TYR", "TRP", "HIS") and atom_name not in ("CA", "C", "CB"):
            return "A"
        return "C"
    if el == "ZN":
        return "Zn"
    if el == "MG":
        return "Mg"
    if el == "CA":
        return "Ca"
    if el == "FE":
        return "Fe"
    if el == "MN":
        return "Mn"
    return el[:2]


def generate_receptor_pdbqt(st: gemmi.Structure, exclude_ligands: bool = True) -> str:
    """Export AutoDock Vina / GNINA receptor PDBQT with polar hydrogens, partial charges, and AD4 atom types."""
    lines = [
        "REMARK  Generated by PDBPrep Studio — AutoDock Vina / GNINA Receptor PDBQT",
        f"REMARK  Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
    ]
    serial = 1
    if len(st) > 0:
        for chain in st[0]:
            for res in chain:
                rn = res.name.strip().upper()
                is_prot = rn in ALL_PROTEIN_RESNAMES
                is_metal = rn in METAL_IONS
                if exclude_ligands and not (is_prot or is_metal):
                    continue
                if rn in WATER_RESNAMES or rn in CRYSTALLIZATION_ARTIFACTS:
                    continue

                for atom in res:
                    el = atom.element.name.upper()
                    an = atom.name.strip()
                    # AutoDock PDBQT keeps polar hydrogens (HD), drops aliphatic C-H
                    if el == "H" and an.startswith(("HA", "HB", "HG2", "HG3", "HD2", "HD3", "HE3")):
                        continue
                    q = _estimate_atom_charge(rn, an, el)
                    ad_type = _autodock_atom_type(rn, an, el)
                    rec = "ATOM  " if is_prot else "HETATM"
                    an_fmt = f" {an:<3}" if len(an) < 4 else f"{an:<4}"
                    lines.append(
                        f"{rec}{serial:5d} {an_fmt} {rn:>3} {chain.name[:1]:1s}{res.seqid.num:4d}    "
                        f"{atom.pos.x:8.3f}{atom.pos.y:8.3f}{atom.pos.z:8.3f}"
                        f"{atom.occ:6.2f}{atom.b_iso:6.2f}    {q:+6.3f} {ad_type:<2s}"
                    )
                    serial += 1
            lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def generate_pqr(st: gemmi.Structure) -> str:
    """Export APBS-compatible PQR file with partial charges and atomic van der Waals radii."""
    lines = [
        "REMARK  Generated by PDBPrep Studio — APBS PQR Format",
    ]
    serial = 1
    if len(st) > 0:
        for chain in st[0]:
            for res in chain:
                rn = res.name.strip().upper()
                rec = "ATOM  " if rn in ALL_PROTEIN_RESNAMES else "HETATM"
                for atom in res:
                    el = atom.element.name.upper()
                    an = atom.name.strip()
                    q = _estimate_atom_charge(rn, an, el)
                    rad = VDW_RADII.get(el, 1.70)
                    lines.append(
                        f"{rec}{serial:5d} {an:<4s} {rn:>3s} {chain.name[:1]:1s}{res.seqid.num:4d}    "
                        f"{atom.pos.x:8.3f}{atom.pos.y:8.3f}{atom.pos.z:8.3f} {q:+7.4f} {rad:6.3f}"
                    )
                    serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def generate_vina_config(docking_boxes: List[Dict[str, Any]], selected_ligand_id: Optional[str] = None, padding: float = 6.0) -> str:
    box = None
    if selected_ligand_id:
        box = next((b for b in docking_boxes if b["ligand_id"] == selected_ligand_id), None)
    if box is None and docking_boxes:
        box = docking_boxes[0]
    cx, cy, cz = box["center"] if box else (0.0, 0.0, 0.0)
    sx, sy, sz = box["size"] if box else (22.0, 22.0, 22.0)
    return (
        "# AutoDock Vina / GNINA Configuration — Generated by PDBPrep Studio\n"
        "receptor = receptor.pdbqt\n"
        "ligand = ligand.pdbqt\n"
        "out = docked_poses.pdbqt\n\n"
        f"center_x = {cx:.2f}\n"
        f"center_y = {cy:.2f}\n"
        f"center_z = {cz:.2f}\n\n"
        f"size_x = {max(14.0, sx + padding - 6.0):.1f}\n"
        f"size_y = {max(14.0, sy + padding - 6.0):.1f}\n"
        f"size_z = {max(14.0, sz + padding - 6.0):.1f}\n\n"
        "exhaustiveness = 16\n"
        "num_modes = 20\n"
        "energy_range = 4\n"
    )


def _adapt_protein_for_forcefield(st: gemmi.Structure, force_field: str) -> gemmi.Structure:
    """
    Clone protein-only chains and adapt protonation residue names cleanly to the target
    GROMACS force field's .rtp nomenclature (CHARMM36m/27 vs AMBER vs OPLS/GROMOS).
    """
    ff_low = force_field.lower()
    is_charmm = "charmm" in ff_low
    is_amber = "amber" in ff_low

    charmm_map = {
        "HID": "HSD", "HIE": "HSE", "HIP": "HSP",
        "CYX": "CYS", "CYM": "CYS",
        "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "TYM": "TYR", "ARN": "ARG",
    }
    generic_map = {
        "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
        "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
        "CYX": "CYS", "CYM": "CYS",
        "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "TYM": "TYR", "ARN": "ARG",
    }
    amber_map = {
        "HSD": "HID", "HSE": "HIE", "HSP": "HIP",
        "CYM": "CYS", "TYM": "TYR", "ARN": "ARG",
    }

    prot_only = st.clone()
    if len(prot_only) > 0:
        for ch in prot_only[0]:
            for ri in range(len(ch) - 1, -1, -1):
                rn = ch[ri].name.strip().upper()
                if rn not in ALL_PROTEIN_RESNAMES:
                    del ch[ri]
                    continue
                if is_charmm and rn in charmm_map:
                    ch[ri].name = charmm_map[rn]
                elif is_amber and rn in amber_map:
                    ch[ri].name = amber_map[rn]
                elif not is_charmm and not is_amber and rn in generic_map:
                    ch[ri].name = generic_map[rn]
    return prot_only


def _mdp_templates(salt_conc: float = 0.15, force_field: str = "amber99sb-ildn") -> Dict[str, str]:
    is_charmm = "charmm" in force_field.lower()
    if is_charmm:
        nonbonded = (
            "cutoff-scheme = Verlet\n"
            "coulombtype = PME\nrcoulomb    = 1.2\n"
            "vdwtype     = cutoff\nvdw-modifier = force-switch\n"
            "rlist       = 1.2\nrvdw        = 1.2\nrvdw-switch = 1.0\n"
            "DispCorr    = no\n"
        )
        ff_note = f"CHARMM36m/CHARMM recommended Force-Switch 1.0-1.2 nm ({force_field})"
    else:
        nonbonded = (
            "cutoff-scheme = Verlet\n"
            "coulombtype = PME\nrcoulomb    = 1.0\n"
            "rvdw        = 1.0\nDispCorr    = EnerPres\n"
        )
        ff_note = f"{force_field} standard PME 1.0 nm cutoffs"

    return {
        "ions.mdp": (
            f"; Ions placement MDP (PDBPrep Studio — {ff_note})\n"
            "integrator  = steep\nemtol       = 1000.0\nemstep      = 0.01\nnsteps      = 5000\n"
            "nstlist     = 1\ncutoff-scheme = Verlet\ncoulombtype = cutoff\nrcoulomb    = 1.2\nrvdw        = 1.2\npbc         = xyz\n"
        ),
        "em.mdp": (
            f"; Energy Minimization MDP (PDBPrep Studio — {ff_note})\n"
            "integrator  = steep\nemtol       = 500.0\nemstep      = 0.01\nnsteps      = 50000\n"
            f"nstlist     = 10\n{nonbonded}pbc         = xyz\n"
        ),
        "nvt.mdp": (
            f"; NVT Equilibration 100 ps at 300 K (PDBPrep Studio — {ff_note})\n"
            "define      = -DPOSRES\nintegrator  = md\nnsteps      = 50000\ndt          = 0.002\n"
            "nstxout-compressed = 5000\ncontinuation = no\nconstraint_algorithm = lincs\nconstraints = h-bonds\n"
            f"{nonbonded}"
            "tcoupl      = V-rescale\ntc-grps     = Protein Non-Protein\ntau_t       = 0.1 0.1\nref_t       = 300 300\n"
            "pcoupl      = no\npbc         = xyz\ngen_vel     = yes\ngen_temp    = 300\ngen_seed    = -1\n"
        ),
        "npt.mdp": (
            f"; NPT Equilibration 100 ps at 300 K, 1 bar (PDBPrep Studio — {ff_note})\n"
            "define      = -DPOSRES\nintegrator  = md\nnsteps      = 50000\ndt          = 0.002\n"
            "nstxout-compressed = 5000\ncontinuation = yes\nconstraint_algorithm = lincs\nconstraints = h-bonds\n"
            f"{nonbonded}"
            "tcoupl      = V-rescale\ntc-grps     = Protein Non-Protein\ntau_t       = 0.1 0.1\nref_t       = 300 300\n"
            "pcoupl      = C-rescale\npcoupltype  = isotropic\ntau_p       = 2.0\nref_p       = 1.0\ncompressibility = 4.5e-5\n"
            "refcoord_scaling = com\npbc         = xyz\ngen_vel     = no\n"
        ),
        "md.mdp": (
            f"; Production MD 100 ns (PDBPrep Studio — {ff_note})\n"
            "integrator  = md\nnsteps      = 50000000\ndt          = 0.002\n"
            "nstxout-compressed = 5000\nnstenergy   = 5000\nnstlog      = 5000\n"
            "continuation = yes\nconstraint_algorithm = lincs\nconstraints = h-bonds\n"
            f"{nonbonded}"
            "tcoupl      = V-rescale\ntc-grps     = Protein Non-Protein\ntau_t       = 0.1 0.1\nref_t       = 300 300\n"
            "pcoupl      = Parrinello-Rahman\npcoupltype  = isotropic\ntau_p       = 2.0\nref_p       = 1.0\ncompressibility = 4.5e-5\n"
            "pbc         = xyz\ngen_vel     = no\n"
        ),
    }


def generate_markdown_report(
    name: str,
    metadata: Dict[str, Any],
    readiness: Dict[str, Any],
    quality: Dict[str, Any],
    propka: Dict[str, Any],
    lig_info: Dict[str, Any],
    history: List[str],
) -> str:
    lines = [
        f"# PDBPrep Studio — Deep Protein Preparation Report",
        f"**Structure:** `{name}`  ",
        f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  ",
        f"**Title:** {metadata.get('title', 'N/A')}  ",
        f"**Method / Resolution:** {metadata.get('method', 'N/A')} / {metadata.get('resolution') or 'N/A'} Å  ",
        "",
        "## 1. Executive Readiness Scores",
        f"| Target Workflow | Score | Status |",
        f"|---|---|---|",
        f"| **Molecular Dynamics (GROMACS / OpenMM)** | **{readiness['md']['score']}/100** | {'READY' if readiness['md']['ready'] else 'NEEDS ATTENTION'} |",
        f"| **Molecular Docking (Vina / GNINA)** | **{readiness['docking']['score']}/100** | {'READY' if readiness['docking']['ready'] else 'NEEDS ATTENTION'} |",
        f"| **Crystallographic Quality (COOT / MolProbity)** | **{readiness['crystallographic']['score']}/100** | {'EXCELLENT' if readiness['crystallographic']['ready'] else 'REVIEW'} |",
        "",
        "## 2. Crystallographic & Geometry Validation (COOT / MolProbity)",
        f"- **Ramachandran Favored / Allowed / Outliers:** {quality['ramachandran']['favored_pct']}% / {quality['ramachandran']['allowed_pct']}% / {quality['ramachandran']['outlier_pct']}% ({quality['ramachandran']['counts']['outlier']} outliers)",
        f"- **Rotamer Favored / Outliers:** {quality['rotamers']['favored_pct']}% / {quality['rotamers']['outlier_pct']}% ({quality['rotamers']['counts']['outlier']} outliers)",
        f"- **All-Atom Clashscore:** {quality['clashes']['clashscore']} ({quality['clashes']['total_clashes']} clashes across {quality['clashes']['total_atoms']} atoms)",
        f"- **Engh & Huber RMSZ:** Bonds = {quality['geometry']['bond_rmsz']}, Angles = {quality['geometry']['angle_rmsz']}",
        f"- **Cis / Twisted Peptides:** {len(quality['geometry']['cis_peptides'])}",
        f"- **Chirality Outliers:** {len(quality['geometry']['chiral_outliers'])}",
        "",
        "## 3. PROPKA 3 Protonation & Electrostatics",
        f"- **Engine:** {propka['engine']}",
        f"- **Target pH:** {propka['ph']} | **Isoelectric Point (pI):** {propka['isoelectric_point']}",
        f"- **Net Protein Charge:** {propka['net_charge_discrete']:+.1f} (discrete) / {propka['net_charge_continuous']:+.2f} (continuous)",
        f"- **Titratable Residues:** {propka['titratable_count']} ({propka['shifted_count']} with |ΔpKa| ≥ 1.5, {propka['metal_coordinated_count']} metal-coordinated)",
        "",
        "## 4. Ligands, Cofactors & Metal Sites",
        f"- **Detected Ligands / Cofactors:** {len(lig_info['ligands'])}",
        f"- **Detected Metal Ions:** {len(lig_info['metals'])}",
        f"- **Remaining Crystallization Buffers:** {len(lig_info['buffers'])}",
        "",
        "## 5. Preparation Audit Trail",
    ]
    for step in history:
        lines.append(f"- {step}")
    return "\n".join(lines) + "\n"


def build_gromacs_zip_bundle(
    st: gemmi.Structure,
    name: str,
    opts: Dict[str, Any],
    metadata: Dict[str, Any],
    readiness: Dict[str, Any],
    quality: Dict[str, Any],
    propka: Dict[str, Any],
    lig_info: Dict[str, Any],
    history: List[str],
) -> bytes:
    """Build a complete production ZIP bundle for GROMACS MD + AutoDock Vina Docking + Reports."""
    ff = opts.get("force_field", "amber99sb-ildn")
    water = opts.get("water_model", "tip3p")
    box_type = opts.get("box_type", "dodecahedron")
    box_dist = float(opts.get("box_distance_nm", 1.0))
    salt = float(opts.get("salt_conc_m", 0.15))

    # Create force-field-adapted protein-only PDB for gmx pdb2gmx
    prot_only = _adapt_protein_for_forcefield(st, ff)

    full_pdb = structure_to_pdb(st)
    prot_pdb = structure_to_pdb(prot_only)
    pdbqt_str = generate_receptor_pdbqt(st)
    pqr_str = generate_pqr(st)
    vina_cfg = generate_vina_config(
        lig_info.get("docking_boxes", []),
        opts.get("selected_ligand_id"),
        float(opts.get("docking_padding", 6.0)),
    )
    md_report = generate_markdown_report(name, metadata, readiness, quality, propka, lig_info, history)

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("structures/complex_prepared.pdb", full_pdb)
        zf.writestr("structures/protein_only.pdb", prot_pdb)
        zf.writestr("structures/complex_prepared.gro", structure_to_gro(st))
        zf.writestr("structures/complex_prepared.cif", structure_to_mmcif(st))
        zf.writestr("structures/receptor.pqr", pqr_str)
        zf.writestr("docking/receptor.pdbqt", pdbqt_str)
        zf.writestr("docking/vina_config.txt", vina_cfg)

        for fname, content in _mdp_templates(salt, ff).items():
            zf.writestr(f"gromacs/mdp/{fname}", content)

        # Try executing real `gmx pdb2gmx` if GROMACS is installed
        tools = check_external_tools()
        gmx_log = "GROMACS not detected on host; ready-to-run script & MDP templates included."
        if tools["gromacs"] and tools["gromacs_path"]:
            with tempfile.TemporaryDirectory() as td:
                in_pdb = os.path.join(td, "protein_only.pdb")
                out_gro = os.path.join(td, "processed.gro")
                out_top = os.path.join(td, "topol.top")
                with open(in_pdb, "w") as f:
                    f.write(prot_pdb)
                cmd = [
                    tools["gromacs_path"], "pdb2gmx",
                    "-f", in_pdb, "-o", out_gro, "-p", out_top,
                    "-ff", ff, "-water", water, "-ignh", "-missing",
                ]
                try:
                    proc = subprocess.run(cmd, cwd=td, capture_output=True, text=True, timeout=25)
                    gmx_log = f"Command: {' '.join(cmd)}\nReturn code: {proc.returncode}\n\nSTDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}"
                    for produced in ("processed.gro", "topol.top", "posre.itp"):
                        p = os.path.join(td, produced)
                        if os.path.exists(p):
                            with open(p, "r") as pf:
                                zf.writestr(f"gromacs/{produced}", pf.read())
                except Exception as exc:
                    gmx_log = f"gmx pdb2gmx execution skipped: {exc}"

        zf.writestr("gromacs/pdb2gmx.log", gmx_log)

        run_sh = f"""#!/usr/bin/env bash
# Production GROMACS MD Pipeline — Generated by PDBPrep Studio
set -euo pipefail

gmx pdb2gmx -f ../structures/protein_only.pdb -o processed.gro -p topol.top -ff {ff} -water {water} -ignh
gmx editconf -f processed.gro -o boxed.gro -c -d {box_dist:.2f} -bt {box_type}
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
gmx grompp -f mdp/ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2
echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc {salt:.2f}
gmx grompp -f mdp/em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
gmx grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
gmx grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt
gmx grompp -f mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -v -deffnm md
"""
        zf.writestr("gromacs/run_gromacs_pipeline.sh", run_sh)
        zf.writestr("reports/preparation_report.md", md_report)
        zf.writestr(
            "reports/preparation_report.json",
            json.dumps(
                {
                    "name": name,
                    "metadata": metadata,
                    "readiness": readiness,
                    "propka_summary": {k: v for k, v in propka.items() if k != "residues"},
                    "history": history,
                },
                indent=2,
            ),
        )

    return buf.getvalue()


def estimate_solvation_system(
    st: gemmi.Structure,
    box_type: str = "dodecahedron",
    box_distance_nm: float = 1.0,
    salt_conc_m: float = 0.15,
    net_charge: float = 0.0,
) -> Dict[str, Any]:
    """Estimate periodic simulation box volume, TIP3P/SPC water molecules, and Na+/Cl- ions."""
    pts = []
    if len(st) > 0:
        for ch in st[0]:
            for res in ch:
                if res.name.strip().upper() in ALL_PROTEIN_RESNAMES:
                    for a in res:
                        pts.append([a.pos.x / 10.0, a.pos.y / 10.0, a.pos.z / 10.0])  # nm

    if not pts:
        return {
            "protein_span_nm": [0.0, 0.0, 0.0],
            "box_diameter_nm": 0.0,
            "box_volume_nm3": 0.0,
            "estimated_waters": 0,
            "na_ions": 0,
            "cl_ions": 0,
            "total_atoms_est": 0,
        }

    arr = np.array(pts, dtype=np.float64)
    span = np.max(arr, axis=0) - np.min(arr, axis=0)
    diam = float(np.linalg.norm(span))
    d_box = max(float(np.max(span)), diam * 0.85) + 2.0 * box_distance_nm

    if box_type == "dodecahedron":
        vol = 0.7071 * (d_box ** 3)
    elif box_type == "octahedron":
        vol = 0.7698 * (d_box ** 3)
    else:  # cubic
        vol = d_box ** 3

    # Protein excluded volume (~0.015 nm³ per heavy atom)
    prot_vol = len(pts) * 0.014
    solvent_vol = max(1.0, vol - prot_vol)
    waters = int(round(solvent_vol * 33.36))

    # Salt pairs for salt_conc_m (mol/L): N = C * V(nm³) * 0.6022
    salt_pairs = int(round(salt_conc_m * solvent_vol * 0.6022))
    q_int = int(round(net_charge))
    na_ions = salt_pairs + max(0, -q_int)
    cl_ions = salt_pairs + max(0, q_int)
    waters_final = max(0, waters - na_ions - cl_ions)
    total_atoms = len(pts) + waters_final * 3 + na_ions + cl_ions

    return {
        "protein_span_nm": [round(float(x), 2) for x in span],
        "box_diameter_nm": round(d_box, 2),
        "box_volume_nm3": round(vol, 1),
        "estimated_waters": waters_final,
        "na_ions": na_ions,
        "cl_ions": cl_ions,
        "net_charge_neutralized": q_int,
        "total_atoms_est": total_atoms,
    }


def run_gromacs_preflight_test(
    st: gemmi.Structure,
    force_field: str = "amber99sb-ildn",
    water_model: str = "tip3p",
) -> Dict[str, Any]:
    """Execute a live `gmx pdb2gmx` pre-flight topology check and parse exact GROMACS diagnostics."""
    tools = check_external_tools()
    if not tools["gromacs"] or not tools["gromacs_path"]:
        return {
            "available": False,
            "passed": False,
            "message": "GROMACS (`gmx`) binary not found on host PATH",
            "force_field": force_field,
            "water_model": water_model,
            "stdout": "",
            "stderr": "",
        }

    # Build clean force-field-adapted protein-only PDB for pdb2gmx topology verification
    prot_only = _adapt_protein_for_forcefield(st, force_field)
    prot_pdb = structure_to_pdb(prot_only)

    with tempfile.TemporaryDirectory() as td:
        in_pdb = os.path.join(td, "protein_only.pdb")
        out_gro = os.path.join(td, "processed.gro")
        out_top = os.path.join(td, "topol.top")
        with open(in_pdb, "w") as f:
            f.write(prot_pdb)

        cmd = [
            tools["gromacs_path"],
            "pdb2gmx",
            "-f", in_pdb,
            "-o", out_gro,
            "-p", out_top,
            "-ff", force_field,
            "-water", water_model,
            "-ignh",
        ]
        try:
            proc = subprocess.run(cmd, cwd=td, capture_output=True, text=True, timeout=30)
            combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
            mass_m = re.search(r"Total mass\s+([-\d.]+)\s+a\.m\.u\.", combined)
            charge_m = re.search(r"Total charge\s+([-\d.]+)\s+e", combined)
            fatal_m = re.search(r"Fatal error:\s*(.+?)(?:\n\n|\Z)", combined, re.DOTALL)
            gro_atoms = None
            if os.path.exists(out_gro):
                with open(out_gro, "r") as gf:
                    glines = gf.readlines()
                    if len(glines) >= 2 and glines[1].strip().isdigit():
                        gro_atoms = int(glines[1].strip())

            passed = proc.returncode == 0 and os.path.exists(out_top)
            return {
                "available": True,
                "passed": passed,
                "return_code": proc.returncode,
                "force_field": force_field,
                "water_model": water_model,
                "total_mass_amu": float(mass_m.group(1)) if mass_m else None,
                "total_charge_e": float(charge_m.group(1)) if charge_m else None,
                "topology_atoms": gro_atoms,
                "fatal_error": fatal_m.group(1).strip() if fatal_m else None,
                "message": (
                    f"GROMACS pdb2gmx passed ({force_field} / {water_model}) — {gro_atoms} topology atoms, charge {charge_m.group(1) if charge_m else '?'} e"
                    if passed
                    else (fatal_m.group(1).strip().splitlines()[0] if fatal_m else f"pdb2gmx exited with code {proc.returncode}")
                ),
                "log": combined[-3500:],
            }
        except Exception as exc:
            return {
                "available": True,
                "passed": False,
                "force_field": force_field,
                "water_model": water_model,
                "message": f"GROMACS test failed: {exc}",
                "log": str(exc),
            }
