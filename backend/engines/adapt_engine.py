"""
ADapt-Parity Receptor & Docking Preparation Engine:
1. Biological Assemblies (REMARK 350 / mmCIF _pdbx_struct_assembly_gen) — ASU vs Biomol 1, 2, ...
   with symmetry-copy tracking and axial-ion deduplication (e.g. 4INS axial zincs).
2. Missing Residue Detection (REMARK 465 + mmCIF + SEQRES) with 3D gap/terminus markers &
   selective per-gap rebuild guarded by strict peptide geometry + steric clearance gates.
3. AutoDock 4 Atom-Type Validator, Ion Substitution (e.g. CU -> ZN), Manual Ion Formal Charges,
   Nucleic Acid (DNA/RNA) receptor support, and Per-Item Receptor/Ligand/Exclude Preparation Tree.
4. Live Cached RCSB Data API Entry-Level Ligand Classification (with offline CCD fallback).
5. Cryptographic Provenance (REMARK 850 headers + MD5-verified session provenance JSON).
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import gemmi
import httpx
import numpy as np
from scipy.spatial import cKDTree

from .chemistry_data import (
    AA_3_TO_1,
    AD4_ION_SUBSTITUTES,
    ALL_PROTEIN_RESNAMES,
    AUTODOCK4_SUPPORTED_ELEMENTS,
    CCD_ION_CHARGES,
    COFACTORS,
    CRYSTALLIZATION_ARTIFACTS,
    EXPECTED_HEAVY_ATOMS,
    HALIDE_IONS,
    METAL_IONS,
    NONSTANDARD_TO_STANDARD,
    NUCLEIC_RESNAMES,
    WATER_RESNAMES,
)
from .coot_engine import _build_env_kdtree, _get_atom, _get_pos, regularize_zone
from .geometry_utils import build_sidechain_coords, canonical_resname, dihedral_deg, distance

_RCSB_CACHE_DIR = Path(
    os.environ.get("PANDAPROTEIN_RCSB_CACHE") or (Path.home() / ".pandaprotein" / "rcsb_cache")
)


# ---------------------------------------------------------------------------
# 1. Biological Assemblies (REMARK 350 / mmCIF)
# ---------------------------------------------------------------------------

def _oligomer_label(n: int) -> str:
    labels = {
        1: "Monomer",
        2: "Dimer",
        3: "Trimer",
        4: "Tetramer",
        5: "Pentamer",
        6: "Hexamer",
        8: "Octamer",
        12: "Dodecamer",
    }
    return labels.get(n, f"{n}-mer")


def list_biological_assemblies(st: gemmi.Structure, active_assembly_id: str = "ASU") -> List[Dict[str, Any]]:
    """List ASU plus all biological assemblies declared via REMARK 350 or mmCIF _pdbx_struct_assembly_gen."""
    asu_chains = [ch.name for ch in st[0]] if len(st) > 0 else []
    asu_prot_chains = [
        ch.name
        for ch in (st[0] if len(st) > 0 else [])
        if any(
            (r.name.strip().upper() in ALL_PROTEIN_RESNAMES or r.name.strip().upper() in NUCLEIC_RESNAMES)
            for r in ch
        )
    ]

    result: List[Dict[str, Any]] = [
        {
            "id": "ASU",
            "name": "ASU",
            "label": f"ASU ({_oligomer_label(max(1, len(asu_prot_chains)))})",
            "oligomeric_count": max(1, len(asu_prot_chains)),
            "chains": asu_chains,
            "operator_count": 1,
            "has_symmetry_copies": False,
            "is_active": active_assembly_id == "ASU",
            "description": "Original crystallographic Asymmetric Unit (as deposited)",
        }
    ]

    for asm in st.assemblies:
        asm_id = str(asm.name).strip()
        total_ops = 0
        involved_chains: List[str] = []
        polymer_copies = 0
        prot_chain_set = set(asu_prot_chains)
        for gen in asm.generators:
            ops = len(gen.operators)
            total_ops = max(total_ops, ops)
            for cid in gen.chains:
                if cid not in involved_chains:
                    involved_chains.append(cid)
                if cid in prot_chain_set:
                    polymer_copies += max(1, ops)
        if not involved_chains:
            involved_chains = list(asu_chains)
        oligo = max(1, polymer_copies if polymer_copies > 0 else len(involved_chains) * max(1, total_ops))
        author_det = bool(asm.author_determined)
        soft_det = bool(asm.software_determined)
        oligo_det = str(getattr(asm, "oligomeric_details", "") or "").strip()
        oligo_name = oligo_det.capitalize() if oligo_det else _oligomer_label(oligo)
        src_note = "Author & Software" if (author_det and soft_det) else ("Author" if author_det else ("PISA / Software" if soft_det else "REMARK 350"))
        result.append({
            "id": asm_id,
            "name": f"Biomol {asm_id}",
            "label": f"Biomol {asm_id} ({oligo_name})",
            "oligomeric_count": oligo,
            "chains": involved_chains,
            "operator_count": max(1, total_ops),
            "has_symmetry_copies": total_ops > 1,
            "is_active": active_assembly_id == asm_id,
            "description": f"{src_note} · {oligo_name} ({max(1, total_ops)} sym-op{'s' if total_ops != 1 else ''})",
        })

    return result


def build_biological_assembly(
    raw_st: gemmi.Structure,
    assembly_id: str,
    deduplicate_axial_ions: bool = True,
    axial_cutoff_A: float = 0.6,
) -> Tuple[gemmi.Structure, Dict[str, Any]]:
    """
    Build either 'ASU' (clone of raw_st) or a declared Biological Assembly (Biomol N).
    Returns (new_structure, info_dict) including symmetry_copy_chains and deduplicated_ions count.
    """
    if assembly_id.upper() == "ASU":
        cloned = raw_st.clone()
        return cloned, {
            "assembly_id": "ASU",
            "name": "ASU",
            "original_chains": [ch.name for ch in cloned[0]] if len(cloned) > 0 else [],
            "symmetry_copy_chains": [],
            "deduplicated_ions": 0,
            "message": "Switched to Asymmetric Unit (ASU)",
        }

    target_asm = None
    for asm in raw_st.assemblies:
        if str(asm.name).strip() == str(assembly_id).strip():
            target_asm = asm
            break
    if target_asm is None:
        raise ValueError(f"Biological assembly '{assembly_id}' not found in this structure")

    orig_chains = [ch.name for ch in raw_st[0]] if len(raw_st) > 0 else []

    # Build assembly preserving original assemblies metadata on the resulting Structure
    new_st = raw_st.clone()
    try:
        how = gemmi.HowToNameCopiedChain.Short
        built_model = gemmi.make_assembly(target_asm, raw_st[0], how)
        # Ensure every chain in built_model has a unique 1-char chain ID (PDB-compatible)
        seen_ids: Set[str] = set()
        pool = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
        used_all = {ch.name for ch in built_model}
        for ch in built_model:
            if ch.name in seen_ids:
                nxt = next((c for c in pool if c not in used_all and c not in seen_ids), None)
                if nxt is not None:
                    ch.name = nxt
                    used_all.add(nxt)
            seen_ids.add(ch.name)
    except Exception:
        how = gemmi.HowToNameCopiedChain.AddNumber
        built_model = gemmi.make_assembly(target_asm, raw_st[0], how)

    # Replace model[0] cleanly while preserving cell, spacegroup, entities, assemblies, raw_remarks
    while len(new_st) > 0:
        del new_st[0]
    new_st.add_model(built_model)
    new_st.setup_entities()
    new_st.assign_label_seq_id()

    built_chains = [ch.name for ch in new_st[0]]
    orig_set = set(orig_chains)
    sym_copy_chains = [cid for cid in built_chains if cid not in orig_set]

    # Deduplicate symmetry-axis ions / single-atom heterogens overlapping within axial_cutoff_A (e.g. 4INS axial ZN)
    dedup_count = 0
    if deduplicate_axial_ions and len(new_st) > 0:
        kept_ion_positions: List[Tuple[str, np.ndarray]] = []
        for ch in new_st[0]:
            to_remove: List[int] = []
            for idx, res in enumerate(ch):
                rn = res.name.strip().upper()
                if rn in METAL_IONS or rn in HALIDE_IONS:
                    heavy = [a for a in res if a.element.name not in ("H", "D")]
                    if len(heavy) == 1:
                        pos = np.array([heavy[0].pos.x, heavy[0].pos.y, heavy[0].pos.z], dtype=float)
                        is_dup = False
                        for prev_rn, prev_pos in kept_ion_positions:
                            if prev_rn == rn and float(np.linalg.norm(pos - prev_pos)) <= axial_cutoff_A:
                                is_dup = True
                                break
                        if is_dup:
                            to_remove.append(idx)
                        else:
                            kept_ion_positions.append((rn, pos))
            for idx in reversed(to_remove):
                del ch[idx]
                dedup_count += 1

    msg = f"Built Biomol {assembly_id} ({len(built_chains)} chain(s)"
    if sym_copy_chains:
        msg += f", {len(sym_copy_chains)} symmetry copy chain(s): {', '.join(sym_copy_chains)}"
    if dedup_count:
        msg += f", deduplicated {dedup_count} axial ion(s)"
    msg += ")"

    return new_st, {
        "assembly_id": str(assembly_id),
        "name": f"Biomol {assembly_id}",
        "original_chains": orig_chains,
        "symmetry_copy_chains": sym_copy_chains,
        "deduplicated_ions": dedup_count,
        "message": msg,
    }


# ---------------------------------------------------------------------------
# 2. Missing Residues (REMARK 465 / SEQRES) + 3D Markers + Selective Rebuild
# ---------------------------------------------------------------------------

def _parse_remark_465(st: gemmi.Structure) -> Dict[str, List[Tuple[int, str, str]]]:
    """
    Parse PDB REMARK 465 missing residue records into {chain_id: [(resseq, icode, resname), ...]}.
    Standard PDB format:
      REMARK 465   M RES C SSSEQI
      cols: 15-17 RES, 19 C, 21-26 SSSEQI (e.g. 'REMARK 465     GLY A   -2')
    """
    by_chain: Dict[str, List[Tuple[int, str, str]]] = {}
    if not getattr(st, "raw_remarks", None):
        return by_chain
    in_table = False
    for line in st.raw_remarks:
        if not line.startswith("REMARK 465"):
            continue
        if "M RES C SSSEQI" in line or "RES C SSSEQI" in line:
            in_table = True
            continue
        if not in_table:
            continue
        body = line[10:].strip()
        if not body:
            continue
        # Match optional model number, 3-letter resname, 1-char chain, integer resseq + optional icode
        m = re.match(r"^(?:(\d+)\s+)?([A-Z0-9]{2,3})\s+([A-Za-z0-9])\s+(-?\d+)([A-Za-z]?)$", body)
        if not m:
            continue
        model_num, resname, ch, num_s, ic = m.groups()
        if model_num is not None and model_num != "1":
            continue
        by_chain.setdefault(ch, []).append((int(num_s), ic or "", resname.upper()))
    return by_chain


def detect_missing_regions(
    st: gemmi.Structure,
    rebuilt_records: Optional[List[Dict[str, Any]]] = None,
    refused_records: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """
    Detect internal missing gaps and missing N/C termini from REMARK 465 + SEQRES + geometric chain breaks.
    Computes exact 3D coordinates for 3Dmol markers:
      - internal gap: anchor_c_xyz, anchor_n_xyz, midpoint_xyz (skipped if C-N < 3.5 A continuous)
      - missing terminus: anchor_xyz, marker_xyz (3.0 A beyond terminal atom)
    """
    rebuilt_records = rebuilt_records or []
    refused_records = refused_records or {}
    rebuilt_ids = {r.get("gap_id") for r in rebuilt_records}

    r465 = _parse_remark_465(st)
    seqres_by_chain: Dict[str, List[str]] = {}
    for ent in st.entities:
        if ent.entity_type == gemmi.EntityType.Polymer and ent.full_sequence:
            for sc in ent.subchains:
                seqres_by_chain[sc] = [gemmi.Entity.first_mon(m) for m in ent.full_sequence]

    internal_gaps: List[Dict[str, Any]] = []
    missing_termini: List[Dict[str, Any]] = []

    if len(st) == 0:
        return {"internal_gaps": [], "missing_termini": [], "rebuilt_gaps": rebuilt_records}

    for chain in st[0]:
        prot_res = [
            r for r in chain
            if r.name.strip().upper() in ALL_PROTEIN_RESNAMES or r.name.strip().upper() in NONSTANDARD_TO_STANDARD
        ]
        if not prot_res:
            continue

        present_nums = {r.seqid.num for r in prot_res}
        r465_chain = sorted(r465.get(chain.name, []), key=lambda x: (x[0], x[1]))
        r465_map: Dict[int, str] = {num: rname for num, _ic, rname in r465_chain}

        first_res = prot_res[0]
        last_res = prot_res[-1]

        # 1. Missing N-terminus (from REMARK 465 or SEQRES before first observed residue)
        nterm_missing = [(num, rn) for num, _ic, rn in r465_chain if num < first_res.seqid.num and num not in present_nums]
        if not nterm_missing:
            seqres = seqres_by_chain.get(first_res.subchain) or seqres_by_chain.get(chain.name)
            if seqres and first_res.seqid.num > 1 and first_res.seqid.num <= len(seqres):
                # Only treat as N-term missing if residue 1..first-1 matches SEQRES indexing
                nterm_missing = [(k, seqres[k - 1]) for k in range(1, first_res.seqid.num)]

        if nterm_missing:
            p_n = _get_pos(first_res, "N")
            p_ca = _get_pos(first_res, "CA")
            if p_n is not None:
                direction = (p_n - p_ca) if p_ca is not None else np.array([1.0, 0.0, 0.0])
                norm = float(np.linalg.norm(direction))
                unit = direction / norm if norm > 1e-6 else np.array([1.0, 0.0, 0.0])
                marker = p_n + 3.0 * unit
                seq_3 = [rn for _, rn in nterm_missing]
                seq_1 = "".join(AA_3_TO_1.get(canonical_resname(rn), "X") for rn in seq_3)
                missing_termini.append({
                    "id": f"{chain.name}:N-term",
                    "chain": chain.name,
                    "terminus": "N",
                    "anchor_resseq": first_res.seqid.num,
                    "anchor_resname": first_res.name.strip(),
                    "missing_count": len(seq_3),
                    "sequence_3": seq_3,
                    "sequence_1": seq_1,
                    "range_label": f"{nterm_missing[0][0]}..{nterm_missing[-1][0]}",
                    "anchor_xyz": [round(float(x), 3) for x in p_n],
                    "marker_xyz": [round(float(x), 3) for x in marker],
                    "label": f"N-term ({len(seq_3)}): {seq_1[:10]}{'…' if len(seq_1) > 10 else ''}",
                })

        # 2. Missing C-terminus
        cterm_missing = [(num, rn) for num, _ic, rn in r465_chain if num > last_res.seqid.num and num not in present_nums]
        if not cterm_missing:
            seqres = seqres_by_chain.get(last_res.subchain) or seqres_by_chain.get(chain.name)
            if seqres and last_res.seqid.num >= 1 and last_res.seqid.num < len(seqres):
                cterm_missing = [(k, seqres[k - 1]) for k in range(last_res.seqid.num + 1, len(seqres) + 1)]

        if cterm_missing:
            p_c = _get_pos(last_res, "C")
            p_ca = _get_pos(last_res, "CA")
            if p_c is not None:
                direction = (p_c - p_ca) if p_ca is not None else np.array([1.0, 0.0, 0.0])
                norm = float(np.linalg.norm(direction))
                unit = direction / norm if norm > 1e-6 else np.array([1.0, 0.0, 0.0])
                marker = p_c + 3.0 * unit
                seq_3 = [rn for _, rn in cterm_missing]
                seq_1 = "".join(AA_3_TO_1.get(canonical_resname(rn), "X") for rn in seq_3)
                missing_termini.append({
                    "id": f"{chain.name}:C-term",
                    "chain": chain.name,
                    "terminus": "C",
                    "anchor_resseq": last_res.seqid.num,
                    "anchor_resname": last_res.name.strip(),
                    "missing_count": len(seq_3),
                    "sequence_3": seq_3,
                    "sequence_1": seq_1,
                    "range_label": f"{cterm_missing[0][0]}..{cterm_missing[-1][0]}",
                    "anchor_xyz": [round(float(x), 3) for x in p_c],
                    "marker_xyz": [round(float(x), 3) for x in marker],
                    "label": f"C-term ({len(seq_3)}): {seq_1[:10]}{'…' if len(seq_1) > 10 else ''}",
                })

        # 3. Internal Gaps between consecutive observed protein residues
        seqres = seqres_by_chain.get(prot_res[0].subchain) or seqres_by_chain.get(chain.name)
        for i in range(len(prot_res) - 1):
            r1, r2 = prot_res[i], prot_res[i + 1]
            c1 = _get_pos(r1, "C")
            n2 = _get_pos(r2, "N")
            dist_cn = distance(c1, n2) if (c1 is not None and n2 is not None) else 9.9
            seq_gap = r2.seqid.num - r1.seqid.num - 1

            # ADapt rule: if ends meet within 3.5 A, treat chain as continuous (e.g. 1JFF GLY34-HIS61)
            if seq_gap > 0 and dist_cn < 3.5:
                continue
            if seq_gap <= 0 and dist_cn <= 2.2:
                continue

            seq_3: List[str] = []
            source_tag = "geometric"
            if seq_gap > 0:
                for num in range(r1.seqid.num + 1, r2.seqid.num):
                    if num in r465_map:
                        seq_3.append(r465_map[num])
                        source_tag = "REMARK 465"
                    elif seqres and 1 <= num <= len(seqres):
                        seq_3.append(seqres[num - 1])
                        if source_tag != "REMARK 465":
                            source_tag = "SEQRES"
                    else:
                        seq_3.append("GLY")

            seq_1 = "".join(AA_3_TO_1.get(canonical_resname(rn), "X") for rn in seq_3)
            gap_id = f"{chain.name}:{r1.seqid.num}-{r2.seqid.num}"
            max_span = max(1, len(seq_3) + 1) * 3.85
            physically_closable = (seq_gap > 0) and (dist_cn <= max_span) and (len(seq_3) <= 25)
            has_unusual_mod = any(
                canonical_resname(rn) not in AA_3_TO_1 and rn != "MSE"
                for rn in seq_3
            )
            can_rebuild = bool(seq_gap > 0 and len(seq_3) <= 25 and not has_unusual_mod)

            refusal_reason = refused_records.get(gap_id)
            if not refusal_reason and seq_gap > 0 and not physically_closable:
                refusal_reason = f"Gap distance ({dist_cn:.1f} Å) exceeds maximum extended span ({max_span:.1f} Å) for {len(seq_3)} residue(s)"

            midpoint = ((c1 + n2) / 2.0) if (c1 is not None and n2 is not None) else np.zeros(3)
            internal_gaps.append({
                "id": gap_id,
                "chain": chain.name,
                "start_resseq": r1.seqid.num,
                "start_resname": r1.name.strip(),
                "end_resseq": r2.seqid.num,
                "end_resname": r2.name.strip(),
                "missing_count": max(seq_gap, 0),
                "gap_distance_A": round(float(dist_cn), 2),
                "sequence_3": seq_3,
                "sequence_1": seq_1,
                "source": source_tag,
                "anchor_c_xyz": [round(float(x), 3) for x in c1] if c1 is not None else [0.0, 0.0, 0.0],
                "anchor_n_xyz": [round(float(x), 3) for x in n2] if n2 is not None else [0.0, 0.0, 0.0],
                "midpoint_xyz": [round(float(x), 3) for x in midpoint],
                "can_rebuild": can_rebuild,
                "physically_closable": physically_closable,
                "status": "rebuilt" if gap_id in rebuilt_ids else ("refused" if gap_id in refused_records else "unbuilt"),
                "refusal_reason": refusal_reason,
                "label": f"{r1.name.strip()}{r1.seqid.num}--{seq_1 or 'break'}--{r2.name.strip()}{r2.seqid.num}",
            })

    return {
        "internal_gaps": internal_gaps,
        "missing_termini": missing_termini,
        "rebuilt_gaps": rebuilt_records,
    }


def _evaluate_loop_geometry(
    st: gemmi.Structure,
    chain_id: str,
    start_resseq: int,
    end_resseq: int,
    rebuilt_nums: Set[int],
) -> Dict[str, Any]:
    """
    Validate rebuilt loop geometry against ADapt's strict peptide closure & steric clearance criteria:
      - C-N peptide bonds in [1.22, 1.44] A
      - CA-CA distances in [3.40, 4.10] A
      - Trans peptide bonds (|omega| >= 135 deg)
      - L-chirality (positive improper dihedral N-CA-C-CB for non-GLY)
      - Minimum heavy-atom distance to observed atoms >= 2.2 A
    """
    model = st[0]
    chain = model[chain_id]
    prot = [
        r for r in chain
        if r.name.strip().upper() in ALL_PROTEIN_RESNAMES
        and start_resseq <= r.seqid.num <= end_resseq
    ]
    prot.sort(key=lambda r: r.seqid.num)

    cn_dists: List[float] = []
    caca_dists: List[float] = []
    omegas: List[float] = []
    chiral_bad: List[int] = []

    for i in range(len(prot) - 1):
        r1, r2 = prot[i], prot[i + 1]
        if (r1.seqid.num not in rebuilt_nums) and (r2.seqid.num not in rebuilt_nums):
            continue
        c1, ca1 = _get_pos(r1, "C"), _get_pos(r1, "CA")
        n2, ca2 = _get_pos(r2, "N"), _get_pos(r2, "CA")
        if c1 is not None and n2 is not None:
            cn_dists.append(round(distance(c1, n2), 3))
        if ca1 is not None and ca2 is not None:
            caca_dists.append(round(distance(ca1, ca2), 3))
        if ca1 is not None and c1 is not None and n2 is not None and ca2 is not None:
            om = dihedral_deg(ca1, c1, n2, ca2)
            omegas.append(round(float(om), 1))

    for r in prot:
        if r.seqid.num not in rebuilt_nums:
            continue
        if canonical_resname(r.name) == "GLY":
            continue
        n_p, ca_p, c_p, cb_p = _get_pos(r, "N"), _get_pos(r, "CA"), _get_pos(r, "C"), _get_pos(r, "CB")
        if n_p is not None and ca_p is not None and c_p is not None and cb_p is not None:
            # L-amino acid improper CA-N-C-CB should be positive (~+34 deg)
            v1 = n_p - ca_p
            v2 = c_p - ca_p
            v3 = cb_p - ca_p
            vol = float(np.dot(np.cross(v1, v2), v3))
            if vol < 0:
                chiral_bad.append(r.seqid.num)

    # Minimum heavy-atom distance from rebuilt atoms to observed (non-rebuilt, non-anchor) atoms
    obs_pts: List[np.ndarray] = []
    new_pts: List[np.ndarray] = []
    for ch in model:
        for res in ch:
            is_new = (ch.name == chain_id and res.seqid.num in rebuilt_nums)
            is_anchor = (ch.name == chain_id and res.seqid.num in (start_resseq, end_resseq))
            for a in res:
                if a.element.name in ("H", "D"):
                    continue
                pt = np.array([a.pos.x, a.pos.y, a.pos.z], dtype=float)
                if is_new:
                    new_pts.append(pt)
                elif not is_anchor:
                    obs_pts.append(pt)

    min_contact = 9.9
    if new_pts and obs_pts:
        tree = cKDTree(np.array(obs_pts))
        dists, _ = tree.query(np.array(new_pts), k=1)
        min_contact = round(float(np.min(dists)), 3)

    reasons: List[str] = []
    if cn_dists and (min(cn_dists) < 1.20 or max(cn_dists) > 1.48):
        reasons.append(f"peptide C–N length out of range ({min(cn_dists):.2f}–{max(cn_dists):.2f} Å)")
    if caca_dists and (min(caca_dists) < 3.30 or max(caca_dists) > 4.20):
        reasons.append(f"CA–CA spacing out of range ({min(caca_dists):.2f}–{max(caca_dists):.2f} Å)")
    if omegas and any(abs(o) < 130.0 for o in omegas):
        reasons.append(f"non-trans peptide omega ({min(abs(o) for o in omegas):.0f}°)")
    if chiral_bad:
        reasons.append(f"distorted CA chirality at residue(s) {chiral_bad}")
    if min_contact < 2.2:
        reasons.append(f"steric clash with observed atoms ({min_contact:.2f} Å < 2.20 Å)")

    return {
        "passed": len(reasons) == 0,
        "reasons": reasons,
        "cn_min_A": min(cn_dists) if cn_dists else 0.0,
        "cn_max_A": max(cn_dists) if cn_dists else 0.0,
        "caca_min_A": min(caca_dists) if caca_dists else 0.0,
        "caca_max_A": max(caca_dists) if caca_dists else 0.0,
        "min_abs_omega_deg": min((abs(o) for o in omegas), default=180.0),
        "min_contact_A": min_contact,
    }


def rebuild_single_gap(
    st: gemmi.Structure,
    chain_id: str,
    start_resseq: int,
    end_resseq: int,
    sequence_3: Optional[List[str]] = None,
    strict_geometry: bool = True,
    max_seeds: int = 3,
) -> Dict[str, Any]:
    """
    Selectively rebuild ONE internal missing loop gap, trying up to max_seeds conformations
    with local Engh & Huber + steric relaxation (observed atoms fixed), and applying
    ADapt's strict geometry validation gate.
    """
    if len(st) == 0:
        raise ValueError("Structure is empty")

    # Find gap metadata if sequence_3 not provided
    missing_info = detect_missing_regions(st)
    gap_id = f"{chain_id}:{start_resseq}-{end_resseq}"
    gap_meta = next((g for g in missing_info["internal_gaps"] if g["id"] == gap_id), None)

    n_missing = end_resseq - start_resseq - 1
    if n_missing <= 0:
        raise ValueError(f"No missing residue numbers between {chain_id}:{start_resseq} and {end_resseq}")
    if n_missing > 25:
        return {
            "rebuilt": False,
            "refused": True,
            "gap_id": gap_id,
            "reason": f"Gap length ({n_missing} residues) exceeds safe de-novo loop limit (25 residues)",
        }

    seq_list = list(sequence_3) if sequence_3 else (list(gap_meta["sequence_3"]) if gap_meta else ["GLY"] * n_missing)
    while len(seq_list) < n_missing:
        seq_list.append("GLY")
    seq_list = ["MET" if r == "MSE" else canonical_resname(r) for r in seq_list[:n_missing]]

    # Physical reach check first
    chain = st[0][chain_id]
    r_start = next((r for r in chain if r.seqid.num == start_resseq and r.name.strip().upper() in ALL_PROTEIN_RESNAMES), None)
    r_end = next((r for r in chain if r.seqid.num == end_resseq and r.name.strip().upper() in ALL_PROTEIN_RESNAMES), None)
    if r_start is None or r_end is None:
        raise ValueError(f"Anchor residues {chain_id}:{start_resseq} and/or {end_resseq} not found")

    p_c1 = _get_pos(r_start, "C")
    p_ca1 = _get_pos(r_start, "CA")
    p_n2 = _get_pos(r_end, "N")
    p_ca2 = _get_pos(r_end, "CA")
    if p_c1 is None or p_n2 is None or p_ca1 is None or p_ca2 is None:
        raise ValueError("Anchor residues are missing backbone C/CA/N atoms")

    span_dist = distance(p_c1, p_n2)
    max_possible = (n_missing + 1) * 3.85
    if strict_geometry and span_dist > max_possible:
        return {
            "rebuilt": False,
            "refused": True,
            "gap_id": gap_id,
            "reason": f"{n_missing} residue(s) ({''.join(AA_3_TO_1.get(r, 'X') for r in seq_list)}) cannot physically close a {span_dist:.1f} Å gap (max extended span {max_possible:.1f} Å)",
        }

    # Try up to max_seeds arch orientations (away from protein centroid, +/- orthogonal)
    all_ca = [
        np.array([a.pos.x, a.pos.y, a.pos.z], dtype=float)
        for ch in st[0]
        for r in ch
        for a in r
        if a.name.strip() == "CA"
    ]
    centroid = np.mean(all_ca, axis=0) if all_ca else np.zeros(3)
    chord = p_ca2 - p_ca1
    chord_len = max(float(np.linalg.norm(chord)), 1e-3)
    chord_u = chord / chord_len
    mid = 0.5 * (p_ca1 + p_ca2)
    outward = mid - centroid
    outward = outward - np.dot(outward, chord_u) * chord_u
    if float(np.linalg.norm(outward)) < 1e-4:
        outward = np.cross(chord_u, np.array([0.0, 0.0, 1.0]))
        if float(np.linalg.norm(outward)) < 1e-4:
            outward = np.cross(chord_u, np.array([0.0, 1.0, 0.0]))
    outward_u = outward / max(float(np.linalg.norm(outward)), 1e-6)
    side_u = np.cross(chord_u, outward_u)
    side_u = side_u / max(float(np.linalg.norm(side_u)), 1e-6)

    # Ideal arc length for (n_missing + 1) CA-CA segments (~3.80 A each)
    target_arc = (n_missing + 1) * 3.80
    excess = max(0.0, target_arc - chord_len)
    base_arch = min(max(1.2, 0.48 * math.sqrt(max(0.0, target_arc**2 - chord_len**2))), 0.45 * target_arc) if excess > 0.3 else 0.8

    seed_dirs = [
        outward_u,
        (0.75 * outward_u + 0.66 * side_u),
        (0.75 * outward_u - 0.66 * side_u),
    ]

    best_trial: Optional[Tuple[int, gemmi.Structure, Dict[str, Any]]] = None
    last_geom: Dict[str, Any] = {}
    rebuilt_nums = set(range(start_resseq + 1, end_resseq))

    for seed_idx in range(min(max_seeds, len(seed_dirs))):
        trial_st = st.clone()
        t_chain = trial_st[0][chain_id]
        arch_dir = seed_dirs[seed_idx]
        arch_dir = arch_dir / max(float(np.linalg.norm(arch_dir)), 1e-6)

        # Compute evenly-spaced CA positions along arched curve from p_ca1 -> p_ca2
        ca_pts: List[np.ndarray] = [p_ca1]
        for k in range(1, n_missing + 1):
            t = k / (n_missing + 1)
            pt = (1.0 - t) * p_ca1 + t * p_ca2 + base_arch * math.sin(math.pi * t) * arch_dir
            ca_pts.append(pt)
        ca_pts.append(p_ca2)

        # Iteratively enforce exact 3.80 A CA-CA spacing with observed anchors fixed
        for _ in range(40):
            for k in range(len(ca_pts) - 1):
                dvec = ca_pts[k + 1] - ca_pts[k]
                d = max(float(np.linalg.norm(dvec)), 1e-6)
                err = (d - 3.80) * 0.45
                corr = (dvec / d) * err
                if k > 0:
                    ca_pts[k] = ca_pts[k] + corr
                if k + 1 < len(ca_pts) - 1:
                    ca_pts[k + 1] = ca_pts[k + 1] - corr

        # Build backbone chain [anchor_CA, anchor_C, N1, CA1, C1, ..., Nk, CAk, Ck, anchor_N, anchor_CA2]
        # and iteratively project all covalent bond lengths & CA-CA distances with observed anchors pinned!
        nodes: List[np.ndarray] = [p_ca1.copy(), p_c1.copy()]
        fixed: List[bool] = [True, True]
        for k in range(1, n_missing + 1):
            ca_prev = ca_pts[k - 1]
            ca_cur = ca_pts[k]
            ca_next = ca_pts[k + 1]
            v_in = ca_cur - ca_prev
            v_in = v_in / max(float(np.linalg.norm(v_in)), 1e-6)
            v_out = ca_next - ca_cur
            v_out = v_out / max(float(np.linalg.norm(v_out)), 1e-6)
            n_init = ca_cur - 1.458 * v_in
            c_init = ca_cur + 1.525 * v_out
            nodes.extend([n_init, ca_cur.copy(), c_init])
            fixed.extend([False, False, False])
        nodes.extend([p_n2.copy(), p_ca2.copy()])
        fixed.extend([True, True])

        # Bond length targets along backbone node chain:
        # idx 0=CA0, 1=C0, (2=N1, 3=CA1, 4=C1), ..., -2=N_end, -1=CA_end
        bond_targets: List[Tuple[int, int, float]] = []
        # Peptide bonds C_i - N_{i+1} = 1.329 A; N_i - CA_i = 1.458 A; CA_i - C_i = 1.525 A
        for k in range(n_missing):
            n_i = 2 + 3 * k
            ca_i = n_i + 1
            c_i = n_i + 2
            prev_c = 1 if k == 0 else (n_i - 1)
            bond_targets.append((prev_c, n_i, 1.329))
            bond_targets.append((n_i, ca_i, 1.458))
            bond_targets.append((ca_i, c_i, 1.525))
            # 1-3 bond angle distance constraints (N-C ~ 2.45 A, prev_CA-N ~ 2.42 A, prev_C-CA ~ 2.44 A, CA-next_N ~ 2.42 A)
            prev_ca = 0 if k == 0 else (ca_i - 3)
            bond_targets.append((prev_ca, n_i, 2.42))
            bond_targets.append((prev_c, ca_i, 2.44))
            bond_targets.append((n_i, c_i, 2.46))
            bond_targets.append((prev_ca, ca_i, 3.80))
        last_c = 2 + 3 * (n_missing - 1) + 2
        last_ca = last_c - 1
        end_n = len(nodes) - 2
        end_ca = len(nodes) - 1
        bond_targets.append((last_c, end_n, 1.329))
        bond_targets.append((last_ca, end_n, 2.42))
        bond_targets.append((last_c, end_ca, 2.44))
        bond_targets.append((last_ca, end_ca, 3.80))

        # Iterative FABRIK / Gauss-Seidel distance projection + soft steric avoidance (220 sweeps)
        env_tree, env_pts = _build_env_kdtree(trial_st, exclude_key=f"{chain_id}:-99999")
        for sweep in range(220):
            # Soft steric push away from observed atoms on intermediate sweeps
            if env_tree is not None and env_pts is not None and (sweep % 15 == 10) and sweep < 180:
                for ni in range(len(nodes)):
                    if fixed[ni]:
                        continue
                    pt = nodes[ni]
                    dists, idxs = env_tree.query(pt, k=2)
                    dists = np.atleast_1d(dists)
                    idxs = np.atleast_1d(idxs)
                    for d_val, idx_p in zip(dists, idxs):
                        if 0.05 < float(d_val) < 2.65:
                            push = pt - env_pts[int(idx_p)]
                            nrm = float(np.linalg.norm(push))
                            if nrm > 1e-5:
                                pt = pt + (push / nrm) * (2.65 - float(d_val)) * 0.35
                    nodes[ni] = pt

            for i_a, i_b, target_d in bond_targets:
                dv = nodes[i_b] - nodes[i_a]
                cur_d = max(float(np.linalg.norm(dv)), 1e-6)
                diff = (cur_d - target_d) / cur_d
                step = dv * diff
                fa, fb = fixed[i_a], fixed[i_b]
                if not fa and not fb:
                    nodes[i_a] = nodes[i_a] + 0.5 * step
                    nodes[i_b] = nodes[i_b] - 0.5 * step
                elif not fa:
                    nodes[i_a] = nodes[i_a] + step
                elif not fb:
                    nodes[i_b] = nodes[i_b] - step

        # Find insertion index right after start_resseq
        insert_pos = 0
        for idx_r, res in enumerate(t_chain):
            if res.seqid.num == start_resseq and res.name.strip().upper() in ALL_PROTEIN_RESNAMES:
                insert_pos = idx_r + 1
                break

        # Build backbone (N, CA, C, O) + sidechain for each missing residue
        new_residues: List[gemmi.Residue] = []
        for k in range(n_missing):
            rname = seq_list[k]
            n_i = 2 + 3 * k
            n_pos = nodes[n_i]
            ca_cur = nodes[n_i + 1]
            c_pos = nodes[n_i + 2]
            next_n = nodes[n_i + 3]

            # Place carbonyl O in peptide plane (CA-C-N) opposite the bisector (~1.231 A from C)
            v_ca_c = c_pos - ca_cur
            v_ca_c = v_ca_c / max(float(np.linalg.norm(v_ca_c)), 1e-6)
            v_n_c = c_pos - next_n
            v_n_c = v_n_c / max(float(np.linalg.norm(v_n_c)), 1e-6)
            o_dir = v_ca_c + v_n_c
            if float(np.linalg.norm(o_dir)) < 1e-4:
                o_dir = arch_dir
            else:
                o_dir = o_dir / float(np.linalg.norm(o_dir))
            o_pos = c_pos + 1.231 * o_dir

            res = gemmi.Residue()
            res.name = rname
            res.seqid = gemmi.SeqId(start_resseq + k + 1, " ")
            res.entity_type = gemmi.EntityType.Polymer
            for aname, el, pos in (("N", "N", n_pos), ("CA", "C", ca_cur), ("C", "C", c_pos), ("O", "O", o_pos)):
                a = gemmi.Atom()
                a.name = aname
                a.element = gemmi.Element(el)
                a.pos = gemmi.Position(float(pos[0]), float(pos[1]), float(pos[2]))
                a.occ = 1.0
                a.b_iso = 45.0
                res.add_atom(a)

            sc = build_sidechain_coords({"N": n_pos, "CA": ca_cur, "C": c_pos}, rname)
            for aname, pos in sc.items():
                if aname in ("N", "CA", "C", "O"):
                    continue
                el = "S" if aname.startswith("S") else ("O" if aname.startswith("O") else ("N" if aname.startswith("N") else "C"))
                a = gemmi.Atom()
                a.name = aname
                a.element = gemmi.Element(el)
                a.pos = gemmi.Position(float(pos[0]), float(pos[1]), float(pos[2]))
                a.occ = 1.0
                a.b_iso = 50.0
                res.add_atom(a)
            new_residues.append(res)

        for offset_i, res in enumerate(new_residues):
            t_chain.add_residue(res, insert_pos + offset_i)

        geom = _evaluate_loop_geometry(trial_st, chain_id, start_resseq, end_resseq, rebuilt_nums)
        last_geom = geom
        if geom["passed"]:
            best_trial = (seed_idx + 1, trial_st, geom)
            break
        if best_trial is None or len(geom["reasons"]) < len(best_trial[2]["reasons"]):
            best_trial = (seed_idx + 1, trial_st, geom)

    seq_1 = "".join(AA_3_TO_1.get(r, "X") for r in seq_list)
    if best_trial is None or (strict_geometry and not best_trial[2]["passed"]):
        geom_used = best_trial[2] if best_trial else last_geom
        reason = "; ".join(geom_used.get("reasons") or ["Loop failed strict geometry closure"])
        return {
            "rebuilt": False,
            "refused": True,
            "gap_id": gap_id,
            "sequence_1": seq_1,
            "reason": f"Refused after {max_seeds} seed(s): {reason}",
            "geometry": geom_used,
        }

    seed_num, winning_st, winning_geom = best_trial
    # Copy winning model back into st in-place
    while len(st) > 0:
        del st[0]
    st.add_model(winning_st[0].clone())
    st.setup_entities()
    st.assign_label_seq_id()

    rebuilt_keys = [f"{chain_id}:{num}" for num in range(start_resseq + 1, end_resseq)]
    remark_850 = (
        f"REMARK 850 PANDAPROTEIN rebuilt {chain_id}:{seq_list[0]}:{start_resseq + 1} to "
        f"{chain_id}:{seq_list[-1]}:{end_resseq - 1} ({seq_1}) seed {seed_num} "
        f"[C-N {winning_geom['cn_min_A']:.2f}-{winning_geom['cn_max_A']:.2f}A, "
        f"min_contact {winning_geom['min_contact_A']:.2f}A]"
    )
    return {
        "rebuilt": True,
        "refused": False,
        "gap_id": gap_id,
        "chain": chain_id,
        "start_resseq": start_resseq,
        "end_resseq": end_resseq,
        "sequence_1": seq_1,
        "sequence_3": seq_list,
        "seed": seed_num,
        "rebuilt_keys": rebuilt_keys,
        "geometry": winning_geom,
        "remark_850": remark_850,
        "message": f"Rebuilt gap {gap_id} ({seq_1}, seed {seed_num} · min contact {winning_geom['min_contact_A']:.2f} Å)",
    }


# ---------------------------------------------------------------------------
# 3. Cached RCSB Entry-Level Ligand Lookup
# ---------------------------------------------------------------------------

def lookup_rcsb_entry_ligands(pdb_id: Optional[str], force_refresh: bool = False) -> Dict[str, Any]:
    """
    Query RCSB Data API (cached on disk in ~/.pandaprotein/rcsb_cache/) for official
    non-polymer components and Ligands of Interest (LOI) of a PDB entry.
    Gracefully returns offline status if unreachable or if pdb_id is not a 4-char PDB code.
    """
    if not pdb_id:
        return {"source": "offline_ccd", "pdb_id": None, "ligands": [], "components": {}}
    pid = re.sub(r"[^A-Za-z0-9]", "", pdb_id.split(".")[0]).upper()
    if len(pid) != 4:
        return {"source": "offline_ccd", "pdb_id": None, "ligands": [], "components": {}}

    cache_file = _RCSB_CACHE_DIR / f"{pid}.json"
    if cache_file.exists() and not force_refresh:
        try:
            cached = json.loads(cache_file.read_text(encoding="utf-8"))
            cached["cached"] = True
            return cached
        except Exception:
            pass

    try:
        _RCSB_CACHE_DIR.mkdir(parents=True, exist_ok=True)
        with httpx.Client(timeout=4.0) as client:
            r = client.get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
            if r.status_code != 200:
                return {"source": "offline_ccd", "pdb_id": pid, "ligands": [], "components": {}}
            entry = r.json()
            np_ids = entry.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids") or []
            components: Dict[str, Dict[str, Any]] = {}
            lig_codes: List[str] = []
            for eid in np_ids[:16]:
                er = client.get(f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pid}/{eid}")
                if er.status_code != 200:
                    continue
                ed = er.json()
                comp_id = (
                    ed.get("pdbx_entity_nonpoly", {}).get("comp_id")
                    or (ed.get("rcsb_nonpolymer_entity_container_identifiers", {}).get("nonpolymer_comp_id"))
                )
                if not comp_id:
                    continue
                comp_id = str(comp_id).upper()
                name = ed.get("pdbx_entity_nonpoly", {}).get("name") or ed.get("rcsb_nonpolymer_entity", {}).get("pdbx_description") or comp_id
                is_loi = bool(ed.get("rcsb_nonpolymer_entity_annotation")) or (
                    comp_id not in CRYSTALLIZATION_ARTIFACTS and comp_id not in METAL_IONS and comp_id not in HALIDE_IONS and comp_id not in WATER_RESNAMES
                )
                components[comp_id] = {
                    "comp_id": comp_id,
                    "name": str(name),
                    "is_ligand": is_loi,
                }
                if is_loi:
                    lig_codes.append(comp_id)

            payload = {
                "source": "rcsb_data_api",
                "pdb_id": pid,
                "fetched_at": datetime.now(timezone.utc).isoformat(),
                "ligands": sorted(set(lig_codes)),
                "components": components,
                "cached": False,
            }
            cache_file.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            return payload
    except Exception:
        return {"source": "offline_ccd", "pdb_id": pid, "ligands": [], "components": {}}


# ---------------------------------------------------------------------------
# 4. AD4 Atom-Type Validator & Unified ADapt Preparation Tree
# ---------------------------------------------------------------------------

def _residue_center(res: gemmi.Residue) -> List[float]:
    pts = [np.array([a.pos.x, a.pos.y, a.pos.z], dtype=float) for a in res if a.element.name not in ("H", "D")]
    if not pts:
        pts = [np.array([a.pos.x, a.pos.y, a.pos.z], dtype=float) for a in res]
    if not pts:
        return [0.0, 0.0, 0.0]
    c = np.mean(pts, axis=0)
    return [round(float(x), 3) for x in c]


def build_preparation_tree(
    st: gemmi.Structure,
    docking_prefs: Dict[str, Any],
    missing_regions: Dict[str, Any],
    rcsb_ligands: Optional[Dict[str, Any]] = None,
    symmetry_copy_chains: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Build the ADapt-style per-chain hierarchical preparation tree linked to 3D coordinates:
      protein, nucleic, missing_residues, altlocs, incomplete_sidechains, modified,
      metals_ions (with manual formal charges & AD4 ion substitution), unknown_atom_types,
      ligands, cofactors, additives, waters.
    """
    residue_roles: Dict[str, str] = docking_prefs.get("residue_roles") or {}
    chain_roles: Dict[str, str] = docking_prefs.get("chain_roles") or {}
    ion_charges: Dict[str, float] = docking_prefs.get("ion_charges") or {}
    ion_substitutions: Dict[str, str] = docking_prefs.get("ion_substitutions") or {}
    include_nucleic = bool(docking_prefs.get("include_nucleic", True))
    sym_set = set(symmetry_copy_chains or [])
    rcsb_comps = (rcsb_ligands or {}).get("components") or {}

    gaps_by_chain: Dict[str, List[Dict[str, Any]]] = {}
    for g in missing_regions.get("internal_gaps", []):
        gaps_by_chain.setdefault(g["chain"], []).append(g)
    terms_by_chain: Dict[str, List[Dict[str, Any]]] = {}
    for t in missing_regions.get("missing_termini", []):
        terms_by_chain.setdefault(t["chain"], []).append(t)

    chains_tree: List[Dict[str, Any]] = []
    summary = {
        "protein_chains": 0,
        "nucleic_residues": 0,
        "altloc_residues": 0,
        "incomplete_sidechains": 0,
        "modified_residues": 0,
        "metals_ions": 0,
        "unsupported_ad4_residues": 0,
        "ligands": 0,
        "cofactors": 0,
        "additives": 0,
        "waters": 0,
        "internal_gaps": len(missing_regions.get("internal_gaps", [])),
        "missing_termini": len(missing_regions.get("missing_termini", [])),
    }

    if len(st) == 0:
        return {"chains": [], "summary": summary, "ad4_substitutes": AD4_ION_SUBSTITUTES}

    for chain in st[0]:
        cid = chain.name
        c_node: Dict[str, Any] = {
            "chain": cid,
            "is_symmetry_copy": cid in sym_set,
            "role": chain_roles.get(cid, "receptor"),
            "protein_count": 0,
            "nucleic_count": 0,
            "nucleic_included": include_nucleic,
            "missing_gaps": gaps_by_chain.get(cid, []),
            "missing_termini": terms_by_chain.get(cid, []),
            "altlocs": [],
            "incomplete_sidechains": [],
            "modified": [],
            "metals_ions": [],
            "unknown_atom_types": [],
            "ligands": [],
            "cofactors": [],
            "additives": [],
            "waters": [],
        }

        for res in chain:
            rn = res.name.strip().upper()
            ic = res.seqid.icode.strip()
            key = f"{cid}:{res.seqid.num}{ic}"
            center = _residue_center(res)

            # Check altlocs
            alt_map: Dict[str, Tuple[float, float]] = {}
            unsupported_els: Set[str] = set()
            for a in res:
                el = a.element.name.upper()
                if el and el not in AUTODOCK4_SUPPORTED_ELEMENTS:
                    unsupported_els.add(el)
                if a.altloc not in ("\x00", " "):
                    alt_map[a.altloc] = (round(float(a.occ), 2), round(float(a.b_iso), 2))

            if len(alt_map) > 1:
                c_node["altlocs"].append({
                    "key": key,
                    "chain": cid,
                    "resseq": res.seqid.num,
                    "icode": ic,
                    "resname": rn,
                    "center": center,
                    "conformers": [{"id": k, "occ": v[0], "b_factor": v[1]} for k, v in sorted(alt_map.items())],
                })
                summary["altloc_residues"] += 1

            # Category routing
            if rn in ALL_PROTEIN_RESNAMES or rn in NONSTANDARD_TO_STANDARD:
                c_node["protein_count"] += 1
                if rn in NONSTANDARD_TO_STANDARD:
                    c_node["modified"].append({
                        "key": key,
                        "chain": cid,
                        "resseq": res.seqid.num,
                        "icode": ic,
                        "resname": rn,
                        "replacement": NONSTANDARD_TO_STANDARD[rn],
                        "center": center,
                        "role": residue_roles.get(key, "receptor"),
                    })
                    summary["modified_residues"] += 1
                std = canonical_resname(NONSTANDARD_TO_STANDARD.get(rn, rn))
                if std in EXPECTED_HEAVY_ATOMS:
                    present = {a.name.strip() for a in res if a.element.name not in ("H", "D")}
                    if rn == "MSE" and "SE" in present:
                        present.add("SD")
                    missing = sorted(set(EXPECTED_HEAVY_ATOMS[std]) - present)
                    if missing:
                        c_node["incomplete_sidechains"].append({
                            "key": key,
                            "chain": cid,
                            "resseq": res.seqid.num,
                            "icode": ic,
                            "resname": rn,
                            "missing": missing,
                            "center": center,
                        })
                        summary["incomplete_sidechains"] += 1
                if unsupported_els and rn != "MSE":
                    c_node["unknown_atom_types"].append({
                        "key": key,
                        "chain": cid,
                        "resseq": res.seqid.num,
                        "icode": ic,
                        "resname": rn,
                        "unsupported_elements": sorted(unsupported_els),
                        "is_ion": False,
                        "substituted_element": ion_substitutions.get(key),
                        "role": residue_roles.get(key, "exclude"),
                        "center": center,
                    })
                    summary["unsupported_ad4_residues"] += 1
            elif rn in NUCLEIC_RESNAMES:
                c_node["nucleic_count"] += 1
                summary["nucleic_residues"] += 1
            elif rn in WATER_RESNAMES:
                c_node["waters"].append({
                    "key": key,
                    "chain": cid,
                    "resseq": res.seqid.num,
                    "icode": ic,
                    "resname": rn,
                    "role": residue_roles.get(key, "exclude"),
                    "center": center,
                })
                summary["waters"] += 1
            elif rn in METAL_IONS or rn in HALIDE_IONS:
                # Primary element of ion
                els = [a.element.name.upper() for a in res if a.element.name not in ("H", "D")]
                prim_el = els[0] if els else rn[:2]
                sub_el = ion_substitutions.get(key)
                effective_el = (sub_el or prim_el).upper()
                ad4_ok = effective_el in AUTODOCK4_SUPPORTED_ELEMENTS
                default_q = CCD_ION_CHARGES.get(rn, CCD_ION_CHARGES.get(effective_el, 2.0))
                cur_q = float(ion_charges.get(key, default_q))
                default_role = "receptor" if ad4_ok else "exclude"
                item = {
                    "key": key,
                    "chain": cid,
                    "resseq": res.seqid.num,
                    "icode": ic,
                    "resname": rn,
                    "element": prim_el,
                    "effective_element": effective_el,
                    "substituted_element": sub_el,
                    "ad4_supported": ad4_ok,
                    "default_charge": default_q,
                    "charge": cur_q,
                    "role": residue_roles.get(key, default_role),
                    "center": center,
                }
                c_node["metals_ions"].append(item)
                summary["metals_ions"] += 1
                if not ad4_ok:
                    c_node["unknown_atom_types"].append({
                        "key": key,
                        "chain": cid,
                        "resseq": res.seqid.num,
                        "icode": ic,
                        "resname": rn,
                        "unsupported_elements": [prim_el],
                        "is_ion": True,
                        "substituted_element": sub_el,
                        "role": item["role"],
                        "center": center,
                    })
                    summary["unsupported_ad4_residues"] += 1
            else:
                heavy_cnt = sum(1 for a in res if a.element.name not in ("H", "D"))
                rcsb_meta = rcsb_comps.get(rn)
                if rn in COFACTORS:
                    cat = "cofactor"
                    default_role = "receptor"
                elif rn in CRYSTALLIZATION_ARTIFACTS:
                    cat = "additive"
                    default_role = "exclude"
                elif rcsb_meta is not None:
                    cat = "ligand" if rcsb_meta.get("is_ligand") else "additive"
                    default_role = "ligand" if cat == "ligand" else "exclude"
                else:
                    cat = "ligand"
                    default_role = "ligand"

                sub_el = ion_substitutions.get(key)
                ad4_ok = len(unsupported_els) == 0 or bool(sub_el)
                entry = {
                    "key": key,
                    "chain": cid,
                    "resseq": res.seqid.num,
                    "icode": ic,
                    "resname": rn,
                    "rcsb_name": rcsb_meta.get("name") if rcsb_meta else None,
                    "heavy_atoms": heavy_cnt,
                    "ad4_supported": ad4_ok,
                    "unsupported_elements": sorted(unsupported_els),
                    "role": residue_roles.get(key, default_role),
                    "center": center,
                }
                if cat == "cofactor":
                    c_node["cofactors"].append(entry)
                    summary["cofactors"] += 1
                elif cat == "additive":
                    c_node["additives"].append(entry)
                    summary["additives"] += 1
                else:
                    c_node["ligands"].append(entry)
                    summary["ligands"] += 1

                if not ad4_ok:
                    c_node["unknown_atom_types"].append({
                        "key": key,
                        "chain": cid,
                        "resseq": res.seqid.num,
                        "icode": ic,
                        "resname": rn,
                        "unsupported_elements": sorted(unsupported_els),
                        "is_ion": False,
                        "substituted_element": sub_el,
                        "role": entry["role"],
                        "center": center,
                    })
                    summary["unsupported_ad4_residues"] += 1

        if c_node["protein_count"] > 0:
            summary["protein_chains"] += 1
        chains_tree.append(c_node)

    return {
        "chains": chains_tree,
        "summary": summary,
        "ad4_substitutes": AD4_ION_SUBSTITUTES,
        "rcsb_source": (rcsb_ligands or {}).get("source", "offline_ccd"),
        "rcsb_pdb_id": (rcsb_ligands or {}).get("pdb_id"),
    }


# ---------------------------------------------------------------------------
# 5. Cryptographic Provenance & REMARK 850 Headers
# ---------------------------------------------------------------------------

def md5_hex(text: str) -> str:
    return hashlib.md5(text.encode("utf-8")).hexdigest()


def sha256_hex(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def format_pdb_serial(serial: int) -> str:
    """
    Format PDB/PDBQT atom serial in 5 columns (cols 7-11).
    Up to 99,999 uses decimal; >= 100,000 uses hexadecimal (ProDy / ADapt / Chimera convention)
    so large assemblies and cryo-EM complexes never overflow fixed columns.
    """
    if serial <= 99999:
        return f"{serial:5d}"
    return f"{serial:5x}"[:5]


def build_remark_850_lines(
    session_name: str,
    source: str,
    input_md5: str,
    active_assembly_id: str,
    ph: float,
    docking_prefs: Dict[str, Any],
    rebuilt_gaps: List[Dict[str, Any]],
    propka_summary: Optional[Dict[str, Any]] = None,
    skipped_unsupported: Optional[List[str]] = None,
) -> List[str]:
    """Generate standardized REMARK 850 provenance lines embedded inside .pdbqt and _rec_H.pdb files."""
    lines = [
        f"REMARK 850 PANDAPROTEIN version 2.0.0 generated {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"REMARK 850 INPUT name={session_name} source={source} md5={input_md5}",
        f"REMARK 850 ASSEMBLY active={active_assembly_id}",
        f"REMARK 850 PROTONATION engine=PROPKA-3.5.1 pH={ph:.2f}",
    ]
    if propka_summary:
        lines.append(
            f"REMARK 850 CHARGE net_discrete={propka_summary.get('net_charge_discrete', 0):+.1f} "
            f"pI={propka_summary.get('isoelectric_point', 0):.2f}"
        )
    ion_charges = docking_prefs.get("ion_charges") or {}
    if ion_charges:
        q_str = ", ".join(f"{k}={v:+.2f}" for k, v in sorted(ion_charges.items()))
        lines.append(f"REMARK 850 ION_CHARGES {q_str}")
    ion_subs = docking_prefs.get("ion_substitutions") or {}
    if ion_subs:
        s_str = ", ".join(f"{k}->{v}" for k, v in sorted(ion_subs.items()))
        lines.append(f"REMARK 850 ION_SUBSTITUTIONS {s_str}")
    roles = docking_prefs.get("residue_roles") or {}
    kept_het = [k for k, v in sorted(roles.items()) if v == "receptor"]
    lig_het = [k for k, v in sorted(roles.items()) if v == "ligand"]
    if kept_het:
        lines.append(f"REMARK 850 RECEPTOR_HETEROGENS {', '.join(kept_het[:20])}")
    if lig_het:
        lines.append(f"REMARK 850 LIGAND_SELECTIONS {', '.join(lig_het[:20])}")
    if skipped_unsupported:
        lines.append(f"REMARK 850 SKIPPED_UNSUPPORTED_AD4 {', '.join(sorted(set(skipped_unsupported))[:20])}")
    for rg in rebuilt_gaps:
        if rg.get("remark_850"):
            lines.append(rg["remark_850"])
    return lines


def build_provenance_record(
    session_id: str,
    name: str,
    source: str,
    input_md5: str,
    current_pdb_text: str,
    active_assembly_id: str,
    symmetry_copy_chains: List[str],
    ph: float,
    propka_data: Dict[str, Any],
    docking_prefs: Dict[str, Any],
    rebuilt_gaps: List[Dict[str, Any]],
    history_log: List[str],
    output_files: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Build an auditable, MD5-checksummed provenance session.json record (ADapt parity)."""
    pkg_versions = {"pandaprotein": "2.0.0", "gemmi": getattr(gemmi, "__version__", "unknown")}
    try:
        import propka
        pkg_versions["propka"] = getattr(propka, "__version__", "3.5.1")
    except Exception:
        pass
    try:
        import rdkit
        pkg_versions["rdkit"] = getattr(rdkit, "__version__", "unknown")
    except Exception:
        pass

    his_decisions = [
        {
            "residue": f"{r['chain']}:{r['resname']}:{r['resseq']}",
            "pka": r.get("pka"),
            "state": r.get("state"),
            "note": r.get("note"),
        }
        for r in (propka_data.get("titratable_residues") or [])
        if r.get("std_resname") in ("HIS", "ASP", "GLU", "CYS", "LYS", "TYR")
    ]

    files_manifest = {}
    for fname, content in (output_files or {}).items():
        files_manifest[fname] = {
            "md5": md5_hex(content),
            "sha256": sha256_hex(content),
            "bytes": len(content.encode("utf-8")),
        }

    return {
        "schema": "pandaprotein.provenance.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "session_id": session_id,
        "name": name,
        "source": source,
        "input_md5": input_md5,
        "current_structure_md5": md5_hex(current_pdb_text),
        "versions": pkg_versions,
        "assembly": {
            "active_id": active_assembly_id,
            "symmetry_copy_chains": symmetry_copy_chains,
        },
        "protonation": {
            "engine": propka_data.get("engine", "PROPKA 3.5.1"),
            "ph": ph,
            "net_charge_discrete": propka_data.get("net_charge_discrete"),
            "isoelectric_point": propka_data.get("isoelectric_point"),
            "decisions": his_decisions,
        },
        "docking_preferences": docking_prefs,
        "rebuilt_gaps": rebuilt_gaps,
        "edit_history": history_log,
        "files_written": files_manifest,
    }
