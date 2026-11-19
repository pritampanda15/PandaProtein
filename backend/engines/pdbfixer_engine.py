"""
Deep PDBFixer Engine:
Full diagnostics & structural repair pipeline for MD & Docking preparation:
- Non-standard residue detection & conversion (MSE->MET, SEP->SER, TPO->THR, etc.)
- Alternate conformation (AltLoc) resolution (highest occupancy / lowest B / A)
- Missing heavy atom reconstruction (NeRF + clash-free rotamer) & C-terminal OXT
- Missing internal loop / chain break detection (SEQRES + C-N geometric break) & kinematic loop builder
- N/C-terminal capping (ACE / NME) for truncated chains
- Disulfide bridge (S-S) detection, CYX assignment & SSBOND records
- Smart water curation (keep buried / catalytic / ligand-bridging waters, strip bulk)
- Heterogen curation (keep ligands/cofactors/metals, remove crystallization buffers)
- Chain filtering, insertion code normalization, and renumbering
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .chemistry_data import (
    ALL_PROTEIN_RESNAMES,
    COFACTORS,
    CRYSTALLIZATION_ARTIFACTS,
    EXPECTED_HEAVY_ATOMS,
    HALIDE_IONS,
    METAL_IONS,
    NONSTANDARD_TO_STANDARD,
    ROTAMER_LIBRARY,
    STANDARD_AA,
    WATER_RESNAMES,
)
from .coot_engine import _build_env_kdtree, _get_atom, _get_pos, _score_sidechain_clashes, regularize_zone
from .geometry_utils import build_sidechain_coords, canonical_resname, distance, nerf_place


def _is_protein_or_mod(res: gemmi.Residue) -> bool:
    r = res.name.strip().upper()
    return r in ALL_PROTEIN_RESNAMES or r in NONSTANDARD_TO_STANDARD


def _classify_heterogen(rname: str) -> str:
    r = rname.strip().upper()
    if r in WATER_RESNAMES:
        return "water"
    if r in METAL_IONS or r in HALIDE_IONS:
        return "metal"
    if r in COFACTORS:
        return "cofactor"
    if r in CRYSTALLIZATION_ARTIFACTS:
        return "buffer"
    return "ligand"


def diagnose_pdbfixer(st: gemmi.Structure, disulfide_cutoff: float = 2.5, water_keep_distance: float = 3.5) -> Dict[str, Any]:
    """Comprehensive pre-fix structural audit (PDBFixer diagnostic report)."""
    if len(st) == 0:
        return {
            "nonstandard_residues": [],
            "missing_atoms": [],
            "missing_loops": [],
            "altlocs": [],
            "disulfides": [],
            "heterogens": [],
            "waters": {"total": 0, "structural": 0, "bulk": 0},
            "insertion_codes": [],
            "chains": [],
            "has_hydrogens": False,
            "model_count": 0,
        }

    model = st[0]
    nonstandard = []
    missing_atoms_list = []
    missing_loops = []
    altlocs_list = []
    insertion_codes = []
    heterogens = []
    cys_sg_atoms: List[Tuple[str, int, str, str, np.ndarray]] = []
    ligand_metal_pts: List[np.ndarray] = []
    protein_polar_pts: List[np.ndarray] = []
    water_entries: List[Tuple[str, int, np.ndarray, float]] = []
    has_hydrogens = False
    chains_info = []

    # Map subchain -> SEQRES full sequence from entities
    seqres_by_chain: Dict[str, List[str]] = {}
    for ent in st.entities:
        if ent.entity_type == gemmi.EntityType.Polymer and ent.full_sequence:
            for sc in ent.subchains:
                seqres_by_chain[sc] = [gemmi.Entity.first_mon(m) for m in ent.full_sequence]

    for chain in model:
        prot_res = []
        for res in chain:
            rname = res.name.strip().upper()
            ic = res.seqid.icode.strip()
            if ic:
                insertion_codes.append({"chain": chain.name, "resseq": res.seqid.num, "icode": ic, "resname": rname})

            # Altloc inspection
            alt_map: Dict[str, Tuple[float, float]] = {}
            for a in res:
                if a.element.name in ("H", "D"):
                    has_hydrogens = True
                if a.altloc not in ("\x00", " "):
                    alt_map[a.altloc] = (round(a.occ, 2), round(a.b_iso, 2))
            if len(alt_map) > 1:
                altlocs_list.append({
                    "chain": chain.name,
                    "resseq": res.seqid.num,
                    "icode": ic,
                    "resname": rname,
                    "conformers": [{"id": k, "occ": v[0], "b_factor": v[1]} for k, v in sorted(alt_map.items())],
                })

            if _is_protein_or_mod(res):
                prot_res.append(res)
                if rname in NONSTANDARD_TO_STANDARD:
                    nonstandard.append({
                        "chain": chain.name,
                        "resseq": res.seqid.num,
                        "icode": ic,
                        "resname": rname,
                        "replacement": NONSTANDARD_TO_STANDARD[rname],
                    })
                std = canonical_resname(NONSTANDARD_TO_STANDARD.get(rname, rname))
                if std in EXPECTED_HEAVY_ATOMS:
                    present = {a.name.strip() for a in res if a.element.name not in ("H", "D")}
                    if rname == "MSE" and "SE" in present:
                        present.add("SD")
                    missing = sorted(set(EXPECTED_HEAVY_ATOMS[std]) - present)
                    if missing:
                        missing_atoms_list.append({
                            "chain": chain.name,
                            "resseq": res.seqid.num,
                            "icode": ic,
                            "resname": rname,
                            "missing": missing,
                        })
                if std == "CYS":
                    sg = _get_pos(res, "SG")
                    if sg is not None:
                        cys_sg_atoms.append((chain.name, res.seqid.num, ic, rname, sg))
                for a in res:
                    if a.element.name in ("N", "O"):
                        protein_polar_pts.append(np.array([a.pos.x, a.pos.y, a.pos.z]))
            else:
                cat = _classify_heterogen(rname)
                if cat == "water":
                    o_p = _get_pos(res, "O")
                    if o_p is None:
                        o_p = _get_pos(res, "OW")
                    if o_p is not None:
                        bval = next((a.b_iso for a in res), 30.0)
                        water_entries.append((chain.name, res.seqid.num, o_p, bval))
                else:
                    heavy_cnt = sum(1 for a in res if a.element.name not in ("H", "D"))
                    heterogens.append({
                        "chain": chain.name,
                        "resseq": res.seqid.num,
                        "icode": ic,
                        "resname": rname,
                        "category": cat,
                        "heavy_atoms": heavy_cnt,
                        "keep_recommended": cat in ("ligand", "cofactor", "metal"),
                    })
                    if cat in ("ligand", "cofactor", "metal"):
                        for a in res:
                            if a.element.name not in ("H", "D"):
                                ligand_metal_pts.append(np.array([a.pos.x, a.pos.y, a.pos.z]))

        # Detect chain breaks / missing internal residues
        chains_info.append({
            "id": chain.name,
            "protein_residues": len(prot_res),
            "first_res": prot_res[0].seqid.num if prot_res else None,
            "last_res": prot_res[-1].seqid.num if prot_res else None,
        })
        if prot_res:
            # Check missing C-terminal OXT
            last_r = prot_res[-1]
            if last_r.name.strip().upper() not in ("NME", "NMA", "NH2") and _get_atom(last_r, "OXT") is None:
                missing_atoms_list.append({
                    "chain": chain.name,
                    "resseq": last_r.seqid.num,
                    "icode": last_r.seqid.icode.strip(),
                    "resname": last_r.name.strip(),
                    "missing": ["OXT"],
                })

            seqres = seqres_by_chain.get(prot_res[0].subchain) or seqres_by_chain.get(chain.name)
            for i in range(len(prot_res) - 1):
                r1, r2 = prot_res[i], prot_res[i + 1]
                c1 = _get_pos(r1, "C")
                n2 = _get_pos(r2, "N")
                seq_gap = r2.seqid.num - r1.seqid.num - 1
                dist_cn = distance(c1, n2) if (c1 is not None and n2 is not None) else 9.9
                if seq_gap > 0 or dist_cn > 2.2:
                    missing_seq = []
                    if seqres and seq_gap > 0:
                        for num in range(r1.seqid.num + 1, r2.seqid.num):
                            if 1 <= num <= len(seqres):
                                missing_seq.append(seqres[num - 1])
                            else:
                                missing_seq.append("GLY")
                    elif seq_gap > 0:
                        missing_seq = ["GLY"] * seq_gap
                    missing_loops.append({
                        "chain": chain.name,
                        "start_resseq": r1.seqid.num,
                        "end_resseq": r2.seqid.num,
                        "missing_count": max(seq_gap, 0),
                        "gap_distance_A": round(dist_cn, 2),
                        "sequence": missing_seq,
                        "type": "missing_residues" if seq_gap > 0 else "peptide_chain_break",
                    })

    # Disulfide detection (SG-SG distance up to 3.2 Å; <= cutoff is active S-S bond)
    disulfides = []
    for i in range(len(cys_sg_atoms)):
        for j in range(i + 1, len(cys_sg_atoms)):
            c1, n1, ic1, rn1, p1 = cys_sg_atoms[i]
            c2, n2, ic2, rn2, p2 = cys_sg_atoms[j]
            d = distance(p1, p2)
            if d <= 3.2:
                disulfides.append({
                    "cys1": f"{c1}:{n1}{ic1}",
                    "cys2": f"{c2}:{n2}{ic2}",
                    "resname1": rn1,
                    "resname2": rn2,
                    "distance_A": round(d, 3),
                    "bonded": d <= disulfide_cutoff,
                })

    # Smart water classification: structural (within water_keep_distance of ligand/metal OR >=3 protein H-bonds) vs bulk
    structural_waters = 0
    lig_tree = cKDTree(np.array(ligand_metal_pts)) if ligand_metal_pts else None
    prot_tree = cKDTree(np.array(protein_polar_pts)) if protein_polar_pts else None
    for _, _, wp, _ in water_entries:
        near_lig = lig_tree is not None and len(lig_tree.query_ball_point(wp, water_keep_distance)) > 0
        hb_prot = prot_tree is not None and len(prot_tree.query_ball_point(wp, 3.3)) >= 3
        if near_lig or hb_prot:
            structural_waters += 1

    return {
        "nonstandard_residues": nonstandard,
        "missing_atoms": missing_atoms_list,
        "missing_loops": missing_loops,
        "altlocs": altlocs_list,
        "disulfides": disulfides,
        "heterogens": heterogens,
        "waters": {
            "total": len(water_entries),
            "structural": structural_waters,
            "bulk": len(water_entries) - structural_waters,
        },
        "insertion_codes": insertion_codes,
        "chains": chains_info,
        "has_hydrogens": has_hydrogens,
        "model_count": len(st),
    }


def _resolve_residue_altlocs(res: gemmi.Residue, mode: str) -> bool:
    if mode == "keep_all":
        return False
    altlocs = sorted({a.altloc for a in res if a.altloc not in ("\x00", " ")})
    if len(altlocs) <= 1:
        # Still clean stray altloc flags
        changed = False
        for a in res:
            if a.altloc not in ("\x00", " "):
                a.altloc = "\x00"
                a.occ = 1.0
                changed = True
        return changed

    if mode == "keep_A":
        chosen = "A" if "A" in altlocs else altlocs[0]
    elif mode == "lowest_bfactor":
        chosen = min(
            altlocs,
            key=lambda al: float(np.mean([a.b_iso for a in res if a.altloc == al] or [999.0])),
        )
    else:  # highest_occupancy
        chosen = max(
            altlocs,
            key=lambda al: float(np.mean([a.occ for a in res if a.altloc == al] or [0.0])),
        )

    kept = []
    seen = set()
    for a in res:
        al = a.altloc if a.altloc not in ("\x00", " ") else ""
        if al == "" or al == chosen:
            nm = a.name.strip()
            if nm not in seen:
                c = a.clone()
                c.altloc = "\x00"
                c.occ = 1.0
                kept.append(c)
                seen.add(nm)
    while len(res) > 0:
        del res[0]
    for a in kept:
        res.add_atom(a)
    return True


def _build_missing_loop_residues(
    r_prev: gemmi.Residue,
    r_next: gemmi.Residue,
    missing_seq: List[str],
    subchain: str,
) -> List[gemmi.Residue]:
    """
    Construct geometrically smooth missing loop residues bridging r_prev.C → r_next.N
    using curved cubic interpolation + NeRF peptide plane placement + ideal side chains.
    """
    c_start = _get_pos(r_prev, "C")
    ca_start = _get_pos(r_prev, "CA")
    n_end = _get_pos(r_next, "N")
    ca_end = _get_pos(r_next, "CA")
    if any(p is None for p in (c_start, ca_start, n_end, ca_end)):
        return []

    k = len(missing_seq)
    if k == 0:
        return []

    # Tangent vectors extending smoothly out from flanking CA->C and CA<-N
    t_start = c_start - ca_start
    t_start = t_start / (np.linalg.norm(t_start) + 1e-8)
    t_end = ca_end - n_end
    t_end = t_end / (np.linalg.norm(t_end) + 1e-8)

    span = distance(ca_start, ca_end)
    arc_scale = max(span * 0.45, k * 1.4)
    # Build smooth Hermite curve for CA positions
    built: List[gemmi.Residue] = []
    ca_pts: List[np.ndarray] = []
    for idx in range(1, k + 1):
        u = idx / (k + 1)
        h00 = 2 * u**3 - 3 * u**2 + 1
        h10 = u**3 - 2 * u**2 + u
        h01 = -2 * u**3 + 3 * u**2
        h11 = u**3 - u**2
        p_ca = h00 * ca_start + h10 * (t_start * arc_scale) + h01 * ca_end + h11 * (t_end * arc_scale)
        ca_pts.append(p_ca)

    all_cas = [ca_start] + ca_pts + [ca_end]
    mean_b = 45.0

    for idx, resname in enumerate(missing_seq):
        std = canonical_resname(NONSTANDARD_TO_STANDARD.get(resname.upper(), resname.upper()))
        if std not in STANDARD_AA:
            std = "GLY"
        ca_prev_p = all_cas[idx]
        ca_curr_p = all_cas[idx + 1]
        ca_next_p = all_cas[idx + 2]

        v_in = ca_curr_p - ca_prev_p
        v_in = v_in / (np.linalg.norm(v_in) + 1e-8)
        v_out = ca_next_p - ca_curr_p
        v_out = v_out / (np.linalg.norm(v_out) + 1e-8)

        # Place N, CA, C, O around CA_curr with standard backbone distances
        n_p = ca_curr_p - v_in * 1.458
        c_p = ca_curr_p + v_out * 1.525
        o_p = nerf_place(n_p, ca_curr_p, c_p, 1.231, 120.8, 180.0)

        bb = {"N": n_p, "CA": ca_curr_p, "C": c_p, "O": o_p}
        coords = build_sidechain_coords(bb, std)

        new_res = gemmi.Residue()
        new_res.name = std
        new_res.seqid = gemmi.SeqId(r_prev.seqid.num + idx + 1, " ")
        new_res.subchain = subchain
        new_res.entity_type = gemmi.EntityType.Polymer
        for aname, pos in coords.items():
            atom = gemmi.Atom()
            atom.name = aname
            atom.element = gemmi.Element(aname[0].upper())
            atom.pos = gemmi.Position(float(pos[0]), float(pos[1]), float(pos[2]))
            atom.occ = 1.0
            atom.b_iso = mean_b
            new_res.add_atom(atom)
        built.append(new_res)

    return built


def _add_terminal_caps(chain: gemmi.Chain) -> Tuple[bool, bool]:
    """Add N-terminal ACE and C-terminal NME capping residues to a protein chain."""
    prot_indices = [i for i, r in enumerate(chain) if _is_protein_or_mod(r)]
    if not prot_indices:
        return False, False

    added_ace = False
    added_nme = False

    # N-terminal ACE
    first_idx = prot_indices[0]
    first_res = chain[first_idx]
    if first_res.name.strip().upper() != "ACE":
        n_p, ca_p, c_p = _get_pos(first_res, "N"), _get_pos(first_res, "CA"), _get_pos(first_res, "C")
        if all(p is not None for p in (n_p, ca_p, c_p)):
            c_ace = nerf_place(c_p, ca_p, n_p, 1.329, 121.7, -60.0)
            o_ace = nerf_place(ca_p, n_p, c_ace, 1.231, 123.0, 0.0)
            ch3_ace = nerf_place(ca_p, n_p, c_ace, 1.520, 116.2, 180.0)
            ace = gemmi.Residue()
            ace.name = "ACE"
            ace.seqid = gemmi.SeqId(max(0, first_res.seqid.num - 1), " ")
            ace.subchain = first_res.subchain
            ace.entity_type = gemmi.EntityType.Polymer
            for aname, pos in [("CH3", ch3_ace), ("C", c_ace), ("O", o_ace)]:
                a = gemmi.Atom()
                a.name = aname
                a.element = gemmi.Element(aname[0])
                a.pos = gemmi.Position(*[float(x) for x in pos])
                a.occ = 1.0
                a.b_iso = 30.0
                ace.add_atom(a)
            chain.add_residue(ace, first_idx)
            added_ace = True

    # Re-find last protein residue
    prot_indices = [i for i, r in enumerate(chain) if _is_protein_or_mod(r)]
    last_idx = prot_indices[-1]
    last_res = chain[last_idx]
    if last_res.name.strip().upper() not in ("NME", "NMA", "NH2"):
        # Remove OXT if adding NME
        for ai in range(len(last_res) - 1, -1, -1):
            if last_res[ai].name.strip() == "OXT":
                del last_res[ai]
        n_p, ca_p, c_p = _get_pos(last_res, "N"), _get_pos(last_res, "CA"), _get_pos(last_res, "C")
        if all(p is not None for p in (n_p, ca_p, c_p)):
            n_nme = nerf_place(n_p, ca_p, c_p, 1.329, 116.2, 180.0)
            c_nme = nerf_place(ca_p, c_p, n_nme, 1.458, 121.7, 180.0)
            nme = gemmi.Residue()
            nme.name = "NME"
            nme.seqid = gemmi.SeqId(last_res.seqid.num + 1, " ")
            nme.subchain = last_res.subchain
            nme.entity_type = gemmi.EntityType.Polymer
            for aname, pos in [("N", n_nme), ("C", c_nme)]:
                a = gemmi.Atom()
                a.name = aname
                a.element = gemmi.Element(aname[0])
                a.pos = gemmi.Position(*[float(x) for x in pos])
                a.occ = 1.0
                a.b_iso = 30.0
                nme.add_atom(a)
            chain.add_residue(nme, last_idx + 1)
            added_nme = True

    return added_ace, added_nme


def run_pdbfixer(st: gemmi.Structure, opts: Dict[str, Any]) -> Dict[str, Any]:
    """Execute the full configurable PDBFixer preparation pipeline in-place on `st`."""
    log: List[str] = []

    if len(st) == 0:
        raise ValueError("Structure contains no models")

    # 1. Keep only Model 1 for MD/Docking preparation
    if len(st) > 1:
        n_models = len(st)
        while len(st) > 1:
            del st[1]
        log.append(f"Retained Model 1 of {n_models} models")

    model = st[0]

    # 2. Chain filtering
    keep_chains = opts.get("keep_chains")
    if keep_chains:
        keep_set = set(keep_chains)
        removed_ch = []
        for ci in range(len(model) - 1, -1, -1):
            if model[ci].name not in keep_set:
                removed_ch.append(model[ci].name)
                del model[ci]
        if removed_ch:
            log.append(f"Removed chain(s): {', '.join(reversed(removed_ch))}")

    # 3. Resolve AltLocs
    altloc_mode = opts.get("resolve_altlocs", "highest_occupancy")
    alt_resolved = 0
    for chain in model:
        for res in chain:
            if _resolve_residue_altlocs(res, altloc_mode):
                alt_resolved += 1
    if alt_resolved:
        log.append(f"Resolved alternate conformations on {alt_resolved} residue(s) ({altloc_mode})")

    # 4. Replace non-standard residues
    if opts.get("replace_nonstandard", True):
        overrides = opts.get("nonstandard_overrides", {}) or {}
        conv_count = 0
        for chain in model:
            for res in chain:
                rname = res.name.strip().upper()
                target = overrides.get(rname) or NONSTANDARD_TO_STANDARD.get(rname)
                if target and target in STANDARD_AA:
                    if rname == "MSE":
                        for a in res:
                            if a.name.strip() == "SE":
                                a.name = "SD"
                                a.element = gemmi.Element("S")
                    # Strip non-canonical modification atoms (e.g., PO3 on SEP/TPO/PTR)
                    allowed = set(EXPECTED_HEAVY_ATOMS[target]) | {"OXT"}
                    for ai in range(len(res) - 1, -1, -1):
                        if res[ai].element.name not in ("H", "D") and res[ai].name.strip() not in allowed:
                            del res[ai]
                    res.name = target
                    res.het_flag = "A"
                    conv_count += 1
        if conv_count:
            log.append(f"Converted {conv_count} non-standard residue(s) to canonical amino acids")

    # 5. Standardize atom names (e.g., ILE CD -> CD1)
    if opts.get("standardize_atom_names", True):
        renamed = 0
        for chain in model:
            for res in chain:
                if canonical_resname(res.name) == "ILE":
                    for a in res:
                        if a.name.strip() == "CD":
                            a.name = "CD1"
                            renamed += 1
        if renamed:
            log.append(f"Standardized {renamed} legacy atom name(s) (IUPAC/AMBER)")

    # 6. Pre-compute ligand/metal & protein polar points for smart water retention
    ligand_metal_pts = []
    protein_polar_pts = []
    for chain in model:
        for res in chain:
            rname = res.name.strip().upper()
            if _is_protein_or_mod(res):
                for a in res:
                    if a.element.name in ("N", "O"):
                        protein_polar_pts.append(np.array([a.pos.x, a.pos.y, a.pos.z]))
            else:
                cat = _classify_heterogen(rname)
                if cat in ("ligand", "cofactor", "metal"):
                    for a in res:
                        if a.element.name not in ("H", "D"):
                            ligand_metal_pts.append(np.array([a.pos.x, a.pos.y, a.pos.z]))

    lig_tree = cKDTree(np.array(ligand_metal_pts)) if ligand_metal_pts else None
    prot_tree = cKDTree(np.array(protein_polar_pts)) if protein_polar_pts else None

    # 7. Curate Heterogens & Waters
    water_mode = opts.get("remove_waters", "bulk")
    water_dist = float(opts.get("water_keep_distance", 3.5))
    remove_buffers = opts.get("remove_crystallization_artifacts", True)
    keep_ligands = opts.get("keep_ligands", True)
    keep_cofactors = opts.get("keep_cofactors", True)
    keep_metals = opts.get("keep_metals", True)

    removed_waters = 0
    kept_waters = 0
    removed_het: List[str] = []

    for chain in model:
        for ri in range(len(chain) - 1, -1, -1):
            res = chain[ri]
            rname = res.name.strip().upper()
            if _is_protein_or_mod(res):
                continue
            cat = _classify_heterogen(rname)
            if cat == "water":
                if water_mode == "all":
                    del chain[ri]
                    removed_waters += 1
                elif water_mode == "bulk":
                    wp = _get_pos(res, "O")
                    if wp is None:
                        wp = _get_pos(res, "OW")
                    keep_w = False
                    if wp is not None:
                        if lig_tree is not None and len(lig_tree.query_ball_point(wp, water_dist)) > 0:
                            keep_w = True
                        elif prot_tree is not None and len(prot_tree.query_ball_point(wp, 3.3)) >= 3:
                            keep_w = True
                    if keep_w:
                        kept_waters += 1
                    else:
                        del chain[ri]
                        removed_waters += 1
                else:
                    kept_waters += 1
            elif cat == "buffer" and remove_buffers:
                removed_het.append(rname)
                del chain[ri]
            elif cat == "ligand" and not keep_ligands:
                removed_het.append(rname)
                del chain[ri]
            elif cat == "cofactor" and not keep_cofactors:
                removed_het.append(rname)
                del chain[ri]
            elif cat == "metal" and not keep_metals:
                removed_het.append(rname)
                del chain[ri]

    if removed_waters:
        log.append(f"Removed {removed_waters} bulk/solvent water(s); retained {kept_waters} structural/active-site water(s)")
    if removed_het:
        log.append(f"Removed {len(removed_het)} unwanted heterogen(s): {', '.join(sorted(set(removed_het)))}")

    # 8. Rebuild Missing Internal Loops / Residues
    if opts.get("add_missing_residues", True):
        max_loop = int(opts.get("max_loop_length", 12))
        seqres_by_chain: Dict[str, List[str]] = {}
        for ent in st.entities:
            if ent.entity_type == gemmi.EntityType.Polymer and ent.full_sequence:
                for sc in ent.subchains:
                    seqres_by_chain[sc] = [gemmi.Entity.first_mon(m) for m in ent.full_sequence]

        total_rebuilt_res = 0
        rebuilt_zones = []
        for chain in model:
            prot_positions = [(idx, r) for idx, r in enumerate(chain) if _is_protein_or_mod(r)]
            if len(prot_positions) < 2:
                continue
            seqres = seqres_by_chain.get(prot_positions[0][1].subchain) or seqres_by_chain.get(chain.name)
            insertions: List[Tuple[int, List[gemmi.Residue], int, int]] = []
            for p_i in range(len(prot_positions) - 1):
                _, r1 = prot_positions[p_i]
                idx2, r2 = prot_positions[p_i + 1]
                gap = r2.seqid.num - r1.seqid.num - 1
                if 1 <= gap <= max_loop:
                    mseq = []
                    for num in range(r1.seqid.num + 1, r2.seqid.num):
                        if seqres and 1 <= num <= len(seqres):
                            mseq.append(seqres[num - 1])
                        else:
                            mseq.append("GLY")
                    new_res_list = _build_missing_loop_residues(r1, r2, mseq, r1.subchain or chain.name)
                    if new_res_list:
                        insertions.append((idx2, new_res_list, r1.seqid.num, r2.seqid.num))

            # Insert from back to front so indices stay valid
            for idx2, new_res_list, s_num, e_num in reversed(insertions):
                for offset, nr in enumerate(new_res_list):
                    chain.add_residue(nr, idx2 + offset)
                total_rebuilt_res += len(new_res_list)
                rebuilt_zones.append((chain.name, s_num, e_num))

        for ch_id, s_num, e_num in rebuilt_zones:
            try:
                regularize_zone(st, ch_id, s_num, e_num, iterations=45)
            except Exception:
                pass
        if total_rebuilt_res:
            log.append(f"Rebuilt {total_rebuilt_res} missing loop residue(s) across {len(rebuilt_zones)} gap(s) with local regularization")

    # 9. Add Missing Heavy Atoms (+ C-terminal OXT)
    if opts.get("add_missing_atoms", True):
        atoms_added = 0
        residues_repaired = 0
        for chain in model:
            prot_res = [r for r in chain if _is_protein_or_mod(r)]
            for res in prot_res:
                std = canonical_resname(res.name)
                if std not in EXPECTED_HEAVY_ATOMS:
                    continue
                present = {a.name.strip() for a in res if a.element.name not in ("H", "D")}
                missing = [a for a in EXPECTED_HEAVY_ATOMS[std] if a not in present]
                if not missing:
                    continue
                bb = {bn: _get_pos(res, bn) for bn in ("N", "CA", "C", "O") if _get_pos(res, bn) is not None}
                if not all(k in bb for k in ("N", "CA", "C")):
                    continue

                # Preserve existing heavy atoms and build missing ones with clash-aware rotamer if needed
                existing_map = {a.name.strip(): np.array([a.pos.x, a.pos.y, a.pos.z]) for a in res if a.element.name not in ("H", "D")}
                lib = ROTAMER_LIBRARY.get(std, [])
                best_coords = build_sidechain_coords(existing_map, std)
                if lib and any(m not in ("N", "CA", "C", "O", "CB") for m in missing):
                    tree, radii = _build_env_kdtree(st, f"{chain.name}:{res.seqid.num}{res.seqid.icode.strip()}")
                    best_score = -1e9
                    for _, prob, chis in lib:
                        cand = build_sidechain_coords(bb, std, chi_targets=chis)
                        # Only score newly added atoms
                        sub = {m: cand[m] for m in missing if m in cand}
                        clashes, energy = _score_sidechain_clashes(sub, tree, radii)
                        score = prob - energy * 4.0 - clashes * 12.0
                        if score > best_score:
                            best_score = score
                            best_coords = cand

                mean_b = float(np.mean([a.b_iso for a in res] or [30.0]))
                added_here = 0
                for m in missing:
                    if m in best_coords:
                        pos = best_coords[m]
                        atom = gemmi.Atom()
                        atom.name = m
                        atom.element = gemmi.Element(m[0].upper())
                        atom.pos = gemmi.Position(float(pos[0]), float(pos[1]), float(pos[2]))
                        atom.occ = 1.0
                        atom.b_iso = round(mean_b, 2)
                        res.add_atom(atom)
                        atoms_added += 1
                        added_here += 1
                if added_here:
                    residues_repaired += 1

            # C-terminal OXT (unless NME capping will be applied)
            if prot_res and opts.get("cap_termini", "none") != "ACE_NME":
                last_r = prot_res[-1]
                if last_r.name.strip().upper() not in ("NME", "NMA", "NH2") and _get_atom(last_r, "OXT") is None:
                    ca_p, c_p, o_p = _get_pos(last_r, "CA"), _get_pos(last_r, "C"), _get_pos(last_r, "O")
                    if all(p is not None for p in (ca_p, c_p, o_p)):
                        oxt_p = nerf_place(o_p, ca_p, c_p, 1.250, 118.4, 180.0)
                        oxt = gemmi.Atom()
                        oxt.name = "OXT"
                        oxt.element = gemmi.Element("O")
                        oxt.pos = gemmi.Position(*[float(x) for x in oxt_p])
                        oxt.occ = 1.0
                        oxt.b_iso = 30.0
                        last_r.add_atom(oxt)
                        atoms_added += 1
        if atoms_added:
            log.append(f"Reconstructed {atoms_added} missing heavy atom(s) across {residues_repaired} residue(s) (+ OXT)")

    # 10. Terminal Capping (ACE / NME)
    if opts.get("cap_termini", "none") == "ACE_NME":
        n_capped = 0
        for chain in model:
            a_ok, n_ok = _add_terminal_caps(chain)
            if a_ok or n_ok:
                n_capped += 1
        if n_capped:
            log.append(f"Added N-terminal ACE and C-terminal NME caps to {n_capped} chain(s)")

    # 11. Assign Disulfide Bridges (CYX + SSBOND connections)
    disulfide_pairs = []
    if opts.get("assign_disulfides", True):
        cutoff = float(opts.get("disulfide_cutoff", 2.5))
        cys_list = []
        for chain in model:
            for res in chain:
                if canonical_resname(res.name) == "CYS":
                    sg = _get_pos(res, "SG")
                    if sg is not None:
                        cys_list.append((chain.name, res, sg))

        used = set()
        for i in range(len(cys_list)):
            if i in used:
                continue
            for j in range(i + 1, len(cys_list)):
                if j in used:
                    continue
                d = distance(cys_list[i][2], cys_list[j][2])
                if d <= cutoff:
                    used.add(i)
                    used.add(j)
                    c1, r1, _ = cys_list[i]
                    c2, r2, _ = cys_list[j]
                    r1.name = "CYX"
                    r2.name = "CYX"
                    # Remove thiol HG if present
                    for r in (r1, r2):
                        for ai in range(len(r) - 1, -1, -1):
                            if r[ai].name.strip() in ("HG", "HG1"):
                                del r[ai]
                    conn = gemmi.Connection()
                    conn.name = f"disulf{len(disulfide_pairs) + 1}"
                    conn.type = gemmi.ConnectionType.Disulf
                    conn.partner1 = gemmi.AtomAddress(c1, r1.seqid, "CYX", "SG")
                    conn.partner2 = gemmi.AtomAddress(c2, r2.seqid, "CYX", "SG")
                    st.connections.append(conn)
                    disulfide_pairs.append(f"{c1}:{r1.seqid.num}–{c2}:{r2.seqid.num} ({d:.2f} Å)")
                    break
        if disulfide_pairs:
            log.append(f"Assigned {len(disulfide_pairs)} disulfide bridge(s) as CYX: {', '.join(disulfide_pairs)}")

    # 12. Fix insertion codes & optional renumbering
    if opts.get("fix_insertion_codes", True) or opts.get("renumber_residues", False):
        renum = opts.get("renumber_residues", False)
        fixed_ic = 0
        for chain in model:
            counter = 1
            for res in chain:
                if renum:
                    res.seqid = gemmi.SeqId(counter, " ")
                    counter += 1
                elif res.seqid.icode.strip():
                    res.seqid = gemmi.SeqId(res.seqid.num, " ")
                    fixed_ic += 1
        if renum:
            log.append("Renumbered all residues sequentially starting from 1 per chain")
        elif fixed_ic:
            log.append(f"Normalized {fixed_ic} insertion code(s)")

    st.setup_entities()
    return {
        "actions": log,
        "disulfides_assigned": disulfide_pairs,
        "diagnostics_after": diagnose_pdbfixer(st),
    }
