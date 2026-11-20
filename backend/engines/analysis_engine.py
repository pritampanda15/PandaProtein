"""
Structural Superposition (Kabsch RMSD), Membrane Hydrophobic-Belt Orientation,
and Non-Covalent Interaction Network Engine (Salt Bridges, Pi-Pi Stacking, Interfaces).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .chemistry_data import ALL_PROTEIN_RESNAMES, NONSTANDARD_TO_STANDARD
from .geometry_utils import canonical_resname, distance
from .structure_io import structure_to_pdb

HYDROPHOBIC_TM = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
AROMATIC_RINGS: Dict[str, List[str]] = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}
ACIDIC_ATOMS: Dict[str, List[str]] = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}
BASIC_ATOMS: Dict[str, List[str]] = {
    "LYS": ["NZ"],
    "ARG": ["NE", "NH1", "NH2"],
    "HIS": ["ND1", "NE2"],
}


def _kabsch_align(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute optimal rotation R and translation t mapping P -> Q minimizing RMSD.
    Returns (R, centroid_P, centroid_Q) so aligned_P = (P - cP) @ R + cQ.
    """
    cP = np.mean(P, axis=0)
    cQ = np.mean(Q, axis=0)
    P0 = P - cP
    Q0 = Q - cQ
    H = P0.T @ Q0
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, -1.0 if d < 0 else 1.0])
    R = U @ D @ Vt
    return R, cP, cQ


def compute_structural_diff(
    st_ref: gemmi.Structure,
    st_cur: gemmi.Structure,
    ref_label: str = "Original Structure",
) -> Dict[str, Any]:
    """
    Align `st_ref` onto `st_cur` via Kabsch C-alpha superposition and compute
    C-alpha RMSD, backbone RMSD, all-heavy-atom RMSD, per-residue shifts, and aligned reference PDB.
    """
    if len(st_ref) == 0 or len(st_cur) == 0:
        raise ValueError("Cannot compare empty structure")

    def _collect_atoms(st: gemmi.Structure) -> Dict[Tuple[str, int, str], np.ndarray]:
        out = {}
        for ch in st[0]:
            for res in ch:
                rn = res.name.strip().upper()
                if rn not in ALL_PROTEIN_RESNAMES and rn not in NONSTANDARD_TO_STANDARD:
                    continue
                for a in res:
                    if a.element.name in ("H", "D"):
                        continue
                    k = (ch.name, res.seqid.num, a.name.strip())
                    if k not in out:
                        out[k] = np.array([a.pos.x, a.pos.y, a.pos.z], dtype=np.float64)
        return out

    ref_atoms = _collect_atoms(st_ref)
    cur_atoms = _collect_atoms(st_cur)

    # Match C-alpha pairs by (chain, resseq); if chains were renamed, fall back to sequential CA order
    ca_keys = sorted(k for k in ref_atoms if k[2] == "CA" and k in cur_atoms)
    used_seq_fallback = False
    if len(ca_keys) < 4:
        ref_cas = [v for k, v in sorted(ref_atoms.items()) if k[2] == "CA"]
        cur_cas = [v for k, v in sorted(cur_atoms.items()) if k[2] == "CA"]
        n_min = min(len(ref_cas), len(cur_cas))
        if n_min < 3:
            raise ValueError("Fewer than 3 matched Cα atoms between reference and current structures")
        P_ca = np.array(ref_cas[:n_min], dtype=np.float64)
        Q_ca = np.array(cur_cas[:n_min], dtype=np.float64)
        used_seq_fallback = True
    else:
        P_ca = np.array([ref_atoms[k] for k in ca_keys], dtype=np.float64)
        Q_ca = np.array([cur_atoms[k] for k in ca_keys], dtype=np.float64)

    R, cP, cQ = _kabsch_align(P_ca, Q_ca)
    P_ca_al = (P_ca - cP) @ R + cQ
    ca_shifts = np.linalg.norm(P_ca_al - Q_ca, axis=1)
    ca_rmsd = round(float(np.sqrt(np.mean(ca_shifts ** 2))), 3)

    # Compute Backbone & All-Heavy-Atom RMSD on matched keys
    bb_diffs = []
    heavy_diffs = []
    res_heavy_shifts: Dict[Tuple[str, int], List[float]] = {}

    if not used_seq_fallback:
        for k, p_ref in ref_atoms.items():
            if k not in cur_atoms:
                continue
            p_al = (p_ref - cP) @ R + cQ
            d = float(np.linalg.norm(p_al - cur_atoms[k]))
            heavy_diffs.append(d * d)
            res_heavy_shifts.setdefault((k[0], k[1]), []).append(d)
            if k[2] in ("N", "CA", "C", "O"):
                bb_diffs.append(d * d)

    bb_rmsd = round(float(math.sqrt(sum(bb_diffs) / max(len(bb_diffs), 1))), 3) if bb_diffs else ca_rmsd
    heavy_rmsd = round(float(math.sqrt(sum(heavy_diffs) / max(len(heavy_diffs), 1))), 3) if heavy_diffs else ca_rmsd

    # Per-residue shift profile
    per_residue = []
    if not used_seq_fallback:
        res_names = {}
        for ch in st_cur[0]:
            for res in ch:
                res_names[(ch.name, res.seqid.num)] = res.name.strip()
        for idx, (ch_id, rseq, _) in enumerate(ca_keys):
            h_list = res_heavy_shifts.get((ch_id, rseq), [float(ca_shifts[idx])])
            per_residue.append({
                "chain": ch_id,
                "resseq": rseq,
                "resname": res_names.get((ch_id, rseq), "UNK"),
                "key": f"{ch_id}:{rseq}",
                "ca_shift_A": round(float(ca_shifts[idx]), 3),
                "max_atom_shift_A": round(float(max(h_list)), 3),
                "mean_heavy_shift_A": round(float(np.mean(h_list)), 3),
            })

    top_shifted = sorted(per_residue, key=lambda x: x["max_atom_shift_A"], reverse=True)[:15]

    # Build aligned copy of st_ref for translucent 3D ghost superposition!
    st_aligned = st_ref.clone()
    for m in st_aligned:
        for ch in m:
            for res in ch:
                for a in res:
                    p = np.array([a.pos.x, a.pos.y, a.pos.z], dtype=np.float64)
                    pa = (p - cP) @ R + cQ
                    a.pos = gemmi.Position(float(pa[0]), float(pa[1]), float(pa[2]))

    return {
        "ref_label": ref_label,
        "matched_ca": int(len(P_ca)),
        "matched_heavy": int(len(heavy_diffs)) if heavy_diffs else int(len(P_ca)),
        "ca_rmsd_A": ca_rmsd,
        "backbone_rmsd_A": bb_rmsd,
        "all_heavy_rmsd_A": heavy_rmsd,
        "max_shift_A": round(float(max((r["max_atom_shift_A"] for r in per_residue), default=0.0)), 3),
        "top_shifted": top_shifted,
        "per_residue": per_residue,
        "aligned_ref_pdb": structure_to_pdb(st_aligned),
    }


def analyze_membrane_and_interactions(st: gemmi.Structure) -> Dict[str, Any]:
    """
    Analyze non-covalent interaction networks (Salt Bridges, Pi-Pi Stacking, Inter-Chain Interface)
    and estimate transmembrane hydrophobic-belt orientation & bilayer planes.
    """
    if len(st) == 0:
        return {
            "salt_bridges": [],
            "pi_stacking": [],
            "interface_contacts": [],
            "membrane": None,
        }

    model = st[0]
    acidic_pts: List[Tuple[str, int, str, str, str, np.ndarray]] = []
    basic_pts: List[Tuple[str, int, str, str, str, np.ndarray]] = []
    ring_entries: List[Tuple[str, int, str, str, np.ndarray, np.ndarray]] = []  # (ch, resseq, ic, rn, centroid, normal)
    ca_pts: List[np.ndarray] = []
    ca_is_hydro: List[bool] = []
    chain_atoms: Dict[str, List[Tuple[int, str, str, str, np.ndarray]]] = {}

    for ch in model:
        for res in ch:
            rn = res.name.strip().upper()
            std = canonical_resname(rn)
            if rn not in ALL_PROTEIN_RESNAMES and rn not in NONSTANDARD_TO_STANDARD:
                continue
            ic = res.seqid.icode.strip()
            atom_map = {}
            for a in res:
                if a.element.name in ("H", "D"):
                    continue
                pos = np.array([a.pos.x, a.pos.y, a.pos.z], dtype=np.float64)
                an = a.name.strip()
                atom_map[an] = pos
                chain_atoms.setdefault(ch.name, []).append((res.seqid.num, ic, rn, an, pos))

            if "CA" in atom_map:
                ca_pts.append(atom_map["CA"])
                ca_is_hydro.append(std in HYDROPHOBIC_TM)

            # Acidic & basic charged atoms (exclude neutral ASH/GLH/LYN)
            if std in ACIDIC_ATOMS and rn not in ("ASH", "ASPH", "GLH", "GLUH"):
                for an in ACIDIC_ATOMS[std]:
                    if an in atom_map:
                        acidic_pts.append((ch.name, res.seqid.num, ic, rn, an, atom_map[an]))
            if std in BASIC_ATOMS and rn not in ("LYN", "LSN", "ARN"):
                # Only count HIS as charged basic if HIP/HSP or default HIS
                if std == "HIS" and rn in ("HID", "HIE", "HSD", "HSE"):
                    pass
                else:
                    for an in BASIC_ATOMS[std]:
                        if an in atom_map:
                            basic_pts.append((ch.name, res.seqid.num, ic, rn, an, atom_map[an]))

            # Aromatic ring centroids + plane normals
            if std in AROMATIC_RINGS:
                r_pts = [atom_map[an] for an in AROMATIC_RINGS[std] if an in atom_map]
                if len(r_pts) >= 5:
                    arr = np.array(r_pts, dtype=np.float64)
                    cent = np.mean(arr, axis=0)
                    centered = arr - cent
                    _, _, vh = np.linalg.svd(centered)
                    normal = vh[-1]
                    ring_entries.append((ch.name, res.seqid.num, ic, rn, cent, normal))

    # 1. Salt Bridges (<= 4.0 Å)
    salt_bridges = []
    seen_sb = set()
    for a_ch, a_num, a_ic, a_rn, a_an, a_pos in acidic_pts:
        for b_ch, b_num, b_ic, b_rn, b_an, b_pos in basic_pts:
            if a_ch == b_ch and abs(a_num - b_num) <= 1:
                continue
            d = distance(a_pos, b_pos)
            if d <= 4.0:
                pair_k = (a_ch, a_num, b_ch, b_num)
                if pair_k not in seen_sb:
                    seen_sb.add(pair_k)
                    salt_bridges.append({
                        "acid": f"{a_rn} {a_ch}:{a_num}{a_ic}",
                        "acid_chain": a_ch,
                        "acid_resseq": a_num,
                        "acid_atom": a_an,
                        "base": f"{b_rn} {b_ch}:{b_num}{b_ic}",
                        "base_chain": b_ch,
                        "base_resseq": b_num,
                        "base_atom": b_an,
                        "distance_A": round(d, 2),
                        "inter_chain": a_ch != b_ch,
                    })
    salt_bridges.sort(key=lambda x: x["distance_A"])

    # 2. Pi-Pi Aromatic Stacking (centroid distance <= 5.6 Å)
    pi_stacking = []
    for i in range(len(ring_entries)):
        for j in range(i + 1, len(ring_entries)):
            c1, n1, ic1, rn1, p1, norm1 = ring_entries[i]
            c2, n2, ic2, rn2, p2, norm2 = ring_entries[j]
            if c1 == c2 and abs(n1 - n2) <= 1:
                continue
            d = distance(p1, p2)
            if d <= 5.6:
                cos_ang = abs(float(np.clip(np.dot(norm1, norm2), -1.0, 1.0)))
                ang = math.degrees(math.acos(cos_ang))
                geom = "Parallel / Face-to-Face" if ang < 30.0 else ("T-Shaped / Edge-to-Face" if ang > 60.0 else "Tilted Offset")
                pi_stacking.append({
                    "res1": f"{rn1} {c1}:{n1}{ic1}",
                    "chain1": c1,
                    "resseq1": n1,
                    "res2": f"{rn2} {c2}:{n2}{ic2}",
                    "chain2": c2,
                    "resseq2": n2,
                    "distance_A": round(d, 2),
                    "angle_deg": round(ang, 1),
                    "geometry": geom,
                    "inter_chain": c1 != c2,
                })
    pi_stacking.sort(key=lambda x: x["distance_A"])

    # 3. Inter-Chain Interface Contacts (<= 4.2 Å across distinct protein chains)
    interface_contacts = []
    chain_ids = sorted(chain_atoms.keys())
    if len(chain_ids) >= 2:
        for ci in range(len(chain_ids)):
            c1 = chain_ids[ci]
            pts1 = np.array([x[4] for x in chain_atoms[c1]], dtype=np.float64)
            tree1 = cKDTree(pts1)
            for cj in range(ci + 1, len(chain_ids)):
                c2 = chain_ids[cj]
                pts2 = np.array([x[4] for x in chain_atoms[c2]], dtype=np.float64)
                seen_pair = set()
                pairs = tree1.query_ball_point(pts2, r=4.2)
                for idx2, hits1 in enumerate(pairs):
                    if not hits1:
                        continue
                    r2_num, r2_ic, r2_rn, a2_name, p2 = chain_atoms[c2][idx2]
                    for idx1 in hits1:
                        r1_num, r1_ic, r1_rn, a1_name, p1 = chain_atoms[c1][idx1]
                        pk = (c1, r1_num, c2, r2_num)
                        if pk in seen_pair:
                            continue
                        seen_pair.add(pk)
                        d = distance(p1, p2)
                        interface_contacts.append({
                            "res1": f"{r1_rn} {c1}:{r1_num}{r1_ic}",
                            "chain1": c1,
                            "resseq1": r1_num,
                            "res2": f"{r2_rn} {c2}:{r2_num}{r2_ic}",
                            "chain2": c2,
                            "resseq2": r2_num,
                            "distance_A": round(d, 2),
                        })
        interface_contacts.sort(key=lambda x: x["distance_A"])

    # 4. Transmembrane Hydrophobic-Belt Orientation Estimator (OPM-inspired)
    membrane_info = None
    if len(ca_pts) >= 20:
        ca_arr = np.array(ca_pts, dtype=np.float64)
        hydro_mask = np.array(ca_is_hydro, dtype=bool)
        center = np.mean(ca_arr, axis=0)
        centered = ca_arr - center
        # Evaluate candidate axis directions (PCA axes + Fibonacci sphere directions) for best 30 Å hydrophobic belt
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        cand_axes = [vh[0], vh[1], vh[2], np.array([0.0, 0.0, 1.0]), np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0])]
        # Add 20 evenly distributed Fibonacci directions
        golden = (1 + math.sqrt(5)) / 2
        for k in range(20):
            theta = 2 * math.pi * k / golden
            phi = math.acos(1 - 2 * (k + 0.5) / 20)
            cand_axes.append(np.array([math.sin(phi) * math.cos(theta), math.sin(phi) * math.sin(theta), math.cos(phi)]))

        best_score = -1e9
        best_axis = np.array([0.0, 0.0, 1.0])
        best_mid = 0.0
        best_hydro_belt_pct = 0.0
        best_hydro_out_pct = 0.0
        half_thick = 15.0  # Standard 30 Å hydrocarbon bilayer core

        # Estimate surface exposure for each Cα (lipids contact outer hydrophobic surface!)
        ca_tree = cKDTree(ca_arr)
        nbr_counts = np.array([len(ca_tree.query_ball_point(pt, r=10.0)) for pt in ca_arr], dtype=np.float64)
        # Surface weight: buried core residues get lower weight than lipid-facing surface residues
        surf_w = np.clip((26.0 - nbr_counts) / 14.0, 0.25, 1.2)

        for ax in cand_axes:
            n_ax = ax / (np.linalg.norm(ax) + 1e-8)
            proj = centered @ n_ax
            p_min, p_max = float(np.min(proj)), float(np.max(proj))
            # Slide the 30 Å membrane slab center along the axis to find the true TM belt
            # even when the protein has a large extracellular or cytoplasmic domain!
            if p_max - p_min > 36.0:
                mids = np.linspace(p_min + half_thick * 0.85, p_max - half_thick * 0.85, 9)
            else:
                mids = np.array([0.0])

            for mid in mids:
                in_slab = np.abs(proj - mid) <= half_thick
                out_slab = ~in_slab
                n_in = int(np.sum(in_slab))
                if n_in < 14:
                    continue
                w_in = surf_w[in_slab]
                h_in_w = float(np.sum(hydro_mask[in_slab] * w_in) / (np.sum(w_in) + 1e-8))
                h_in_raw = float(np.mean(hydro_mask[in_slab]))
                if np.sum(out_slab) >= 8:
                    w_out = surf_w[out_slab]
                    h_out_w = float(np.sum(hydro_mask[out_slab] * w_out) / (np.sum(w_out) + 1e-8))
                    h_out_raw = float(np.mean(hydro_mask[out_slab]))
                else:
                    h_out_w = 0.35
                    h_out_raw = 0.35
                # Reward high surface-exposed hydrophobic enrichment inside the 30 Å slab vs polar caps
                score = (h_in_w - h_out_w) * 110.0 + h_in_w * 30.0 + min(n_in, 180) * 0.04
                if score > best_score:
                    best_score = score
                    best_axis = n_ax
                    best_mid = float(mid)
                    best_hydro_belt_pct = round(h_in_raw * 100.0, 1)
                    best_hydro_out_pct = round(h_out_raw * 100.0, 1)

        # Ensure normal points generally toward +Z
        if best_axis[2] < 0:
            best_axis = -best_axis
            best_mid = -best_mid

        slab_center = center + best_axis * best_mid
        top_center = slab_center + best_axis * half_thick
        bot_center = slab_center - best_axis * half_thick
        # Compute radial extent of residues specifically inside the 30 Å membrane belt
        proj_best = (ca_arr - slab_center) @ best_axis
        in_best = np.abs(proj_best) <= half_thick + 3.0
        if np.any(in_best):
            perp = (ca_arr[in_best] - slab_center) - np.outer(proj_best[in_best], best_axis)
            belt_r = float(np.percentile(np.linalg.norm(perp, axis=1), 95))
        else:
            belt_r = float(np.max(np.linalg.norm(centered, axis=1))) * 0.65
        disc_radius = round(max(22.0, min(52.0, belt_r + 12.0)), 1)

        # Tilt angle relative to laboratory Z axis
        tilt_deg = round(math.degrees(math.acos(float(np.clip(best_axis[2], -1.0, 1.0)))), 1)
        likely_tm = best_hydro_belt_pct >= 52.0 and (best_hydro_belt_pct - best_hydro_out_pct) >= 8.0

        membrane_info = {
            "normal": [round(float(x), 4) for x in best_axis],
            "center": [round(float(x), 2) for x in slab_center],
            "top_leaflet_center": [round(float(x), 2) for x in top_center],
            "bottom_leaflet_center": [round(float(x), 2) for x in bot_center],
            "thickness_A": 30.0,
            "disc_radius_A": disc_radius,
            "tilt_from_z_deg": tilt_deg,
            "belt_hydrophobic_pct": best_hydro_belt_pct,
            "flank_hydrophobic_pct": best_hydro_out_pct,
            "classification": "Likely Transmembrane Domain" if likely_tm else "Soluble / Globular Protein",
        }

    return {
        "salt_bridges": salt_bridges[:80],
        "pi_stacking": pi_stacking[:80],
        "interface_contacts": interface_contacts[:80],
        "counts": {
            "salt_bridges": len(salt_bridges),
            "inter_chain_salt_bridges": sum(1 for s in salt_bridges if s["inter_chain"]),
            "pi_stacking": len(pi_stacking),
            "interface_contacts": len(interface_contacts),
        },
        "membrane": membrane_info,
    }


def orient_structure_along_z(st: gemmi.Structure, custom_axis: Optional[List[float]] = None) -> Dict[str, Any]:
    """
    Rotate & center the structure so its membrane normal (or pore axis) aligns cleanly
    with the laboratory +Z axis [0, 0, 1] and the bilayer slab sits centered at Z = 0.
    """
    info = analyze_membrane_and_interactions(st)
    mem = info.get("membrane")
    if custom_axis and len(custom_axis) == 3:
        n_vec = np.array(custom_axis, dtype=np.float64)
        n_vec = n_vec / (np.linalg.norm(n_vec) + 1e-8)
        center = np.array(mem["center"], dtype=np.float64) if mem else np.zeros(3)
    elif mem is not None:
        n_vec = np.array(mem["normal"], dtype=np.float64)
        center = np.array(mem["center"], dtype=np.float64)
    else:
        raise ValueError("Structure has too few Cα atoms to orient")

    z_axis = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    v = np.cross(n_vec, z_axis)
    s = float(np.linalg.norm(v))
    c = float(np.dot(n_vec, z_axis))

    if s < 1e-6:
        R = np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], dtype=np.float64)
        R = np.eye(3) + vx + (vx @ vx) * ((1.0 - c) / (s * s))

    for m in st:
        for ch in m:
            for res in ch:
                for a in res:
                    p = np.array([a.pos.x, a.pos.y, a.pos.z], dtype=np.float64) - center
                    pr = R @ p
                    a.pos = gemmi.Position(float(pr[0]), float(pr[1]), float(pr[2]))

    st.setup_entities()
    return {
        "changed": True,
        "tilt_corrected_deg": mem["tilt_from_z_deg"] if mem else 0.0,
        "message": f"Oriented channel/membrane axis along +Z and centered at (0, 0, 0) (corrected {mem['tilt_from_z_deg'] if mem else 0.0}° tilt)",
    }
