"""
HOLE Ion-Channel Pore Radius & Conductance Engine (Smart et al., J. Mol. Graph. 1996).

Computes the exact maximum-inscribed-sphere pore radius profile R(s) along a channel axis
using vectorized Monte Carlo + local Nelder-Mead optimization on perpendicular slices,
with classic AMBER/Bondi van der Waals radii, Smart et al. resistance/conductance integration,
red/green/blue water-accessibility coloring (<1.15 Å / 1.15–2.30 Å / >2.30 Å),
pore-lining residue identification, and CSV / PyMOL PDB exports.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import gemmi
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from .analysis_engine import analyze_membrane_and_interactions
from .chemistry_data import ALL_PROTEIN_RESNAMES, NONSTANDARD_TO_STANDARD
from .geometry_utils import canonical_resname

# Classic HOLE (Smart et al. 1996, simple.rad / AMBER) van der Waals radii (Å)
HOLE_VDW_RADII: Dict[str, float] = {
    "C": 1.85,
    "N": 1.75,
    "O": 1.65,
    "S": 2.00,
    "P": 2.10,
    "H": 1.00,
    "D": 1.00,
    "F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "ZN": 1.39,
    "MG": 1.73,
    "CA": 1.97,
    "NA": 2.27,
    "K": 2.75,
}


def _extract_hole_atoms(
    st: gemmi.Structure,
    include_hetatms: bool = False,
) -> Tuple[np.ndarray, np.ndarray, List[Dict[str, Any]]]:
    """Extract atomic coordinates, HOLE vdW radii, and residue metadata."""
    pts: List[List[float]] = []
    rads: List[float] = []
    meta: List[Dict[str, Any]] = []
    if len(st) == 0:
        return np.zeros((0, 3)), np.zeros((0,)), []

    water_names = {"HOH", "WAT", "SOL", "DOD"}
    for ch in st[0]:
        for res in ch:
            rn = res.name.strip().upper()
            if rn in water_names:
                continue
            is_prot = rn in ALL_PROTEIN_RESNAMES or rn in NONSTANDARD_TO_STANDARD
            if not is_prot and not include_hetatms:
                continue
            std = canonical_resname(rn)
            ic = res.seqid.icode.strip()
            key = f"{ch.name}:{res.seqid.num}{ic}"
            for a in res:
                elem = a.element.name.strip().upper()
                pts.append([a.pos.x, a.pos.y, a.pos.z])
                rads.append(HOLE_VDW_RADII.get(elem, 1.80))
                meta.append({
                    "chain": ch.name,
                    "resseq": res.seqid.num,
                    "icode": ic,
                    "resname": rn,
                    "std_resname": std,
                    "atom": a.name.strip(),
                    "key": key,
                    "is_backbone": a.name.strip() in ("N", "CA", "C", "O", "H", "HA"),
                })

    if not pts:
        return np.zeros((0, 3)), np.zeros((0,)), []
    return np.array(pts, dtype=np.float64), np.array(rads, dtype=np.float64), meta


def _detect_channel_axis_and_center(
    st: gemmi.Structure,
    axis_mode: str = "auto",
    custom_point: Optional[List[float]] = None,
) -> Tuple[np.ndarray, np.ndarray, str]:
    """
    Determine channel axis vector (cvect) and pore seed center (cpoint).
    For symmetric oligomeric channels (>=3 homologous-sized chains), detects the exact
    symmetry axis perpendicular to the subunit ring; otherwise uses the transmembrane
    hydrophobic-belt normal or explicit X/Y/Z axis.
    """
    model = st[0]
    chain_cas: Dict[str, np.ndarray] = {}
    all_cas: List[np.ndarray] = []
    for ch in model:
        c_pts = []
        for res in ch:
            rn = res.name.strip().upper()
            if rn in ALL_PROTEIN_RESNAMES or rn in NONSTANDARD_TO_STANDARD:
                ca = res.find_atom("CA", "\0")
                if ca is not None:
                    p = np.array([ca.pos.x, ca.pos.y, ca.pos.z], dtype=np.float64)
                    c_pts.append(p)
                    all_cas.append(p)
        if len(c_pts) >= 20:
            chain_cas[ch.name] = np.array(c_pts, dtype=np.float64)

    all_ca_arr = np.array(all_cas, dtype=np.float64) if all_cas else np.zeros((1, 3))
    global_center = np.mean(all_ca_arr, axis=0)

    # Explicit lab axes
    mode = axis_mode.strip().lower()
    if mode in ("z", "y", "x"):
        axis_map = {
            "z": np.array([0.0, 0.0, 1.0]),
            "y": np.array([0.0, 1.0, 0.0]),
            "x": np.array([1.0, 0.0, 0.0]),
        }
        cvect = axis_map[mode]
        cpoint = np.array(custom_point, dtype=np.float64) if (custom_point and len(custom_point) == 3) else global_center
        return cvect, cpoint, f"Manual +{mode.upper()} Axis"

    # Auto mode: check for radial oligomeric channel (>= 3 major chains arranged around central pore)
    if len(chain_cas) >= 3:
        sizes = sorted([len(v) for v in chain_cas.values()], reverse=True)
        # Keep chains within 40% of the largest subunit (ignore tiny peptides/antibodies if present)
        major = [pts for pts in chain_cas.values() if len(pts) >= sizes[0] * 0.45]
        if len(major) >= 3:
            centroids = np.array([np.mean(pts, axis=0) for pts in major], dtype=np.float64)
            ring_center = np.mean(centroids, axis=0)
            rel = centroids - ring_center
            _, svals, vh = np.linalg.svd(rel, full_matrices=False)
            # If centroids are roughly coplanar (ring of subunits), vh[-1] is the exact pore symmetry axis!
            if len(svals) == 3 and svals[0] > 4.0 and (svals[2] / (svals[0] + 1e-6)) < 0.45:
                cvect = vh[-1] / (np.linalg.norm(vh[-1]) + 1e-8)
                if cvect[2] < 0:
                    cvect = -cvect
                cpoint = np.array(custom_point, dtype=np.float64) if (custom_point and len(custom_point) == 3) else ring_center
                return cvect, cpoint, f"Oligomeric {len(major)}-Subunit Symmetry Axis"

    # Fallback: Transmembrane Hydrophobic-Belt normal from analysis_engine
    mem = analyze_membrane_and_interactions(st).get("membrane")
    if mem is not None:
        cvect = np.array(mem["normal"], dtype=np.float64)
        cvect = cvect / (np.linalg.norm(cvect) + 1e-8)
        cpoint = (
            np.array(custom_point, dtype=np.float64)
            if (custom_point and len(custom_point) == 3)
            else np.array(mem["center"], dtype=np.float64)
        )
        return cvect, cpoint, "Transmembrane Hydrophobic-Belt Axis"

    cvect = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    cpoint = np.array(custom_point, dtype=np.float64) if (custom_point and len(custom_point) == 3) else global_center
    return cvect, cpoint, "Principal +Z Axis"


def _orthonormal_basis(n: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    ref = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.85 else np.array([0.0, 1.0, 0.0])
    u = np.cross(ref, n)
    u = u / (np.linalg.norm(u) + 1e-12)
    v = np.cross(n, u)
    v = v / (np.linalg.norm(v) + 1e-12)
    return u, v


def _max_sphere_on_plane(
    plane_origin: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    init_uv: Tuple[float, float],
    atom_pts: np.ndarray,
    atom_rads: np.ndarray,
    atom_tree: cKDTree,
    search_radius: float = 5.5,
    max_drift_from_init: float = 3.5,
    end_radius: float = 10.0,
) -> Tuple[float, Tuple[float, float], np.ndarray, int]:
    """
    HOLE slice optimization: finds (u, v) on the perpendicular slice plane maximizing
    the largest inscribed sphere clearance R = min_i (||p - r_i|| - vdw_i).
    Uses fast coarse polar search (Monte Carlo-like multi-start) followed by Nelder-Mead polish.
    """
    seed_pt = plane_origin + init_uv[0] * u + init_uv[1] * v
    nbr_idx = atom_tree.query_ball_point(seed_pt, r=end_radius + search_radius + 4.5)
    if len(nbr_idx) < 4:
        return end_radius + 1.0, init_uv, seed_pt, -1

    sub_pts = atom_pts[nbr_idx]
    sub_rads = atom_rads[nbr_idx]
    nbr_arr = np.array(nbr_idx, dtype=np.int64)

    def eval_uv(uv: np.ndarray) -> Tuple[float, int]:
        p = plane_origin + uv[0] * u + uv[1] * v
        dists = np.linalg.norm(sub_pts - p, axis=1) - sub_rads
        min_i = int(np.argmin(dists))
        return float(dists[min_i]), int(nbr_arr[min_i])

    # Step 1: Vectorized 2D grid/polar multi-start around init_uv (fast & escapes local traps)
    cand_offsets = [(0.0, 0.0)]
    for rad in (0.6, 1.4, 2.4, min(3.4, max_drift_from_init)):
        n_ang = 8 if rad < 2.0 else 12
        for k in range(n_ang):
            ang = 2.0 * math.pi * k / n_ang
            cand_offsets.append((rad * math.cos(ang), rad * math.sin(ang)))
    cands_uv = np.array([[init_uv[0] + du, init_uv[1] + dv] for du, dv in cand_offsets], dtype=np.float64)
    cands_3d = plane_origin[None, :] + cands_uv[:, 0:1] * u[None, :] + cands_uv[:, 1:2] * v[None, :]

    # Vectorized clearance matrix: (M_cands, N_atoms)
    diff = cands_3d[:, None, :] - sub_pts[None, :, :]
    dist_mat = np.sqrt(np.sum(diff * diff, axis=2)) - sub_rads[None, :]
    r_cands = np.min(dist_mat, axis=1)
    best_k = int(np.argmax(r_cands))
    best_uv = cands_uv[best_k].copy()
    best_r = float(r_cands[best_k])

    # Step 2: Nelder-Mead local refinement constrained within max_drift_from_init
    def objective(x: np.ndarray) -> float:
        du = x[0] - init_uv[0]
        dv = x[1] - init_uv[1]
        drift = math.hypot(du, dv)
        p = plane_origin + x[0] * u + x[1] * v
        r_val = float(np.min(np.linalg.norm(sub_pts - p, axis=1) - sub_rads))
        if drift > max_drift_from_init:
            r_val -= (drift - max_drift_from_init) * 3.0
        return -r_val

    res = minimize(
        objective,
        best_uv,
        method="Nelder-Mead",
        options={"xatol": 0.04, "fatol": 0.015, "maxiter": 60},
    )
    opt_uv = res.x
    if math.hypot(opt_uv[0] - init_uv[0], opt_uv[1] - init_uv[1]) <= max_drift_from_init + 0.25:
        r_opt, closest_idx = eval_uv(opt_uv)
        if r_opt >= best_r:
            best_r = r_opt
            best_uv = opt_uv
        else:
            best_r, closest_idx = eval_uv(best_uv)
    else:
        best_r, closest_idx = eval_uv(best_uv)

    best_pt = plane_origin + best_uv[0] * u + best_uv[1] * v
    return best_r, (float(best_uv[0]), float(best_uv[1])), best_pt, closest_idx


def _pore_color(radius: float) -> Tuple[str, str]:
    """Classic HOLE (Smart et al. 1996) 3-zone water/ion permeability coloring."""
    if radius < 1.15:
        return "#FF3B30", "closed"      # Red: < 1 water molecule (gate / constriction)
    if radius <= 2.30:
        return "#34C759", "single_file" # Green: single-file water / selectivity filter
    return "#007AFF", "hydrated"        # Blue: wide hydrated vestibule


def analyze_hole_pore(
    st: gemmi.Structure,
    axis_mode: str = "auto",
    step: float = 0.5,
    end_radius: float = 10.0,
    custom_point: Optional[List[float]] = None,
    include_hetatms: bool = False,
) -> Dict[str, Any]:
    """
    Run HOLE-algorithm ion-channel pore radius analysis along the channel axis.
    Walks bidirectionally (+s and -s) from the pore center until reaching bulk solvent
    (radius > end_radius) or the protein boundary.
    """
    step = float(max(0.25, min(2.0, step)))
    end_radius = float(max(4.5, min(18.0, end_radius)))

    atom_pts, atom_rads, atom_meta = _extract_hole_atoms(st, include_hetatms=include_hetatms)
    if len(atom_pts) < 30:
        return {
            "axis_label": "Unavailable",
            "cvect": [0.0, 0.0, 1.0],
            "cpoint": [0.0, 0.0, 0.0],
            "profile": [],
            "lining_residues": [],
            "summary": {
                "length_A": 0.0,
                "min_radius_A": 0.0,
                "max_radius_A": 0.0,
                "mean_radius_A": 0.0,
                "constriction_coord_A": 0.0,
                "constriction_point": [0.0, 0.0, 0.0],
                "constriction_residues": [],
                "gate_status": "No pore detected",
                "conductance_pS": 0.0,
                "pore_resistance_GOhm": 0.0,
                "red_fraction_pct": 0.0,
                "green_fraction_pct": 0.0,
                "blue_fraction_pct": 0.0,
            },
        }

    cvect, cpoint, axis_label = _detect_channel_axis_and_center(st, axis_mode=axis_mode, custom_point=custom_point)
    u, v = _orthonormal_basis(cvect)
    atom_tree = cKDTree(atom_pts)

    # Determine total protein span along cvect relative to cpoint
    projs = (atom_pts - cpoint) @ cvect
    s_min_bound = float(np.min(projs)) - 2.0
    s_max_bound = float(np.max(projs)) + 2.0

    # Find best starting center on the s=0 plane (wider search radius initially so even slightly off-center seeds snap into the pore)
    r0, uv0, pt0, cidx0 = _max_sphere_on_plane(
        cpoint, u, v, (0.0, 0.0), atom_pts, atom_rads, atom_tree,
        search_radius=7.0, max_drift_from_init=6.0, end_radius=end_radius,
    )

    def walk_direction(sign: float) -> List[Tuple[float, float, np.ndarray, int]]:
        collected: List[Tuple[float, float, np.ndarray, int]] = []
        cur_uv = uv0
        s = sign * step
        wide_streak = 0
        while (s <= s_max_bound if sign > 0 else s >= s_min_bound):
            plane_org = cpoint + s * cvect
            r_s, next_uv, pt_s, c_idx = _max_sphere_on_plane(
                plane_org, u, v, cur_uv, atom_pts, atom_rads, atom_tree,
                search_radius=5.0, max_drift_from_init=2.2, end_radius=end_radius,
            )
            if r_s >= end_radius:
                wide_streak += 1
                if wide_streak >= 2:
                    break
            else:
                wide_streak = 0
                collected.append((s, r_s, pt_s, c_idx))
                cur_uv = next_uv
            s += sign * step
        return collected

    neg_pts = walk_direction(-1.0)
    neg_pts.reverse()
    pos_pts = walk_direction(+1.0)

    raw_slices: List[Tuple[float, float, np.ndarray, int]] = []
    raw_slices.extend(neg_pts)
    if r0 < end_radius:
        raw_slices.append((0.0, r0, pt0, cidx0))
    raw_slices.extend(pos_pts)

    if not raw_slices:
        return {
            "axis_label": axis_label,
            "cvect": [round(float(x), 4) for x in cvect],
            "cpoint": [round(float(x), 2) for x in cpoint],
            "profile": [],
            "lining_residues": [],
            "summary": {
                "length_A": 0.0,
                "min_radius_A": 0.0,
                "max_radius_A": 0.0,
                "mean_radius_A": 0.0,
                "constriction_coord_A": 0.0,
                "constriction_point": [round(float(x), 2) for x in cpoint],
                "constriction_residues": [],
                "gate_status": "No enclosed pore along this axis",
                "conductance_pS": 0.0,
                "pore_resistance_GOhm": 0.0,
                "red_fraction_pct": 0.0,
                "green_fraction_pct": 0.0,
                "blue_fraction_pct": 0.0,
            },
        }

    # Build full profile + lining residue inventory
    profile: List[Dict[str, Any]] = []
    res_stats: Dict[str, Dict[str, Any]] = {}

    for s_val, r_raw, pt, c_idx in raw_slices:
        r_clamped = max(0.0, r_raw)
        col_hex, state = _pore_color(r_clamped)

        # Find all residues lining this pore slice (within r_clamped + vdw + 2.0 Å)
        nbrs = atom_tree.query_ball_point(pt, r=max(3.5, r_clamped + 4.2))
        slice_res_keys: List[str] = []
        seen_sr = set()
        for ni in nbrs:
            m = atom_meta[ni]
            clr = float(np.linalg.norm(atom_pts[ni] - pt) - atom_rads[ni])
            if clr <= r_clamped + 1.8:
                rk = m["key"]
                lbl = f"{m['resname']} {rk}"
                if lbl not in seen_sr:
                    seen_sr.add(lbl)
                    slice_res_keys.append(lbl)
                if rk not in res_stats:
                    res_stats[rk] = {
                        "chain": m["chain"],
                        "resseq": m["resseq"],
                        "icode": m["icode"],
                        "resname": m["resname"],
                        "key": rk,
                        "min_pore_radius_A": r_clamped,
                        "axis_min_A": s_val,
                        "axis_max_A": s_val,
                        "closest_atom": m["atom"],
                        "min_clearance_A": clr,
                        "backbone_contact": m["is_backbone"],
                    }
                else:
                    st_r = res_stats[rk]
                    if r_clamped < st_r["min_pore_radius_A"]:
                        st_r["min_pore_radius_A"] = r_clamped
                    if clr < st_r["min_clearance_A"]:
                        st_r["min_clearance_A"] = clr
                        st_r["closest_atom"] = m["atom"]
                    st_r["axis_min_A"] = min(st_r["axis_min_A"], s_val)
                    st_r["axis_max_A"] = max(st_r["axis_max_A"], s_val)
                    st_r["backbone_contact"] = st_r["backbone_contact"] or m["is_backbone"]

        closest_m = atom_meta[c_idx] if (0 <= c_idx < len(atom_meta)) else None
        profile.append({
            "axis_coord_A": round(float(s_val), 2),
            "radius_A": round(float(r_clamped), 2),
            "x": round(float(pt[0]), 2),
            "y": round(float(pt[1]), 2),
            "z": round(float(pt[2]), 2),
            "color": col_hex,
            "state": state,
            "closest_residue": f"{closest_m['resname']} {closest_m['key']}" if closest_m else "",
            "closest_atom": closest_m["atom"] if closest_m else "",
            "lining_residues": slice_res_keys[:8],
        })

    # Constriction gate (minimum pore radius along the channel)
    min_pt = min(profile, key=lambda p: p["radius_A"])
    min_r = min_pt["radius_A"]
    max_r = max(p["radius_A"] for p in profile)
    mean_r = sum(p["radius_A"] for p in profile) / len(profile)
    pore_len = round(profile[-1]["axis_coord_A"] - profile[0]["axis_coord_A"], 2)

    # Identify constriction-gate residues (within 2.5 Å of the narrowest constriction)
    constriction_res = []
    c_s = min_pt["axis_coord_A"]
    lining_list = []
    for rinfo in res_stats.values():
        is_gate = (rinfo["axis_min_A"] - 1.5 <= c_s <= rinfo["axis_max_A"] + 1.5) and (
            rinfo["min_pore_radius_A"] <= min_r + 0.85
        )
        r_entry = {
            "chain": rinfo["chain"],
            "resseq": rinfo["resseq"],
            "icode": rinfo["icode"],
            "resname": rinfo["resname"],
            "key": rinfo["key"],
            "min_pore_radius_A": round(rinfo["min_pore_radius_A"], 2),
            "axis_min_A": round(rinfo["axis_min_A"], 1),
            "axis_max_A": round(rinfo["axis_max_A"], 1),
            "closest_atom": rinfo["closest_atom"],
            "backbone_contact": rinfo["backbone_contact"],
            "is_constriction": is_gate,
        }
        lining_list.append(r_entry)
        if is_gate:
            constriction_res.append(f"{rinfo['resname']} {rinfo['key']}")

    lining_list.sort(key=lambda x: (not x["is_constriction"], x["min_pore_radius_A"], x["axis_min_A"]))

    # Classic HOLE Conductance Estimate (Smart et al. 1997, ~1M KCl rho = 0.12 Ohm·m)
    # R_pore = rho * sum(ds / (pi * r^2)) + access resistance; note 1 Å = 1e-10 m -> ds/(pi*r^2) has units 1e10 m^-1
    rho = 0.12  # Ohm·m (physiological ~1M KCl standard HOLE reference; 0.65 for 150 mM)
    if min_r < 0.85:
        conductance_pS = 0.0
        resistance_GOhm = 999.0
        gate_status = "Closed Gate (< 1.15 Å — Water-Impermeable)"
    else:
        integral = sum(step / (math.pi * max(p["radius_A"], 0.5) ** 2) for p in profile)
        r_end1 = max(profile[0]["radius_A"], 1.0)
        r_end2 = max(profile[-1]["radius_A"], 1.0)
        r_ohm = rho * 1e10 * (integral + 0.25 / r_end1 + 0.25 / r_end2)
        # Empirical HOLE correction factor (~5.0x for non-bulk nano-confined ion mobility, Smart et al. 1997)
        r_ohm_corr = r_ohm * 5.0
        resistance_GOhm = round(r_ohm_corr / 1e9, 2)
        conductance_pS = round((1.0 / r_ohm_corr) * 1e12, 1) if r_ohm_corr > 0 else 0.0
        if min_r < 1.15:
            gate_status = "Tight Constriction (< 1.15 Å — Dehydrated / Closed)"
        elif min_r <= 2.30:
            gate_status = "Single-File Water / Selectivity Filter (1.15–2.30 Å)"
        else:
            gate_status = "Open Hydrated Channel (> 2.30 Å)"

    n_tot = len(profile)
    red_pct = round(100.0 * sum(1 for p in profile if p["state"] == "closed") / n_tot, 1)
    green_pct = round(100.0 * sum(1 for p in profile if p["state"] == "single_file") / n_tot, 1)
    blue_pct = round(100.0 * sum(1 for p in profile if p["state"] == "hydrated") / n_tot, 1)

    return {
        "axis_label": axis_label,
        "cvect": [round(float(x), 4) for x in cvect],
        "cpoint": [round(float(x), 2) for x in cpoint],
        "params": {
            "axis_mode": axis_mode,
            "step": step,
            "end_radius": end_radius,
            "include_hetatms": include_hetatms,
        },
        "profile": profile,
        "lining_residues": lining_list[:100],
        "summary": {
            "length_A": pore_len,
            "min_radius_A": round(min_r, 2),
            "max_radius_A": round(max_r, 2),
            "mean_radius_A": round(mean_r, 2),
            "constriction_coord_A": min_pt["axis_coord_A"],
            "constriction_point": [min_pt["x"], min_pt["y"], min_pt["z"]],
            "constriction_residues": constriction_res[:12],
            "gate_status": gate_status,
            "conductance_pS": conductance_pS,
            "pore_resistance_GOhm": resistance_GOhm,
            "red_fraction_pct": red_pct,
            "green_fraction_pct": green_pct,
            "blue_fraction_pct": blue_pct,
        },
    }


def export_hole_csv(hole_data: Dict[str, Any]) -> str:
    """Export HOLE pore radius profile with per-slice closest & lining residues as CSV."""
    lines = [
        "axis_coord_A,pore_radius_A,x,y,z,water_accessibility,closest_residue,closest_atom,lining_residues"
    ]
    for pt in hole_data.get("profile", []):
        lining = "; ".join(pt.get("lining_residues", []))
        lines.append(
            f"{pt['axis_coord_A']},{pt['radius_A']},{pt['x']},{pt['y']},{pt['z']},"
            f"{pt['state']},\"{pt['closest_residue']}\",{pt['closest_atom']},\"{lining}\""
        )
    return "\n".join(lines) + "\n"


def export_hole_pdb(hole_data: Dict[str, Any]) -> str:
    """Export HOLE pore spheres as a PyMOL/ChimeraX/VMD-ready PDB (radius in B-factor, zone in chain R/G/B)."""
    lines = [
        "REMARK   1 Generated by PandaProtein HOLE Pore Radius Engine (Smart et al. 1996)",
        f"REMARK   2 Axis: {hole_data.get('axis_label', '')}  MinRadius: {hole_data.get('summary', {}).get('min_radius_A', 0)} A",
        "REMARK   3 Chain R (<1.15 A closed) | Chain G (1.15-2.30 A single-file) | Chain B (>2.30 A hydrated)",
        "REMARK   4 B-factor column stores exact HOLE pore radius in Angstroms (set sphere_scale in PyMOL).",
    ]
    serial = 1
    for idx, pt in enumerate(hole_data.get("profile", []), start=1):
        ch = "R" if pt["state"] == "closed" else ("G" if pt["state"] == "single_file" else "B")
        lines.append(
            f"HETATM{serial:5d}  SPH HOL {ch}{idx:4d}    "
            f"{pt['x']:8.3f}{pt['y']:8.3f}{pt['z']:8.3f}  1.00{pt['radius_A']:6.2f}           H"
        )
        serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"
