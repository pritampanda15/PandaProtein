"""
MOLE 2.5 / CAVER-Inspired Tunnels, Channels & Pores Engine:
- 3D Delaunay triangulation / Voronoi medial-axis void-space skeleton (SciPy)
- Exact clearance radius calculation with Bondi van der Waals radii
- Weighted Dijkstra pathfinding (cost = L / (r + eps)^2) with tunnel diversification
- Supports:
    1) Tunnels: from active-site / ligand / metal / pocket / custom residue origin to bulk solvent surface
    2) Pores / Channels: through-protein pathways connecting opposite solvent openings
- Computes smooth centerline profiles (radius, bottleneck, tortuosity, throughput)
- Computes MOLE 2.5 physicochemical profiles (Kyte-Doolittle hydropathy, Cid hydrophobicity,
  Zimmerman polarity, net charge, Jones mutability) and full lining + bottleneck residue lists
- Exports tunnel geometry to JSON, CSV, and PyMOL/ChimeraX-compatible pseudo-atom PDB
"""

from __future__ import annotations

import heapq
import math
from typing import Any, Dict, List, Optional, Tuple

import gemmi
import numpy as np
from scipy.spatial import Delaunay, cKDTree

from .chemistry_data import (
    ALL_PROTEIN_RESNAMES,
    NONSTANDARD_TO_STANDARD,
    VDW_RADII,
    WATER_RESNAMES,
)
from .geometry_utils import canonical_resname
from .pocket_engine import detect_binding_pockets

# MOLE 2.5 standard residue physicochemical scales
# Kyte-Doolittle hydropathy (-4.5 hydrophilic .. +4.5 hydrophobic)
KYTE_DOOLITTLE: Dict[str, float] = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
    "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
    "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
    "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}

# Zimmerman polarity (0 nonpolar .. 52 highly polar)
POLARITY_SCALE: Dict[str, float] = {
    "ARG": 52.0, "ASP": 49.7, "GLU": 49.9, "HIS": 51.6, "LYS": 49.5,
    "ASN": 3.38, "GLN": 3.53, "SER": 1.67, "THR": 1.66, "TYR": 1.61,
    "TRP": 2.10, "CYS": 1.48, "PRO": 1.58, "MET": 1.43, "PHE": 0.35,
    "VAL": 0.13, "LEU": 0.13, "ILE": 0.13, "ALA": 0.0, "GLY": 0.0,
}

# Formal charge at physiological pH
CHARGE_SCALE: Dict[str, int] = {
    "ARG": 1, "LYS": 1, "HIP": 1,
    "ASP": -1, "GLU": -1, "CYM": -1, "TYM": -1,
}

# Jones relative mutability (MOLE 2.5 mutability metric, ~18-134)
MUTABILITY_SCALE: Dict[str, float] = {
    "ASN": 134, "SER": 120, "ASP": 106, "GLU": 102, "ALA": 100,
    "THR": 97, "ILE": 96, "MET": 94, "GLN": 93, "VAL": 74,
    "HIS": 66, "ARG": 65, "LYS": 56, "PRO": 56, "GLY": 49,
    "PHE": 41, "TYR": 41, "LEU": 40, "CYS": 20, "TRP": 18,
}

TUNNEL_PALETTE = [
    "#007AFF",  # Apple Blue
    "#34C759",  # Apple Green
    "#FF9500",  # Apple Orange
    "#AF52DE",  # Apple Purple
    "#5AC8FA",  # Apple Teal
    "#FF2D55",  # Apple Pink
    "#5856D6",  # Apple Indigo
    "#FFCC00",  # Apple Yellow
]


def _circumcenter_batch(pts: np.ndarray, simplices: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Vectorized exact circumcenter & circumradius for all 3D Delaunay tetrahedra.
    Returns (centers [M, 3], radii [M]) with NaN for degenerate coplanar simplices.
    """
    p0 = pts[simplices[:, 0]]
    p1 = pts[simplices[:, 1]]
    p2 = pts[simplices[:, 2]]
    p3 = pts[simplices[:, 3]]

    a = p1 - p0
    b = p2 - p0
    c = p3 - p0

    # Solve 3x3 linear system: 2 * [a; b; c] @ u = [|a|^2, |b|^2, |c|^2]
    A = np.stack([a, b, c], axis=1)  # (M, 3, 3)
    rhs = 0.5 * np.stack([
        np.sum(a * a, axis=1),
        np.sum(b * b, axis=1),
        np.sum(c * c, axis=1),
    ], axis=1)  # (M, 3)

    dets = np.linalg.det(A)
    valid = np.abs(dets) > 1e-4

    centers = np.full((len(simplices), 3), np.nan, dtype=np.float64)
    radii = np.full(len(simplices), np.nan, dtype=np.float64)

    if np.any(valid):
        u = np.linalg.solve(A[valid], rhs[valid][..., None])[..., 0]
        centers[valid] = p0[valid] + u
        radii[valid] = np.linalg.norm(u, axis=1)

    return centers, radii


def _extract_atoms(
    st: gemmi.Structure, include_hetatms: bool = False
) -> Tuple[np.ndarray, np.ndarray, List[Dict[str, Any]]]:
    """Extract heavy atoms (coords, vdW radii, metadata) for Voronoi channel computation."""
    coords: List[List[float]] = []
    radii: List[float] = []
    meta: List[Dict[str, Any]] = []

    if len(st) == 0:
        return np.zeros((0, 3)), np.zeros((0,)), []

    for ch in st[0]:
        for res in ch:
            rn = res.name.strip().upper()
            if rn in WATER_RESNAMES:
                continue
            is_prot = rn in ALL_PROTEIN_RESNAMES or rn in NONSTANDARD_TO_STANDARD
            if not is_prot and not include_hetatms:
                continue
            ic = res.seqid.icode.strip()
            key = f"{ch.name}:{res.seqid.num}{ic}"
            for a in res:
                if a.element.name in ("H", "D"):
                    continue
                el = a.element.name.upper()
                r_vdw = VDW_RADII.get(el, 1.70)
                coords.append([a.pos.x, a.pos.y, a.pos.z])
                radii.append(r_vdw)
                meta.append({
                    "chain": ch.name,
                    "resseq": res.seqid.num,
                    "icode": ic,
                    "resname": rn,
                    "std_resname": canonical_resname(rn),
                    "key": key,
                    "atom": a.name.strip(),
                    "is_backbone": a.name.strip() in ("N", "CA", "C", "O", "OXT"),
                    "is_protein": is_prot,
                })

    if not coords:
        return np.zeros((0, 3)), np.zeros((0,)), []
    return np.array(coords, dtype=np.float64), np.array(radii, dtype=np.float64), meta


def _resolve_start_points(
    st: gemmi.Structure,
    start_mode: str,
    start_residues: Optional[List[str]] = None,
    start_point: Optional[List[float]] = None,
) -> List[Tuple[str, np.ndarray]]:
    """Resolve one or more starting origins (label, 3D point) for tunnel calculation."""
    if start_mode == "coordinates" and start_point and len(start_point) == 3:
        return [("Custom Point", np.array(start_point, dtype=np.float64))]

    if start_mode == "residues" and start_residues:
        # Parse specs like "A:41", "A:145", or "41"
        wanted = set()
        for tok in start_residues:
            t = tok.strip().upper()
            if not t:
                continue
            if ":" in t:
                ch_p, num_p = t.split(":", 1)
                wanted.add((ch_p.strip(), num_p.strip()))
            else:
                wanted.add(("", t))
        pts = []
        if len(st) > 0:
            for ch in st[0]:
                for res in ch:
                    num_s = f"{res.seqid.num}{res.seqid.icode.strip()}".upper()
                    if (ch.name.upper(), num_s) in wanted or ("", num_s) in wanted:
                        for a in res:
                            if a.element.name not in ("H", "D"):
                                pts.append([a.pos.x, a.pos.y, a.pos.z])
        if pts:
            lbl = ", ".join(start_residues[:3])
            return [(f"Residues ({lbl})", np.mean(np.array(pts), axis=0))]

    origins: List[Tuple[str, np.ndarray]] = []

    # Auto mode: prioritize co-crystallized ligand / cofactor / catalytic metal centroids, then deepest pockets
    if start_mode in ("auto", "ligand") and len(st) > 0:
        for ch in st[0]:
            for res in ch:
                rn = res.name.strip().upper()
                if rn in ALL_PROTEIN_RESNAMES or rn in NONSTANDARD_TO_STANDARD or rn in WATER_RESNAMES:
                    continue
                # Skip tiny crystallization ions/buffers if bigger ligands exist
                h_pts = [[a.pos.x, a.pos.y, a.pos.z] for a in res if a.element.name not in ("H", "D")]
                if len(h_pts) >= 5:
                    origins.append((f"Ligand {rn} ({ch.name}:{res.seqid.num})", np.mean(np.array(h_pts), axis=0)))

    if start_mode in ("auto", "pocket") or not origins:
        pockets = detect_binding_pockets(st, max_pockets=3)
        for pk in pockets:
            origins.append((f"Cavity #{pk['rank']} ({pk['volume_A3']} Å³)", np.array(pk["center"], dtype=np.float64)))

    if not origins and len(st) > 0:
        # Fallback to protein geometric centroid
        all_p = [
            [a.pos.x, a.pos.y, a.pos.z]
            for ch in st[0] for res in ch
            if res.name.strip().upper() in ALL_PROTEIN_RESNAMES
            for a in res if a.element.name not in ("H", "D")
        ]
        if all_p:
            origins.append(("Protein Core Centroid", np.mean(np.array(all_p), axis=0)))

    return origins[:3]


def _build_voronoi_graph(
    atom_pts: np.ndarray,
    atom_radii: np.ndarray,
    bottleneck_min: float = 0.9,
    probe_surface: float = 4.2,
) -> Dict[str, Any]:
    """
    Construct the Voronoi medial-axis void graph from 3D Delaunay triangulation.
    Each node is a Delaunay tetrahedron circumcenter with exact clearance radius `r_free`.
    """
    n_atoms = len(atom_pts)
    # Add 8 distant bounding box points so surface exterior tetrahedra exist cleanly
    pmin = np.min(atom_pts, axis=0)
    pmax = np.max(atom_pts, axis=0)
    pad = max(12.0, probe_surface * 3.0)
    bmin = pmin - pad
    bmax = pmax + pad

    corners = np.array([
        [bmin[0], bmin[1], bmin[2]], [bmax[0], bmin[1], bmin[2]],
        [bmin[0], bmax[1], bmin[2]], [bmax[0], bmax[1], bmin[2]],
        [bmin[0], bmin[1], bmax[2]], [bmax[0], bmin[1], bmax[2]],
        [bmin[0], bmax[1], bmax[2]], [bmax[0], bmax[1], bmax[2]],
    ], dtype=np.float64)

    all_pts = np.vstack([atom_pts, corners])
    tri = Delaunay(all_pts)
    simplices = tri.simplices
    neighbors = tri.neighbors

    centers, _ = _circumcenter_batch(all_pts, simplices)
    atom_tree = cKDTree(atom_pts)

    # Filter out NaN or wildly distant circumcenters (outside padded protein box)
    valid_mask = (
        ~np.isnan(centers[:, 0])
        & (centers[:, 0] >= pmin[0] - 8.0) & (centers[:, 0] <= pmax[0] + 8.0)
        & (centers[:, 1] >= pmin[1] - 8.0) & (centers[:, 1] <= pmax[1] + 8.0)
        & (centers[:, 2] >= pmin[2] - 8.0) & (centers[:, 2] <= pmax[2] + 8.0)
    )
    valid_indices = np.where(valid_mask)[0]
    if len(valid_indices) == 0:
        return {"nodes": [], "adj": {}, "exits": set(), "atom_tree": atom_tree}

    v_centers = centers[valid_indices]
    # Exact clearance to nearest protein atom surface (accounting for vdW radius)
    dists, nearest_idx = atom_tree.query(v_centers, k=1)
    r_free_arr = dists - atom_radii[nearest_idx]

    # Map original simplex index -> compact node index
    pass_mask = r_free_arr >= bottleneck_min
    kept_orig = valid_indices[pass_mask]
    kept_centers = v_centers[pass_mask]
    kept_rfree = r_free_arr[pass_mask]

    if len(kept_orig) == 0:
        return {"nodes": [], "adj": {}, "exits": set(), "atom_tree": atom_tree}

    orig_to_compact = {int(orig): i for i, orig in enumerate(kept_orig)}

    # Identify exterior / bulk solvent exit nodes:
    # A node is an exit if its clearance >= probe_surface OR any vertex of its simplex is a bounding-box corner
    # OR it sits on the outer solvent boundary (few protein atoms within 9 Å)
    exits = set()
    nbr_counts = [len(x) for x in atom_tree.query_ball_point(kept_centers, r=8.5)]
    for idx, orig_s in enumerate(kept_orig):
        simp = simplices[orig_s]
        has_corner = bool(np.any(simp >= n_atoms))
        if has_corner or kept_rfree[idx] >= probe_surface or nbr_counts[idx] < 10:
            exits.add(idx)

    # Build weighted adjacency list between adjacent Delaunay tetrahedra
    adj: Dict[int, List[Tuple[int, float, float]]] = {i: [] for i in range(len(kept_orig))}
    for i, orig_s in enumerate(kept_orig):
        c_i = kept_centers[i]
        r_i = kept_rfree[i]
        for nb_orig in neighbors[orig_s]:
            j = orig_to_compact.get(int(nb_orig))
            if j is None or j <= i:
                continue
            seg_len = float(np.linalg.norm(kept_centers[j] - c_i))
            if seg_len > 9.5:
                continue
            r_edge = min(r_i, kept_rfree[j])
            # MOLE cost function: rewards wider radius strongly and shorter length
            base_cost = seg_len / ((max(r_edge, 0.35) + 0.15) ** 2)
            adj[i].append((j, base_cost, seg_len))
            adj[j].append((i, base_cost, seg_len))

    return {
        "centers": kept_centers,
        "r_free": kept_rfree,
        "adj": adj,
        "exits": exits,
        "nbr_counts": nbr_counts,
        "atom_tree": atom_tree,
    }


def _dijkstra_to_exits(
    start_node: int,
    adj: Dict[int, List[Tuple[int, float, float]]],
    exits: set,
    node_penalties: Dict[int, float],
    blocked_exits: set,
) -> Optional[Tuple[List[int], float]]:
    """Find minimum-cost path from start_node to any allowed exterior exit node."""
    dist: Dict[int, float] = {start_node: 0.0}
    parent: Dict[int, int] = {}
    pq: List[Tuple[float, int]] = [(0.0, start_node)]

    while pq:
        d_u, u = heapq.heappop(pq)
        if d_u > dist.get(u, 1e18):
            continue
        if u != start_node and u in exits and u not in blocked_exits:
            # Reconstruct path
            path = [u]
            cur = u
            while cur in parent:
                cur = parent[cur]
                path.append(cur)
            path.reverse()
            return path, d_u

        for v, base_w, _ in adj.get(u, ()):
            # Do not traverse through other exterior nodes once outside
            if u in exits and u != start_node:
                continue
            pen = node_penalties.get(v, 1.0)
            nd = d_u + base_w * pen
            if nd < dist.get(v, 1e18):
                dist[v] = nd
                parent[v] = u
                heapq.heappush(pq, (nd, v))
    return None


def _smooth_and_sample_path(
    raw_pts: np.ndarray,
    atom_tree: cKDTree,
    atom_pts: np.ndarray,
    atom_radii: np.ndarray,
    atom_meta: List[Dict[str, Any]],
    step_A: float = 0.65,
) -> Dict[str, Any]:
    """
    Smooth the raw Voronoi centerline (moving-average preserving endpoints), resample
    at uniform arc-length intervals, recompute exact clearance radii & MOLE physicochemical profiles.
    """
    pts = raw_pts.copy()
    # Gentle Laplacian smoothing (2 passes, endpoints pinned)
    if len(pts) >= 4:
        for _ in range(2):
            smoothed = pts.copy()
            smoothed[1:-1] = 0.25 * pts[:-2] + 0.50 * pts[1:-1] + 0.25 * pts[2:]
            pts = smoothed

    # Arc-length resampling
    seg_lens = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg_lens)])
    total_len = float(cum[-1])
    if total_len < 0.5:
        return {}

    n_samples = max(6, min(90, int(round(total_len / step_A)) + 1))
    s_grid = np.linspace(0.0, total_len, n_samples)
    sx = np.interp(s_grid, cum, pts[:, 0])
    sy = np.interp(s_grid, cum, pts[:, 1])
    sz = np.interp(s_grid, cum, pts[:, 2])
    sampled = np.column_stack([sx, sy, sz])

    # Compute exact clearance radius & local lining residues at every sample point
    profile = []
    lining_residues: Dict[str, Dict[str, Any]] = {}
    bottleneck_radius = 999.0
    bottleneck_idx = 0
    cost_sum = 0.0

    for idx, p in enumerate(sampled):
        # Find nearby atoms within 6.5 Å to compute exact free radius & lining properties
        nbr_idxs = atom_tree.query_ball_point(p, r=7.5)
        if not nbr_idxs:
            d1, i1 = atom_tree.query(p, k=1)
            r_free = float(d1 - atom_radii[int(i1)])
            nbr_idxs = [int(i1)]
        else:
            dists = np.linalg.norm(atom_pts[nbr_idxs] - p, axis=1) - atom_radii[nbr_idxs]
            r_free = float(np.min(dists))

        r_clamped = max(0.35, min(6.5, r_free))
        if r_clamped < bottleneck_radius:
            bottleneck_radius = r_clamped
            bottleneck_idx = idx

        if idx > 0:
            ds = float(s_grid[idx] - s_grid[idx - 1])
            cost_sum += ds / ((r_clamped + 0.2) ** 2)

        # Lining atoms: within r_free + 2.8 Å (direct wall contact)
        wall_cut = r_clamped + 2.8
        local_res_seen = set()
        slice_res_labels: List[str] = []
        slice_hydro = []
        slice_polar = []
        slice_charge = 0
        slice_mut = []

        for ai in nbr_idxs:
            d_surf = float(np.linalg.norm(atom_pts[ai] - p) - atom_radii[ai])
            if d_surf > wall_cut:
                continue
            m = atom_meta[ai]
            rk = m["key"]
            std = m["std_resname"]
            if rk not in local_res_seen:
                local_res_seen.add(rk)
                slice_res_labels.append(f"{m['resname']} {rk}")
                if std in KYTE_DOOLITTLE:
                    slice_hydro.append(KYTE_DOOLITTLE[std])
                    slice_polar.append(POLARITY_SCALE.get(std, 0.0))
                    slice_charge += CHARGE_SCALE.get(m["resname"], CHARGE_SCALE.get(std, 0))
                    slice_mut.append(MUTABILITY_SCALE.get(std, 70.0))

            if rk not in lining_residues:
                lining_residues[rk] = {
                    "chain": m["chain"],
                    "resseq": m["resseq"],
                    "icode": m["icode"],
                    "resname": m["resname"],
                    "std_resname": std,
                    "key": rk,
                    "min_dist_A": round(d_surf, 2),
                    "start_dist": round(float(s_grid[idx]), 1),
                    "end_dist": round(float(s_grid[idx]), 1),
                    "backbone_contact": m["is_backbone"],
                    "hydropathy": KYTE_DOOLITTLE.get(std, 0.0),
                    "charge": CHARGE_SCALE.get(m["resname"], CHARGE_SCALE.get(std, 0)),
                    "is_bottleneck": False,
                }
            else:
                entry = lining_residues[rk]
                entry["end_dist"] = round(float(s_grid[idx]), 1)
                if d_surf < entry["min_dist_A"]:
                    entry["min_dist_A"] = round(d_surf, 2)
                if m["is_backbone"]:
                    entry["backbone_contact"] = True

        profile.append({
            "distance": round(float(s_grid[idx]), 2),
            "x": round(float(p[0]), 2),
            "y": round(float(p[1]), 2),
            "z": round(float(p[2]), 2),
            "radius": round(r_clamped, 2),
            "hydropathy": round(float(np.mean(slice_hydro)) if slice_hydro else 0.0, 2),
            "polarity": round(float(np.mean(slice_polar)) if slice_polar else 0.0, 2),
            "charge": int(slice_charge),
            "mutability": round(float(np.mean(slice_mut)) if slice_mut else 70.0, 1),
            "residues": sorted(slice_res_labels),
            "is_bottleneck": False,
        })

    # Mark exact bottleneck slice in profile
    if profile:
        profile[bottleneck_idx]["is_bottleneck"] = True

    # Tag bottleneck-forming residues (atoms closest to bottleneck sphere)
    bn_pt = sampled[bottleneck_idx]
    bn_res_keys = []
    for ai in atom_tree.query_ball_point(bn_pt, r=bottleneck_radius + 3.0):
        rk = atom_meta[ai]["key"]
        if rk in lining_residues and not lining_residues[rk]["is_bottleneck"]:
            lining_residues[rk]["is_bottleneck"] = True
            bn_res_keys.append(f"{lining_residues[rk]['resname']} {rk}")

    euclid = float(np.linalg.norm(sampled[-1] - sampled[0]))
    tortuosity = round(total_len / max(euclid, 0.5), 2)
    throughput = round(math.exp(-min(cost_sum * 0.12, 12.0)), 3)

    lining_list = sorted(lining_residues.values(), key=lambda x: (x["start_dist"], x["min_dist_A"]))
    avg_hydro = round(float(np.mean([p["hydropathy"] for p in profile])), 2)
    avg_polar = round(float(np.mean([p["polarity"] for p in profile])), 2)
    avg_mut = round(float(np.mean([p["mutability"] for p in profile])), 1)
    net_charge = sum(r["charge"] for r in lining_list)
    pos_count = sum(1 for r in lining_list if r["charge"] > 0)
    neg_count = sum(1 for r in lining_list if r["charge"] < 0)

    return {
        "length_A": round(total_len, 2),
        "bottleneck_radius_A": round(bottleneck_radius, 2),
        "bottleneck_distance_A": profile[bottleneck_idx]["distance"],
        "bottleneck_point": [profile[bottleneck_idx]["x"], profile[bottleneck_idx]["y"], profile[bottleneck_idx]["z"]],
        "bottleneck_residues": bn_res_keys[:8],
        "max_radius_A": round(max(p["radius"] for p in profile), 2),
        "tortuosity": tortuosity,
        "throughput": throughput,
        "properties": {
            "hydropathy": avg_hydro,
            "polarity": avg_polar,
            "mutability": avg_mut,
            "net_charge": net_charge,
            "positive_residues": pos_count,
            "negative_residues": neg_count,
            "character": (
                "Hydrophobic" if avg_hydro > 0.6 else ("Hydrophilic / Polar" if avg_hydro < -0.6 else "Amphipathic")
            ),
        },
        "profile": profile,
        "lining_residues": lining_list,
    }


def analyze_tunnels_and_pores(
    st: gemmi.Structure,
    mode: str = "both",
    bottleneck_radius: float = 0.9,
    probe_radius: float = 4.2,
    min_length: float = 5.0,
    max_tunnels: int = 6,
    start_mode: str = "auto",
    start_residues: Optional[List[str]] = None,
    start_point: Optional[List[float]] = None,
    include_hetatms: bool = False,
) -> Dict[str, Any]:
    """
    Run MOLE 2.5-style Voronoi tunnel & pore detection on a Gemmi Structure.
    """
    atom_pts, atom_radii, atom_meta = _extract_atoms(st, include_hetatms=include_hetatms)
    if len(atom_pts) < 20:
        return {
            "tunnels": [],
            "pores": [],
            "origins": [],
            "params": {
                "mode": mode,
                "bottleneck_radius": bottleneck_radius,
                "probe_radius": probe_radius,
                "min_length": min_length,
            },
            "summary": {"total_tunnels": 0, "total_pores": 0, "voronoi_nodes": 0},
        }

    graph = _build_voronoi_graph(
        atom_pts, atom_radii, bottleneck_min=bottleneck_radius, probe_surface=probe_radius
    )
    centers = graph.get("centers")
    if centers is None or len(centers) < 6:
        return {
            "tunnels": [],
            "pores": [],
            "origins": [],
            "params": {
                "mode": mode,
                "bottleneck_radius": bottleneck_radius,
                "probe_radius": probe_radius,
                "min_length": min_length,
            },
            "summary": {"total_tunnels": 0, "total_pores": 0, "voronoi_nodes": 0},
        }

    node_tree = cKDTree(centers)
    adj = graph["adj"]
    exits = graph["exits"]
    atom_tree = graph["atom_tree"]

    origins = _resolve_start_points(st, start_mode, start_residues, start_point)
    origins_meta = [{"label": lbl, "point": [round(float(x), 2) for x in pt]} for lbl, pt in origins]

    tunnels: List[Dict[str, Any]] = []
    pores: List[Dict[str, Any]] = []

    # Precompute which Voronoi nodes are reachable from bulk solvent exits (MOLE origin snapping)
    reachable_from_exit = set(exits)
    bfs_q = list(exits)
    qi = 0
    while qi < len(bfs_q):
        u = bfs_q[qi]
        qi += 1
        for v, _, _ in adj.get(u, ()):
            if v not in reachable_from_exit:
                reachable_from_exit.add(v)
                bfs_q.append(v)

    # 1. TUNNEL DISCOVERY (Origin -> Distinct Solvent Exits)
    if mode in ("tunnels", "both") and origins:
        seen_exit_pts: List[np.ndarray] = []
        per_origin_limit = max(2, math.ceil(max_tunnels / len(origins)) + 1)

        for o_label, o_pt in origins:
            # Snap starting origin to the nearest interior Voronoi vertex that connects to bulk solvent
            k_cand = min(80, len(centers))
            _, cand_ids = node_tree.query(o_pt, k=k_cand)
            if np.isscalar(cand_ids):
                cand_ids = [int(cand_ids)]
            start_node = next(
                (int(ci) for ci in cand_ids if int(ci) in reachable_from_exit and int(ci) not in exits),
                next((int(ci) for ci in cand_ids if int(ci) in reachable_from_exit), int(cand_ids[0])),
            )

            node_penalties: Dict[int, float] = {}
            blocked_exits: set = set()

            for _ in range(per_origin_limit * 2):
                if len(tunnels) >= max_tunnels:
                    break
                found = _dijkstra_to_exits(start_node, adj, exits, node_penalties, blocked_exits)
                if found is None:
                    break
                path_nodes, _ = found
                exit_node = path_nodes[-1]
                exit_pt = centers[exit_node]

                # Block a 5.5 Å neighborhood around this exit so next search finds a genuinely different opening
                for near_e in node_tree.query_ball_point(exit_pt, r=5.5):
                    if near_e in exits:
                        blocked_exits.add(int(near_e))

                # Penalize intermediate nodes to diversify tunnel routes
                for pn in path_nodes[1:-1]:
                    node_penalties[pn] = node_penalties.get(pn, 1.0) * 2.6

                # Skip duplicate exit openings too close to an already accepted tunnel
                if any(float(np.linalg.norm(exit_pt - ep)) < 5.0 for ep in seen_exit_pts):
                    continue

                raw_pts = centers[path_nodes]
                analyzed = _smooth_and_sample_path(raw_pts, atom_tree, atom_pts, atom_radii, atom_meta)
                if not analyzed or analyzed["length_A"] < min_length:
                    continue
                if analyzed["bottleneck_radius_A"] < bottleneck_radius * 0.85:
                    continue

                seen_exit_pts.append(exit_pt)
                t_rank = len(tunnels) + 1
                tunnels.append({
                    "id": f"T{t_rank}",
                    "kind": "tunnel",
                    "name": f"Tunnel {t_rank}",
                    "origin_label": o_label,
                    "color": TUNNEL_PALETTE[(t_rank - 1) % len(TUNNEL_PALETTE)],
                    **analyzed,
                })

        tunnels.sort(key=lambda t: (t["throughput"], t["bottleneck_radius_A"]), reverse=True)
        for i, t in enumerate(tunnels, start=1):
            t["id"] = f"T{i}"
            t["name"] = f"Tunnel {i}"
            t["color"] = TUNNEL_PALETTE[(i - 1) % len(TUNNEL_PALETTE)]

    # 2. PORE / CHANNEL DISCOVERY (Through-Protein Surface Exit A -> Interior -> Opposite Exit B)
    if mode in ("pores", "both"):
        # Pick deepest interior Voronoi nodes as candidate pore midpoints (like MOLE Pore mode!)
        interior_nodes = [
            i for i in range(len(centers))
            if i not in exits and graph["nbr_counts"][i] >= 16 and graph["r_free"][i] >= bottleneck_radius
        ]
        # Rank interior seeds by burial * clearance radius
        interior_nodes.sort(key=lambda i: graph["r_free"][i] * min(graph["nbr_counts"][i], 45), reverse=True)

        used_pore_pairs: List[Tuple[np.ndarray, np.ndarray]] = []
        for seed in interior_nodes[:25]:
            if len(pores) >= 3:
                break
            # Find first exit A
            res_a = _dijkstra_to_exits(seed, adj, exits, {}, set())
            if res_a is None:
                continue
            path_a, _ = res_a
            exit_a = centers[path_a[-1]]

            # Find second exit B on the OPPOSITE side of the protein (dot product < -0.25 & distance >= 12 Å)
            v_a = exit_a - centers[seed]
            va_norm = float(np.linalg.norm(v_a))
            if va_norm < 3.0:
                continue
            u_a = v_a / va_norm

            blocked_b = set()
            for e_idx in exits:
                v_e = centers[e_idx] - centers[seed]
                ve_n = float(np.linalg.norm(v_e))
                if ve_n < 3.0 or float(np.dot(v_e / ve_n, u_a)) > -0.20:
                    blocked_b.add(e_idx)
            pen_b = {n: 6.0 for n in path_a[1:]}
            res_b = _dijkstra_to_exits(seed, adj, exits, pen_b, blocked_b)
            if res_b is None:
                continue
            path_b, _ = res_b
            exit_b = centers[path_b[-1]]

            if float(np.linalg.norm(exit_a - exit_b)) < max(10.0, min_length * 1.6):
                continue

            # Deduplicate pores sharing nearly identical entrance/exit openings
            dup = False
            for ea, eb in used_pore_pairs:
                if (
                    (np.linalg.norm(exit_a - ea) < 7.0 and np.linalg.norm(exit_b - eb) < 7.0)
                    or (np.linalg.norm(exit_a - eb) < 7.0 and np.linalg.norm(exit_b - ea) < 7.0)
                ):
                    dup = True
                    break
            if dup:
                continue

            full_nodes = list(reversed(path_a)) + path_b[1:]
            raw_pts = centers[full_nodes]
            analyzed = _smooth_and_sample_path(raw_pts, atom_tree, atom_pts, atom_radii, atom_meta)
            if not analyzed or analyzed["length_A"] < max(10.0, min_length * 1.5):
                continue
            if analyzed["bottleneck_radius_A"] < bottleneck_radius * 0.85:
                continue

            used_pore_pairs.append((exit_a, exit_b))
            p_rank = len(pores) + 1
            pores.append({
                "id": f"P{p_rank}",
                "kind": "pore",
                "name": f"Pore {p_rank}",
                "origin_label": "Trans-Protein Channel",
                "color": "#AF52DE" if p_rank == 1 else ("#5AC8FA" if p_rank == 2 else "#FF2D55"),
                **analyzed,
            })

        pores.sort(key=lambda p: (p["bottleneck_radius_A"], p["length_A"]), reverse=True)
        for i, p in enumerate(pores, start=1):
            p["id"] = f"P{i}"
            p["name"] = f"Pore {i}"

    return {
        "tunnels": tunnels,
        "pores": pores,
        "origins": origins_meta,
        "params": {
            "mode": mode,
            "bottleneck_radius": bottleneck_radius,
            "probe_radius": probe_radius,
            "min_length": min_length,
            "max_tunnels": max_tunnels,
            "start_mode": start_mode,
            "include_hetatms": include_hetatms,
        },
        "summary": {
            "total_tunnels": len(tunnels),
            "total_pores": len(pores),
            "voronoi_nodes": len(centers),
            "widest_bottleneck_A": round(
                max([t["bottleneck_radius_A"] for t in tunnels + pores], default=0.0), 2
            ),
        },
    }


def _csv_escape(val: Any) -> str:
    s = str(val)
    if any(ch in s for ch in (",", '"', "\n")):
        return '"' + s.replace('"', '""') + '"'
    return s


def export_mole_csv(mole_data: Dict[str, Any]) -> str:
    """
    Export all detected tunnel & pore centerline profiles as CSV, including
    per-slice lining residues, bottleneck flags, bottleneck residues, and full channel residue lists.
    """
    header = [
        "id", "kind", "name", "origin", "length_A", "bottleneck_radius_A", "throughput",
        "distance_A", "x", "y", "z", "radius_A", "is_bottleneck",
        "hydropathy", "polarity", "charge", "mutability",
        "slice_residues", "bottleneck_residues", "all_lining_residues",
    ]
    lines = [",".join(header)]
    for item in mole_data.get("tunnels", []) + mole_data.get("pores", []):
        bn_res_str = "; ".join(item.get("bottleneck_residues", []))
        all_res_str = "; ".join(
            f"{r['resname']} {r['key']}" for r in item.get("lining_residues", [])
        )
        for pt in item.get("profile", []):
            slice_res_str = "; ".join(pt.get("residues", []))
            row = [
                item["id"],
                item["kind"],
                item["name"],
                item.get("origin_label", ""),
                item.get("length_A", 0.0),
                item.get("bottleneck_radius_A", 0.0),
                item.get("throughput", 0.0),
                pt["distance"],
                pt["x"],
                pt["y"],
                pt["z"],
                pt["radius"],
                "yes" if pt.get("is_bottleneck") else "no",
                pt["hydropathy"],
                pt["polarity"],
                pt["charge"],
                pt["mutability"],
                slice_res_str,
                bn_res_str,
                all_res_str,
            ]
            lines.append(",".join(_csv_escape(v) for v in row))
    return "\n".join(lines) + "\n"


def export_mole_residues_csv(mole_data: Dict[str, Any]) -> str:
    """Export one row per lining residue per tunnel/pore (ideal for mutation design & spreadsheets)."""
    header = [
        "channel_id", "kind", "channel_name", "origin", "bottleneck_radius_A",
        "chain", "resseq", "icode", "resname", "residue_key",
        "is_bottleneck", "start_dist_A", "end_dist_A", "min_dist_A",
        "contact_type", "hydropathy", "charge",
    ]
    lines = [",".join(header)]
    for item in mole_data.get("tunnels", []) + mole_data.get("pores", []):
        for lr in item.get("lining_residues", []):
            row = [
                item["id"],
                item["kind"],
                item["name"],
                item.get("origin_label", ""),
                item.get("bottleneck_radius_A", 0.0),
                lr["chain"],
                lr["resseq"],
                lr.get("icode", ""),
                lr["resname"],
                lr["key"],
                "yes" if lr.get("is_bottleneck") else "no",
                lr["start_dist"],
                lr["end_dist"],
                lr["min_dist_A"],
                "backbone+sidechain" if lr.get("backbone_contact") else "sidechain",
                lr["hydropathy"],
                lr["charge"],
            ]
            lines.append(",".join(_csv_escape(v) for v in row))
    return "\n".join(lines) + "\n"


def export_mole_pdb(mole_data: Dict[str, Any]) -> str:
    """
    Export MOLE tunnels & pores as PyMOL/ChimeraX-loadable pseudo-atom PDB
    (radius stored in B-factor column so `alter all, vdw=b; show spheres` renders exact MOLE geometry!).
    """
    lines = [
        "REMARK  Generated by PandaProtein — MOLE Tunnels & Pores Geometry",
        "REMARK  PyMOL tip: alter resn TUN+POR, vdw=b; show spheres, resn TUN+POR",
    ]
    serial = 1
    chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    all_items = mole_data.get("tunnels", []) + mole_data.get("pores", [])
    for idx, item in enumerate(all_items):
        ch = chain_letters[idx % len(chain_letters)]
        rname = "TUN" if item["kind"] == "tunnel" else "POR"
        for ri, pt in enumerate(item.get("profile", []), start=1):
            lines.append(
                f"HETATM{serial:5d}  C   {rname} {ch}{ri:4d}    "
                f"{pt['x']:8.3f}{pt['y']:8.3f}{pt['z']:8.3f}  1.00{pt['radius']:6.2f}           C"
            )
            serial += 1
        lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def scan_bottleneck_mutations(
    st: gemmi.Structure,
    bottleneck_residues: List[str],
    bottleneck_point: List[float],
    baseline_bottleneck_A: float,
    target_aa: str = "ALA",
    include_hetatms: bool = False,
) -> Dict[str, Any]:
    """
    In-silico MOLE Tunnel Engineering (Bottleneck Mutation Scanner):
    Tests mutating each bottleneck-forming residue on an isolated clone to `target_aa`
    (e.g. ALA or GLY) via COOT Mutate & Auto-Fit, recomputes exact local channel clearance
    at the bottleneck and along the local neighborhood, and ranks mutations by ΔBottleneck (Å)!
    """
    from .coot_engine import mutate_residue

    target = target_aa.strip().upper()
    if target not in ("ALA", "GLY", "VAL", "SER"):
        target = "ALA"

    bn_pt = np.array(bottleneck_point, dtype=np.float64)
    results: List[Dict[str, Any]] = []
    seen = set()

    for label in bottleneck_residues[:10]:
        # Parse formats like "LEU A:141" or "A:141"
        parts = label.strip().split()
        key = parts[-1] if parts else ""
        if ":" not in key or key in seen:
            continue
        seen.add(key)
        ch_id, num_str = key.split(":", 1)
        digits = "".join(c for c in num_str if c.isdigit() or c == "-")
        icode = "".join(c for c in num_str if c.isalpha())
        if not digits:
            continue
        resseq = int(digits)

        # Find current residue name
        wt_resname = parts[0] if len(parts) >= 2 else ""
        std_wt = canonical_resname(wt_resname)
        mut_target = "GLY" if (std_wt == target and target == "ALA") else target
        if std_wt == mut_target:
            continue

        try:
            clone = st.clone()
            mutate_residue(clone, ch_id, resseq, mut_target, icode)
            # Recompute exact clearance radius at the bottleneck sphere (and nearby 1.2 Å probe grid)
            pts, rads, _ = _extract_atoms(clone, include_hetatms=include_hetatms)
            if len(pts) < 10:
                continue
            tree = cKDTree(pts)
            # Evaluate max achievable clearance within 1.5 Å of the original bottleneck center
            offsets = np.array([
                [0, 0, 0], [0.8, 0, 0], [-0.8, 0, 0], [0, 0.8, 0], [0, -0.8, 0], [0, 0, 0.8], [0, 0, -0.8],
                [0.6, 0.6, 0], [-0.6, 0.6, 0], [0.6, -0.6, 0], [-0.6, -0.6, 0],
            ], dtype=np.float64)
            best_r = 0.0
            for off in offsets:
                q_pt = bn_pt + off
                nbrs = tree.query_ball_point(q_pt, r=6.5)
                if nbrs:
                    r_here = float(np.min(np.linalg.norm(pts[nbrs] - q_pt, axis=1) - rads[nbrs]))
                    if r_here > best_r:
                        best_r = r_here

            new_r = round(max(0.35, best_r), 2)
            delta_r = round(new_r - baseline_bottleneck_A, 2)
            results.append({
                "chain": ch_id,
                "resseq": resseq,
                "icode": icode,
                "key": key,
                "wt_resname": wt_resname or std_wt,
                "mut_resname": mut_target,
                "mutation": f"{wt_resname or std_wt}{resseq}{mut_target} ({ch_id})",
                "baseline_bottleneck_A": round(baseline_bottleneck_A, 2),
                "new_bottleneck_A": new_r,
                "delta_bottleneck_A": delta_r,
                "effect": "Widens Tunnel" if delta_r >= 0.15 else ("Neutral" if delta_r >= -0.05 else "Narrows"),
            })
        except Exception:
            continue

    results.sort(key=lambda x: x["delta_bottleneck_A"], reverse=True)
    return {
        "target_aa": target,
        "baseline_bottleneck_A": round(baseline_bottleneck_A, 2),
        "candidates": results,
    }
