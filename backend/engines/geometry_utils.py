"""
High-precision 3D molecular geometry utilities:
- Bond lengths, angles, dihedrals, signed chiral volumes
- NeRF (Natural Extension Reference Frame) internal-coordinate placement
- Rodrigues rotation for rotamers, pepflips, cis-trans flips, and NQH flips
- Ramachandran region classifier (General, Glycine, Proline, Pre-Proline)
- Ideal side-chain builder for missing-atom reconstruction and point mutations
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from .chemistry_data import CHI_ATOMS, PROTONATION_SYNONYMS, SIDECHAIN_IC


def canonical_resname(resname: str) -> str:
    r = resname.strip().upper()
    return PROTONATION_SYNONYMS.get(r, r)


def distance(p0: np.ndarray, p1: np.ndarray) -> float:
    return float(np.linalg.norm(p1 - p0))


def angle_deg(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> float:
    v1 = p0 - p1
    v2 = p2 - p1
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < 1e-8 or n2 < 1e-8:
        return 0.0
    cos_a = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
    return float(math.degrees(math.acos(cos_a)))


def dihedral_deg(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1_norm = np.linalg.norm(b1)
    if b1_norm < 1e-8:
        return 0.0
    b1 = b1 / b1_norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(math.degrees(math.atan2(y, x)))


def wrap_angle(deg: float) -> float:
    return ((deg + 180.0) % 360.0) - 180.0


def angular_diff(a: float, b: float) -> float:
    return abs(wrap_angle(a - b))


def chiral_volume(center: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Signed triple product (a - center) . ((b - center) x (c - center))."""
    return float(np.dot(a - center, np.cross(b - center, c - center)))


def nerf_place(
    p_a: np.ndarray,
    p_b: np.ndarray,
    p_c: np.ndarray,
    bond_len: float,
    bond_angle_deg: float,
    torsion_deg: float,
) -> np.ndarray:
    """
    Place atom D given A-B-C coordinates, bond length |C-D|, angle (B-C-D),
    and dihedral (A-B-C-D) using the Natural Extension Reference Frame (NeRF).
    """
    theta = math.radians(180.0 - bond_angle_deg)
    phi = math.radians(torsion_deg)

    bc = p_c - p_b
    bc_norm = np.linalg.norm(bc)
    if bc_norm < 1e-8:
        bc_hat = np.array([1.0, 0.0, 0.0])
    else:
        bc_hat = bc / bc_norm

    ab = p_b - p_a
    n = np.cross(ab, bc_hat)
    n_norm = np.linalg.norm(n)
    if n_norm < 1e-8:
        # Degenerate collinear fallback
        ref = np.array([0.0, 1.0, 0.0]) if abs(bc_hat[0]) > 0.9 else np.array([1.0, 0.0, 0.0])
        n = np.cross(ref, bc_hat)
        n_norm = np.linalg.norm(n)
    n_hat = n / n_norm
    m_hat = np.cross(n_hat, bc_hat)

    d_local = np.array([
        bond_len * math.cos(theta),
        bond_len * math.sin(theta) * math.cos(phi),
        bond_len * math.sin(theta) * math.sin(phi),
    ])
    return p_c + d_local[0] * bc_hat + d_local[1] * m_hat + d_local[2] * n_hat


def rotate_point_around_axis(
    point: np.ndarray, origin: np.ndarray, axis: np.ndarray, angle_deg_val: float
) -> np.ndarray:
    """Rotate a 3D point around an axis passing through origin by angle_deg_val (Rodrigues)."""
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-8:
        return point.copy()
    k = axis / axis_norm
    theta = math.radians(angle_deg_val)
    v = point - origin
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)
    rotated = v * cos_t + np.cross(k, v) * sin_t + k * np.dot(k, v) * (1.0 - cos_t)
    return origin + rotated


def ideal_cbeta_position(n_pos: np.ndarray, ca_pos: np.ndarray, c_pos: np.ndarray) -> np.ndarray:
    """Ideal L-amino acid C-beta position from backbone (N, CA, C) for C-beta deviation check."""
    return nerf_place(c_pos, n_pos, ca_pos, 1.530, 110.5, -122.6)


def build_sidechain_coords(
    backbone: Dict[str, np.ndarray],
    resname: str,
    chi_targets: Optional[List[float]] = None,
) -> Dict[str, np.ndarray]:
    """
    Build all ideal side-chain heavy atoms from backbone N, CA, C coordinates.
    Optionally applies target chi angles [chi1, chi2, chi3, chi4].
    """
    std = canonical_resname(resname)
    coords: Dict[str, np.ndarray] = {k: v.copy() for k, v in backbone.items()}
    if not all(k in coords for k in ("N", "CA", "C")):
        return coords

    specs = SIDECHAIN_IC.get(std, [])
    for atom_name, p1, p2, p3, bl, ba, dih in specs:
        if atom_name in coords:
            continue
        if p1 in coords and p2 in coords and p3 in coords:
            coords[atom_name] = nerf_place(coords[p3], coords[p2], coords[p1], bl, ba, dih)

    if chi_targets and std in CHI_ATOMS:
        apply_chi_angles_to_dict(coords, std, chi_targets)

    return coords


def apply_chi_angles_to_dict(
    coords: Dict[str, np.ndarray], std_resname: str, chi_targets: List[float]
) -> None:
    """Rotate side-chain atoms in-place to match desired chi angles."""
    chi_defs = CHI_ATOMS.get(std_resname, [])
    if not chi_defs:
        return

    # Downstream dependency order based on SIDECHAIN_IC order
    specs = SIDECHAIN_IC.get(std_resname, [])
    atom_order = [s[0] for s in specs]

    for idx, target_chi in enumerate(chi_targets):
        if idx >= len(chi_defs):
            break
        a0, a1, a2, a3 = chi_defs[idx]
        if not all(a in coords for a in (a0, a1, a2, a3)):
            continue
        current_chi = dihedral_deg(coords[a0], coords[a1], coords[a2], coords[a3])
        delta = wrap_angle(target_chi - current_chi)
        if abs(delta) < 1e-3:
            continue

        origin = coords[a1]
        axis = coords[a2] - coords[a1]
        # Atoms downstream of bond a1-a2: a3 and any sidechain atom after a2 in build order
        if a2 in atom_order:
            cut = atom_order.index(a2) + 1
            downstream = atom_order[cut:]
        else:
            downstream = [a for a in atom_order if a not in ("N", "CA", "C", "O", "CB")]
        if a3 not in downstream:
            downstream = [a3] + downstream

        for aname in downstream:
            if aname in coords and aname not in (a0, a1, a2):
                coords[aname] = rotate_point_around_axis(coords[aname], origin, axis, delta)


def rama_category(resname: str, next_resname: Optional[str]) -> str:
    std = canonical_resname(resname)
    nxt = canonical_resname(next_resname) if next_resname else ""
    if std == "GLY":
        return "glycine"
    if std == "PRO":
        return "proline"
    if nxt == "PRO":
        return "pre_proline"
    return "general"


def _in_ellipse(phi: float, psi: float, cx: float, cy: float, rx: float, ry: float, rot_deg: float = 0.0) -> bool:
    dx = wrap_angle(phi - cx)
    dy = wrap_angle(psi - cy)
    if rot_deg != 0.0:
        rad = math.radians(rot_deg)
        c, s = math.cos(rad), math.sin(rad)
        rx_p = dx * c + dy * s
        ry_p = -dx * s + dy * c
        dx, dy = rx_p, ry_p
    return (dx / rx) ** 2 + (dy / ry) ** 2 <= 1.0


def evaluate_ramachandran(phi: float, psi: float, category: str) -> str:
    """
    Classify (phi, psi) into 'favored', 'allowed', or 'outlier' using
    Top8000/MolProbity-calibrated analytical contours for General, Glycine,
    Proline, and Pre-Proline residues.
    """
    if category == "glycine":
        # Symmetric alpha-R, alpha-L, beta, and extended polyproline regions
        favored_regions = [
            (-75, -35, 42, 42, 0),
            (75, 35, 42, 42, 0),
            (-115, 145, 55, 48, 0),
            (115, -145, 55, 48, 0),
            (85, -170, 50, 35, 0),
            (-85, 170, 50, 35, 0),
            (175, 175, 45, 45, 0),
            (-175, -175, 45, 45, 0),
        ]
        allowed_regions = [
            (-75, -30, 68, 70, 0),
            (75, 30, 68, 70, 0),
            (-120, 140, 80, 75, 0),
            (120, -140, 80, 75, 0),
            (180, 180, 75, 75, 0),
            (180, 0, 55, 55, 0),
            (0, 180, 55, 55, 0),
        ]
    elif category == "proline":
        favored_regions = [
            (-63, -38, 24, 30, 0),
            (-65, 145, 26, 36, 0),
            (-70, -10, 22, 26, 0),
        ]
        allowed_regions = [
            (-65, -35, 42, 54, 0),
            (-65, 142, 42, 56, 0),
            (-72, 20, 38, 46, 0),
            (-65, 180, 38, 34, 0),
        ]
    elif category == "pre_proline":
        favored_regions = [
            (-115, 138, 54, 42, 0),
            (-72, -32, 32, 30, 0),
            (-72, 145, 36, 36, 0),
            (55, 42, 24, 26, 0),
        ]
        allowed_regions = [
            (-112, 135, 76, 62, 0),
            (-75, -30, 56, 54, 0),
            (-115, 22, 60, 54, 0),
            (58, 42, 42, 46, 0),
        ]
    else:  # general (Richardson Top8000 98% favored / 99.8% allowed contours)
        favored_regions = [
            (-64, -41, 34, 28, -15),   # Right-handed alpha helix
            (-92, -8, 38, 28, -10),    # Bridge / 3_10 helix
            (-118, 132, 56, 44, -15),  # Beta sheet
            (-70, 142, 32, 34, 0),     # PPII
            (-140, 165, 44, 30, 0),    # Extended beta
            (58, 38, 24, 28, 0),       # Left-handed alpha
        ]
        allowed_regions = [
            (-75, -32, 62, 58, -12),
            (-115, 132, 78, 68, -12),
            (-112, 18, 64, 56, 0),     # Top8000 alpha-beta bridge region
            (-135, -55, 42, 44, 0),    # Extended lower-left helical fringe
            (-120, -168, 64, 36, 0),
            (60, 35, 44, 52, 0),
            (65, -118, 30, 34, 0),     # Classic gamma-turn / type-II' turn allowed basin
            (72, -45, 30, 34, 0),
        ]

    for cx, cy, rx, ry, rot in favored_regions:
        if _in_ellipse(phi, psi, cx, cy, rx, ry, rot):
            return "favored"
    for cx, cy, rx, ry, rot in allowed_regions:
        if _in_ellipse(phi, psi, cx, cy, rx, ry, rot):
            return "allowed"
    return "outlier"
