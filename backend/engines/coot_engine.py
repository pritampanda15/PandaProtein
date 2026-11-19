"""
COOT-Inspired Crystallographic Validation & Interactive Model-Building Engine:
- Ramachandran (phi/psi) validation across General, Glycine, Proline, Pre-Proline
- Dunbrack/Richardson Rotamer analysis, ranked rotamer library explorer, Auto-Fit Rotamer
- Engh & Huber bond-length / bond-angle RMSZ & C-beta deviation validation
- Cis/Twisted peptide detection + COOT 180° Pepflip & Cis↔Trans conversion
- Chirality validation (C-alpha L/D & Ile/Thr C-beta) + inversion repair
- Fast KDTree MolProbity-style all-atom steric clashscore & bump list
- Asn / Gln / His 180° side-chain flip scoring & execution
- COOT Mutate & Auto-Fit Residue
- COOT Zone Regularization (local harmonic + Ramachandran + steric minimization)
- Automated COOT Polish pipeline
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .chemistry_data import (
    ALL_PROTEIN_RESNAMES,
    CHI_ATOMS,
    EXPECTED_HEAVY_ATOMS,
    IDEAL_BACKBONE_ANGLES,
    IDEAL_BACKBONE_BONDS,
    NONSTANDARD_TO_STANDARD,
    ROTAMER_LIBRARY,
    STANDARD_AA,
    VDW_RADII,
)
from .geometry_utils import (
    angle_deg,
    angular_diff,
    apply_chi_angles_to_dict,
    build_sidechain_coords,
    canonical_resname,
    chiral_volume,
    dihedral_deg,
    distance,
    evaluate_ramachandran,
    ideal_cbeta_position,
    nerf_place,
    rama_category,
    rotate_point_around_axis,
    wrap_angle,
)


def _atom_pos(atom: gemmi.Atom) -> np.ndarray:
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=np.float64)


def _get_atom(res: gemmi.Residue, name: str) -> Optional[gemmi.Atom]:
    """Return primary (highest occupancy / first) atom matching name."""
    best = None
    for a in res:
        if a.name.strip() == name:
            if best is None or a.occ > best.occ:
                best = a
    return best


def _get_pos(res: gemmi.Residue, name: str) -> Optional[np.ndarray]:
    a = _get_atom(res, name)
    return _atom_pos(a) if a is not None else None


def _is_protein_res(res: gemmi.Residue) -> bool:
    r = res.name.strip().upper()
    return r in ALL_PROTEIN_RESNAMES or r in NONSTANDARD_TO_STANDARD


def _find_residue(st: gemmi.Structure, chain_id: str, resseq: int, icode: str = "") -> Tuple[Optional[gemmi.Chain], Optional[gemmi.Residue], int]:
    if len(st) == 0:
        return None, None, -1
    model = st[0]
    for chain in model:
        if chain.name == chain_id:
            for idx, res in enumerate(chain):
                if res.seqid.num == resseq and res.seqid.icode.strip() == icode.strip():
                    return chain, res, idx
    return None, None, -1


def _res_key(chain_id: str, res: gemmi.Residue) -> str:
    ic = res.seqid.icode.strip()
    return f"{chain_id}:{res.seqid.num}{ic}"


def _extract_protein_chains(st: gemmi.Structure) -> List[Tuple[str, List[gemmi.Residue]]]:
    if len(st) == 0:
        return []
    result = []
    for chain in st[0]:
        pres = [r for r in chain if _is_protein_res(r) and r.name.strip().upper() not in ("ACE", "NME", "NMA", "NH2")]
        if pres:
            result.append((chain.name, pres))
    return result


def _compute_clashes(st: gemmi.Structure, overlap_cutoff: float = 0.40) -> Dict[str, Any]:
    """Fast KDTree MolProbity-like all-atom steric clash calculation."""
    if len(st) == 0:
        return {"clashscore": 0.0, "total_clashes": 0, "total_atoms": 0, "clashes": [], "residue_clash_counts": {}}

    coords = []
    meta = []  # (chain_id, resseq, icode, resname, atom_name, elem, res_idx_global, altloc, is_polar)
    res_global_idx = -1
    prev_res_uid = None

    for chain in st[0]:
        for res in chain:
            res_uid = (chain.name, res.seqid.num, res.seqid.icode.strip())
            if res_uid != prev_res_uid:
                res_global_idx += 1
                prev_res_uid = res_uid
            rname = res.name.strip().upper()
            for atom in res:
                el = atom.element.name.upper()
                if el not in VDW_RADII:
                    el = atom.name.strip()[0].upper() if atom.name.strip() else "C"
                coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                meta.append((
                    chain.name,
                    res.seqid.num,
                    res.seqid.icode.strip(),
                    rname,
                    atom.name.strip(),
                    el,
                    res_global_idx,
                    atom.altloc if atom.altloc != "\x00" else "",
                    el in ("N", "O", "F"),
                ))

    n_atoms = len(coords)
    if n_atoms < 2:
        return {"clashscore": 0.0, "total_clashes": 0, "total_atoms": n_atoms, "clashes": [], "residue_clash_counts": {}}

    pts = np.array(coords, dtype=np.float64)
    tree = cKDTree(pts)
    pairs = tree.query_pairs(r=3.3)

    clashes = []
    res_clash_counts: Dict[str, int] = {}

    for i, j in pairs:
        m1 = meta[i]
        m2 = meta[j]

        # Skip same residue (intra-residue bonded/1-3/1-4)
        if m1[6] == m2[6]:
            continue
        # Skip incompatible altlocs (e.g. A vs B)
        if m1[7] and m2[7] and m1[7] != m2[7]:
            continue

        # Skip sequential backbone peptide-bond & adjacent Cbeta/H pairs (i and i+1)
        if m1[0] == m2[0] and abs(m1[1] - m2[1]) == 1:
            bb = {"N", "CA", "C", "O", "H", "HN", "HA", "HA2", "HA3", "CB", "CD"}
            if m1[4] in bb and m2[4] in bb:
                continue
            # Also skip adjacent-residue H-to-backbone/sidechain 1-4 contacts
            if m1[5] == "H" or m2[5] == "H":
                continue

        # Skip disulfide SG-SG bonds
        if m1[4] == "SG" and m2[4] == "SG":
            continue

        # Skip metal coordination bonds
        if m1[5] in ("ZN", "MG", "CA", "FE", "MN", "CU", "CO", "NI") or m2[5] in ("ZN", "MG", "CA", "FE", "MN", "CU", "CO", "NI"):
            continue

        d = float(np.linalg.norm(pts[i] - pts[j]))

        # Skip genuine inter-residue covalent bonds (glycans, covalent ligands, terminal caps)
        if d < 1.80 and m1[5] != "H" and m2[5] != "H":
            continue

        # Skip hydrogen-to-polar H-bonds (H...O or H...N down to 1.55 Å are H-bonds, not clashes)
        if (m1[5] == "H" and m2[5] in ("N", "O", "S")) or (m2[5] == "H" and m1[5] in ("N", "O", "S")):
            if d >= 1.50:
                continue

        # MolProbity calibrated radii (polar H = 1.00 Å, aliphatic H = 1.10 Å)
        r1 = 1.05 if m1[5] == "H" else VDW_RADII.get(m1[5], 1.70)
        r2 = 1.05 if m2[5] == "H" else VDW_RADII.get(m2[5], 1.70)

        # H-bond forgiveness for polar heavy atoms (N/O - N/O distance down to ~2.45 Å is H-bond, not clash)
        hb_allowance = 0.65 if (m1[8] and m2[8]) else 0.0
        overlap = (r1 + r2) - d - hb_allowance

        if overlap > overlap_cutoff:
            key1 = f"{m1[0]}:{m1[1]}{m1[2]}"
            key2 = f"{m2[0]}:{m2[1]}{m2[2]}"
            res_clash_counts[key1] = res_clash_counts.get(key1, 0) + 1
            res_clash_counts[key2] = res_clash_counts.get(key2, 0) + 1
            clashes.append({
                "atom1": {"chain": m1[0], "resseq": m1[1], "icode": m1[2], "resname": m1[3], "atom": m1[4], "key": key1},
                "atom2": {"chain": m2[0], "resseq": m2[1], "icode": m2[2], "resname": m2[3], "atom": m2[4], "key": key2},
                "distance": round(d, 3),
                "overlap": round(overlap, 3),
                "severity": "severe" if overlap >= 0.8 else "moderate",
                "midpoint": [round(float((pts[i][k] + pts[j][k]) / 2.0), 3) for k in range(3)],
            })

    clashes.sort(key=lambda x: x["overlap"], reverse=True)
    clashscore = round((len(clashes) / max(n_atoms, 1)) * 1000.0, 2)
    return {
        "clashscore": clashscore,
        "total_clashes": len(clashes),
        "total_atoms": n_atoms,
        "clashes": clashes[:120],
        "residue_clash_counts": res_clash_counts,
    }


def _evaluate_rotamer(res: gemmi.Residue) -> Optional[Dict[str, Any]]:
    std = canonical_resname(res.name)
    if std not in CHI_ATOMS or std not in ROTAMER_LIBRARY:
        return None
    chi_defs = CHI_ATOMS[std]
    chis: List[float] = []
    for a0, a1, a2, a3 in chi_defs:
        p0, p1, p2, p3 = _get_pos(res, a0), _get_pos(res, a1), _get_pos(res, a2), _get_pos(res, a3)
        if any(p is None for p in (p0, p1, p2, p3)):
            return None
        chis.append(round(dihedral_deg(p0, p1, p2, p3), 1))

    lib = ROTAMER_LIBRARY[std]
    best_dist = 999.0
    best_name = ""
    best_prob = 0.0
    for rot_name, prob, ref_chis in lib:
        n_cmp = min(len(chis), len(ref_chis))
        diffs = [angular_diff(chis[k], ref_chis[k]) for k in range(n_cmp)]
        # Symmetric aromatic / carboxylate chi2/chi3 180° periodicity
        if std in ("PHE", "TYR", "ASP") and n_cmp >= 2:
            diffs[1] = min(diffs[1], abs(180.0 - diffs[1]))
        if std == "GLU" and n_cmp >= 3:
            diffs[2] = min(diffs[2], abs(180.0 - diffs[2]))
        rmsd = math.sqrt(sum(d * d for d in diffs) / n_cmp)
        if rmsd < best_dist:
            best_dist = rmsd
            best_name = rot_name
            best_prob = prob

    # For flexible sp3 side chains (MET, LYS, ARG, GLU, GLN, ILE, LEU), all staggered wells
    # (-65, +62, 180) are physically populated in Richardson Top8000; only eclipsed angles are outliers.
    sp3_chi_counts = {"MET": 3, "LYS": 4, "ARG": 3, "GLU": 2, "GLN": 2, "ILE": 2, "LEU": 2}
    staggered_wells = (-65.0, 62.0, 180.0)
    staggered_rmsd = best_dist
    if std in sp3_chi_counts:
        n_sp3 = min(len(chis), sp3_chi_counts[std])
        stag_diffs = [min(angular_diff(chis[k], w) for w in staggered_wells) for k in range(n_sp3)]
        if stag_diffs:
            staggered_rmsd = math.sqrt(sum(d * d for d in stag_diffs) / len(stag_diffs))

    eff_dist = min(best_dist, staggered_rmsd)

    if eff_dist <= 30.0:
        status = "favored"
    elif eff_dist <= 45.0:
        status = "allowed"
    else:
        status = "outlier"

    return {
        "chis": chis,
        "nearest_rotamer": best_name,
        "library_prob": best_prob,
        "rmsd_deg": round(eff_dist, 1),
        "status": status,
    }


def _detect_nqh_flips(st: gemmi.Structure) -> List[Dict[str, Any]]:
    """
    Evaluate ASN, GLN, HIS residues for 180° COOT/Reduce side-chain flip improvement
    by comparing H-bond donor/acceptor geometry against nearby polar atoms.
    """
    if len(st) == 0:
        return []
    polar_pts = []
    polar_meta = []  # (res_uid, role: 'donor'|'acceptor'|'both', atom_name)
    for chain in st[0]:
        for res in chain:
            uid = (chain.name, res.seqid.num, res.seqid.icode.strip())
            rname = canonical_resname(res.name)
            for atom in res:
                an = atom.name.strip()
                el = atom.element.name.upper()
                if el not in ("N", "O"):
                    continue
                if an == "O" or an in ("OD1", "OD2", "OE1", "OE2"):
                    role = "acceptor"
                elif an == "N" or an in ("ND2", "NE2", "NE", "NH1", "NH2", "NZ", "NE1"):
                    role = "donor"
                else:
                    role = "both"
                if rname == "HOH":
                    role = "both"
                polar_pts.append([atom.pos.x, atom.pos.y, atom.pos.z])
                polar_meta.append((uid, role, an))

    if len(polar_pts) < 4:
        return []
    pts = np.array(polar_pts, dtype=np.float64)
    tree = cKDTree(pts)

    candidates = []
    for chain in st[0]:
        for res in chain:
            std = canonical_resname(res.name)
            if std not in ("ASN", "GLN", "HIS"):
                continue
            uid = (chain.name, res.seqid.num, res.seqid.icode.strip())
            if std == "ASN":
                a_acc, a_don = _get_pos(res, "OD1"), _get_pos(res, "ND2")
            elif std == "GLN":
                a_acc, a_don = _get_pos(res, "OE1"), _get_pos(res, "NE2")
            else:
                a_acc, a_don = _get_pos(res, "ND1"), _get_pos(res, "NE2")

            if a_acc is None or a_don is None:
                continue

            def _score_pair(p_acc: np.ndarray, p_don: np.ndarray) -> float:
                score = 0.0
                for idx in tree.query_ball_point(p_acc, 3.5):
                    m_uid, m_role, _ = polar_meta[idx]
                    if m_uid == uid:
                        continue
                    d = float(np.linalg.norm(pts[idx] - p_acc))
                    if d < 2.2:
                        score -= 3.0
                    elif m_role in ("donor", "both"):
                        score += 2.0
                    else:
                        score -= 1.2  # acceptor-acceptor repulsion
                for idx in tree.query_ball_point(p_don, 3.5):
                    m_uid, m_role, _ = polar_meta[idx]
                    if m_uid == uid:
                        continue
                    d = float(np.linalg.norm(pts[idx] - p_don))
                    if d < 2.2:
                        score -= 3.0
                    elif m_role in ("acceptor", "both"):
                        score += 2.0
                    else:
                        score -= 1.2  # donor-donor repulsion
                return score

            curr_score = _score_pair(a_acc, a_don)
            flipped_score = _score_pair(a_don, a_acc)
            gain = round(flipped_score - curr_score, 2)
            if gain >= 1.5:
                candidates.append({
                    "chain": chain.name,
                    "resseq": res.seqid.num,
                    "icode": res.seqid.icode.strip(),
                    "resname": res.name.strip(),
                    "key": _res_key(chain.name, res),
                    "current_score": round(curr_score, 2),
                    "flipped_score": round(flipped_score, 2),
                    "hbond_gain": gain,
                    "reason": "180° amide/imidazole flip improves H-bond donor/acceptor complementarity",
                })

    candidates.sort(key=lambda x: x["hbond_gain"], reverse=True)
    return candidates


def analyze_structure_quality(st: gemmi.Structure) -> Dict[str, Any]:
    """Run complete COOT + MolProbity crystallographic validation on a Gemmi Structure."""
    chains = _extract_protein_chains(st)
    clash_data = _compute_clashes(st)
    res_clashes = clash_data["residue_clash_counts"]

    rama_points: List[Dict[str, Any]] = []
    rama_counts = {"favored": 0, "allowed": 0, "outlier": 0}

    rotamer_points: List[Dict[str, Any]] = []
    rotamer_counts = {"favored": 0, "allowed": 0, "outlier": 0}

    bond_outliers: List[Dict[str, Any]] = []
    angle_outliers: List[Dict[str, Any]] = []
    bond_z_sq: List[float] = []
    angle_z_sq: List[float] = []

    cis_peptides: List[Dict[str, Any]] = []
    chiral_outliers: List[Dict[str, Any]] = []
    cbeta_outliers: List[Dict[str, Any]] = []
    bfactor_profile: List[Dict[str, Any]] = []
    residues_index: List[Dict[str, Any]] = []

    for chain_id, pres in chains:
        for i, res in enumerate(pres):
            rname = res.name.strip().upper()
            std = canonical_resname(rname)
            key = _res_key(chain_id, res)
            issues: List[str] = []

            n_p = _get_pos(res, "N")
            ca_p = _get_pos(res, "CA")
            c_p = _get_pos(res, "C")
            o_p = _get_pos(res, "O")
            cb_p = _get_pos(res, "CB")

            # Missing heavy atoms check
            present_heavy = {a.name.strip() for a in res if a.element.name not in ("H", "D")}
            expected = set(EXPECTED_HEAVY_ATOMS.get(std, []))
            missing_atoms = sorted(expected - present_heavy) if expected else []
            if missing_atoms:
                issues.append("missing_atoms")

            # Altloc check
            altlocs = sorted({a.altloc for a in res if a.altloc not in ("\x00", " ")})
            if len(altlocs) > 1:
                issues.append("altloc")

            # B-factor & occupancy
            bvals = [a.b_iso for a in res if a.element.name not in ("H", "D")]
            occs = [a.occ for a in res if a.element.name not in ("H", "D")]
            mean_b = round(float(np.mean(bvals)), 2) if bvals else 0.0
            min_occ = round(float(min(occs)), 2) if occs else 1.0
            bfactor_profile.append({
                "chain": chain_id,
                "resseq": res.seqid.num,
                "icode": res.seqid.icode.strip(),
                "resname": rname,
                "key": key,
                "mean_b": mean_b,
                "min_occ": min_occ,
            })

            # Intra-residue backbone bond lengths
            for (a1, a2), (ideal, sigma) in IDEAL_BACKBONE_BONDS.items():
                if a1 == "C_prev":
                    continue
                p1, p2 = _get_pos(res, a1), _get_pos(res, a2)
                if p1 is not None and p2 is not None:
                    d = distance(p1, p2)
                    z = (d - ideal) / sigma
                    bond_z_sq.append(z * z)
                    if abs(z) > 4.0:
                        issues.append("bond_outlier")
                        bond_outliers.append({
                            "chain": chain_id,
                            "resseq": res.seqid.num,
                            "icode": res.seqid.icode.strip(),
                            "resname": rname,
                            "key": key,
                            "atoms": f"{a1}-{a2}",
                            "observed": round(d, 3),
                            "ideal": ideal,
                            "z_score": round(z, 2),
                        })

            # Intra-residue backbone angles (Engh & Huber Proline N-CA-CB ideal is 103.0° due to pyrrolidine ring!)
            for (a1, a2, a3), (ideal, sigma) in IDEAL_BACKBONE_ANGLES.items():
                if "_prev" in a1 or "_prev" in a2:
                    continue
                if std == "PRO" and (a1, a2, a3) == ("N", "CA", "CB"):
                    ideal, sigma = 103.0, 2.2
                elif std == "PRO" and (a1, a2, a3) == ("C", "CA", "CB"):
                    ideal, sigma = 111.5, 2.4
                p1, p2, p3 = _get_pos(res, a1), _get_pos(res, a2), _get_pos(res, a3)
                if p1 is not None and p2 is not None and p3 is not None:
                    ang = angle_deg(p1, p2, p3)
                    z = (ang - ideal) / sigma
                    angle_z_sq.append(z * z)
                    if abs(z) > 4.0:
                        issues.append("angle_outlier")
                        angle_outliers.append({
                            "chain": chain_id,
                            "resseq": res.seqid.num,
                            "icode": res.seqid.icode.strip(),
                            "resname": rname,
                            "key": key,
                            "atoms": f"{a1}-{a2}-{a3}",
                            "observed": round(ang, 2),
                            "ideal": ideal,
                            "z_score": round(z, 2),
                        })

            # Peptide bond & omega (from prev residue i-1 -> i)
            phi_val: Optional[float] = None
            psi_val: Optional[float] = None
            omega_val: Optional[float] = None

            if i > 0:
                prev_res = pres[i - 1]
                prev_c = _get_pos(prev_res, "C")
                prev_ca = _get_pos(prev_res, "CA")
                prev_o = _get_pos(prev_res, "O")
                if prev_c is not None and n_p is not None and distance(prev_c, n_p) < 2.2:
                    # Peptide bond length
                    d_pep = distance(prev_c, n_p)
                    ideal_p, sig_p = IDEAL_BACKBONE_BONDS[("C_prev", "N")]
                    zp = (d_pep - ideal_p) / sig_p
                    bond_z_sq.append(zp * zp)
                    if abs(zp) > 4.0:
                        issues.append("bond_outlier")
                        bond_outliers.append({
                            "chain": chain_id,
                            "resseq": res.seqid.num,
                            "icode": res.seqid.icode.strip(),
                            "resname": rname,
                            "key": key,
                            "atoms": "C(prev)-N",
                            "observed": round(d_pep, 3),
                            "ideal": ideal_p,
                            "z_score": round(zp, 2),
                        })

                    if prev_ca is not None and ca_p is not None:
                        omega_val = round(dihedral_deg(prev_ca, prev_c, n_p, ca_p), 1)
                        abs_om = abs(omega_val)
                        if abs_om < 30.0:
                            is_pro = std == "PRO"
                            if not is_pro:
                                issues.append("cis_nonpro")
                            cis_peptides.append({
                                "chain": chain_id,
                                "resseq": res.seqid.num,
                                "icode": res.seqid.icode.strip(),
                                "prev_resseq": prev_res.seqid.num,
                                "prev_resname": prev_res.name.strip(),
                                "resname": rname,
                                "key": key,
                                "omega": omega_val,
                                "type": "cis-proline" if is_pro else "cis-non-proline",
                                "severity": "info" if is_pro else "outlier",
                            })
                        elif abs_om < 150.0:
                            issues.append("twisted_peptide")
                            cis_peptides.append({
                                "chain": chain_id,
                                "resseq": res.seqid.num,
                                "icode": res.seqid.icode.strip(),
                                "prev_resseq": prev_res.seqid.num,
                                "prev_resname": prev_res.name.strip(),
                                "resname": rname,
                                "key": key,
                                "omega": omega_val,
                                "type": "twisted-peptide",
                                "severity": "outlier",
                            })

                    if ca_p is not None and c_p is not None:
                        phi_val = round(dihedral_deg(prev_c, n_p, ca_p, c_p), 1)

            if i + 1 < len(pres):
                next_res = pres[i + 1]
                next_n = _get_pos(next_res, "N")
                if c_p is not None and next_n is not None and distance(c_p, next_n) < 2.2:
                    if n_p is not None and ca_p is not None:
                        psi_val = round(dihedral_deg(n_p, ca_p, c_p, next_n), 1)

            # Ramachandran classification
            rama_status = None
            if phi_val is not None and psi_val is not None:
                next_rname = pres[i + 1].name.strip() if i + 1 < len(pres) else None
                cat = rama_category(rname, next_rname)
                rama_status = evaluate_ramachandran(phi_val, psi_val, cat)
                rama_counts[rama_status] += 1
                if rama_status == "outlier":
                    issues.append("rama_outlier")
                rama_points.append({
                    "chain": chain_id,
                    "resseq": res.seqid.num,
                    "icode": res.seqid.icode.strip(),
                    "resname": rname,
                    "key": key,
                    "phi": phi_val,
                    "psi": psi_val,
                    "omega": omega_val,
                    "category": cat,
                    "status": rama_status,
                })

            # Chirality & C-beta deviation
            if std != "GLY" and n_p is not None and ca_p is not None and c_p is not None and cb_p is not None:
                # L-amino acid chiral volume: (N - CA) . ((C - CA) x (CB - CA)) is ~ +2.45 Å^3 for L-aa
                vol = chiral_volume(ca_p, n_p, c_p, cb_p)
                if vol < 1.2:
                    issues.append("chiral_outlier")
                    chiral_outliers.append({
                        "chain": chain_id,
                        "resseq": res.seqid.num,
                        "icode": res.seqid.icode.strip(),
                        "resname": rname,
                        "key": key,
                        "chiral_volume": round(vol, 2),
                        "expected": "+2.45 (L-amino acid)",
                        "diagnosis": "D-amino acid inversion" if vol < 0 else "Planarized C-alpha center",
                    })
                ideal_cb = ideal_cbeta_position(n_p, ca_p, c_p)
                cb_dev = distance(cb_p, ideal_cb)
                if cb_dev > 0.25:
                    issues.append("cbeta_deviation")
                    cbeta_outliers.append({
                        "chain": chain_id,
                        "resseq": res.seqid.num,
                        "icode": res.seqid.icode.strip(),
                        "resname": rname,
                        "key": key,
                        "deviation": round(cb_dev, 3),
                    })

            # Rotamer analysis
            rot_info = _evaluate_rotamer(res)
            if rot_info is not None:
                rotamer_counts[rot_info["status"]] += 1
                if rot_info["status"] == "outlier":
                    issues.append("rotamer_outlier")
                rotamer_points.append({
                    "chain": chain_id,
                    "resseq": res.seqid.num,
                    "icode": res.seqid.icode.strip(),
                    "resname": rname,
                    "key": key,
                    **rot_info,
                })

            if res_clashes.get(key, 0) > 0:
                issues.append("clash")
            if rname in NONSTANDARD_TO_STANDARD:
                issues.append("nonstandard")

            residues_index.append({
                "chain": chain_id,
                "resseq": res.seqid.num,
                "icode": res.seqid.icode.strip(),
                "resname": rname,
                "std_resname": std,
                "key": key,
                "phi": phi_val,
                "psi": psi_val,
                "omega": omega_val,
                "rama_status": rama_status,
                "rotamer_status": rot_info["status"] if rot_info else None,
                "chis": rot_info["chis"] if rot_info else [],
                "clash_count": res_clashes.get(key, 0),
                "missing_atoms": missing_atoms,
                "altlocs": altlocs,
                "mean_b": mean_b,
                "min_occ": min_occ,
                "issues": sorted(set(issues)),
            })

    total_rama = max(sum(rama_counts.values()), 1)
    total_rot = max(sum(rotamer_counts.values()), 1)
    bond_rmsz = round(math.sqrt(sum(bond_z_sq) / max(len(bond_z_sq), 1)), 3) if bond_z_sq else 0.0
    angle_rmsz = round(math.sqrt(sum(angle_z_sq) / max(len(angle_z_sq), 1)), 3) if angle_z_sq else 0.0
    nqh_flips = _detect_nqh_flips(st)

    # MolProbity-like composite crystallographic quality score (0-100)
    rama_fav_pct = round(100.0 * rama_counts["favored"] / total_rama, 1)
    rama_out_pct = round(100.0 * rama_counts["outlier"] / total_rama, 1)
    rot_fav_pct = round(100.0 * rotamer_counts["favored"] / total_rot, 1)
    rot_out_pct = round(100.0 * rotamer_counts["outlier"] / total_rot, 1)

    quality_score = 100.0
    quality_score -= min(28.0, clash_data["clashscore"] * 1.1)
    quality_score -= min(22.0, rama_out_pct * 3.0)
    quality_score -= min(18.0, rot_out_pct * 2.0)
    quality_score -= min(12.0, len(bond_outliers) * 1.5 + len(angle_outliers) * 0.8)
    quality_score -= min(10.0, sum(1 for c in cis_peptides if c["severity"] == "outlier") * 4.0)
    quality_score -= min(10.0, len(chiral_outliers) * 5.0)
    quality_score = round(max(0.0, min(100.0, quality_score)), 1)

    return {
        "quality_score": quality_score,
        "ramachandran": {
            "favored_pct": rama_fav_pct,
            "allowed_pct": round(100.0 * rama_counts["allowed"] / total_rama, 1),
            "outlier_pct": rama_out_pct,
            "counts": rama_counts,
            "points": rama_points,
        },
        "rotamers": {
            "favored_pct": rot_fav_pct,
            "allowed_pct": round(100.0 * rotamer_counts["allowed"] / total_rot, 1),
            "outlier_pct": rot_out_pct,
            "counts": rotamer_counts,
            "outliers": [r for r in rotamer_points if r["status"] == "outlier"],
            "all": rotamer_points,
        },
        "geometry": {
            "bond_rmsz": bond_rmsz,
            "angle_rmsz": angle_rmsz,
            "bond_outliers": sorted(bond_outliers, key=lambda x: abs(x["z_score"]), reverse=True)[:80],
            "angle_outliers": sorted(angle_outliers, key=lambda x: abs(x["z_score"]), reverse=True)[:80],
            "cbeta_outliers": cbeta_outliers,
            "cis_peptides": cis_peptides,
            "chiral_outliers": chiral_outliers,
        },
        "clashes": clash_data,
        "nqh_flips": nqh_flips,
        "bfactor_profile": bfactor_profile,
        "residues": residues_index,
    }


def _build_env_kdtree(st: gemmi.Structure, exclude_key: str) -> Tuple[Optional[cKDTree], List[float]]:
    pts = []
    radii = []
    for chain in st[0]:
        for res in chain:
            if _res_key(chain.name, res) == exclude_key:
                continue
            for atom in res:
                if atom.element.name in ("H", "D"):
                    continue
                pts.append([atom.pos.x, atom.pos.y, atom.pos.z])
                radii.append(VDW_RADII.get(atom.element.name.upper(), 1.70))
    if not pts:
        return None, []
    return cKDTree(np.array(pts, dtype=np.float64)), radii


def _score_sidechain_clashes(
    sc_coords: Dict[str, np.ndarray], tree: Optional[cKDTree], env_radii: List[float]
) -> Tuple[int, float]:
    if tree is None:
        return 0, 0.0
    n_clashes = 0
    penalty = 0.0
    for aname, pos in sc_coords.items():
        if aname in ("N", "CA", "C", "O", "CB"):
            continue
        r1 = VDW_RADII.get(aname[0].upper(), 1.70)
        neighbors = tree.query_ball_point(pos, r1 + 1.9)
        for ni in neighbors:
            d = float(np.linalg.norm(tree.data[ni] - pos))
            overlap = (r1 + env_radii[ni]) - d - 0.35
            if overlap > 0.25:
                n_clashes += 1
                penalty += overlap * overlap * 10.0
    return n_clashes, round(penalty, 2)


def get_residue_rotamers(st: gemmi.Structure, chain_id: str, resseq: int, icode: str = "") -> Dict[str, Any]:
    """Return ranked Dunbrack/Richardson rotamers for a residue with live clash scoring."""
    _, res, _ = _find_residue(st, chain_id, resseq, icode)
    if res is None:
        raise ValueError(f"Residue {chain_id}:{resseq}{icode} not found")

    std = canonical_resname(res.name)
    lib = ROTAMER_LIBRARY.get(std, [])
    bb = {}
    for bn in ("N", "CA", "C", "O"):
        p = _get_pos(res, bn)
        if p is not None:
            bb[bn] = p
    if not all(k in bb for k in ("N", "CA", "C")):
        raise ValueError("Residue lacks complete backbone N/CA/C atoms for rotamer modeling")

    tree, radii = _build_env_kdtree(st, _res_key(chain_id, res))
    current = _evaluate_rotamer(res)
    options = []
    for idx, (rname, prob, chis) in enumerate(lib):
        sc = build_sidechain_coords(bb, std, chi_targets=chis)
        clashes, energy = _score_sidechain_clashes(sc, tree, radii)
        options.append({
            "index": idx,
            "name": rname,
            "probability": prob,
            "chis": chis,
            "clashes": clashes,
            "clash_energy": energy,
            "score": round(prob - energy * 4.0 - clashes * 12.0, 2),
        })
    options.sort(key=lambda x: x["score"], reverse=True)
    return {
        "key": _res_key(chain_id, res),
        "resname": res.name.strip(),
        "std_resname": std,
        "current": current,
        "rotamers": options,
    }


def _replace_sidechain(res: gemmi.Residue, target_resname: str, new_coords: Dict[str, np.ndarray]) -> None:
    """Replace side-chain atoms of a Gemmi Residue while preserving backbone & B-factor."""
    bb_atoms = {}
    mean_b = 25.0
    bvals = [a.b_iso for a in res]
    if bvals:
        mean_b = float(np.mean(bvals))

    for a in list(res):
        if a.name.strip() in ("N", "CA", "C", "O", "OXT") and a.altloc in ("\x00", " ", "A"):
            if a.name.strip() not in bb_atoms:
                a.altloc = "\x00"
                a.occ = 1.0
                bb_atoms[a.name.strip()] = a.clone()

    # Clear all atoms and re-add backbone + built sidechain
    while len(res) > 0:
        del res[0]
    res.name = target_resname.upper()

    for bn in ("N", "CA", "C", "O", "OXT"):
        if bn in bb_atoms:
            res.add_atom(bb_atoms[bn])

    for aname, pos in new_coords.items():
        if aname in ("N", "CA", "C", "O", "OXT"):
            continue
        atom = gemmi.Atom()
        atom.name = aname
        atom.element = gemmi.Element(aname[0].upper())
        atom.pos = gemmi.Position(float(pos[0]), float(pos[1]), float(pos[2]))
        atom.occ = 1.0
        atom.b_iso = round(mean_b, 2)
        res.add_atom(atom)


def apply_rotamer(
    st: gemmi.Structure,
    chain_id: str,
    resseq: int,
    icode: str = "",
    rotamer_index: Optional[int] = None,
    chi_angles: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Apply a library rotamer (or auto-best if index is None) or explicit chi angles."""
    _, res, _ = _find_residue(st, chain_id, resseq, icode)
    if res is None:
        raise ValueError(f"Residue {chain_id}:{resseq}{icode} not found")

    std = canonical_resname(res.name)
    lib = ROTAMER_LIBRARY.get(std, [])
    bb = {bn: _get_pos(res, bn) for bn in ("N", "CA", "C", "O") if _get_pos(res, bn) is not None}
    if not all(k in bb for k in ("N", "CA", "C")):
        raise ValueError("Missing backbone atoms N/CA/C")

    if chi_angles is not None:
        target_chis = chi_angles
        chosen_name = "custom"
    elif rotamer_index is not None and 0 <= rotamer_index < len(lib):
        chosen_name, _, target_chis = lib[rotamer_index]
    else:
        # COOT Auto-Fit Rotamer: choose highest composite score
        ranked = get_residue_rotamers(st, chain_id, resseq, icode)["rotamers"]
        if not ranked:
            return {"message": f"{std} has no rotatable chi angles", "changed": False}
        best = ranked[0]
        chosen_name = best["name"]
        target_chis = best["chis"]

    new_coords = build_sidechain_coords(bb, std, chi_targets=target_chis)
    _replace_sidechain(res, res.name.strip(), new_coords)
    return {
        "changed": True,
        "key": _res_key(chain_id, res),
        "rotamer": chosen_name,
        "chis": target_chis,
        "message": f"Applied rotamer '{chosen_name}' ({target_chis}) to {res.name.strip()} {chain_id}:{resseq}",
    }


def mutate_residue(
    st: gemmi.Structure, chain_id: str, resseq: int, target_resname: str, icode: str = ""
) -> Dict[str, Any]:
    """COOT Simple Mutate & Auto-Fit Rotamer to any of the 20 standard amino acids."""
    target = target_resname.strip().upper()
    if target not in STANDARD_AA:
        raise ValueError(f"Target '{target}' is not a standard 20 amino acid")
    _, res, _ = _find_residue(st, chain_id, resseq, icode)
    if res is None:
        raise ValueError(f"Residue {chain_id}:{resseq}{icode} not found")

    old_name = res.name.strip()
    bb = {bn: _get_pos(res, bn) for bn in ("N", "CA", "C", "O") if _get_pos(res, bn) is not None}
    if not all(k in bb for k in ("N", "CA", "C")):
        raise ValueError("Residue lacks complete backbone atoms")

    lib = ROTAMER_LIBRARY.get(target, [])
    tree, radii = _build_env_kdtree(st, _res_key(chain_id, res))
    best_coords = build_sidechain_coords(bb, target)
    best_rot = "default"

    if lib:
        best_score = -1e9
        for rname, prob, chis in lib:
            sc = build_sidechain_coords(bb, target, chi_targets=chis)
            clashes, energy = _score_sidechain_clashes(sc, tree, radii)
            score = prob - energy * 4.0 - clashes * 12.0
            if score > best_score:
                best_score = score
                best_coords = sc
                best_rot = rname

    _replace_sidechain(res, target, best_coords)
    return {
        "changed": True,
        "key": _res_key(chain_id, res),
        "from_resname": old_name,
        "to_resname": target,
        "rotamer": best_rot,
        "message": f"Mutated {chain_id}:{resseq} {old_name} → {target} (auto-fit rotamer '{best_rot}')",
    }


def pepflip_residue(st: gemmi.Structure, chain_id: str, resseq: int, icode: str = "") -> Dict[str, Any]:
    """COOT Pepflip: rotate peptide carbonyl C=O and next N by 180° about the CA(i)-CA(i+1) axis."""
    chain, res, idx = _find_residue(st, chain_id, resseq, icode)
    if chain is None or res is None:
        raise ValueError(f"Residue {chain_id}:{resseq}{icode} not found")
    if idx + 1 >= len(chain) or not _is_protein_res(chain[idx + 1]):
        raise ValueError("Pepflip requires a covalently bonded following residue (i+1)")

    next_res = chain[idx + 1]
    ca1 = _get_pos(res, "CA")
    ca2 = _get_pos(next_res, "CA")
    if ca1 is None or ca2 is None or distance(ca1, ca2) > 4.5:
        raise ValueError("Adjacent CA atoms are missing or across a chain break")

    axis = ca2 - ca1
    for a in res:
        if a.name.strip() in ("C", "O"):
            p = rotate_point_around_axis(_atom_pos(a), ca1, axis, 180.0)
            a.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))
    for a in next_res:
        if a.name.strip() in ("N", "H", "HN"):
            p = rotate_point_around_axis(_atom_pos(a), ca1, axis, 180.0)
            a.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))

    return {
        "changed": True,
        "key": _res_key(chain_id, res),
        "message": f"Executed 180° COOT Pepflip on peptide {chain_id}:{resseq}-{next_res.seqid.num}",
    }


def flip_cis_trans(st: gemmi.Structure, chain_id: str, resseq: int, icode: str = "") -> Dict[str, Any]:
    """
    COOT Cis ↔ Trans conversion on peptide bond (i-1 → i):
    Rebuilds C(i-1), O(i-1), N(i) cleanly to invert omega (~0° ↔ ~180°) while keeping CA anchors fixed.
    """
    chain, res, idx = _find_residue(st, chain_id, resseq, icode)
    if chain is None or res is None or idx == 0:
        raise ValueError(f"Residue {chain_id}:{resseq}{icode} has no preceding peptide bond")
    prev_res = chain[idx - 1]
    if not _is_protein_res(prev_res):
        raise ValueError("Preceding residue is not a peptide residue")

    ca_prev = _get_pos(prev_res, "CA")
    c_prev = _get_pos(prev_res, "C")
    n_curr = _get_pos(res, "N")
    ca_curr = _get_pos(res, "CA")
    if any(p is None for p in (ca_prev, c_prev, n_curr, ca_curr)):
        raise ValueError("Missing backbone atoms around peptide bond")

    old_omega = dihedral_deg(ca_prev, c_prev, n_curr, ca_curr)
    is_cis = abs(old_omega) < 90.0

    # Place C_prev and N_curr symmetrically on opposite sides (for trans) or same side (for cis) of CA-CA vector
    v_ca = ca_curr - ca_prev
    dist_ca = np.linalg.norm(v_ca)
    if dist_ca < 1e-3:
        raise ValueError("Degenerate CA-CA distance")
    u_x = v_ca / dist_ca

    # Compute current plane normal from CA_prev -> C_prev
    v_c = c_prev - ca_prev
    perp = v_c - np.dot(v_c, u_x) * u_x
    if np.linalg.norm(perp) < 1e-3:
        perp = np.array([0.0, 1.0, 0.0])
    u_y = perp / np.linalg.norm(perp)

    if is_cis:
        # Convert CIS -> TRANS: place C_prev above axis, N_curr below axis
        new_c = ca_prev + u_x * (dist_ca * 0.36) + u_y * 0.65
        new_n = ca_prev + u_x * (dist_ca * 0.64) - u_y * 0.65
        mid = 0.5 * (new_c + new_n)
        cn = new_n - new_c
        cn_hat = cn / (np.linalg.norm(cn) + 1e-8)
        new_c = mid - cn_hat * 0.665
        new_n = mid + cn_hat * 0.665
        new_o = nerf_place(new_n, ca_prev, new_c, 1.231, 120.8, 180.0)
        target_label = "trans"
    else:
        # Convert TRANS -> CIS: place both C_prev and N_curr on the same side of CA-CA
        new_c = ca_prev + u_x * (dist_ca * 0.34) + u_y * 0.95
        new_n = ca_prev + u_x * (dist_ca * 0.66) + u_y * 0.95
        mid = 0.5 * (new_c + new_n)
        cn = new_n - new_c
        cn_hat = cn / (np.linalg.norm(cn) + 1e-8)
        new_c = mid - cn_hat * 0.665
        new_n = mid + cn_hat * 0.665
        new_o = nerf_place(new_n, ca_prev, new_c, 1.231, 120.8, 180.0)
        target_label = "cis"

    for a in prev_res:
        if a.name.strip() == "C":
            a.pos = gemmi.Position(*[float(x) for x in new_c])
        elif a.name.strip() == "O":
            a.pos = gemmi.Position(*[float(x) for x in new_o])
    for a in res:
        if a.name.strip() == "N":
            a.pos = gemmi.Position(*[float(x) for x in new_n])

    new_omega = round(dihedral_deg(ca_prev, new_c, new_n, ca_curr), 1)
    return {
        "changed": True,
        "key": _res_key(chain_id, res),
        "old_omega": round(old_omega, 1),
        "new_omega": new_omega,
        "message": f"Converted peptide {prev_res.seqid.num}→{resseq} to {target_label} (ω: {round(old_omega, 1)}° → {new_omega}°)",
    }


def flip_nqh_residue(st: gemmi.Structure, chain_id: str, resseq: int, icode: str = "") -> Dict[str, Any]:
    """COOT / Reduce 180° flip for ASN (OD1↔ND2), GLN (OE1↔NE2), or HIS (imidazole ring 180°)."""
    _, res, _ = _find_residue(st, chain_id, resseq, icode)
    if res is None:
        raise ValueError(f"Residue {chain_id}:{resseq}{icode} not found")
    std = canonical_resname(res.name)
    if std == "ASN":
        a1, a2 = _get_atom(res, "OD1"), _get_atom(res, "ND2")
        if a1 is None or a2 is None:
            raise ValueError("ASN missing OD1 or ND2")
        p1, p2 = a1.pos, a2.pos
        a1.pos, a2.pos = gemmi.Position(p2.x, p2.y, p2.z), gemmi.Position(p1.x, p1.y, p1.z)
    elif std == "GLN":
        a1, a2 = _get_atom(res, "OE1"), _get_atom(res, "NE2")
        if a1 is None or a2 is None:
            raise ValueError("GLN missing OE1 or NE2")
        p1, p2 = a1.pos, a2.pos
        a1.pos, a2.pos = gemmi.Position(p2.x, p2.y, p2.z), gemmi.Position(p1.x, p1.y, p1.z)
    elif std == "HIS":
        cb = _get_pos(res, "CB")
        cg = _get_pos(res, "CG")
        if cb is None or cg is None:
            raise ValueError("HIS missing CB or CG")
        axis = cg - cb
        for a in res:
            if a.name.strip() in ("ND1", "CD2", "CE1", "NE2", "HD1", "HD2", "HE1", "HE2"):
                p = rotate_point_around_axis(_atom_pos(a), cg, axis, 180.0)
                a.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))
    else:
        raise ValueError(f"NQH 180° flip applies only to ASN, GLN, or HIS (got {std})")

    return {
        "changed": True,
        "key": _res_key(chain_id, res),
        "message": f"Flipped {std} {chain_id}:{resseq} terminal group by 180°",
    }


def pick_altloc(st: gemmi.Structure, chain_id: str, resseq: int, altloc: str, icode: str = "") -> Dict[str, Any]:
    """Resolve alternate conformations on a residue to a single chosen conformer."""
    _, res, _ = _find_residue(st, chain_id, resseq, icode)
    if res is None:
        raise ValueError(f"Residue {chain_id}:{resseq}{icode} not found")

    target_alt = altloc.strip().upper()
    kept_names = set()
    to_keep = []
    for a in res:
        al = a.altloc if a.altloc not in ("\x00", " ") else ""
        if al == "" or al.upper() == target_alt:
            if a.name.strip() not in kept_names:
                c = a.clone()
                c.altloc = "\x00"
                c.occ = 1.0
                to_keep.append(c)
                kept_names.add(a.name.strip())

    while len(res) > 0:
        del res[0]
    for a in to_keep:
        res.add_atom(a)

    return {
        "changed": True,
        "key": _res_key(chain_id, res),
        "message": f"Resolved {chain_id}:{resseq} alternate conformations → Conformer '{target_alt}' (occ=1.0)",
    }


def fix_chirality(st: gemmi.Structure, chain_id: str, resseq: int, icode: str = "") -> Dict[str, Any]:
    """Restore ideal L-amino acid C-alpha chirality and rebuild side chain with clash-free rotamer."""
    _, res, _ = _find_residue(st, chain_id, resseq, icode)
    if res is None:
        raise ValueError(f"Residue {chain_id}:{resseq}{icode} not found")
    std = canonical_resname(res.name)
    if std == "GLY":
        return {"changed": False, "message": "Glycine is achiral"}
    return mutate_residue(st, chain_id, resseq, std, icode)


def regularize_zone(
    st: gemmi.Structure, chain_id: str, start_resseq: int, end_resseq: int, iterations: int = 60
) -> Dict[str, Any]:
    """
    COOT 'Regularize Zone' / Local Geometry Idealization:
    Performs damped gradient descent on Engh & Huber backbone bond lengths (N-CA, CA-C, C-O, C-N),
    backbone bond angles, soft-sphere steric repulsion, and Ramachandran outlier steering,
    then re-anchors side chains cleanly onto refined backbones.
    """
    if len(st) == 0:
        raise ValueError("Empty structure")

    target_chain = next((c for c in st[0] if c.name == chain_id), None)
    if target_chain is None:
        raise ValueError(f"Chain {chain_id} not found")

    zone_res = [
        r for r in target_chain
        if _is_protein_res(r) and start_resseq <= r.seqid.num <= end_resseq
    ]
    if not zone_res:
        raise ValueError(f"No protein residues found in {chain_id}:{start_resseq}-{end_resseq}")

    # Collect backbone coordinates for the zone
    bb_names = ("N", "CA", "C", "O")
    atom_refs: List[Tuple[int, str, gemmi.Atom]] = []
    pos_list: List[np.ndarray] = []
    idx_map: Dict[Tuple[int, str], int] = {}

    for ri, res in enumerate(zone_res):
        for an in bb_names:
            a = _get_atom(res, an)
            if a is not None:
                idx = len(pos_list)
                atom_refs.append((ri, an, a))
                pos_list.append(_atom_pos(a))
                idx_map[(ri, an)] = idx

    if not pos_list:
        return {"changed": False, "message": "No backbone atoms to regularize"}

    X = np.array(pos_list, dtype=np.float64)
    X0 = X.copy()

    # Build bond restraint list: (i, j, ideal_len, weight)
    bonds: List[Tuple[int, int, float, float]] = []
    for ri in range(len(zone_res)):
        for a1, a2, d0, w in [
            ("N", "CA", 1.458, 12.0),
            ("CA", "C", 1.525, 12.0),
            ("C", "O", 1.231, 16.0),
        ]:
            if (ri, a1) in idx_map and (ri, a2) in idx_map:
                bonds.append((idx_map[(ri, a1)], idx_map[(ri, a2)], d0, w))
        if ri + 1 < len(zone_res):
            if (ri, "C") in idx_map and (ri + 1, "N") in idx_map:
                i1, i2 = idx_map[(ri, "C")], idx_map[(ri + 1, "N")]
                if distance(X[i1], X[i2]) < 2.2:
                    bonds.append((i1, i2, 1.329, 14.0))
        # Pseudo-angle restraints via 1-3 ideal distances (fast, unconditionally stable)
        for a1, a2, d13, w in [
            ("N", "C", 2.462, 5.0),
            ("CA", "O", 2.400, 6.0),
        ]:
            if (ri, a1) in idx_map and (ri, a2) in idx_map:
                bonds.append((idx_map[(ri, a1)], idx_map[(ri, a2)], d13, w))
        if ri + 1 < len(zone_res):
            if (ri, "CA") in idx_map and (ri + 1, "N") in idx_map and (ri, "C") in idx_map:
                if distance(X[idx_map[(ri, "C")]], X[idx_map[(ri + 1, "N")]]) < 2.2:
                    bonds.append((idx_map[(ri, "CA")], idx_map[(ri + 1, "N")], 2.425, 5.0))
                    if (ri, "O") in idx_map:
                        bonds.append((idx_map[(ri, "O")], idx_map[(ri + 1, "N")], 2.250, 5.0))
                    if (ri + 1, "CA") in idx_map:
                        bonds.append((idx_map[(ri, "C")], idx_map[(ri + 1, "CA")], 2.435, 5.0))

    # Iterative harmonic relaxation with light anchor tether so >4σ bond/angle outliers actually resolve
    step = 0.018
    for _ in range(iterations):
        grad = np.zeros_like(X)
        # Light tether preserves global fold while allowing Engh & Huber outliers to relax cleanly
        grad += 0.06 * (X - X0)

        for i1, i2, d0, w in bonds:
            diff = X[i2] - X[i1]
            d = float(np.linalg.norm(diff))
            if d < 1e-6:
                continue
            f = (w * 1.6) * (d - d0) / d * diff
            grad[i1] -= f
            grad[i2] += f

        # Clamp max displacement per iteration for stability
        norms = np.linalg.norm(grad, axis=1, keepdims=True)
        scale = np.minimum(1.0, 0.10 / (norms * step + 1e-9))
        X -= step * grad * scale

    # Write regularized backbone coordinates back and shift side-chain atoms with CA displacement
    ca_shifts: Dict[int, np.ndarray] = {}
    for idx, (ri, an, atom) in enumerate(atom_refs):
        if an == "CA":
            ca_shifts[ri] = X[idx] - X0[idx]
        atom.pos = gemmi.Position(float(X[idx, 0]), float(X[idx, 1]), float(X[idx, 2]))

    for ri, res in enumerate(zone_res):
        shift = ca_shifts.get(ri)
        if shift is not None and np.linalg.norm(shift) > 1e-5:
            for a in res:
                if a.name.strip() not in bb_names:
                    p = _atom_pos(a) + shift
                    a.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))
        # Also snap CB back to ideal L-amino acid position if C-beta deviation > 0.15 Å or CA-CB bond > 4σ
        std = canonical_resname(res.name)
        if std != "GLY":
            n_p, ca_p, c_p, cb_a = _get_pos(res, "N"), _get_pos(res, "CA"), _get_pos(res, "C"), _get_atom(res, "CB")
            if all(p is not None for p in (n_p, ca_p, c_p)) and cb_a is not None:
                ideal_cb = ideal_cbeta_position(n_p, ca_p, c_p)
                cb_cur = _atom_pos(cb_a)
                if distance(cb_cur, ideal_cb) > 0.15 or abs(distance(ca_p, cb_cur) - 1.530) > 0.06:
                    cb_shift = ideal_cb - cb_cur
                    for a in res:
                        if a.name.strip() not in ("N", "CA", "C", "O", "OXT", "H", "HN", "HA"):
                            p = _atom_pos(a) + cb_shift
                            a.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))

    rmsd_moved = round(float(np.sqrt(np.mean(np.sum((X - X0) ** 2, axis=1)))), 3)
    return {
        "changed": True,
        "zone": f"{chain_id}:{start_resseq}-{end_resseq}",
        "residues_regularized": len(zone_res),
        "rmsd_shift_A": rmsd_moved,
        "message": f"Regularized {len(zone_res)} residues in {chain_id}:{start_resseq}-{end_resseq} (RMSD shift {rmsd_moved} Å)",
    }


def _fix_rama_outlier_residue(st: gemmi.Structure, chain_id: str, resseq: int, icode: str = "") -> bool:
    """
    COOT peptide-plane rotation to move a Ramachandran outlier into a favored/allowed basin:
    Rotates peptide plane (C_prev, O_prev, N_i) around CA_prev-CA_i and/or
    (C_i, O_i, N_next) around CA_i-CA_next to achieve valid (phi, psi) with minimum shift.
    """
    chain, res, idx = _find_residue(st, chain_id, resseq, icode)
    if chain is None or res is None or idx <= 0 or idx + 1 >= len(chain):
        return False
    prev_res = chain[idx - 1]
    next_res = chain[idx + 1]
    if not (_is_protein_res(prev_res) and _is_protein_res(next_res)):
        return False

    ca_prev = _get_pos(prev_res, "CA")
    c_prev = _get_pos(prev_res, "C")
    o_prev = _get_pos(prev_res, "O")
    n_i = _get_pos(res, "N")
    ca_i = _get_pos(res, "CA")
    c_i = _get_pos(res, "C")
    o_i = _get_pos(res, "O")
    n_next = _get_pos(next_res, "N")
    ca_next = _get_pos(next_res, "CA")
    if any(p is None for p in (ca_prev, c_prev, n_i, ca_i, c_i, n_next, ca_next)):
        return False
    if distance(c_prev, n_i) > 2.2 or distance(c_i, n_next) > 2.2:
        return False

    cat = rama_category(res.name, next_res.name)
    cat_prev = rama_category(prev_res.name, res.name)
    next2_name = chain[idx + 2].name if (idx + 2 < len(chain) and _is_protein_res(chain[idx + 2])) else None
    cat_next = rama_category(next_res.name, next2_name)
    prev2_c = _get_pos(chain[idx - 2], "C") if (idx - 2 >= 0 and _is_protein_res(chain[idx - 2])) else None
    prev_n = _get_pos(prev_res, "N")
    next_c = _get_pos(next_res, "C")
    next2_n = _get_pos(chain[idx + 2], "N") if (idx + 2 < len(chain) and _is_protein_res(chain[idx + 2])) else None

    axis1 = ca_i - ca_prev
    axis2 = ca_next - ca_i

    best_angles = None
    best_cost = 1e9
    # Coarse-to-fine search over peptide plane rotations keeping CA fixed and protecting neighbor (i-1, i+1) Rama & angles
    candidates = [0, 15, -15, 30, -30, 50, -50, 75, -75, 105, -105, 140, -140, 180]
    for a1 in candidates:
        c_p_rot = rotate_point_around_axis(c_prev, ca_prev, axis1, a1)
        n_i_rot = rotate_point_around_axis(n_i, ca_prev, axis1, a1)
        # Protect prev_res (i-1) Ramachandran & N-CA-C angle
        if prev_n is not None:
            if abs(angle_deg(prev_n, ca_prev, c_p_rot) - 111.2) > 9.0:
                continue
            if prev2_c is not None:
                phi_prev = dihedral_deg(prev2_c, prev_n, ca_prev, c_p_rot)
                psi_prev = dihedral_deg(prev_n, ca_prev, c_p_rot, n_i_rot)
                if evaluate_ramachandran(phi_prev, psi_prev, cat_prev) == "outlier":
                    continue
        for a2 in candidates:
            c_i_rot = rotate_point_around_axis(c_i, ca_i, axis2, a2)
            n_n_rot = rotate_point_around_axis(n_next, ca_i, axis2, a2)

            phi = dihedral_deg(c_p_rot, n_i_rot, ca_i, c_i_rot)
            psi = dihedral_deg(n_i_rot, ca_i, c_i_rot, n_n_rot)
            status = evaluate_ramachandran(phi, psi, cat)
            if status == "outlier":
                continue
            # Keep N-CA-C angle strictly within Engh & Huber (< 3.5σ)
            ang_nac = angle_deg(n_i_rot, ca_i, c_i_rot)
            if abs(ang_nac - 111.2) > 9.0:
                continue
            # Protect next_res (i+1) Ramachandran & N-CA-C angle
            if next_c is not None:
                if abs(angle_deg(n_n_rot, ca_next, next_c) - 111.2) > 9.0:
                    continue
                if next2_n is not None:
                    phi_nxt = dihedral_deg(c_i_rot, n_n_rot, ca_next, next_c)
                    psi_nxt = dihedral_deg(n_n_rot, ca_next, next_c, next2_n)
                    if evaluate_ramachandran(phi_nxt, psi_nxt, cat_next) == "outlier":
                        continue
            cost = abs(a1) + abs(a2) + (0.0 if status == "favored" else 45.0) + abs(ang_nac - 111.2) * 3.0
            if cost < best_cost:
                best_cost = cost
                best_angles = (a1, a2)

    if best_angles is None:
        return False

    a1, a2 = best_angles
    if a1 != 0:
        for a in prev_res:
            if a.name.strip() in ("C", "O"):
                p = rotate_point_around_axis(_atom_pos(a), ca_prev, axis1, a1)
                a.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))
        for a in res:
            if a.name.strip() in ("N", "H", "HN"):
                p = rotate_point_around_axis(_atom_pos(a), ca_prev, axis1, a1)
                a.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))
    if a2 != 0:
        for a in res:
            if a.name.strip() in ("C", "O"):
                p = rotate_point_around_axis(_atom_pos(a), ca_i, axis2, a2)
                a.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))
        for a in next_res:
            if a.name.strip() in ("N", "H", "HN"):
                p = rotate_point_around_axis(_atom_pos(a), ca_i, axis2, a2)
                a.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))

    # Re-anchor CB/sidechain to ideal tetrahedral geometry after backbone plane rotation
    std = canonical_resname(res.name)
    if std != "GLY":
        n_new, ca_new, c_new, cb_a = _get_pos(res, "N"), _get_pos(res, "CA"), _get_pos(res, "C"), _get_atom(res, "CB")
        if all(p is not None for p in (n_new, ca_new, c_new)) and cb_a is not None:
            ideal_cb = ideal_cbeta_position(n_new, ca_new, c_new)
            shift = ideal_cb - _atom_pos(cb_a)
            for a in res:
                if a.name.strip() not in ("N", "CA", "C", "O", "OXT", "H", "HN", "HA"):
                    p = _atom_pos(a) + shift
                    a.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))
    return True


def _relieve_steric_clashes(st: gemmi.Structure, max_passes: int = 18) -> int:
    """
    Fast vectorized KDTree steric clash relief:
    Gently separates overlapping non-bonded atom pairs (preferring side-chain / outer atoms
    over backbone anchors) until steric overlaps drop below the 0.40 Å MolProbity threshold.
    """
    if len(st) == 0:
        return 0

    atom_objs: List[gemmi.Atom] = []
    coords: List[List[float]] = []
    meta = []
    res_global = -1
    prev_uid = None
    bb_set = {"N", "CA", "C", "O", "OXT", "CB"}

    for ch in st[0]:
        for res in ch:
            uid = (ch.name, res.seqid.num, res.seqid.icode.strip())
            if uid != prev_uid:
                res_global += 1
                prev_uid = uid
            for a in res:
                el = a.element.name.upper()
                if el not in VDW_RADII:
                    el = a.name.strip()[0].upper() if a.name.strip() else "C"
                atom_objs.append(a)
                coords.append([a.pos.x, a.pos.y, a.pos.z])
                an = a.name.strip()
                meta.append((
                    ch.name,
                    res.seqid.num,
                    an,
                    el,
                    res_global,
                    a.altloc if a.altloc != "\x00" else "",
                    el in ("N", "O", "F"),
                    an in bb_set,
                ))

    n_atoms = len(coords)
    if n_atoms < 2:
        return 0

    X = np.array(coords, dtype=np.float64)
    total_moved_pairs = 0

    for _ in range(max_passes):
        tree = cKDTree(X)
        pairs = tree.query_pairs(r=3.2)
        disp = np.zeros_like(X)
        counts = np.zeros(n_atoms, dtype=np.float64)
        n_clashing = 0

        for i, j in pairs:
            m1, m2 = meta[i], meta[j]
            if m1[4] == m2[4]:
                continue
            if m1[5] and m2[5] and m1[5] != m2[5]:
                continue
            if m1[0] == m2[0] and abs(m1[1] - m2[1]) == 1:
                bb_ext = {"N", "CA", "C", "O", "H", "HN", "HA", "HA2", "HA3", "CB", "CD"}
                if (m1[2] in bb_ext and m2[2] in bb_ext) or m1[3] == "H" or m2[3] == "H":
                    continue
            if m1[2] == "SG" and m2[2] == "SG":
                continue
            if m1[3] in ("ZN", "MG", "CA", "FE", "MN", "CU", "CO", "NI") or m2[3] in ("ZN", "MG", "CA", "FE", "MN", "CU", "CO", "NI"):
                continue

            vec = X[j] - X[i]
            d = float(np.linalg.norm(vec))
            if d < 1.80 and m1[3] != "H" and m2[3] != "H":
                continue
            if (m1[3] == "H" and m2[3] in ("N", "O", "S")) or (m2[3] == "H" and m1[3] in ("N", "O", "S")):
                if d >= 1.50:
                    continue

            r1 = 1.05 if m1[3] == "H" else VDW_RADII.get(m1[3], 1.70)
            r2 = 1.05 if m2[3] == "H" else VDW_RADII.get(m2[3], 1.70)
            hb_allow = 0.65 if (m1[6] and m2[6]) else 0.0
            overlap = (r1 + r2) - d - hb_allow

            if overlap > 0.36:
                n_clashing += 1
                u = vec / (d + 1e-6) if d > 1e-5 else np.array([1.0, 0.0, 0.0])
                push = min(0.35, (overlap - 0.30) * 0.65)
                # Keep rigid backbone + CB intact whenever the partner atom can absorb the push
                if m1[7] and not m2[7]:
                    w1, w2 = 0.0, 1.0
                elif m2[7] and not m1[7]:
                    w1, w2 = 1.0, 0.0
                else:
                    w1, w2 = 0.5, 0.5
                wsum = w1 + w2
                disp[i] -= u * push * (w1 / wsum)
                disp[j] += u * push * (w2 / wsum)
                counts[i] += 1.0
                counts[j] += 1.0

        if n_clashing == 0:
            break
        total_moved_pairs += n_clashing
        mask = counts > 0
        X[mask] += disp[mask] / np.maximum(counts[mask, None], 1.0)

    if total_moved_pairs > 0:
        for idx, a in enumerate(atom_objs):
            a.pos = gemmi.Position(float(X[idx, 0]), float(X[idx, 1]), float(X[idx, 2]))
    return total_moved_pairs


def delete_residue(st: gemmi.Structure, chain_id: str, resseq: int, icode: str = "") -> Dict[str, Any]:
    chain, res, idx = _find_residue(st, chain_id, resseq, icode)
    if chain is None or res is None:
        raise ValueError(f"Residue {chain_id}:{resseq}{icode} not found")
    rname = res.name.strip()
    del chain[idx]
    return {
        "changed": True,
        "message": f"Deleted residue {rname} {chain_id}:{resseq}{icode}",
    }


def _fast_autofit_rotamers(st: gemmi.Structure, targets: List[Tuple[str, int, str]]) -> int:
    """Fast batch rotamer auto-fit sharing a single KDTree across all target residues."""
    if not targets or len(st) == 0:
        return 0

    # Build one shared heavy-atom coordinate tree for fast O(log N) clash scoring
    pts = []
    radii = []
    res_keys = []
    for ch in st[0]:
        for res in ch:
            rk = _res_key(ch.name, res)
            for a in res:
                if a.element.name in ("H", "D"):
                    continue
                pts.append([a.pos.x, a.pos.y, a.pos.z])
                radii.append(VDW_RADII.get(a.element.name.upper(), 1.70))
                res_keys.append(rk)
    if not pts:
        return 0
    pts_arr = np.array(pts, dtype=np.float64)
    tree = cKDTree(pts_arr)

    fixed = 0
    for ch_id, rseq, ic in targets:
        _, res, _ = _find_residue(st, ch_id, rseq, ic)
        if res is None:
            continue
        std = canonical_resname(res.name)
        lib = ROTAMER_LIBRARY.get(std, [])
        if not lib:
            continue
        bb = {bn: _get_pos(res, bn) for bn in ("N", "CA", "C", "O") if _get_pos(res, bn) is not None}
        if not all(k in bb for k in ("N", "CA", "C")):
            continue

        self_key = _res_key(ch_id, res)
        best_score = -1e9
        best_coords = None

        for _, prob, chis in lib:
            sc = build_sidechain_coords(bb, std, chi_targets=chis)
            n_cl = 0
            pen = 0.0
            for aname, pos in sc.items():
                if aname in ("N", "CA", "C", "O", "CB"):
                    continue
                r1 = VDW_RADII.get(aname[0].upper(), 1.70)
                for ni in tree.query_ball_point(pos, r1 + 1.85):
                    if res_keys[ni] == self_key:
                        continue
                    d = float(np.linalg.norm(pts_arr[ni] - pos))
                    ov = (r1 + radii[ni]) - d - 0.35
                    if ov > 0.25:
                        n_cl += 1
                        pen += ov * ov * 12.0
            score = prob - pen * 4.0 - n_cl * 14.0
            if score > best_score:
                best_score = score
                best_coords = sc

        if best_coords is not None:
            _replace_sidechain(res, res.name.strip(), best_coords)
            fixed += 1
    return fixed


def auto_fix_coot(
    st: gemmi.Structure,
    fix_clashes: bool = True,
    fix_rotamers: bool = True,
    fix_nqh_flips: bool = True,
    fix_cis_nonpro: bool = True,
    fix_ramachandran: bool = True,
    regularize_outliers: bool = True,
) -> Dict[str, Any]:
    """
    Comprehensive One-Click COOT Auto-Refinement pipeline:
    Resolves structural blockers (AltLocs, non-standard residues, missing heavy atoms, chirality,
    C-beta deviations, cis-non-Pro/twisted peptides, NQH flips, ALL rotamer outliers,
    Ramachandran outliers, Engh & Huber bond/angle outliers, and all-atom steric clashes),
    then re-validates so resolved issues vanish from the COOT issue tracker.
    """
    actions: List[str] = []
    before = analyze_structure_quality(st)
    before_issues = sum(1 for r in before["residues"] if r["issues"])

    # 1. Resolve AltLocs, non-standard residues, and missing heavy atoms so those issue badges clear
    alt_resolved = 0
    nonstd_fixed = 0
    atoms_rebuilt = 0
    if len(st) > 0:
        for ch in st[0]:
            for res in ch:
                # Resolve multi-conformer AltLocs (keep highest occupancy)
                altlocs = sorted({a.altloc for a in res if a.altloc not in ("\x00", " ")})
                if len(altlocs) > 1:
                    chosen = max(
                        altlocs,
                        key=lambda al: float(np.mean([a.occ for a in res if a.altloc == al] or [0.0])),
                    )
                    pick_altloc(st, ch.name, res.seqid.num, chosen, res.seqid.icode.strip())
                    alt_resolved += 1

                # Convert non-standard residues (e.g. MSE -> MET)
                rn = res.name.strip().upper()
                if rn in NONSTANDARD_TO_STANDARD:
                    target = NONSTANDARD_TO_STANDARD[rn]
                    if rn == "MSE":
                        for a in res:
                            if a.name.strip() == "SE":
                                a.name = "SD"
                                a.element = gemmi.Element("S")
                    allowed = set(EXPECTED_HEAVY_ATOMS[target]) | {"OXT"}
                    for ai in range(len(res) - 1, -1, -1):
                        if res[ai].element.name not in ("H", "D") and res[ai].name.strip() not in allowed:
                            del res[ai]
                    res.name = target
                    res.het_flag = "A"
                    nonstd_fixed += 1

                # Rebuild missing heavy atoms on standard residues
                std = canonical_resname(res.name)
                if std in EXPECTED_HEAVY_ATOMS:
                    present = {a.name.strip() for a in res if a.element.name not in ("H", "D")}
                    missing = [a for a in EXPECTED_HEAVY_ATOMS[std] if a not in present]
                    if missing:
                        existing_map = {
                            a.name.strip(): _atom_pos(a)
                            for a in res
                            if a.element.name not in ("H", "D")
                        }
                        if all(k in existing_map for k in ("N", "CA", "C")):
                            built = build_sidechain_coords(existing_map, std)
                            mean_b = float(np.mean([a.b_iso for a in res] or [25.0]))
                            for m in missing:
                                if m in built:
                                    atom = gemmi.Atom()
                                    atom.name = m
                                    atom.element = gemmi.Element(m[0].upper())
                                    atom.pos = gemmi.Position(*[float(x) for x in built[m]])
                                    atom.occ = 1.0
                                    atom.b_iso = round(mean_b, 2)
                                    res.add_atom(atom)
                                    atoms_rebuilt += 1

    if alt_resolved:
        actions.append(f"Resolved AltLoc conformers on {alt_resolved} residue(s)")
    if nonstd_fixed:
        actions.append(f"Converted {nonstd_fixed} non-standard residue(s) to canonical amino acids")
    if atoms_rebuilt:
        actions.append(f"Rebuilt {atoms_rebuilt} missing heavy atom(s)")

    # 2. Fix cis non-proline and twisted peptides
    if fix_cis_nonpro:
        cis_fixed = 0
        for cp in before["geometry"]["cis_peptides"]:
            if cp["severity"] == "outlier":
                try:
                    flip_cis_trans(st, cp["chain"], cp["resseq"], cp["icode"])
                    cis_fixed += 1
                except Exception:
                    pass
        if cis_fixed:
            actions.append(f"Converted {cis_fixed} cis-non-Pro / twisted peptide(s) to planar trans")

    # 3. Fix chirality & C-beta outliers
    chiral_fixed = 0
    for ch_o in before["geometry"]["chiral_outliers"]:
        try:
            fix_chirality(st, ch_o["chain"], ch_o["resseq"], ch_o["icode"])
            chiral_fixed += 1
        except Exception:
            pass
    for cb_o in before["geometry"]["cbeta_outliers"]:
        try:
            fix_chirality(st, cb_o["chain"], cb_o["resseq"], cb_o["icode"])
            chiral_fixed += 1
        except Exception:
            pass
    if chiral_fixed:
        actions.append(f"Restored L-chirality / ideal Cβ geometry on {chiral_fixed} residue(s)")

    # 4. Execute beneficial Asn/Gln/His 180° flips
    if fix_nqh_flips:
        nqh_fixed = 0
        for cand in before["nqh_flips"]:
            try:
                flip_nqh_residue(st, cand["chain"], cand["resseq"], cand["icode"])
                nqh_fixed += 1
            except Exception:
                pass
        if nqh_fixed:
            actions.append(f"Executed 180° H-bond flip on {nqh_fixed} Asn/Gln/His residue(s)")

    # 5. Fix Ramachandran outliers via neighbor-safe peptide-plane rotation (monotonic check!)
    if fix_ramachandran:
        rama_fixed = 0
        for rp in before["ramachandran"]["points"]:
            if rp["status"] == "outlier":
                try:
                    if _fix_rama_outlier_residue(st, rp["chain"], rp["resseq"], rp["icode"]):
                        rama_fixed += 1
                except Exception:
                    pass
        if rama_fixed:
            actions.append(f"Steered {rama_fixed} Ramachandran outlier(s) into favored/allowed basins")

    # 6. Auto-Fit ALL outlier rotamers & ALL clashing side chains (no 35-item cap!)
    targets_to_rotamer = set()
    if fix_rotamers:
        for r in before["rotamers"]["outliers"]:
            targets_to_rotamer.add((r["chain"], r["resseq"], r["icode"]))
    if fix_clashes:
        for r in before["residues"]:
            if r["clash_count"] > 0:
                targets_to_rotamer.add((r["chain"], r["resseq"], r["icode"]))

    rot_fixed = _fast_autofit_rotamers(st, sorted(targets_to_rotamer))
    if rot_fixed:
        actions.append(f"Auto-fit clash-free rotamers on {rot_fixed} residue(s)")

    # 7. Fast vectorized all-atom steric clash relief BEFORE final Engh & Huber geometry idealization
    if fix_clashes:
        relieved = _relieve_steric_clashes(st, max_passes=22)
        if relieved > 0:
            actions.append("Relieved residual all-atom steric overlaps via local relaxation")

    # 8. Final Engh & Huber bond/angle/C-beta zone regularization (with monotonic rollback guard per zone!)
    if regularize_outliers:
        mid_check = analyze_structure_quality(st)
        outlier_zones: Dict[str, List[int]] = {}
        for bo in (
            mid_check["geometry"]["bond_outliers"]
            + mid_check["geometry"]["angle_outliers"]
            + mid_check["geometry"]["cbeta_outliers"]
        ):
            outlier_zones.setdefault(bo["chain"], []).append(bo["resseq"])
        zones_count = 0
        for ch_id, rseqs in outlier_zones.items():
            uniq = sorted(set(rseqs))
            start = uniq[0]
            end = uniq[0]
            for num in uniq[1:] + [None]:
                if num is not None and num <= end + 3:
                    end = num
                else:
                    try:
                        snap = st.clone()
                        before_z = analyze_structure_quality(st)
                        before_z_cnt = sum(1 for r in before_z["residues"] if r["issues"])
                        regularize_zone(st, ch_id, max(1, start - 1), end + 1, iterations=120)
                        after_z = analyze_structure_quality(st)
                        after_z_cnt = sum(1 for r in after_z["residues"] if r["issues"])
                        if after_z_cnt > before_z_cnt and after_z["quality_score"] < before_z["quality_score"]:
                            # Roll back zone move that worsened geometry
                            while len(st) > 0:
                                del st[0]
                            for m in snap:
                                st.add_model(m.clone())
                        else:
                            zones_count += 1
                    except Exception:
                        pass
                    if num is not None:
                        start = end = num
        if zones_count:
            actions.append(f"Regularized {zones_count} Engh & Huber bond/angle/Cβ outlier zone(s)")

    # 9. Final convergence pass: fix any residual rotamers, C-beta/chirality, or clashes introduced
    post_check = analyze_structure_quality(st)
    if fix_rotamers and post_check["rotamers"]["outliers"]:
        extra_rot = [(r["chain"], r["resseq"], r["icode"]) for r in post_check["rotamers"]["outliers"]]
        _fast_autofit_rotamers(st, extra_rot)
    if fix_clashes and post_check["clashes"]["total_clashes"] > 0:
        _relieve_steric_clashes(st, max_passes=15)

    # Final validation pass after all refinements
    after = analyze_structure_quality(st)
    after_issues = sum(1 for r in after["residues"] if r["issues"])

    return {
        "actions": actions,
        "before_score": before["quality_score"],
        "after_score": after["quality_score"],
        "before_clashscore": before["clashes"]["clashscore"],
        "after_clashscore": after["clashes"]["clashscore"],
        "before_issues": before_issues,
        "after_issues": after_issues,
        "before_rotamers": before["rotamers"]["counts"]["outlier"],
        "after_rotamers": after["rotamers"]["counts"]["outlier"],
        "before_rama": before["ramachandran"]["counts"]["outlier"],
        "after_rama": after["ramachandran"]["counts"]["outlier"],
        "validation": after,
    }
