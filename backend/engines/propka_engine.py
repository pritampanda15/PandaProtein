"""
PROPKA 3 & Protonation / H-Bond Network Engine:
- Executes real PROPKA 3.5.1 pKa predictions on the current structure
- Computes full titration curve (pH 0.0 - 14.0), net charge, and Isoelectric Point (pI)
- Detects metal-coordinating residues and disulfide bridges automatically
- Intelligent Histidine tautomer selection (HID / HIE / HIP) via H-bond microenvironment
- Supports per-residue manual protonation overrides
- Builds accurate polar & nonpolar hydrogens respecting exact assigned protonation states
"""

from __future__ import annotations

import os
import tempfile
from typing import Any, Dict, List, Optional, Tuple

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .chemistry_data import ALL_PROTEIN_RESNAMES, METAL_IONS, MODEL_PKA, TITRATABLE_CHARGES
from .coot_engine import _get_pos, flip_nqh_residue, _detect_nqh_flips
from .geometry_utils import canonical_resname, distance, nerf_place
from .structure_io import structure_to_pdb

try:
    import logging
    import propka.run as propka_run
    logging.getLogger("propka").setLevel(logging.CRITICAL)
    HAS_PROPKA = True
except Exception:
    HAS_PROPKA = False


def _detect_metal_coordinating_residues(st: gemmi.Structure, cutoff: float = 2.85) -> Dict[str, List[str]]:
    """Map residue key 'A:18' -> list of coordinated metal descriptions (e.g. 'ZN:301 (2.08 Å)')."""
    if len(st) == 0:
        return {}
    metals: List[Tuple[str, np.ndarray]] = []
    for chain in st[0]:
        for res in chain:
            rname = res.name.strip().upper()
            if rname in METAL_IONS:
                for a in res:
                    metals.append((f"{rname} {chain.name}:{res.seqid.num}", np.array([a.pos.x, a.pos.y, a.pos.z])))
    if not metals:
        return {}

    coord_map: Dict[str, List[str]] = {}
    donor_atoms = {"SG", "ND1", "NE2", "OD1", "OD2", "OE1", "OE2", "OH"}
    for chain in st[0]:
        for res in chain:
            std = canonical_resname(res.name)
            if std not in ("CYS", "HIS", "ASP", "GLU", "TYR"):
                continue
            key = f"{chain.name}:{res.seqid.num}{res.seqid.icode.strip()}"
            for a in res:
                if a.name.strip() not in donor_atoms:
                    continue
                pos = np.array([a.pos.x, a.pos.y, a.pos.z])
                for m_label, m_pos in metals:
                    d = distance(pos, m_pos)
                    if d <= cutoff:
                        coord_map.setdefault(key, []).append(f"{m_label} via {a.name.strip()} ({d:.2f} Å)")
    return coord_map


def _detect_disulfide_keys(st: gemmi.Structure, cutoff: float = 2.55) -> set:
    if len(st) == 0:
        return set()
    cys = []
    for chain in st[0]:
        for res in chain:
            if canonical_resname(res.name) == "CYS":
                sg = _get_pos(res, "SG")
                if sg is not None:
                    cys.append((f"{chain.name}:{res.seqid.num}{res.seqid.icode.strip()}", res.name.strip().upper(), sg))
    bonded = set()
    for k, rn, _ in cys:
        if rn in ("CYX", "CYS2"):
            bonded.add(k)
    for i in range(len(cys)):
        for j in range(i + 1, len(cys)):
            if distance(cys[i][2], cys[j][2]) <= cutoff:
                bonded.add(cys[i][0])
                bonded.add(cys[j][0])
    return bonded


def _choose_his_tautomer(st: gemmi.Structure, chain_id: str, resseq: int, icode: str, metal_info: List[str]) -> Tuple[str, str]:
    """Select HID vs HIE based on metal coordination or H-bond acceptor proximity."""
    if metal_info:
        # If metal binds NE2, ND1 should be protonated (HID); if metal binds ND1, protonate NE2 (HIE)
        if any("via NE2" in m for m in metal_info):
            return "HID", "Metal coordinates NE2 → δ-protonated (HID)"
        if any("via ND1" in m for m in metal_info):
            return "HIE", "Metal coordinates ND1 → ε-protonated (HIE)"

    if len(st) == 0:
        return "HIE", "Default neutral ε-tautomer (HIE)"

    target_res = None
    acc_pts = []
    for chain in st[0]:
        for res in chain:
            is_self = (chain.name == chain_id and res.seqid.num == resseq and res.seqid.icode.strip() == icode)
            if is_self:
                target_res = res
                continue
            for a in res:
                if a.element.name == "O":
                    acc_pts.append([a.pos.x, a.pos.y, a.pos.z])

    if target_res is None or not acc_pts:
        return "HIE", "Standard physiological neutral ε-tautomer (HIE)"

    nd1 = _get_pos(target_res, "ND1")
    ne2 = _get_pos(target_res, "NE2")
    if nd1 is None or ne2 is None:
        return "HIE", "Standard neutral ε-tautomer (HIE)"

    tree = cKDTree(np.array(acc_pts))
    n_nd1 = len(tree.query_ball_point(nd1, 3.4))
    n_ne2 = len(tree.query_ball_point(ne2, 3.4))
    if n_nd1 > n_ne2:
        return "HID", f"ND1 has {n_nd1} H-bond acceptor(s) vs {n_ne2} at NE2 → HID"
    if n_ne2 > n_nd1:
        return "HIE", f"NE2 has {n_ne2} H-bond acceptor(s) vs {n_nd1} at ND1 → HIE"
    return "HIE", "Balanced microenvironment → physiological ε-tautomer (HIE)"


def _assign_state(
    resname: str,
    pka: float,
    ph: float,
    his_default: str,
    is_disulfide: bool,
    metal_info: List[str],
    respect_metal: bool,
    respect_disulfide: bool,
    tautomer_hint: Tuple[str, str],
) -> Tuple[str, float, str]:
    """Determine protonation state label, fractional charge at pH, and rationale."""
    std = canonical_resname(resname)
    if std in ("N+", "C-"):
        frac_prot = 1.0 / (1.0 + 10.0 ** (ph - pka))
        q_p, q_d = TITRATABLE_CHARGES[std]
        charge = q_p * frac_prot + q_d * (1.0 - frac_prot)
        state = ("NH3+" if ph < pka else "NH2") if std == "N+" else ("COOH" if ph < pka else "COO-")
        return state, round(charge, 3), f"Terminus pKa {pka:.2f} vs pH {ph:.1f}"

    if std == "CYS":
        if is_disulfide and respect_disulfide:
            return "CYX", 0.0, "Disulfide bridge (S–S covalent bond)"
        if metal_info and respect_metal:
            return "CYM", -1.0, f"Deprotonated thiolate coordinating metal ({metal_info[0]})"
        if ph > pka:
            return "CYM", -1.0, f"pH ({ph:.1f}) > pKa ({pka:.2f}) → deprotonated thiolate"
        return "CYS", 0.0, f"pH ({ph:.1f}) ≤ pKa ({pka:.2f}) → neutral thiol"

    if std == "HIS":
        if metal_info and respect_metal:
            return tautomer_hint[0], 0.0, tautomer_hint[1]
        if his_default in ("HID", "HIE", "HIP"):
            return his_default, (1.0 if his_default == "HIP" else 0.0), f"Global HIS preference ({his_default})"
        if ph < pka:
            return "HIP", 1.0, f"pH ({ph:.1f}) < pKa ({pka:.2f}) → doubly protonated (+1)"
        return tautomer_hint[0], 0.0, tautomer_hint[1]

    if std == "ASP":
        if ph < pka:
            return "ASH", 0.0, f"Shifted pKa ({pka:.2f}) > pH ({ph:.1f}) → protonated neutral ASH"
        return "ASP", -1.0, f"pH ({ph:.1f}) ≥ pKa ({pka:.2f}) → deprotonated carboxylate (-1)"

    if std == "GLU":
        if ph < pka:
            return "GLH", 0.0, f"Shifted pKa ({pka:.2f}) > pH ({ph:.1f}) → protonated neutral GLH"
        return "GLU", -1.0, f"pH ({ph:.1f}) ≥ pKa ({pka:.2f}) → deprotonated carboxylate (-1)"

    if std == "LYS":
        if ph > pka:
            return "LYN", 0.0, f"pH ({ph:.1f}) > pKa ({pka:.2f}) → neutral deprotonated LYN"
        return "LYS", 1.0, f"pH ({ph:.1f}) ≤ pKa ({pka:.2f}) → protonated ammonium (+1)"

    if std == "TYR":
        if ph > pka:
            return "TYM", -1.0, f"pH ({ph:.1f}) > pKa ({pka:.2f}) → phenolate (-1)"
        return "TYR", 0.0, f"pH ({ph:.1f}) ≤ pKa ({pka:.2f}) → neutral phenol"

    if std == "ARG":
        if ph > pka:
            return "ARN", 0.0, f"pH ({ph:.1f}) > pKa ({pka:.2f}) → neutral guanidine"
        return "ARG", 1.0, f"pH ({ph:.1f}) ≤ pKa ({pka:.2f}) → protonated guanidinium (+1)"

    return std, 0.0, "Non-titratable"


def _run_propka_raw(st: gemmi.Structure) -> Tuple[List[Dict[str, Any]], str]:
    """Run real PROPKA 3.5.1 on a clean copy of the structure; fallback to microenvironment model if needed."""
    # Prepare clean temporary PDB with canonical names so PROPKA recognizes all residues
    st_copy = st.clone()
    if len(st_copy) > 0:
        for ch in st_copy[0]:
            for res in ch:
                std = canonical_resname(res.name)
                if std in ALL_PROTEIN_RESNAMES:
                    res.name = std

    pdb_str = structure_to_pdb(st_copy)

    if HAS_PROPKA:
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tf:
                tf.write(pdb_str)
                tmp_path = tf.name
            mol = propka_run.single(tmp_path, optargs=["--quiet"], write_pka=False)
            conf_name = "AVR" if "AVR" in mol.conformations else next(iter(mol.conformations.keys()))
            conf = mol.conformations[conf_name]

            groups = []
            for g in conf.groups:
                if not getattr(g, "titratable", True):
                    continue
                res_type = g.residue_type.strip().upper()
                if res_type not in MODEL_PKA:
                    continue
                atom = g.atom
                chain_id = atom.chain_id.strip() or "A"
                resseq = int(atom.res_num)
                icode = (atom.icode or "").strip()
                pka_val = float(g.pka_value)
                model_pka = float(g.model_pka)
                buried_pct = round(float(getattr(g, "buried", 0.0)) * 100.0, 1)
                groups.append({
                    "chain": chain_id,
                    "resseq": resseq,
                    "icode": icode,
                    "resname": res_type,
                    "key": f"{chain_id}:{resseq}{icode}" if res_type not in ("N+", "C-") else f"{chain_id}:{resseq}{icode}:{res_type}",
                    "model_pka": round(model_pka, 2),
                    "pka": round(pka_val, 2),
                    "delta_pka": round(pka_val - model_pka, 2),
                    "buried_pct": min(100.0, max(0.0, buried_pct)),
                })
            if groups:
                return groups, "PROPKA 3.5.1"
        except Exception:
            pass
        finally:
            if tmp_path and os.path.exists(tmp_path):
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass

    # Microenvironment empirical fallback (burial desolvation + local Coulombic shift)
    groups = []
    if len(st) == 0:
        return groups, "Empirical Microenvironment Engine"
    all_pts = []
    for ch in st[0]:
        for r in ch:
            for a in r:
                if a.element.name not in ("H", "D"):
                    all_pts.append([a.pos.x, a.pos.y, a.pos.z])
    tree = cKDTree(np.array(all_pts)) if all_pts else None

    for ch in st[0]:
        for res in ch:
            std = canonical_resname(res.name)
            if std not in ("ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "ARG"):
                continue
            ref_atom = {"ASP": "CG", "GLU": "CD", "HIS": "CG", "CYS": "SG", "TYR": "OH", "LYS": "NZ", "ARG": "CZ"}[std]
            p = _get_pos(res, ref_atom)
            if p is None:
                p = _get_pos(res, "CA")
            buried = 0.0
            if p is not None and tree is not None:
                n_nbr = len(tree.query_ball_point(p, 8.0))
                buried = min(100.0, max(0.0, (n_nbr - 15) * 1.6))
            mpka = MODEL_PKA[std]
            # Burial favors neutral state: raises acid pKa, lowers base pKa
            is_acid = std in ("ASP", "GLU", "CYS", "TYR")
            shift = (buried / 100.0) * (1.8 if is_acid else -1.5)
            pka = round(mpka + shift, 2)
            groups.append({
                "chain": ch.name,
                "resseq": res.seqid.num,
                "icode": res.seqid.icode.strip(),
                "resname": std,
                "key": f"{ch.name}:{res.seqid.num}{res.seqid.icode.strip()}",
                "model_pka": mpka,
                "pka": pka,
                "delta_pka": round(pka - mpka, 2),
                "buried_pct": round(buried, 1),
            })
    return groups, "Empirical Microenvironment Engine"


def analyze_propka(
    st: gemmi.Structure,
    ph: float = 7.4,
    respect_metal_coordination: bool = True,
    respect_disulfides: bool = True,
    his_default: str = "auto",
    overrides: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Run full PROPKA analysis, titration curve, pI, and protonation state predictions."""
    overrides = overrides or {}
    raw_groups, engine_name = _run_propka_raw(st)
    metal_map = _detect_metal_coordinating_residues(st)
    disulfide_keys = _detect_disulfide_keys(st)

    enriched = []
    for g in raw_groups:
        key = g["key"]
        base_key = f"{g['chain']}:{g['resseq']}{g['icode']}"
        is_ss = base_key in disulfide_keys and g["resname"] == "CYS"
        m_info = metal_map.get(base_key, [])
        taut_hint = (
            _choose_his_tautomer(st, g["chain"], g["resseq"], g["icode"], m_info)
            if g["resname"] == "HIS"
            else (g["resname"], "")
        )

        if base_key in overrides and g["resname"] not in ("N+", "C-"):
            state = overrides[base_key].upper()
            q_map = {"HIP": 1.0, "HID": 0.0, "HIE": 0.0, "ASP": -1.0, "ASH": 0.0, "GLU": -1.0, "GLH": 0.0, "CYS": 0.0, "CYM": -1.0, "CYX": 0.0, "LYS": 1.0, "LYN": 0.0, "TYR": 0.0, "TYM": -1.0, "ARG": 1.0, "ARN": 0.0}
            charge = q_map.get(state, 0.0)
            reason = "Manual user override"
            is_override = True
        else:
            state, charge, reason = _assign_state(
                g["resname"], g["pka"], ph, his_default, is_ss, m_info,
                respect_metal_coordination, respect_disulfides, taut_hint,
            )
            is_override = False

        # Available allowed states for UI dropdown
        choices_map = {
            "HIS": ["HID", "HIE", "HIP"],
            "ASP": ["ASP", "ASH"],
            "GLU": ["GLU", "GLH"],
            "CYS": ["CYS", "CYM", "CYX"],
            "LYS": ["LYS", "LYN"],
            "TYR": ["TYR", "TYM"],
            "ARG": ["ARG", "ARN"],
        }

        enriched.append({
            **g,
            "base_key": base_key,
            "assigned_state": state,
            "charge_at_ph": charge,
            "rationale": reason,
            "is_disulfide": is_ss,
            "metal_coordination": m_info,
            "is_override": is_override,
            "shifted_significantly": abs(g["delta_pka"]) >= 1.5,
            "available_states": choices_map.get(g["resname"], [state]),
        })

    # Compute titration curve (pH 0.0 to 14.0 in 0.5 increments) and exact pI
    def _net_charge_at(test_ph: float) -> float:
        q_total = 0.0
        for g in enriched:
            if g["is_disulfide"] and respect_disulfides:
                continue
            rname = g["resname"]
            if rname not in TITRATABLE_CHARGES:
                continue
            qp, qd = TITRATABLE_CHARGES[rname]
            fp = 1.0 / (1.0 + 10.0 ** (test_ph - g["pka"]))
            q_total += qp * fp + qd * (1.0 - fp)
        return q_total

    titration_curve = []
    pi_est = 7.0
    prev_ph, prev_q = 0.0, _net_charge_at(0.0)
    for step in range(29):
        p_val = round(step * 0.5, 1)
        q_val = _net_charge_at(p_val)
        titration_curve.append({"ph": p_val, "net_charge": round(q_val, 2)})
        if prev_q >= 0.0 >= q_val and abs(prev_q - q_val) > 1e-6:
            pi_est = round(prev_ph + (p_val - prev_ph) * (prev_q / (prev_q - q_val)), 2)
        prev_ph, prev_q = p_val, q_val

    net_charge_discrete = round(sum(g["charge_at_ph"] for g in enriched), 2)
    net_charge_exact = round(_net_charge_at(ph), 2)

    return {
        "engine": engine_name,
        "ph": ph,
        "isoelectric_point": pi_est,
        "net_charge_discrete": net_charge_discrete,
        "net_charge_continuous": net_charge_exact,
        "titratable_count": len(enriched),
        "shifted_count": sum(1 for g in enriched if g["shifted_significantly"]),
        "metal_coordinated_count": sum(1 for g in enriched if g["metal_coordination"]),
        "residues": enriched,
        "titration_curve": titration_curve,
    }


def _add_hydrogens_to_structure(st: gemmi.Structure, assigned_states: Dict[str, str]) -> int:
    """
    Build accurate polar (amide NH, side-chain OH/SH/NH/NH2/NH3+/His HD1/HE2) and
    C-alpha hydrogens respecting assigned protonation states (HID/HIE/HIP/ASH/GLH/CYM/CYX/LYN).
    """
    if len(st) == 0:
        return 0
    added = 0

    for chain in st[0]:
        prot_res = [r for r in chain if canonical_resname(r.name) in ALL_PROTEIN_RESNAMES]
        # Strip existing hydrogens first to rebuild cleanly without duplicates
        for res in prot_res:
            for ai in range(len(res) - 1, -1, -1):
                if res[ai].element.name in ("H", "D"):
                    del res[ai]

        for idx, res in enumerate(prot_res):
            key = f"{chain.name}:{res.seqid.num}{res.seqid.icode.strip()}"
            state = assigned_states.get(key, res.name.strip().upper())
            std = canonical_resname(state)
            b_iso = float(np.mean([a.b_iso for a in res] or [25.0]))

            def _put_h(hname: str, pos: np.ndarray) -> None:
                nonlocal added
                h = gemmi.Atom()
                h.name = hname
                h.element = gemmi.Element("H")
                h.pos = gemmi.Position(float(pos[0]), float(pos[1]), float(pos[2]))
                h.occ = 1.0
                h.b_iso = round(b_iso, 2)
                res.add_atom(h)
                added += 1

            n_p = _get_pos(res, "N")
            ca_p = _get_pos(res, "CA")
            c_p = _get_pos(res, "C")

            # Backbone amide H (non-Pro, non-N-terminal)
            if idx > 0 and std != "PRO" and n_p is not None and ca_p is not None:
                prev_c = _get_pos(prot_res[idx - 1], "C")
                if prev_c is not None and distance(prev_c, n_p) < 2.2:
                    bis = (prev_c - n_p) / (np.linalg.norm(prev_c - n_p) + 1e-8) + (ca_p - n_p) / (np.linalg.norm(ca_p - n_p) + 1e-8)
                    if np.linalg.norm(bis) > 1e-4:
                        h_pos = n_p - (bis / np.linalg.norm(bis)) * 1.01
                        _put_h("H", h_pos)
            elif idx == 0 and std not in ("PRO", "ACE") and all(p is not None for p in (n_p, ca_p, c_p)):
                # N-terminal NH3+
                for hi, tors in enumerate([60.0, 180.0, -60.0], start=1):
                    _put_h(f"H{hi}", nerf_place(c_p, ca_p, n_p, 1.01, 109.5, tors))

            # C-alpha HA (except GLY)
            cb_p = _get_pos(res, "CB")
            if std != "GLY" and all(p is not None for p in (n_p, ca_p, c_p, cb_p)):
                v = (n_p - ca_p) + (c_p - ca_p) + (cb_p - ca_p)
                if np.linalg.norm(v) > 1e-4:
                    _put_h("HA", ca_p - (v / np.linalg.norm(v)) * 1.09)

            # Side-chain polar hydrogens according to protonation state
            if state in ("HID", "HIP"):
                cg, nd1, ce1 = _get_pos(res, "CG"), _get_pos(res, "ND1"), _get_pos(res, "CE1")
                if all(p is not None for p in (cg, nd1, ce1)):
                    v = (cg - nd1) / (np.linalg.norm(cg - nd1) + 1e-8) + (ce1 - nd1) / (np.linalg.norm(ce1 - nd1) + 1e-8)
                    if np.linalg.norm(v) > 1e-4:
                        _put_h("HD1", nd1 - (v / np.linalg.norm(v)) * 1.01)
            if state in ("HIE", "HIP"):
                cd2, ne2, ce1 = _get_pos(res, "CD2"), _get_pos(res, "NE2"), _get_pos(res, "CE1")
                if all(p is not None for p in (cd2, ne2, ce1)):
                    v = (cd2 - ne2) / (np.linalg.norm(cd2 - ne2) + 1e-8) + (ce1 - ne2) / (np.linalg.norm(ce1 - ne2) + 1e-8)
                    if np.linalg.norm(v) > 1e-4:
                        _put_h("HE2", ne2 - (v / np.linalg.norm(v)) * 1.01)

            if state == "ASH":
                cb, cg, od2 = _get_pos(res, "CB"), _get_pos(res, "CG"), _get_pos(res, "OD2")
                if all(p is not None for p in (cb, cg, od2)):
                    _put_h("HD2", nerf_place(cb, cg, od2, 0.97, 109.5, 180.0))
            elif state == "GLH":
                cg, cd, oe2 = _get_pos(res, "CG"), _get_pos(res, "CD"), _get_pos(res, "OE2")
                if all(p is not None for p in (cg, cd, oe2)):
                    _put_h("HE2", nerf_place(cg, cd, oe2, 0.97, 109.5, 180.0))
            elif state == "CYS":
                ca, cb, sg = _get_pos(res, "CA"), _get_pos(res, "CB"), _get_pos(res, "SG")
                if all(p is not None for p in (ca, cb, sg)):
                    _put_h("HG", nerf_place(ca, cb, sg, 1.33, 96.0, 180.0))
            elif state == "SER":
                ca, cb, og = _get_pos(res, "CA"), _get_pos(res, "CB"), _get_pos(res, "OG")
                if all(p is not None for p in (ca, cb, og)):
                    _put_h("HG", nerf_place(ca, cb, og, 0.96, 109.5, 180.0))
            elif state == "THR":
                ca, cb, og1 = _get_pos(res, "CA"), _get_pos(res, "CB"), _get_pos(res, "OG1")
                if all(p is not None for p in (ca, cb, og1)):
                    _put_h("HG1", nerf_place(ca, cb, og1, 0.96, 109.5, 180.0))
            elif state == "TYR":
                ce1, cz, oh = _get_pos(res, "CE1"), _get_pos(res, "CZ"), _get_pos(res, "OH")
                if all(p is not None for p in (ce1, cz, oh)):
                    _put_h("HH", nerf_place(ce1, cz, oh, 0.96, 109.5, 0.0))
            elif state in ("LYS", "LYN"):
                cd, ce, nz = _get_pos(res, "CD"), _get_pos(res, "CE"), _get_pos(res, "NZ")
                if all(p is not None for p in (cd, ce, nz)):
                    angles = [60.0, 180.0, -60.0] if state == "LYS" else [60.0, -60.0]
                    for hi, tors in enumerate(angles, start=1):
                        _put_h(f"HZ{hi}", nerf_place(cd, ce, nz, 1.01, 109.5, tors))
            elif state == "ASN":
                cb, cg, nd2 = _get_pos(res, "CB"), _get_pos(res, "CG"), _get_pos(res, "ND2")
                if all(p is not None for p in (cb, cg, nd2)):
                    _put_h("HD21", nerf_place(cb, cg, nd2, 1.01, 120.0, 0.0))
                    _put_h("HD22", nerf_place(cb, cg, nd2, 1.01, 120.0, 180.0))
            elif state == "GLN":
                cg, cd, ne2 = _get_pos(res, "CG"), _get_pos(res, "CD"), _get_pos(res, "NE2")
                if all(p is not None for p in (cg, cd, ne2)):
                    _put_h("HE21", nerf_place(cg, cd, ne2, 1.01, 120.0, 0.0))
                    _put_h("HE22", nerf_place(cg, cd, ne2, 1.01, 120.0, 180.0))
            elif state == "TRP":
                cd1, ne1, ce2 = _get_pos(res, "CD1"), _get_pos(res, "NE1"), _get_pos(res, "CE2")
                if all(p is not None for p in (cd1, ne1, ce2)):
                    v = (cd1 - ne1) / (np.linalg.norm(cd1 - ne1) + 1e-8) + (ce2 - ne1) / (np.linalg.norm(ce2 - ne1) + 1e-8)
                    if np.linalg.norm(v) > 1e-4:
                        _put_h("HE1", ne1 - (v / np.linalg.norm(v)) * 1.01)
            elif state in ("ARG", "ARN"):
                cd, ne, cz = _get_pos(res, "CD"), _get_pos(res, "NE"), _get_pos(res, "CZ")
                nh1, nh2 = _get_pos(res, "NH1"), _get_pos(res, "NH2")
                if all(p is not None for p in (cd, ne, cz)):
                    v = (cd - ne) / (np.linalg.norm(cd - ne) + 1e-8) + (cz - ne) / (np.linalg.norm(cz - ne) + 1e-8)
                    if np.linalg.norm(v) > 1e-4:
                        _put_h("HE", ne - (v / np.linalg.norm(v)) * 1.01)
                if all(p is not None for p in (ne, cz, nh1)):
                    _put_h("HH11", nerf_place(ne, cz, nh1, 1.01, 120.0, 0.0))
                    _put_h("HH12", nerf_place(ne, cz, nh1, 1.01, 120.0, 180.0))
                if all(p is not None for p in (ne, cz, nh2)):
                    _put_h("HH21", nerf_place(ne, cz, nh2, 1.01, 120.0, 0.0))
                    if state == "ARG":
                        _put_h("HH22", nerf_place(ne, cz, nh2, 1.01, 120.0, 180.0))

    return added


def apply_protonation(st: gemmi.Structure, opts: Dict[str, Any]) -> Dict[str, Any]:
    """Apply PROPKA-driven protonation states, H-bond NQH flips, and hydrogen placement."""
    actions: List[str] = []

    if opts.get("optimize_hbond_network", True):
        flips = _detect_nqh_flips(st)
        for cand in flips:
            try:
                res = flip_nqh_residue(st, cand["chain"], cand["resseq"], cand["icode"])
                actions.append(res["message"])
            except Exception:
                pass

    analysis = analyze_propka(
        st,
        ph=float(opts.get("ph", 7.4)),
        respect_metal_coordination=bool(opts.get("respect_metal_coordination", True)),
        respect_disulfides=bool(opts.get("respect_disulfides", True)),
        his_default=str(opts.get("his_default", "auto")),
        overrides=opts.get("overrides") or {},
    )

    state_map: Dict[str, str] = {}
    for g in analysis["residues"]:
        if g["resname"] not in ("N+", "C-"):
            state_map[g["base_key"]] = g["assigned_state"]

    renamed_counts: Dict[str, int] = {}
    if len(st) > 0:
        for chain in st[0]:
            for res in chain:
                key = f"{chain.name}:{res.seqid.num}{res.seqid.icode.strip()}"
                if key in state_map:
                    new_state = state_map[key]
                    if res.name.strip().upper() != new_state:
                        res.name = new_state
                        renamed_counts[new_state] = renamed_counts.get(new_state, 0) + 1

    if renamed_counts:
        summary = ", ".join(f"{k}×{v}" for k, v in sorted(renamed_counts.items()))
        actions.append(f"Assigned protonation states at pH {analysis['ph']}: {summary}")

    h_added = 0
    if opts.get("add_hydrogens", True):
        h_added = _add_hydrogens_to_structure(st, state_map)
        actions.append(f"Added {h_added} state-specific hydrogen atom(s) (polar + Cα)")

    return {
        "actions": actions,
        "hydrogens_added": h_added,
        "propka": analysis,
    }
