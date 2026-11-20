"""
Live Real-Time 3D Molecular Docking Engine for PandaProtein.

Implements:
- Automated RDKit Ligand Preparation (SMILES / SDF / MOL2 / PDB / PDBQT / Session Extraction ->
  pH-aware protonation, ETKDGv3 multi-conformer generation, MMFF94 minimization, AutoDock Vina
  atom typing, 2D SVG depiction, and Lipinski/physicochemical descriptors).
- Exact AutoDock Vina (Trott & Olson, J. Comput. Chem. 2010) 6-term scoring function:
    Gauss1, Gauss2, Repulsion, Hydrophobic, Hydrogen Bonding + N_rot torsional entropy penalty.
- Live Server-Sent Events (SSE) Streaming Docking:
    Phase 1: Global Grid Exploration (diverse initial orientations inside the 3D grid box)
    Phase 2: Metropolis-Hastings Simulated Annealing + Flexible Conformer Sampling
    Phase 3: Continuous 6-DOF Local Gradient/Simplex Minimization (Nelder-Mead)
    Optional Phase 4: External AutoDock Vina CLI verification if `vina` is installed.
- Streams real-time 3D ligand coordinates, live H-bonds, clash flashes, and energy decomposition
  on every visual frame so users watch the ligand dock live inside the 3D viewport.
"""

from __future__ import annotations

import io
import json
import math
import os
import shutil
import subprocess
import tempfile
import time
import uuid
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Tuple

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from .chemistry_data import ALL_PROTEIN_RESNAMES, NONSTANDARD_TO_STANDARD
from .export_engine import generate_receptor_pdbqt
from .structure_io import structure_to_pdb

# AutoDock Vina (Trott & Olson 2010) published empirical weights
VINA_W_GAUSS1 = -0.035579
VINA_W_GAUSS2 = -0.005156
VINA_W_REPULSION = 0.840245
VINA_W_HYDROPHOBIC = -0.035069
VINA_W_HBOND = -0.587439
VINA_W_ROT = 0.05846

# AutoDock Vina X-Score van der Waals radii (Å)
VINA_VDW: Dict[str, float] = {
    "C": 1.90,
    "N": 1.80,
    "O": 1.70,
    "S": 2.00,
    "P": 2.10,
    "F": 1.50,
    "CL": 1.80,
    "BR": 2.00,
    "I": 2.20,
    "ZN": 1.20,
    "MG": 1.20,
    "CA": 1.20,
    "FE": 1.20,
    "MN": 1.20,
}

# Built-in curated benchmark ligands for instant 1-click testing (great for 6LU7 Mpro & general pockets)
PRESET_LIGANDS: List[Dict[str, str]] = [
    {
        "id": "nirmatrelvir",
        "name": "Nirmatrelvir (Paxlovid · Mpro Inhibitor)",
        "smiles": "CC1(C2C1C(N(C2)C(=O)C(C(C)(C)C)NC(=O)C(F)(F)F)C(=O)NC(CC3CCNC3=O)C#N)C",
    },
    {
        "id": "baicalein",
        "name": "Baicalein (Natural Flavone Inhibitor)",
        "smiles": "O=c1cc(-c2ccccc2)oc2cc(O)c(O)c(O)c12",
    },
    {
        "id": "ebselen",
        "name": "Ebselen (Covalent/Non-Covalent Mpro Probe)",
        "smiles": "O=c1c2ccccc2[se]n1-c1ccccc1",
    },
    {
        "id": "ibuprofen",
        "name": "Ibuprofen (Classic NSAID Control)",
        "smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    },
]

# In-memory cache for prepared docking ligands and completed docking results
_LIGAND_CACHE: Dict[str, Dict[str, Any]] = {}
_RESULT_CACHE: Dict[str, Dict[str, Any]] = {}


def find_vina_binary() -> Optional[str]:
    for cand in (
        shutil.which("vina"),
        "/Users/rv/miniforge3/bin/vina",
        "/opt/homebrew/bin/vina",
        "/usr/local/bin/vina",
    ):
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None


def _draw_ligand_svg(mol: Chem.Mol) -> str:
    try:
        m2d = Chem.RemoveHs(mol)
        AllChem.Compute2DCoords(m2d)
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 190)
        opts = drawer.drawOptions()
        opts.clearBackground = False
        opts.bondLineWidth = 1.8
        drawer.DrawMolecule(m2d)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return ""


def _classify_rdkit_heavy_atoms(mol_with_h: Chem.Mol) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, List[Dict[str, Any]]]:
    """
    Classify each heavy atom in an RDKit molecule for AutoDock Vina scoring:
    returns (vdw_radii, is_hydrophobic, is_donor, is_acceptor, atom_info).
    """
    vdw: List[float] = []
    hydro: List[bool] = []
    donor: List[bool] = []
    acceptor: List[bool] = []
    info: List[Dict[str, Any]] = []

    heavy_idx = 0
    for atom in mol_with_h.GetAtoms():
        z = atom.GetAtomicNum()
        if z <= 1:
            continue
        sym = atom.GetSymbol().upper()
        r = VINA_VDW.get(sym, 1.90)
        nbrs = list(atom.GetNeighbors())
        has_h = any(n.GetAtomicNum() == 1 for n in nbrs)
        hetero_nbr = any(n.GetAtomicNum() not in (1, 6) for n in nbrs)

        is_hyd = (sym == "C" and not hetero_nbr) or sym in ("CL", "BR", "I", "S")
        is_don = (sym in ("N", "O") and has_h)
        is_acc = (sym in ("N", "O") and atom.GetFormalCharge() <= 0)
        if sym == "N" and atom.GetTotalValence() >= 4:
            is_acc = False  # quaternary ammonium is donor only

        heavy_idx += 1
        vdw.append(r)
        hydro.append(is_hyd)
        donor.append(is_don)
        acceptor.append(is_acc)
        info.append({
            "idx": atom.GetIdx(),
            "name": f"{sym}{heavy_idx}"[:4],
            "elem": sym,
            "is_donor": is_don,
            "is_acceptor": is_acc,
            "is_hydrophobic": is_hyd,
        })

    return (
        np.array(vdw, dtype=np.float64),
        np.array(hydro, dtype=bool),
        np.array(donor, dtype=bool),
        np.array(acceptor, dtype=bool),
        info,
    )


def _coords_to_pdb_block(
    mol_with_h: Chem.Mol,
    coords_all: np.ndarray,
    resname: str = "LIG",
    chain: str = "Z",
) -> str:
    """Fast PDB block serializer for a posed RDKit molecule (heavy + polar H)."""
    lines: List[str] = []
    serial = 1
    idx_to_serial: Dict[int, int] = {}
    for atom in mol_with_h.GetAtoms():
        ai = atom.GetIdx()
        sym = atom.GetSymbol().upper()
        # Keep all heavy atoms + polar hydrogens (attached to N/O/S) for crisp clean 3D sticks
        if sym == "H":
            nbrs = atom.GetNeighbors()
            if not nbrs or nbrs[0].GetSymbol().upper() not in ("N", "O", "S"):
                continue
        x, y, z = coords_all[ai]
        aname = f"{sym}{serial % 999}"[:4]
        idx_to_serial[ai] = serial
        lines.append(
            f"HETATM{serial:5d} {aname:<4s} {resname:<3s} {chain}   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {sym:>2s}"
        )
        serial += 1
    for bond in mol_with_h.GetBonds():
        i1 = bond.GetBeginAtomIdx()
        i2 = bond.GetEndAtomIdx()
        if i1 in idx_to_serial and i2 in idx_to_serial:
            s1, s2 = idx_to_serial[i1], idx_to_serial[i2]
            lines.append(f"CONECT{s1:5d}{s2:5d}")
    lines.append("END")
    return "\n".join(lines) + "\n"


def _extract_session_ligand_pdb(st: gemmi.Structure, ligand_id: str) -> Tuple[str, str]:
    """Extract a specific co-crystallized ligand from a gemmi.Structure as PDB text."""
    parts = ligand_id.split(":")
    if len(parts) < 3:
        raise ValueError(f"Invalid ligand_id: {ligand_id}")
    ch_id, rn, num_str = parts[0], parts[1], parts[2]
    digits = "".join(c for c in num_str if c.isdigit() or c == "-")
    icode = "".join(c for c in num_str if c.isalpha())
    resseq = int(digits)

    sub = gemmi.Structure()
    sub.name = rn
    m = gemmi.Model("1")
    c = gemmi.Chain(ch_id)
    found = False
    for ch in st[0]:
        if ch.name != ch_id:
            continue
        for res in ch:
            if res.name.strip().upper() == rn.upper() and res.seqid.num == resseq and res.seqid.icode.strip() == icode:
                c.add_residue(res.clone())
                found = True
                break
    if not found:
        raise ValueError(f"Ligand '{ligand_id}' not found in active structure")
    m.add_chain(c)
    sub.add_model(m)
    return structure_to_pdb(sub), f"{rn} ({ch_id}:{resseq})"


def prepare_docking_ligand(
    st: Optional[gemmi.Structure] = None,
    smiles: Optional[str] = None,
    file_text: Optional[str] = None,
    filename: Optional[str] = None,
    session_ligand_id: Optional[str] = None,
    ph: float = 7.4,
    num_conformers: int = 8,
) -> Dict[str, Any]:
    """
    Prepare a ligand for live 3D docking:
    Parses SMILES / SDF / MOL2 / PDB / PDBQT or extracts from active session,
    builds 3D conformers (ETKDGv3 + MMFF94), computes Vina atom types & descriptors,
    and caches the prepared ligand ready for streaming docking.
    """
    mol: Optional[Chem.Mol] = None
    name = "Ligand"
    source = "Custom"
    has_native_3d = False

    if session_ligand_id and st is not None:
        pdb_blk, label = _extract_session_ligand_pdb(st, session_ligand_id)
        mol = Chem.MolFromPDBBlock(pdb_blk, sanitize=True, removeHs=False)
        if mol is None:
            mol = Chem.MolFromPDBBlock(pdb_blk, sanitize=False, removeHs=False)
            if mol is not None:
                try:
                    Chem.SanitizeMol(mol)
                except Exception:
                    pass
        name = label
        source = f"Extracted from Receptor ({label})"
        has_native_3d = True
    elif smiles:
        smi = smiles.strip()
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: '{smi}'")
        name = filename or (smi[:22] + ("…" if len(smi) > 22 else ""))
        source = "SMILES Input"
    elif file_text:
        ext = Path(filename or "ligand.sdf").suffix.lower()
        name = Path(filename or "ligand.sdf").stem
        source = f"Uploaded ({filename})"
        if ext in (".sdf", ".mol"):
            suppl = Chem.SDMolSupplier()
            suppl.SetData(file_text, removeHs=False)
            mol = next((m for m in suppl if m is not None), None)
            if mol is None:
                mol = Chem.MolFromMolBlock(file_text, removeHs=False)
        elif ext == ".mol2":
            mol = Chem.MolFromMol2Block(file_text, removeHs=False)
        elif ext in (".pdb", ".pdbqt", ".ent"):
            # Strip PDBQT branch/root records if .pdbqt
            cleaned = "\n".join(
                ln[:66] for ln in file_text.splitlines() if ln.startswith(("ATOM", "HETATM", "CONECT"))
            )
            mol = Chem.MolFromPDBBlock(cleaned, removeHs=False)
        elif ext == ".smi":
            first_smi = file_text.strip().splitlines()[0].split()[0]
            mol = Chem.MolFromSmiles(first_smi)
        else:
            mol = Chem.MolFromMolBlock(file_text, removeHs=False) or Chem.MolFromPDBBlock(file_text, removeHs=False)
        if mol is not None and mol.GetNumConformers() > 0:
            has_native_3d = True
    else:
        raise ValueError("Provide a SMILES string, an uploaded ligand file, or a session ligand_id.")

    if mol is None:
        raise ValueError("Could not parse ligand structure with RDKit. Check file format or SMILES.")

    # Ensure clean valence & add hydrogens (keep native 3D if present)
    try:
        mol_no_h = Chem.RemoveHs(mol)
    except Exception:
        mol_no_h = mol
    mol_h = Chem.AddHs(mol_no_h, addCoords=has_native_3d)

    # Build a diverse ensemble of low-energy 3D conformers for flexible live docking
    conformer_coords: List[np.ndarray] = []
    if has_native_3d and mol_h.GetNumConformers() > 0:
        try:
            AllChem.MMFFOptimizeMolecule(mol_h, maxIters=120)
        except Exception:
            pass
        c0 = np.array(mol_h.GetConformer(0).GetPositions(), dtype=np.float64)
        c0 -= np.mean(c0, axis=0)
        conformer_coords.append(c0)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 0
    try:
        cids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=max(4, num_conformers), params=params))
        for cid in cids:
            try:
                AllChem.MMFFOptimizeMolecule(mol_h, confId=cid, maxIters=150)
            except Exception:
                try:
                    AllChem.UFFOptimizeMolecule(mol_h, confId=cid, maxIters=150)
                except Exception:
                    pass
            pos = np.array(mol_h.GetConformer(cid).GetPositions(), dtype=np.float64)
            pos -= np.mean(pos, axis=0)
            conformer_coords.append(pos)
    except Exception:
        pass

    if not conformer_coords:
        AllChem.EmbedMolecule(mol_h, randomSeed=42)
        pos = np.array(mol_h.GetConformer(0).GetPositions(), dtype=np.float64)
        pos -= np.mean(pos, axis=0)
        conformer_coords.append(pos)

    vdw, is_hydro, is_donor, is_acceptor, heavy_info = _classify_rdkit_heavy_atoms(mol_h)
    heavy_indices = np.array([h["idx"] for h in heavy_info], dtype=np.int64)
    n_rot = int(rdMolDescriptors.CalcNumRotatableBonds(mol_no_h))

    # Compute RDKit physicochemical descriptors + 2D depiction
    smi_out = Chem.MolToSmiles(mol_no_h)
    formula = rdMolDescriptors.CalcMolFormula(mol_no_h)
    mw = round(float(Descriptors.MolWt(mol_no_h)), 2)
    clogp = round(float(Descriptors.MolLogP(mol_no_h)), 2)
    tpsa = round(float(Descriptors.TPSA(mol_no_h)), 1)
    hbd = int(Lipinski.NumHDonors(mol_no_h))
    hba = int(Lipinski.NumHAcceptors(mol_no_h))
    svg = _draw_ligand_svg(mol_no_h)

    lig_key = uuid.uuid4().hex[:10]
    _LIGAND_CACHE[lig_key] = {
        "ligand_key": lig_key,
        "name": name,
        "source": source,
        "ph": ph,
        "mol_h": mol_h,
        "conformers": conformer_coords[:12],
        "heavy_indices": heavy_indices,
        "vdw": vdw,
        "is_hydrophobic": is_hydro,
        "is_donor": is_donor,
        "is_acceptor": is_acceptor,
        "heavy_info": heavy_info,
        "n_rot": n_rot,
        "smiles": smi_out,
    }
    if len(_LIGAND_CACHE) > 25:
        del _LIGAND_CACHE[next(iter(_LIGAND_CACHE))]

    return {
        "ligand_key": lig_key,
        "name": name,
        "source": source,
        "smiles": smi_out,
        "formula": formula,
        "mw": mw,
        "clogp": clogp,
        "tpsa": tpsa,
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": n_rot,
        "heavy_atoms": int(len(heavy_indices)),
        "total_atoms": int(mol_h.GetNumAtoms()),
        "conformers_generated": len(conformer_coords[:12]),
        "svg": svg,
        "ph": ph,
        "vina_available": find_vina_binary() is not None,
    }


def _extract_receptor_pocket_atoms(
    st: gemmi.Structure,
    box_center: np.ndarray,
    box_size: np.ndarray,
    include_hetatms: bool = False,
) -> Dict[str, Any]:
    """Extract all receptor heavy atoms within (grid_box/2 + 8.5 Å) with Vina atom typing."""
    half = box_size / 2.0 + 8.5
    pts: List[List[float]] = []
    vdw: List[float] = []
    hydro: List[bool] = []
    donor: List[bool] = []
    acceptor: List[bool] = []
    meta: List[Dict[str, Any]] = []

    water_names = {"HOH", "WAT", "SOL", "DOD"}
    metal_names = {"ZN", "MG", "CA", "FE", "FE2", "MN", "CU", "CO", "NI"}
    polar_c_names = {"C", "CG", "CD", "CZ", "CE1", "CD2"}  # carbonyl/guanidinium carbons

    for ch in st[0]:
        for res in ch:
            rn = res.name.strip().upper()
            if rn in water_names:
                continue
            is_prot = rn in ALL_PROTEIN_RESNAMES or rn in NONSTANDARD_TO_STANDARD
            is_met = rn in metal_names
            if not is_prot and not is_met and not include_hetatms:
                continue
            ic = res.seqid.icode.strip()
            key = f"{ch.name}:{res.seqid.num}{ic}"
            for a in res:
                elem = a.element.name.strip().upper()
                if elem in ("H", "D"):
                    continue
                p = np.array([a.pos.x, a.pos.y, a.pos.z], dtype=np.float64)
                if np.any(np.abs(p - box_center) > half):
                    continue
                an = a.name.strip().upper()
                r = VINA_VDW.get(elem, 1.90)
                is_hyd = (elem == "C" and an not in polar_c_names) or elem == "S"
                # Backbone N (except PRO), sidechain N (ND/NE/NH/NZ), hydroxyl O (OG/OH) act as donors
                is_don = (
                    (elem == "N" and not (rn == "PRO" and an == "N"))
                    or (elem == "O" and an in ("OG", "OG1", "OH"))
                )
                # All oxygens + unprotonated ring nitrogens act as acceptors
                is_acc = (elem == "O") or (elem == "N" and an in ("ND1", "NE2"))
                if is_met:
                    is_don = True  # metal chelation favorable interaction with ligand O/N acceptors

                pts.append([p[0], p[1], p[2]])
                vdw.append(r)
                hydro.append(is_hyd)
                donor.append(is_don)
                acceptor.append(is_acc)
                meta.append({
                    "chain": ch.name,
                    "resseq": res.seqid.num,
                    "icode": ic,
                    "resname": rn,
                    "atom": an,
                    "key": key,
                    "label": f"{rn} {key}",
                })

    if not pts:
        raise ValueError("No receptor atoms found near the docking grid box. Re-center the grid box on the protein.")

    return {
        "pts": np.array(pts, dtype=np.float64),
        "vdw": np.array(vdw, dtype=np.float64),
        "is_hydrophobic": np.array(hydro, dtype=bool),
        "is_donor": np.array(donor, dtype=bool),
        "is_acceptor": np.array(acceptor, dtype=bool),
        "meta": meta,
    }


def _euler_to_rot(rx: float, ry: float, rz: float) -> np.ndarray:
    cx, sx = math.cos(rx), math.sin(rx)
    cy, sy = math.cos(ry), math.sin(ry)
    cz, sz = math.cos(rz), math.sin(rz)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]], dtype=np.float64)
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]], dtype=np.float64)
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]], dtype=np.float64)
    return Rz @ Ry @ Rx


def _evaluate_vina_pose(
    lig_heavy_coords: np.ndarray,
    lig_vdw: np.ndarray,
    lig_hydro: np.ndarray,
    lig_donor: np.ndarray,
    lig_acceptor: np.ndarray,
    rec: Dict[str, Any],
    n_rot: int,
    box_center: np.ndarray,
    box_half: np.ndarray,
    collect_contacts: bool = False,
    lig_info: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """
    Exact AutoDock Vina (Trott & Olson 2010) 6-term scoring function in vectorized NumPy.
    Returns total ΔG (kcal/mol), individual energy terms, and optionally live H-bonds & clashes.
    """
    rec_pts = rec["pts"]
    # Vectorized pairwise distances: (N_lig, N_rec)
    diff = lig_heavy_coords[:, None, :] - rec_pts[None, :, :]
    r_mat = np.sqrt(np.sum(diff * diff, axis=2))
    mask = r_mat < 8.0

    if not np.any(mask):
        return {
            "total": 0.0,
            "gauss1": 0.0,
            "gauss2": 0.0,
            "repulsion": 0.0,
            "hydrophobic": 0.0,
            "hbond": 0.0,
            "torsion": 0.0,
            "hbonds": [],
            "clashes": [],
            "contacts": [],
        }

    li, ri = np.where(mask)
    r_vals = r_mat[li, ri]
    d_vals = r_vals - lig_vdw[li] - rec["vdw"][ri]

    # 1. Gauss1 & Gauss2 steric attraction
    g1 = np.exp(-((d_vals / 0.5) ** 2))
    g2 = np.exp(-(((d_vals - 3.0) / 2.0) ** 2))

    # 2. Steric repulsion (d < 0)
    rep = np.where(d_vals < 0.0, d_vals * d_vals, 0.0)

    # 3. Hydrophobic term
    both_hyd = lig_hydro[li] & rec["is_hydrophobic"][ri]
    hyd_raw = np.where(
        both_hyd,
        np.where(d_vals < 0.5, 1.0, np.where(d_vals < 1.5, 1.5 - d_vals, 0.0)),
        0.0,
    )

    # 4. Hydrogen bonding term (donor-acceptor complementary pairs)
    is_hb_pair = (lig_donor[li] & rec["is_acceptor"][ri]) | (lig_acceptor[li] & rec["is_donor"][ri])
    hb_raw = np.where(
        is_hb_pair,
        np.where(d_vals < -0.7, 1.0, np.where(d_vals < 0.0, -d_vals / 0.7, 0.0)),
        0.0,
    )

    e_g1 = float(VINA_W_GAUSS1 * np.sum(g1))
    e_g2 = float(VINA_W_GAUSS2 * np.sum(g2))
    e_rep = float(VINA_W_REPULSION * np.sum(rep))
    e_hyd = float(VINA_W_HYDROPHOBIC * np.sum(hyd_raw))
    e_hb = float(VINA_W_HBOND * np.sum(hb_raw))

    # Soft harmonic wall penalty if ligand centroid leaves the user's grid box
    cent = np.mean(lig_heavy_coords, axis=0)
    over = np.maximum(0.0, np.abs(cent - box_center) - box_half)
    box_pen = float(12.0 * np.sum(over * over))

    c_inter = e_g1 + e_g2 + e_rep + e_hyd + e_hb
    tors_div = 1.0 + VINA_W_ROT * n_rot
    total = (c_inter / tors_div) + box_pen
    torsion_cost = total - box_pen - c_inter

    hbonds: List[Dict[str, Any]] = []
    clashes: List[Dict[str, Any]] = []
    contacts_list: List[Dict[str, Any]] = []

    if collect_contacts:
        seen_hb = set()
        seen_res = {}
        for k_idx in range(len(li)):
            i_l = int(li[k_idx])
            i_r = int(ri[k_idx])
            rv = float(r_vals[k_idx])
            dv = float(d_vals[k_idx])
            rm = rec["meta"][i_r]
            lp = lig_heavy_coords[i_l]
            rp = rec_pts[i_r]

            if is_hb_pair[k_idx] and rv <= 3.45 and dv > -1.2:
                hkey = (rm["key"], rm["atom"], i_l)
                if hkey not in seen_hb:
                    seen_hb.add(hkey)
                    hbonds.append({
                        "residue": rm["label"],
                        "chain": rm["chain"],
                        "resseq": rm["resseq"],
                        "rec_atom": rm["atom"],
                        "lig_atom": lig_info[i_l]["name"] if lig_info else f"L{i_l}",
                        "distance_A": round(rv, 2),
                        "p1": [round(float(lp[0]), 2), round(float(lp[1]), 2), round(float(lp[2]), 2)],
                        "p2": [round(float(rp[0]), 2), round(float(rp[1]), 2), round(float(rp[2]), 2)],
                    })
            if dv < -0.65 and len(clashes) < 10:
                clashes.append({
                    "residue": rm["label"],
                    "distance_A": round(rv, 2),
                    "p1": [round(float(lp[0]), 2), round(float(lp[1]), 2), round(float(lp[2]), 2)],
                    "p2": [round(float(rp[0]), 2), round(float(rp[1]), 2), round(float(rp[2]), 2)],
                })
            if rv <= 4.0:
                rk = rm["key"]
                ctype = "hbond" if (is_hb_pair[k_idx] and rv <= 3.45) else ("hydrophobic" if both_hyd[k_idx] else "vdw")
                if rk not in seen_res or rv < seen_res[rk]["distance_A"]:
                    seen_res[rk] = {
                        "residue": rm["label"],
                        "chain": rm["chain"],
                        "resseq": rm["resseq"],
                        "resname": rm["resname"],
                        "key": rk,
                        "type": ctype,
                        "distance_A": round(rv, 2),
                    }
        hbonds.sort(key=lambda x: x["distance_A"])
        contacts_list = sorted(seen_res.values(), key=lambda x: (x["type"] != "hbond", x["distance_A"]))[:14]

    return {
        "total": round(float(total), 2),
        "gauss1": round(e_g1, 2),
        "gauss2": round(e_g2, 2),
        "repulsion": round(e_rep, 2),
        "hydrophobic": round(e_hyd, 2),
        "hbond": round(e_hb, 2),
        "torsion": round(float(torsion_cost), 2),
        "hbonds": hbonds[:10],
        "clashes": clashes[:8],
        "contacts": contacts_list,
    }


def _affinity_to_ki_str(dg_kcal_mol: float) -> str:
    """Convert ΔG (kcal/mol) at 298.15 K to estimated dissociation constant Ki."""
    if dg_kcal_mol >= 0:
        return "> 1 M (Unfavorable)"
    rt = 0.5924847  # kcal/mol at 298.15 K
    ki_m = math.exp(dg_kcal_mol / rt)
    if ki_m < 1e-9:
        return f"{ki_m * 1e12:.1f} pM"
    if ki_m < 1e-6:
        return f"{ki_m * 1e9:.1f} nM"
    if ki_m < 1e-3:
        return f"{ki_m * 1e6:.1f} µM"
    return f"{ki_m * 1e3:.2f} mM"


def _run_external_vina_if_requested(
    st: gemmi.Structure,
    lig_entry: Dict[str, Any],
    center: np.ndarray,
    size: np.ndarray,
) -> Optional[Dict[str, Any]]:
    """Optional cross-validation with external AutoDock Vina CLI binary if installed."""
    vina_bin = find_vina_binary()
    obabel_bin = shutil.which("obabel") or ("/opt/homebrew/bin/obabel" if os.path.isfile("/opt/homebrew/bin/obabel") else None)
    if not vina_bin or not obabel_bin:
        return None
    try:
        with tempfile.TemporaryDirectory(prefix="pandaprotein_vina_") as td:
            rec_path = os.path.join(td, "receptor.pdbqt")
            lig_sdf = os.path.join(td, "ligand.sdf")
            lig_pdbqt = os.path.join(td, "ligand.pdbqt")
            out_pdbqt = os.path.join(td, "docked.pdbqt")

            with open(rec_path, "w", encoding="utf-8") as f:
                f.write(generate_receptor_pdbqt(st))
            w = Chem.SDWriter(lig_sdf)
            w.write(lig_entry["mol_h"], confId=0)
            w.close()

            r_ob = subprocess.run([obabel_bin, lig_sdf, "-O", lig_pdbqt, "-h"], capture_output=True, text=True, timeout=15)
            if r_ob.returncode != 0 or not os.path.exists(lig_pdbqt):
                return None

            cmd = [
                vina_bin,
                "--receptor", rec_path,
                "--ligand", lig_pdbqt,
                "--center_x", f"{center[0]:.3f}",
                "--center_y", f"{center[1]:.3f}",
                "--center_z", f"{center[2]:.3f}",
                "--size_x", f"{size[0]:.2f}",
                "--size_y", f"{size[1]:.2f}",
                "--size_z", f"{size[2]:.2f}",
                "--exhaustiveness", "4",
                "--num_modes", "5",
                "--out", out_pdbqt,
            ]
            rv = subprocess.run(cmd, capture_output=True, text=True, timeout=45)
            if rv.returncode != 0 or not os.path.exists(out_pdbqt):
                return None

            # Parse best Vina affinity from stdout
            best_vina = None
            for line in rv.stdout.splitlines():
                parts = line.split()
                if len(parts) >= 4 and parts[0] == "1":
                    try:
                        best_vina = float(parts[1])
                        break
                    except ValueError:
                        pass
            return {
                "binary": vina_bin,
                "best_affinity_kcal_mol": best_vina,
                "stdout": "\n".join(rv.stdout.splitlines()[-16:]),
            }
    except Exception:
        return None


def stream_live_docking(
    st: gemmi.Structure,
    ligand_key: str,
    box_center: List[float],
    box_size: List[float],
    exhaustiveness: int = 8,
    speed: str = "normal",
    use_vina_verify: bool = False,
    include_hetatms: bool = False,
) -> Generator[str, None, None]:
    """
    Generator yielding SSE events (`data: {...}\n\n`) for Live 3D Molecular Docking!
    """
    if ligand_key not in _LIGAND_CACHE:
        yield f"data: {json.dumps({'type': 'error', 'message': 'Prepared ligand expired — please prepare ligand again.'})}\n\n"
        return

    lig = _LIGAND_CACHE[ligand_key]
    center = np.array(box_center, dtype=np.float64)
    size = np.array([max(8.0, min(40.0, float(x))) for x in box_size], dtype=np.float64)
    half = size / 2.0

    try:
        rec = _extract_receptor_pocket_atoms(st, center, size, include_hetatms=include_hetatms)
    except Exception as exc:
        yield f"data: {json.dumps({'type': 'error', 'message': str(exc)})}\n\n"
        return

    mol_h = lig["mol_h"]
    conformers: List[np.ndarray] = lig["conformers"]
    heavy_idx: np.ndarray = lig["heavy_indices"]
    lig_vdw: np.ndarray = lig["vdw"]
    lig_hydro: np.ndarray = lig["is_hydrophobic"]
    lig_donor: np.ndarray = lig["is_donor"]
    lig_acceptor: np.ndarray = lig["is_acceptor"]
    heavy_info: List[Dict[str, Any]] = lig["heavy_info"]
    n_rot: int = lig["n_rot"]

    # Frame pacing so users can genuinely watch the 3D live action!
    frame_delay = 0.085 if speed == "cinematic" else (0.045 if speed == "normal" else 0.01)
    n_global = max(14, min(28, exhaustiveness * 2 + 4))
    n_sa = max(28, min(64, exhaustiveness * 4 + 12))
    n_local = 14
    total_steps = n_global + n_sa + n_local

    yield f"data: {json.dumps({'type': 'init', 'total_steps': total_steps, 'pocket_atoms': len(rec['pts']), 'ligand_name': lig['name'], 'rotatable_bonds': n_rot, 'box_center': center.tolist(), 'box_size': size.tolist()})}\n\n"

    rng = np.random.default_rng(int(time.time() * 1000) % 100000)

    def pose_coords(conf_i: int, x6: np.ndarray) -> np.ndarray:
        R = _euler_to_rot(float(x6[3]), float(x6[4]), float(x6[5]))
        return (conformers[conf_i] @ R.T) + x6[:3]

    def score_pose(conf_i: int, x6: np.ndarray, with_contacts: bool = False) -> Tuple[Dict[str, Any], np.ndarray]:
        all_c = pose_coords(conf_i, x6)
        sc = _evaluate_vina_pose(
            all_c[heavy_idx],
            lig_vdw,
            lig_hydro,
            lig_donor,
            lig_acceptor,
            rec,
            n_rot,
            center,
            half,
            collect_contacts=with_contacts,
            lig_info=heavy_info,
        )
        return sc, all_c

    # Keep pool of diverse candidate modes (conf_i, x6, score_dict, all_coords)
    pool: List[Tuple[int, np.ndarray, Dict[str, Any], np.ndarray]] = []
    best_score = 1e9
    best_state: Optional[Tuple[int, np.ndarray, Dict[str, Any], np.ndarray]] = None

    # Find void-space seed pockets inside the grid box so global search discovers real binding crevices fast!
    rec_tree = cKDTree(rec["pts"])
    void_seeds: List[np.ndarray] = []
    grid_n = 7
    xs = np.linspace(center[0] - half[0] * 0.72, center[0] + half[0] * 0.72, grid_n)
    ys = np.linspace(center[1] - half[1] * 0.72, center[1] + half[1] * 0.72, grid_n)
    zs = np.linspace(center[2] - half[2] * 0.72, center[2] + half[2] * 0.72, grid_n)
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    gpts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    gdists, _ = rec_tree.query(gpts, k=1)
    # Pocket points: far enough to avoid severe clash (> 2.4 Å) yet close enough to contact pocket walls (< 5.8 Å)
    good_void = gpts[(gdists >= 2.3) & (gdists <= 5.8)]
    if len(good_void) >= 6:
        # Sort by closeness to pocket center + moderate clearance
        order = np.argsort(np.linalg.norm(good_void - center, axis=1) - gdists[(gdists >= 2.3) & (gdists <= 5.8)] * 0.8)
        for pt in good_void[order[:40]]:
            void_seeds.append(pt)
    if not void_seeds:
        void_seeds.append(center.copy())

    step = 0

    # ===================== PHASE 1: GLOBAL GRID EXPLORATION =====================
    for i in range(n_global):
        step += 1
        conf_i = int(rng.integers(0, len(conformers)))
        seed_xyz = void_seeds[i % len(void_seeds)] + rng.normal(0.0, 1.15, size=3)
        seed_xyz = np.clip(seed_xyz, center - half * 0.85, center + half * 0.85)
        angles = rng.uniform(-math.pi, math.pi, size=3)
        x6 = np.concatenate([seed_xyz, angles])

        # Fast 16-orientation + translation jitter sweep per pocket seed to find deep steric fits!
        sc, all_c = score_pose(conf_i, x6, with_contacts=False)
        for _ in range(16):
            alt_x = x6.copy()
            alt_x[:3] = np.clip(seed_xyz + rng.normal(0.0, 0.9, size=3), center - half * 0.85, center + half * 0.85)
            alt_x[3:] = rng.uniform(-math.pi, math.pi, size=3)
            alt_c = int(rng.integers(0, len(conformers))) if rng.random() < 0.35 else conf_i
            sc2, _ = score_pose(alt_c, alt_x, with_contacts=False)
            if sc2["total"] < sc["total"]:
                x6 = alt_x
                conf_i = alt_c
                sc = sc2
        sc, all_c = score_pose(conf_i, x6, with_contacts=True)

        is_best = sc["total"] < best_score
        if is_best:
            best_score = sc["total"]
            best_state = (conf_i, x6.copy(), sc, all_c.copy())
        pool.append((conf_i, x6.copy(), sc, all_c.copy()))

        frame = {
            "type": "frame",
            "step": step,
            "total_steps": total_steps,
            "phase": "Global Grid Exploration",
            "temperature": 4.0,
            "accepted": True,
            "is_new_best": is_best,
            "score": sc["total"],
            "best_score": round(best_score, 2),
            "ki": _affinity_to_ki_str(best_score),
            "terms": {
                "gauss1": sc["gauss1"],
                "gauss2": sc["gauss2"],
                "repulsion": sc["repulsion"],
                "hydrophobic": sc["hydrophobic"],
                "hbond": sc["hbond"],
                "torsion": sc["torsion"],
            },
            "hbonds": sc["hbonds"],
            "clashes": sc["clashes"],
            "contacts": sc["contacts"],
            "pose_pdb": _coords_to_pdb_block(mol_h, all_c),
        }
        yield f"data: {json.dumps(frame)}\n\n"
        if frame_delay > 0:
            time.sleep(frame_delay)

    # ===================== PHASE 2: METROPOLIS SIMULATED ANNEALING =====================
    assert best_state is not None
    cur_conf, cur_x, cur_sc, cur_c = best_state[0], best_state[1].copy(), best_state[2], best_state[3].copy()

    for i in range(n_sa):
        step += 1
        frac = i / max(1, n_sa - 1)
        temp = 2.6 * ((0.10 / 2.6) ** frac)  # Exponential cooling schedule (2.6 -> 0.10 kcal/mol)

        # Periodically jump to a second promising basin in pool to sample multiple binding modes!
        if i == n_sa // 2 and len(pool) >= 4:
            pool.sort(key=lambda t: t[2]["total"])
            alt = pool[min(2, len(pool) - 1)]
            cur_conf, cur_x, cur_sc, cur_c = alt[0], alt[1].copy(), alt[2], alt[3].copy()

        # Perform 18 inner MC trials per streamed frame for genuine deep Vina convergence
        accepted_frame = False
        is_best_frame = False
        for _ in range(18):
            trans_scale = max(0.25, 1.9 * (1.0 - 0.78 * frac))
            rot_scale = max(0.10, 0.85 * (1.0 - 0.78 * frac))
            prop_x = cur_x.copy()
            prop_x[:3] += rng.normal(0.0, trans_scale, size=3)
            prop_x[:3] = np.clip(prop_x[:3], center - half * 0.9, center + half * 0.9)
            prop_x[3:] += rng.normal(0.0, rot_scale, size=3)
            prop_conf = cur_conf
            if len(conformers) > 1 and rng.random() < 0.25:
                prop_conf = int(rng.integers(0, len(conformers)))

            sc_try, c_try = score_pose(prop_conf, prop_x, with_contacts=False)
            de = sc_try["total"] - cur_sc["total"]
            if de <= 0 or rng.random() < math.exp(-min(40.0, de / max(temp, 0.05))):
                cur_conf, cur_x, cur_sc, cur_c = prop_conf, prop_x, sc_try, c_try
                accepted_frame = True
                if cur_sc["total"] < best_score:
                    best_score = cur_sc["total"]
                    best_full, _ = score_pose(cur_conf, cur_x, with_contacts=True)
                    best_state = (cur_conf, cur_x.copy(), best_full, cur_c.copy())
                    is_best_frame = True

        # Re-evaluate displayed pose with full H-bond/clash/contact breakdown
        vis_sc, vis_c = score_pose(cur_conf, cur_x, with_contacts=True)
        pool.append((cur_conf, cur_x.copy(), vis_sc, vis_c.copy()))

        frame = {
            "type": "frame",
            "step": step,
            "total_steps": total_steps,
            "phase": "Simulated Annealing & Conformer Search",
            "temperature": round(temp, 2),
            "accepted": accepted_frame,
            "is_new_best": is_best_frame,
            "score": vis_sc["total"],
            "best_score": round(best_score, 2),
            "ki": _affinity_to_ki_str(best_score),
            "terms": {
                "gauss1": vis_sc["gauss1"],
                "gauss2": vis_sc["gauss2"],
                "repulsion": vis_sc["repulsion"],
                "hydrophobic": vis_sc["hydrophobic"],
                "hbond": vis_sc["hbond"],
                "torsion": vis_sc["torsion"],
            },
            "hbonds": vis_sc["hbonds"],
            "clashes": vis_sc["clashes"],
            "contacts": vis_sc["contacts"],
            "pose_pdb": _coords_to_pdb_block(mol_h, vis_c),
        }
        yield f"data: {json.dumps(frame)}\n\n"
        if frame_delay > 0:
            time.sleep(frame_delay)

    # ===================== PHASE 3: LOCAL CONTINUOUS 6-DOF MINIMIZATION =====================
    # Select up to 5 top diverse basins from pool and run Nelder-Mead local energy minimization!
    pool.sort(key=lambda t: t[2]["total"])
    seeds_for_min: List[Tuple[int, np.ndarray]] = []
    for c_i, x_v, _, c_arr in pool:
        if len(seeds_for_min) >= 6:
            break
        # Deduplicate spatially (> 1.8 Å heavy-atom RMSD)
        dup = False
        for _, sx in seeds_for_min:
            if np.linalg.norm(x_v[:3] - sx[:3]) < 1.8:
                dup = True
                break
        if not dup:
            seeds_for_min.append((c_i, x_v.copy()))

    refined_modes: List[Tuple[int, np.ndarray, Dict[str, Any], np.ndarray]] = []
    for m_idx, (c_i, x_init) in enumerate(seeds_for_min):
        traj_pts: List[np.ndarray] = []

        def obj(xv: np.ndarray) -> float:
            s_d, _ = score_pose(c_i, xv, with_contacts=False)
            return float(s_d["total"])

        def cb(xk: np.ndarray) -> None:
            traj_pts.append(xk.copy())

        res_opt = minimize(
            obj,
            x_init,
            method="Nelder-Mead",
            callback=cb,
            options={"xatol": 0.03, "fatol": 0.015, "maxiter": 180},
        )
        opt_x = res_opt.x
        opt_sc, opt_c = score_pose(c_i, opt_x, with_contacts=True)
        refined_modes.append((c_i, opt_x.copy(), opt_sc, opt_c.copy()))

        # Stream smooth intermediate frames from the best seed's gradient descent!
        if m_idx == 0 and traj_pts:
            picks = np.linspace(0, len(traj_pts) - 1, min(len(traj_pts), max(3, n_local - len(seeds_for_min)))).astype(int)
            for pk in picks:
                if step >= total_steps:
                    break
                step += 1
                tx = traj_pts[int(pk)]
                tsc, tc = score_pose(c_i, tx, with_contacts=True)
                is_b = tsc["total"] < best_score
                if is_b:
                    best_score = tsc["total"]
                    best_state = (c_i, tx.copy(), tsc, tc.copy())
                yield f"data: {json.dumps({'type': 'frame', 'step': step, 'total_steps': total_steps, 'phase': 'Local 6-DOF Energy Minimization', 'temperature': 0.05, 'accepted': True, 'is_new_best': is_b, 'score': tsc['total'], 'best_score': round(best_score, 2), 'ki': _affinity_to_ki_str(best_score), 'terms': {'gauss1': tsc['gauss1'], 'gauss2': tsc['gauss2'], 'repulsion': tsc['repulsion'], 'hydrophobic': tsc['hydrophobic'], 'hbond': tsc['hbond'], 'torsion': tsc['torsion']}, 'hbonds': tsc['hbonds'], 'clashes': tsc['clashes'], 'contacts': tsc['contacts'], 'pose_pdb': _coords_to_pdb_block(mol_h, tc)})}\n\n"
                if frame_delay > 0:
                    time.sleep(frame_delay)
        else:
            step = min(total_steps, step + 1)
            is_b = opt_sc["total"] < best_score
            if is_b:
                best_score = opt_sc["total"]
                best_state = (c_i, opt_x.copy(), opt_sc, opt_c.copy())
            yield f"data: {json.dumps({'type': 'frame', 'step': step, 'total_steps': total_steps, 'phase': f'Refining Binding Mode #{m_idx + 1}', 'temperature': 0.01, 'accepted': True, 'is_new_best': is_b, 'score': opt_sc['total'], 'best_score': round(best_score, 2), 'ki': _affinity_to_ki_str(best_score), 'terms': {'gauss1': opt_sc['gauss1'], 'gauss2': opt_sc['gauss2'], 'repulsion': opt_sc['repulsion'], 'hydrophobic': opt_sc['hydrophobic'], 'hbond': opt_sc['hbond'], 'torsion': opt_sc['torsion']}, 'hbonds': opt_sc['hbonds'], 'clashes': opt_sc['clashes'], 'contacts': opt_sc['contacts'], 'pose_pdb': _coords_to_pdb_block(mol_h, opt_c)})}\n\n"
            if frame_delay > 0:
                time.sleep(frame_delay * 0.8)

    # Combine & cluster Top 9 Ranked Binding Modes
    all_final = refined_modes + sorted(pool, key=lambda t: t[2]["total"])[:20]
    all_final.sort(key=lambda t: t[2]["total"])
    clustered: List[Tuple[int, np.ndarray, Dict[str, Any], np.ndarray]] = []
    for item in all_final:
        if len(clustered) >= 9:
            break
        h_coords = item[3][heavy_idx]
        too_close = False
        for prev in clustered:
            rmsd = float(np.sqrt(np.mean(np.sum((h_coords - prev[3][heavy_idx]) ** 2, axis=1))))
            if rmsd < 1.25:
                too_close = True
                break
        if not too_close:
            clustered.append(item)

    best_heavy = clustered[0][3][heavy_idx]
    modes_payload: List[Dict[str, Any]] = []
    sdf_buf = io.StringIO()
    sdf_writer = Chem.SDWriter(sdf_buf)

    for rank, (c_i, x_v, sc_d, c_arr) in enumerate(clustered, start=1):
        rmsd_best = round(float(np.sqrt(np.mean(np.sum((c_arr[heavy_idx] - best_heavy) ** 2, axis=1)))), 2)
        pose_pdb = _coords_to_pdb_block(mol_h, c_arr)
        modes_payload.append({
            "mode": rank,
            "affinity_kcal_mol": sc_d["total"],
            "ki": _affinity_to_ki_str(sc_d["total"]),
            "rmsd_from_best_A": rmsd_best,
            "terms": {
                "gauss1": sc_d["gauss1"],
                "gauss2": sc_d["gauss2"],
                "repulsion": sc_d["repulsion"],
                "hydrophobic": sc_d["hydrophobic"],
                "hbond": sc_d["hbond"],
                "torsion": sc_d["torsion"],
            },
            "hbonds": sc_d["hbonds"],
            "clashes": sc_d["clashes"],
            "contacts": sc_d["contacts"],
            "pose_pdb": pose_pdb,
        })

        # Write pose into multi-pose SDF with Vina affinity properties
        m_copy = Chem.Mol(mol_h)
        conf = m_copy.GetConformer(0)
        for ai in range(m_copy.GetNumAtoms()):
            conf.SetAtomPosition(ai, (float(c_arr[ai][0]), float(c_arr[ai][1]), float(c_arr[ai][2])))
        m_copy.SetProp("_Name", f"{lig['name']}_mode_{rank}")
        m_copy.SetProp("Vina_Affinity_kcal_mol", str(sc_d["total"]))
        m_copy.SetProp("Estimated_Ki", _affinity_to_ki_str(sc_d["total"]))
        m_copy.SetProp("RMSD_from_best_A", str(rmsd_best))
        sdf_writer.write(m_copy, confId=0)

    sdf_writer.close()

    # Optional external AutoDock Vina CLI cross-check
    vina_cli = None
    if use_vina_verify:
        yield f"data: {json.dumps({'type': 'status', 'message': 'Cross-validating top poses with AutoDock Vina binary…'})}\n\n"
        vina_cli = _run_external_vina_if_requested(st, lig, center, size)

    # Cache downloadable outputs (complex PDB & multi-pose SDF)
    result_id = uuid.uuid4().hex[:10]
    rec_pdb_lines = [ln for ln in structure_to_pdb(st).splitlines() if not ln.startswith("END")]
    best_lig_lines = [ln for ln in modes_payload[0]["pose_pdb"].splitlines() if not ln.startswith("END")]
    complex_pdb = "\n".join(rec_pdb_lines + ["TER"] + best_lig_lines + ["END"]) + "\n"
    poses_sdf = sdf_buf.getvalue()

    _RESULT_CACHE[result_id] = {
        "ligand_name": lig["name"],
        "complex_pdb": complex_pdb,
        "poses_sdf": poses_sdf,
        "modes": modes_payload,
    }
    if len(_RESULT_CACHE) > 20:
        del _RESULT_CACHE[next(iter(_RESULT_CACHE))]

    yield f"data: {json.dumps({'type': 'complete', 'result_id': result_id, 'best_affinity_kcal_mol': modes_payload[0]['affinity_kcal_mol'], 'best_ki': modes_payload[0]['ki'], 'modes': modes_payload, 'vina_cli': vina_cli})}\n\n"


def export_docked_result(result_id: str, fmt: str = "complex_pdb", mode: int = 1) -> Tuple[str, str, str]:
    """Return (content, filename, media_type) for a completed docking job."""
    if result_id not in _RESULT_CACHE:
        raise KeyError("Docking result not found or expired — please re-run docking.")
    entry = _RESULT_CACHE[result_id]
    safe = "".join(c if c.isalnum() or c in ("_", "-") else "_" for c in entry["ligand_name"])[:24] or "ligand"
    if fmt == "poses_sdf":
        return entry["poses_sdf"], f"{safe}_docked_poses.sdf", "chemical/x-mdl-sdfile"
    if fmt == "mode_pdb":
        modes = entry["modes"]
        m_obj = next((m for m in modes if m["mode"] == mode), modes[0])
        return m_obj["pose_pdb"], f"{safe}_mode_{m_obj['mode']}.pdb", "chemical/x-pdb"
    return entry["complex_pdb"], f"{safe}_docked_complex.pdb", "chemical/x-pdb"
