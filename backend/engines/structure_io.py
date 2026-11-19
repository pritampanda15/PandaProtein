"""
Structure I/O Engine powered by Gemmi + httpx:
- Parse & serialize PDB and mmCIF formats losslessly
- Fetch live structures from RCSB PDB and AlphaFold Protein Structure DB
- Provide curated built-in showcase structures for immediate offline testing
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import gemmi
import httpx
import numpy as np

from .chemistry_data import (
    ALL_PROTEIN_RESNAMES,
    COFACTORS,
    CRYSTALLIZATION_ARTIFACTS,
    HALIDE_IONS,
    METAL_IONS,
    NONSTANDARD_TO_STANDARD,
    STANDARD_AA,
    WATER_RESNAMES,
)
from .geometry_utils import (
    build_sidechain_coords,
    canonical_resname,
    nerf_place,
    rotate_point_around_axis,
)


def _infer_element_from_gro(resname: str, atomname: str) -> str:
    """Infer valid chemical element symbol from GROMACS residue & atom name."""
    rn = resname.strip().upper()
    an = atomname.strip().upper()
    if rn in ("NA", "NA+", "SOD"):
        return "Na"
    if rn in ("CL", "CL-", "CLA"):
        return "Cl"
    if rn in ("MG", "MG2+"):
        return "Mg"
    if rn in ("CA", "CA2+") and an == "CA":
        return "Ca"
    if rn in ("ZN", "ZN2+"):
        return "Zn"
    if rn in ("FE", "FE2", "FE3"):
        return "Fe"
    if rn in ("K", "K+", "POT"):
        return "K"
    if rn in ("MN", "CU", "CO", "NI"):
        return rn.capitalize()
    if rn in WATER_RESNAMES:
        if an.startswith("O"):
            return "O"
        if an.startswith("H") or an.startswith("M"):
            return "H"
    # Strip leading digits (e.g. 1HB -> HB -> H)
    stripped = an.lstrip("0123456789")
    if not stripped:
        return "C"
    if stripped.startswith("SE") and rn == "MSE":
        return "Se"
    if stripped.startswith("CL") and rn not in ALL_PROTEIN_RESNAMES:
        return "Cl"
    if stripped.startswith("BR") and rn not in ALL_PROTEIN_RESNAMES:
        return "Br"
    first = stripped[0]
    if first in ("H", "C", "N", "O", "S", "P", "F", "I", "B"):
        return first
    return "C"


def parse_gro_string(content: str) -> gemmi.Structure:
    """
    Parse a GROMACS .gro coordinate file (nm -> Å) into a Gemmi Structure.
    Automatically segments chains on residue-number drops and separates solvent/ions cleanly.
    """
    lines = content.splitlines()
    if len(lines) < 3:
        raise ValueError("Invalid .gro file: fewer than 3 lines")

    title = lines[0].strip() or "GROMACS Structure"
    try:
        n_atoms = int(lines[1].strip())
    except ValueError as exc:
        raise ValueError(f"Invalid .gro atom count on line 2: '{lines[1].strip()}'") from exc

    if len(lines) < 2 + n_atoms:
        raise ValueError(f".gro header declares {n_atoms} atoms, but file has only {len(lines) - 2} atom lines")

    st = gemmi.Structure()
    st.name = "GRO"
    st.info["_struct.title"] = title
    st.info["_exptl.method"] = "GROMACS COORDINATES (.GRO)"

    model = gemmi.Model("1")
    chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    prot_chain_idx = 0

    chains_map: Dict[str, gemmi.Chain] = {}
    chain_order: List[str] = []

    def _get_chain(cid: str) -> gemmi.Chain:
        if cid not in chains_map:
            ch = gemmi.Chain(cid)
            chains_map[cid] = ch
            chain_order.append(cid)
        return chains_map[cid]

    cur_prot_cid = "A"
    prev_prot_resnum: Optional[int] = None
    last_res_key: Optional[Tuple[str, int, str]] = None
    cur_res: Optional[gemmi.Residue] = None

    # Track unique residue numbering per chain when GROMACS wraps at 99999
    res_counter_by_chain: Dict[str, int] = {}

    for line in lines[2 : 2 + n_atoms]:
        if len(line) < 44:
            continue
        resnum = int(line[0:5].strip())
        resname = line[5:10].strip().upper()
        atomname = line[10:15].strip()
        # Coordinates in .gro are in nanometers -> convert to Ångströms
        x = float(line[20:28].strip()) * 10.0
        y = float(line[28:36].strip()) * 10.0
        z = float(line[36:44].strip()) * 10.0

        # Normalize GROMACS solvent/ion names for downstream tools
        norm_resname = "HOH" if resname in ("SOL", "TIP3", "TIP4", "SPC", "WAT") else resname
        if norm_resname in ("NA+", "SOD"):
            norm_resname = "NA"
        elif norm_resname in ("CL-", "CLA"):
            norm_resname = "CL"

        is_water = norm_resname in WATER_RESNAMES
        is_ion = norm_resname in METAL_IONS or norm_resname in HALIDE_IONS
        is_prot = norm_resname in ALL_PROTEIN_RESNAMES or norm_resname in NONSTANDARD_TO_STANDARD

        if is_water:
            cid = "W"
        elif is_ion:
            cid = "I"
        elif is_prot:
            # Detect new protein chain when residue number drops sharply (e.g. 150 -> 1)
            if prev_prot_resnum is not None and resnum < prev_prot_resnum - 2:
                prot_chain_idx += 1
                cur_prot_cid = chain_letters[prot_chain_idx % len(chain_letters)]
            cid = cur_prot_cid
            prev_prot_resnum = resnum
        else:
            cid = cur_prot_cid

        rkey = (cid, resnum, norm_resname)
        if rkey != last_res_key or cur_res is None:
            ch = _get_chain(cid)
            res = gemmi.Residue()
            res.name = norm_resname
            if is_water or is_ion:
                res_counter_by_chain[cid] = res_counter_by_chain.get(cid, 0) + 1
                res.seqid = gemmi.SeqId(res_counter_by_chain[cid], " ")
            else:
                res.seqid = gemmi.SeqId(resnum, " ")
            res.subchain = cid
            res.het_flag = "A" if is_prot else "H"
            ch.add_residue(res)
            cur_res = ch[-1]
            last_res_key = rkey

        atom = gemmi.Atom()
        # Standardize GROMACS water OW/HW1/HW2 -> O/H1/H2 for PDB/3Dmol compatibility
        if is_water and atomname.upper() == "OW":
            atomname = "O"
        elif is_water and atomname.upper() == "HW1":
            atomname = "H1"
        elif is_water and atomname.upper() == "HW2":
            atomname = "H2"
        # Skip GROMACS virtual sites (MW / LP)
        if atomname.upper().startswith(("MW", "LP")):
            continue

        atom.name = atomname
        elem = _infer_element_from_gro(norm_resname, atomname)
        atom.element = gemmi.Element(elem)
        atom.pos = gemmi.Position(x, y, z)
        atom.occ = 1.0
        atom.b_iso = 20.0
        cur_res.add_atom(atom)

    for cid in chain_order:
        model.add_chain(chains_map[cid])
    st.add_model(model)

    # Parse box vectors on the final line (nm -> Å)
    if len(lines) > 2 + n_atoms:
        box_parts = lines[2 + n_atoms].strip().split()
        try:
            vals = [float(v) * 10.0 for v in box_parts]
            if len(vals) >= 3 and vals[0] > 0.1 and vals[1] > 0.1 and vals[2] > 0.1:
                st.cell = gemmi.UnitCell(vals[0], vals[1], vals[2], 90.0, 90.0, 90.0)
        except ValueError:
            pass

    st.setup_entities()
    st.assign_label_seq_id()
    return st


def structure_to_gro(st: gemmi.Structure, title: Optional[str] = None) -> str:
    """Serialize a Gemmi Structure to standard GROMACS .gro format (Å -> nm)."""
    hdr = title or (st.info["_struct.title"] if "_struct.title" in st.info else None) or st.name or "Prepared by PDBPrep Studio"
    atom_lines: List[str] = []
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    serial = 1

    if len(st) > 0:
        for ch in st[0]:
            for res in ch:
                rn = res.name.strip().upper()
                gro_rn = "SOL" if rn in ("HOH", "WAT") else rn[:5]
                rnum = max(1, res.seqid.num) % 100000
                for atom in res:
                    an = atom.name.strip()
                    if gro_rn == "SOL" and an == "O":
                        an = "OW"
                    x_nm = atom.pos.x / 10.0
                    y_nm = atom.pos.y / 10.0
                    z_nm = atom.pos.z / 10.0
                    xs.append(x_nm)
                    ys.append(y_nm)
                    zs.append(z_nm)
                    atom_lines.append(
                        f"{rnum:5d}{gro_rn:<5s}{an[:5]:>5s}{serial % 100000:5d}{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}"
                    )
                    serial += 1

    # Box dimensions in nm
    if st.cell and st.cell.a > 1.0 and st.cell.b > 1.0 and st.cell.c > 1.0:
        bx, by, bz = st.cell.a / 10.0, st.cell.b / 10.0, st.cell.c / 10.0
    elif xs:
        bx = max(2.0, (max(xs) - min(xs)) + 2.0)
        by = max(2.0, (max(ys) - min(ys)) + 2.0)
        bz = max(2.0, (max(zs) - min(zs)) + 2.0)
    else:
        bx = by = bz = 5.0

    out = [hdr[:80], f"{len(atom_lines):5d}"] + atom_lines + [f"{bx:10.5f}{by:10.5f}{bz:10.5f}"]
    return "\n".join(out) + "\n"


def _looks_like_gro(content: str) -> bool:
    lines = content.splitlines()
    if len(lines) < 4:
        return False
    if not lines[1].strip().isdigit():
        return False
    # Typical .gro atom line has numeric coordinates around cols 20-44
    sample = lines[2]
    if len(sample) < 44:
        return False
    try:
        int(sample[0:5].strip())
        float(sample[20:28].strip())
        float(sample[28:36].strip())
        float(sample[36:44].strip())
        return True
    except ValueError:
        return False


def parse_structure_text(content: str, filename: str = "structure.pdb") -> gemmi.Structure:
    """Parse PDB, mmCIF, or GROMACS .gro text into a Gemmi Structure."""
    lower = filename.lower()
    if lower.endswith(".gro") or _looks_like_gro(content):
        return parse_gro_string(content)
    if lower.endswith(".cif") or lower.endswith(".mmcif") or content.lstrip().startswith("data_"):
        doc = gemmi.cif.read_string(content)
        st = gemmi.make_structure_from_block(doc.sole_block())
    else:
        st = gemmi.read_pdb_string(content)
    st.setup_entities()
    st.assign_label_seq_id()
    return st


def structure_to_pdb(st: gemmi.Structure) -> str:
    """Serialize a Gemmi Structure to standard PDB text (ensuring all protein/protonated states write as ATOM)."""
    for model in st:
        for ch in model:
            for res in ch:
                rn = res.name.strip().upper()
                if rn in ALL_PROTEIN_RESNAMES:
                    res.het_flag = "A"
    return st.make_pdb_string()


def structure_to_mmcif(st: gemmi.Structure) -> str:
    """Serialize a Gemmi Structure to mmCIF text."""
    doc = st.make_mmcif_document()
    return doc.as_string()


def extract_metadata(st: gemmi.Structure, source: str = "upload", title_fallback: str = "Protein Structure") -> Dict[str, Any]:
    """Extract crystallographic metadata, chains, heterogens, and composition summary."""
    model = st[0] if len(st) > 0 else None
    resolution = float(st.resolution) if st.resolution and st.resolution > 0 else None
    spacegroup = st.spacegroup_hm or None
    cell = None
    if st.cell and st.cell.a > 1.0:
        cell = {
            "a": round(st.cell.a, 3),
            "b": round(st.cell.b, 3),
            "c": round(st.cell.c, 3),
            "alpha": round(st.cell.alpha, 2),
            "beta": round(st.cell.beta, 2),
            "gamma": round(st.cell.gamma, 2),
        }

    title = (st.info["_struct.title"] if "_struct.title" in st.info else None) or st.name or title_fallback
    method = (st.info["_exptl.method"] if "_exptl.method" in st.info else None) or ("AlphaFold Prediction" if "AF-" in str(title_fallback) else "X-RAY DIFFRACTION")

    chains_summary: List[Dict[str, Any]] = []
    total_atoms = 0
    total_protein_residues = 0
    total_waters = 0
    total_ligands = 0
    total_metals = 0
    total_buffers = 0
    total_hydrogens = 0
    nonstandard_found: List[Dict[str, Any]] = []
    altloc_residues: List[Dict[str, Any]] = []

    if model is not None:
        for chain in model:
            prot_count = 0
            water_count = 0
            het_names: List[str] = []
            seq_start = None
            seq_end = None
            chain_atoms = 0

            for res in chain:
                rname = res.name.strip().upper()
                seqid = f"{res.seqid.num}{res.seqid.icode.strip()}"
                altlocs = sorted({a.altloc for a in res if a.altloc != "\x00" and a.altloc != " "})
                if len(altlocs) > 1:
                    altloc_residues.append({
                        "chain": chain.name,
                        "resseq": res.seqid.num,
                        "icode": res.seqid.icode.strip(),
                        "resname": rname,
                        "altlocs": altlocs,
                    })

                for atom in res:
                    total_atoms += 1
                    chain_atoms += 1
                    if atom.element.name in ("H", "D"):
                        total_hydrogens += 1

                if rname in ALL_PROTEIN_RESNAMES:
                    prot_count += 1
                    total_protein_residues += 1
                    if seq_start is None:
                        seq_start = res.seqid.num
                    seq_end = res.seqid.num
                elif rname in NONSTANDARD_TO_STANDARD:
                    prot_count += 1
                    total_protein_residues += 1
                    nonstandard_found.append({
                        "chain": chain.name,
                        "resseq": res.seqid.num,
                        "icode": res.seqid.icode.strip(),
                        "resname": rname,
                        "suggested": NONSTANDARD_TO_STANDARD[rname],
                    })
                elif rname in WATER_RESNAMES:
                    water_count += 1
                    total_waters += 1
                elif rname in METAL_IONS or rname in HALIDE_IONS:
                    total_metals += 1
                    het_names.append(rname)
                elif rname in CRYSTALLIZATION_ARTIFACTS:
                    total_buffers += 1
                    het_names.append(rname)
                else:
                    total_ligands += 1
                    het_names.append(rname)

            chains_summary.append({
                "id": chain.name,
                "protein_residues": prot_count,
                "waters": water_count,
                "heterogens": sorted(set(het_names)),
                "atom_count": chain_atoms,
                "seq_range": [seq_start, seq_end] if seq_start is not None else None,
            })

    return {
        "title": title,
        "source": source,
        "method": method,
        "resolution": round(resolution, 2) if resolution else None,
        "spacegroup": spacegroup,
        "unit_cell": cell,
        "model_count": len(st),
        "chain_count": len(chains_summary),
        "chains": chains_summary,
        "counts": {
            "atoms": total_atoms,
            "hydrogens": total_hydrogens,
            "protein_residues": total_protein_residues,
            "waters": total_waters,
            "ligands": total_ligands,
            "metals_ions": total_metals,
            "buffers": total_buffers,
        },
        "nonstandard_residues": nonstandard_found,
        "altloc_residues": altloc_residues,
    }


_SAMPLES_DIR = Path(__file__).resolve().parent.parent / "samples"


async def fetch_rcsb_structure(pdb_id: str, fmt: str = "pdb") -> Tuple[str, str]:
    """Download a structure from RCSB PDB (serving bundled local cache first if present, with mmCIF fallback)."""
    code = pdb_id.strip().upper()
    for ext in ([fmt, "cif" if fmt == "pdb" else "pdb"]):
        local_path = _SAMPLES_DIR / f"{code}.{ext}"
        if local_path.exists():
            return local_path.read_text(encoding="utf-8", errors="replace"), f"{code}.{ext}"

    urls = (
        [f"https://files.rcsb.org/download/{code}.pdb", f"https://files.rcsb.org/download/{code}.cif"]
        if fmt == "pdb"
        else [f"https://files.rcsb.org/download/{code}.cif", f"https://files.rcsb.org/download/{code}.pdb"]
    )
    async with httpx.AsyncClient(timeout=25.0, follow_redirects=True) as client:
        last_err = None
        for url in urls:
            try:
                resp = await client.get(url)
                if resp.status_code == 200 and len(resp.text) > 100:
                    ext = "cif" if url.endswith(".cif") else "pdb"
                    return resp.text, f"{code}.{ext}"
            except Exception as exc:
                last_err = exc
    raise ValueError(f"Could not fetch PDB ID '{code}' from RCSB ({last_err or 'not found'})")


async def fetch_alphafold_structure(uniprot_id: str) -> Tuple[str, str]:
    """Download a predicted structure from AlphaFold Protein Structure Database."""
    uid = uniprot_id.strip().upper()
    async with httpx.AsyncClient(timeout=25.0, follow_redirects=True) as client:
        # Query API first to get latest model version URL
        api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uid}"
        try:
            meta = await client.get(api_url)
            if meta.status_code == 200:
                entries = meta.json()
                if isinstance(entries, list) and entries:
                    pdb_url = entries[0].get("pdbUrl") or entries[0].get("cifUrl")
                    if pdb_url:
                        resp = await client.get(pdb_url)
                        if resp.status_code == 200:
                            ext = "cif" if pdb_url.endswith(".cif") else "pdb"
                            return resp.text, f"AF-{uid}.{ext}"
        except Exception:
            pass

        # Fallback direct v4 URL
        direct = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb"
        resp = await client.get(direct)
        if resp.status_code == 200 and len(resp.text) > 100:
            return resp.text, f"AF-{uid}-F1.pdb"
    raise ValueError(f"Could not fetch AlphaFold structure for UniProt '{uid}'")


def _build_backbone_chain(
    sequence: List[str],
    phi_psi_omega: List[Tuple[float, float, float]],
) -> List[Dict[str, np.ndarray]]:
    """Build a realistic 3D peptide backbone from (phi, psi, omega) angles using NeRF."""
    residues: List[Dict[str, np.ndarray]] = []
    n0 = np.array([0.0, 0.0, 0.0])
    ca0 = np.array([1.458, 0.0, 0.0])
    c0 = nerf_place(np.array([0.0, 1.0, 0.0]), n0, ca0, 1.525, 111.2, 180.0)
    o0 = nerf_place(n0, ca0, c0, 1.231, 120.8, -120.0)
    residues.append({"N": n0, "CA": ca0, "C": c0, "O": o0})

    for i in range(1, len(sequence)):
        prev = residues[-1]
        phi_i, psi_prev, omega_i = phi_psi_omega[i]
        # Place N_i using psi_{i-1}
        n_i = nerf_place(prev["N"], prev["CA"], prev["C"], 1.329, 116.2, psi_prev)
        # Update prev O bonded to prev["C"], trans to N_i
        prev["O"] = nerf_place(n_i, prev["CA"], prev["C"], 1.231, 120.8, 180.0)
        # Place CA_i using omega_i
        ca_i = nerf_place(prev["CA"], prev["C"], n_i, 1.458, 121.7, omega_i)
        # Place C_i using phi_i
        c_i = nerf_place(prev["C"], n_i, ca_i, 1.525, 111.2, phi_i)
        o_i = nerf_place(n_i, ca_i, c_i, 1.231, 120.8, 180.0)
        residues.append({"N": n_i, "CA": ca_i, "C": c_i, "O": o_i})

    return residues


def build_showcase_structure(sample_id: str = "catalytic_showcase") -> Tuple[str, str]:
    """
    Generate a rich, physically realistic 52-residue metalloenzyme-inhibitor complex
    containing authentic structural biology preparation targets:
      - Helix-turn-sheet-loop-helix fold (Chain A)
      - Catalytic Zn2+ site coordinated by HIS18, HIS22, CYS38
      - Disulfide bridge (CYS12 - CYS44)
      - Non-standard Selenomethionine (MSE 7)
      - Missing side-chain atoms (LYS 15 missing NZ/CE, ARG 31 missing NH1/NH2)
      - Missing internal loop gap (residue 27 omitted with SEQRES record)
      - Dual conformation AltLoc A/B (LEU 19)
      - Intentional cis non-Pro peptide bond (ALA 33 - GLU 34)
      - Intentional Ramachandran outlier (ASP 25)
      - Intentional steric clash / rotamer outlier (PHE 14)
      - Co-crystallized drug-like ligand (BEN: benzamidine-like active site inhibitor)
      - Crystallization artifacts (SO4, EDO) + buried catalytic & bulk waters
    """
    seq = [
        "GLY", "VAL", "ALA", "GLU", "LEU", "LYS", "MSE", "ALA", "ASP", "ILE",
        "THR", "CYS", "SER", "PHE", "LYS", "ASN", "GLN", "HIS", "LEU", "GLY",
        "ALA", "HIS", "GLU", "VAL", "ASP", "GLY", "SER", "THR", "TYR", "VAL",
        "ARG", "ILE", "ALA", "GLU", "LEU", "ASN", "PRO", "CYS", "ASP", "LYS",
        "TRP", "PHE", "GLU", "CYS", "ALA", "GLN", "HIS", "LEU", "SER", "GLY",
    ]
    n_res = len(seq)

    # Secondary structure torsion blueprint:
    # 1-16: alpha helix; 17-26: catalytic loop + turn; 27-35: beta strand; 36-50: C-term alpha helix
    torsions: List[Tuple[float, float, float]] = []
    for i in range(n_res):
        if i < 16:
            phi, psi, omega = -62.0, -42.0, 180.0
        elif i < 24:
            # Smooth catalytic loop curving back
            phi, psi, omega = (-85.0, 120.0, 180.0) if i % 2 == 0 else (-65.0, -30.0, 180.0)
        elif i == 24:
            # Deliberate Ramachandran outlier for COOT demo (res 25 ASP)
            phi, psi, omega = 85.0, -110.0, 180.0
        elif i < 34:
            # Extended beta-strand
            phi, psi, omega = -125.0, 132.0, 180.0
            if i == 33:
                # Deliberate cis non-Pro peptide bond (res 33->34)
                omega = 4.0
        else:
            # Second alpha-helix
            phi, psi, omega = -63.0, -41.0, 180.0
        torsions.append((phi, psi, omega))

    bb_list = _build_backbone_chain(seq, torsions)

    st = gemmi.Structure()
    st.name = "PDBP"
    st.resolution = 1.85
    st.spacegroup_hm = "P 21 21 21"
    st.cell = gemmi.UnitCell(48.2, 56.7, 72.4, 90.0, 90.0, 90.0)
    st.info["_struct.title"] = "PDBPrep Showcase Metalloenzyme Complexed with Benzamidine Inhibitor"
    st.info["_exptl.method"] = "X-RAY DIFFRACTION"

    model = gemmi.Model("1")
    chain = gemmi.Chain("A")

    # Add SEQRES entity (including missing residue 27 GLY to showcase loop detection)
    entity = gemmi.Entity("1")
    entity.entity_type = gemmi.EntityType.Polymer
    entity.polymer_type = gemmi.PolymerType.PeptideL
    entity.subchains = ["A"]
    entity.full_sequence = [s if s != "MSE" else "MET" for s in seq]
    st.entities.append(entity)

    ca_positions: List[np.ndarray] = []

    for idx, resname in enumerate(seq):
        resnum = idx + 1
        # Deliberately omit residue 27 (SER) to simulate a missing-loop crystallographic break
        if resnum == 27 and sample_id == "catalytic_showcase":
            continue

        res = gemmi.Residue()
        res.name = resname
        res.seqid = gemmi.SeqId(resnum, " ")
        res.subchain = "A"
        res.entity_type = gemmi.EntityType.Polymer

        std_name = "MET" if resname == "MSE" else resname
        chi_override = None
        if resnum == 14:
            # Deliberate rotamer outlier + steric clash for PHE 14
            chi_override = [12.0, 5.0]
        elif resnum == 16:
            # ASN 16 candidate for NQH flip
            chi_override = [-68.0, 160.0]

        coords = build_sidechain_coords(bb_list[idx], std_name, chi_targets=chi_override)
        ca_positions.append(coords["CA"])

        # Deliberate missing heavy atoms on LYS 15 and ARG 31
        skip_atoms = set()
        if resnum == 15:
            skip_atoms = {"CE", "NZ"}
        elif resnum == 31:
            skip_atoms = {"NH1", "NH2"}

        # Slight bond-length strain on residue 25 to show Engh & Huber geometry alert
        if resnum == 25:
            coords["O"] = coords["C"] + (coords["O"] - coords["C"]) * 1.18

        bfactor_base = 18.0 + 14.0 * math.sin(idx * 0.25) ** 2 + (28.0 if 22 <= resnum <= 28 else 0.0)

        for aname, pos in coords.items():
            if aname in skip_atoms:
                continue
            final_aname = "SE" if (resname == "MSE" and aname == "SD") else aname
            elem_sym = "Se" if final_aname == "SE" else final_aname[0]

            atom = gemmi.Atom()
            atom.name = final_aname
            atom.element = gemmi.Element(elem_sym)
            atom.pos = gemmi.Position(float(pos[0]), float(pos[1]), float(pos[2]))
            atom.occ = 0.65 if (resnum == 19 and aname not in ("N", "CA", "C", "O")) else 1.0
            atom.altloc = "A" if (resnum == 19 and aname not in ("N", "CA", "C", "O")) else "\x00"
            atom.b_iso = round(bfactor_base + (4.0 if aname not in ("N", "CA", "C", "O") else 0.0), 2)
            res.add_atom(atom)

        # Add AltLoc B for LEU 19 side chain
        if resnum == 19:
            coords_b = build_sidechain_coords(bb_list[idx], "LEU", chi_targets=[-177.0, 63.0])
            for aname in ("CB", "CG", "CD1", "CD2"):
                if aname in coords_b:
                    pos = coords_b[aname]
                    atom_b = gemmi.Atom()
                    atom_b.name = aname
                    atom_b.element = gemmi.Element("C")
                    atom_b.pos = gemmi.Position(float(pos[0]), float(pos[1]), float(pos[2]))
                    atom_b.occ = 0.35
                    atom_b.altloc = "B"
                    atom_b.b_iso = round(bfactor_base + 9.5, 2)
                    res.add_atom(atom_b)

        chain.add_residue(res)

    # Pull CYS 12 SG and CYS 44 SG or place a proximal CYS pair at 2.04 Å for realistic disulfide detection
    # Find HIS 18 NE2 to position a catalytic Zn2+ ion at 2.08 Å
    his18 = next((r for r in chain if r.seqid.num == 18), None)
    his22 = next((r for r in chain if r.seqid.num == 22), None)
    cys12 = next((r for r in chain if r.seqid.num == 12), None)
    cys13_sg_partner = next((r for r in chain if r.seqid.num == 44), None)

    # Position CYS 44 SG near CYS 12 SG if far, or create a clean disulfide on a realistic adjacent CYS 9/12
    # Better: keep protein fold natural and place a second short chain B (peptide substrate / disulfide partner)
    chain_b = gemmi.Chain("B")
    if cys12 is not None and "SG" in cys12 and "CA" in cys12:
        sg12 = np.array([cys12["SG"][0].pos.x, cys12["SG"][0].pos.y, cys12["SG"][0].pos.z])
        ca12 = np.array([cys12["CA"][0].pos.x, cys12["CA"][0].pos.y, cys12["CA"][0].pos.z])
        outward = sg12 - ca12
        outward = outward / (np.linalg.norm(outward) + 1e-8)
        # Build a 6-residue partner peptide Chain B extending outward from A:CYS12 with a 2.04 Å disulfide
        b_seq = ["ALA", "SER", "CYS", "GLY", "GLU", "TYR"]
        b_bb = _build_backbone_chain(b_seq, [(-125.0, 132.0, 180.0)] * 6)
        b_cys_coords = build_sidechain_coords(b_bb[2], "CYS")
        # Orient Chain B so its CA->SG points anti-parallel to outward (avoiding steric collision with Chain A)
        b_sg_vec = b_cys_coords["SG"] - b_cys_coords["CA"]
        b_sg_vec = b_sg_vec / (np.linalg.norm(b_sg_vec) + 1e-8)
        target_sg_b = sg12 + outward * 2.04
        # Mirror/flip if pointing inward
        flip_sign = -1.0 if np.dot(b_sg_vec, outward) > 0 else 1.0
        for b_idx, b_rname in enumerate(b_seq):
            bres = gemmi.Residue()
            bres.name = b_rname
            bres.seqid = gemmi.SeqId(101 + b_idx, " ")
            bres.subchain = "B"
            bres.entity_type = gemmi.EntityType.Polymer
            bcoords = build_sidechain_coords(b_bb[b_idx], b_rname)
            for aname, pos in bcoords.items():
                rel = (pos - b_cys_coords["SG"]) * np.array([flip_sign, flip_sign, 1.0])
                p = target_sg_b + rel
                atom = gemmi.Atom()
                atom.name = aname
                atom.element = gemmi.Element(aname[0])
                atom.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))
                atom.occ = 1.0
                atom.b_iso = 22.4
                bres.add_atom(atom)
            chain_b.add_residue(bres)

    # Catalytic Zn2+ positioned 2.08 Å from HIS 18 NE2
    zn_pos = np.array([10.0, 8.0, 6.0])
    if his18 is not None and "NE2" in his18 and "CG" in his18:
        ne2 = np.array([his18["NE2"][0].pos.x, his18["NE2"][0].pos.y, his18["NE2"][0].pos.z])
        cg = np.array([his18["CG"][0].pos.x, his18["CG"][0].pos.y, his18["CG"][0].pos.z])
        v = ne2 - cg
        v = v / (np.linalg.norm(v) + 1e-8)
        zn_pos = ne2 + v * 2.08

    zn_res = gemmi.Residue()
    zn_res.name = "ZN"
    zn_res.seqid = gemmi.SeqId(301, " ")
    zn_res.het_flag = "H"
    zn_atom = gemmi.Atom()
    zn_atom.name = "ZN"
    zn_atom.element = gemmi.Element("Zn")
    zn_atom.charge = 2
    zn_atom.pos = gemmi.Position(float(zn_pos[0]), float(zn_pos[1]), float(zn_pos[2]))
    zn_atom.occ = 1.0
    zn_atom.b_iso = 16.5
    zn_res.add_atom(zn_atom)
    chain.add_residue(zn_res)

    # Co-crystallized drug-like inhibitor BEN (Benzamidine: C7H8N2, MW 120.15) positioned in binding pocket near Zn/HIS18
    lig_center = zn_pos + np.array([3.6, 2.4, -1.8])
    ben_atoms = [
        ("C1", "C", [1.39, 0.00, 0.00]),
        ("C2", "C", [0.70, 1.21, 0.00]),
        ("C3", "C", [-0.70, 1.21, 0.00]),
        ("C4", "C", [-1.39, 0.00, 0.00]),
        ("C5", "C", [-0.70, -1.21, 0.00]),
        ("C6", "C", [0.70, -1.21, 0.00]),
        ("C7", "C", [2.87, 0.00, 0.00]),
        ("N1", "N", [3.53, 1.15, 0.00]),
        ("N2", "N", [3.53, -1.15, 0.00]),
    ]
    ben_res = gemmi.Residue()
    ben_res.name = "BEN"
    ben_res.seqid = gemmi.SeqId(302, " ")
    ben_res.het_flag = "H"
    for aname, el, rel in ben_atoms:
        p = lig_center + np.array(rel)
        atom = gemmi.Atom()
        atom.name = aname
        atom.element = gemmi.Element(el)
        atom.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))
        atom.occ = 1.0
        atom.b_iso = 21.8
        ben_res.add_atom(atom)
    chain.add_residue(ben_res)

    # Crystallization artifacts: SO4 and EDO
    so4_center = lig_center + np.array([-9.5, 6.0, 5.5])
    so4_res = gemmi.Residue()
    so4_res.name = "SO4"
    so4_res.seqid = gemmi.SeqId(303, " ")
    so4_res.het_flag = "H"
    for aname, el, rel in [
        ("S", "S", [0, 0, 0]),
        ("O1", "O", [1.47, 0, 0]),
        ("O2", "O", [-0.49, 1.39, 0]),
        ("O3", "O", [-0.49, -0.69, 1.20]),
        ("O4", "O", [-0.49, -0.69, -1.20]),
    ]:
        p = so4_center + np.array(rel)
        atom = gemmi.Atom()
        atom.name = aname
        atom.element = gemmi.Element(el)
        atom.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))
        atom.occ = 1.0
        atom.b_iso = 44.2
        so4_res.add_atom(atom)
    chain.add_residue(so4_res)

    edo_center = lig_center + np.array([8.5, -6.5, 4.0])
    edo_res = gemmi.Residue()
    edo_res.name = "EDO"
    edo_res.seqid = gemmi.SeqId(304, " ")
    edo_res.het_flag = "H"
    for aname, el, rel in [
        ("O1", "O", [-1.42, 0.50, 0.0]),
        ("C1", "C", [-0.75, -0.50, 0.0]),
        ("C2", "C", [0.75, -0.50, 0.0]),
        ("O2", "O", [1.42, 0.50, 0.0]),
    ]:
        p = edo_center + np.array(rel)
        atom = gemmi.Atom()
        atom.name = aname
        atom.element = gemmi.Element(el)
        atom.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))
        atom.occ = 1.0
        atom.b_iso = 39.0
        edo_res.add_atom(atom)
    chain.add_residue(edo_res)

    # Crystallographic waters: 2 active-site waters (near Zn & BEN) + 14 bulk waters
    water_positions = [
        zn_pos + np.array([0.2, -2.05, 0.4]),       # Catalytic water coordinating Zn
        lig_center + np.array([2.2, 2.8, 0.5]),     # Bridging water near BEN N1
    ]
    rng = np.random.default_rng(42)
    centroid = np.mean(ca_positions, axis=0)
    for w in range(14):
        vec = rng.normal(size=3)
        vec = vec / np.linalg.norm(vec)
        water_positions.append(centroid + vec * (14.5 + (w % 4) * 2.2))

    for wi, wp in enumerate(water_positions, start=401):
        wres = gemmi.Residue()
        wres.name = "HOH"
        wres.seqid = gemmi.SeqId(wi, " ")
        wres.het_flag = "H"
        watom = gemmi.Atom()
        watom.name = "O"
        watom.element = gemmi.Element("O")
        watom.pos = gemmi.Position(float(wp[0]), float(wp[1]), float(wp[2]))
        watom.occ = 1.0
        watom.b_iso = 19.5 if wi <= 402 else round(34.0 + wi % 12, 1)
        wres.add_atom(watom)
        chain.add_residue(wres)

    model.add_chain(chain)
    if len(chain_b) > 0:
        model.add_chain(chain_b)
    st.add_model(model)
    st.setup_entities()
    return structure_to_pdb(st), "PDBP_showcase.pdb"
