"""
Ligand, Cofactor & Metal Site Analysis Engine powered by RDKit + Gemmi:
- Extracts & classifies all heterogens (Ligand, Cofactor, Metal, Crystallization Buffer)
- Computes RDKit physicochemical & Lipinski Rule-of-5 properties + crisp 2D SVG depiction
- Maps 3D binding pocket interactions (H-bonds, hydrophobic contacts, metal coordination)
- Computes metal coordination sphere geometry
- Generates 3D ETKDGv3 conformers & SDF export
- Computes automatic docking grid box around ligands or active pocket
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import gemmi
import numpy as np

from .chemistry_data import (
    ALL_PROTEIN_RESNAMES,
    COFACTORS,
    CRYSTALLIZATION_ARTIFACTS,
    HALIDE_IONS,
    METAL_IONS,
    NONSTANDARD_TO_STANDARD,
    WATER_RESNAMES,
)
from .geometry_utils import distance

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
except Exception:
    HAS_RDKIT = False


def _is_protein_or_mod(rname: str) -> bool:
    return rname in ALL_PROTEIN_RESNAMES or rname in NONSTANDARD_TO_STANDARD


def _residue_to_pdb_block(chain_id: str, res: gemmi.Residue) -> str:
    sub = gemmi.Structure()
    m = gemmi.Model("1")
    c = gemmi.Chain(chain_id)
    c.add_residue(res.clone())
    m.add_chain(c)
    sub.add_model(m)
    return sub.make_pdb_string()


def _render_mol_svg(mol: "Chem.Mol", width: int = 280, height: int = 190) -> str:
    try:
        m2d = Chem.Mol(mol)
        AllChem.Compute2DCoords(m2d)
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = False
        opts.bondLineWidth = 1.6
        opts.padding = 0.12
        drawer.DrawMolecule(m2d)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return ""


def _analyze_rdkit_ligand(pdb_block: str) -> Dict[str, Any]:
    if not HAS_RDKIT:
        return {}
    try:
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=True, removeHs=False)
        if mol is None:
            mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
            if mol is not None:
                try:
                    Chem.SanitizeMol(mol)
                except Exception:
                    pass
        if mol is None:
            return {}
        mol_no_h = Chem.RemoveHs(mol, sanitize=False)
        mw = round(float(Descriptors.MolWt(mol)), 2)
        clogp = round(float(Descriptors.MolLogP(mol)), 2)
        tpsa = round(float(Descriptors.TPSA(mol)), 2)
        hbd = int(Lipinski.NumHDonors(mol))
        hba = int(Lipinski.NumHAcceptors(mol))
        rot = int(Lipinski.NumRotatableBonds(mol))
        rings = int(rdMolDescriptors.CalcNumRings(mol))
        formula = rdMolDescriptors.CalcMolFormula(mol)
        smiles = Chem.MolToSmiles(mol_no_h)
        violations = int(mw > 500) + int(clogp > 5) + int(hbd > 5) + int(hba > 10)
        svg = _render_mol_svg(mol_no_h)
        return {
            "formula": formula,
            "smiles": smiles,
            "mw": mw,
            "clogp": clogp,
            "tpsa": tpsa,
            "hbd": hbd,
            "hba": hba,
            "rotatable_bonds": rot,
            "rings": rings,
            "lipinski_violations": violations,
            "drug_like": violations <= 1,
            "svg": svg,
        }
    except Exception:
        return {}


def analyze_ligands_and_metals(st: gemmi.Structure, contact_cutoff: float = 4.2) -> Dict[str, Any]:
    """Detect and analyze all ligands, cofactors, metals, buffers, and binding-site contacts."""
    if len(st) == 0:
        return {"ligands": [], "metals": [], "buffers": [], "docking_boxes": []}

    model = st[0]
    protein_atoms: List[Tuple[str, int, str, str, str, np.ndarray]] = []
    for chain in model:
        for res in chain:
            rn = res.name.strip().upper()
            if _is_protein_or_mod(rn):
                for a in res:
                    if a.element.name not in ("H", "D"):
                        protein_atoms.append((
                            chain.name,
                            res.seqid.num,
                            res.seqid.icode.strip(),
                            rn,
                            a.name.strip(),
                            np.array([a.pos.x, a.pos.y, a.pos.z]),
                        ))

    ligands = []
    metals = []
    buffers = []
    docking_boxes = []

    for chain in model:
        for res in chain:
            rn = res.name.strip().upper()
            if _is_protein_or_mod(rn) or rn in WATER_RESNAMES:
                continue

            ic = res.seqid.icode.strip()
            lig_id = f"{chain.name}:{rn}:{res.seqid.num}{ic}"
            heavy_atoms = [a for a in res if a.element.name not in ("H", "D")]
            if not heavy_atoms:
                continue
            pts = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in heavy_atoms], dtype=np.float64)
            center = np.mean(pts, axis=0)
            span = np.max(pts, axis=0) - np.min(pts, axis=0)

            if rn in METAL_IONS or rn in HALIDE_IONS:
                # Coordination sphere analysis
                coordinators = []
                for pch, prs, pic, prn, pan, ppos in protein_atoms:
                    if pan[0] in ("N", "O", "S"):
                        d = distance(center, ppos)
                        if d <= 2.9:
                            coordinators.append({
                                "chain": pch,
                                "resseq": prs,
                                "resname": prn,
                                "atom": pan,
                                "key": f"{pch}:{prs}{pic}",
                                "distance_A": round(d, 2),
                            })
                coordinators.sort(key=lambda x: x["distance_A"])
                cn = len(coordinators)
                geom = {1: "Monodentate", 2: "Linear/Bent", 3: "Trigonal", 4: "Tetrahedral/Square Planar", 5: "Trigonal Bipyramidal", 6: "Octahedral"}.get(cn, f"CN={cn}")
                metals.append({
                    "id": lig_id,
                    "chain": chain.name,
                    "resseq": res.seqid.num,
                    "resname": rn,
                    "element": heavy_atoms[0].element.name,
                    "center": [round(float(x), 3) for x in center],
                    "coordination_number": cn,
                    "geometry": geom,
                    "coordinators": coordinators,
                })
                continue

            if rn in CRYSTALLIZATION_ARTIFACTS:
                buffers.append({
                    "id": lig_id,
                    "chain": chain.name,
                    "resseq": res.seqid.num,
                    "resname": rn,
                    "heavy_atoms": len(heavy_atoms),
                    "center": [round(float(x), 3) for x in center],
                })
                continue

            # Drug-like Ligand or Cofactor
            cat = "cofactor" if rn in COFACTORS else "ligand"
            pdb_block = _residue_to_pdb_block(chain.name, res)
            rd_props = _analyze_rdkit_ligand(pdb_block)

            # Binding pocket contacts
            contacts: List[Dict[str, Any]] = []
            seen_pair = set()
            for a in heavy_atoms:
                la_pos = np.array([a.pos.x, a.pos.y, a.pos.z])
                la_el = a.element.name.upper()
                for pch, prs, pic, prn, pan, ppos in protein_atoms:
                    d = distance(la_pos, ppos)
                    if d > contact_cutoff:
                        continue
                    pa_el = pan[0].upper()
                    if la_el in ("N", "O") and pa_el in ("N", "O") and d <= 3.55:
                        itype = "hbond"
                    elif la_el == "C" and pa_el == "C" and d <= 4.0:
                        itype = "hydrophobic"
                    else:
                        itype = "vdw"
                    pair_k = (pch, prs, itype)
                    if pair_k not in seen_pair:
                        seen_pair.add(pair_k)
                        contacts.append({
                            "chain": pch,
                            "resseq": prs,
                            "icode": pic,
                            "resname": prn,
                            "key": f"{pch}:{prs}{pic}",
                            "protein_atom": pan,
                            "ligand_atom": a.name.strip(),
                            "distance_A": round(d, 2),
                            "type": itype,
                        })
            contacts.sort(key=lambda x: x["distance_A"])

            box = {
                "ligand_id": lig_id,
                "label": f"{rn} ({chain.name}:{res.seqid.num})",
                "center": [round(float(x), 2) for x in center],
                "size": [round(float(max(16.0, span[k] + 10.0)), 1) for k in range(3)],
            }
            docking_boxes.append(box)

            ligands.append({
                "id": lig_id,
                "chain": chain.name,
                "resseq": res.seqid.num,
                "icode": ic,
                "resname": rn,
                "category": cat,
                "heavy_atoms": len(heavy_atoms),
                "mean_b": round(float(np.mean([a.b_iso for a in heavy_atoms])), 2),
                "center": [round(float(x), 3) for x in center],
                "docking_box": box,
                "contacts": contacts[:40],
                **rd_props,
            })

    # Fallback global docking box if no co-crystallized ligand present
    if not docking_boxes and protein_atoms:
        p_pts = np.array([x[5] for x in protein_atoms], dtype=np.float64)
        p_center = np.mean(p_pts, axis=0)
        docking_boxes.append({
            "ligand_id": "global_protein",
            "label": "Protein Centroid (Blind Docking)",
            "center": [round(float(x), 2) for x in p_center],
            "size": [24.0, 24.0, 24.0],
        })

    return {
        "ligands": ligands,
        "metals": metals,
        "buffers": buffers,
        "docking_boxes": docking_boxes,
    }


def generate_ligand_conformers_sdf(st: gemmi.Structure, ligand_id: str, num_confs: int = 5) -> Tuple[str, List[Dict[str, Any]]]:
    """Generate ETKDGv3 3D conformers for a specified ligand and return SDF + energy ranking."""
    if not HAS_RDKIT:
        raise RuntimeError("RDKit is required for conformer generation")
    if len(st) == 0:
        raise ValueError("Empty structure")

    target_res = None
    target_chain = ""
    for chain in st[0]:
        for res in chain:
            lid = f"{chain.name}:{res.name.strip().upper()}:{res.seqid.num}{res.seqid.icode.strip()}"
            if lid == ligand_id:
                target_res = res
                target_chain = chain.name
                break
    if target_res is None:
        raise ValueError(f"Ligand '{ligand_id}' not found")

    pdb_block = _residue_to_pdb_block(target_chain, target_res)
    mol = Chem.MolFromPDBBlock(pdb_block, sanitize=True, removeHs=False)
    if mol is None:
        raise ValueError("Could not sanitize ligand structure with RDKit")

    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    cids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=max(1, min(20, num_confs)), params=params))
    if not cids:
        raise ValueError("ETKDGv3 could not embed conformers for this ligand")

    energies = []
    for cid in cids:
        try:
            props = AllChem.MMFFGetMoleculeProperties(mol_h)
            if props is not None:
                ff = AllChem.MMFFGetMoleculeForceField(mol_h, props, confId=cid)
                ff.Minimize(maxIts=200)
                e = float(ff.CalcEnergy())
            else:
                ff = AllChem.UFFGetMoleculeForceField(mol_h, confId=cid)
                ff.Minimize(maxIts=200)
                e = float(ff.CalcEnergy())
        except Exception:
            e = 0.0
        energies.append((cid, round(e, 2)))

    energies.sort(key=lambda x: x[1])
    min_e = energies[0][1] if energies else 0.0

    sdf_blocks = []
    summary = []
    for rank, (cid, energy) in enumerate(energies, start=1):
        mol_h.SetProp("_Name", f"{target_res.name.strip()}_conf_{rank}")
        mol_h.SetProp("MMFF94_Energy_kcal_mol", f"{energy:.2f}")
        mol_h.SetProp("Delta_E_kcal_mol", f"{energy - min_e:.2f}")
        sdf_blocks.append(Chem.MolToMolBlock(mol_h, confId=cid) + "\n$$$$\n")
        summary.append({
            "rank": rank,
            "conformer_id": cid,
            "energy_kcal_mol": energy,
            "delta_e": round(energy - min_e, 2),
        })

    return "".join(sdf_blocks), summary
