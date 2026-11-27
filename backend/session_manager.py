"""
Multi-Structure Session Manager with COOT-style Undo / Redo Snapshot Stack.
Every structural modification pushes a checkpoint with a human-readable audit log.
"""

from __future__ import annotations

import json
import os
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import gemmi

_SESSIONS_DIR = Path(
    os.environ.get("PANDAPROTEIN_SESSION_DIR") or (Path.home() / ".pandaprotein" / "sessions")
)
_SESSIONS_DIR.mkdir(parents=True, exist_ok=True)

from .engines.adapt_engine import (
    build_preparation_tree,
    build_provenance_record,
    build_remark_850_lines,
    detect_missing_regions,
    list_biological_assemblies,
    lookup_rcsb_entry_ligands,
    md5_hex,
)
from .engines.chemistry_data import AA_3_TO_1
from .engines.coot_engine import _get_pos, analyze_structure_quality
from .engines.export_engine import compute_readiness_audit, estimate_solvation_system
from .engines.geometry_utils import canonical_resname, distance
from .engines.ligand_engine import analyze_ligands_and_metals
from .engines.pdbfixer_engine import diagnose_pdbfixer
from .engines.pocket_engine import detect_binding_pockets
from .engines.propka_engine import analyze_propka
from .engines.structure_io import extract_metadata, structure_to_pdb


def _default_docking_prefs() -> Dict[str, Any]:
    return {
        "residue_roles": {},
        "chain_roles": {},
        "ion_charges": {},
        "ion_substitutions": {},
        "include_nucleic": True,
    }


@dataclass
class Checkpoint:
    pdb_text: str
    action: str
    timestamp: str


@dataclass
class StructureSession:
    session_id: str
    name: str
    source: str
    structure: gemmi.Structure
    ph: float = 7.4
    propka_overrides: Dict[str, str] = field(default_factory=dict)
    his_default: str = "auto"
    respect_metal: bool = True
    respect_disulfides: bool = True
    undo_stack: List[Checkpoint] = field(default_factory=list)
    redo_stack: List[Checkpoint] = field(default_factory=list)
    history_log: List[str] = field(default_factory=list)
    baseline: Optional[Dict[str, Any]] = field(default=None)
    original_pdb: Optional[str] = field(default=None)
    raw_input_structure: Optional[gemmi.Structure] = field(default=None, repr=False)
    input_md5: str = ""
    active_assembly_id: str = "ASU"
    symmetry_copy_chains: List[str] = field(default_factory=list)
    docking_prefs: Dict[str, Any] = field(default_factory=_default_docking_prefs)
    rebuilt_gaps: List[Dict[str, Any]] = field(default_factory=list)
    refused_gaps: Dict[str, str] = field(default_factory=dict)
    rcsb_ligands: Optional[Dict[str, Any]] = field(default=None)
    _cached_state: Optional[Dict[str, Any]] = field(default=None, repr=False)

    def invalidate_cache(self) -> None:
        self._cached_state = None
        self.save_to_disk()

    def save_to_disk(self) -> None:
        """Persist session coordinates & state to disk so uvicorn reloads never lose active sessions."""
        try:
            _SESSIONS_DIR.mkdir(parents=True, exist_ok=True)
            pdb_path = _SESSIONS_DIR / f"{self.session_id}.pdb"
            orig_path = _SESSIONS_DIR / f"{self.session_id}_orig.pdb"
            meta_path = _SESSIONS_DIR / f"{self.session_id}.json"
            pdb_path.write_text(structure_to_pdb(self.structure), encoding="utf-8")
            if self.original_pdb and not orig_path.exists():
                orig_path.write_text(self.original_pdb, encoding="utf-8")
            payload = {
                "session_id": self.session_id,
                "name": self.name,
                "source": self.source,
                "ph": self.ph,
                "propka_overrides": self.propka_overrides,
                "his_default": self.his_default,
                "respect_metal": self.respect_metal,
                "respect_disulfides": self.respect_disulfides,
                "history_log": self.history_log[-80:],
                "baseline": self.baseline,
                "input_md5": self.input_md5,
                "active_assembly_id": self.active_assembly_id,
                "symmetry_copy_chains": self.symmetry_copy_chains,
                "docking_prefs": self.docking_prefs,
                "rebuilt_gaps": self.rebuilt_gaps,
                "refused_gaps": self.refused_gaps,
                "rcsb_ligands": self.rcsb_ligands,
                "undo_stack": [
                    {"pdb_text": c.pdb_text, "action": c.action, "timestamp": c.timestamp}
                    for c in self.undo_stack[-15:]
                ],
                "redo_stack": [
                    {"pdb_text": c.pdb_text, "action": c.action, "timestamp": c.timestamp}
                    for c in self.redo_stack[-15:]
                ],
            }
            meta_path.write_text(json.dumps(payload), encoding="utf-8")
        except Exception:
            pass

    @classmethod
    def load_from_disk(cls, sid: str) -> Optional["StructureSession"]:
        pdb_path = _SESSIONS_DIR / f"{sid}.pdb"
        meta_path = _SESSIONS_DIR / f"{sid}.json"
        if not (pdb_path.exists() and meta_path.exists()):
            return None
        try:
            data = json.loads(meta_path.read_text(encoding="utf-8"))
            st = gemmi.read_pdb_string(pdb_path.read_text(encoding="utf-8"))
            st.setup_entities()
            st.assign_label_seq_id()
            orig_path = _SESSIONS_DIR / f"{sid}_orig.pdb"
            orig_pdb = orig_path.read_text(encoding="utf-8") if orig_path.exists() else pdb_path.read_text(encoding="utf-8")
            raw_st = gemmi.read_pdb_string(orig_pdb)
            raw_st.setup_entities()
            raw_st.assign_label_seq_id()
            dp = _default_docking_prefs()
            dp.update(data.get("docking_prefs") or {})
            sess = cls(
                session_id=data["session_id"],
                name=data["name"],
                source=data["source"],
                structure=st,
                ph=float(data.get("ph", 7.4)),
                propka_overrides=dict(data.get("propka_overrides") or {}),
                his_default=str(data.get("his_default", "auto")),
                respect_metal=bool(data.get("respect_metal", True)),
                respect_disulfides=bool(data.get("respect_disulfides", True)),
                history_log=list(data.get("history_log") or []),
                baseline=data.get("baseline"),
                original_pdb=orig_pdb,
                raw_input_structure=raw_st,
                input_md5=str(data.get("input_md5") or md5_hex(orig_pdb)),
                active_assembly_id=str(data.get("active_assembly_id", "ASU")),
                symmetry_copy_chains=list(data.get("symmetry_copy_chains") or []),
                docking_prefs=dp,
                rebuilt_gaps=list(data.get("rebuilt_gaps") or []),
                refused_gaps=dict(data.get("refused_gaps") or {}),
                rcsb_ligands=data.get("rcsb_ligands"),
                undo_stack=[Checkpoint(**c) for c in (data.get("undo_stack") or [])],
                redo_stack=[Checkpoint(**c) for c in (data.get("redo_stack") or [])],
            )
            return sess
        except Exception:
            return None

    def push_undo(self, action_label: str) -> None:
        ts = datetime.now().strftime("%H:%M:%S")
        self.undo_stack.append(Checkpoint(pdb_text=structure_to_pdb(self.structure), action=action_label, timestamp=ts))
        if len(self.undo_stack) > 40:
            self.undo_stack.pop(0)
        self.redo_stack.clear()

    def record_history(self, message: str) -> None:
        ts = datetime.now().strftime("%H:%M:%S")
        self.history_log.append(f"[{ts}] {message}")

    def undo(self) -> Optional[str]:
        if not self.undo_stack:
            return None
        cp = self.undo_stack.pop()
        ts = datetime.now().strftime("%H:%M:%S")
        self.redo_stack.append(Checkpoint(pdb_text=structure_to_pdb(self.structure), action=cp.action, timestamp=ts))
        self.structure = gemmi.read_pdb_string(cp.pdb_text)
        self.structure.setup_entities()
        self.invalidate_cache()
        self.record_history(f"Undo: {cp.action}")
        return cp.action

    def redo(self) -> Optional[str]:
        if not self.redo_stack:
            return None
        cp = self.redo_stack.pop()
        ts = datetime.now().strftime("%H:%M:%S")
        self.undo_stack.append(Checkpoint(pdb_text=structure_to_pdb(self.structure), action=cp.action, timestamp=ts))
        self.structure = gemmi.read_pdb_string(cp.pdb_text)
        self.structure.setup_entities()
        self.invalidate_cache()
        self.record_history(f"Redo: {cp.action}")
        return cp.action

    def get_full_state(self, force_refresh: bool = False) -> Dict[str, Any]:
        if self._cached_state is not None and not force_refresh:
            return self._cached_state

        metadata = extract_metadata(self.structure, source=self.source, title_fallback=self.name)
        fixer_diag = diagnose_pdbfixer(self.structure)
        quality = analyze_structure_quality(self.structure)
        propka = analyze_propka(
            self.structure,
            ph=self.ph,
            respect_metal_coordination=self.respect_metal,
            respect_disulfides=self.respect_disulfides,
            his_default=self.his_default,
            overrides=self.propka_overrides,
        )
        lig_info = analyze_ligands_and_metals(self.structure)
        pockets = detect_binding_pockets(self.structure, max_pockets=5)

        # Merge detected druggable cavities into docking_boxes so apo proteins have targeted Vina boxes
        existing_box_ids = {b["ligand_id"] for b in lig_info["docking_boxes"]}
        for pk in pockets:
            if pk["id"] not in existing_box_ids:
                lig_info["docking_boxes"].append({
                    "ligand_id": pk["id"],
                    "label": pk["label"],
                    "center": pk["center"],
                    "size": pk["size"],
                })

        readiness = compute_readiness_audit(fixer_diag, quality, propka, lig_info)
        solvation_est = estimate_solvation_system(
            self.structure,
            net_charge=propka["net_charge_discrete"],
        )

        # Compact 1-letter sequence strip per chain with QA status badges
        seq_by_chain: Dict[str, List[Dict[str, Any]]] = {}
        lig_contact_keys = {
            ct["key"] for lig in lig_info["ligands"] for ct in lig.get("contacts", [])
        }
        metal_coord_keys = {
            co["key"] for m in lig_info["metals"] for co in m.get("coordinators", [])
        }
        for r in quality["residues"]:
            aa1 = AA_3_TO_1.get(r["std_resname"], "X")
            sev = "clean"
            if any(x in r["issues"] for x in ("clash", "rama_outlier", "cis_nonpro", "chiral_outlier")):
                sev = "error"
            elif r["issues"]:
                sev = "warn"
            elif r["key"] in metal_coord_keys:
                sev = "metal"
            elif r["key"] in lig_contact_keys:
                sev = "pocket"
            seq_by_chain.setdefault(r["chain"], []).append({
                "chain": r["chain"],
                "resseq": r["resseq"],
                "icode": r["icode"],
                "resname": r["resname"],
                "aa1": aa1,
                "key": r["key"],
                "severity": sev,
                "issues": r["issues"],
            })

        # 3D structural overlays: disulfide cylinders, metal coordination lines, ligand H-bonds
        disulfide_bonds = []
        if len(self.structure) > 0:
            cys_sgs = []
            for ch in self.structure[0]:
                for res in ch:
                    if canonical_resname(res.name) == "CYS":
                        sg = _get_pos(res, "SG")
                        if sg is not None:
                            cys_sgs.append((f"{ch.name}:{res.seqid.num}", sg))
            for i in range(len(cys_sgs)):
                for j in range(i + 1, len(cys_sgs)):
                    d = distance(cys_sgs[i][1], cys_sgs[j][1])
                    if d <= 2.6:
                        disulfide_bonds.append({
                            "label": f"{cys_sgs[i][0]}–{cys_sgs[j][0]} ({d:.2f} Å)",
                            "p1": [round(float(x), 3) for x in cys_sgs[i][1]],
                            "p2": [round(float(x), 3) for x in cys_sgs[j][1]],
                        })

        metal_bonds = []
        if len(self.structure) > 0:
            res_lookup = {}
            for ch in self.structure[0]:
                for res in ch:
                    res_lookup[f"{ch.name}:{res.seqid.num}{res.seqid.icode.strip()}"] = res
            for m in lig_info["metals"]:
                for co in m["coordinators"]:
                    r_obj = res_lookup.get(co["key"])
                    if r_obj is not None:
                        ap = _get_pos(r_obj, co["atom"])
                        if ap is not None:
                            metal_bonds.append({
                                "label": f"{m['resname']}–{co['resname']}{co['resseq']}.{co['atom']} ({co['distance_A']} Å)",
                                "p1": m["center"],
                                "p2": [round(float(x), 3) for x in ap],
                            })

        ligand_hbonds = []
        if len(self.structure) > 0:
            res_lookup = {}
            for ch in self.structure[0]:
                for res in ch:
                    res_lookup[f"{ch.name}:{res.seqid.num}{res.seqid.icode.strip()}"] = res
            for lig in lig_info["ligands"]:
                lig_key = f"{lig['chain']}:{lig['resseq']}{lig['icode']}"
                lig_res = res_lookup.get(lig_key)
                if lig_res is None:
                    continue
                for ct in lig.get("contacts", []):
                    if ct["type"] != "hbond":
                        continue
                    pr_res = res_lookup.get(ct["key"])
                    if pr_res is None:
                        continue
                    p_lig = _get_pos(lig_res, ct["ligand_atom"])
                    p_pro = _get_pos(pr_res, ct["protein_atom"])
                    if p_lig is not None and p_pro is not None:
                        ligand_hbonds.append({
                            "label": f"{lig['resname']}.{ct['ligand_atom']}…{ct['resname']}{ct['resseq']}.{ct['protein_atom']} ({ct['distance_A']} Å)",
                            "p1": [round(float(x), 3) for x in p_lig],
                            "p2": [round(float(x), 3) for x in p_pro],
                        })

        # Capture immutable baseline metrics on initial load for Before -> After comparison
        issues_count = sum(1 for r in quality["residues"] if r["issues"])
        if self.baseline is None:
            self.baseline = {
                "overall_score": readiness["overall_score"],
                "md_score": readiness["md"]["score"],
                "docking_score": readiness["docking"]["score"],
                "quality_score": readiness["crystallographic"]["score"],
                "clashscore": quality["clashes"]["clashscore"],
                "issues_count": issues_count,
                "atoms": metadata["counts"]["atoms"],
                "hydrogens": metadata["counts"]["hydrogens"],
                "rama_favored_pct": quality["ramachandran"]["favored_pct"],
                "rotamer_outliers": quality["rotamers"]["counts"]["outlier"],
            }

        # Assemblies, REMARK 465 missing regions, and ADapt unified preparation tree
        asm_source = self.raw_input_structure if self.raw_input_structure is not None else self.structure
        assemblies = list_biological_assemblies(asm_source, active_assembly_id=self.active_assembly_id)
        missing_regions = detect_missing_regions(
            self.structure,
            rebuilt_records=self.rebuilt_gaps,
            refused_records=self.refused_gaps,
        )
        prep_tree = build_preparation_tree(
            self.structure,
            docking_prefs=self.docking_prefs,
            missing_regions=missing_regions,
            rcsb_ligands=self.rcsb_ligands,
            symmetry_copy_chains=self.symmetry_copy_chains,
        )
        rebuilt_keys = [
            k
            for rg in self.rebuilt_gaps
            for k in (rg.get("rebuilt_keys") or [])
        ]

        state = {
            "session_id": self.session_id,
            "name": self.name,
            "source": self.source,
            "pdb_text": structure_to_pdb(self.structure),
            "input_md5": self.input_md5,
            "can_undo": len(self.undo_stack) > 0,
            "can_redo": len(self.redo_stack) > 0,
            "undo_count": len(self.undo_stack),
            "redo_count": len(self.redo_stack),
            "history": list(self.history_log[-80:]),
            "baseline": self.baseline,
            "metadata": metadata,
            "pdbfixer": fixer_diag,
            "validation": quality,
            "propka": propka,
            "ligands_info": lig_info,
            "pockets": pockets,
            "sequence_chains": seq_by_chain,
            "assemblies": assemblies,
            "active_assembly_id": self.active_assembly_id,
            "symmetry_copy_chains": list(self.symmetry_copy_chains),
            "missing_regions": missing_regions,
            "prep_tree": prep_tree,
            "docking_prefs": self.docking_prefs,
            "rebuilt_keys": rebuilt_keys,
            "overlays": {
                "disulfides": disulfide_bonds,
                "metal_bonds": metal_bonds,
                "ligand_hbonds": ligand_hbonds,
            },
            "solvation_estimate": solvation_est,
            "readiness": readiness,
        }
        self._cached_state = state
        return state

    def get_remark_850_lines(self) -> List[str]:
        st = self.get_full_state()
        return build_remark_850_lines(
            session_name=self.name,
            source=self.source,
            input_md5=self.input_md5,
            active_assembly_id=self.active_assembly_id,
            ph=self.ph,
            docking_prefs=self.docking_prefs,
            rebuilt_gaps=self.rebuilt_gaps,
            propka_summary=st.get("propka"),
        )

    def get_provenance_record(self, output_files: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
        st = self.get_full_state()
        return build_provenance_record(
            session_id=self.session_id,
            name=self.name,
            source=self.source,
            input_md5=self.input_md5,
            current_pdb_text=st["pdb_text"],
            active_assembly_id=self.active_assembly_id,
            symmetry_copy_chains=self.symmetry_copy_chains,
            ph=self.ph,
            propka_data=st["propka"],
            docking_prefs=self.docking_prefs,
            rebuilt_gaps=self.rebuilt_gaps,
            history_log=list(self.history_log),
            output_files=output_files,
        )


class SessionStore:
    def __init__(self) -> None:
        self._sessions: Dict[str, StructureSession] = {}

    def create(self, name: str, source: str, st: gemmi.Structure) -> StructureSession:
        sid = uuid.uuid4().hex[:10]
        orig_pdb = structure_to_pdb(st)
        raw_clone = st.clone()
        # Fast non-blocking check if RCSB ligand cache already exists for this PDB ID
        rcsb_cached = None
        base_id = name.split(".")[0].upper()
        if len(base_id) == 4 and base_id.isalnum() and source.startswith(("RCSB", "Bundled")):
            try:
                rcsb_cached = lookup_rcsb_entry_ligands(base_id, force_refresh=False)
            except Exception:
                rcsb_cached = None

        sess = StructureSession(
            session_id=sid,
            name=name,
            source=source,
            structure=st,
            original_pdb=orig_pdb,
            raw_input_structure=raw_clone,
            input_md5=md5_hex(orig_pdb),
            rcsb_ligands=rcsb_cached,
        )
        sess.record_history(f"Loaded structure '{name}' ({source}) [md5:{sess.input_md5[:10]}]")
        self._sessions[sid] = sess
        sess.save_to_disk()
        return sess

    def get(self, sid: str) -> StructureSession:
        if sid not in self._sessions:
            loaded = StructureSession.load_from_disk(sid)
            if loaded is not None:
                self._sessions[sid] = loaded
                return loaded
            raise KeyError(f"Session '{sid}' not found")
        return self._sessions[sid]

    def list_sessions(self) -> List[Dict[str, Any]]:
        return [
            {
                "session_id": s.session_id,
                "name": s.name,
                "source": s.source,
                "undo_count": len(s.undo_stack),
            }
            for s in self._sessions.values()
        ]


STORE = SessionStore()
