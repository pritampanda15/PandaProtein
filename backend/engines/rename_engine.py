"""
Structure, Chain, Residue, Ligand & File Renaming Engine (Gemmi-powered):
- Rename structure filename & PDB TITLE header
- Rename single or batch chain IDs (collision-safe two-pass swap, e.g. A<->B)
- Merge chains (with optional sequential renumbering)
- Reorder chains, delete chains, auto-assign sequential chain IDs (A, B, C...)
- Per-chain residue renumbering (custom start number or +/- offset)
- Rename ligand / heterogen 3-letter codes and move heterogens between chains
- Batch residue-name renaming (e.g. HIS->HIE, or canonicalize HID/HIE/HIP/CYX/ASH/GLH -> standard PDB names)
- Normalize occupancies (1.0), B-factors, and segment IDs
"""

from __future__ import annotations

import string
from typing import Any, Dict, List, Optional

import gemmi

from .chemistry_data import ALL_PROTEIN_RESNAMES, PROTONATION_SYNONYMS, WATER_RESNAMES
from .geometry_utils import canonical_resname


def rename_chains_batch(st: gemmi.Structure, mapping: Dict[str, str]) -> Dict[str, Any]:
    """
    Collision-safe batch chain renaming (supports swaps like {'A': 'B', 'B': 'A'}).
    Updates both model chains and SSBOND/Connection partner chain IDs.
    """
    if len(st) == 0:
        raise ValueError("Empty structure")

    cleaned: Dict[str, str] = {}
    for old_id, new_id in mapping.items():
        o = old_id.strip()
        n = new_id.strip()
        if not o or not n:
            continue
        if len(n) > 4:
            raise ValueError(f"Chain ID '{n}' is too long (max 4 chars; 1 char recommended for PDB)")
        if o != n:
            cleaned[o] = n

    if not cleaned:
        return {"changed": False, "message": "No chain ID changes specified"}

    # Validate no two chains collapse onto the same target name unintentionally
    model = st[0]
    existing_names = [c.name for c in model]
    final_names = [cleaned.get(nm, nm) for nm in existing_names]
    if len(set(final_names)) != len(final_names):
        raise ValueError(
            "Duplicate target chain ID detected — use 'Merge Chains' if you want to combine two chains."
        )

    # Two-pass rename to avoid collisions during swaps (e.g., A->B, B->A)
    temp_map = {old: f"__TMP_{idx}__" for idx, old in enumerate(cleaned.keys())}
    for m in st:
        for ch in m:
            if ch.name in temp_map:
                ch.name = temp_map[ch.name]
    inv_temp = {tmp: cleaned[old] for old, tmp in temp_map.items()}
    for m in st:
        for ch in m:
            if ch.name in inv_temp:
                new_nm = inv_temp[ch.name]
                ch.name = new_nm
                for res in ch:
                    res.subchain = new_nm

    # Update connection / SSBOND chain references
    for conn in st.connections:
        if conn.partner1.chain_name in cleaned:
            conn.partner1.chain_name = cleaned[conn.partner1.chain_name]
        if conn.partner2.chain_name in cleaned:
            conn.partner2.chain_name = cleaned[conn.partner2.chain_name]

    st.setup_entities()
    summary = ", ".join(f"{k}→{v}" for k, v in cleaned.items())
    return {"changed": True, "mapping": cleaned, "message": f"Renamed chain(s): {summary}"}


def auto_standardize_chain_ids(st: gemmi.Structure) -> Dict[str, Any]:
    """Assign clean sequential single-letter chain IDs A, B, C, ... in current order."""
    if len(st) == 0:
        raise ValueError("Empty structure")
    alphabet = list(string.ascii_uppercase) + list(string.digits)
    mapping = {}
    for idx, ch in enumerate(st[0]):
        target = alphabet[idx] if idx < len(alphabet) else f"C{idx}"
        if ch.name != target:
            mapping[ch.name] = target
    if not mapping:
        return {"changed": False, "message": "Chain IDs are already standardized (A, B, C…)"}
    return rename_chains_batch(st, mapping)


def reorder_chains(st: gemmi.Structure, chain_order: List[str]) -> Dict[str, Any]:
    """Reorder chains in Model 1 according to `chain_order`."""
    if len(st) == 0:
        raise ValueError("Empty structure")
    model = st[0]
    by_name = {ch.name: ch.clone() for ch in model}
    ordered = [c for c in chain_order if c in by_name]
    for ch in model:
        if ch.name not in ordered:
            ordered.append(ch.name)

    while len(model) > 0:
        del model[0]
    for cid in ordered:
        model.add_chain(by_name[cid])
    st.setup_entities()
    return {"changed": True, "order": ordered, "message": f"Reordered chains → {' → '.join(ordered)}"}


def merge_chains(
    st: gemmi.Structure,
    source_chain: str,
    target_chain: str,
    renumber_merged: bool = True,
) -> Dict[str, Any]:
    """Merge all residues from `source_chain` into `target_chain` and remove `source_chain`."""
    if len(st) == 0:
        raise ValueError("Empty structure")
    if source_chain == target_chain:
        raise ValueError("Source and target chains must be different")
    model = st[0]
    src_idx = next((i for i, c in enumerate(model) if c.name == source_chain), None)
    tgt = next((c for c in model if c.name == target_chain), None)
    if src_idx is None or tgt is None:
        raise ValueError(f"Chain '{source_chain}' or '{target_chain}' not found")

    src = model[src_idx]
    max_num = max((r.seqid.num for r in tgt), default=0)
    next_num = max_num + 1
    moved = 0
    for res in src:
        cloned = res.clone()
        cloned.subchain = target_chain
        if renumber_merged:
            cloned.seqid = gemmi.SeqId(next_num, " ")
            next_num += 1
        tgt.add_residue(cloned)
        moved += 1

    del model[src_idx]
    st.setup_entities()
    return {
        "changed": True,
        "message": f"Merged chain {source_chain} ({moved} residues) into chain {target_chain}",
    }


def delete_chains_batch(st: gemmi.Structure, chain_ids: List[str]) -> Dict[str, Any]:
    """Delete one or more chains across all models and clean up dangling SSBOND/LINK connections."""
    if len(st) == 0:
        raise ValueError("Empty structure")
    targets = {c.strip() for c in chain_ids if c and c.strip()}
    if not targets:
        raise ValueError("No chain IDs specified for deletion")

    model0 = st[0]
    existing = [c.name for c in model0]
    found = [c for c in existing if c in targets]
    if not found:
        raise ValueError(f"Chain(s) {', '.join(sorted(targets))} not found")
    if len(found) >= len(existing):
        raise ValueError("Cannot delete all chains — at least one chain must remain in the structure")

    for m in st:
        for idx in range(len(m) - 1, -1, -1):
            if m[idx].name in targets:
                del m[idx]

    # Clean up dangling SSBOND / LINK connections pointing to deleted chains
    for ci in range(len(st.connections) - 1, -1, -1):
        conn = st.connections[ci]
        if conn.partner1.chain_name in targets or conn.partner2.chain_name in targets:
            del st.connections[ci]

    st.setup_entities()
    return {
        "changed": True,
        "deleted": found,
        "message": f"Deleted chain(s): {', '.join(found)}",
    }


def delete_chain(st: gemmi.Structure, chain_id: str) -> Dict[str, Any]:
    """Delete a single chain from the structure."""
    return delete_chains_batch(st, [chain_id])


def renumber_chain(
    st: gemmi.Structure,
    chain_id: str,
    mode: str = "start",
    start_num: int = 1,
    offset: int = 0,
    include_heterogens: bool = False,
    clear_icodes: bool = True,
) -> Dict[str, Any]:
    """
    Renumber residues in a chain:
      - mode='start': sequential numbering starting at `start_num`
      - mode='offset': shift existing residue numbers by `offset` (e.g. +10 or -5)
    """
    if len(st) == 0:
        raise ValueError("Empty structure")
    chain = next((c for c in st[0] if c.name == chain_id), None)
    if chain is None:
        raise ValueError(f"Chain '{chain_id}' not found")

    count = 0
    cur = start_num
    het_cur = max(start_num + 500, 501)
    for res in chain:
        rname = res.name.strip().upper()
        is_prot = rname in ALL_PROTEIN_RESNAMES
        if not is_prot and not include_heterogens:
            continue
        ic = " " if clear_icodes else res.seqid.icode
        if mode == "offset":
            new_n = max(1, res.seqid.num + offset)
            res.seqid = gemmi.SeqId(new_n, ic)
        else:
            if is_prot:
                res.seqid = gemmi.SeqId(cur, ic)
                cur += 1
            else:
                res.seqid = gemmi.SeqId(het_cur, ic)
                het_cur += 1
        count += 1

    st.setup_entities()
    desc = f"shifted by {offset:+d}" if mode == "offset" else f"sequentially from {start_num}"
    return {
        "changed": True,
        "message": f"Renumbered {count} residue(s) in chain {chain_id} ({desc})",
    }


def rename_residue_code(
    st: gemmi.Structure,
    old_resname: str,
    new_resname: str,
    chain_id: Optional[str] = None,
    resseq: Optional[int] = None,
    move_to_chain: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Rename a specific residue/ligand/heterogen 3-letter code (e.g. BEN -> LIG, UNL -> INH, HIS -> HIE),
    optionally moving it to another chain (creates target chain if it doesn't exist).
    """
    if len(st) == 0:
        raise ValueError("Empty structure")
    old_r = old_resname.strip().upper()
    new_r = new_resname.strip().upper()
    if not (1 <= len(new_r) <= 4):
        raise ValueError("New residue/ligand code must be 1–4 characters (3 standard for PDB)")

    model = st[0]
    renamed_count = 0
    to_move: List[tuple] = []

    for ch in model:
        if chain_id and ch.name != chain_id:
            continue
        for ri, res in enumerate(ch):
            if res.name.strip().upper() != old_r:
                continue
            if resseq is not None and res.seqid.num != resseq:
                continue
            res.name = new_r
            renamed_count += 1
            if move_to_chain and move_to_chain.strip() and move_to_chain.strip() != ch.name:
                to_move.append((ch.name, ri, res.clone()))

    if renamed_count == 0:
        raise ValueError(f"No residue matching '{old_r}' found")

    if to_move and move_to_chain:
        tgt_name = move_to_chain.strip()
        tgt_chain = next((c for c in model if c.name == tgt_name), None)
        if tgt_chain is None:
            tgt_chain = gemmi.Chain(tgt_name)
            model.add_chain(tgt_chain)
            tgt_chain = next(c for c in model if c.name == tgt_name)
        # Remove from source chains back-to-front
        for src_name, ri, cloned in reversed(to_move):
            src_ch = next((c for c in model if c.name == src_name), None)
            if src_ch is not None and ri < len(src_ch):
                del src_ch[ri]
            cloned.subchain = tgt_name
            tgt_chain.add_residue(cloned)

    st.setup_entities()
    move_msg = f" and moved to chain {move_to_chain}" if (to_move and move_to_chain) else ""
    return {
        "changed": True,
        "count": renamed_count,
        "message": f"Renamed {renamed_count} occurrence(s) of {old_r} → {new_r}{move_msg}",
    }


def canonicalize_residue_names(st: gemmi.Structure) -> Dict[str, Any]:
    """Reset AMBER/GROMACS protonation aliases (HID/HIE/HIP/CYX/CYM/ASH/GLH/LYN) back to standard PDB 20-AA names."""
    if len(st) == 0:
        raise ValueError("Empty structure")
    changed = 0
    for ch in st[0]:
        for res in ch:
            rn = res.name.strip().upper()
            if rn in PROTONATION_SYNONYMS:
                res.name = PROTONATION_SYNONYMS[rn]
                changed += 1
    st.setup_entities()
    return {
        "changed": changed > 0,
        "message": f"Canonicalized {changed} force-field residue name(s) back to standard PDB names"
        if changed
        else "All residue names are already canonical PDB names",
    }


def normalize_structure_fields(
    st: gemmi.Structure,
    reset_occupancy: bool = False,
    reset_bfactor: Optional[float] = None,
    strip_waters_to_chain: Optional[str] = None,
) -> Dict[str, Any]:
    """Normalize occupancies to 1.00, optionally set uniform B-factor, or consolidate waters to a single chain."""
    if len(st) == 0:
        raise ValueError("Empty structure")
    actions = []
    model = st[0]

    if reset_occupancy:
        for ch in model:
            for res in ch:
                for a in res:
                    a.occ = 1.0
        actions.append("Reset all atomic occupancies to 1.00")

    if reset_bfactor is not None:
        val = float(reset_bfactor)
        for ch in model:
            for res in ch:
                for a in res:
                    a.b_iso = round(val, 2)
        actions.append(f"Set uniform B-factor to {val:.2f} Å²")

    if strip_waters_to_chain:
        w_chain_id = strip_waters_to_chain.strip()
        waters = []
        for ch in model:
            for ri in range(len(ch) - 1, -1, -1):
                if ch[ri].name.strip().upper() in WATER_RESNAMES:
                    waters.append(ch[ri].clone())
                    del ch[ri]
        if waters:
            waters.reverse()
            tgt = next((c for c in model if c.name == w_chain_id), None)
            if tgt is None:
                model.add_chain(gemmi.Chain(w_chain_id))
                tgt = next(c for c in model if c.name == w_chain_id)
            for wi, w in enumerate(waters, start=1):
                w.subchain = w_chain_id
                w.seqid = gemmi.SeqId(wi, " ")
                tgt.add_residue(w)
            actions.append(f"Consolidated {len(waters)} water(s) into chain {w_chain_id} (1..{len(waters)})")

    if not actions:
        return {"changed": False, "message": "No normalization options selected"}
    st.setup_entities()
    return {"changed": True, "message": "; ".join(actions)}
