"""
Authoritative structural biology & crystallographic reference data:
- Standard amino acid heavy-atom sets & internal-coordinate templates
- Engh & Huber (1991/2001) ideal backbone bond lengths & angles (with sigma)
- Dunbrack/Richardson rotamer library (chi1-chi4 canonical states + probabilities)
- MolProbity / Bondi van der Waals radii
- Non-standard residue conversions, cofactors, metals, crystallization buffers
- Protonation state nomenclature (AMBER/GROMACS/CHARMM)
"""

from __future__ import annotations

from typing import Dict, List, Set, Tuple

STANDARD_AA: Set[str] = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}

PROTONATION_SYNONYMS: Dict[str, str] = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "ASH": "ASP", "ASPH": "ASP",
    "GLH": "GLU", "GLUH": "GLU",
    "CYX": "CYS", "CYM": "CYS", "CYS2": "CYS", "CYSH": "CYS",
    "LYN": "LYS", "LSN": "LYS",
    "TYM": "TYR",
    "ARN": "ARG",
}

ALL_PROTEIN_RESNAMES: Set[str] = STANDARD_AA | set(PROTONATION_SYNONYMS.keys()) | {"ACE", "NME", "NMA", "NH2"}

AA_3_TO_1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

NONSTANDARD_TO_STANDARD: Dict[str, str] = {
    "MSE": "MET",  # Selenomethionine
    "SEP": "SER",  # Phosphoserine
    "TPO": "THR",  # Phosphothreonine
    "PTR": "TYR",  # Phosphotyrosine
    "CSO": "CYS",  # S-hydroxycysteine
    "CSD": "CYS",  # 3-sulfinoalanine
    "CME": "CYS",  # S,S-(2-hydroxyethyl)thiocysteine
    "OCS": "CYS",  # Cysteinesulfonic acid
    "CSS": "CYS",  # S-mercaptocysteine
    "KCX": "LYS",  # Lysine NZ-carboxylic acid
    "MLY": "LYS",  # N-dimethyl-lysine
    "M3L": "LYS",  # N-trimethyllysine
    "ALY": "LYS",  # N(6)-acetyllysine
    "HYP": "PRO",  # 4-hydroxyproline
    "PCA": "GLU",  # Pyroglutamic acid
    "CGU": "GLU",  # Gamma-carboxyglutamic acid
    "TYS": "TYR",  # O-sulfo-L-tyrosine
    "LLP": "LYS",  # PLP-lysine
    "DAL": "ALA", "DVA": "VAL", "DLE": "LEU", "DPR": "PRO",
    "MEN": "ASN", "SMC": "CYS", "SCH": "CYS", "SEC": "CYS",
}

WATER_RESNAMES: Set[str] = {"HOH", "WAT", "DOD", "H2O", "SOL", "TIP", "TIP3", "SPC"}

METAL_IONS: Set[str] = {
    "ZN", "MG", "CA", "FE", "FE2", "MN", "CU", "CU1", "CO", "NI",
    "NA", "K", "LI", "CD", "HG", "PB", "SR", "BA", "CS", "RB", "MO", "W",
}

HALIDE_IONS: Set[str] = {"CL", "BR", "I", "F"}

CRYSTALLIZATION_ARTIFACTS: Set[str] = {
    "EDO", "GOL", "PEG", "PGE", "1PE", "P6G", "PG4", "PE4", "PE8",
    "SO4", "PO4", "SUL", "ACT", "ACY", "FMT", "DMS", "DMSO",
    "MPD", "TRS", "MES", "EPE", "HEP", "CIT", "FLC", "TAR", "TLA",
    "BME", "DTT", "IPA", "EOH", "MOH", "ACN", "ACE", "NH4", "NO3",
    "SCN", "AZI", "CAC", "MLI", "IMD", "BCT", "UNX", "UNL", "BOG",
    "LDA", "LMT", "OLC", "OLA", "CLR", "POP", "MYR",
}

COFACTORS: Set[str] = {
    "NAD", "NAI", "NAP", "NDP", "FAD", "FMN", "HEM", "HEC", "HEA", "HEB",
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "ANP", "GNP", "AGS", "ACP",
    "PLP", "PMP", "SAM", "SAH", "COA", "ACO", "TPP", "THP", "BTN", "B12",
    "CNC", "SF4", "FES", "F3S", "MQ7", "UQ1", "UQ2", "PQQ", "MGD", "GSH",
}

# Expected heavy atoms per standard amino acid (excluding OXT)
EXPECTED_HEAVY_ATOMS: Dict[str, List[str]] = {
    "ALA": ["N", "CA", "C", "O", "CB"],
    "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
    "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
    "CYS": ["N", "CA", "C", "O", "CB", "SG"],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "GLY": ["N", "CA", "C", "O"],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"],
    "SER": ["N", "CA", "C", "O", "CB", "OG"],
    "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2"],
    "TRP": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
}

# Internal coordinate definitions for sidechain heavy atoms:
# atom_name: (parent1, parent2, parent3, bond_length_A, bond_angle_deg, default_dihedral_deg)
# Placed using NeRF (Natural Extension Reference Frame) from (parent3, parent2, parent1).
SIDECHAIN_IC: Dict[str, List[Tuple[str, str, str, str, float, float, float]]] = {
    "ALA": [
        ("CB", "CA", "N", "C", 1.521, 110.4, -122.6),
    ],
    "GLY": [],
    "SER": [
        ("CB", "CA", "N", "C", 1.530, 110.4, -122.6),
        ("OG", "CB", "CA", "N", 1.417, 110.8, 62.0),
    ],
    "CYS": [
        ("CB", "CA", "N", "C", 1.530, 110.8, -122.6),
        ("SG", "CB", "CA", "N", 1.808, 114.4, -65.0),
    ],
    "VAL": [
        ("CB", "CA", "N", "C", 1.540, 111.5, -122.6),
        ("CG1", "CB", "CA", "N", 1.521, 110.5, 177.0),
        ("CG2", "CB", "CA", "N", 1.521, 110.5, -63.0),
    ],
    "THR": [
        ("CB", "CA", "N", "C", 1.540, 111.5, -122.6),
        ("OG1", "CB", "CA", "N", 1.433, 109.6, 60.0),
        ("CG2", "CB", "CA", "N", 1.521, 110.5, -60.0),
    ],
    "ILE": [
        ("CB", "CA", "N", "C", 1.540, 111.5, -122.6),
        ("CG1", "CB", "CA", "N", 1.530, 110.4, -60.0),
        ("CG2", "CB", "CA", "N", 1.521, 110.5, 62.0),
        ("CD1", "CG1", "CB", "CA", 1.513, 113.8, 170.0),
    ],
    "LEU": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.530, 116.3, -65.0),
        ("CD1", "CG", "CB", "CA", 1.521, 110.7, 175.0),
        ("CD2", "CG", "CB", "CA", 1.521, 110.7, -65.0),
    ],
    "ASP": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.516, 112.6, -70.0),
        ("OD1", "CG", "CB", "CA", 1.249, 118.4, -15.0),
        ("OD2", "CG", "CB", "CA", 1.249, 118.4, 165.0),
    ],
    "ASN": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.516, 112.6, -70.0),
        ("OD1", "CG", "CB", "CA", 1.231, 120.8, -20.0),
        ("ND2", "CG", "CB", "CA", 1.328, 116.4, 160.0),
    ],
    "GLU": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.520, 114.1, -67.0),
        ("CD", "CG", "CB", "CA", 1.516, 112.6, 177.0),
        ("OE1", "CD", "CG", "CB", 1.249, 118.4, 0.0),
        ("OE2", "CD", "CG", "CB", 1.249, 118.4, 180.0),
    ],
    "GLN": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.520, 114.1, -67.0),
        ("CD", "CG", "CB", "CA", 1.516, 112.6, 177.0),
        ("OE1", "CD", "CG", "CB", 1.231, 120.8, 0.0),
        ("NE2", "CD", "CG", "CB", 1.328, 116.4, 180.0),
    ],
    "MET": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.520, 113.7, -67.0),
        ("SD", "CG", "CB", "CA", 1.803, 112.7, -65.0),
        ("CE", "SD", "CG", "CB", 1.791, 100.9, -70.0),
    ],
    "LYS": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.520, 114.1, -67.0),
        ("CD", "CG", "CB", "CA", 1.520, 111.3, 180.0),
        ("CE", "CD", "CG", "CB", 1.520, 111.3, 180.0),
        ("NZ", "CE", "CD", "CG", 1.489, 111.9, 180.0),
    ],
    "ARG": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.520, 114.1, -67.0),
        ("CD", "CG", "CB", "CA", 1.520, 111.3, 180.0),
        ("NE", "CD", "CG", "CB", 1.460, 112.0, 180.0),
        ("CZ", "NE", "CD", "CG", 1.329, 124.2, 180.0),
        ("NH1", "CZ", "NE", "CD", 1.326, 120.0, 0.0),
        ("NH2", "CZ", "NE", "CD", 1.326, 120.0, 180.0),
    ],
    "HIS": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.497, 113.8, -65.0),
        ("ND1", "CG", "CB", "CA", 1.378, 122.7, -75.0),
        ("CD2", "CG", "CB", "CA", 1.354, 131.2, 105.0),
        ("CE1", "ND1", "CG", "CB", 1.320, 107.2, 180.0),
        ("NE2", "CD2", "CG", "CB", 1.374, 107.2, 180.0),
    ],
    "PHE": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.502, 113.8, -65.0),
        ("CD1", "CG", "CB", "CA", 1.384, 120.7, 90.0),
        ("CD2", "CG", "CB", "CA", 1.384, 120.7, -90.0),
        ("CE1", "CD1", "CG", "CB", 1.382, 120.7, 180.0),
        ("CE2", "CD2", "CG", "CB", 1.382, 120.7, 180.0),
        ("CZ", "CE1", "CD1", "CG", 1.382, 120.0, 0.0),
    ],
    "TYR": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.512, 113.9, -65.0),
        ("CD1", "CG", "CB", "CA", 1.389, 120.8, 90.0),
        ("CD2", "CG", "CB", "CA", 1.389, 120.8, -90.0),
        ("CE1", "CD1", "CG", "CB", 1.382, 121.0, 180.0),
        ("CE2", "CD2", "CG", "CB", 1.382, 121.0, 180.0),
        ("CZ", "CE1", "CD1", "CG", 1.378, 119.6, 0.0),
        ("OH", "CZ", "CE1", "CD1", 1.376, 119.7, 180.0),
    ],
    "TRP": [
        ("CB", "CA", "N", "C", 1.530, 110.5, -122.6),
        ("CG", "CB", "CA", "N", 1.498, 114.1, -65.0),
        ("CD1", "CG", "CB", "CA", 1.365, 126.9, 95.0),
        ("CD2", "CG", "CB", "CA", 1.433, 126.8, -85.0),
        ("NE1", "CD1", "CG", "CB", 1.374, 110.2, 180.0),
        ("CE2", "CD2", "CG", "CB", 1.409, 107.2, 180.0),
        ("CE3", "CD2", "CG", "CB", 1.398, 133.9, 0.0),
        ("CZ2", "CE2", "CD2", "CG", 1.394, 122.4, 180.0),
        ("CZ3", "CE3", "CD2", "CG", 1.382, 118.8, 180.0),
        ("CH2", "CZ2", "CE2", "CD2", 1.368, 117.5, 0.0),
    ],
    "PRO": [
        ("CB", "CA", "N", "C", 1.530, 103.2, -118.0),
        ("CG", "CB", "CA", "N", 1.492, 104.5, 25.0),
        ("CD", "CG", "CB", "CA", 1.503, 105.0, -35.0),
    ],
}

# Chi angle atom quadruplets per amino acid
CHI_ATOMS: Dict[str, List[Tuple[str, str, str, str]]] = {
    "ARG": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "NE"), ("CG", "CD", "NE", "CZ")],
    "ASN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "ASP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "CYS": [("N", "CA", "CB", "SG")],
    "GLN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
    "GLU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
    "HIS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1")],
    "ILE": [("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD1")],
    "LEU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "LYS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "CE"), ("CG", "CD", "CE", "NZ")],
    "MET": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "SD"), ("CB", "CG", "SD", "CE")],
    "PHE": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "PRO": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")],
    "SER": [("N", "CA", "CB", "OG")],
    "THR": [("N", "CA", "CB", "OG1")],
    "TRP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "TYR": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "VAL": [("N", "CA", "CB", "CG1")],
}

# Richardson / Dunbrack Top-Ranked Rotamer Library:
# resname -> list of (rotamer_name, probability_pct, [chi1, chi2, ...])
ROTAMER_LIBRARY: Dict[str, List[Tuple[str, float, List[float]]]] = {
    "SER": [("p", 48.0, [62.0]), ("m", 29.0, [-65.0]), ("t", 22.0, [180.0])],
    "CYS": [("m", 56.0, [-65.0]), ("t", 26.0, [-177.0]), ("p", 17.0, [62.0])],
    "VAL": [("t", 73.0, [175.0]), ("m", 20.0, [-60.0]), ("p", 6.0, [63.0])],
    "THR": [("p", 49.0, [60.0]), ("m", 43.0, [-60.0]), ("t", 7.0, [-175.0])],
    "ILE": [
        ("mt", 59.0, [-63.0, 169.0]),
        ("mm", 16.0, [-60.0, -61.0]),
        ("pt", 12.0, [62.0, 170.0]),
        ("tt", 7.0, [168.0, 165.0]),
        ("tp", 4.0, [170.0, 66.0]),
    ],
    "LEU": [
        ("mt", 59.0, [-65.0, 175.0]),
        ("tp", 29.0, [-177.0, 63.0]),
        ("tt", 4.0, [-172.0, 155.0]),
        ("mp", 3.0, [-78.0, 66.0]),
        ("pt", 2.0, [62.0, 170.0]),
    ],
    "ASP": [
        ("m-20", 48.0, [-70.0, -15.0]),
        ("t0", 24.0, [-174.0, 0.0]),
        ("p-10", 16.0, [62.0, -10.0]),
        ("t70", 7.0, [-175.0, 65.0]),
        ("p30", 4.0, [62.0, 30.0]),
    ],
    "ASN": [
        ("m-20", 34.0, [-70.0, -20.0]),
        ("m-80", 16.0, [-70.0, -75.0]),
        ("t-20", 18.0, [-174.0, -20.0]),
        ("t30", 12.0, [-174.0, 30.0]),
        ("p-10", 11.0, [62.0, -10.0]),
        ("m120", 6.0, [-65.0, 120.0]),
    ],
    "GLU": [
        ("mt-10", 33.0, [-67.0, 177.0, -10.0]),
        ("tt0", 24.0, [-177.0, 178.0, 0.0]),
        ("mm-40", 14.0, [-65.0, -65.0, -40.0]),
        ("tp10", 10.0, [-177.0, 65.0, 10.0]),
        ("pt-20", 8.0, [62.0, 180.0, -20.0]),
        ("mp0", 6.0, [-65.0, 80.0, 0.0]),
    ],
    "GLN": [
        ("mt-30", 35.0, [-67.0, 177.0, -25.0]),
        ("tt0", 20.0, [-177.0, 177.0, 0.0]),
        ("mm-40", 14.0, [-66.0, -65.0, -40.0]),
        ("tp-100", 11.0, [-177.0, 65.0, -100.0]),
        ("pt20", 9.0, [64.0, 180.0, 20.0]),
        ("mp0", 6.0, [-67.0, 75.0, 0.0]),
    ],
    "MET": [
        ("mmm", 23.0, [-67.0, -65.0, -70.0]),
        ("mtt", 17.0, [-67.0, 178.0, 175.0]),
        ("mtm", 15.0, [-67.0, 178.0, -75.0]),
        ("tpt", 14.0, [-177.0, 65.0, 175.0]),
        ("ttt", 12.0, [-177.0, 178.0, 175.0]),
        ("mtp", 10.0, [-67.0, 178.0, 75.0]),
    ],
    "LYS": [
        ("mttt", 27.0, [-67.0, 180.0, 180.0, 180.0]),
        ("tttt", 21.0, [-177.0, 180.0, 180.0, 180.0]),
        ("mmtt", 12.0, [-65.0, -68.0, 180.0, 180.0]),
        ("tptt", 10.0, [-177.0, 68.0, 180.0, 180.0]),
        ("mttm", 9.0, [-67.0, 180.0, 180.0, -65.0]),
        ("pttt", 8.0, [62.0, 180.0, 180.0, 180.0]),
    ],
    "ARG": [
        ("mtt180", 20.0, [-67.0, 180.0, 180.0, 180.0]),
        ("ttt180", 16.0, [-177.0, 180.0, 180.0, 180.0]),
        ("mtm-85", 13.0, [-67.0, 180.0, -65.0, -85.0]),
        ("mtt85", 12.0, [-67.0, 180.0, 180.0, 85.0]),
        ("mmt180", 10.0, [-65.0, -68.0, 180.0, 180.0]),
        ("tpt180", 9.0, [-177.0, 68.0, 180.0, 180.0]),
        ("ptt180", 7.0, [62.0, 180.0, 180.0, 180.0]),
    ],
    "HIS": [
        ("m-70", 32.0, [-65.0, -75.0]),
        ("t60", 17.0, [-177.0, 65.0]),
        ("t-80", 16.0, [-177.0, -80.0]),
        ("m80", 13.0, [-65.0, 80.0]),
        ("p-80", 11.0, [62.0, -75.0]),
        ("p80", 7.0, [62.0, 80.0]),
    ],
    "PHE": [
        ("m-85", 44.0, [-65.0, -85.0]),
        ("t80", 33.0, [-177.0, 75.0]),
        ("p90", 13.0, [62.0, 90.0]),
        ("m-30", 7.0, [-65.0, -30.0]),
    ],
    "TYR": [
        ("m-85", 45.0, [-65.0, -85.0]),
        ("t80", 34.0, [-177.0, 75.0]),
        ("p90", 12.0, [62.0, 90.0]),
        ("m-30", 6.0, [-65.0, -30.0]),
    ],
    "TRP": [
        ("m95", 32.0, [-67.0, 95.0]),
        ("t-105", 18.0, [-177.0, -105.0]),
        ("t90", 16.0, [-177.0, 90.0]),
        ("m0", 12.0, [-67.0, -5.0]),
        ("p-90", 11.0, [62.0, -90.0]),
        ("m-90", 7.0, [-67.0, -90.0]),
    ],
    "PRO": [
        ("Cgamma-endo", 52.0, [28.0, -34.0]),
        ("Cgamma-exo", 48.0, [-25.0, 32.0]),
    ],
}

# Engh & Huber (2001) ideal backbone bond lengths (Å) and sigma (Å)
IDEAL_BACKBONE_BONDS: Dict[Tuple[str, str], Tuple[float, float]] = {
    ("N", "CA"): (1.458, 0.019),
    ("CA", "C"): (1.525, 0.021),
    ("C", "O"): (1.231, 0.020),
    ("CA", "CB"): (1.530, 0.020),
    ("C_prev", "N"): (1.329, 0.014),
}

# Engh & Huber (2001) ideal backbone bond angles (deg) and sigma (deg)
IDEAL_BACKBONE_ANGLES: Dict[Tuple[str, str, str], Tuple[float, float]] = {
    ("N", "CA", "C"): (111.2, 2.8),
    ("CA", "C", "O"): (120.8, 1.7),
    ("N", "CA", "CB"): (110.5, 1.7),
    ("C", "CA", "CB"): (110.1, 1.9),
    ("C_prev", "N", "CA"): (121.7, 1.8),
    ("CA_prev", "C_prev", "N"): (116.2, 2.0),
    ("O_prev", "C_prev", "N"): (123.0, 1.6),
}

# MolProbity / Bondi van der Waals radii (Å)
VDW_RADII: Dict[str, float] = {
    "H": 1.17, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "SE": 1.90, "ZN": 1.39, "MG": 1.73, "CA": 1.94, "FE": 1.40,
    "MN": 1.40, "CU": 1.40, "NA": 2.27, "K": 2.75,
}

# Model pKa values for titratable protein groups (PROPKA 3 reference)
MODEL_PKA: Dict[str, float] = {
    "ASP": 3.80,
    "GLU": 4.50,
    "HIS": 6.50,
    "CYS": 9.00,
    "TYR": 10.00,
    "LYS": 10.50,
    "ARG": 12.50,
    "N+": 8.00,
    "C-": 3.20,
}

# Formal charges of protonated vs deprotonated forms
TITRATABLE_CHARGES: Dict[str, Tuple[int, int]] = {
    # resname: (charge_when_protonated, charge_when_deprotonated)
    "ASP": (0, -1),
    "GLU": (0, -1),
    "HIS": (+1, 0),
    "CYS": (0, -1),
    "TYR": (0, -1),
    "LYS": (+1, 0),
    "ARG": (+1, 0),
    "N+": (+1, 0),
    "C-": (0, -1),
}

# Approximate empirical atomic partial charges for PQR / PDBQT fallback (AMBER-like)
BACKBONE_CHARGES: Dict[str, float] = {
    "N": -0.4157, "H": 0.2719, "HN": 0.2719, "CA": 0.0337, "HA": 0.0823,
    "C": 0.5973, "O": -0.5679, "OXT": -0.8055, "CB": -0.0500,
}
