# PandaProtein

**Deep Protein Preparation & Crystallographic Workbench for Molecular Dynamics & Docking**

[![License: MIT](https://img.shields.io/badge/License-MIT-007AFF.svg)](LICENSE)
[![Frontend: React + Tailwind CSS](https://img.shields.io/badge/Frontend-React%20%2B%20Tailwind%20CSS-38BDF8.svg)]()
[![Backend: FastAPI + Python](https://img.shields.io/badge/Backend-FastAPI%20%2B%20Gemmi%20%2B%20PROPKA%203-34C759.svg)]()
[![MD & Docking: GROMACS · CHARMM36m · AutoDock Vina](https://img.shields.io/badge/Ready%20for-GROMACS%20%C2%B7%20CHARMM36m%20%C2%B7%20Vina%2FGNINA-AF52DE.svg)]()

**PandaProtein** is a modern full-stack protein preparation studio built with **React, TypeScript, and Tailwind CSS** (styled with Apple Human Interface Guidelines) and a high-performance **FastAPI** structural biology backend. It unifies **COOT crystallographic validation & interactive model building**, **PDBFixer structural repair**, **PROPKA 3.5.1 protonation**, **RDKit ligand/pocket analysis**, and **GROMACS (CHARMM36m / AMBER) + AutoDock Vina / GNINA** export pipelines inside a single live 3D web workbench.

> **Repository:** [https://github.com/pritampanda15/PandaProtein](https://github.com/pritampanda15/PandaProtein)

---

## Core Capabilities

### 1. COOT-Inspired Crystallographic Validation & Interactive Model Building
- **Ramachandran Plot & Refinement**: Interactive clickable SVG Ramachandran plot calibrated to Richardson Top8000 contours (General, Glycine, Proline, Pre-Proline). Clicking any outlier focuses the 3D viewer directly on that residue; neighbor-safe peptide-plane rotation steers outliers cleanly into favored/allowed basins without disrupting adjacent residues.
- **Dunbrack / Richardson Rotamer Library & Auto-Fit**: Computes $\chi_1$–$\chi_4$ dihedrals, evaluates staggered conformational wells, scores every library rotamer against the 3D environment for clashes, and applies clash-free rotamers in one click.
- **Engh & Huber Geometry, Cis-Peptides & Chirality**: Bond length & angle RMSZ ($>4\sigma$ flags with Proline pyrrolidine geometry), $C_\beta$ deviations, cis-Pro vs cis-non-Pro / twisted peptide detection, and $C_\alpha$ L/D chirality validation + inversion repair.
- **Interactive COOT Model-Building Palette**:
  - **Mutate & Auto-Fit** (any of the 20 standard amino acids)
  - **180° Pepflip** (peptide carbonyl orientation flip)
  - **Cis ↔ Trans Peptide Conversion**
  - **180° NQH Flip** (Asn/Gln/His H-bond donor–acceptor flip)
  - **AltLoc Conformer Picker**
  - **COOT Regularize Zone** (local harmonic Engh & Huber + $C_\beta$ idealization)
  - **COOT Auto-Refine All** with monotonic rollback guards, automatic re-validation, and full **Undo / Redo** checkpoint history (`⌘Z` / `⌘⇧Z`).
- **All-Atom Steric Clashscore**: Fast KDTree MolProbity-calibrated overlap detection with 3D red clash spheres and automated side-chain + local relaxation clash relief.

### 2. Deep PDBFixer Pipeline
- Convert non-standard/modified residues (`MSE→MET`, `SEP→SER`, `TPO→THR`, `PTR→TYR`, `CSO/CSD→CYS`, etc.)
- Reconstruct missing heavy side-chain atoms (NeRF internal coordinates + clash-free rotamers) and C-terminal `OXT`
- Detect and rebuild missing internal loops (`SEQRES` + $C\text{–}N$ peptide breaks) with local regularization
- Resolve alternate conformations (`highest_occupancy`, `lowest_bfactor`, `keep_A`)
- Optional N/C-terminal capping (`ACE` / `NME`)
- Detect and assign disulfide bridges (`SG–SG` cutoff → `CYX` + `SSBOND`)
- Smart water & heterogen curation: retain active-site/bridging/buried waters, strip bulk waters and crystallization buffers (`SO4`, `EDO`, `GOL`, `PEG`, `DMS`, etc.) while keeping ligands, cofactors, and catalytic metals.

### 3. PROPKA 3.5.1 & Hydrogen-Bond Network Optimization
- Real **PROPKA 3.5.1** empirical $pK_a$ predictions with a live debounced **pH slider (0.0 – 14.0)**, interactive **SVG titration curve**, net charge, and **Isoelectric Point ($pI$)**.
- Automatic **metal-coordination locking** (`CYM`, `HID`/`HIE`) and **disulfide protection** (`CYX`).
- Per-residue manual overrides (`HID`, `HIE`, `HIP`, `ASH`, `GLH`, `CYM`, `CYX`, `LYN`, `TYM`, `ARN`), H-bond network NQH flips, and state-specific hydrogen building.

### 4. Chains, File & Residue Renaming Studio (with 3D Chain Spotlight)
- **Live 3D Chain Spotlight**: Click any chain to illuminate it in Apple Blue and dim other chains; check chains for deletion to highlight them in Red before deleting.
- Collision-safe **batch chain renaming & swapping** (`A ↔ B`), **Auto-Standardize (`A, B, C…`)**, **Reorder**, **Merge**, **Keep Only**, and **Single/Batch Chain Deletion**.
- Per-chain sequential residue renumbering or `+/-` offset shifting, ligand/heterogen 3-letter code renaming & chain reassignment (`BEN → LIG`, move to chain `L`), PDB filename & `TITLE` editing, canonicalizing protonation names, normalizing occupancies/B-factors, and a live PDB coordinate preview.

### 5. RDKit Ligands, Metals & 3D Druggable Binding Pocket Finder
- RDKit 2D vector SVG depictions, SMILES, Lipinski Rule-of-5 properties (MW, cLogP, TPSA, HBD/HBA, rotatable bonds), and **ETKDGv3 + MMFF94** 3D conformer generation (`.sdf`).
- 3D binding pocket contact mapping (H-bonds, hydrophobic contacts, metal coordination spheres).
- Built-in **3D Druggable Cavity Detector** (volume in $\text{Å}^3$, druggability score, hydrophobic %, lining residues, translucent 3D pocket cloud overlay, and automatic AutoDock Vina grid boxes for apo proteins).

### 6. Multi-Target Export Studio & Live GROMACS Verification
- Supports **PDB** (`.pdb`), **mmCIF** (`.cif`), **GROMACS `.gro`** (both **import & export** with nm conversion and box vectors), **RCSB PDB ID fetch**, and **AlphaFold DB UniProt fetch** — ships with **SARS-CoV-2 $\text{M}^\text{pro}$ (`6LU7`)** bundled locally as the instant default structure.
- **Force-Field-Aware GROMACS Pipeline**:
  - **CHARMM36m (`charmm36-jul2022` · Recommended)** with automatic `HID/HIE/HIP → HSD/HSE/HSP` nomenclature mapping and official CHARMM MDP nonbonded physics (`vdw-modifier = force-switch`, `rvdw-switch = 1.0`, `rvdw = 1.2`, `rcoulomb = 1.2`, `DispCorr = no`).
  - Also supports **AMBER99SB-ILDN**, **CHARMM27**, **OPLS-AA/L**, **GROMOS96 54a7**, and **AMBER03**.
  - **Live `gmx pdb2gmx` Pre-Flight Test** right from the browser + **Periodic Box, Solvation & Ion Calculator** (box volume, TIP3P waters, $\text{Na}^+/\text{Cl}^-$ counter-ions, total MD atom count).
- **Molecular Docking Export**: AutoDock Vina / GNINA receptor **`.pdbqt`** (polar H, partial charges, AD4 atom types) + auto-computed **`vina_config.txt`** with 3D wireframe Grid Box overlay, APBS **`.pqr`**, complete **GROMACS Production `.zip` Bundle**, and publication-grade Markdown/JSON reports.

### 7. Apple-Crafted 3D Workbench UX
- Persistent **3Dmol.js** WebGL viewport (`keepH: true`) with Cartoon, Sticks, Ball & Stick, and Surface modes, QA Issues Heatmap, B-factor/pLDDT coloring, 3D Disulfide/Metal/H-Bond cylinders (`S–S/HB`), 3D Clash spheres, Pocket clouds, and PNG snapshot export.
- Interactive **1-Letter Sequence Ribbon Strip** docked above the 3D viewport.
- **Apple Spotlight Command Palette (`⌘K`)** to jump straight to any residue (`HIS 18`, `A:41`), chain, ligand, or action; **`⌘Z` / `⌘⇧Z`** Undo/Redo; **Light & Dark Mode**; and disk-persisted sessions.

---

## Quick Start

### Prerequisites
- **Python 3.10+** (tested on Python 3.13)
- **Node.js 18+** & **npm** (for building or developing the frontend)
- **GROMACS** *(optional, for live `gmx pdb2gmx` verification & topology generation)*

### Installation

```bash
# 1. Clone the repository
git clone https://github.com/pritampanda15/PandaProtein.git
cd PandaProtein

# 2. Set up the Python virtual environment and install backend engines
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

# 3. Install frontend dependencies and build the production React UI
npm --prefix frontend install
npm --prefix frontend run build
```

### Run PandaProtein

```bash
chmod +x start.sh
./start.sh
```

Open **[http://127.0.0.1:8000](http://127.0.0.1:8000)** in your browser.

For live frontend hot-reload development (FastAPI `:8000` + Vite HMR `:5173`):

```bash
./start.sh --dev
```

---

## Project Architecture

```text
PandaProtein/
├── backend/
│   ├── main.py                      # FastAPI REST API + static UI serving
│   ├── models.py                    # Pydantic request/response schemas
│   ├── session_manager.py           # Disk-persisted sessions, Baseline diff & COOT Undo/Redo stack
│   ├── samples/
│   │   └── 6LU7.pdb                 # Bundled default SARS-CoV-2 Mpro + N3 benchmark structure
│   └── engines/
│       ├── chemistry_data.py        # Engh & Huber ideals, Dunbrack/Richardson rotamers, vdW radii
│       ├── geometry_utils.py        # NeRF builder, Rodrigues rotation, Top8000 Ramachandran contours
│       ├── structure_io.py          # PDB / mmCIF / GROMACS .gro I/O + RCSB & AlphaFold fetchers
│       ├── coot_engine.py           # COOT validation, model-building, & monotonic Auto-Refine
│       ├── pdbfixer_engine.py       # Deep PDBFixer diagnostics & structural repair pipeline
│       ├── propka_engine.py         # PROPKA 3.5.1 titration, protonation & H-bond optimization
│       ├── ligand_engine.py         # RDKit 2D depictions, Ro5, pocket contacts & ETKDGv3 conformers
│       ├── pocket_engine.py         # 3D concave druggable binding pocket detector
│       ├── rename_engine.py         # File, chain swap/delete/merge, residue & ligand renaming
│       └── export_engine.py         # CHARMM36m/AMBER GROMACS bundle, live gmx test, PDBQT, PQR
├── frontend/
│   └── src/
│       ├── App.tsx                  # Split-pane studio, ⌘K / ⌘Z / ⌘⇧Z shortcuts, theme manager
│       ├── api/client.ts            # Typed FastAPI client
│       ├── types/pdbprep.ts         # TypeScript structural biology models
│       └── components/
│           ├── TopBar.tsx           # Translucent macOS toolbar + Readiness indicator
│           ├── SidebarNav.tsx       # Studio navigation & live structure telemetry
│           ├── SequenceStrip.tsx    # Interactive 1-letter amino acid ribbon bar
│           ├── Viewer3D.tsx         # 3Dmol.js WebGL viewport (keepH, spotlight, pockets, S-S/HB)
│           ├── OverviewPanel.tsx    # Activity readiness rings & Before → After progress diff
│           ├── CootWorkbench.tsx    # Interactive COOT residue editor & rotamer explorer
│           ├── RenameChainsPanel.tsx# PDB file, chain spotlight/swap/delete & ligand renamer
│           ├── PdbFixerPanel.tsx    # Granular PDBFixer repair controls
│           ├── PropkaPanel.tsx      # Interactive pH slider, SVG titration curve & overrides
│           ├── ValidationPanel.tsx  # Interactive SVG Ramachandran plot & Engh-Huber charts
│           ├── ClashesPanel.tsx     # MolProbity clashscore, rotamer outliers & NQH flips
│           ├── LigandPanel.tsx      # RDKit 2D SVG, contacts, metals, pockets & conformers
│           ├── ExportPanel.tsx      # CHARMM36m/AMBER GROMACS test, solvation calc & downloads
│           ├── CommandPalette.tsx   # Apple Spotlight ⌘K quick jump & action launcher
│           └── ImportModal.tsx      # Upload (.pdb/.cif/.gro), RCSB PDB & AlphaFold DB fetch
├── start.sh                         # One-command launcher (with backend auto-reload)
├── requirements.txt                 # Python backend dependencies
└── LICENSE                          # MIT License
```

---

## Keyboard Shortcuts

| Shortcut | Action |
|---|---|
| `⌘K` / `Ctrl+K` | Open Spotlight Quick Jump & Command Palette |
| `⌘Z` / `Ctrl+Z` | Undo last COOT / PDBFixer / PROPKA / Chain edit |
| `⌘⇧Z` / `Ctrl+Y` | Redo edit |

---

## Author & Contact

- **Author**: Pritam Kumar Panda
- **Email**: pritam@stanford.edu
- **GitHub**: [https://github.com/pritampanda15/PandaProtein](https://github.com/pritampanda15/PandaProtein)

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
