# PDBPrep Usage Guide

This guide provides detailed instructions for using PDBPrep to prepare protein structures for molecular dynamics simulations.

## Table of Contents

1. [Getting Started](#getting-started)
2. [Main Interface](#main-interface)
3. [Step-by-Step Workflow](#step-by-step-workflow)
4. [Advanced Features](#advanced-features)
5. [Output Files Explained](#output-files-explained)
6. [Common Use Cases](#common-use-cases)
7. [Tips and Best Practices](#tips-and-best-practices)

## Getting Started

### Launching PDBPrep

```bash
python pdbprep_gui.py
```

Or if installed via setup.py:
```bash
pdbprep
```

### First Time Setup

Before using PDBPrep, ensure:
1. GROMACS is installed and accessible (`gmx --version`)
2. All Python dependencies are installed
3. You have PDB files ready for processing

## Main Interface

The PDBPrep interface consists of several tabs:

### 1. Input Tab
- **PDB File Selection**: Browse and select your input structure
- **Output Directory**: Choose where processed files will be saved
- **Quick Info**: View basic information about loaded structure

### 2. Configuration Tab
- **Force Field Selection**: Choose from AMBER, CHARMM, OPLS
- **Water Model**: Select TIP3P, TIP4P, or SPC/E
- **Protonation Settings**: Configure pH and protonation states
- **Processing Options**: Enable/disable various features

### 3. Processing Tab
- **Progress Monitor**: Real-time processing status
- **Log Output**: Detailed processing information
- **Control Buttons**: Start, pause, or cancel processing

### 4. Results Tab
- **Output Files**: Browse generated files
- **Analysis Reports**: View processing summaries
- **Visualization**: Preview structures (if PyMOL is available)

## Step-by-Step Workflow

### Basic Protein Preparation

#### Step 1: Load Your Structure
1. Click the **"Browse"** button in the Input tab
2. Navigate to your PDB file
3. Select the file (e.g., `6LU7.pdb`)
4. The structure information will be displayed

#### Step 2: Configure Settings
1. Switch to the **Configuration** tab
2. Select force field:
   - **AMBER99SB-ILDN**: Recommended for most proteins
   - **CHARMM27**: Good for membrane proteins
   - **OPLS-AA/L**: Alternative all-atom force field

3. Choose water model:
   - **TIP3P**: Most common, fastest
   - **TIP4P**: Better water properties
   - **SPC/E**: Good for specific conditions

4. Set pH value (typically 7.0-7.4 for physiological conditions)

#### Step 3: Advanced Options (Optional)
- **Remove Heteroatoms**: Clean up non-protein molecules
- **Add Missing Atoms**: Reconstruct missing heavy atoms
- **Generate Protonation States**: Automatic pH-based protonation
- **Process Ligands**: Extract and parameterize small molecules
- **Generate Conformers**: Create multiple ligand geometries

#### Step 4: Process
1. Switch to the **Processing** tab
2. Click **"Start Processing"**
3. Monitor progress in real-time
4. Wait for completion message

#### Step 5: Review Results
1. Switch to the **Results** tab
2. Browse output files
3. Read summary reports
4. Check for any warnings or errors

## Advanced Features

### Protonation State Management

#### Automatic Protonation
```
pH Settings:
- Target pH: 7.4
- Tolerance: ±0.5
- Method: PropKa or H++
```

#### Manual Residue Selection
For specific residues (e.g., HIS):
1. Open **Advanced Protonation Dialog**
2. Select residue from list
3. Choose protonation state:
   - HID (δ-protonated)
   - HIE (ε-protonated)
   - HIP (doubly protonated)

### Ligand Processing

When your PDB contains ligands:

1. **Automatic Detection**: Ligands are identified automatically
2. **Extraction**: Each ligand is saved separately
3. **Topology Generation**: GROMACS-compatible `.itp` files created
4. **Conformer Generation**: Multiple 3D conformations generated
5. **Property Calculation**:
   - Molecular weight
   - LogP (lipophilicity)
   - H-bond donors/acceptors
   - Lipinski's rule of five compliance

### Batch Processing

Process multiple structures:

```python
# Example batch script
import subprocess

pdb_files = ['6LU7.pdb', '2XNI.pdb', '1ABC.pdb']

for pdb in pdb_files:
    # Launch PDBPrep with specific settings
    # (Command-line mode - future feature)
    pass
```

## Output Files Explained

### Directory Structure
```
pdbprep_output/
└── pdbprep_20250918_141713/
    ├── structures/
    │   ├── prepared_6LU7_protein.pdb    # Clean protein structure
    │   ├── prepared_6LU7.gro            # GROMACS coordinates
    │   ├── prepared_6LU7.top            # Main topology file
    │   ├── topol_Protein_chain_A.itp    # Chain-specific topology
    │   ├── posre_Protein_chain_A.itp    # Position restraints
    │   └── processed.gro                # Final processed structure
    │
    ├── ligands/
    │   ├── ligand_02J_C_1.pdb           # Extracted ligand
    │   ├── ligand_02J_C_1.itp           # Ligand topology
    │   └── ligand_02J_C_1_conformers/   # Multiple conformations
    │       ├── ligand_02J_C_1_conf_1.pdb
    │       ├── ligand_02J_C_1_conf_2.pdb
    │       └── ligand_02J_C_1_conf_3.pdb
    │
    └── reports/
        ├── summary.txt                   # Processing summary
        ├── processing_report.json        # Detailed processing log
        ├── ligand_analysis.json          # Ligand properties
        └── ligand_summary.txt            # Ligand summary
```

### Key Files

#### 1. prepared_*.pdb
Clean protein structure with:
- Corrected atom names
- Standard residue names
- Proper chain identifiers
- Missing atoms added (if requested)

#### 2. prepared_*.gro
GROMACS coordinate file:
- Binary-compatible format
- Used for MD simulations
- Contains position data

#### 3. prepared_*.top
Main topology file:
- Force field parameters
- Molecular connectivity
- Include statements for chains and ligands

#### 4. ligand_*.itp
Ligand topology files:
- Atom types and charges
- Bond parameters
- Angle and dihedral terms
- Compatible with selected force field

#### 5. Reports
- **summary.txt**: Human-readable overview
- **processing_report.json**: Machine-readable log
- **ligand_analysis.json**: Chemical properties and drug-likeness

## Common Use Cases

### Use Case 1: Protein-Ligand Docking Preparation

**Objective**: Prepare a protein-ligand complex for docking studies

**Steps**:
1. Load protein-ligand complex PDB
2. Enable "Process Ligands" option
3. Enable "Generate Conformers" (set to 10-20 conformers)
4. Process structure
5. Use generated ligand conformers for docking
6. Use prepared protein for receptor preparation

**Output**: Multiple ligand conformers ready for AutoDock, Vina, or Glide

### Use Case 2: MD Simulation Setup

**Objective**: Prepare structure for GROMACS MD simulation

**Steps**:
1. Load protein PDB
2. Select AMBER99SB-ILDN force field
3. Choose TIP3P water model
4. Set pH to 7.4
5. Enable "Add Missing Atoms"
6. Process structure
7. Use generated .gro and .top files directly in GROMACS

**Next steps** (in GROMACS):
```bash
# Add to your GROMACS workflow
gmx editconf -f prepared_protein.gro -o boxed.gro -bt cubic -d 1.0
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
# ... continue with standard GROMACS workflow
```

### Use Case 3: Drug-Likeness Screening

**Objective**: Analyze ligands from PDB structures for drug-likeness

**Steps**:
1. Load PDB with multiple ligands
2. Enable "Process Ligands"
3. Process structure
4. Review `ligand_analysis.json` for:
   - Molecular weight
   - LogP
   - H-bond donors/acceptors
   - Lipinski violations

**Interpretation**:
- Lipinski's Rule of Five compliance
- Identify most drug-like candidates
- Flag potential ADME issues

## Tips and Best Practices

### Before Processing

1. **Validate Your PDB**:
   - Check for missing residues
   - Verify chain identifiers
   - Look for alternate conformations

2. **Check GROMACS Installation**:
   ```bash
   gmx --version
   ```

3. **Backup Original Files**: Keep pristine copies of input structures

### During Processing

1. **Monitor the Log**: Watch for warnings about:
   - Missing atoms
   - Unrecognized residues
   - Protonation issues

2. **Don't Interrupt**: Let processing complete fully

### After Processing

1. **Validate Output**:
   - Check atom counts match expectations
   - Verify topology file completeness
   - Review any warnings in reports

2. **Test in GROMACS**:
   ```bash
   gmx check -f prepared_protein.gro
   gmx check -s prepared_protein.top
   ```

3. **Visualize Results**:
   - Use PyMOL, VMD, or Chimera
   - Check for structural artifacts
   - Verify ligand positions

### Common Pitfalls

1. **Non-standard Residues**: Manually add parameters if needed
2. **Metal Ions**: May require special topology files
3. **Disulfide Bonds**: Verify correct connectivity
4. **Membrane Proteins**: Consider using CHARMM27 force field

## Troubleshooting

### Problem: "GROMACS not found"
**Solution**: Add GROMACS to PATH or reinstall

### Problem: "Ligand topology generation failed"
**Solution**:
- Check if ligand has valid 3D coordinates
- Try different force field
- Generate topology manually with ACPYPE or CGenFF

### Problem: "Missing heavy atoms"
**Solution**: Enable "Add Missing Atoms" option

### Problem: "Unusual protonation state"
**Solution**: Manually set protonation in Advanced Settings

## Additional Resources

- [GROMACS Documentation](http://manual.gromacs.org/)
- [Biopython Tutorial](https://biopython.org/wiki/Documentation)
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [PDB Format Specification](https://www.wwpdb.org/documentation/file-format)

## Getting Help

If you encounter issues:
1. Check this documentation
2. Review log files in output directory
3. Open an issue on GitHub with:
   - PDBPrep version
   - Input file (if possible)
   - Error messages
   - System information

## Appendix: Command Reference

### Force Fields
- **AMBER99SB-ILDN**: General-purpose protein force field
- **CHARMM27**: Versatile, good for membranes
- **OPLS-AA/L**: All-atom optimized potentials

### Water Models
- **TIP3P**: 3-site, fast, most common
- **TIP4P**: 4-site, better properties
- **SPC/E**: Extended simple point charge

### File Extensions
- `.pdb`: Protein Data Bank format
- `.gro`: GROMACS coordinate format
- `.top`: GROMACS topology file
- `.itp`: GROMACS include topology
- `.mdp`: GROMACS parameter file (not generated by PDBPrep)
