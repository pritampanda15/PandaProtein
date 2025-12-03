# PDBPrep

**Professional PDB Preparation Tool for GROMACS**

A comprehensive GUI application for preparing protein structures for molecular dynamics (MD) simulations with special focus on drug design applications.

## Features

- **User-Friendly GUI**: Built with PySide6 for a modern, intuitive interface
- **GROMACS Integration**: Seamless preparation of structures for GROMACS MD simulations
- **Ligand Processing**: Advanced ligand extraction and analysis
- **Protonation Control**: Sophisticated pH-dependent protonation state management
- **Multiple File Formats**: Support for PDB, GRO, and topology files
- **Batch Processing**: Handle multiple structures efficiently
- **Detailed Reports**: Comprehensive analysis and processing reports
- **Conformer Generation**: Generate multiple conformers for ligands

## Prerequisites

Before installing PDBPrep, ensure you have the following installed:

- **Python 3.8+**
- **GROMACS** (for MD preparation features)
- **Reduce** (optional, for advanced protonation)
- **Open Babel** (optional, for format conversions)

### Installing GROMACS

```bash
# Ubuntu/Debian
sudo apt-get install gromacs

# macOS (via Homebrew)
brew install gromacs

# Or build from source: http://www.gromacs.org/
```

## Installation

### Method 1: Basic Installation

1. Clone the repository:
```bash
git clone https://github.com/pritampanda15/pdbprep.git
cd pdbprep
```

2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

3. Run the application:
```bash
python pdbprep_gui_latest.py
```

### Method 2: Development Installation

```bash
git clone https://github.com/pritampanda15/pdbprep.git
cd pdbprep
pip install -e .
```

Then run:
```bash
pdbprep
```

### Method 3: Using Virtual Environment (Recommended)

```bash
# Create virtual environment
python -m venv pdbprep_env

# Activate virtual environment
# On macOS/Linux:
source pdbprep_env/bin/activate
# On Windows:
# pdbprep_env\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run application
python pdbprep_gui.py
```

## Quick Start

1. **Launch the application**:
   ```bash
   python pdbprep_gui.py
   ```

2. **Load a PDB file**: Click "Browse" and select your protein structure file

3. **Configure settings**:
   - Choose force field (AMBER, CHARMM, OPLS)
   - Set water model
   - Configure protonation settings

4. **Process**: Click "Start Processing" to prepare your structure

5. **Review outputs**: Find processed files in `pdbprep_output/` directory

## Usage Examples

### Basic Protein Preparation

1. Load your PDB file (e.g., `6LU7.pdb`)
2. Select AMBER99SB-ILDN force field
3. Choose TIP3P water model
4. Set pH to 7.4
5. Click "Start Processing"

Output files will be generated in `pdbprep_output/pdbprep_YYYYMMDD_HHMMSS/`:
- `structures/` - Prepared protein structures
- `ligands/` - Extracted and processed ligands
- `reports/` - Analysis and summary reports

### Processing with Ligands

When your PDB contains ligands (small molecules):
- Ligands are automatically extracted
- Individual topology files (`.itp`) are generated
- Multiple conformers can be generated for docking studies
- Drug-likeness properties are calculated (Lipinski's rule of five)

## Project Structure

```
pdbprep/
├── pdbprep_gui.py    # Main application
├── requirements.txt          # Python dependencies
├── README.md                 # This file
├── LICENSE                   # License information
├── setup.py                  # Installation script
├── .gitignore               # Git ignore rules
├── docs/                    # Additional documentation
│   └── USAGE.md            # Detailed usage guide
└── pdbprep_output/         # Output directory (created at runtime)
    └── pdbprep_YYYYMMDD_HHMMSS/
        ├── structures/      # Prepared protein structures
        ├── ligands/        # Extracted ligands and conformers
        └── reports/        # Analysis reports
```

## Configuration Options

### Force Fields
- AMBER99SB-ILDN (recommended for proteins)
- CHARMM27
- OPLS-AA/L

### Water Models
- TIP3P (default)
- TIP4P
- SPC/E

### Protonation Settings
- Target pH (0.0 - 14.0)
- pH tolerance
- Residue-specific protonation states

## Output Files

### Structure Files
- `prepared_*.pdb` - Cleaned protein structure
- `prepared_*.gro` - GROMACS coordinate file
- `prepared_*.top` - GROMACS topology file

### Ligand Files
- `ligand_*.pdb` - Individual ligand structures
- `ligand_*.itp` - Ligand topology files
- `ligand_*_conformers/` - Generated conformers

### Reports
- `summary.txt` - Processing summary
- `ligand_analysis.json` - Ligand properties
- `processing_report.json` - Detailed processing log

## Troubleshooting

### Common Issues

**"GROMACS not found"**
- Ensure GROMACS is installed and in your PATH
- Test with: `gmx --version`

**"Biopython not found"**
- Install with: `pip install biopython`

**"RDKit not found"**
- Install with: `pip install rdkit`
- Or via conda: `conda install -c conda-forge rdkit`

**GUI doesn't launch**
- Check PySide6 installation: `pip install --upgrade PySide6`
- Verify Python version: `python --version` (must be 3.8+)

### Getting Help

- Check the [documentation](docs/USAGE.md)
- Open an issue on GitHub
- Contact the maintainers

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- GROMACS development team
- Biopython contributors
- RDKit developers
- PySide6/Qt framework

## Contact

- **Author**: [Your Name]
- **Email**: [your.email@example.com]
- **GitHub**: [https://github.com/yourusername/pdbprep](https://github.com/yourusername/pdbprep)

## Roadmap

- [ ] Support for more force fields
- [ ] Integration with PyMOL for visualization
- [ ] Automated MD simulation setup
- [ ] Web-based interface
- [ ] Cloud processing support

## Version History

- **v1.0.0** - Initial release
  - Basic protein preparation
  - Ligand extraction
  - GROMACS integration
  - GUI interface
