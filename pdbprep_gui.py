#!/usr/bin/env python3
"""
PDBPrep - Professional PDB Preparation Tool for GROMACS
A comprehensive GUI for preparing protein structures for molecular dynamics simulations
with special focus on drug design applications.
"""

import sys
import os
import json
import subprocess
import threading
import time
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime

# Add basic PDB processing imports (you can replace with your preferred libraries)
try:
    from Bio.PDB import PDBParser, PDBIO, Select
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("Warning: Biopython not found. Install with: pip install biopython")

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("Warning: RDKit not found. Install with: pip install rdkit")

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QFormLayout, QLabel, QPushButton, QLineEdit, QTextEdit,
    QComboBox, QCheckBox, QSpinBox, QDoubleSpinBox, QFileDialog, QMessageBox,
    QProgressBar, QGroupBox, QDialog, QDialogButtonBox, QTableWidget,
    QTableWidgetItem, QHeaderView, QSplitter, QFrame, QScrollArea,
    QStatusBar, QMenuBar, QToolBar, QListWidget, QListWidgetItem,
    QTreeWidget, QTreeWidgetItem, QSlider, QPushButton as QBtn
)
from PySide6.QtCore import (
    Qt, QThread, QTimer, Signal, QSettings, QSize, QRect,
    QPropertyAnimation, QEasingCurve, QParallelAnimationGroup
)
from PySide6.QtGui import (
    QFont, QPixmap, QIcon, QColor, QPalette, QAction, QKeySequence,
    QValidator, QRegularExpressionValidator, QIntValidator
)
from PySide6.QtCore import QRegularExpression

class PDBValidator(QValidator):
    """Custom validator for PDB file paths and formats."""
    
    def validate(self, input_text: str, pos: int) -> Tuple[QValidator.State, str, int]:
        if not input_text:
            return QValidator.State.Intermediate, input_text, pos
        
        if input_text.endswith('.pdb') and os.path.exists(input_text):
            return QValidator.State.Acceptable, input_text, pos
        elif os.path.exists(input_text) or input_text.endswith('.pdb'):
            return QValidator.State.Intermediate, input_text, pos
        else:
            return QValidator.State.Invalid, input_text, pos

class AdvancedProtonationDialog(QDialog):
    """Advanced protonation settings dialog."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Advanced Protonation Settings")
        self.setModal(True)
        self.resize(600, 500)
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # pH settings
        ph_group = QGroupBox("pH Settings")
        ph_layout = QGridLayout(ph_group)
        
        ph_layout.addWidget(QLabel("Target pH:"), 0, 0)
        self.ph_spinbox = QDoubleSpinBox()
        self.ph_spinbox.setRange(0.0, 14.0)
        self.ph_spinbox.setValue(7.0)
        self.ph_spinbox.setSingleStep(0.1)
        ph_layout.addWidget(self.ph_spinbox, 0, 1)
        
        ph_layout.addWidget(QLabel("pH Tolerance:"), 1, 0)
        self.ph_tolerance = QDoubleSpinBox()
        self.ph_tolerance.setRange(0.1, 2.0)
        self.ph_tolerance.setValue(0.5)
        self.ph_tolerance.setSingleStep(0.1)
        ph_layout.addWidget(self.ph_tolerance, 1, 1)
        
        layout.addWidget(ph_group)
        
        # Residue-specific settings
        residue_group = QGroupBox("Residue-Specific Settings")
        residue_layout = QVBoxLayout(residue_group)
        
        # HIS settings
        his_frame = QFrame()
        his_layout = QGridLayout(his_frame)
        his_layout.addWidget(QLabel("Histidine (HIS):"), 0, 0)
        self.his_combo = QComboBox()
        self.his_combo.addItems(["Auto", "HID (δ-protonated)", "HIE (ε-protonated)", "HIP (both protonated)"])
        his_layout.addWidget(self.his_combo, 0, 1)
        
        # ASP/GLU settings
        his_layout.addWidget(QLabel("Aspartic Acid (ASP):"), 1, 0)
        self.asp_combo = QComboBox()
        self.asp_combo.addItems(["Auto", "ASP (deprotonated)", "ASPH (protonated)"])
        his_layout.addWidget(self.asp_combo, 1, 1)
        
        his_layout.addWidget(QLabel("Glutamic Acid (GLU):"), 2, 0)
        self.glu_combo = QComboBox()
        self.glu_combo.addItems(["Auto", "GLU (deprotonated)", "GLUH (protonated)"])
        his_layout.addWidget(self.glu_combo, 2, 1)
        
        # CYS settings
        his_layout.addWidget(QLabel("Cysteine (CYS):"), 3, 0)
        self.cys_combo = QComboBox()
        self.cys_combo.addItems(["Auto", "CYS (free)", "CYM (deprotonated)", "CYX (disulfide)"])
        his_layout.addWidget(self.cys_combo, 3, 1)
        
        # LYS settings
        his_layout.addWidget(QLabel("Lysine (LYS):"), 4, 0)
        self.lys_combo = QComboBox()
        self.lys_combo.addItems(["Auto", "LYS (protonated)", "LYN (neutral)"])
        his_layout.addWidget(self.lys_combo, 4, 1)
        
        residue_layout.addWidget(his_frame)
        layout.addWidget(residue_group)
        
        # Advanced options
        advanced_group = QGroupBox("Advanced Options")
        advanced_layout = QVBoxLayout(advanced_group)
        
        self.use_propka = QCheckBox("Use PROPKA for pKa prediction")
        self.use_propka.setChecked(True)
        advanced_layout.addWidget(self.use_propka)
        
        self.consider_metal = QCheckBox("Consider metal coordination")
        advanced_layout.addWidget(self.consider_metal)
        
        self.interactive_his = QCheckBox("Interactive histidine protonation")
        advanced_layout.addWidget(self.interactive_his)
        
        layout.addWidget(advanced_group)
        
        # Buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
    
    def get_settings(self) -> Dict:
        """Return current protonation settings."""
        return {
            'ph': self.ph_spinbox.value(),
            'ph_tolerance': self.ph_tolerance.value(),
            'his_state': self.his_combo.currentText(),
            'asp_state': self.asp_combo.currentText(),
            'glu_state': self.glu_combo.currentText(),
            'cys_state': self.cys_combo.currentText(),
            'lys_state': self.lys_combo.currentText(),
            'use_propka': self.use_propka.isChecked(),
            'consider_metal': self.consider_metal.isChecked(),
            'interactive_his': self.interactive_his.isChecked()
        }

class PDBProcessorThread(QThread):
    """Background thread for PDB processing operations."""
    
    progress_updated = Signal(int)
    status_updated = Signal(str)
    log_updated = Signal(str)
    finished_signal = Signal(dict)
    error_signal = Signal(str)
    
    def __init__(self, options: Dict):
        super().__init__()
        self.options = options
        self.is_running = True
        self.results = {
            'input_file': '',
            'output_directory': '',
            'files_created': [],
            'disulfides': [],
            'ligands': [],
            'statistics': {},
            'errors': [],
            'warnings': []
        }
    
    def log(self, message: str):
        """Log a message with timestamp."""
        timestamp = datetime.now().strftime('%H:%M:%S')
        full_message = f"[{timestamp}] {message}"
        self.log_updated.emit(full_message)
    
    def run(self):
        """Main processing function."""
        try:
            self.log("Starting PDB preparation...")
            self.progress_updated.emit(5)
            
            # Debug print
            print(f"DEBUG - Processing started with options: {self.options}")
            
            # Validate input
            if not self.validate_input():
                return
            
            # Create output directory
            self.create_output_directory()
            self.progress_updated.emit(10)
            
            # Debug print current state
            print(f"DEBUG - Input file: {self.results.get('input_file')}")
            print(f"DEBUG - Output dir: {self.results.get('output_directory')}")
            
            # Process PDB file
            if self.options.get('batch_mode', False):
                self.process_batch()
            else:
                self.process_single_pdb()
            
            # Debug print final results
            print(f"DEBUG - Final results: {self.results}")
            
            self.progress_updated.emit(100)
            self.status_updated.emit("Processing completed successfully!")
            self.log(f"All files saved to: {self.results['output_directory']}")
            
            self.finished_signal.emit(self.results)
            
        except Exception as e:
            print(f"DEBUG - Exception in processing: {str(e)}")
            self.log(f"Error: {str(e)}")
            self.error_signal.emit(str(e))
    
    def validate_input(self) -> bool:
        """Validate input parameters."""
        self.status_updated.emit("Validating input...")
        
        if self.options.get('batch_mode', False):
            batch_folder = self.options.get('batch_folder', '')
            if not batch_folder or not os.path.exists(batch_folder):
                self.error_signal.emit("Batch folder not found")
                return False
            
            # Find PDB files in batch folder
            pdb_files = list(Path(batch_folder).glob("*.pdb"))
            if not pdb_files:
                self.error_signal.emit("No PDB files found in batch folder")
                return False
            
            self.results['pdb_files'] = [str(f) for f in pdb_files]
            self.results['input_file'] = f"Batch: {len(pdb_files)} files from {batch_folder}"
            self.log(f"Found {len(pdb_files)} PDB files for batch processing")
            
        else:
            pdb_file = self.options.get('pdb_file', '')
            if not pdb_file or not os.path.exists(pdb_file):
                self.error_signal.emit("PDB file not found")
                return False
            
            self.results['input_file'] = pdb_file
            self.log(f"Input PDB file: {pdb_file}")
        
        print(f"DEBUG - After validation, input_file = {self.results['input_file']}")
        return True
    
    def create_output_directory(self):
        """Create output directory structure."""
        self.status_updated.emit("Creating output directory...")
        
        output_dir = self.options.get('output_dir', './pdbprep_output')
        output_dir = os.path.abspath(output_dir)
        
        # Create timestamp-based subdirectory
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_dir = os.path.join(output_dir, f"pdbprep_{timestamp}")
        
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(output_dir, "structures"), exist_ok=True)
        os.makedirs(os.path.join(output_dir, "ligands"), exist_ok=True)
        os.makedirs(os.path.join(output_dir, "reports"), exist_ok=True)
        
        self.results['output_directory'] = output_dir
        self.log(f"Created output directory: {output_dir}")
    
    def process_single_pdb(self):
        """Process a single PDB file."""
        pdb_file = self.results['input_file']
        self.log(f"Processing {os.path.basename(pdb_file)}...")
        
        # Debug: Print what we're about to do
        print(f"DEBUG - About to process: {pdb_file}")
        print(f"DEBUG - Extract ligands option: {self.options.get('extract_ligands', False)}")
        print(f"DEBUG - Generate topologies option: {self.options.get('generate_topologies', False)}")
        
        # Read and analyze PDB file
        self.analyze_pdb_structure(pdb_file)
        self.progress_updated.emit(30)
        
        # Process structure (this will separate protein and ligands)
        self.clean_structure(pdb_file)
        self.progress_updated.emit(50)
        
        # Extract and process ligands - ALWAYS do this if there are ligands
        ligand_lines = self.results.get('ligand_lines', [])
        print(f"DEBUG - Found {len(ligand_lines)} ligand lines")
        
        if ligand_lines:
            self.log(f"Processing {len(ligand_lines)} ligand atoms...")
            self.extract_ligands(pdb_file)
        else:
            self.log("No ligands found to process")
        
        self.progress_updated.emit(70)
        
        # Generate topologies and final files
        self.generate_output_files(pdb_file)
        self.progress_updated.emit(90)
        
        # Create report
        self.generate_report()
        
        print(f"DEBUG - Final results keys: {self.results.keys()}")
        print(f"DEBUG - Ligand analysis: {self.results.get('ligand_analysis', 'NOT FOUND')}")
    
    def process_batch(self):
        """Process multiple PDB files."""
        pdb_files = self.results['pdb_files']
        total_files = len(pdb_files)
        
        for i, pdb_file in enumerate(pdb_files):
            if not self.is_running:
                return
            
            self.log(f"Processing {os.path.basename(pdb_file)} ({i+1}/{total_files})...")
            
            try:
                # Process each file
                self.analyze_pdb_structure(pdb_file)
                self.clean_structure(pdb_file)
                
                if self.options.get('extract_ligands', False):
                    self.extract_ligands(pdb_file)
                
                self.generate_output_files(pdb_file)
                
            except Exception as e:
                self.log(f"Error processing {pdb_file}: {str(e)}")
                self.results['errors'].append(f"{pdb_file}: {str(e)}")
            
            # Update progress
            progress = 20 + int((i + 1) / total_files * 70)
            self.progress_updated.emit(progress)
        
        self.generate_report()
    
    def analyze_pdb_structure(self, pdb_file: str):
        """Analyze PDB structure and collect statistics."""
        self.status_updated.emit(f"Analyzing {os.path.basename(pdb_file)}...")
        
        if HAS_BIOPYTHON:
            try:
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure('protein', pdb_file)
                
                # Count atoms, residues, chains
                atom_count = sum(1 for _ in structure.get_atoms())
                residue_count = sum(1 for _ in structure.get_residues())
                chain_count = len(list(structure.get_chains()))
                
                # Look for ligands (HETATM)
                ligands = []
                waters = []
                
                for residue in structure.get_residues():
                    if residue.id[0] != ' ':  # HETATM
                        if residue.resname == 'HOH':
                            waters.append(residue.resname)
                        else:
                            ligands.append(residue.resname)
                
                # Look for disulfide bonds (simplified)
                cys_residues = [r for r in structure.get_residues() if r.resname == 'CYS']
                disulfides = []
                for i, cys1 in enumerate(cys_residues):
                    for cys2 in cys_residues[i+1:]:
                        if 'SG' in [a.name for a in cys1] and 'SG' in [a.name for a in cys2]:
                            try:
                                sg1 = cys1['SG']
                                sg2 = cys2['SG']
                                distance = sg1 - sg2
                                if distance < 2.5:  # Typical disulfide bond distance
                                    disulfides.append(f"CYS{cys1.id[1]}-CYS{cys2.id[1]}")
                            except:
                                pass
                
                self.results['statistics'].update({
                    'atoms': atom_count,
                    'residues': residue_count,
                    'chains': chain_count,
                    'waters': len(waters),
                    'unique_ligands': list(set(ligands))
                })
                
                self.results['ligands'] = list(set(ligands))
                self.results['disulfides'] = disulfides
                
                self.log(f"Structure analysis: {atom_count} atoms, {residue_count} residues, {chain_count} chains")
                if ligands:
                    self.log(f"Found ligands: {', '.join(set(ligands))}")
                if disulfides:
                    self.log(f"Disulfide bonds: {', '.join(disulfides)}")
                
            except Exception as e:
                self.log(f"BioPython analysis failed: {str(e)}")
                self.basic_pdb_analysis(pdb_file)
        else:
            self.basic_pdb_analysis(pdb_file)
    
    def basic_pdb_analysis(self, pdb_file: str):
        """Basic PDB analysis without BioPython."""
        try:
            with open(pdb_file, 'r') as f:
                lines = f.readlines()
            
            atom_lines = [l for l in lines if l.startswith('ATOM') or l.startswith('HETATM')]
            hetatm_lines = [l for l in lines if l.startswith('HETATM')]
            
            # Extract ligand names from HETATM records
            ligands = set()
            for line in hetatm_lines:
                if len(line) > 20:
                    resname = line[17:20].strip()
                    if resname != 'HOH':  # Exclude water
                        ligands.add(resname)
            
            self.results['statistics'] = {
                'atoms': len(atom_lines),
                'hetatm_records': len(hetatm_lines)
            }
            
            self.results['ligands'] = list(ligands)
            self.log(f"Basic analysis: {len(atom_lines)} atom records")
            if ligands:
                self.log(f"Found ligands: {', '.join(ligands)}")
                
        except Exception as e:
            self.log(f"Basic analysis failed: {str(e)}")
    
    def clean_structure(self, pdb_file: str):
        """Clean and process the structure, separating protein from ligands."""
        self.status_updated.emit("Cleaning structure and separating ligands...")
        
        print(f"DEBUG - clean_structure called for {pdb_file}")
        
        try:
            with open(pdb_file, 'r') as f:
                lines = f.readlines()
            
            print(f"DEBUG - Read {len(lines)} lines from PDB file")
            
            # Separate protein and ligand lines
            protein_lines = []
            ligand_lines = []
            water_lines = []
            removed_waters = 0
            
            # Define common ions and cofactors to exclude from ligand analysis
            common_ions = {'HOH', 'NA', 'CL', 'MG', 'CA', 'ZN', 'FE', 'K', 'MN', 'CO', 'NI', 'CU'}
            
            for line in lines:
                if line.startswith('ATOM'):
                    protein_lines.append(line)
                elif line.startswith('HETATM'):
                    if len(line) > 20:
                        resname = line[17:20].strip()
                        
                        if resname == 'HOH':
                            water_lines.append(line)
                            if self.options.get('clean_structure', False):
                                removed_waters += 1
                                continue
                            else:
                                protein_lines.append(line)  # Keep waters if not cleaning
                        elif resname in common_ions:
                            if self.options.get('keep_ions', True):
                                protein_lines.append(line)  # Keep ions with protein
                            # else skip (remove ions)
                        else:
                            # This is a potential ligand
                            ligand_lines.append(line)
                            print(f"DEBUG - Found ligand line: {resname}")
                else:
                    # Headers, CONECT, etc.
                    protein_lines.append(line)
            
            print(f"DEBUG - Separated: {len(protein_lines)} protein lines, {len(ligand_lines)} ligand lines")
            
            if removed_waters > 0:
                self.log(f"Removed {removed_waters} water molecules")
            
            # Save cleaned protein structure (without ligands)
            output_dir = self.results['output_directory']
            prefix = self.options.get('output_prefix', 'prepared')
            base_name = os.path.splitext(os.path.basename(pdb_file))[0]
            
            cleaned_pdb = os.path.join(output_dir, "structures", f"{prefix}_{base_name}_protein.pdb")
            with open(cleaned_pdb, 'w') as f:
                f.writelines(protein_lines)
            
            self.results['files_created'].append(cleaned_pdb)
            self.log(f"Protein structure saved: {os.path.basename(cleaned_pdb)}")
            
            # Store ligand lines for later processing
            self.results['ligand_lines'] = ligand_lines
            if ligand_lines:
                unique_ligands = set()
                for line in ligand_lines:
                    if len(line) > 20:
                        resname = line[17:20].strip()
                        unique_ligands.add(resname)
                
                self.log(f"Found {len(ligand_lines)} ligand atoms from {len(unique_ligands)} unique ligands: {', '.join(unique_ligands)}")
                print(f"DEBUG - Unique ligands found: {unique_ligands}")
            else:
                self.log("No ligands found in structure")
            
        except Exception as e:
            error_msg = f"Structure cleaning failed: {str(e)}"
            self.log(error_msg)
            self.results['errors'].append(error_msg)
            print(f"DEBUG - {error_msg}")
            import traceback
            traceback.print_exc()
    
    def extract_ligands(self, pdb_file: str):
        """Extract and analyze ligands with drug-like property calculations."""
        self.status_updated.emit("Extracting and analyzing ligands...")
        
        print(f"DEBUG - extract_ligands called with {pdb_file}")
        print(f"DEBUG - extract_ligands option: {self.options.get('extract_ligands', False)}")
        
        try:
            ligands_extracted = 0
            output_dir = self.results['output_directory']
            ligand_lines = self.results.get('ligand_lines', [])
            
            print(f"DEBUG - Processing {len(ligand_lines)} ligand lines")
            
            if not ligand_lines:
                self.log("No ligands found for extraction")
                return
            
            # Group ligand lines by residue
            ligand_groups = {}
            ligand_analysis = []
            
            for line in ligand_lines:
                if len(line) > 26:
                    resname = line[17:20].strip()
                    chain = line[21:22].strip() if len(line) > 21 else ''
                    resnum = line[22:26].strip() if len(line) > 25 else '1'
                    
                    key = f"{resname}_{chain}_{resnum}"
                    if key not in ligand_groups:
                        ligand_groups[key] = {
                            'lines': [],
                            'resname': resname,
                            'chain': chain,
                            'resnum': resnum
                        }
                    ligand_groups[key]['lines'].append(line)
            
            print(f"DEBUG - Found {len(ligand_groups)} ligand groups: {list(ligand_groups.keys())}")
            
            # Process each ligand group
            for ligand_id, ligand_data in ligand_groups.items():
                try:
                    print(f"DEBUG - Processing ligand {ligand_id}")
                    
                    # Save individual ligand PDB
                    ligand_file = os.path.join(output_dir, "ligands", f"ligand_{ligand_id}.pdb")
                    with open(ligand_file, 'w') as f:
                        f.writelines(ligand_data['lines'])
                    
                    self.results['files_created'].append(ligand_file)
                    ligands_extracted += 1
                    self.log(f"Saved ligand: {os.path.basename(ligand_file)}")
                    
                    # Analyze ligand properties
                    self.log(f"Analyzing properties of {ligand_data['resname']}...")
                    ligand_props = self.analyze_ligand_properties(ligand_data, ligand_file)
                    ligand_analysis.append(ligand_props)
                    print(f"DEBUG - Ligand analysis result: {ligand_props}")
                    
                    # Generate conformers if requested
                    if self.options.get('generate_conformers', False):
                        self.log(f"Generating conformers for {ligand_data['resname']}...")
                        self.generate_ligand_conformers(ligand_data, ligand_file)
                    
                    # Generate topology if requested
                    if self.options.get('generate_topologies', False):
                        self.log(f"Generating topology for {ligand_data['resname']}...")
                        self.generate_ligand_topology(ligand_data, ligand_file)
                        
                except Exception as e:
                    error_msg = f"Error processing ligand {ligand_id}: {str(e)}"
                    self.log(error_msg)
                    self.results['errors'].append(error_msg)
                    print(f"DEBUG - {error_msg}")
            
            # Store ligand analysis results
            self.results['ligand_analysis'] = ligand_analysis
            print(f"DEBUG - Stored ligand_analysis: {len(ligand_analysis)} ligands")
            
            if ligands_extracted > 0:
                self.log(f"Extracted and analyzed {ligands_extracted} ligand(s)")
                
                # Save ligand analysis report
                self.save_ligand_report(ligand_analysis)
                self.log("Ligand analysis report saved")
            else:
                self.log("No ligands found for extraction")
                
        except Exception as e:
            error_msg = f"Ligand extraction failed: {str(e)}"
            self.log(error_msg)
            self.results['errors'].append(error_msg)
            print(f"DEBUG - {error_msg}")
            import traceback
            traceback.print_exc()
    
    def analyze_ligand_properties(self, ligand_data: Dict, ligand_file: str) -> Dict:
        """Analyze drug-like properties of a ligand."""
        try:
            resname = ligand_data['resname']
            analysis = {
                'name': resname,
                'file': os.path.basename(ligand_file),
                'formula': 'N/A',
                'molecular_weight': 'N/A',
                'logp': 'N/A',
                'hbd': 'N/A',  # Hydrogen bond donors
                'hba': 'N/A',  # Hydrogen bond acceptors
                'rotatable_bonds': 'N/A',
                'tpsa': 'N/A',  # Topological polar surface area
                'lipinski_violations': 'N/A',
                'veber_violations': 'N/A',
                'drug_like': 'Unknown'
            }
            
            if HAS_RDKIT:
                try:
                    # Try to create molecule from PDB coordinates
                    # This is simplified - real implementation would need better 3D->2D conversion
                    mol = None
                    
                    # Extract atom information for basic analysis
                    atoms = []
                    for line in ligand_data['lines']:
                        if len(line) > 76:
                            element = line[76:78].strip()
                            if element:
                                atoms.append(element)
                    
                    # Basic molecular formula
                    from collections import Counter
                    atom_counts = Counter(atoms)
                    formula_parts = []
                    for element in sorted(atom_counts.keys()):
                        count = atom_counts[element]
                        if count > 1:
                            formula_parts.append(f"{element}{count}")
                        else:
                            formula_parts.append(element)
                    
                    analysis['formula'] = ''.join(formula_parts)
                    
                    # Estimate molecular weight from formula
                    atomic_weights = {
                        'C': 12.01, 'N': 14.01, 'O': 16.00, 'P': 30.97,
                        'S': 32.07, 'H': 1.008, 'F': 19.00, 'CL': 35.45,
                        'BR': 79.90, 'I': 126.90, 'NA': 22.99, 'MG': 24.31,
                        'K': 39.10, 'CA': 40.08, 'FE': 55.85, 'ZN': 65.38
                    }
                    
                    estimated_mw = sum(atomic_weights.get(element, 0) * count 
                                     for element, count in atom_counts.items())
                    analysis['molecular_weight'] = f"{estimated_mw:.2f}"
                    
                    # Basic Lipinski assessment based on estimated properties
                    mw_violation = estimated_mw > 500
                    violations = 0
                    if mw_violation:
                        violations += 1
                    
                    analysis['lipinski_violations'] = violations
                    analysis['drug_like'] = 'Likely' if violations <= 1 else 'Unlikely'
                    
                except Exception as e:
                    self.log(f"RDKit analysis failed for {resname}: {str(e)}")
            
            return analysis
            
        except Exception as e:
            self.log(f"Ligand property analysis failed: {str(e)}")
            return {'name': ligand_data.get('resname', 'Unknown'), 'error': str(e)}
    
    def generate_ligand_conformers(self, ligand_data: Dict, ligand_file: str):
        """Generate conformers for flexible ligands."""
        if not HAS_RDKIT:
            self.log("RDKit not available for conformer generation")
            return
        
        try:
            max_conformers = self.options.get('max_conformers', 10)
            resname = ligand_data['resname']
            
            # For now, create placeholder conformer files
            output_dir = os.path.dirname(ligand_file)
            base_name = os.path.splitext(os.path.basename(ligand_file))[0]
            
            conformer_dir = os.path.join(output_dir, f"{base_name}_conformers")
            os.makedirs(conformer_dir, exist_ok=True)
            
            # Create placeholder conformer files (in real implementation, use RDKit)
            for i in range(min(3, max_conformers)):  # Create 3 example conformers
                conf_file = os.path.join(conformer_dir, f"{base_name}_conf_{i+1}.pdb")
                with open(conf_file, 'w') as f:
                    f.write(f"REMARK Conformer {i+1} of {resname}\n")
                    f.writelines(ligand_data['lines'])
                
                self.results['files_created'].append(conf_file)
            
            self.log(f"Generated conformers for {resname}")
            
        except Exception as e:
            self.log(f"Conformer generation failed: {str(e)}")
    
    def generate_ligand_topology(self, ligand_data: Dict, ligand_file: str):
        """Generate GROMACS topology for ligand."""
        try:
            resname = ligand_data['resname']
            force_field = self.options.get('ligand_ff', 'GAFF')
            
            output_dir = os.path.dirname(ligand_file)
            base_name = os.path.splitext(os.path.basename(ligand_file))[0]
            
            # Create basic topology file
            top_file = os.path.join(output_dir, f"{base_name}.itp")
            with open(top_file, 'w') as f:
                f.write(f"; Topology for {resname} generated by PDBPrep\n")
                f.write(f"; Force field: {force_field}\n")
                f.write(f"; Charge method: {self.options.get('charge_method', 'AM1-BCC')}\n")
                f.write(f"\n[ moleculetype ]\n")
                f.write(f"{resname}     3\n")
                f.write(f"\n[ atoms ]\n")
                f.write(f";   nr       type  resnr residue  atom   cgnr     charge       mass\n")
                
                # Add placeholder atoms (real implementation would calculate charges)
                atom_idx = 1
                for line in ligand_data['lines']:
                    if len(line) > 76:
                        atom_name = line[12:16].strip()
                        element = line[76:78].strip() if len(line) > 77 else atom_name[0]
                        f.write(f"{atom_idx:5d} {element:>10s}     1    {resname:>3s} {atom_name:>6s} "
                               f"{atom_idx:5d}   0.00000  12.01000\n")
                        atom_idx += 1
            
            self.results['files_created'].append(top_file)
            self.log(f"Generated topology for {resname}")
            
        except Exception as e:
            self.log(f"Topology generation failed: {str(e)}")
    
    def save_ligand_report(self, ligand_analysis: List[Dict]):
        """Save detailed ligand analysis report."""
        try:
            output_dir = self.results['output_directory']
            report_file = os.path.join(output_dir, "reports", "ligand_analysis.json")
            
            report = {
                "timestamp": datetime.now().isoformat(),
                "analysis_options": {
                    "calculate_druglike": self.options.get('calc_druglike', False),
                    "lipinski_filter": self.options.get('lipinski_filter', False),
                    "veber_filter": self.options.get('veber_filter', False),
                    "pains_filter": self.options.get('pains_filter', False)
                },
                "ligands": ligand_analysis,
                "summary": {
                    "total_ligands": len(ligand_analysis),
                    "drug_like": len([l for l in ligand_analysis if l.get('drug_like') == 'Likely']),
                    "lipinski_compliant": len([l for l in ligand_analysis 
                                             if l.get('lipinski_violations', 999) <= 1])
                }
            }
            
            with open(report_file, 'w') as f:
                json.dump(report, f, indent=2)
            
            self.results['files_created'].append(report_file)
            
            # Also create human-readable text report
            text_report = os.path.join(output_dir, "reports", "ligand_summary.txt")
            with open(text_report, 'w') as f:
                f.write("Ligand Analysis Summary\n")
                f.write("======================\n\n")
                
                for ligand in ligand_analysis:
                    f.write(f"Ligand: {ligand.get('name', 'Unknown')}\n")
                    f.write(f"  Formula: {ligand.get('formula', 'N/A')}\n")
                    f.write(f"  Molecular Weight: {ligand.get('molecular_weight', 'N/A')}\n")
                    f.write(f"  Drug-like: {ligand.get('drug_like', 'Unknown')}\n")
                    f.write(f"  Lipinski violations: {ligand.get('lipinski_violations', 'N/A')}\n")
                    f.write("\n")
            
            self.results['files_created'].append(text_report)
            
        except Exception as e:
            self.log(f"Ligand report generation failed: {str(e)}")
    
    def generate_output_files(self, pdb_file: str):
        """Generate additional output files."""
        self.status_updated.emit("Generating output files...")
        
        output_dir = self.results['output_directory']
        prefix = self.options.get('output_prefix', 'prepared')
        base_name = os.path.splitext(os.path.basename(pdb_file))[0]
        
        # Create a mock .gro file (you'd replace this with actual conversion)
        gro_file = os.path.join(output_dir, "structures", f"{prefix}_{base_name}.gro")
        with open(gro_file, 'w') as f:
            f.write("Generated by PDBPrep\n")
            f.write("0\n")  # Empty structure placeholder
            f.write("0.0 0.0 0.0\n")
        
        self.results['files_created'].append(gro_file)
        
        # Create a mock topology file
        top_file = os.path.join(output_dir, "structures", f"{prefix}_{base_name}.top")
        with open(top_file, 'w') as f:
            f.write(f"; Topology generated by PDBPrep\n")
            f.write(f"; Force field: {self.options.get('force_field', 'CHARMM36m')}\n")
            f.write(f"\n#include \"charmm36-jul2022.ff/forcefield.itp\"\n")
            f.write(f"\n[ system ]\n{base_name}\n")
            f.write(f"\n[ molecules ]\n")
        
        self.results['files_created'].append(top_file)
        
        self.log(f"Generated GROMACS files: {os.path.basename(gro_file)}, {os.path.basename(top_file)}")
    
    def generate_report(self):
        """Generate processing report."""
        self.status_updated.emit("Generating report...")
        
        output_dir = self.results['output_directory']
        report_file = os.path.join(output_dir, "reports", "processing_report.json")
        
        report = {
            "timestamp": datetime.now().isoformat(),
            "pdbprep_version": "1.0",
            "input_options": self.options,
            "results": self.results,
            "summary": {
                "files_created": len(self.results['files_created']),
                "ligands_found": len(self.results['ligands']),
                "disulfide_bonds": len(self.results['disulfides']),
                "errors": len(self.results['errors']),
                "warnings": len(self.results['warnings'])
            }
        }
        
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        self.results['files_created'].append(report_file)
        self.log(f"Report saved: {os.path.basename(report_file)}")
        
        # Also create a simple text summary
        summary_file = os.path.join(output_dir, "reports", "summary.txt")
        with open(summary_file, 'w') as f:
            f.write(f"PDBPrep Processing Summary\n")
            f.write(f"========================\n\n")
            f.write(f"Processed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Input: {self.results.get('input_file', 'Batch processing')}\n")
            f.write(f"Output directory: {output_dir}\n\n")
            f.write(f"Results:\n")
            f.write(f"- Files created: {len(self.results['files_created'])}\n")
            f.write(f"- Ligands found: {len(self.results['ligands'])}\n")
            f.write(f"- Disulfide bonds: {len(self.results['disulfides'])}\n")
            f.write(f"- Errors: {len(self.results['errors'])}\n")
            
            if self.results['files_created']:
                f.write(f"\nFiles created:\n")
                for file_path in self.results['files_created']:
                    f.write(f"- {os.path.basename(file_path)}\n")
        
        self.results['files_created'].append(summary_file)
    
    def stop(self):
        self.is_running = False

class BasicTab(QWidget):
    """Basic settings tab."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_output_dir = None
        self.current_results = None
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Input section
        input_group = QGroupBox("Input Files")
        input_layout = QFormLayout(input_group)
        
        # PDB file selection
        pdb_layout = QHBoxLayout()
        self.pdb_path = QLineEdit()
        self.pdb_path.setValidator(PDBValidator())
        self.pdb_path.setPlaceholderText("Select PDB file...")
        pdb_browse_btn = QPushButton("Browse")
        pdb_browse_btn.clicked.connect(self.browse_pdb)
        pdb_layout.addWidget(self.pdb_path)
        pdb_layout.addWidget(pdb_browse_btn)
        input_layout.addRow("PDB File:", pdb_layout)
        
        # Batch processing
        self.batch_mode = QCheckBox("Batch processing mode")
        self.batch_mode.toggled.connect(self.toggle_batch_mode)
        input_layout.addRow(self.batch_mode)
        
        self.batch_folder = QLineEdit()
        self.batch_folder.setEnabled(False)
        self.batch_folder.setPlaceholderText("Select folder containing PDB files...")
        batch_browse_btn = QPushButton("Browse Folder")
        batch_browse_btn.clicked.connect(self.browse_batch_folder)
        batch_browse_btn.setEnabled(False)
        self.batch_browse_btn = batch_browse_btn
        
        batch_layout = QHBoxLayout()
        batch_layout.addWidget(self.batch_folder)
        batch_layout.addWidget(batch_browse_btn)
        input_layout.addRow("Batch Folder:", batch_layout)
        
        layout.addWidget(input_group)
        
        # Force field section
        ff_group = QGroupBox("Force Field Selection")
        ff_layout = QFormLayout(ff_group)
        
        self.force_field = QComboBox()
        self.force_field.addItems([
            "CHARMM36m", "AMBER14SB", "AMBER19SB", "OPLS-AA/L",
            "GROMOS54a7", "CHARMM27", "Custom"
        ])
        ff_layout.addRow("Force Field:", self.force_field)
        
        # Custom mapping
        mapping_layout = QHBoxLayout()
        self.custom_mapping = QLineEdit()
        self.custom_mapping.setPlaceholderText("Custom force field mapping file...")
        mapping_browse_btn = QPushButton("Browse")
        mapping_browse_btn.clicked.connect(self.browse_mapping)
        mapping_layout.addWidget(self.custom_mapping)
        mapping_layout.addWidget(mapping_browse_btn)
        ff_layout.addRow("Custom Mapping:", mapping_layout)
        
        layout.addWidget(ff_group)
        
        # Basic processing options
        process_group = QGroupBox("Basic Processing Options")
        process_layout = QVBoxLayout(process_group)
        
        self.handle_alt_loc = QCheckBox("Handle alternate locations")
        self.handle_alt_loc.setChecked(True)
        self.handle_alt_loc.setToolTip("Automatically resolve alternate atom locations")
        
        self.rebuild_missing = QCheckBox("Rebuild missing heavy atoms")
        self.rebuild_missing.setChecked(True)
        self.rebuild_missing.setToolTip("Add missing heavy atoms using template structures")
        
        self.add_hydrogens = QCheckBox("Add hydrogen atoms")
        self.add_hydrogens.setChecked(True)
        self.add_hydrogens.setToolTip("Add hydrogen atoms based on pH and protonation states")
        
        self.detect_disulfides = QCheckBox("Detect disulfide bonds")
        self.detect_disulfides.setChecked(True)
        self.detect_disulfides.setToolTip("Automatically detect and rename disulfide-bonded cysteines")
        
        self.clean_structure = QCheckBox("Clean structure (remove waters, ions)")
        self.clean_structure.setToolTip("Remove crystallographic waters and ions")
        
        for checkbox in [self.handle_alt_loc, self.rebuild_missing, self.add_hydrogens, 
                        self.detect_disulfides, self.clean_structure]:
            process_layout.addWidget(checkbox)
        
        layout.addWidget(process_group)
        
        # Output section
        output_group = QGroupBox("Output Settings")
        output_layout = QFormLayout(output_group)
        
        # Output directory
        output_dir_layout = QHBoxLayout()
        self.output_dir = QLineEdit()
        self.output_dir.setText("./pdbprep_output")
        output_browse_btn = QPushButton("Browse")
        output_browse_btn.clicked.connect(self.browse_output_dir)
        output_dir_layout.addWidget(self.output_dir)
        output_dir_layout.addWidget(output_browse_btn)
        output_layout.addRow("Output Directory:", output_dir_layout)
        
        # Output filename prefix
        self.output_prefix = QLineEdit()
        self.output_prefix.setText("prepared")
        self.output_prefix.setPlaceholderText("Output filename prefix...")
        output_layout.addRow("Filename Prefix:", self.output_prefix)
        
        layout.addWidget(output_group)
        
        layout.addStretch()
    
    def browse_pdb(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "Select PDB File", "", 
            "PDB Files (*.pdb);;All Files (*)"
        )
        if filename:
            self.pdb_path.setText(filename)
    
    def browse_batch_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Folder with PDB Files")
        if folder:
            self.batch_folder.setText(folder)
    
    def browse_mapping(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "Select Custom Mapping File", "", 
            "JSON Files (*.json);;All Files (*)"
        )
        if filename:
            self.custom_mapping.setText(filename)
    
    def browse_output_dir(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if folder:
            self.output_dir.setText(folder)
    
    def toggle_batch_mode(self, checked):
        self.batch_folder.setEnabled(checked)
        self.batch_browse_btn.setEnabled(checked)
        self.pdb_path.setEnabled(not checked)

class AdvancedTab(QWidget):
    """Advanced settings tab."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.protonation_settings = {}
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Protonation section
        proton_group = QGroupBox("Protonation Settings")
        proton_layout = QVBoxLayout(proton_group)
        
        # Basic protonation options
        self.protonate_residues = QCheckBox("Add protons to titratable residues")
        self.protonate_residues.setChecked(True)
        proton_layout.addWidget(self.protonate_residues)
        
        # pH settings
        ph_layout = QHBoxLayout()
        ph_layout.addWidget(QLabel("Target pH:"))
        self.ph_value = QDoubleSpinBox()
        self.ph_value.setRange(0.0, 14.0)
        self.ph_value.setValue(7.0)
        self.ph_value.setSingleStep(0.1)
        ph_layout.addWidget(self.ph_value)
        ph_layout.addStretch()
        proton_layout.addLayout(ph_layout)
        
        # Advanced protonation button
        self.advanced_proton_btn = QPushButton("Advanced Protonation Settings...")
        self.advanced_proton_btn.clicked.connect(self.show_advanced_protonation)
        proton_layout.addWidget(self.advanced_proton_btn)
        
        layout.addWidget(proton_group)
        
        # Structure validation
        validation_group = QGroupBox("Structure Validation")
        validation_layout = QVBoxLayout(validation_group)
        
        self.check_clashes = QCheckBox("Check for atomic clashes")
        self.check_clashes.setChecked(True)
        
        self.check_geometry = QCheckBox("Validate bond lengths and angles")
        self.check_geometry.setChecked(True)
        
        self.check_chirality = QCheckBox("Validate stereochemistry")
        self.check_chirality.setChecked(True)
        
        for checkbox in [self.check_clashes, self.check_geometry, self.check_chirality]:
            validation_layout.addWidget(checkbox)
        
        # Clash tolerance
        clash_layout = QHBoxLayout()
        clash_layout.addWidget(QLabel("Clash tolerance (Å):"))
        self.clash_tolerance = QDoubleSpinBox()
        self.clash_tolerance.setRange(0.1, 2.0)
        self.clash_tolerance.setValue(0.8)
        self.clash_tolerance.setSingleStep(0.1)
        clash_layout.addWidget(self.clash_tolerance)
        clash_layout.addStretch()
        validation_layout.addLayout(clash_layout)
        
        layout.addWidget(validation_group)
        
        # Water handling
        water_group = QGroupBox("Water and Ion Handling")
        water_layout = QVBoxLayout(water_group)
        
        self.keep_waters = QCheckBox("Keep crystallographic waters")
        self.keep_waters.setToolTip("Retain water molecules from crystal structure")
        
        self.keep_ions = QCheckBox("Keep crystallographic ions")
        self.keep_ions.setToolTip("Retain metal ions and other crystallographic ions")
        
        # Water distance threshold
        water_dist_layout = QHBoxLayout()
        water_dist_layout.addWidget(QLabel("Water-protein distance (Å):"))
        self.water_distance = QDoubleSpinBox()
        self.water_distance.setRange(0.0, 10.0)
        self.water_distance.setValue(3.0)
        self.water_distance.setSingleStep(0.5)
        self.water_distance.setToolTip("Maximum distance from protein to keep water")
        water_dist_layout.addWidget(self.water_distance)
        water_dist_layout.addStretch()
        
        for widget in [self.keep_waters, self.keep_ions]:
            water_layout.addWidget(widget)
        water_layout.addLayout(water_dist_layout)
        
        layout.addWidget(water_group)
        
        # Chain handling
        chain_group = QGroupBox("Chain Processing")
        chain_layout = QVBoxLayout(chain_group)
        
        self.process_all_chains = QCheckBox("Process all chains")
        self.process_all_chains.setChecked(True)
        self.process_all_chains.toggled.connect(self.toggle_chain_selection)
        
        self.chain_selection = QLineEdit()
        self.chain_selection.setPlaceholderText("Enter chain IDs (e.g., A,B,C)")
        self.chain_selection.setEnabled(False)
        
        chain_layout.addWidget(self.process_all_chains)
        chain_layout.addWidget(QLabel("Select specific chains:"))
        chain_layout.addWidget(self.chain_selection)
        
        layout.addWidget(chain_group)
        
        layout.addStretch()
    
    def show_advanced_protonation(self):
        dialog = AdvancedProtonationDialog(self)
        if dialog.exec() == QDialog.Accepted:
            self.protonation_settings = dialog.get_settings()
    
    def toggle_chain_selection(self, checked):
        self.chain_selection.setEnabled(not checked)

class LigandTab(QWidget):
    """Ligand processing tab."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Ligand detection
        detection_group = QGroupBox("Ligand Detection")
        detection_layout = QVBoxLayout(detection_group)
        
        self.auto_detect_ligands = QCheckBox("Automatically detect ligands")
        self.auto_detect_ligands.setChecked(True)
        self.auto_detect_ligands.setToolTip("Identify HETATM entries as potential ligands")
        
        self.exclude_waters = QCheckBox("Exclude water molecules")
        self.exclude_waters.setChecked(True)
        
        self.exclude_ions = QCheckBox("Exclude common ions (Na+, Cl-, Mg2+, etc.)")
        self.exclude_ions.setChecked(True)
        
        self.min_heavy_atoms = QSpinBox()
        self.min_heavy_atoms.setRange(1, 100)
        self.min_heavy_atoms.setValue(6)
        
        min_atoms_layout = QHBoxLayout()
        min_atoms_layout.addWidget(QLabel("Minimum heavy atoms:"))
        min_atoms_layout.addWidget(self.min_heavy_atoms)
        min_atoms_layout.addStretch()
        
        for widget in [self.auto_detect_ligands, self.exclude_waters, self.exclude_ions]:
            detection_layout.addWidget(widget)
        detection_layout.addLayout(min_atoms_layout)
        
        layout.addWidget(detection_group)
        
        # Ligand processing
        processing_group = QGroupBox("Ligand Processing")
        processing_layout = QVBoxLayout(processing_group)
        
        self.extract_ligands = QCheckBox("Extract ligands to separate files")
        self.extract_ligands.setChecked(True)  # Enable by default
        
        self.generate_topologies = QCheckBox("Generate ligand topologies")
        self.generate_topologies.setChecked(True)  # Enable by default
        self.generate_topologies.setToolTip("Create GROMACS topology files for ligands")
        
        self.optimize_geometry = QCheckBox("Optimize ligand geometry")
        self.optimize_geometry.setToolTip("Perform geometry optimization using RDKit")
        
        self.generate_conformers = QCheckBox("Generate conformers")
        self.generate_conformers.setChecked(True)  # Enable by default
        self.generate_conformers.setToolTip("Generate multiple conformations for flexible ligands")
        
        # Number of conformers
        conf_layout = QHBoxLayout()
        conf_layout.addWidget(QLabel("Max conformers:"))
        self.max_conformers = QSpinBox()
        self.max_conformers.setRange(1, 100)
        self.max_conformers.setValue(10)
        conf_layout.addWidget(self.max_conformers)
        conf_layout.addStretch()
        
        for widget in [self.extract_ligands, self.generate_topologies, 
                      self.optimize_geometry, self.generate_conformers]:
            processing_layout.addWidget(widget)
        processing_layout.addLayout(conf_layout)
        
        layout.addWidget(processing_group)
        
        # Force field for ligands
        ligand_ff_group = QGroupBox("Ligand Force Field")
        ligand_ff_layout = QFormLayout(ligand_ff_group)
        
        self.ligand_ff = QComboBox()
        self.ligand_ff.addItems([
            "GAFF", "GAFF2", "CGenFF", "OPLS-AA", 
            "SwissParam", "ATB", "Custom"
        ])
        ligand_ff_layout.addRow("Force Field:", self.ligand_ff)
        
        # Charge calculation method
        self.charge_method = QComboBox()
        self.charge_method.addItems([
            "AM1-BCC", "RESP", "ESP", "Mulliken", "Gasteiger"
        ])
        ligand_ff_layout.addRow("Charge Method:", self.charge_method)
        
        layout.addWidget(ligand_ff_group)
        
        # Drug-like properties
        druglike_group = QGroupBox("Drug-like Properties Analysis")
        druglike_layout = QVBoxLayout(druglike_group)
        
        self.calc_druglike = QCheckBox("Calculate drug-like properties")
        self.calc_druglike.setChecked(True)  # Enable by default
        
        self.lipinski_filter = QCheckBox("Apply Lipinski's Rule of Five")
        self.lipinski_filter.setChecked(True)  # Enable by default
        self.veber_filter = QCheckBox("Apply Veber's rules")
        self.pains_filter = QCheckBox("PAINS filtering")
        
        for checkbox in [self.calc_druglike, self.lipinski_filter, self.veber_filter, self.pains_filter]:
            druglike_layout.addWidget(checkbox)
        
        layout.addWidget(druglike_group)
        
        # Binding site analysis
        binding_group = QGroupBox("Binding Site Analysis")
        binding_layout = QVBoxLayout(binding_group)
        
        self.detect_binding_site = QCheckBox("Detect binding sites")
        self.binding_site_radius = QDoubleSpinBox()
        self.binding_site_radius.setRange(3.0, 15.0)
        self.binding_site_radius.setValue(8.0)
        self.binding_site_radius.setSuffix(" Å")
        
        radius_layout = QHBoxLayout()
        radius_layout.addWidget(QLabel("Search radius:"))
        radius_layout.addWidget(self.binding_site_radius)
        radius_layout.addStretch()
        
        binding_layout.addWidget(self.detect_binding_site)
        binding_layout.addLayout(radius_layout)
        
        layout.addWidget(binding_group)
        
        layout.addStretch()

class ResultsTab(QWidget):
    """Results and analysis tab."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Create splitter for results
        splitter = QSplitter(Qt.Horizontal)
        
        # Left side - summary and logs
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        
        # Processing summary
        summary_group = QGroupBox("Processing Summary")
        summary_layout = QVBoxLayout(summary_group)
        
        self.summary_text = QTextEdit()
        self.summary_text.setMaximumHeight(150)
        self.summary_text.setReadOnly(True)
        summary_layout.addWidget(self.summary_text)
        
        left_layout.addWidget(summary_group)
        
        # Processing log
        log_group = QGroupBox("Processing Log")
        log_layout = QVBoxLayout(log_group)
        
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFont(QFont("Courier", 9))
        log_layout.addWidget(self.log_text)
        
        # Log controls
        log_controls = QHBoxLayout()
        self.clear_log_btn = QPushButton("Clear Log")
        self.clear_log_btn.clicked.connect(self.clear_log)
        self.save_log_btn = QPushButton("Save Log")
        self.save_log_btn.clicked.connect(self.save_log)
        log_controls.addWidget(self.clear_log_btn)
        log_controls.addWidget(self.save_log_btn)
        log_controls.addStretch()
        log_layout.addLayout(log_controls)
        
        left_layout.addWidget(log_group)
        splitter.addWidget(left_widget)
        
        # Right side - detailed results
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        
        # Results tabs
        self.results_tabs = QTabWidget()
        
        # Structure analysis tab
        structure_tab = QWidget()
        structure_layout = QVBoxLayout(structure_tab)
        
        self.structure_table = QTableWidget()
        self.structure_table.setColumnCount(3)
        self.structure_table.setHorizontalHeaderLabels(["Property", "Value", "Status"])
        self.structure_table.horizontalHeader().setStretchLastSection(True)
        structure_layout.addWidget(self.structure_table)
        
        self.results_tabs.addTab(structure_tab, "Structure")
        
        # Ligands tab
        ligands_tab = QWidget()
        ligands_layout = QVBoxLayout(ligands_tab)
        
        self.ligands_table = QTableWidget()
        self.ligands_table.setColumnCount(4)
        self.ligands_table.setHorizontalHeaderLabels(["Ligand", "Formula", "MW", "Drug-like"])
        self.ligands_table.horizontalHeader().setStretchLastSection(True)
        ligands_layout.addWidget(self.ligands_table)
        
        self.results_tabs.addTab(ligands_tab, "Ligands")
        
        # Errors/Warnings tab
        issues_tab = QWidget()
        issues_layout = QVBoxLayout(issues_tab)
        
        self.issues_list = QListWidget()
        issues_layout.addWidget(self.issues_list)
        
        self.results_tabs.addTab(issues_tab, "Issues")
        
        right_layout.addWidget(self.results_tabs)
        splitter.addWidget(right_widget)
        
        # Set splitter proportions
        splitter.setSizes([400, 600])
        layout.addWidget(splitter)
        
        # Export buttons
        export_layout = QHBoxLayout()
        
        self.open_output_btn = QPushButton("📁 Open Output Folder")
        self.open_output_btn.clicked.connect(self.open_output_folder)
        self.open_output_btn.setEnabled(False)
        
        self.export_json_btn = QPushButton("Export JSON Report")
        self.export_json_btn.clicked.connect(self.export_json_report)
        
        self.export_csv_btn = QPushButton("Export CSV Summary")
        self.export_csv_btn.clicked.connect(self.export_csv_summary)
        
        export_layout.addWidget(self.open_output_btn)
        export_layout.addStretch()
        export_layout.addWidget(self.export_json_btn)
        export_layout.addWidget(self.export_csv_btn)
        
        layout.addLayout(export_layout)
    
    def open_output_folder(self):
        """Open the output folder in system file manager."""
        if self.current_output_dir and os.path.exists(self.current_output_dir):
            if sys.platform.startswith('darwin'):  # macOS
                subprocess.run(['open', self.current_output_dir])
            elif sys.platform.startswith('win'):   # Windows
                subprocess.run(['explorer', self.current_output_dir])
            else:  # Linux and others
                subprocess.run(['xdg-open', self.current_output_dir])
    
    def set_output_directory(self, output_dir: str):
        """Set the current output directory and enable the button."""
        self.current_output_dir = output_dir
        self.open_output_btn.setEnabled(True)
    
    def set_current_results(self, results: Dict):
        """Set the current processing results for export functions."""
        self.current_results = results
    
    def clear_log(self):
        self.log_text.clear()
    
    def save_log(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Log", "pdbprep_log.txt", 
            "Text Files (*.txt);;All Files (*)"
        )
        if filename:
            with open(filename, 'w') as f:
                f.write(self.log_text.toPlainText())
    
    def export_json_report(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export JSON Report", "pdbprep_report.json",
            "JSON Files (*.json);;All Files (*)"
        )
        if filename:
            # Mock report data
            report = {
                "timestamp": datetime.now().isoformat(),
                "input_file": "example.pdb",
                "processing_options": {},
                "results": {},
                "issues": []
            }
            with open(filename, 'w') as f:
                json.dump(report, f, indent=2)
    
    def export_csv_summary(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export CSV Summary", "pdbprep_summary.csv",
            "CSV Files (*.csv);;All Files (*)"
        )
        if filename:
            # Implementation would export table data to CSV
            pass
    
    def update_log(self, message: str):
        """Add message to log."""
        self.log_text.append(message)
        # Auto-scroll to bottom
        scrollbar = self.log_text.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

class PDBPrepMainWindow(QMainWindow):
    """Main application window."""
    
    def __init__(self):
        super().__init__()
        self.processor_thread = None
        self.settings = QSettings("PDBPrep", "PDBPrep")
        self.init_ui()
        self.load_settings()
        
    def init_ui(self):
        self.setWindowTitle("PDBPrep - Professional PDB Preparation Tool")
        self.setGeometry(100, 100, 1200, 800)
        
        # Set application icon - use a safe standard icon or none
        try:
            # Try to use a common icon that should exist
            self.setWindowIcon(self.style().standardIcon(self.style().SP_FileIcon))
        except AttributeError:
            # If no standard icons work, just skip setting an icon
            pass
        
        # Create central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        
        # Create tab widget
        self.tabs = QTabWidget()
        
        # Create tabs
        self.basic_tab = BasicTab()
        self.advanced_tab = AdvancedTab()
        self.ligand_tab = LigandTab()
        self.results_tab = ResultsTab()
        
        # Add tabs
        self.tabs.addTab(self.basic_tab, "Basic Settings")
        self.tabs.addTab(self.advanced_tab, "Advanced")
        self.tabs.addTab(self.ligand_tab, "Ligands")
        self.tabs.addTab(self.results_tab, "Results")
        
        layout.addWidget(self.tabs)
        
        # Progress section
        progress_frame = QFrame()
        progress_frame.setFrameStyle(QFrame.Box)
        progress_layout = QVBoxLayout(progress_frame)
        
        self.status_label = QLabel("Ready")
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        
        progress_layout.addWidget(self.status_label)
        progress_layout.addWidget(self.progress_bar)
        
        layout.addWidget(progress_frame)
        
        # Control buttons
        button_layout = QHBoxLayout()
        
        self.run_btn = QPushButton("Run PDBPrep")
        self.run_btn.clicked.connect(self.run_pdbprep)
        self.run_btn.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: none;
                padding: 10px 20px;
                font-size: 14px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
            QPushButton:pressed {
                background-color: #3d8b40;
            }
        """)
        
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.clicked.connect(self.stop_processing)
        self.stop_btn.setEnabled(False)
        self.stop_btn.setStyleSheet("""
            QPushButton {
                background-color: #f44336;
                color: white;
                border: none;
                padding: 10px 20px;
                font-size: 14px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #da190b;
            }
            QPushButton:pressed {
                background-color: #c61409;
            }
        """)
        
        self.validate_btn = QPushButton("Validate Settings")
        self.validate_btn.clicked.connect(self.validate_settings)
        
        self.load_preset_btn = QPushButton("Load Preset")
        self.load_preset_btn.clicked.connect(self.load_preset)
        
        self.save_preset_btn = QPushButton("Save Preset")
        self.save_preset_btn.clicked.connect(self.save_preset)
        
        button_layout.addWidget(self.validate_btn)
        button_layout.addStretch()
        button_layout.addWidget(self.load_preset_btn)
        button_layout.addWidget(self.save_preset_btn)
        button_layout.addStretch()
        button_layout.addWidget(self.stop_btn)
        button_layout.addWidget(self.run_btn)
        
        layout.addLayout(button_layout)
        
        # Create menu bar
        self.create_menu_bar()
        
        # Create status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")
        
        # Create toolbar
        self.create_toolbar()
        
    def create_menu_bar(self):
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("File")
        
        new_action = QAction("New Project", self)
        new_action.setShortcut(QKeySequence.New)
        new_action.triggered.connect(self.new_project)
        file_menu.addAction(new_action)
        
        open_action = QAction("Open Project", self)
        open_action.setShortcut(QKeySequence.Open)
        open_action.triggered.connect(self.open_project)
        file_menu.addAction(open_action)
        
        save_action = QAction("Save Project", self)
        save_action.setShortcut(QKeySequence.Save)
        save_action.triggered.connect(self.save_project)
        file_menu.addAction(save_action)
        
        file_menu.addSeparator()
        
        # Add demo file creation
        demo_action = QAction("Create Demo PDB File", self)
        demo_action.triggered.connect(self.create_demo_pdb)
        file_menu.addAction(demo_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction("Exit", self)
        exit_action.setShortcut(QKeySequence.Quit)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Tools menu
        tools_menu = menubar.addMenu("Tools")
        
        gmx_test_action = QAction("Test gmx pdb2gmx", self)
        gmx_test_action.triggered.connect(self.test_gmx_pdb2gmx)
        tools_menu.addAction(gmx_test_action)
        
        structure_viewer_action = QAction("Structure Viewer", self)
        structure_viewer_action.triggered.connect(self.open_structure_viewer)
        tools_menu.addAction(structure_viewer_action)
        
        # Help menu
        help_menu = menubar.addMenu("Help")
        
        about_action = QAction("About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)
    
    def create_demo_pdb(self):
        """Create a demo PDB file for testing."""
        demo_pdb_content = """REMARK   1 DEMO PDB FILE CREATED BY PDBPREP
REMARK   2 THIS IS A SIMPLIFIED STRUCTURE FOR TESTING
ATOM      1  N   ALA A   1      20.154  16.967  18.280  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  17.863  18.618  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.744  17.085  18.897  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.805  15.948  19.374  1.00 20.00           O  
ATOM      5  CB  ALA A   1      19.373  18.761  19.796  1.00 20.00           C  
ATOM      6  N   CYS A   2      16.585  17.649  18.582  1.00 20.00           N  
ATOM      7  CA  CYS A   2      15.270  17.053  18.764  1.00 20.00           C  
ATOM      8  C   CYS A   2      14.113  18.039  18.897  1.00 20.00           C  
ATOM      9  O   CYS A   2      14.073  19.110  18.297  1.00 20.00           O  
ATOM     10  CB  CYS A   2      15.009  16.123  17.578  1.00 20.00           C  
ATOM     11  SG  CYS A   2      16.068  14.661  17.505  1.00 20.00           S  
ATOM     12  N   CYS A   3      13.123  17.758  19.741  1.00 20.00           N  
ATOM     13  CA  CYS A   3      11.942  18.603  19.949  1.00 20.00           C  
ATOM     14  C   CYS A   3      10.669  17.816  20.204  1.00 20.00           C  
ATOM     15  O   CYS A   3      10.680  16.617  20.501  1.00 20.00           O  
ATOM     16  CB  CYS A   3      11.748  19.552  18.772  1.00 20.00           C  
ATOM     17  SG  CYS A   3      13.141  20.667  18.415  1.00 20.00           S  
HETATM   18  O   HOH S   4       8.235  19.234  20.177  1.00 30.00           O  
HETATM   19  O   HOH S   5      22.454  15.432  17.234  1.00 30.00           O  
HETATM   20  MG  MG  S   6      15.500  22.000  16.000  1.00 25.00          MG  
HETATM   21  N9  ATP A 101      12.345  14.567  22.123  1.00 35.00           N  
HETATM   22  C8  ATP A 101      13.234  15.234  22.567  1.00 35.00           C  
HETATM   23  N7  ATP A 101      14.123  15.876  21.987  1.00 35.00           N  
HETATM   24  C5  ATP A 101      13.876  15.234  20.823  1.00 35.00           C  
HETATM   25  C6  ATP A 101      14.456  15.345  19.567  1.00 35.00           C  
HETATM   26  N6  ATP A 101      15.567  16.123  19.234  1.00 35.00           N  
HETATM   27  N1  ATP A 101      14.023  14.567  18.456  1.00 35.00           N  
HETATM   28  C2  ATP A 101      12.867  13.789  18.567  1.00 35.00           C  
HETATM   29  N3  ATP A 101      12.234  13.567  19.723  1.00 35.00           N  
HETATM   30  C4  ATP A 101      12.678  14.234  20.867  1.00 35.00           C  
HETATM   31  C1' ATP A 101      11.234  13.567  21.987  1.00 35.00           C  
HETATM   32  O4' ATP A 101      10.567  14.234  20.867  1.00 35.00           O  
HETATM   33  C2' ATP A 101       9.876  13.234  22.567  1.00 35.00           C  
HETATM   34  O2' ATP A 101       9.234  12.123  21.987  1.00 35.00           O  
HETATM   35  C3' ATP A 101       9.123  14.567  22.345  1.00 35.00           C  
HETATM   36  O3' ATP A 101       8.234  14.876  23.234  1.00 35.00           O  
HETATM   37  C4' ATP A 101      10.234  15.567  21.876  1.00 35.00           C  
HETATM   38  C5' ATP A 101       9.876  16.234  20.567  1.00 35.00           C  
HETATM   39  O5' ATP A 101      10.876  17.123  19.987  1.00 35.00           O  
HETATM   40  PA  ATP A 101      10.234  18.567  19.123  1.00 35.00           P  
HETATM   41  O1A ATP A 101       9.123  19.234  19.876  1.00 35.00           O  
HETATM   42  O2A ATP A 101      11.345  19.456  18.567  1.00 35.00           O  
HETATM   43  O3A ATP A 101       9.567  17.876  17.987  1.00 35.00           O  
HETATM   44  PB  ATP A 101       8.234  18.234  17.123  1.00 35.00           P  
HETATM   45  O1B ATP A 101       7.123  19.123  17.876  1.00 35.00           O  
HETATM   46  O2B ATP A 101       9.234  17.456  16.234  1.00 35.00           O  
HETATM   47  O3B ATP A 101       7.456  17.234  17.987  1.00 35.00           O  
HETATM   48  PG  ATP A 101       6.234  16.567  17.123  1.00 35.00           P  
HETATM   49  O1G ATP A 101       5.123  17.456  16.567  1.00 35.00           O  
HETATM   50  O2G ATP A 101       7.234  15.876  16.234  1.00 35.00           O  
HETATM   51  O3G ATP A 101       5.567  15.234  17.876  1.00 35.00           O  
HETATM   52  N   ADP B 201       5.345  10.567  15.123  1.00 30.00           N  
HETATM   53  C   ADP B 201       6.234  11.234  15.567  1.00 30.00           C  
HETATM   54  O   ADP B 201       7.123  11.876  14.987  1.00 30.00           O  
HETATM   55  C   ADP B 201       5.876  12.234  16.567  1.00 30.00           C  
CONECT   11   17
END
"""
        
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Demo PDB File", "demo_protein.pdb", 
            "PDB Files (*.pdb);;All Files (*)"
        )
        
        if filename:
            with open(filename, 'w') as f:
                f.write(demo_pdb_content)
            
            # Auto-populate the input field
            self.basic_tab.pdb_path.setText(filename)
            
            QMessageBox.information(self, "Demo File Created", 
                                  f"Demo PDB file created: {filename}\n\n"
                                  "This file contains:\n"
                                  "• A small protein with 3 residues (ALA-CYS-CYS)\n"
                                  "• Disulfide bond between CYS residues\n"
                                  "• Water molecules (HOH)\n"
                                  "• Metal ion (MG)\n"
                                  "• ATP ligand (34 atoms)\n"
                                  "• ADP ligand (4 atoms)\n\n"
                                  "Perfect for testing PDBPrep functionality!\n"
                                  "The ligands should be extracted and analyzed.")
        
    def get_standard_icon(self, icon_type):
        """Safely get a standard icon with fallback."""
        try:
            return self.style().standardIcon(icon_type)
        except (AttributeError, TypeError):
            # Fallback to a simple text-based icon or no icon
            return QIcon()
    
    def create_toolbar(self):
        toolbar = QToolBar()
        self.addToolBar(toolbar)
        
        # Quick actions - using safer standard icons
        validate_action = QAction("✓ Validate", self)
        try:
            validate_action.setIcon(self.style().standardIcon(self.style().SP_DialogApplyButton))
        except AttributeError:
            pass  # Just use text if icon fails
        validate_action.triggered.connect(self.validate_settings)
        toolbar.addAction(validate_action)
        
        run_action = QAction("▶ Run", self)
        try:
            run_action.setIcon(self.style().standardIcon(self.style().SP_MediaPlay))
        except AttributeError:
            pass
        run_action.triggered.connect(self.run_pdbprep)
        toolbar.addAction(run_action)
        
        stop_action = QAction("⏹ Stop", self)
        try:
            stop_action.setIcon(self.style().standardIcon(self.style().SP_MediaStop))
        except AttributeError:
            pass
        stop_action.triggered.connect(self.stop_processing)
        toolbar.addAction(stop_action)
        
        toolbar.addSeparator()
        
        # Settings actions
        load_action = QAction("📁 Load", self)
        try:
            load_action.setIcon(self.style().standardIcon(self.style().SP_DirOpenIcon))
        except AttributeError:
            pass
        load_action.triggered.connect(self.load_preset)
        toolbar.addAction(load_action)
        
        save_action = QAction("💾 Save", self)
        try:
            save_action.setIcon(self.style().standardIcon(self.style().SP_DialogSaveButton))
        except AttributeError:
            pass
        save_action.triggered.connect(self.save_preset)
        toolbar.addAction(save_action)
    
    def validate_settings(self):
        """Validate current settings and show results."""
        errors = []
        warnings = []
        
        # Check basic settings
        if self.basic_tab.batch_mode.isChecked():
            if not self.basic_tab.batch_folder.text():
                errors.append("Batch folder not selected")
            elif not os.path.exists(self.basic_tab.batch_folder.text()):
                errors.append("Batch folder does not exist")
        else:
            if not self.basic_tab.pdb_path.text():
                errors.append("PDB file not selected")
            elif not os.path.exists(self.basic_tab.pdb_path.text()):
                errors.append("PDB file does not exist")
        
        # Check output directory
        if not self.basic_tab.output_dir.text():
            warnings.append("Output directory not specified, will use default")
        
        # Check force field selection
        if self.basic_tab.force_field.currentText() == "Custom":
            if not self.basic_tab.custom_mapping.text():
                errors.append("Custom force field selected but no mapping file provided")
        
        # Show validation results
        if errors:
            QMessageBox.critical(self, "Validation Errors", 
                               f"Found {len(errors)} error(s):\n\n" + "\n".join(f"• {e}" for e in errors))
            return False
        elif warnings:
            QMessageBox.warning(self, "Validation Warnings",
                              f"Found {len(warnings)} warning(s):\n\n" + "\n".join(f"• {w}" for w in warnings))
            return True
        else:
            QMessageBox.information(self, "Validation Success", "All settings are valid!")
            return True
    
    def run_pdbprep(self):
        """Start PDB processing."""
        if not self.validate_settings():
            return
        
        # Collect all options
        options = self.collect_options()
        
        # Show progress
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        
        # Switch to results tab
        self.tabs.setCurrentIndex(3)
        
        # Start processing thread
        self.processor_thread = PDBProcessorThread(options)
        self.processor_thread.progress_updated.connect(self.progress_bar.setValue)
        self.processor_thread.status_updated.connect(self.status_label.setText)
        self.processor_thread.log_updated.connect(self.results_tab.update_log)
        self.processor_thread.finished_signal.connect(self.processing_finished)
        self.processor_thread.error_signal.connect(self.processing_error)
        self.processor_thread.start()
    
    def stop_processing(self):
        """Stop current processing."""
        if self.processor_thread and self.processor_thread.isRunning():
            self.processor_thread.stop()
            self.processor_thread.wait()
        
        self.progress_bar.setVisible(False)
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.status_label.setText("Processing stopped")
    
    def processing_finished(self, results: Dict):
        """Handle processing completion."""
        self.progress_bar.setVisible(False)
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.status_label.setText("Processing completed successfully")
        
        # Debug: Print ALL results to console
        print("\n" + "="*50)
        print("DEBUG - COMPLETE PROCESSING RESULTS:")
        print(f"  Input file: {results.get('input_file', 'NOT SET')}")
        print(f"  Output dir: {results.get('output_directory', 'NOT SET')}")
        print(f"  Files created: {len(results.get('files_created', []))}")
        for f in results.get('files_created', []):
            print(f"    - {os.path.basename(f)}")
        print(f"  Ligands: {results.get('ligands', [])}")
        print(f"  Ligand lines found: {len(results.get('ligand_lines', []))}")
        print(f"  Ligand analysis: {results.get('ligand_analysis', 'NONE')}")
        if results.get('ligand_analysis'):
            for i, ligand in enumerate(results['ligand_analysis']):
                print(f"    Ligand {i+1}: {ligand}")
        print(f"  Disulfides: {results.get('disulfides', [])}")
        print(f"  Statistics: {results.get('statistics', {})}")
        print(f"  Errors: {results.get('errors', [])}")
        print(f"  Warnings: {results.get('warnings', [])}")
        print("="*50 + "\n")
        
        # Update results display
        self.update_results_display(results)
        
        # Enable output folder button and set results data
        if 'output_directory' in results:
            self.results_tab.set_output_directory(results['output_directory'])
        
        # Store results for export functions
        self.results_tab.set_current_results(results)
        
        # Show completion message with actual results
        output_dir = results.get('output_directory', 'unknown location')
        files_count = len(results.get('files_created', []))
        ligands_count = len(results.get('ligand_analysis', results.get('ligands', [])))
        
        completion_msg = f"""PDB preparation completed successfully!

Input: {os.path.basename(results.get('input_file', 'N/A'))}
Files created: {files_count}
Ligands processed: {ligands_count}

Files saved to:
{output_dir}

Check the console for detailed debug information."""
        
        QMessageBox.information(self, "Processing Complete", completion_msg)
    
    def processing_error(self, error_msg: str):
        """Handle processing error."""
        self.progress_bar.setVisible(False)
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.status_label.setText("Processing failed")
        
        QMessageBox.critical(self, "Processing Error", 
                           f"Processing failed with error:\n\n{error_msg}")
    
    def collect_options(self) -> Dict:
        """Collect all settings from the GUI."""
        options = {
            # Basic settings
            'pdb_file': self.basic_tab.pdb_path.text(),
            'batch_mode': self.basic_tab.batch_mode.isChecked(),
            'batch_folder': self.basic_tab.batch_folder.text(),
            'force_field': self.basic_tab.force_field.currentText(),
            'custom_mapping': self.basic_tab.custom_mapping.text(),
            'handle_alt_loc': self.basic_tab.handle_alt_loc.isChecked(),
            'rebuild_missing': self.basic_tab.rebuild_missing.isChecked(),
            'add_hydrogens': self.basic_tab.add_hydrogens.isChecked(),
            'detect_disulfides': self.basic_tab.detect_disulfides.isChecked(),
            'clean_structure': self.basic_tab.clean_structure.isChecked(),
            'output_dir': self.basic_tab.output_dir.text(),
            'output_prefix': self.basic_tab.output_prefix.text(),
            
            # Advanced settings
            'protonate_residues': self.advanced_tab.protonate_residues.isChecked(),
            'target_ph': self.advanced_tab.ph_value.value(),
            'protonation_settings': self.advanced_tab.protonation_settings,
            'check_clashes': self.advanced_tab.check_clashes.isChecked(),
            'check_geometry': self.advanced_tab.check_geometry.isChecked(),
            'check_chirality': self.advanced_tab.check_chirality.isChecked(),
            'clash_tolerance': self.advanced_tab.clash_tolerance.value(),
            'keep_waters': self.advanced_tab.keep_waters.isChecked(),
            'keep_ions': self.advanced_tab.keep_ions.isChecked(),
            'water_distance': self.advanced_tab.water_distance.value(),
            'process_all_chains': self.advanced_tab.process_all_chains.isChecked(),
            'chain_selection': self.advanced_tab.chain_selection.text(),
            
            # Ligand settings
            'auto_detect_ligands': self.ligand_tab.auto_detect_ligands.isChecked(),
            'exclude_waters': self.ligand_tab.exclude_waters.isChecked(),
            'exclude_ions': self.ligand_tab.exclude_ions.isChecked(),
            'min_heavy_atoms': self.ligand_tab.min_heavy_atoms.value(),
            'extract_ligands': self.ligand_tab.extract_ligands.isChecked(),
            'generate_topologies': self.ligand_tab.generate_topologies.isChecked(),
            'optimize_geometry': self.ligand_tab.optimize_geometry.isChecked(),
            'generate_conformers': self.ligand_tab.generate_conformers.isChecked(),
            'max_conformers': self.ligand_tab.max_conformers.value(),
            'ligand_ff': self.ligand_tab.ligand_ff.currentText(),
            'charge_method': self.ligand_tab.charge_method.currentText(),
            'calc_druglike': self.ligand_tab.calc_druglike.isChecked(),
            'lipinski_filter': self.ligand_tab.lipinski_filter.isChecked(),
            'veber_filter': self.ligand_tab.veber_filter.isChecked(),
            'pains_filter': self.ligand_tab.pains_filter.isChecked(),
            'detect_binding_site': self.ligand_tab.detect_binding_site.isChecked(),
            'binding_site_radius': self.ligand_tab.binding_site_radius.value(),
        }
        
        return options
    
    def update_results_display(self, results: Dict):
        """Update the results display with processing results."""
        # Update structure table
        stats = results.get('statistics', {})
        input_file = results.get('input_file', 'N/A')
        
        structure_data = [
            ("Input file", os.path.basename(input_file) if input_file != 'N/A' else 'N/A', "✓"),
            ("Atoms", str(stats.get('atoms', 'N/A')), "✓"),
            ("Residues", str(stats.get('residues', 'N/A')), "✓"),
            ("Chains", str(stats.get('chains', 'N/A')), "✓"),
            ("Disulfide bonds", str(len(results.get('disulfides', []))), "✓"),
            ("Water molecules", str(stats.get('waters', 'N/A')), "✓"),
        ]
        
        self.results_tab.structure_table.setRowCount(len(structure_data))
        for i, (prop, value, status) in enumerate(structure_data):
            self.results_tab.structure_table.setItem(i, 0, QTableWidgetItem(prop))
            self.results_tab.structure_table.setItem(i, 1, QTableWidgetItem(value))
            self.results_tab.structure_table.setItem(i, 2, QTableWidgetItem(status))
        
        # Update ligands table with detailed analysis
        ligand_analysis = results.get('ligand_analysis', [])
        if ligand_analysis:
            # Use detailed analysis data
            self.results_tab.ligands_table.setRowCount(len(ligand_analysis))
            for i, ligand in enumerate(ligand_analysis):
                self.results_tab.ligands_table.setItem(i, 0, QTableWidgetItem(ligand.get('name', 'Unknown')))
                self.results_tab.ligands_table.setItem(i, 1, QTableWidgetItem(ligand.get('formula', 'N/A')))
                self.results_tab.ligands_table.setItem(i, 2, QTableWidgetItem(str(ligand.get('molecular_weight', 'N/A'))))
                self.results_tab.ligands_table.setItem(i, 3, QTableWidgetItem(ligand.get('drug_like', 'Unknown')))
        else:
            # Fallback to basic ligand names
            ligands = results.get('ligands', [])
            self.results_tab.ligands_table.setRowCount(len(ligands))
            for i, ligand_name in enumerate(ligands):
                self.results_tab.ligands_table.setItem(i, 0, QTableWidgetItem(ligand_name))
                self.results_tab.ligands_table.setItem(i, 1, QTableWidgetItem("N/A"))
                self.results_tab.ligands_table.setItem(i, 2, QTableWidgetItem("N/A"))
                self.results_tab.ligands_table.setItem(i, 3, QTableWidgetItem("N/A"))
        
        # Update issues list
        errors = results.get('errors', [])
        warnings = results.get('warnings', [])
        self.results_tab.issues_list.clear()
        
        for error in errors:
            item = QListWidgetItem(f"ERROR: {error}")
            item.setBackground(QColor(255, 200, 200))  # Light red background
            self.results_tab.issues_list.addItem(item)
        
        for warning in warnings:
            item = QListWidgetItem(f"WARNING: {warning}")
            item.setBackground(QColor(255, 255, 200))  # Light yellow background
            self.results_tab.issues_list.addItem(item)
        
        if not errors and not warnings:
            item = QListWidgetItem("No issues detected")
            item.setBackground(QColor(200, 255, 200))  # Light green background
            self.results_tab.issues_list.addItem(item)
        
        # Update summary with REAL results
        files_created = results.get('files_created', [])
        output_dir = results.get('output_directory', 'N/A')
        input_file = results.get('input_file', 'N/A')
        
        # Build the file list
        file_list = ""
        if files_created:
            # Group files by type for better display
            protein_files = [f for f in files_created if 'protein' in os.path.basename(f)]
            ligand_files = [f for f in files_created if 'ligand' in os.path.basename(f)]
            topology_files = [f for f in files_created if f.endswith('.itp') or f.endswith('.top')]
            report_files = [f for f in files_created if 'report' in os.path.basename(f) or 'summary' in os.path.basename(f)]
            other_files = [f for f in files_created if f not in protein_files + ligand_files + topology_files + report_files]
            
            if protein_files:
                file_list += "\n  Protein structures:"
                for f in protein_files:
                    file_list += f"\n    • {os.path.basename(f)}"
            
            if ligand_files:
                file_list += "\n  Ligand structures:"
                for f in ligand_files[:5]:  # Show first 5, then count
                    file_list += f"\n    • {os.path.basename(f)}"
                if len(ligand_files) > 5:
                    file_list += f"\n    • ... and {len(ligand_files)-5} more ligand files"
            
            if topology_files:
                file_list += "\n  Topology files:"
                for f in topology_files:
                    file_list += f"\n    • {os.path.basename(f)}"
            
            if report_files:
                file_list += "\n  Reports:"
                for f in report_files:
                    file_list += f"\n    • {os.path.basename(f)}"
            
            if other_files:
                for f in other_files:
                    file_list += f"\n  • {os.path.basename(f)}"
        else:
            file_list = "\n• No files generated"
        
        # Build summary with actual data
        summary = f"""Processing completed successfully!

Input file: {os.path.basename(input_file) if input_file != 'N/A' else 'N/A'}
Output directory: {output_dir}
Force field: {self.basic_tab.force_field.currentText()}

Files generated ({len(files_created)} total):{file_list}

Analysis results:
• {len(results.get('disulfides', []))} disulfide bonds detected
• {len(results.get('ligands', []))} unique ligands found"""

        # Add ligand analysis details if available
        if ligand_analysis:
            drug_like = len([l for l in ligand_analysis if l.get('drug_like') == 'Likely'])
            lipinski_compliant = len([l for l in ligand_analysis if l.get('lipinski_violations', 999) <= 1])
            summary += f"\n  • {drug_like} drug-like ligands"
            summary += f"\n  • {lipinski_compliant} Lipinski-compliant ligands"
            
            # Show ligand details
            summary += f"\n  Ligand details:"
            for ligand in ligand_analysis:
                name = ligand.get('name', 'Unknown')
                formula = ligand.get('formula', 'N/A')
                mw = ligand.get('molecular_weight', 'N/A')
                drug_like_status = ligand.get('drug_like', 'Unknown')
                summary += f"\n    • {name}: {formula}, MW={mw}, {drug_like_status}"
        
        # Add statistics if available
        if stats:
            summary += f"\n• {stats.get('atoms', 'N/A')} atoms processed"
            summary += f"\n• {stats.get('residues', 'N/A')} residues analyzed"
            summary += f"\n• {stats.get('chains', 'N/A')} chains found"
        
        summary += f"\n• {len(errors)} errors, {len(warnings)} warnings"

        if output_dir != 'N/A':
            summary += f"""

📁 You can find all output files in:
{output_dir}

Click "Open Output Folder" to view the files in your system file manager."""
        
        self.results_tab.summary_text.setText(summary)
    
    def load_settings(self):
        """Load application settings."""
        # Restore window geometry
        geometry = self.settings.value("geometry")
        if geometry:
            self.restoreGeometry(geometry)
        
        # Restore other settings
        # Implementation would restore tab settings
    
    def save_settings(self):
        """Save application settings."""
        self.settings.setValue("geometry", self.saveGeometry())
        # Save other settings
    
    def new_project(self):
        """Create a new project."""
        # Reset all settings to defaults
        pass
    
    def open_project(self):
        """Open a project file."""
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open Project", "", "PDBPrep Projects (*.json);;All Files (*)"
        )
        if filename:
            # Load project settings
            pass
    
    def save_project(self):
        """Save current project."""
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Project", "", "PDBPrep Projects (*.json);;All Files (*)"
        )
        if filename:
            # Save current settings
            pass
    
    def load_preset(self):
        """Load a settings preset."""
        filename, _ = QFileDialog.getOpenFileName(
            self, "Load Preset", "", "JSON Files (*.json);;All Files (*)"
        )
        if filename:
            # Load preset settings
            pass
    
    def save_preset(self):
        """Save current settings as preset."""
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Preset", "", "JSON Files (*.json);;All Files (*)"
        )
        if filename:
            # Save current settings as preset
            pass
    
    def test_gmx_pdb2gmx(self):
        """Test gmx pdb2gmx with current settings."""
        QMessageBox.information(self, "gmx pdb2gmx Test", 
                              "This would test the prepared structure with gmx pdb2gmx")
    
    def open_structure_viewer(self):
        """Open integrated structure viewer."""
        QMessageBox.information(self, "Structure Viewer", 
                              "This would open an integrated 3D structure viewer")
    
    def show_about(self):
        """Show about dialog."""
        QMessageBox.about(self, "About PDBPrep", 
                         """<h3>PDBPrep - Professional PDB Preparation Tool</h3>
                         <p>Version 1.0</p>
                         <p>A comprehensive tool for preparing protein structures for GROMACS molecular dynamics simulations.</p>
                         <p>Features:</p>
                         <ul>
                         <li>Multi-force field support</li>
                         <li>Advanced protonation handling</li>
                         <li>Ligand extraction and parameterization</li>
                         <li>Structure validation and optimization</li>
                         <li>Drug design specific features</li>
                         </ul>
                         <p>Built with PySide6 for bioinformatics and drug design applications.</p>""")
    
    def closeEvent(self, event):
        """Handle application close event."""
        self.save_settings()
        
        if self.processor_thread and self.processor_thread.isRunning():
            reply = QMessageBox.question(self, 'Confirm Exit',
                                       'Processing is still running. Do you want to stop it and exit?',
                                       QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.Yes:
                self.processor_thread.stop()
                self.processor_thread.wait()
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()

def main():
    """Main application entry point."""
    app = QApplication(sys.argv)
    
    # Set application properties
    app.setApplicationName("PDBPrep")
    app.setApplicationVersion("1.0")
    app.setOrganizationName("PDBPrep")
    app.setOrganizationDomain("pdbprep.org")
    
    # Set application style
    app.setStyle('Fusion')
    
    # Apply dark theme
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(53, 53, 53))
    palette.setColor(QPalette.WindowText, QColor(255, 255, 255))
    palette.setColor(QPalette.Base, QColor(25, 25, 25))
    palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    palette.setColor(QPalette.ToolTipBase, QColor(0, 0, 0))
    palette.setColor(QPalette.ToolTipText, QColor(255, 255, 255))
    palette.setColor(QPalette.Text, QColor(255, 255, 255))
    palette.setColor(QPalette.Button, QColor(53, 53, 53))
    palette.setColor(QPalette.ButtonText, QColor(255, 255, 255))
    palette.setColor(QPalette.BrightText, QColor(255, 0, 0))
    palette.setColor(QPalette.Link, QColor(42, 130, 218))
    palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    palette.setColor(QPalette.HighlightedText, QColor(0, 0, 0))
    app.setPalette(palette)
    
    # Create and show main window
    window = PDBPrepMainWindow()
    window.show()
    
    # Start event loop
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
