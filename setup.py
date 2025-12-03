#!/usr/bin/env python3
"""Setup script for PDBPrep."""

from setuptools import setup, find_packages
from pathlib import Path

# Read the contents of README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

# Read requirements
requirements = []
with open('requirements.txt') as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith('#'):
            requirements.append(line)

setup(
    name='pdbprep',
    version='1.0.0',
    description='Professional PDB Preparation Tool for GROMACS',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='PDBPrep Contributors',
    author_email='your.email@example.com',
    url='https://github.com/pritampanda15/pdbprep',
    license='MIT',
    py_modules=['pdbprep_gui_latest'],
    install_requires=requirements,
    python_requires='>=3.8',
    entry_points={
        'console_scripts': [
            'pdbprep=pdbprep_gui:main',
        ],
        'gui_scripts': [
            'pdbprep-gui=pdbprep_gui:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Operating System :: OS Independent',
        'Environment :: X11 Applications :: Qt',
    ],
    keywords='pdb gromacs molecular-dynamics drug-design protein-preparation',
    project_urls={
        'Bug Reports': 'https://github.com/pritampanda15/pdbprep/issues',
        'Source': 'https://github.com/pritampanda15/pdbprep',
        'Documentation': 'https://github.com/pritampanda15/pdbprep/blob/main/docs/USAGE.md',
    },
)
