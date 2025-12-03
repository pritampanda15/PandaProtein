# Contributing to PDBPrep

Thank you for your interest in contributing to PDBPrep! This document provides guidelines for contributing to the project.

## Code of Conduct

By participating in this project, you agree to maintain a respectful and inclusive environment for all contributors.

## How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check existing issues to avoid duplicates. When creating a bug report, include:

- **Clear title and description**
- **Steps to reproduce**
- **Expected behavior**
- **Actual behavior**
- **System information** (OS, Python version, dependency versions)
- **Sample PDB file** (if applicable and not confidential)
- **Error messages and logs**

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion, include:

- **Clear title and description**
- **Use case**: Why is this enhancement needed?
- **Proposed solution**: How should it work?
- **Alternatives considered**: What other approaches did you think about?

### Pull Requests

1. **Fork the repository** and create your branch from `main`
2. **Make your changes** following the coding standards below
3. **Add tests** if applicable
4. **Update documentation** if needed
5. **Ensure tests pass**
6. **Submit your pull request**

## Development Setup

### Prerequisites
- Python 3.8+
- GROMACS
- Git

### Setup Steps

```bash
# Clone your fork
git clone https://github.com/pritampanda15/pdbprep.git
cd pdbprep

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e .
pip install -r requirements.txt

# Install development dependencies
pip install pytest pytest-cov black flake8
```

## Coding Standards

### Python Style

- Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/)
- Use 4 spaces for indentation
- Maximum line length: 100 characters
- Use meaningful variable names

### Code Formatting

Use `black` for automatic formatting:
```bash
black pdbprep_gui.py
```

### Linting

Use `flake8` for linting:
```bash
flake8 pdbprep_gui.py --max-line-length=100
```

### Documentation

- Add docstrings to all functions and classes
- Use Google-style docstrings
- Update README.md if adding new features
- Update docs/USAGE.md for user-facing changes

Example:
```python
def process_structure(pdb_file: str, force_field: str) -> dict:
    """Process a PDB structure for MD simulation.

    Args:
        pdb_file: Path to input PDB file
        force_field: Force field to use (e.g., 'amber99sb-ildn')

    Returns:
        Dictionary containing processed file paths and metadata

    Raises:
        FileNotFoundError: If PDB file doesn't exist
        ValueError: If force field is not supported
    """
    pass
```

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=pdbprep_gui --cov-report=html
```

### Writing Tests

- Add tests for new functionality
- Ensure existing tests still pass
- Aim for high code coverage (>80%)

Example test:
```python
def test_pdb_loading():
    """Test that PDB files are loaded correctly."""
    # Test implementation
    pass
```

## Commit Messages

Write clear, meaningful commit messages:

```
Add feature: ligand conformer generation

- Implement conformer generation using RDKit
- Add user option to specify number of conformers
- Update documentation with new feature

Fixes #123
```

Format:
- Use present tense ("Add feature" not "Added feature")
- First line: brief summary (50 chars or less)
- Blank line
- Detailed description if needed
- Reference issues/PRs

## Branch Naming

Use descriptive branch names:
- `feature/add-ligand-parameterization`
- `fix/topology-generation-bug`
- `docs/update-installation-guide`
- `refactor/clean-gui-code`

## Pull Request Process

1. **Update documentation** as needed
2. **Add tests** for new functionality
3. **Ensure all tests pass**
4. **Update CHANGELOG** (if exists)
5. **Request review** from maintainers
6. **Address review comments**
7. **Squash commits** if requested

## Project Structure

```
pdbprep/
├── pdbprep_gui_latest.py    # Main application code
├── tests/                    # Test files (to be added)
├── docs/                     # Documentation
├── examples/                 # Example files (to be added)
└── scripts/                  # Utility scripts (to be added)
```

## Areas for Contribution

We especially welcome contributions in these areas:

### High Priority
- [ ] Unit tests and integration tests
- [ ] Support for additional force fields
- [ ] Command-line interface (CLI mode)
- [ ] Improved error handling and validation

### Medium Priority
- [ ] PyMOL integration for visualization
- [ ] Support for more file formats (MOL2, SDF)
- [ ] Batch processing mode
- [ ] Configuration file support

### Nice to Have
- [ ] Web interface
- [ ] Plugin system
- [ ] Cloud processing support
- [ ] Automated MD setup

## Questions?

If you have questions:
- Open an issue for discussion
- Contact the maintainers
- Check existing documentation

## Recognition

Contributors will be recognized in:
- README.md contributors section
- Release notes
- Git history

Thank you for contributing to PDBPrep!
