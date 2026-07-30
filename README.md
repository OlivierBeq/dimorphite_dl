# ⚗️ Dimorphite-DL

<!-- Badges -->
<div align="center">

[![PyPI version](https://img.shields.io/pypi/v/dimorphite-ojmb.svg)](https://pypi.org/project/dimorphite-ojmb/)
[![Supported Python versions](https://img.shields.io/pypi/pyversions/dimorphite-ojmb.svg)](https://pypi.org/project/dimorphite-ojmb/)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Tests](https://github.com/OlivierBeq/dimorphite_dl/actions/workflows/ci.yml/badge.svg)](https://github.com/OlivierBeq/dimorphite_dl/actions/workflows/ci.yml)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
<br>
</div>

A fast, accurate, and modular open-source program for enumerating the ionization states of drug-like small molecules over a user-specified pH range. 

⚠️ This repository is a modernized fork of [Durrant Lab](https://durrantlab.pitt.edu/)'s original [Dimorphite-DL](https://git.durrantlab.pitt.edu/jdurrant/dimorphite_dl) ([original article](https://www.doi.org/10.1186/s13321-019-0336-9)), reworked for easier installation and integration into bigger cheminformatics projects — see [CHANGES.md](CHANGES.md) for the full history.

## ✨ Features

- ⚡ **Fast pH-range enumeration** — adds/removes hydrogens for any pH window you specify, with a tunable `--pka_precision` (number of standard deviations around each predicted pKa).
- 🧬 **Broad functional-group coverage** — alcohols, amides, amines/anilines, azides, carboxyls, guanidines, imides, nitro groups, phosphates/phosphonates (including ATP/NAD-style polyphosphate chains), sulfates/sulfonates/sulfonamides, thiols, indoles/pyrroles, aromatic nitrogens, and more.
- 🖥️ **CLI and 🐍 Python API** — run it as a `dimorphite` command, or call `dimorphite_dl.run(...)` directly from your own scripts.
- 📄 **Flexible input** — SMILES from the command line, an `.smi` file, or a Python list of SMILES strings / RDKit `Mol` objects.
- 🏷️ **Optional state labeling** — tag each output SMILES as `PROTONATED`, `DEPROTONATED`, or `BOTH`.
- 🛡️ **Bounded output** — `--max_variants` caps the number of protonation variants generated per input compound.
- ✅ **Modern tooling** — `pyproject.toml` packaging, a `pytest` test suite, `ruff` linting/formatting, and CI across Linux/macOS/Windows × Python 3.11-3.13.

## ✍️ Copyright and Citation Notice

Olivier J. M. Béquignon is **neither** the copyright holder of the original Dimorphite-DL algorithm **nor** responsible for it. This repository packages and modernizes Jacob D. Durrant's original work (Apache-2.0 licensed) — see [CONTRIBUTORS.md](CONTRIBUTORS.md) for the original authors.

### Citing

If you use Dimorphite-DL in your research, please cite the original publication:

> Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An open-source program for enumerating the ionization states of drug-like small molecules. *J Cheminform* 11:14.
> [DOI: 10.1186/s13321-019-0336-9](https://www.doi.org/10.1186/s13321-019-0336-9)

## 📦 Installation

```bash
pip install dimorphite-ojmb
```

Or from source:

```bash
git clone https://github.com/OlivierBeq/dimorphite_dl.git
pip install ./dimorphite_dl
```

## 🛠️ Requirements

- Python 3.11+
- [RDKit](https://www.rdkit.org/docs/Install.html)
- [Click](https://click.palletsprojects.com/) (installed automatically)

## 💡 Usage

### Command line

```
usage: dimorphite [-h] [--min_ph MIN] [--max_ph MAX] [--pka_precision PRE]
                  [--smiles SMI] [--smiles_file FILE] [--output_file FILE]
                  [--max_variants MXV] [--label_states] [--silent]

Dimorphite 1.2.4: Creates models of appropriately protonated small moleucles.
Apache 2.0 License. Copyright 2020 Jacob D. Durrant.

Options:
  -h, --help           Show this message and exit.
  --min_ph MIN         minimum pH to consider (default: 6.4)
  --max_ph MAX         maximum pH to consider (default: 8.4)
  --pka_precision PRE  pKa precision factor (number of standard devations,
                       default: 1.0)
  --smiles SMI         SMILES string to protonate NOTE: This argument is
                       mutually exclusive with smiles_file.  [required]
  --smiles_file FILE   file that contains SMILES strings to protonate NOTE:
                       This argument is mutually exclusive with smiles.
                       [required]
  --output_file FILE   output file to write protonated SMILES (optional)
  --max_variants MXV   limit number of variants per input compound (default:
                       128)
  --label_states       label protonated SMILES with target state (i.e.,
                       "DEPROTONATED", "PROTONATED", or "BOTH").
  --silent             do not print any messages to the screen
```

The default pH range is 6.4 to 8.4, considered biologically relevant pH.

```bash
dimorphite --smiles_file sample_molecules.smi
dimorphite --smiles "CCC(=O)O" --min_ph -3.0 --max_ph -2.0
dimorphite --smiles "CCCN" --min_ph -3.0 --max_ph -2.0 --output_file output.smi
dimorphite --smiles_file sample_molecules.smi --pka_precision 2.0 --label_states
```

### Python API

```python
import dimorphite_dl

# Using the dimorphite_dl.run() function, you can run Dimorphite-DL exactly as
# you would from the command line. Here's an example:
dimorphite_dl.run(smiles="CCCN", min_ph=-3.0, max_ph=-2.0)

# One can also provide multiple SMILES at once.
dimorphite_dl.run(["C[C@](F)(Br)CC(O)=O", "CCCCCN"], min_ph=5.0, max_ph=9.0, silent=True)
```

RDKit `Mol` objects are accepted too, either on their own or in a list:

```python
from rdkit import Chem

mol = Chem.MolFromSmiles("CCCN")
dimorphite_dl.run(mol, min_ph=-3.0, max_ph=-2.0, silent=True)
dimorphite_dl.run([mol], min_ph=-3.0, max_ph=-2.0, silent=True)
```

### Running the tests

From a source checkout:

```bash
pip install -e ".[test]"
pytest
```

## ⚠️ Caveats

Dimorphite-DL deprotonates indoles and pyrroles around pH 14.5. But these substructures can also be protonated around pH -3.5. Dimorphite does not perform the protonation.

## 📄 License

This project is licensed under the Apache 2.0 License - see the [LICENSE.txt](https://github.com/OlivierBeq/dimorphite_dl/blob/master/LICENSE.txt) file for details.

## 📚 API Documentation

```python
def run(smiles, min_ph=6.4, max_ph=8.4, pka_precision=1.0, max_variants=128, label_states=False, silent=False):
```

Runs Dimorphite-DL protonation on one or more molecules.

#### Parameters

- ***smiles  : str | Chem.Mol | list[str] | list[Chem.Mol]***
  A SMILES string, an RDKit `Mol` object, or a list of either, to protonate. `None` returns `[None]`.
- ***min_ph  : float***
  Minimum pH to consider (default: `6.4`).
- ***max_ph  : float***
  Maximum pH to consider (default: `8.4`).
- ***pka_precision  : float***
  pKa precision factor — number of standard deviations around each predicted pKa to consider (default: `1.0`).
- ***max_variants  : int***
  Limit on the number of protonation variants generated per input compound (default: `128`).
- ***label_states  : bool***
  If `True`, each result becomes a `(smiles_list, state)` tuple, where `state` is `"PROTONATED"`, `"DEPROTONATED"`, `"BOTH"`, a list thereof (one per site), or `None` if the molecule has no ionizable site.
- ***silent  : bool***
  Suppress the citation banner and warning messages.
- ***return_type  : list***
  One entry per input molecule, in the same order: a SMILES string (or `(smiles_list, state)` tuple if `label_states=True`), or `None` if the molecule could not be parsed.

## 👥 Authors and Contacts

See the [CONTRIBUTORS.md](CONTRIBUTORS.md) file for a full list of original contributors. Please contact Jacob Durrant (durrantj@pitt.edu) with any questions about the original algorithm, or open an [issue](https://github.com/OlivierBeq/dimorphite_dl/issues) for questions about this fork.
