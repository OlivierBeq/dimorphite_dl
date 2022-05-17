Dimorphite-DL
=============

This repository is a modified version of [Durrant Lab](https://durrantlab.pitt.edu/)'s original [Dimorphite-DL](https://git.durrantlab.pitt.edu/jdurrant/dimorphite_dl) implementation ([original article](https://www.doi.org/10.1186/s13321-019-0336-9)).

Changes were brought to facilitate its installation and integration in bigger cheminformatics projects (see [CHANGES.md](CHANGES.md)).

What is it?
-----------

As stated on Durrant's Lab original GitLab page: 
<table>
  <tr>
    <td align="justify">
Dimorphite-DL adds hydrogen atoms to molecular representations, as appropriate
for a user-specified pH range. 
It is a fast, accurate, accessible, and modular open-source program for enumerating small-molecule ionization states.  
Users can provide SMILES strings from the command line or via an .smi file.
</td>
  </tr>
</table>

Citation
--------

If you use Dimorphite-DL in your research, please cite:

Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An
open-source program for enumerating the ionization states of drug-like small
molecules. J Cheminform 11:14. doi:[10.1186/s13321-019-0336-9](https://www.doi.org/10.1186/s13321-019-0336-9).

Licensing
---------

Dimorphite-DL is released under the Apache 2.0 license. See LICENCE.txt for
details.

Installation
---------

This modified version of Dimorphite-DL can be installed using [pip](https://pip.pypa.io/en/stable/getting-started/):
```bash
pip install https://github.com/OlivierBeq/dimorphite_dl/master/tarball
```
If the above does not work, you can run the installation like so:
```bash
python -m pip install https://github.com/OlivierBeq/dimorphite_dl/master/tarball
```



Usage
-----

```
usage: dimorphite [-h] [--min_ph MIN] [--max_ph MAX] [--pka_precision PRE]
                  [--smiles SMI] [--smiles_file FILE] [--output_file FILE]
                  [--max_variants MXV] [--label_states] [--silent] [--test]

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
  --test               run unit tests (for debugging)
```

The default pH range is 6.4 to 8.4, considered biologically relevant pH.

CLI usage examples
--------

```
  dimorphite_dl --smiles_file sample_molecules.smi
  dimorphite_dl --smiles "CCC(=O)O" --min_ph -3.0 --max_ph -2.0
  dimorphite_dl --smiles "CCCN" --min_ph -3.0 --max_ph -2.0 --output_file output.smi
  dimorphite_dl --smiles_file sample_molecules.smi --pka_precision 2.0 --label_states
  dimorphite_dl --test
```

Advanced Usage
--------------

It is also possible to access Dimorphite-DL from another Python script, rather
than from the command line. Here's an example:

```python
from rdkit import Chem
import dimorphite_dl

# Using the dimorphite_dl.run() function, you can run Dimorphite-DL exactly as
# you would from the command line. Here's an example:
dimorphite_dl.run(smiles="CCCN",
                  min_ph=-3.0,
                  max_ph=-2.0)


# One can also provide multiple SMILES at once.
dimorphite_dl.run(["C[C@](F)(Br)CC(O)=O", "CCCCCN"],
                  min_ph=5.0,
                  max_ph=9.0,
                  silent=True)
```

Caveats
-------

Dimorphite-DL deprotonates indoles and pyrroles around pH 14.5. But these
substructures can also be protonated around pH -3.5. Dimorphite does not
perform the protonation.

Authors and Contacts
--------------------

See the `CONTRIBUTORS.md` file for a full list of contributors. Please contact
Jacob Durrant (durrantj@pitt.edu) with any questions.
