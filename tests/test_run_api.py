# Copyright 2026 Olivier J. M. Béquignon
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest
from rdkit import Chem

from dimorphite_dl.dimorphite_dl import run


def test_run_with_none():
    assert run(None) == [None]


def test_run_rejects_non_str_non_list():
    with pytest.raises(ValueError, match="smiles must be a string"):
        run(123)


def test_run_rejects_bare_mol_object():
    # A single (non-list) Chem.Mol is not accepted, despite the type hint
    # suggesting it should be -- this documents existing behavior.
    mol = Chem.MolFromSmiles("CCCN")
    with pytest.raises(ValueError, match="smiles must be a string"):
        run(mol, silent=True)


def test_run_single_smiles_string():
    assert run("CCC(=O)O", min_ph=-3.0, max_ph=-2.0, silent=True) == ["CCC(=O)O"]


def test_run_list_of_smiles_strings():
    result = run(["CCC(=O)O", "CCCN"], min_ph=-3.0, max_ph=-2.0, silent=True)
    assert result == ["CCC(=O)O", "CCC[NH3+]"]


def test_run_list_of_mol_objects():
    mols = [Chem.MolFromSmiles("CCCN")]
    result = run(mols, min_ph=-3.0, max_ph=-2.0, silent=True)
    assert result == ["CCC[NH3+]"]


def test_run_malformed_smiles_yields_none_without_raising():
    result = run(["CCCN", "not_a_smiles"], min_ph=-3.0, max_ph=-2.0, silent=True)
    assert result == ["CCC[NH3+]", None]


def test_run_molecule_with_no_protonation_sites():
    result = run("C", label_states=True, silent=True)
    assert result == [(["C"], None)]


def test_run_max_variants_limit_warning(capsys):
    smi = (
        "CCCC[C@@H](C(=O)N)NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](C(C)C)"
        "NC(=O)[C@@H](NC(=O)[C@H](Cc1c[nH]c2c1cccc2)NC(=O)[C@@H](NC(=O)[C@@H]"
        "(Cc1ccc(cc1)O)N)CCC(=O)N)C)C)Cc1nc[nH]c1)Cc1ccccc1"
    )
    result = run(smi, max_variants=2, silent=False)
    assert len(result[0]) == 2
    assert "Limited number of variants" in capsys.readouterr().out


def test_run_prints_citation_header_when_not_silent(capsys):
    run("CCCN", min_ph=-3.0, max_ph=-2.0, silent=False)
    assert "Dimorphite-DL" in capsys.readouterr().out
