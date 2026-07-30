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

from io import StringIO

import pytest

from dimorphite_dl.dimorphite_dl import LoadSMIFile, Protonator, UtilFuncs


def test_convert_smiles_str_to_mol_rejects_none():
    assert UtilFuncs.convert_smiles_str_to_mol(None) is None


def test_convert_smiles_str_to_mol_rejects_non_str():
    assert UtilFuncs.convert_smiles_str_to_mol(123) is None


def test_load_smi_file_skips_blank_lines():
    loader = LoadSMIFile(StringIO("\nCCCN\n"), silent=True)
    entries = list(loader)
    assert [e["smiles"] for e in entries] == ["CCCN"]


def test_load_smi_file_warns_on_malformed_line(capsys):
    loader = LoadSMIFile(StringIO("not_a_smiles\nCCCN\n"), silent=False)
    entries = list(loader)
    assert [e["smiles"] for e in entries] == ["CCCN"]
    assert "WARNING" in capsys.readouterr().out


def test_protonate_rejects_non_str_smiles():
    with pytest.raises(TypeError, match="only SMILES can be prepared"):
        Protonator(silent=True).protonate(123)
