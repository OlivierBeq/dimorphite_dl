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

from click.testing import CliRunner

from dimorphite_dl.cli import main


def run_cli(args):
    return CliRunner().invoke(main, args)


def test_help():
    result = run_cli(["--help"])
    assert result.exit_code == 0
    assert "Options:" in result.output


def test_missing_required_option():
    result = run_cli([])
    assert result.exit_code != 0
    assert "Missing option" in result.output


def test_mutually_exclusive_options():
    result = run_cli(["--smiles", "CCCN", "--smiles_file", "unused.smi"])
    assert result.exit_code != 0
    assert "mutually exclusive" in result.output


def test_single_smiles_single_variant():
    result = run_cli(["--smiles", "CCC(=O)O", "--min_ph", "-3.0", "--max_ph", "-2.0", "--silent"])
    assert result.exit_code == 0
    assert result.output.strip() == "CCC(=O)O"


def test_single_smiles_multiple_variants():
    result = run_cli(
        [
            "--smiles",
            "C#CCO",
            "--min_ph",
            "14.78",
            "--max_ph",
            "14.78",
            "--pka_precision",
            "0.5",
            "--silent",
        ]
    )
    assert result.exit_code == 0
    lines = result.output.strip().splitlines()
    assert set(lines) == {"C#CC[O-]", "C#CCO"}


def test_multiple_smiles_options():
    result = run_cli(
        [
            "--smiles",
            "CCC(=O)O",
            "--smiles",
            "CCCN",
            "--min_ph",
            "-3.0",
            "--max_ph",
            "-2.0",
            "--silent",
        ]
    )
    assert result.exit_code == 0
    lines = result.output.strip().splitlines()
    assert set(lines) == {"CCC(=O)O", "CCC[NH3+]"}


def test_label_states_single_label():
    result = run_cli(
        [
            "--smiles",
            "C#CCO",
            "--min_ph",
            "14.78",
            "--max_ph",
            "14.78",
            "--pka_precision",
            "0.5",
            "--label_states",
            "--silent",
        ]
    )
    assert result.exit_code == 0
    lines = result.output.strip().splitlines()
    assert set(lines) == {"C#CC[O-]\tBOTH", "C#CCO\tBOTH"}


def test_label_states_multi_site_label_list():
    result = run_cli(
        [
            "--smiles",
            "O=P(O)(O)OCCCC",
            "--min_ph",
            "2.0",
            "--max_ph",
            "9.5",
            "--pka_precision",
            "5",
            "--label_states",
            "--silent",
        ]
    )
    assert result.exit_code == 0
    assert "BOTH\tBOTH" in result.output


def test_smiles_file_with_data_column(tmp_path):
    smi_file = tmp_path / "input.smi"
    smi_file.write_text("CCC(=O)O mycompound\n")
    result = run_cli(
        ["--smiles_file", str(smi_file), "--min_ph", "-3.0", "--max_ph", "-2.0", "--silent"]
    )
    assert result.exit_code == 0
    assert result.output.strip() == "CCC(=O)O\tmycompound"


def test_smiles_file_with_data_column_and_label_states(tmp_path):
    smi_file = tmp_path / "input.smi"
    smi_file.write_text("CCC(=O)O mycompound\n")
    result = run_cli(
        [
            "--smiles_file",
            str(smi_file),
            "--min_ph",
            "-3.0",
            "--max_ph",
            "-2.0",
            "--label_states",
            "--silent",
        ]
    )
    assert result.exit_code == 0
    assert result.output.strip() == "CCC(=O)O\tmycompound\tPROTONATED"


def test_output_file(tmp_path):
    out_file = tmp_path / "output.smi"
    result = run_cli(
        [
            "--smiles",
            "CCC(=O)O",
            "--min_ph",
            "-3.0",
            "--max_ph",
            "-2.0",
            "--silent",
            "--output_file",
            str(out_file),
        ]
    )
    assert result.exit_code == 0
    assert result.output == ""
    assert out_file.read_text().strip() == "CCC(=O)O"


def test_output_file_with_label_states(tmp_path):
    out_file = tmp_path / "output.smi"
    result = run_cli(
        [
            "--smiles",
            "CCC(=O)O",
            "--min_ph",
            "-3.0",
            "--max_ph",
            "-2.0",
            "--label_states",
            "--silent",
            "--output_file",
            str(out_file),
        ]
    )
    assert result.exit_code == 0
    assert result.output == ""
    assert out_file.read_text().strip() == "CCC(=O)O\tPROTONATED"


def test_malformed_smiles_is_silently_skipped():
    result = run_cli(["--smiles", "not_a_smiles", "--silent"])
    assert result.exit_code == 0
    assert result.output == ""


def test_malformed_smiles_skipped_with_label_states():
    result = run_cli(["--smiles", "not_a_smiles", "--label_states", "--silent"])
    assert result.exit_code == 0
    assert result.output == ""
