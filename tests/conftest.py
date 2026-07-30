# Copyright 2020 Jacob D. Durrant
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
#
# This file is an adaptation of test data and logic originally written for
# src/dimorphite_dl/dimorphite_dl.py, adapted by Olivier J. M. Béquignon,
# 2026, as part of the dimorphite-ojmb fork
# (https://github.com/OlivierBeq/dimorphite_dl). See CHANGES.md for a
# summary of changes.

import functools

import pytest

from dimorphite_dl.dimorphite_dl import Protonator, ProtSubstructFuncs

from .data import DIPROTIC_GROUPS


@functools.lru_cache(maxsize=1)
def _load_average_pkas():
    cats_with_two_prot_sites = [group[4] for group in DIPROTIC_GROUPS]
    lines = ProtSubstructFuncs.load_substructre_smarts_file()
    average_pkas = {
        line.split()[0].replace("*", ""): float(line.split()[3])
        for line in lines
        if line.split()[0] not in cats_with_two_prot_sites
    }
    average_pkas_phos = {
        line.split()[0].replace("*", ""): [float(line.split()[3]), float(line.split()[6])]
        for line in lines
        if line.split()[0] in cats_with_two_prot_sites
    }
    return average_pkas, average_pkas_phos


@pytest.fixture(scope="session")
def average_pkas():
    """Average pKa per monoprotic category, loaded from site_substructures.smarts."""
    return _load_average_pkas()[0]


@pytest.fixture(scope="session")
def average_pkas_phos():
    """[pKa1, pKa2] per diprotic (phosphate/phosphonate) category."""
    return _load_average_pkas()[1]


@pytest.fixture
def assert_protonation():
    """Factory fixture: run a Protonator and assert its output matches expectations."""

    def _assert_protonation(args, smiles, expected_output, labels):
        output = list(Protonator(**args).protonate(smiles))
        output = [
            [value] + output[1] if isinstance(output[1], list) else [value, output[1]]
            for value in output[0]
        ]

        num_states = len(expected_output)
        assert len(output) == num_states, (
            f"{smiles} should have {num_states} states at pH {args['min_ph']}: {output}"
        )

        got_smiles = {entry[0] for entry in output}
        assert not got_smiles - set(expected_output), (
            f"{smiles} is not {' AND '.join(expected_output)} at pH "
            f"{args['min_ph']} - {args['max_ph']}; it is {' AND '.join(got_smiles)}"
        )

        got_labels = {entry[1] for entry in output}
        assert not got_labels - set(labels), (
            f"{smiles} not labeled as {' AND '.join(labels)}; it is {' AND '.join(got_labels)}"
        )

    return _assert_protonation
