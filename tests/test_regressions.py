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

import pytest

from dimorphite_dl.dimorphite_dl import Protonator

from .data import SPECIFIC_EXAMPLES


def test_no_carbanion():
    """A molecule that used to trigger a carbanion-producing bug."""
    smi = "Cc1nc2cc(-c3[nH]c4cc5ccccc5c5c4c3CCN(C(=O)O)[C@@H]5O)cc3c(=O)[nH][nH]c(n1)c23"
    output = Protonator(silent=True).protonate(smi)
    assert "[C-]" not in "".join(output).upper()


def test_max_variants_limit():
    """A molecule that used to produce more than max_variants outputs."""
    smi = (
        "CCCC[C@@H](C(=O)N)NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](C(C)C)"
        "NC(=O)[C@@H](NC(=O)[C@H](Cc1c[nH]c2c1cccc2)NC(=O)[C@@H](NC(=O)[C@@H]"
        "(Cc1ccc(cc1)O)N)CCC(=O)N)C)C)Cc1nc[nH]c1)Cc1ccccc1"
    )
    output = Protonator(silent=True).protonate(smi)
    assert len(output) == 128


@pytest.mark.parametrize(
    "example",
    SPECIFIC_EXAMPLES,
    ids=["ATP", "NAD"],
)
def test_specific_examples_at_given_ph(example):
    """ATP and NAD have Internal_phosphate_polyphos_chain / Initial_phosphate_like_in_ATP_ADP
    sites that can't be exercised by monoprotic test compounds, so they're checked here at
    specific, hand-picked pH values instead."""
    smi = example[0]
    for ph, expected_output in example[1:]:
        output = Protonator(min_ph=ph, max_ph=ph, pka_precision=0, silent=True).protonate(smi)
        assert output[0].strip() == expected_output, (
            f"{smi} at pH {ph} should be {expected_output}, but it is {output[0].strip()}"
        )
