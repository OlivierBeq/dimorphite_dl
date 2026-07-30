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

from .data import DIPROTIC_GROUPS

_IDS = [group[4] for group in DIPROTIC_GROUPS]


@pytest.mark.parametrize("smi, protonated, mix, deprotonated, category", DIPROTIC_GROUPS, ids=_IDS)
def test_very_acidic(assert_protonation, smi, protonated, mix, deprotonated, category):
    args = {
        "min_ph": -10000000,
        "max_ph": -10000000,
        "pka_precision": 0.5,
        "label_states": True,
        "silent": True,
    }
    assert_protonation(args, smi, [protonated], ["PROTONATED"])


@pytest.mark.parametrize("smi, protonated, mix, deprotonated, category", DIPROTIC_GROUPS, ids=_IDS)
def test_very_basic(assert_protonation, smi, protonated, mix, deprotonated, category):
    args = {
        "min_ph": 10000000,
        "max_ph": 10000000,
        "pka_precision": 0.5,
        "label_states": True,
        "silent": True,
    }
    assert_protonation(args, smi, [deprotonated], ["DEPROTONATED"])


@pytest.mark.parametrize("smi, protonated, mix, deprotonated, category", DIPROTIC_GROUPS, ids=_IDS)
def test_at_first_pka(
    assert_protonation, average_pkas_phos, smi, protonated, mix, deprotonated, category
):
    avg_pka = average_pkas_phos[category][0]
    args = {
        "min_ph": avg_pka,
        "max_ph": avg_pka,
        "pka_precision": 0.5,
        "label_states": True,
        "silent": True,
    }
    assert_protonation(args, smi, [mix, protonated], ["BOTH"])


@pytest.mark.parametrize("smi, protonated, mix, deprotonated, category", DIPROTIC_GROUPS, ids=_IDS)
def test_at_second_pka(
    assert_protonation, average_pkas_phos, smi, protonated, mix, deprotonated, category
):
    avg_pka = average_pkas_phos[category][1]
    args = {
        "min_ph": avg_pka,
        "max_ph": avg_pka,
        "pka_precision": 0.5,
        "label_states": True,
        "silent": True,
    }
    assert_protonation(args, smi, [mix, deprotonated], ["DEPROTONATED", "DEPROTONATED"])


@pytest.mark.parametrize("smi, protonated, mix, deprotonated, category", DIPROTIC_GROUPS, ids=_IDS)
def test_at_pka_midpoint(
    assert_protonation, average_pkas_phos, smi, protonated, mix, deprotonated, category
):
    pka1, pka2 = average_pkas_phos[category]
    avg_pka = 0.5 * (pka1 + pka2)
    args = {
        "min_ph": avg_pka,
        "max_ph": avg_pka,
        "pka_precision": 5,  # wide enough to produce all three states
        "label_states": True,
        "silent": True,
    }
    assert_protonation(args, smi, [mix, deprotonated, protonated], ["BOTH", "BOTH"])
