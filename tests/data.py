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

"""Protonation test data, ported verbatim from the legacy TestFuncs test suite."""

# fmt: off
MONOPROTIC_GROUPS = [
    # input smiles,            protonated,                  deprotonated,               category
    ["C#CCO",                  "C#CCO",                     "C#CC[O-]",                 "Alcohol"],
    ["C(=O)N",                 "NC=O",                      "[NH-]C=O",                 "Amide"],
    ["CC(=O)NOC(C)=O",         "CC(=O)NOC(C)=O",            "CC(=O)[N-]OC(C)=O",        "Amide_electronegative"],
    ["COC(=N)N",               "COC(N)=[NH2+]",             "COC(=N)N",                 "AmidineGuanidine2"],
    ["Brc1ccc(C2NCCS2)cc1",    "Brc1ccc(C2[NH2+]CCS2)cc1",  "Brc1ccc(C2NCCS2)cc1",      "Amines_primary_secondary_tertiary"],
    ["CC(=O)[n+]1ccc(N)cc1",   "CC(=O)[n+]1ccc([NH3+])cc1", "CC(=O)[n+]1ccc(N)cc1",     "Anilines_primary"],
    ["CCNc1ccccc1",            "CC[NH2+]c1ccccc1",          "CCNc1ccccc1",              "Anilines_secondary"],
    ["Cc1ccccc1N(C)C",         "Cc1ccccc1[NH+](C)C",        "Cc1ccccc1N(C)C",           "Anilines_tertiary"],
    ["BrC1=CC2=C(C=C1)NC=C2",  "Brc1ccc2[nH]ccc2c1",        "Brc1ccc2[n-]ccc2c1",       "Indole_pyrrole"],
    ["O=c1cc[nH]cc1",          "O=c1cc[nH]cc1",             "O=c1cc[n-]cc1",            "Aromatic_nitrogen_protonated"],
    ["C-N=[N+]=[N@H]",         "CN=[N+]=N",                 "CN=[N+]=[N-]",             "Azide"],
    ["BrC(C(O)=O)CBr",         "O=C(O)C(Br)CBr",            "O=C([O-])C(Br)CBr",        "Carboxyl"],
    ["NC(NN=O)=N",             "NC(=[NH2+])NN=O",           "N=C(N)NN=O",               "AmidineGuanidine1"],
    ["C(F)(F)(F)C(=O)NC(=O)C", "CC(=O)NC(=O)C(F)(F)F",      "CC(=O)[N-]C(=O)C(F)(F)F",  "Imide"],
    ["O=C(C)NC(C)=O",          "CC(=O)NC(C)=O",             "CC(=O)[N-]C(C)=O",         "Imide2"],
    ["CC(C)(C)C(N(C)O)=O",     "CN(O)C(=O)C(C)(C)C",        "CN([O-])C(=O)C(C)(C)C",    "N-hydroxyamide"],
    ["C[N+](O)=O",             "C[N+](=O)O",                "C[N+](=O)[O-]",            "Nitro"],
    ["O=C1C=C(O)CC1",          "O=C1C=C(O)CC1",             "O=C1C=C([O-])CC1",         "O=C-C=C-OH"],
    ["C1CC1OO",                "OOC1CC1",                   "[O-]OC1CC1",               "Peroxide2"],
    ["C(=O)OO",                "O=COO",                     "O=CO[O-]",                 "Peroxide1"],
    ["Brc1cc(O)cc(Br)c1",      "Oc1cc(Br)cc(Br)c1",         "[O-]c1cc(Br)cc(Br)c1",     "Phenol"],
    ["CC(=O)c1ccc(S)cc1",      "CC(=O)c1ccc(S)cc1",         "CC(=O)c1ccc([S-])cc1",     "Phenyl_Thiol"],
    ["C=CCOc1ccc(C(=O)O)cc1",  "C=CCOc1ccc(C(=O)O)cc1",     "C=CCOc1ccc(C(=O)[O-])cc1", "Phenyl_carboxyl"],
    ["COP(=O)(O)OC",           "COP(=O)(O)OC",              "COP(=O)([O-])OC",          "Phosphate_diester"],
    ["CP(C)(=O)O",             "CP(C)(=O)O",                "CP(C)(=O)[O-]",            "Phosphinic_acid"],
    ["CC(C)OP(C)(=O)O",        "CC(C)OP(C)(=O)O",           "CC(C)OP(C)(=O)[O-]",       "Phosphonate_ester"],
    ["CC1(C)OC(=O)NC1=O",      "CC1(C)OC(=O)NC1=O",         "CC1(C)OC(=O)[N-]C1=O",     "Ringed_imide1"],
    ["O=C(N1)C=CC1=O",         "O=C1C=CC(=O)N1",            "O=C1C=CC(=O)[N-]1",        "Ringed_imide2"],
    ["O=S(OC)(O)=O",           "COS(=O)(=O)O",              "COS(=O)(=O)[O-]",          "Sulfate"],
    ["COc1ccc(S(=O)O)cc1",     "COc1ccc(S(=O)O)cc1",        "COc1ccc(S(=O)[O-])cc1",    "Sulfinic_acid"],
    ["CS(N)(=O)=O",            "CS(N)(=O)=O",               "CS([NH-])(=O)=O",          "Sulfonamide"],
    ["CC(=O)CSCCS(O)(=O)=O",   "CC(=O)CSCCS(=O)(=O)O",      "CC(=O)CSCCS(=O)(=O)[O-]",  "Sulfonate"],
    ["CC(=O)S",                "CC(=O)S",                   "CC(=O)[S-]",               "Thioic_acid"],
    ["C(C)(C)(C)(S)",          "CC(C)(C)S",                 "CC(C)(C)[S-]",             "Thiol"],
    ["Brc1cc[nH+]cc1",         "Brc1cc[nH+]cc1",            "Brc1ccncc1",               "Aromatic_nitrogen_unprotonated"],
    ["C=C(O)c1c(C)cc(C)cc1C",  "C=C(O)c1c(C)cc(C)cc1C",     "C=C([O-])c1c(C)cc(C)cc1C", "Vinyl_alcohol"],
    ["CC(=O)ON",               "CC(=O)O[NH3+]",             "CC(=O)ON",                 "Primary_hydroxyl_amine"],
    # Note: Internal_phosphate_polyphos_chain and Initial_phosphate_like_in_ATP_ADP
    # are not tested here because there is no way to generate monoprotic compounds
    # for them. See SPECIFIC_EXAMPLES (ATP/NAD) instead.
]

DIPROTIC_GROUPS = [
    # input smiles,      protonated,       mix (one deprotonation), deprotonated,           category
    ["O=P(O)(O)OCCCC", "CCCCOP(=O)(O)O", "CCCCOP(=O)([O-])O", "CCCCOP(=O)([O-])[O-]", "Phosphate"],
    ["CC(P(O)(O)=O)C", "CC(C)P(=O)(O)O", "CC(C)P(=O)([O-])O", "CC(C)P(=O)([O-])[O-]", "Phosphonate"],
]
# fmt: on

# Each entry: [input smiles, *(pH, expected protonated SMILES at that pH)]
SPECIFIC_EXAMPLES = [
    [
        "O=P(O)(OP(O)(OP(O)(OCC1OC(C(C1O)O)N2C=NC3=C2N=CN=C3N)=O)=O)O",  # ATP
        (
            0.5,
            "[NH3+]c1[nH+]c[nH+]c2c1[nH+]cn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O",
        ),
        (
            1.0,
            "[NH3+]c1[nH+]c[nH+]c2c1[nH+]cn2C1OC(COP(=O)(O)OP(=O)([O-])OP(=O)(O)O)C(O)C1O",
        ),
        (
            2.6,
            "[NH3+]c1[nH+]c[nH+]c2c1[nH+]cn2C1OC(COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])O)C(O)C1O",
        ),
        (
            7.0,
            "Nc1ncnc2c1ncn2C1OC(COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])C(O)C1O",
        ),
    ],
    [
        "O=P(O)(OP(O)(OCC1C(O)C(O)C(N2C=NC3=C(N)N=CN=C32)O1)=O)OCC(O4)C(O)C(O)C4[N+]5=CC=CC(C(N)=O)=C5",  # NAD
        (
            0.5,
            "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c([NH3+])ncnc54)C(O)C3O)C(O)C2O)c1",
        ),
        (
            2.5,
            "NC(=O)c1ccc[n+](C2OC(COP(=O)([O-])OP(=O)([O-])OCC3OC(n4cnc5c([NH3+])ncnc54)C(O)C3O)C(O)C2O)c1",
        ),
        (
            7.4,
            "NC(=O)c1ccc[n+](C2OC(COP(=O)([O-])OP(=O)([O-])OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1",
        ),
    ],
]
