#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""
List of hybrid functionals
"""

import copy

pbe0 = {
    "name": "PBE0",
    "description": '    PBE0 Hyb-GGA Exchange-Correlation Functional\n',
    "citation": '    J.P. Perdew et. al., J. Chem. Phys., 105(22), 9982-9985, 1996\n    C. Adamo et. a., J. Chem Phys., 110(13), 6158-6170, 1999\n',
    "x_functionals": {"GGA_X_PBE": {"alpha": 0.75}, "X_HF": {"alpha": 0.25}},
    "c_functionals": {"GGA_C_PBE": {}},
}

b5050lyp = {
    "name": "B5050LYP",
    "description": '    B5050LYP Hyb-GGA Exchange-Correlation Functional\n',
    "citation": '    Y. Shao et. al., J. Chem. Phys., 188, 4807-4818, 2003\n',
    "x_functionals": {"LDA_X": {"alpha": 0.08}, "GGA_X_B88": {"alpha": 0.42}, "X_HF": {"alpha": 0.50}},
    "c_functionals": {"LDA_C_VWN": {"alpha": 0.19}, "GGA_C_LYP": {"alpha": 0.81}},
}

wpbe0 = {
    "name": "WPBE0",
    "description": '    PBE0 SR-XC Functional (HJS Model)\n',
    "citation": '    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754, 2009\n',
    "x_functionals": {"GGA_X_HJS_PBE": {"omega": 0.3, "alpha": 0.75}, "X_HF": {"alpha": 0.25, "beta": 1.0, "omega": 0.3}},
    "c_functionals": {"GGA_C_PBE": {}},
}

wpbe = {
    "name": "WPBE",
    "description": '    PBE SR-XC Functional (HJS Model)\n',
    "citation": '    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754, 2009\n',
    "x_functionals": {"GGA_X_HJS_PBE": {"omega": 0.4}, "X_HF": {"alpha": 0.0, "beta": 1.0, "omega": 0.4}},
    "c_functionals": {"GGA_C_PBE": {}},
}

wb97x_d = {
    "name": "WB97X-D",
    "description": '    Parameterized Hybrid LRC B97 GGA XC Functional with Dispersion\n',
    "citation": '    J.-D. Chai and M. Head-Gordon, Phys. Chem. Chem. Phys., 10, 6615-6620, 2008\n',
    "xc_functionals": {"HYB_GGA_XC_WB97X_D": {}},
    "dispersion": {"type": "chg", "params": {"s6": 1.0}}
}

hf = {
    "name": "HF",
    "x_functionals": {"X_HF": {"alpha": 1.0}},
    "c_functionals": {},
}

hfd = {
    "name": "HF+D",
    "x_functionals": {"X_HF": {"alpha": 1.0}},
    "c_functionals": {},
    "dispersion": {"type": "das2010", "params": {"s6": 1.0}}
}


hf3c = {
    "name": "HF3C",
    "description": '    Hartree Fock as Roothaan prescribed plus 3C\n',
    "citation": '    Sure et al., J. Comput. Chem., 34, 1672-1685, 2013\n',
    "x_functionals": {"X_HF": {"alpha": 1.0}},
    "c_functionals": {},
    "dispersion": {"type": "d3bj", "params": {'s6': 1.000, 's8':  0.8777, 'a1':  0.4171, 'a2': 2.9149}},
}

pbeh3c = {
    "name": "PBEH3C",
    "description": '    PBEH-3C Hybrid GGA Exchange-Correlation Functional plus 3C\n',
    "citation": '    Grimme et. al., J. Chem. Phys., 143, 054107, 2015\n',
    "x_functionals": {"GGA_X_PBE": {"tweak": [1.0245, 0.12345679], "alpha": 0.58}, "X_HF": {"alpha": 0.42}},
    "c_functionals": {"GGA_C_PBE": {"tweak": [0.03]}},
    "dispersion": {"type": "d3bj", "params": {'s6': 1.000, 's8':  0.0000, 'a1':  0.4860, 'a2': 4.5000}},
}

sogga11_x = {
    "name": "SOGGA11-X",
    "description": '   SOGGA11-X Hybrid Exchange-Correlation Functional\n',
    "citation": '    R. Peverati and D. G. Truhlar, J. Chem. Phys. 135, 191102, 2011\n',
    "x_functionals": {"HYB_GGA_X_SOGGA11_X": {"use_libxc": True}},
    "c_functionals": {"GGA_C_SOGGA11_X": {}},
}

mn12_sx = {
    "name": "MN12_SX",
    "description": '   MN12-SX Meta-GGA Hybrid Screened Exchange-Correlation Functional\n',
    "citation": '    R. Peverati, D. G. Truhlar, Phys. Chem. Chem. Phys 14, 16187, 2012\n',
    "x_functionals": {"HYB_MGGA_X_MN12_SX": {"use_libxc": True}},
    "c_functionals": {"MGGA_C_MN12_SX": {}},
}

mn15 = {
    "name": "MN15",
    "description": '   MN15 Hybrid Meta-GGA Exchange-Correlation Functional\n',
    "citation": '    H. S. Yu, X. He, S. L. Li, and D. G. Truhlar, Chem. Sci. 7, 5032-5051, 2016\n',
    "x_functionals": {"HYB_MGGA_X_MN15": {"use_libxc": True}},
    "c_functionals": {"MGGA_C_MN15": {}},
}

beta = 0.0018903811666999256  # 5.0*(36.0*math.pi)**(-5.0/3.0)
X2S = 0.1282782438530421943003109254455883701296
X_FACTOR_C = 0.9305257363491000250020102180716672510262  #    /* 3/8*cur(3/pi)*4^(2/3) */
bt = 0.00538  # paper values
c_pw = 1.7382  # paper values
expo_pw6 = 3.8901  # paperl values
alpha_pw6 = c_pw / X2S / X2S
a_pw6 = 6.0 * bt / X2S
b_pw6 = 1.0 / X2S
c_pw6 = bt / (X_FACTOR_C * X2S * X2S)
d_pw6 = -(bt - beta) / (X_FACTOR_C * X2S * X2S)
f_pw6 = 1.0e-6 / (X_FACTOR_C * X2S**expo_pw6)
copp = 0.00262
css = 0.03668

pw6b95 = {
    "name": "PW6B95",
    "description": '    PW6B95 Hybrid Meta-GGA XC Functional\n',
    "citation": '  Y. Zhao and D. Truhlar, J. Phys. Chem. A., 109,5656-5667, 2005\n',
    "x_functionals": {"GGA_X_PW91": {"tweak": [a_pw6, b_pw6, c_pw6, d_pw6, f_pw6, alpha_pw6, expo_pw6], "alpha": 0.72}, "X_HF": {"alpha": 0.28}},
    "c_functionals": {"MGGA_C_BC95": {"tweak": [css, copp]}},
}


functional_list = {
    "TEST-PBEH3C": pbeh3c,
    "TEST-PBE0": pbe0,
    "TEST-WPBE": wpbe,
    "TEST-WPBE0": wpbe0,
    "TEST-B5050LYP": b5050lyp,
    "TEST-WB97X-D": wb97x_d,
    "TEST-HF-D": hfd,
    "TEST-HF": hf,
    "TEST-SCF": hf,
    "TEST-HF3C": hf3c,
    "TEST-SOGGA11-X": sogga11_x,
    "TEST-MN15": mn15,
    "TEST-MN12-SX": mn12_sx,
    "TEST-PW6B95": pw6b95,
    "TEST-B97-1P":          {"name": "B97-1p",          "xc_functionals": {"HYB_GGA_XC_B97_1p":          {}}},
    "TEST-B3PW91":          {"name": "B3PW91",          "xc_functionals": {"HYB_GGA_XC_B3PW91":          {}}},
    "TEST-B3LYP":           {"name": "B3LYP",           "xc_functionals": {"HYB_GGA_XC_B3LYP":           {}}},
    "TEST-B3P86":           {"name": "B3P86",           "xc_functionals": {"HYB_GGA_XC_B3P86":           {}}},
    "TEST-O3LYP":           {"name": "O3LYP",           "xc_functionals": {"HYB_GGA_XC_O3LYP":           {}}},
    "TEST-MPW1K":           {"name": "mPW1K",           "xc_functionals": {"HYB_GGA_XC_mPW1K":           {}}},
    "TEST-PBEH":            {"name": "PBEH",            "xc_functionals": {"HYB_GGA_XC_PBEH":            {}}},
    "TEST-B97":             {"name": "B97",             "xc_functionals": {"HYB_GGA_XC_B97":             {}}},
    "TEST-B97-1":           {"name": "B97-1",           "xc_functionals": {"HYB_GGA_XC_B97_1":           {}}},
    "TEST-B97-2":           {"name": "B97-2",           "xc_functionals": {"HYB_GGA_XC_B97_2":           {}}},
    "TEST-X3LYP":           {"name": "X3LYP",           "xc_functionals": {"HYB_GGA_XC_X3LYP":           {}}},
    "TEST-B1WC":            {"name": "B1WC",            "xc_functionals": {"HYB_GGA_XC_B1WC":            {}}},
    "TEST-B97-K":           {"name": "B97-K",           "xc_functionals": {"HYB_GGA_XC_B97_K":           {}}},
    "TEST-B97-3":           {"name": "B97-3",           "xc_functionals": {"HYB_GGA_XC_B97_3":           {}}},
    "TEST-MPW3PW":          {"name": "MPW3PW",          "xc_functionals": {"HYB_GGA_XC_MPW3PW":          {}}},
    "TEST-B1LYP":           {"name": "B1LYP",           "xc_functionals": {"HYB_GGA_XC_B1LYP":           {}}},
    "TEST-B1PW91":          {"name": "B1PW91",          "xc_functionals": {"HYB_GGA_XC_B1PW91":          {}}},
    "TEST-MPW1PW":          {"name": "mPW1PW",          "xc_functionals": {"HYB_GGA_XC_mPW1PW":          {}}},
    "TEST-MPW3LYP":         {"name": "MPW3LYP",         "xc_functionals": {"HYB_GGA_XC_MPW3LYP":         {}}},
    "TEST-SB98-1A":         {"name": "SB98-1a",         "xc_functionals": {"HYB_GGA_XC_SB98_1a":         {}}},
    "TEST-SB98-1B":         {"name": "SB98-1b",         "xc_functionals": {"HYB_GGA_XC_SB98_1b":         {}}},
    "TEST-SB98-1C":         {"name": "SB98-1c",         "xc_functionals": {"HYB_GGA_XC_SB98_1c":         {}}},
    "TEST-SB98-2A":         {"name": "SB98-2a",         "xc_functionals": {"HYB_GGA_XC_SB98_2a":         {}}},
    "TEST-SB98-2B":         {"name": "SB98-2b",         "xc_functionals": {"HYB_GGA_XC_SB98_2b":         {}}},
    "TEST-SB98-2C":         {"name": "SB98-2c",         "xc_functionals": {"HYB_GGA_XC_SB98_2c":         {}}},
    "TEST-HSE03":           {"name": "HSE03",           "xc_functionals": {"HYB_GGA_XC_HSE03":           {}}},
    "TEST-HSE06":           {"name": "HSE06",           "xc_functionals": {"HYB_GGA_XC_HSE06":           {}}},
    "TEST-HJS-PBE":         {"name": "HJS-PBE",         "xc_functionals": {"HYB_GGA_XC_HJS_PBE":         {}}},
    "TEST-HJS-PBE-SOL":     {"name": "HJS-PBE_SOL",     "xc_functionals": {"HYB_GGA_XC_HJS_PBE_SOL":     {}}},
    "TEST-HJS-B88":         {"name": "HJS-B88",         "xc_functionals": {"HYB_GGA_XC_HJS_B88":         {}}},
    "TEST-HJS-B97X":        {"name": "HJS-B97X",        "xc_functionals": {"HYB_GGA_XC_HJS_B97X":        {}}},
    "TEST-CAM-B3LYP":       {"name": "CAM-B3LYP",       "xc_functionals": {"HYB_GGA_XC_CAM_B3LYP":       {}}},
    "TEST-TUNED-CAM-B3LYP": {"name": "TUNED-CAM-B3LYP", "xc_functionals": {"HYB_GGA_XC_TUNED_CAM_B3LYP": {}}},
    "TEST-BHANDH":          {"name": "BHANDH",          "xc_functionals": {"HYB_GGA_XC_BHANDH":          {}}},
    "TEST-BHANDHLYP":       {"name": "BHANDHLYP",       "xc_functionals": {"HYB_GGA_XC_BHANDHLYP":       {}}},
    "TEST-MB3LYP-RC04":     {"name": "MB3LYP-RC04",     "xc_functionals": {"HYB_GGA_XC_MB3LYP_RC04":     {}}},
    "TEST-MPWLYP1M":        {"name": "MPWLYP1M",        "xc_functionals": {"HYB_GGA_XC_MPWLYP1M":        {}}},
    "TEST-REVB3LYP":        {"name": "REVB3LYP",        "xc_functionals": {"HYB_GGA_XC_REVB3LYP":        {}}},
    "TEST-CAMY-BLYP":       {"name": "CAMY-BLYP",       "xc_functionals": {"HYB_GGA_XC_CAMY_BLYP":       {}}},
    "TEST-PBE0-13":         {"name": "PBE0-13",         "xc_functionals": {"HYB_GGA_XC_PBE0_13":         {}}},
    "TEST-B3LYPS":          {"name": "B3LYPs",          "xc_functionals": {"HYB_GGA_XC_B3LYPs":          {}}},
    "TEST-WB97":            {"name": "WB97",            "xc_functionals": {"HYB_GGA_XC_WB97":            {}}},
    "TEST-WB97X":           {"name": "WB97X",           "xc_functionals": {"HYB_GGA_XC_WB97X":           {}}},
    "TEST-LRC-WPBEH":       {"name": "LRC-WPBEH",       "xc_functionals": {"HYB_GGA_XC_LRC_WPBEH":       {}}},
    "TEST-WB97X-V":         {"name": "WB97X-V",         "xc_functionals": {"HYB_GGA_XC_WB97X_V":         {}}},
    "TEST-LCY-PBE":         {"name": "LCY-PBE",         "xc_functionals": {"HYB_GGA_XC_LCY_PBE":         {}}},
    "TEST-LCY-BLYP":        {"name": "LCY-BLYP",        "xc_functionals": {"HYB_GGA_XC_LCY_BLYP":        {}}},
    "TEST-LC-VV10":         {"name": "LC-VV10",         "xc_functionals": {"HYB_GGA_XC_LC_VV10":         {}}},
    "TEST-CAMY-B3LYP":      {"name": "CAMY-B3LYP",      "xc_functionals": {"HYB_GGA_XC_CAMY_B3LYP":      {}}},
    "TEST-HPBEINT":         {"name": "HPBEINT",         "xc_functionals": {"HYB_GGA_XC_HPBEINT":         {}}},
    "TEST-LRC-WPBE":        {"name": "LRC-WPBE",        "xc_functionals": {"HYB_GGA_XC_LRC_WPBE":        {}}},
    "TEST-B3LYP5":          {"name": "B3LYP5",          "xc_functionals": {"HYB_GGA_XC_B3LYP5":          {}}},
    "TEST-EDF2":            {"name": "EDF2",            "xc_functionals": {"HYB_GGA_XC_EDF2":            {}}},
    "TEST-CAP0":            {"name": "CAP0",            "xc_functionals": {"HYB_GGA_XC_CAP0":            {}}},
    "TEST-M05":     {"name": "M05",      "xc_functionals": {"HYB_MGGA_XC_M05":      {}}},
    "TEST-M05-2X":  {"name": "M05-2X",   "xc_functionals": {"HYB_MGGA_XC_M05_2X":   {}}},
    "TEST-B88B95":  {"name": "B88B95",   "xc_functionals": {"HYB_MGGA_XC_B88B95":   {}}},
    "TEST-B86B95":  {"name": "B86B95",   "xc_functionals": {"HYB_MGGA_XC_B86B95":   {}}},
    "TEST-PW86B95": {"name": "PW86B95",  "xc_functionals": {"HYB_MGGA_XC_PW86B95":  {}}},
    "TEST-BB1K":    {"name": "BB1K",     "xc_functionals": {"HYB_MGGA_XC_BB1K":     {}}},
    "TEST-M06-HF":  {"name": "M06_HF",   "xc_functionals": {"HYB_MGGA_XC_M06_HF":   {}}},
    "TEST-MPW1B95": {"name": "MPW1B95",  "xc_functionals": {"HYB_MGGA_XC_MPW1B95":  {}}},
    "TEST-MPWB1K":  {"name": "MPWB1K",   "xc_functionals": {"HYB_MGGA_XC_MPWB1K":   {}}},
    "TEST-X1B95":   {"name": "X1B95",    "xc_functionals": {"HYB_MGGA_XC_X1B95":    {}}},
    "TEST-XB1K":    {"name": "XB1K",     "xc_functionals": {"HYB_MGGA_XC_XB1K":     {}}},
    "TEST-M06":     {"name": "M06",      "xc_functionals": {"HYB_MGGA_XC_M06":      {}}},
    "TEST-M06-2X":  {"name": "M06-2X",   "xc_functionals": {"HYB_MGGA_XC_M06_2X":   {}}},
    "TEST-PW6B95":  {"name": "PW6B95",   "xc_functionals": {"HYB_MGGA_XC_PW6B95":   {}}},
    "TEST-PWB6K":   {"name": "PWB6K",    "xc_functionals": {"HYB_MGGA_XC_PWB6K":    {}}},
    "TEST-TPSSH":   {"name": "TPSSH",    "xc_functionals": {"HYB_MGGA_XC_TPSSH":    {}}},
    "TEST-REVTPSSH":{"name": "REVTPSSH", "xc_functionals": {"HYB_MGGA_XC_REVTPSSH": {}}},
    "TEST-M08-HX":  {"name": "M08-HX",   "xc_functionals": {"HYB_MGGA_XC_M08_HX":   {}}},
    "TEST-M08-SO":  {"name": "M08-SO",   "xc_functionals": {"HYB_MGGA_XC_M08_SO":   {}}},
    "TEST-M11":     {"name": "M11",      "xc_functionals": {"HYB_MGGA_XC_M11":      {}}},
    "TEST-WB97M-V": {"name": "WB97M-V",  "xc_functionals": {"HYB_MGGA_XC_WB97M_V":  {}}},
}
