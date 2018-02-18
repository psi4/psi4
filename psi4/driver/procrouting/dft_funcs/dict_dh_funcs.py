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
List of XC functionals
"""


b2plyp = {
    "name": "B2PLYP",
    "x_functionals": {"GGA_X_B88": {"alpha": 0.47}, "X_HF": {"alpha": 0.53}},
    "c_functionals": {"GGA_C_LYP": {"alpha": 0.73}, "C_MP2": {"alpha": 0.27}},
    "citation": '    S. Grimme, J. Chem. Phys., 124, 034108, 2006\n',
    "description": '    B2PLYP Double Hybrid Exchange-Correlation Functional\n',
}

dsd_blyp = {
    "name": "DSD-BLYP",
    "x_functionals": {"GGA_X_B88": {"alpha": 0.30}, "X_HF": {"alpha": 0.70}},
    "c_functionals": {"GGA_C_LYP": {"alpha": 0.56}, "C_MP2": {"alpha": 1.0, "os": 0.46, "ss": 0.40}},
    "citation": '    S. Kozuch, Phys. Chem. Chem. Phys., 13, 20104, 2011\n',
    "description": '    DSD-BLYP Dispersion-corrected SCS Double Hybrid XC Functional\n',
}

core_dsd_blyp = {
    "name": "CORE-DSD-BLYP",
    "x_functionals": {"GGA_X_B88": {"alpha": 0.31}, "X_HF": {"alpha": 0.69}},
    "c_functionals": {"GGA_C_LYP": {"alpha": 0.54}, "C_MP2": {"alpha": 1.0, "os": 0.46, "ss": 0.37}},
    "citation": '    S. Kozuch, Phys. Chem. Chem. Phys., 13, 20104, 2011\n',
    "description": '    DSD-BLYP Dispersion-corrected SCS Double Hybrid XC Functional\n' + \
                   '    (full-core parameterization) \n'  
}

pbe0_2 = {
    "name": "PBE0-2",
    "x_functionals": {"GGA_X_PBE": {"alpha": 0.206299}, "X_HF": {"alpha": 0.793701}},
    "c_functionals": {"GGA_C_PBE": {"alpha": 0.5}, "C_MP2": {"alpha": 0.5}},
    "citation": '    J. Chai, Chem. Phys. Lett., 538, 121-125, 2012\n',
    "description": '    PBE0-2 Double Hybrid Exchange-Correlation Functional\n',
}

pbe0_dh = {
    "name": "PBE0-DH",
    "x_functionals": {"GGA_X_PBE": {"alpha": 0.5}, "X_HF": {"alpha": 0.5}},
    "c_functionals": {"GGA_C_PBE": {"alpha": 0.875}, "C_MP2": {"alpha": 0.125}},
    "citation": '    E. Bremond, C. Adamo, J. Chem. Phys., 135, 024106, 2011\n',
    "description": '    PBE0-DH Double Hybrid Exchange-Correlation Functional\n',
}

dsd_pbep86 = {
    "name": "DSD-PBEP86",
    "x_functionals": {"GGA_X_PBE": {"alpha": 0.32}, "X_HF": {"alpha": 0.68}},
    "c_functionals": {"GGA_C_P86": {"alpha": 0.45}, "C_MP2": {"alpha": 1.0, "os": 0.51, "ss": 0.23}},
    "citation": '    S. Kozuch, Phys. Chem. Chem. Phys., 13, 20104, 2011\n',
    "description": '    DSD-PBEP86 Dispersion-corrected SCS Double Hybrid XC Functional (opt. for -D2)\n',
}

dsd_pbepbe = {
    "name": "DSD-PBEPBE",
    "x_functionals": {"GGA_X_PBE": {"alpha": 0.34}, "X_HF": {"alpha": 0.66}},
    "c_functionals": {"GGA_C_PBE": {"alpha": 0.51}, "C_MP2": {"alpha": 1.0, "os": 0.53, "ss": 0.12}},
    "citation": '    S. Kozuch, Phys. Chem. Chem. Phys., 13, 20104, 2011\n',
    "description": '    DSD-PBEPBE Dispersion-corrected SCS Double Hybrid XC Functional\n',
}

b2gpplyp = {
    "name": "B2GPPLYP",
    "x_functionals": {"GGA_X_B88": {"alpha": 0.35}, "X_HF": {"alpha": 0.65}},
    "c_functionals": {"GGA_C_LYP": {"alpha": 0.64}, "C_MP2": {"alpha": 0.36}},
    "citation": '    A. Karton, A. Tarnopolsky, J.-F. Lamere, G. C. Schatz, J.M. L. Martin, J.Phys. Chem. A, 112, 12868-12886,2008   \n',
    "description": '    B2GPPLYP Double Hybrid Exchange-Correlation Functional\n',
}

beta = 5.0 * (36.0 * 3.141592653589793)**(-5.0 / 3.0)
X2S = 0.1282782438530421943003109254455883701296
X_FACTOR_C = 0.9305257363491000250020102180716672510262  #    /* 3/8*cur(3/pi)*4^(2/3) */
bt = 0.004440  # paper values
c_pw = 0.32620  # paper values
expo_pw6 = 3.7868  # paper values
alpha_pw6 = c_pw / X2S / X2S
a_pw6 = 6.0 * bt / X2S
b_pw6 = 1.0 / X2S
c_pw6 = bt / (X_FACTOR_C * X2S * X2S)
d_pw6 = -(bt - beta) / (X_FACTOR_C * X2S * X2S)
f_pw6 = 1.0e-6 / (X_FACTOR_C * X2S**expo_pw6)
copp = 0.00250
css = 0.03241

pwpb95 = {
    "name": "PWPB95",
    "x_functionals": {"GGA_X_PW91": {"tweak": [a_pw6, b_pw6, c_pw6, d_pw6, f_pw6, alpha_pw6, expo_pw6],"alpha": 0.50}, "X_HF": {"alpha": 0.50}},
    "c_functionals": {"MGGA_C_BC95": {"tweak": [css, copp], "alpha": 0.731}, "C_MP2": {"alpha": 1.0, "ss": 0.0, "os": 0.269}},
    "citation": '    L. Goerigk, S.Grimme, J.Chem. Theory Compt. 7, 291-309, 2011 \n',
    "description": '    PWPB95 SOS Double Hybrid XC Functional\n',
}

ptpss = {
    "name": "PTPSS",
    "x_functionals": {"MGGA_X_TPSS": {"tweak": [0.15, 0.88491, 0.047, 0.872, 0.16952], "alpha": 0.50}, "X_HF": {"alpha": 0.50}},
    "c_functionals": {"MGGA_C_TPSS": {"tweak": [0.06080, 6.3, 0.53, 0.87, 0.50, 2.26], "alpha": 0.625}, "C_MP2": {"alpha": 1.0, "ss": 0.0, "os": 0.375}},
    "citation": '    L. Goerigk, S.Grimme, J. Chem. Theory Comput., 7, 291-309, 2011 \n',
    "description": '    PTPSS SOS Double Hybrid XC Functional\n',
}



functional_list = {
    "TEST-B2PLYP": b2plyp,
    "TEST-PBE0-2": pbe0_2,
    "TEST-PBE0-DH": pbe0_dh,
    "TEST-DSD-BLYP": dsd_blyp,
    "TEST-CORE-DSD-BLYP": core_dsd_blyp,
    "TEST-DSD-PBEP86": dsd_pbep86,
    "TEST-DSD-PBEPBE": dsd_pbepbe,
    "TEST-B2GPPLYP": b2gpplyp,
    "TEST-PWPB95": pwpb95,
    "TEST-PTPSS": ptpss,
}
