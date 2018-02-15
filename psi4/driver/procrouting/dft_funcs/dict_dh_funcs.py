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
    "x_functionals": {"GGA_X_B88": {"alpha": 0.47}},
    "c_functionals": {"GGA_C_LYP": {"alpha": 0.73}},
    "citation": '    S. Grimme, J. Chem. Phys., 124, 034108, 2006\n',
    "description": '    B2PLYP Double Hybrid Exchange-Correlation Functional\n',
}

dsd_blyp = {
    "name": "DSD-BLYP",
    "x_functionals": {"GGA_X_B88": {"alpha": 0.30}},
    "c_functionals": {"GGA_C_LYP": {"alpha": 0.56}, "C_MP2": {"alpha": 1.0, "os": 0.46, "ss": 0.40}},
    "citation": '    S. Kozuch, Phys. Chem. Chem. Phys., 13, 20104, 2011\n',
    "description": '    DSD-BLYP Dispersion-corrected SCS Double Hybrid XC Functional\n',
}

functional_list = {
    "TEST-B2PLYP": b2plyp,
    "TEST-DSD-BLYP": dsd_blyp,
}
