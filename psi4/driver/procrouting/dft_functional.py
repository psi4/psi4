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
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""
Module to provide lightweight definitions of functionals and
SuperFunctionals
"""
import re
import os
import math
from psi4 import core
from psi4.driver.qcdb import interface_dftd3 as dftd3
from psi4.driver.p4util.exceptions import *

## ==> Functionals <== ##


def build_s_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('S_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('S_X')
    # Tab in, trailing newlines
    fun.set_description('    Slater LSDA Exchange\n')
    # Tab in, trailing newlines
    fun.set_citation('    J.C. Slater, Phys. Rev., 81(3):385-390, 1951\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(False)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters

    # => End User-Customization <= #

    return fun


def build_b88_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('B88_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('B88_X')
    # Tab in, trailing newlines
    fun.set_description('    Becke88 GGA Exchange\n')
    # Tab in, trailing newlines
    fun.set_citation('    A.D. Becke, Phys. Rev. A, 38(6):3098-3100, 1988\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('B88_d', 0.0042)
    fun.set_parameter('B88_a', 1.0000)

    # => End User-Customization <= #

    return fun

def build_b86b_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('B86B_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('B86B_X')
    # Tab in, trailing newlines
    fun.set_description('    Becke86B GGA Exchange\n')
    # Tab in, trailing newlines
    fun.set_citation('    A. D. Becke, J. Chem. Phys. 85:7184, 1986.\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # => End User-Customization <= #

    return fun

def build_pw86_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('PW86_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('PW86_X')
    # Tab in, trailing newlines
    fun.set_description('    Perdew-Wang 1986 (PW86) GGA Exchange\n')
    # Tab in, trailing newlines
    fun.set_citation('    J. P. Perdew and W. Yue, Phys. Rev. B 33:8800(R), 1986.\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # => End User-Customization <= #

    return fun

def build_b3_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('B88_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('B3_X')
    # Tab in, trailing newlines
    fun.set_description('    Becke88 GGA Exchange (B3LYP weighting)\n')
    # Tab in, trailing newlines
    fun.set_citation('    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(0.8)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('B88_d', 0.0042)
    fun.set_parameter('B88_a', 0.9000)

    # => End User-Customization <= #

    return fun


def build_pbe_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('PBE_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('PBE_X')
    # Tab in, trailing newlines
    fun.set_description('    PBE GGA Exchange Hole (Parameter Free)\n')
    # Tab in, trailing newlines
    fun.set_citation('    J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('PBE_kp', 0.804)
    fun.set_parameter('PBE_mu', 0.2195149727645171)

    # => End User-Customization <= #

    return fun


def build_revpbe_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('PBE_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('revPBE_X')
    # Tab in, trailing newlines
    fun.set_description('    Revised PBE GGA Exchange Hole (Parameter Free)\n')
    # Tab in, trailing newlines
    fun.set_citation('    Zhang et. al., Phys. Rev. Lett., 80(4), 890, 1998\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('PBE_kp', 1.245)
    fun.set_parameter('PBE_mu', 0.2195149727645171)

    # => End User-Customization <= #

    return fun


def build_rpbe_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('RPBE_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('RPBE_X')
    # Tab in, trailing newlines
    fun.set_description('    RPBE GGA Exchange Hole (Parameter Free)\n')
    # Tab in, trailing newlines
    fun.set_citation('    Hammer et. al. Phys. Rev. B, 59(2), 7413-7421, 1999\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('PBE_kp', 0.804)
    fun.set_parameter('PBE_mu', 0.2195149727645171)

    # => End User-Customization <= #

    return fun


def build_sogga_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('SOGGA_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('SOGGA_X')
    # Tab in, trailing newlines
    fun.set_description('    Second Order GGA Exchange Hole (Parameter Free)\n')
    # Tab in, trailing newlines
    fun.set_citation('    Zhao et. al., J. Chem. Phys., 128(18), 184109, 2008\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('PBE_kp', 0.55208138)
    fun.set_parameter('PBE_mu', 10.0 / 81.0)

    # => End User-Customization <= #

    return fun


def build_pbesol_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('PBE_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('PBEsol_X')
    # Tab in, trailing newlines
    fun.set_description('    PBEsol GGA Exchange Hole (Parameter Free)\n')
    # Tab in, trailing newlines
    fun.set_citation('    J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('PBE_kp', 0.804)
    fun.set_parameter('PBE_mu', 10.0 / 81.0)

    # => End User-Customization <= #

    return fun


def build_pw91_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('PW91_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('PW91_X')
    # Tab in, trailing newlines
    fun.set_description('    PW91 Parameterized GGA Exchange\n')
    # Tab in, trailing newlines
    fun.set_citation('    J.P. Perdew et. al., Phys. Rev. B., 46(11), 6671-6687, 1992\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    k01 = math.pow(6.0 * math.pi * math.pi, 1.0 / 3.0)
    k02 = k01 * k01
    k04 = k02 * k02
    fun.set_parameter('PW91_a1', 0.19645 / (2.0 * k01))
    fun.set_parameter('PW91_a2', 7.79560 / (2.0 * k01))
    fun.set_parameter('PW91_a3', 0.27430 / (4.0 * k02))
    fun.set_parameter('PW91_a4', 0.15080 / (4.0 * k02))
    fun.set_parameter('PW91_a5', 100.000 / (4.0 * k02))
    fun.set_parameter('PW91_a6', 0.00400 / (16.0 * k04))

    # => End User-Customization <= #

    return fun


def build_b97_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('B97_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('B97_X')
    # Tab in, trailing newlines
    fun.set_description('    B97 Parameterized GGA Exchange\n')
    # Tab in, trailing newlines
    fun.set_citation('    A.D. Becke, J. Chem. Phys., 107(20), 8554-8560, 1997\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('B97_gamma', 0.004)

    # => End User-Customization <= #

    return fun


def build_vwn5_c_functional(name):

    # Call this first
    fun = core.Functional.build_base('VWN5_C')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('VWN5_C')
    # Tab in, trailing newlines
    fun.set_description('    VWN5 LSDA Correlation, QMC Parameters, VWN5 Spin Polarization\n')
    # Tab in, trailing newlines
    fun.set_citation('    S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200-1211, 1980\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(False)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('EcP_2', -0.10498)
    fun.set_parameter('EcP_3', 3.72744)
    fun.set_parameter('EcP_4', 12.9352)
    fun.set_parameter('EcF_2', -0.32500)
    fun.set_parameter('EcF_3', 7.06042)
    fun.set_parameter('EcF_4', 18.0578)
    fun.set_parameter('Ac_2', -0.00475840)
    fun.set_parameter('Ac_3', 1.13107)
    fun.set_parameter('Ac_4', 13.0045)

    # => End User-Customization <= #

    return fun


def build_vwn5rpa_c_functional(name):

    # Call this first
    fun = core.Functional.build_base('VWN5_C')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('VWN5RPA_C')
    # Tab in, trailing newlines
    fun.set_description('    VWN5 LSDA Correlation, RPA Parameters, VWN5 Spin Polarization\n')
    # Tab in, trailing newlines
    fun.set_citation('    S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200-1211, 1980\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(False)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('EcP_2', -0.409286)
    fun.set_parameter('EcP_3', 13.0720)
    fun.set_parameter('EcP_4', 42.7198)
    fun.set_parameter('EcF_2', -0.743294)
    fun.set_parameter('EcF_3', 20.1231)
    fun.set_parameter('EcF_4', 101.578)
    fun.set_parameter('Ac_2', -0.228344)
    fun.set_parameter('Ac_3', 1.06835)
    fun.set_parameter('Ac_4', 11.4813)

    # => End User-Customization <= #

    return fun


def build_vwn3_c_functional(name):

    # Call this first
    fun = core.Functional.build_base('VWN3_C')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('VWN3_C')
    # Tab in, trailing newlines
    fun.set_description('    VWN3 LSDA Correlation, QMC Parameters, VWN1 Spin Polarization\n')
    # Tab in, trailing newlines
    fun.set_citation('    S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200-1211, 1980\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(False)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('EcP_2', -0.10498)
    fun.set_parameter('EcP_3', 3.72744)
    fun.set_parameter('EcP_4', 12.9352)
    fun.set_parameter('EcF_2', -0.32500)
    fun.set_parameter('EcF_3', 7.06042)
    fun.set_parameter('EcF_4', 18.0578)

    # => End User-Customization <= #

    return fun


def build_vwn3rpa_c_functional(name):

    # Call this first
    fun = core.Functional.build_base('VWN3_C')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('VWN3RPA_C')
    # Tab in, trailing newlines
    fun.set_description('    VWN3 LSDA Correlation, RPA Parameters, VWN1 Spin Polarization\n')
    # Tab in, trailing newlines
    fun.set_citation('    S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200-1211, 1980\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(False)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('EcP_2', -0.409286)
    fun.set_parameter('EcP_3', 13.0720)
    fun.set_parameter('EcP_4', 42.7198)
    fun.set_parameter('EcF_2', -0.743294)
    fun.set_parameter('EcF_3', 20.1231)
    fun.set_parameter('EcF_4', 101.578)

    # => End User-Customization <= #

    return fun


def build_ws_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('wS_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('wS_X')
    # Tab in, trailing newlines
    fun.set_description('    Slater Short-Range LSDA Exchange\n')
    # Tab in, trailing newlines
    fun.set_citation('    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(False)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.3)

    # Custom parameters

    # => End User-Customization <= #

    return fun


def build_wpbe_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('wPBE_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('wPBE_X')
    # Tab in, trailing newlines
    fun.set_description('    PBE Short-Range GGA Exchange (HJS Formalism)\n')
    # Tab in, trailing newlines
    fun.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.3)

    # Custom parameters
    fun.set_parameter('A', 0.7572110)
    fun.set_parameter('B', -0.1063640)
    fun.set_parameter('C', -0.1186490)
    fun.set_parameter('D', 0.6096500)
    fun.set_parameter('E', -0.0477963)

    fun.set_parameter('Ha0', 0.0000000)
    fun.set_parameter('Ha1', 0.0000000)
    fun.set_parameter('Ha2', 0.0159941)
    fun.set_parameter('Ha3', 0.0852995)
    fun.set_parameter('Ha4', -0.1603680)
    fun.set_parameter('Ha5', 0.1526450)
    fun.set_parameter('Ha6', -0.0971263)
    fun.set_parameter('Ha7', 0.0422061)

    fun.set_parameter('Hb0', 1.0000000)
    fun.set_parameter('Hb1', 5.3331900)
    fun.set_parameter('Hb2', -12.478000)
    fun.set_parameter('Hb3', 11.098800)
    fun.set_parameter('Hb4', -5.1101300)
    fun.set_parameter('Hb5', 1.7146800)
    fun.set_parameter('Hb6', -0.6103800)
    fun.set_parameter('Hb7', 0.3075550)
    fun.set_parameter('Hb8', -0.0770547)
    fun.set_parameter('Hb9', 0.0334840)

    # => End User-Customization <= #

    return fun


def build_wpbesol_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('wPBE_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('wPBEsol_X')
    # Tab in, trailing newlines
    fun.set_description('    PBEsol Short-Range GGA Exchange (HJS Formalism)\n')
    # Tab in, trailing newlines
    fun.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.3)

    # Custom parameters
    fun.set_parameter('A', 0.7572110)
    fun.set_parameter('B', -0.1063640)
    fun.set_parameter('C', -0.1186490)
    fun.set_parameter('D', 0.6096500)
    fun.set_parameter('E', -0.0477963)

    fun.set_parameter('Ha0', 0.0000000)
    fun.set_parameter('Ha1', 0.0000000)
    fun.set_parameter('Ha2', 0.0047333)
    fun.set_parameter('Ha3', 0.0403304)
    fun.set_parameter('Ha4', -0.0574615)
    fun.set_parameter('Ha5', 0.0435395)
    fun.set_parameter('Ha6', -0.0216251)
    fun.set_parameter('Ha7', 0.0063721)

    fun.set_parameter('Hb0', 1.00000)
    fun.set_parameter('Hb1', 8.52056)
    fun.set_parameter('Hb2', -13.9885)
    fun.set_parameter('Hb3', 9.28583)
    fun.set_parameter('Hb4', -3.27287)
    fun.set_parameter('Hb5', 0.843499)
    fun.set_parameter('Hb6', -0.235543)
    fun.set_parameter('Hb7', 0.0847074)
    fun.set_parameter('Hb8', -0.0171561)
    fun.set_parameter('Hb9', 0.0050552)

    # => End User-Customization <= #

    return fun


def build_wb88_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('wB88_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('wB88_X')
    # Tab in, trailing newlines
    fun.set_description('    B88 Short-Range GGA Exchange (HJS Formalism)\n')
    # Tab in, trailing newlines
    fun.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n')

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.3)

    # Custom parameters
    fun.set_parameter('A', 0.7572110)
    fun.set_parameter('B', -0.1063640)
    fun.set_parameter('C', -0.1186490)
    fun.set_parameter('D', 0.6096500)
    fun.set_parameter('E', -0.0477963)

    fun.set_parameter('Ha0', 0.0000000)
    fun.set_parameter('Ha1', 0.0000000)
    fun.set_parameter('Ha2', 0.0253933)
    fun.set_parameter('Ha3', -0.0673075)
    fun.set_parameter('Ha4', 0.0891476)
    fun.set_parameter('Ha5', -0.0454168)
    fun.set_parameter('Ha6', -0.0076581)
    fun.set_parameter('Ha7', 0.0142506)

    fun.set_parameter('Hb0', 1.00000)
    fun.set_parameter('Hb1', -2.65060)
    fun.set_parameter('Hb2', 3.91108)
    fun.set_parameter('Hb3', -3.31509)
    fun.set_parameter('Hb4', 1.54485)
    fun.set_parameter('Hb5', -0.198386)
    fun.set_parameter('Hb6', -0.136112)
    fun.set_parameter('Hb7', 0.0647862)
    fun.set_parameter('Hb8', 0.0159586)
    fun.set_parameter('Hb9', -2.45066E-4)

    # => End User-Customization <= #

    return fun

def build_hf_x_functional(name):

    # Call this first
    fun = core.Functional.build_base('HF_X')

    # => End User-Customization <= #

    return fun

def build_primitive_functional(name):

    # Call this first
    key = name.upper()
    if (key[0] == 'W'):
        key = 'w' + key[1:]
    fun = core.Functional.build_base(key)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name(key)
    # Tab in, trailing newlines
    fun.set_description(fun.description())
    # Tab in, trailing newlines
    fun.set_citation(fun.citation())

    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(fun.is_gga())
    fun.set_meta(fun.is_meta())
    fun.set_alpha(fun.alpha())
    fun.set_omega(fun.omega())

    # Custom parameters
    # Always built-in for this functional

    # => End User-Customization <= #

    return fun

# Functional lookup table
functionals = {
        's_x'         : build_s_x_functional,
        'b88_x'       : build_b88_x_functional,
        'b86b_x'      : build_b86b_x_functional,
        'pw86_x'      : build_pw86_x_functional,
        'b3_x'        : build_b3_x_functional,
        'pbe_x'       : build_pbe_x_functional,
        'revpbe_x'    : build_revpbe_x_functional,
        'rpbe_x'      : build_rpbe_x_functional,
        'sogga_x'     : build_sogga_x_functional,
        'pbesol_x'    : build_pbesol_x_functional,
        'pw91_x'      : build_pw91_x_functional,
        'b97_x'       : build_b97_x_functional,
        'ws_x'        : build_ws_x_functional,
        'wb97_x'      : build_primitive_functional,
        'wpbe_x'      : build_wpbe_x_functional,
        'wpbesol_x'   : build_wpbesol_x_functional,
        'wb88_x'      : build_wb88_x_functional,
        'ft97b_x'     : build_primitive_functional,
        'm_x'         : build_primitive_functional,
        'lyp_c'       : build_primitive_functional,
        'pz81_c'      : build_primitive_functional,
        'p86_c'       : build_primitive_functional,
        'vwn5rpa_c'   : build_vwn5rpa_c_functional,
        'vwn5_c'      : build_vwn5_c_functional,
        'vwn3rpa_c'   : build_vwn3rpa_c_functional,
        'vwn3_c'      : build_vwn3_c_functional,
        'pw91_c'      : build_primitive_functional,
        'pw92_c'      : build_primitive_functional,
        'pbe_c'       : build_primitive_functional,
        'ft97_c'      : build_primitive_functional,
        'b_c'         : build_primitive_functional,
        'm_c'         : build_primitive_functional,
        'pbea_c'      : build_primitive_functional,
        'pw92a_c'     : build_primitive_functional,
        'wpbe_c'      : build_primitive_functional,
        'wpw92_c'     : build_primitive_functional,
        'hf_x'        : build_hf_x_functional,
    }


def build_functional(alias):
    name = alias.lower()
    return functionals[name](name)


def functional_list():
    val = []
    for key in functionals.keys():
        val.append(functionals[key](key))
    return val

## ==> SuperFunctionals <== ##


def build_ws_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wS_X')
    # Tab in, trailing newlines
    sup.set_description('    Slater Short-Range LSDA Exchange\n')
    # Tab in, trailing newlines
    sup.set_citation('    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n')

    # Add member functionals
    sup.add_x_functional(build_functional('wS_X'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_wpbe_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wPBE_X')
    # Tab in, trailing newlines
    sup.set_description('    PBE Short-Range GGA Exchange (HJS Model)\n')
    # Tab in, trailing newlines
    sup.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n')

    # Add member functionals
    sup.add_x_functional(build_functional('wPBE_X'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_wpbesol_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wPBEsol_X')
    # Tab in, trailing newlines
    sup.set_description('    PBEsol Short-Range GGA Exchange (HJS Model)\n')
    # Tab in, trailing newlines
    sup.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n')

    # Add member functionals
    sup.add_x_functional(build_functional('wPBEsol_X'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)

def build_wpw92_c_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wPW92_C')
    # Tab in, trailing newlines
    sup.set_description('    Short-Range PW92 Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    TODO\n')

    # Add member functionals
    sup.add_c_functional(build_functional('wPW92_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.3)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)

def build_wpbe_c_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wPBE_C')
    # Tab in, trailing newlines
    sup.set_description('    Short-Range PBE Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    TODO\n')

    # Add member functionals
    sup.add_c_functional(build_functional('wPBE_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.5)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)

def build_wpbe2_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wPBE2')
    # Tab in, trailing newlines
    sup.set_description('    Double-Hybrid PBE LRC Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    TODO\n')

    # Add member functionals
    sup.add_x_functional(build_functional('wPBE_X'))
    sup.add_c_functional(build_functional('wPBE_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.5)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)

def build_wb88_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wB88_X')
    # Tab in, trailing newlines
    sup.set_description('    B88 Short-Range GGA Exchange (HJS Model)\n')
    # Tab in, trailing newlines
    sup.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n')

    # Add member functionals
    sup.add_x_functional(build_functional('wB88_X'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_svwn_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('SVWN')
    # Tab in, trailing newlines
    sup.set_description('    SVWN3 (RPA) LSDA Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n')

    # Add member functionals
    sup.add_x_functional(build_functional('S_X'))
    sup.add_c_functional(build_functional('VWN3RPA_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_blyp_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('BLYP')
    # Tab in, trailing newlines
    sup.set_description('    BLYP GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n    B. Miehlich et. al., Chem. Phys. Lett., 157(3), 200-206 1989\n')

    # Add member functionals
    sup.add_x_functional(build_functional('B88_X'))
    sup.add_c_functional(build_functional('LYP_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)

def build_b86bpbe_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('B86BPBE')
    # Tab in, trailing newlines
    sup.set_description('    B86BPBE GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    A. D. Becke, J. Chem. Phys. 85:7184, 1986.\n')

    # Add member functionals
    sup.add_x_functional(build_functional('B86B_X'))
    sup.add_c_functional(build_functional('PBE_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)

def build_pw86pbe_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('PW86PBE')
    # Tab in, trailing newlines
    sup.set_description('    PW86PBE GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    J. P. Perdew and W. Yue, Phys. Rev. B 33:8800(R), 1986.\n')

    # Add member functionals
    sup.add_x_functional(build_functional('PW86_X'))
    sup.add_c_functional(build_functional('PBE_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)

def build_pw91_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('PW91')
    # Tab in, trailing newlines
    sup.set_description('    PW91 GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    J.P. Perdew et. al., Phys. Rev. B., 46(11), 6671-6687, 1992\n')

    # Add member functionals
    sup.add_x_functional(build_functional('PW91_X'))
    sup.add_c_functional(build_functional('PW91_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_bp86_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('BP86')
    # Tab in, trailing newlines
    sup.set_description('    BP86 GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    Null\n')

    # Add member functionals
    sup.add_x_functional(build_functional('B88_X'))
    sup.add_c_functional(build_functional('P86_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_ft97_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('FT97')
    # Tab in, trailing newlines
    sup.set_description('   FT97 GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    M. Filatov and W. Theil, Int. J. Quant. Chem., 62, 603-616, 1997\n')

    # Add member functionals
    sup.add_x_functional(build_functional('FT97B_X'))
    sup.add_c_functional(build_functional('FT97_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_pbe_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('PBE')
    # Tab in, trailing newlines
    sup.set_description('    PBE GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996\n')

    # Add member functionals
    sup.add_x_functional(build_functional('PBE_X'))
    sup.add_c_functional(build_functional('PBE_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_pbe0_superfunctional(name, npoints, deriv):

    sup = build_pbe_superfunctional(name, npoints, deriv)[0]
    sup.set_name('PBE0')
    sup.set_description('    PBE0 Hybrid GGA Exchange-Correlation Functional\n')
    sup.set_citation('    Adamo et. al., J. Chem. Phys., 110(13), 6158, 1999\n')
    sup.set_x_alpha(0.25)
    return (sup, False)


def build_pbeh3c_superfunctional(name, npoints, deriv):

    # call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => user-customization <= #

    # no spaces, keep it short and according to convention
    sup.set_name('PBEH3C')
    # tab in, trailing newlines
    sup.set_description('    PBE-3C Hybrid GGA Exchange-Correlation Functional\n')
    # tab in, trailing newlines
    sup.set_citation('    Grimme et. al., J. Chem. Phys., 143, 054107, 2015\n')

    # add member functionals
    pbe_x3c = build_functional('PBE_X')
    pbe_x3c.set_parameter('PBE_kp', 1.0245)
    pbe_x3c.set_parameter('PBE_mu', 0.12345679)
    sup.add_x_functional(pbe_x3c)
    
    pbe_c3c = build_functional('PBE_C')
    pbe_c3c.set_parameter('bet', 0.03)
    sup.add_c_functional(pbe_c3c)

    # set gks up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.42)
    sup.set_c_alpha(0.0)

    # => end user-customization <= #

    # call this last
    sup.allocate()
    return sup, False


def build_sogga_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('SOGGA')
    # Tab in, trailing newlines
    sup.set_description('    Second Order GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    Zhao et. al., J. Chem. Phys., 128(18), 184109, 2008\n')

    # Add member functionals
    sup.add_x_functional(build_functional('SOGGA_X'))

    C = build_functional('PBE_C')
    C.set_parameter('bet', 0.037526)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_b3lyp_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('B3LYP')
    # Tab in, trailing newlines
    sup.set_description('    B3LYP Hybrid-GGA Exchange-Correlation Functional (VWN1-RPA)\n')
    # Tab in, trailing newlines
    sup.set_citation('    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n')

    # Add member functionals
    b3 = build_functional('B3_X')
    b3.set_alpha(1.0)
    sup.add_x_functional(b3)
    lyp = build_functional('LYP_C')
    lyp.set_alpha(0.81)
    vwn = build_functional('VWN3RPA_C')
    vwn.set_alpha(0.19)
    sup.add_c_functional(vwn)
    sup.add_c_functional(lyp)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.2)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_hf_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('HF_X')
    # Tab in, trailing newlines
    sup.set_description('    Hartree-Fock Exchange Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    \n')

    # Add member functionals
    hf_x = build_functional('hf_x')
    hf_x.set_alpha(1.0)
    sup.add_x_functional(hf_x)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(1.0)
    sup.set_c_alpha(0.0)

    sup.allocate()
    return (sup, False)


def build_b3lyp5_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('B3LYP5')
    # Tab in, trailing newlines
    sup.set_description('    B3LYP5 Hybrid-GGA Exchange-Correlation Functional (VWN5)\n')
    # Tab in, trailing newlines
    sup.set_citation('    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n')

    # Add member functionals
    b3 = build_functional('B3_X')
    b3.set_alpha(1.0)
    sup.add_x_functional(b3)
    lyp = build_functional('LYP_C')
    lyp.set_alpha(0.81)
    vwn = build_functional('VWN5_C')
    vwn.set_alpha(0.19)
    sup.add_c_functional(lyp)
    sup.add_c_functional(vwn)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.2)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_b970_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('B97-0')
    # Tab in, trailing newlines
    sup.set_description('    B97-0 Hybrid-GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    A.D. Becke, J. Chem. Phys., 107(20), 8554-8560, 1997\n')

    # Add member functionals
    X = build_functional('B97_X')
    X.set_name('B97-0_X')
    X.set_alpha(1.0 / 0.8057)

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 0.8094)
    X.set_parameter('B97_a1', 0.5073)
    X.set_parameter('B97_a2', 0.7481)

    C = build_functional('B_C')
    C.set_name('B97-0_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 0.9454)
    C.set_parameter('B97_os_a1', 0.7471)
    C.set_parameter('B97_os_a2', -4.5961)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 0.1737)
    C.set_parameter('B97_ss_a1', 2.3487)
    C.set_parameter('B97_ss_a2', -2.4868)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.1943)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_b971_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('B97-1')
    # Tab in, trailing newlines
    sup.set_description('    B97-1 Hybrid-GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    F.A. Hamprecht et. al., J. Chem. Phys., 109(15), 6264-6271, 1998\n')

    # Add member functionals
    X = build_functional('B97_X')
    X.set_name('B97-1_X')
    X.set_alpha(1.0 / 0.79)

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 0.789518)
    X.set_parameter('B97_a1', 0.573805)
    X.set_parameter('B97_a2', 0.660975)

    C = build_functional('B_C')
    C.set_name('B97-1_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 0.955689)
    C.set_parameter('B97_os_a1', 0.788552)
    C.set_parameter('B97_os_a2', -5.47869)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 0.0820011)
    C.set_parameter('B97_ss_a1', 2.71681)
    C.set_parameter('B97_ss_a2', -2.87103)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.21)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_b972_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('B97-2')
    # Tab in, trailing newlines
    sup.set_description('    B97-2 Hybrid-GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    P.J. Wilson et. al., J. Chem. Phys., 115(20), 9233-9242, 2001\n')

    # Add member functionals
    X = build_functional('B97_X')
    X.set_name('B97-2_X')
    X.set_alpha(1.0 / 0.79)

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 0.827642)
    X.set_parameter('B97_a1', 0.047840)
    X.set_parameter('B97_a2', 1.761250)

    C = build_functional('B_C')
    C.set_name('B97-2_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 0.999849)
    C.set_parameter('B97_os_a1', 1.40626)
    C.set_parameter('B97_os_a2', -7.44060)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 0.585808)
    C.set_parameter('B97_ss_a1', -0.691682)
    C.set_parameter('B97_ss_a2', 0.394796)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.21)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_b97d_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('B97-D2P4')
    # Tab in, trailing newlines
    sup.set_description('    B97-D Pure-GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    S. Grimme, J. Comput. Chem., 27, 1787-1799, 2006\n')

    # Add member functionals
    X = build_functional('B97_X')
    X.set_name('B97-D_X')
    X.set_alpha(1.0)

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 1.08662)
    X.set_parameter('B97_a1', -0.52127)
    X.set_parameter('B97_a2', 3.25429)

    C = build_functional('B_C')
    C.set_name('B97-D_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 0.69041)
    C.set_parameter('B97_os_a1', 6.30270)
    C.set_parameter('B97_os_a2', -14.9712)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 0.22340)
    C.set_parameter('B97_ss_a1', -1.56208)
    C.set_parameter('B97_ss_a2', 1.94293)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, ('b97-d', 'd2p4'))


def build_hcth_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('HCTH')
    # Tab in, trailing newlines
    sup.set_description('    HCTH Pure-GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    F.A. Hamprecht et. al., J. Chem. Phys., 109(15), 6264-6271\n')

    # Add member functionals
    X = build_functional('B97_X')
    X.set_name('HCTH_X')
    X.set_alpha(1.0)

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 1.09320)
    X.set_parameter('B97_a1', -0.744056)
    X.set_parameter('B97_a2', 5.59920)
    X.set_parameter('B97_a3', -6.78549)
    X.set_parameter('B97_a4', 4.49357)

    C = build_functional('B_C')
    C.set_name('HCTH_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 0.729974)
    C.set_parameter('B97_os_a1', 3.35287)
    C.set_parameter('B97_os_a2', -11.5430)
    C.set_parameter('B97_os_a3', 8.08564)
    C.set_parameter('B97_os_a4', -4.47857)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 0.222601)
    C.set_parameter('B97_ss_a1', -0.0338622)
    C.set_parameter('B97_ss_a2', -0.0125170)
    C.set_parameter('B97_ss_a3', -0.802496)
    C.set_parameter('B97_ss_a4', 1.55396)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_hcth120_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('HCTH120')
    # Tab in, trailing newlines
    sup.set_description('    HCTH120 Pure-GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    A.D. Boese, et. al., J. Chem. Phys., 112(4), 1670-1678, 2000\n')

    # Add member functionals
    X = build_functional('B97_X')
    X.set_name('HCTH120_X')
    X.set_alpha(1.0)

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 1.09163)
    X.set_parameter('B97_a1', -0.747215)
    X.set_parameter('B97_a2', 5.07833)
    X.set_parameter('B97_a3', -4.10746)
    X.set_parameter('B97_a4', 1.17173)

    C = build_functional('B_C')
    C.set_name('HCTH120_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 0.514730)
    C.set_parameter('B97_os_a1', 6.92982)
    C.set_parameter('B97_os_a2', -24.7073)
    C.set_parameter('B97_os_a3', 23.1098)
    C.set_parameter('B97_os_a4', -11.3234)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 0.489508)
    C.set_parameter('B97_ss_a1', -0.260699)
    C.set_parameter('B97_ss_a2', 0.432917)
    C.set_parameter('B97_ss_a3', -1.99247)
    C.set_parameter('B97_ss_a4', 2.48531)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_hcth147_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('HCTH147')
    # Tab in, trailing newlines
    sup.set_description('    HCTH147 Pure-GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    A.D. Boese, et. al., J. Chem. Phys., 112(4), 1670-1678, 2000\n')

    # Add member functionals
    X = build_functional('B97_X')
    X.set_name('HCTH147_X')
    X.set_alpha(1.0)

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 1.09025)
    X.set_parameter('B97_a1', -0.799194)
    X.set_parameter('B97_a2', 5.57212)
    X.set_parameter('B97_a3', -5.86760)
    X.set_parameter('B97_a4', 3.04544)

    C = build_functional('B_C')
    C.set_name('HCTH147_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 0.542352)
    C.set_parameter('B97_os_a1', 7.01464)
    C.set_parameter('B97_os_a2', -28.3822)
    C.set_parameter('B97_os_a3', 35.0329)
    C.set_parameter('B97_os_a4', -20.4284)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 0.562576)
    C.set_parameter('B97_ss_a1', 0.0171436)
    C.set_parameter('B97_ss_a2', -1.30636)
    C.set_parameter('B97_ss_a3', 1.05747)
    C.set_parameter('B97_ss_a4', 0.885429)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_hcth407_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('HCTH407')
    # Tab in, trailing newlines
    sup.set_description('    HCTH407 Pure-GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    A.D. Boese and N.C. Handy, J. Chem. Phys., 114(13), 5497-5503, 2001\n')

    # Add member functionals
    X = build_functional('B97_X')
    X.set_name('HCTH407_X')
    X.set_alpha(1.0)

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 1.08184)
    X.set_parameter('B97_a1', -0.518339)
    X.set_parameter('B97_a2', 3.42562)
    X.set_parameter('B97_a3', -2.62901)
    X.set_parameter('B97_a4', 2.28855)

    C = build_functional('B_C')
    C.set_name('HCTH407_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 0.589076)
    C.set_parameter('B97_os_a1', 4.42374)
    C.set_parameter('B97_os_a2', -19.2218)
    C.set_parameter('B97_os_a3', 42.5721)
    C.set_parameter('B97_os_a4', -42.0052)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 1.18777)
    C.set_parameter('B97_ss_a1', -2.40292)
    C.set_parameter('B97_ss_a2', 5.61741)
    C.set_parameter('B97_ss_a3', -9.17923)
    C.set_parameter('B97_ss_a4', 6.24798)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_wsvwn_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wSVWN')
    # Tab in, trailing newlines
    sup.set_description('    LSDA SR-XC Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n')

    # Add member functionals
    sup.add_x_functional(build_functional('wS_X'))
    sup.add_c_functional(build_functional('VWN3RPA_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_wpbe_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wPBE')
    # Tab in, trailing newlines
    sup.set_description('    PBE SR-XC Functional (HJS Model)\n')
    # Tab in, trailing newlines
    sup.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n')

    # Add member functionals
    sup.add_x_functional(build_functional('wPBE_X'))
    sup.add_c_functional(build_functional('PBE_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.4)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_wpbe0_superfunctional(name, npoints, deriv):

    sup = build_wpbe_superfunctional(name, npoints, deriv)[0]
    sup.set_name('wPBE0')
    sup.set_description('    PBE0 SR-XC Functional (HJS Model)\n')
    sup.set_x_omega(0.3)
    sup.set_x_alpha(0.25)
    return (sup, False)


def build_wpbesol_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wPBEsol')
    # Tab in, trailing newlines
    sup.set_description('    PBEsol SR-XC Functional (HJS Model)\n')
    # Tab in, trailing newlines
    sup.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n')

    # Add member functionals
    sup.add_x_functional(build_functional('wPBEsol_X'))
    sup.add_c_functional(build_functional('PBE_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.4)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_wpbesol0_superfunctional(name, npoints, deriv):

    sup = build_wpbesol_superfunctional(name, npoints, deriv)[0]
    sup.set_name('wPBEsol0')
    sup.set_description('    PBEsol0 SR-XC Functional (HJS Model)\n')
    sup.set_x_omega(0.3)
    sup.set_x_alpha(0.25)
    return (sup, False)


def build_wblyp_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wBLYP')
    # Tab in, trailing newlines
    sup.set_description('    BLYP SR-XC Functional (HJS Model)\n')
    # Tab in, trailing newlines
    sup.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n')

    # Add member functionals
    sup.add_x_functional(build_functional('wB88_X'))
    sup.add_c_functional(build_functional('LYP_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_wb97_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wB97')
    # Tab in, trailing newlines
    sup.set_description('    Parameterized LRC B97 GGA XC Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    J.-D. Chai and M. Head-Gordon, J. Chem. Phys., 128, 084106, 2008\n')

    # Add member functionals
    X = build_functional('wB97_X')
    X.set_name('wB97_X')
    X.set_alpha(1.0)

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 1.0)
    X.set_parameter('B97_a1', 1.13116E0)
    X.set_parameter('B97_a2', -2.74915E0)
    X.set_parameter('B97_a3', 1.20900E1)
    X.set_parameter('B97_a4', -5.71642E0)

    C = build_functional('B_C')
    C.set_name('wB97_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 1.0)
    C.set_parameter('B97_os_a1', 3.99051E0)
    C.set_parameter('B97_os_a2', -1.70066E1)
    C.set_parameter('B97_os_a3', 1.07292E0)
    C.set_parameter('B97_os_a4', 8.88211E0)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 1.0)
    C.set_parameter('B97_ss_a1', -2.55352E0)
    C.set_parameter('B97_ss_a2', 1.18926E1)
    C.set_parameter('B97_ss_a3', -2.69452E1)
    C.set_parameter('B97_ss_a4', 1.70927E1)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.4)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_wb97x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wB97X')
    # Tab in, trailing newlines
    sup.set_description('    Parameterized Hybrid LRC B97 GGA XC Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    J.-D. Chai and M. Head-Gordon, J. Chem. Phys., 128, 084106, 2008\n')

    # Add member functionals
    X = build_functional('wB97_X')
    X.set_name('wB97X_X')
    X.set_alpha(1.0 / (1.0 - 0.157706))

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 8.42294E-1)
    X.set_parameter('B97_a1', 7.26479E-1)
    X.set_parameter('B97_a2', 1.04760E0)
    X.set_parameter('B97_a3', -5.70635E0)
    X.set_parameter('B97_a4', 1.32794E1)

    C = build_functional('B_C')
    C.set_name('wB97X_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 1.0)
    C.set_parameter('B97_os_a1', 2.37031E0)
    C.set_parameter('B97_os_a2', -1.13995E1)
    C.set_parameter('B97_os_a3', 6.58405E0)
    C.set_parameter('B97_os_a4', -3.78132E0)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 1.0)
    C.set_parameter('B97_ss_a1', -4.33879E0)
    C.set_parameter('B97_ss_a2', 1.82308E1)
    C.set_parameter('B97_ss_a3', -3.17430E1)
    C.set_parameter('B97_ss_a4', 1.72901E1)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.4)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.157706)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_wb97xd_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wB97X-D')
    # Tab in, trailing newlines
    sup.set_description('    Parameterized Hybrid LRC B97 GGA XC Functional with Dispersion\n')
    # Tab in, trailing newlines
    sup.set_citation('    J.-D. Chai and M. Head-Gordon, Phys. Chem. Chem. Phys., 10, 6615-6620, 2008\n')

    # Add member functionals
    alpha = 2.22036E-1
    omega = 0.2
    X = build_functional('wB97_X')
    X.set_name('wB97X_X')
    X.set_alpha(1.0 / (1.0 - alpha))

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 7.77964E-1)     # Table 1: c_{x\sigma,0}
    X.set_parameter('B97_a1', 6.61160E-1)     # Table 1: c_{x\sigma,1}
    X.set_parameter('B97_a2', 5.74541E-1)     # Table 1: c_{x\sigma,2}
    X.set_parameter('B97_a3', -5.25671E0)     # Table 1: c_{x\sigma,3}
    X.set_parameter('B97_a4', 1.16386E1)      # Table 1: c_{x\sigma,4}

    C = build_functional('B_C')
    C.set_name('wB97X_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 1.0)         # Table 1: c_{c\alpha\beta,0}
    C.set_parameter('B97_os_a1', 1.79413E0)   # Table 1: c_{c\alpha\beta,1}
    C.set_parameter('B97_os_a2', -1.20477E1)  # Table 1: c_{c\alpha\beta,2}
    C.set_parameter('B97_os_a3', 1.40847E1)   # Table 1: c_{c\alpha\beta,3}
    C.set_parameter('B97_os_a4', -8.50809E0)  # Table 1: c_{c\alpha\beta,4}

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 1.0)         # Table 1: c_{c\sigma\sigma,0}
    C.set_parameter('B97_ss_a1', -6.90539E0)  # Table 1: c_{c\sigma\sigma,1}
    C.set_parameter('B97_ss_a2', 3.13343E1)   # Table 1: c_{c\sigma\sigma,2}
    C.set_parameter('B97_ss_a3', -5.10533E1)  # Table 1: c_{c\sigma\sigma,3}
    C.set_parameter('B97_ss_a4', 2.64423E1)   # Table 1: c_{c\sigma\sigma,4}

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(omega)   # Table 1: omega
    sup.set_c_omega(0.0)
    sup.set_x_alpha(alpha)   # Table 1: c_x
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, ('wB97', '-CHG'))


def build_m05_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('M05')
    # Tab in, trailing newlines
    sup.set_description('    Heavily Parameterized Hybrid Meta-GGA XC Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    Zhao et. al., J. Chem. Phys., 123, 161103, 2005\n')

    # Add member functionals
    X = build_functional('M_X')
    X.set_name('M05_X')
    X.set_alpha(1.0)

    # LSDA Exchange type is Slater, no parameters

    # GGA Exchange type is PBE, special parameters because Truhlar is lazy
    C1 = 3.36116E-3
    C2 = 4.49267E-3
    K0 = 3.0 / 2.0 * math.pow(3.0 / (math.pi * 4.0), 1.0 / 3.0)
    k0 = math.pow(6.0 * math.pi * math.pi, 1.0 / 3.0)
    kp = C1 / (C2 * K0)
    mu = 4.0 * k0 * k0 * kp * C2
    X.set_parameter('PBE_kp', kp)  # Different effective kp
    X.set_parameter('PBE_mu', mu)  # Different effective mu

    # Meta Exchange type is insane mess of w power series expansion
    X.set_parameter('Meta_a0', 1.0)
    X.set_parameter('Meta_a1', 0.08151)
    X.set_parameter('Meta_a2', -0.43956)
    X.set_parameter('Meta_a3', -3.22422)
    X.set_parameter('Meta_a4', 2.01819)
    X.set_parameter('Meta_a5', 8.79431)
    X.set_parameter('Meta_a6', -0.00295)
    X.set_parameter('Meta_a7', 9.82029)
    X.set_parameter('Meta_a8', -4.82351)
    X.set_parameter('Meta_a9', -48.17574)
    X.set_parameter('Meta_a10', 3.64802)
    X.set_parameter('Meta_a11', 34.02248)

    C = build_functional('M_C')
    C.set_name('M05_C')

    # LSDA Correlation type is PW92, no parameters

    # GGA Correlation type is B97
    C.set_parameter('B97_os_gamma', 0.0031 * 2.0)
    C.set_parameter('B97_os_a0', 1.0)
    C.set_parameter('B97_os_a1', 3.78569)
    C.set_parameter('B97_os_a2', -14.15261)
    C.set_parameter('B97_os_a3', -7.46589)
    C.set_parameter('B97_os_a4', 17.94491)

    C.set_parameter('B97_ss_gamma', 0.06)
    C.set_parameter('B97_ss_a0', 1.0)
    C.set_parameter('B97_ss_a1', 3.77344)
    C.set_parameter('B97_ss_a2', -26.04463)
    C.set_parameter('B97_ss_a3', 30.69913)
    C.set_parameter('B97_ss_a4', -9.22695)

    # Meta Correlation type is Becke metric, no parameters

    # Add the functionals in
    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.28)  # Hartree-Fock exact exchange
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_m05_2x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('M05-2X')
    # Tab in, trailing newlines
    sup.set_description('    Heavily Parameterized Hybrid Meta-GGA XC Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    Zhao et. al., J. Chem. Theory Comput., 2, 364, 2006\n')

    # Add member functionals
    X = build_functional('M_X')
    X.set_name('M05_2X_X')
    X.set_alpha(1.0)

    # LSDA Exchange type is Slater, no parameters

    # GGA Exchange type is PBE, special parameters because Truhlar is lazy
    C1 = 3.36116E-3
    C2 = 4.49267E-3
    K0 = 3.0 / 2.0 * math.pow(3.0 / (math.pi * 4.0), 1.0 / 3.0)
    k0 = math.pow(6.0 * math.pi * math.pi, 1.0 / 3.0)
    kp = C1 / (C2 * K0)
    mu = 4.0 * k0 * k0 * kp * C2
    X.set_parameter('PBE_kp', kp)
    X.set_parameter('PBE_mu', mu)

# Meta Exchange type is insane mess of w power series expansion
    X.set_parameter('Meta_a0', 1.0)
    X.set_parameter('Meta_a1', -0.56833)
    X.set_parameter('Meta_a2', -1.30057)
    X.set_parameter('Meta_a3', 5.50070)
    X.set_parameter('Meta_a4', 9.06402)
    X.set_parameter('Meta_a5', -32.21075)
    X.set_parameter('Meta_a6', -23.73298)
    X.set_parameter('Meta_a7', 70.22996)
    X.set_parameter('Meta_a8', 29.88614)
    X.set_parameter('Meta_a9', -60.25778)
    X.set_parameter('Meta_a10', -13.22205)
    X.set_parameter('Meta_a11', 15.23694)

    C = build_functional('M_C')
    C.set_name('M05_2X_C')

    # LSDA Correlation type is PW92, no parameters

    # GGA Correlation type is B97
    C.set_parameter('B97_os_gamma', 0.0031 * 2.0)
    C.set_parameter('B97_os_a0', 1.00000)
    C.set_parameter('B97_os_a1', 1.09297)
    C.set_parameter('B97_os_a2', -3.79171)
    C.set_parameter('B97_os_a3', 2.82810)
    C.set_parameter('B97_os_a4', -10.58909)

    C.set_parameter('B97_ss_gamma', 0.06)
    C.set_parameter('B97_ss_a0', 1.00000)
    C.set_parameter('B97_ss_a1', -3.05430)
    C.set_parameter('B97_ss_a2', 7.61854)
    C.set_parameter('B97_ss_a3', 1.47665)
    C.set_parameter('B97_ss_a4', -11.92365)

    # Meta Correlation type is Becke metric, no parameters

    # Add the functionals in
    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.56)  # Hartree-Fock exact exchange
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_dldf_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('dlDF')
    # Tab in, trailing newlines
    sup.set_description('    Dispersionless Hybrid Meta-GGA XC Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    Pernal et. al., Phys. Rev. Lett., 103, 263201, 2009\n')

    # Add member functionals
    X = build_functional('M_X')
    X.set_name('dlDF_X')
    X.set_alpha(1.0)

    # LSDA Exchange type is Slater, no parameters

    # GGA Exchange type is PBE
    kp = 4.8827323
    mu = 0.3511128
    X.set_parameter('PBE_kp', kp)
    X.set_parameter('PBE_mu', mu)

# Meta Exchange is a reparametrized truncation of Truhlar's functional
    X.set_parameter('Meta_a0', 1.0)
    X.set_parameter('Meta_a1', -0.1637571)
    X.set_parameter('Meta_a2', -0.1880028)
    X.set_parameter('Meta_a3', -0.4490609)
    X.set_parameter('Meta_a4', -0.0082359)

    C = build_functional('M_C')
    C.set_name('dlDF_C')

    # LSDA Correlation type is PW92, no parameters

    # GGA Correlation type is B97
    C.set_parameter('B97_os_gamma', 0.0031 * 2.0)
    C.set_parameter('B97_os_a0', 1.00000)
    C.set_parameter('B97_os_a1', 5.9515308)
    C.set_parameter('B97_os_a2', -11.1602877)

    C.set_parameter('B97_ss_gamma', 0.06)
    C.set_parameter('B97_ss_a0', 1.00000)
    C.set_parameter('B97_ss_a1', -2.5960897)
    C.set_parameter('B97_ss_a2', 2.2233793)

    # Meta Correlation type is Becke metric, no parameters

    # Add the functionals in
    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.6144129)  # Hartree-Fock exact exchange
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_dldfd09_superfunctional(name, npoints, deriv):
    sup, disp = build_dldf_superfunctional(name, npoints, deriv)
    sup.set_name('dlDF+D09')

    return (sup, ('dlDF', '-DAS2009'))

def build_dldfd10_superfunctional(name, npoints, deriv):
    sup, disp = build_dldf_superfunctional(name, npoints, deriv)
    sup.set_name('dlDF+D')

    return (sup, ('dlDF', '-DAS2010'))


def build_hfd_superfunctional(name, npoints, deriv):

    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)
    sup.set_name('HF+D')
    sup.set_x_alpha(1.0)

    sup.allocate()
    return (sup, ('HF', '-DAS2010'))


def build_b2plyp_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('B2PLYP')
    # Tab in, trailing newlines
    sup.set_description('    B2PLYP Double Hybrid Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    S. Grimme, J. Chem. Phys., 124, 034108, 2006\n')

    # Add member functionals
    becke = build_functional('B88_X')
    becke.set_alpha(1.0)
    sup.add_x_functional(becke)
    lyp = build_functional('LYP_C')
    lyp.set_alpha(1.0)
    sup.add_c_functional(lyp)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.53)
    sup.set_c_alpha(0.27)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_wb97x_2tqz_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wB97X-2(TQZ)')
    # Tab in, trailing newlines
    sup.set_description('    Double Hybrid LRC B97 GGA XC Functional (TQZ parametrization)\n')
    # Tab in, trailing newlines
    sup.set_citation('    J.-D. Chai and M. Head-Gordon, J. Chem. Phys., 131, 174105, 2009\n')

    # Add member functionals
    X = build_functional('wB97_X')
    X.set_name('wB97X_X')
    X.set_alpha(1.0 / (1.0 - 0.636158))

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 3.15503E-1)
    X.set_parameter('B97_a1', 1.04772E0)
    X.set_parameter('B97_a2', -2.33506E0)
    X.set_parameter('B97_a3', 3.19909E0)

    C = build_functional('B_C')
    C.set_name('wB97X_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 5.18198E-1)
    C.set_parameter('B97_os_a1', -5.85956E-1)
    C.set_parameter('B97_os_a2', 4.27080E0)
    C.set_parameter('B97_os_a3', -6.48897E0)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 9.08460E-1)
    C.set_parameter('B97_ss_a1', -2.80936E0)
    C.set_parameter('B97_ss_a2', 6.02676E0)
    C.set_parameter('B97_ss_a3', -4.56981E0)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.636158)
    sup.set_c_alpha(1.0)
    sup.set_c_os_alpha(0.447105)
    sup.set_c_ss_alpha(0.529319)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_wb97x_2lp_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wB97X-2(LP)')
    # Tab in, trailing newlines
    sup.set_description('    Double Hybrid LRC B97 GGA XC Functional (Large Pople parametrization)\n')
    # Tab in, trailing newlines
    sup.set_citation('    J.-D. Chai and M. Head-Gordon, J. Chem. Phys., 131, 174105, 2009\n')

    # Add member functionals
    X = build_functional('wB97_X')
    X.set_name('wB97X_X')
    X.set_alpha(1.0 / (1.0 - 0.678792))

    X.set_parameter('B97_gamma', 0.004)
    X.set_parameter('B97_a0', 2.51767E-1)
    X.set_parameter('B97_a1', 1.57375E0)
    X.set_parameter('B97_a2', -5.26624E0)
    X.set_parameter('B97_a3', 6.74313E0)

    C = build_functional('B_C')
    C.set_name('wB97X_C')

    C.set_parameter('B97_os_gamma', 0.006)
    C.set_parameter('B97_os_a0', 5.53261E-1)
    C.set_parameter('B97_os_a1', -1.16626E0)
    C.set_parameter('B97_os_a2', 6.84409E0)
    C.set_parameter('B97_os_a3', -8.90640E0)

    C.set_parameter('B97_ss_gamma', 0.2)
    C.set_parameter('B97_ss_a0', 1.15698E0)
    C.set_parameter('B97_ss_a1', -3.31669E0)
    C.set_parameter('B97_ss_a2', 6.27265E0)
    C.set_parameter('B97_ss_a3', -4.51464E0)

    sup.add_x_functional(X)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.678792)
    sup.set_c_alpha(1.0)
    sup.set_c_os_alpha(0.477992)
    sup.set_c_ss_alpha(0.581569)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_pbe0_2_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('PBE0-2')
    # Tab in, trailing newlines
    sup.set_description('    PBE0-2 Double Hydrid Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    J. Chai, Chem. Phys. Lett., 538, 121-125, 2012\n')

    # Add member functionals
    X = build_functional('PBE_X')
    X.set_alpha(1.0)
    sup.add_x_functional(X)
    C = build_functional('PBE_C')
    C.set_alpha(1.0)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.793701)
    sup.set_c_alpha(0.5)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_dsd_blyp_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('DSD-BLYP')
    # Tab in, trailing newlines
    sup.set_description('    DSD-BLYP Dispersion-corrected SCS Double Hybrid XC Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    S. Kozuch, Phys. Chem. Chem. Phys., 13, 20104, 2011\n')

    # Add member functionals
    X = build_functional('B88_X')
    X.set_alpha(1.0)
    sup.add_x_functional(X)
    C = build_functional('LYP_C')
    C.set_alpha(0.55)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.71)
    sup.set_c_alpha(1.0)
    sup.set_c_os_alpha(0.46)
    sup.set_c_ss_alpha(0.43)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_dsd_pbep86_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('DSD-PBEP86')
    # Tab in, trailing newlines
    sup.set_description('    DSD-PBEP86 Dispersion-corrected SCS Double Hybrid XC Functional (opt. for -D2)\n')
    # Tab in, trailing newlines
    sup.set_citation('    S. Kozuch, Phys. Chem. Chem. Phys., 13, 20104, 2011\n')

    # Add member functionals
    X = build_functional('PBE_X')
    X.set_alpha(1.0)
    sup.add_x_functional(X)
    C = build_functional('P86_C')
    C.set_alpha(0.45)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.68)
    sup.set_c_alpha(1.0)
    sup.set_c_ss_alpha(0.23)
    sup.set_c_os_alpha(0.51)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_dsd_pbepbe_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('DSD-PBEPBE')
    # Tab in, trailing newlines
    sup.set_description('    DSD-PBEPBE Dispersion-corrected SCS Double Hybrid XC Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    S. Kozuch, Phys. Chem. Chem. Phys., 13, 20104, 2011\n')

    # Add member functionals
    X = build_functional('PBE_X')
    X.set_alpha(1.0)
    sup.add_x_functional(X)
    C = build_functional('PBE_C')
    C.set_alpha(0.51)
    sup.add_c_functional(C)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.66)
    sup.set_c_alpha(1.0)
    sup.set_c_ss_alpha(0.12)
    sup.set_c_os_alpha(0.53)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_primitive_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    key = name.upper()
    fun = build_functional(key)

    # No spaces, keep it short and according to convention
    sup.set_name(key)
    # Tab in, trailing newlines
    sup.set_description(fun.description())
    # Tab in, trailing newlines
    sup.set_citation(fun.citation())

    # Add member functionals

    if (key[-1] == 'X'):
        sup.add_x_functional(fun)
    else:
        sup.add_c_functional(fun)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


def build_hf_superfunctional(name, npoints, deriv):
    # Special "functional" that is simply Hartree Fock

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('HF')
    # Tab in, trailing newlines
    sup.set_description('    Hartree Fock as Roothan prescribed\n')
    # Tab in, trailing newlines
    sup.set_citation('    \n')

    # 100% exact exchange
    sup.set_x_alpha(1.0)

    # Zero out other GKS
    sup.set_c_omega(0.0)
    sup.set_x_omega(0.0)
    sup.set_c_alpha(0.0)

    # Dont allocate, no functionals
    return (sup, False)


def build_hf3c_superfunctional(name, npoints, deriv):

    # call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => user-customization <= #

    # no spaces, keep it short and according to convention
    sup.set_name('HF3C')
    # tab in, trailing newlines
    sup.set_description('    Hartree Fock as Roothan prescribed plus 3C\n')
    # tab in, trailing newlines
    sup.set_citation('    Sure et al., J. Comput. Chem., 34, 1672-1685, 2013\n')

    # 100% exact exchange
    sup.set_x_alpha(1.0)

    # Zero out other GKS
    sup.set_c_omega(0.0)
    sup.set_x_omega(0.0)
    sup.set_c_alpha(0.0)

    # Dont allocate, no functionals
    return (sup, False)


# Superfunctional lookup table
superfunctionals = {
        'hf'              : build_hf_superfunctional,
        'hf+d'            : build_hfd_superfunctional,
        's_x'             : build_primitive_superfunctional,
        'b88_x'           : build_primitive_superfunctional,
        'b3_x'            : build_primitive_superfunctional,
        'pbe_x'           : build_primitive_superfunctional,
        'rpbe_x'          : build_primitive_superfunctional,
        'sogga_x'         : build_primitive_superfunctional,
        'pbesol_x'        : build_primitive_superfunctional,
        'pw91_x'          : build_primitive_superfunctional,
        'ws_x'            : build_ws_x_superfunctional,
        'wpbe_x'          : build_wpbe_x_superfunctional,
        'wpbesol_x'       : build_wpbesol_x_superfunctional,
        'wb88_x'          : build_wb88_x_superfunctional,
        'lyp_c'           : build_primitive_superfunctional,
        'ft97b_x'         : build_primitive_superfunctional,
        'pz81_c'          : build_primitive_superfunctional,
        'p86_c'           : build_primitive_superfunctional,
        'pw91_c'          : build_primitive_superfunctional,
        'pw92_c'          : build_primitive_superfunctional,
        'pbe_c'           : build_primitive_superfunctional,
        'ft97_c'          : build_primitive_superfunctional,
        'vwn5rpa_c'       : build_primitive_superfunctional,
        'vwn5_c'          : build_primitive_superfunctional,
        'vwn3rpa_c'       : build_primitive_superfunctional,
        'vwn3_c'          : build_primitive_superfunctional,
        'svwn'            : build_svwn_superfunctional,
        'blyp'            : build_blyp_superfunctional,
        'b86bpbe'         : build_b86bpbe_superfunctional,
        'pw86pbe'         : build_pw86pbe_superfunctional,
        'bp86'            : build_bp86_superfunctional,
        'pw91'            : build_pw91_superfunctional,
        'pbe'             : build_pbe_superfunctional,
        'hf3c'            : build_hf3c_superfunctional,
        'pbeh3c'          : build_pbeh3c_superfunctional,
        'ft97'            : build_ft97_superfunctional,
        'b3lyp'           : build_b3lyp_superfunctional,
        'b3lyp5'          : build_b3lyp5_superfunctional,
        'hf_x'            : build_hf_x_superfunctional,
        'pbe0'            : build_pbe0_superfunctional,
        'b97-0'           : build_b970_superfunctional,
        'b97-1'           : build_b971_superfunctional,
        'b97-2'           : build_b972_superfunctional,
        'b97-d'           : build_b97d_superfunctional,
        'hcth'            : build_hcth_superfunctional,
        'hcth120'         : build_hcth120_superfunctional,
        'hcth147'         : build_hcth147_superfunctional,
        'hcth407'         : build_hcth407_superfunctional,
        'wsvwn'           : build_wsvwn_superfunctional,
        'wpbe'            : build_wpbe_superfunctional,
        'wpbe0'           : build_wpbe0_superfunctional,
        'wpbesol'         : build_wpbesol_superfunctional,
        'wpbesol0'        : build_wpbesol0_superfunctional,
        'wblyp'           : build_wblyp_superfunctional,
        'wb97'            : build_wb97_superfunctional,
        'wb97x'           : build_wb97x_superfunctional,
        'wb97x-d'         : build_wb97xd_superfunctional,
        'm05'             : build_m05_superfunctional,
        'm05-2x'          : build_m05_2x_superfunctional,
        'dldf'            : build_dldf_superfunctional,
        'dldf+d09'        : build_dldfd09_superfunctional,
        'dldf+d'          : build_dldfd10_superfunctional,
        'sogga'           : build_sogga_superfunctional,
        'b2plyp'          : build_b2plyp_superfunctional,
        #'wb97x-2(tqz)'    : build_wb97x_2tqz_superfunctional,  # removed 26 Feb 2014 pending better handling of SS/OS DH coeff
        #'wb97x-2(lp)'     : build_wb97x_2lp_superfunctional,  # removed 26 Feb 2014 pending better handling of SS/OS DH coeff
        'pbe0-2'          : build_pbe0_2_superfunctional,
        #'dsd-blyp'        : build_dsd_blyp_superfunctional,  # -D variants still need to be added  # removed 26 Feb 2014 pending better handling of SS/OS DH coeff
        #'dsd-pbep86'      : build_dsd_pbep86_superfunctional,  # removed 26 Feb 2014 pending better handling of SS/OS DH coeff
        #'dsd-pbepbe'      : build_dsd_pbepbe_superfunctional,  # removed 26 Feb 2014 pending better handling of SS/OS DH coeff
        'pbea_c'          : build_primitive_superfunctional,
        'pw92a_c'         : build_primitive_superfunctional,
        'wpbe_c'          : build_wpbe_c_superfunctional,
        'wpw92_c'         : build_wpw92_c_superfunctional,
        'wpbe2'           : build_wpbe2_superfunctional,
    }

## Build up the lost of functionals we can compute
# Add in plain values
superfunctional_list = []
for key in superfunctionals.keys():
    sup = superfunctionals[key](key, 1, 1)[0]
    superfunctional_list.append(sup)

# Figure out what Grimme functionals we have
p4_funcs = set(superfunctionals.keys())
p4_funcs -= set(['b97-d'])
p4_funcs -= set(['hf3c'])
p4_funcs -= set(['pbeh3c'])
for dashlvl, superfunctional_listues in dftd3.dashcoeff.items():
    func_list = (set(superfunctional_listues.keys()) & p4_funcs)
    for func in func_list:
        sup = superfunctionals[func](func, 1, 1)[0]
        sup.set_name(sup.name() + '-' + dashlvl.upper())
        superfunctional_list.append(sup)

        if dashlvl == 'd2p4':
            # -D2 overide
            sup = superfunctionals[func](func, 1, 1)[0]
            sup.set_name(sup.name() + '-D2')
            superfunctional_list.append(sup)

            # -D overide
            sup = superfunctionals[func](func, 1, 1)[0]
            sup.set_name(sup.name() + '-D')
            superfunctional_list.append(sup)

        if dashlvl == 'd3zero':
            sup = superfunctionals[func](func, 1, 1)[0]
            sup.set_name(sup.name() + '-D3')
            superfunctional_list.append(sup)

        if dashlvl == 'd3mzero':
            sup = superfunctionals[func](func, 1, 1)[0]
            sup.set_name(sup.name() + '-D3M')
            superfunctional_list.append(sup)

# B97D is an odd one
for dashlvl in dftd3.full_dash_keys:
    if dashlvl == 'd2p4': continue

    sup = superfunctionals['b97-d']('b97-d', 1, 1)[0]
    sup.set_name('B97-' + dashlvl.upper())
    superfunctional_list.append(sup)

# wPBE, grr need a new scheme
for dashlvl in ['d3', 'd3m', 'd3zero', 'd3mzero', 'd3bj', 'd3mbj']:
    sup = superfunctionals['wpbe']('wpbe', 1, 1)[0]
    sup.set_name(sup.name() + '-' + dashlvl.upper())
    superfunctional_list.append(sup)




def build_superfunctional(alias):
    name = alias.lower()

    npoints = core.get_option("SCF", "DFT_BLOCK_MAX_POINTS");
    deriv = 1 # Default depth for now

    # Grab out superfunctional
    if name in ["gen", ""]:
        sup = (core.get_option("DFT_CUSTOM_FUNCTIONAL"), False)
        if not isinstance(sup[0], core.SuperFunctional):
            raise KeyError("SCF: Custom Functional requested, but nothing provided in DFT_CUSTOM_FUNCTIONAL")

    elif name in ['hf3c', 'pbeh3c']:
        func = superfunctionals[name](name, npoints, deriv)[0]
        sup = (func, (name, 'd3bj'))

    elif name in superfunctionals.keys():
        sup = superfunctionals[name](name, npoints, deriv)


    elif any(name.endswith(al) for al in dftd3.full_dash_keys):

        # Odd hack for b97-d
        if 'b97-d' in name:
            name = name.replace('b97', 'b97-d')

        dashparam = [x for x in dftd3.full_dash_keys if name.endswith(x)]
        if len(dashparam) > 1:
            raise Exception("Dashparam %s is ambiguous.")
        else:
            dashparam = dashparam[0]

        base_name = name.replace('-' + dashparam, '')

        if dashparam in ['d2', 'd']:
            dashparam = 'd2p4'

        if dashparam == 'd3':
            dashparam = 'd3zero'

        if dashparam == 'd3m':
            dashparam = 'd3mzero'

        if base_name not in superfunctionals.keys():
            raise KeyError("SCF: Functional (%s) with base (%s) not found!" % (alias, base_name))

        func = superfunctionals[base_name](base_name, npoints, deriv)[0]

        base_name = base_name.replace('wpbe', 'lcwpbe')
        sup = (func, (base_name, dashparam))

    else:
        raise KeyError("SCF: Functional (%s) not found!" % alias)

    # Set options
    if core.has_option_changed("SCF", "DFT_OMEGA") and sup[0].is_x_lrc():
        sup[0].set_x_omega(core.get_option("SCF", "DFT_OMEGA"))
    if core.has_option_changed("SCF", "DFT_OMEGA_C") and sup[0].is_c_lrc():
        sup[0].set_c_omega(core.get_option("SCF", "DFT_OMEGA_C"))

    if core.has_option_changed("SCF", "DFT_ALPHA"):
        sup[0].set_x_alpha(core.get_option("SCF", "DFT_ALPHA"))
    if core.has_option_changed("SCF", "DFT_ALPHA_C"):
        sup[0].set_c_alpha(core.get_option("SCF", "DFT_ALPHA_C"))

    # Check SCF_TYPE
    if sup[0].is_x_lrc() and (core.get_option("SCF", "SCF_TYPE") not in ["DIRECT", "DF", "OUT_OF_CORE", "PK"]):
        raise KeyError("SCF: SCF_TYPE (%s) not supported for range-seperated functionals."
                        % core.get_option("SCF", "SCF_TYPE"))

    if (core.get_global_option('INTEGRAL_PACKAGE') == 'ERD') and (sup[0].is_x_lrc()):
        raise ValidationError('INTEGRAL_PACKAGE ERD does not play nicely with LRC DFT functionals, so stopping.')

    return sup

def test_ccl_functional(functional, ccl_functional):

    check = True

    if (not os.path.exists('data_pt_%s.html' % (ccl_functional))):
        os.system('wget ftp://ftp.dl.ac.uk/qcg/dft_library/data_pt_%s.html' % ccl_functional)
    fh = open('data_pt_%s.html' % (ccl_functional))
    lines = fh.readlines()
    fh.close()

    points = []
    point = {}

    rho_line = re.compile(r'^\s*rhoa=\s*(-?\d+\.\d+E[+-]\d+)\s*rhob=\s*(-?\d+\.\d+E[+-]\d+)\s*sigmaaa=\s*(-?\d+\.\d+E[+-]\d+)\s*sigmaab=\s*(-?\d+\.\d+E[+-]\d+)\s*sigmabb=\s*(-?\d+\.\d+E[+-]\d+)\s*')
    val_line = re.compile(r'^\s*(\w*)\s*=\s*(-?\d+\.\d+E[+-]\d+)')

    aliases = { 'zk'            : 'v',
                'vrhoa'         : 'v_rho_a',
                'vrhob'         : 'v_rho_b',
                'vsigmaaa'      : 'v_gamma_aa',
                'vsigmaab'      : 'v_gamma_ab',
                'vsigmabb'      : 'v_gamma_bb',
                'v2rhoa2'       : 'v_rho_a_rho_a',
                'v2rhoab'       : 'v_rho_a_rho_b',
                'v2rhob2'       : 'v_rho_b_rho_b',
                'v2rhoasigmaaa' : 'v_rho_a_gamma_aa',
                'v2rhoasigmaab' : 'v_rho_a_gamma_ab',
                'v2rhoasigmabb' : 'v_rho_a_gamma_bb',
                'v2rhobsigmaaa' : 'v_rho_b_gamma_aa',
                'v2rhobsigmaab' : 'v_rho_b_gamma_ab',
                'v2rhobsigmabb' : 'v_rho_b_gamma_bb',
                'v2sigmaaa2'    : 'v_gamma_aa_gamma_aa',
                'v2sigmaaaab'   : 'v_gamma_aa_gamma_ab',
                'v2sigmaaabb'   : 'v_gamma_aa_gamma_bb',
                'v2sigmaab2'    : 'v_gamma_ab_gamma_ab',
                'v2sigmaabbb'   : 'v_gamma_ab_gamma_bb',
                'v2sigmabb2'    : 'v_gamma_bb_gamma_bb',
              }

    for line in lines:

        mobj = re.match(rho_line, line)
        if (mobj):

            if len(point):
                points.append(point)
                point = {}

            point['rho_a'] = float(mobj.group(1))
            point['rho_b'] = float(mobj.group(2))
            point['gamma_aa'] = float(mobj.group(3))
            point['gamma_ab'] = float(mobj.group(4))
            point['gamma_bb'] = float(mobj.group(5))

            continue

        mobj = re.match(val_line, line)
        if (mobj):
            point[aliases[mobj.group(1)]] = float(mobj.group(2))

    points.append(point)

    N = len(points)
    rho_a = core.Vector(N)
    rho_b = core.Vector(N)
    gamma_aa = core.Vector(N)
    gamma_ab = core.Vector(N)
    gamma_bb = core.Vector(N)
    tau_a = core.Vector(N)
    tau_b = core.Vector(N)

    index = 0
    for point in points:
        rho_a[index] = point['rho_a']
        rho_b[index] = point['rho_b']
        gamma_aa[index] = point['gamma_aa']
        gamma_ab[index] = point['gamma_ab']
        gamma_bb[index] = point['gamma_bb']
        index = index + 1

    super = build_superfunctional(functional, N, 1)
    super.test_functional(rho_a, rho_b, gamma_aa, gamma_ab, gamma_bb, tau_a, tau_b)

    v = super.value('V')
    v_rho_a = super.value('V_RHO_A')
    v_rho_b = super.value('V_RHO_B')
    v_gamma_aa = super.value('V_GAMMA_AA')
    v_gamma_ab = super.value('V_GAMMA_AB')
    v_gamma_bb = super.value('V_GAMMA_BB')

    if not v_gamma_aa:
        v_gamma_aa = tau_a
        v_gamma_ab = tau_a
        v_gamma_bb = tau_a

    tasks = ['v', 'v_rho_a', 'v_rho_b', 'v_gamma_aa', 'v_gamma_ab', 'v_gamma_bb']
    mapping = {
            'v': v,
            'v_rho_a': v_rho_a,
            'v_rho_b': v_rho_b,
            'v_gamma_aa': v_gamma_aa,
            'v_gamma_ab': v_gamma_ab,
            'v_gamma_bb': v_gamma_bb,
        }

    super.print_detail(3)
    index = 0
    for point in points:
        core.print_out('rho_a= %11.3E, rho_b= %11.3E, gamma_aa= %11.3E, gamma_ab= %11.3E, gamma_bb= %11.3E\n' % (rho_a[index], rho_b[index], gamma_aa[index], gamma_ab[index], gamma_bb[index]))

        for task in tasks:
            v_ref = point[task]
            v_obs = mapping[task][index]
            delta = v_obs - v_ref
            if (v_ref == 0.0):
                epsilon = 0.0
            else:
                epsilon = abs(delta / v_ref)
            if (epsilon < 1.0E-11):
                passed = 'PASSED'
            else:
                passed = 'FAILED'
                check = False

            core.print_out('\t%-15s %24.16E %24.16E %24.16E %24.16E %6s\n' % (task, v_ref, v_obs, delta, epsilon, passed))

        index = index + 1

    core.print_out('\n')
    return check
