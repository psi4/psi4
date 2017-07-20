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
List of GGA SuperFunctionals built from LibXC primitives.
"""

from psi4 import core

def build_pbe0_superfunctional(name, npoints, deriv, restricted):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('PBE0')
    # Tab in, trailing newlines
    sup.set_description('    PBE0 Hyb-GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    J.P. Perdew et. al., J. Chem. Phys., 105(22), 9982-9985, 1996\n    C. Adamo et. a., J. Chem Phys., 110(13), 6158-6170, 1999\n')

    # Add member functionals
    pbe_x = core.LibXCFunctional('XC_GGA_X_PBE', restricted)
    pbe_x.set_alpha(0.75)
    sup.add_x_functional(pbe_x)
    sup.add_c_functional(core.LibXCFunctional('XC_GGA_C_PBE', restricted))

    sup.set_x_alpha(0.25)

    # Call this last
    sup.allocate()
    return (sup, False)

def build_b5050lyp_superfunctional(name, npoints, deriv, restricted):
    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('B5050LYP')
    # Tab in, trailing newlines
    sup.set_description('    B5050LYP Hyb-GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    Y. Shao et. al., J. Chem. Phys., 188, 4807-4818, 2003\n')

    # Add member functionals
    slater = core.LibXCFunctional("XC_LDA_X", restricted)
    slater.set_alpha(0.08)
    sup.add_x_functional(slater)
    becke = core.LibXCFunctional("XC_GGA_X_B88", restricted)
    becke.set_alpha(0.42)
    sup.add_x_functional(becke)
    sup.set_x_alpha(0.50)

    vwn = core.LibXCFunctional("XC_LDA_C_VWN", restricted)
    vwn.set_alpha(0.19)
    sup.add_c_functional(vwn)
    lyp = core.LibXCFunctional("XC_GGA_C_LYP", restricted)
    lyp.set_alpha(0.81)
    sup.add_c_functional(lyp)

    # Call this last
    sup.allocate()
    return (sup, False)

def build_wpbe_superfunctional(name, npoints, deriv, restricted):

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
    pbe_x = core.LibXCFunctional('XC_GGA_X_HJS_PBE', restricted)
    pbe_x.set_omega(0.4)
    sup.add_x_functional(pbe_x)
    sup.add_c_functional(core.LibXCFunctional('XC_GGA_C_PBE', restricted))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.4)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)

def build_wpbe0_superfunctional(name, npoints, deriv, restricted):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wPBE0')
    # Tab in, trailing newlines
    sup.set_description('    PBE0 SR-XC Functional (HJS Model)\n')
    # Tab in, trailing newlines
    sup.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n')

    # Add member functionals
    pbe_x = core.LibXCFunctional('XC_GGA_X_HJS_PBE', restricted)
    pbe_x.set_omega(0.3)
    pbe_x.set_alpha(0.75)
    sup.add_x_functional(pbe_x)
    sup.add_c_functional(core.LibXCFunctional('XC_GGA_C_PBE', restricted))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.25)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)

def build_wb97xd_superfunctional(name, npoints, deriv, restricted):

    # Call this first
    sup = core.SuperFunctional.XC_build("XC_HYB_GGA_XC_WB97X_D", restricted)
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wB97X-D')
    # Tab in, trailing newlines
    sup.set_description('    Parameterized Hybrid LRC B97 GGA XC Functional with Dispersion\n')
    # Tab in, trailing newlines
    sup.set_citation('    J.-D. Chai and M. Head-Gordon, Phys. Chem. Chem. Phys., 10, 6615-6620, 2008\n')

    # Call this last
    sup.allocate()
    return (sup, ('wB97', '-CHG'))


def build_hfd_superfunctional(name, npoints, deriv, restricted):

    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)
    sup.set_name('HF+D')
    sup.set_x_alpha(1.0)

    sup.allocate()
    return (sup, ('HF', '-DAS2010'))

def build_hf_superfunctional(name, npoints, deriv, restricted):

    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)
    sup.set_name('HF')
    sup.set_x_alpha(1.0)

    sup.allocate()
    return (sup, False)

def build_hf3c_superfunctional(name, npoints, deriv, restricted):

    sup = core.SuperFunctional.blank()

    # => user-customization <= #

    # no spaces, keep it short and according to convention
    sup.set_name('HF3C')
    # tab in, trailing newlines
    sup.set_description('    Hartree Fock as Roothaan prescribed plus 3C\n')
    # tab in, trailing newlines
    sup.set_citation('    Sure et al., J. Comput. Chem., 34, 1672-1685, 2013\n')

    sup.set_max_points(npoints)
    sup.set_deriv(deriv)
    sup.set_x_alpha(1.0)

    sup.allocate()
    return (sup, ('HF3C', '-d3bj'))

def build_pbeh3c_superfunctional(name, npoints, deriv, restricted):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('PBEH3C')
    # tab in, trailing newlines
    sup.set_description('    PBEH-3C Hybrid GGA Exchange-Correlation Functional plus 3C\n')
    # tab in, trailing newlines
    sup.set_citation('    Grimme et. al., J. Chem. Phys., 143, 054107, 2015\n')

    # Add member functionals
    pbe_x3c = core.LibXCFunctional('XC_GGA_X_PBE', restricted)
    pbe_x3c.set_alpha(0.58)
    pbe_x3c.set_tweak([1.0245, 0.12345679])
    sup.add_x_functional(pbe_x3c)

    pbe_c3c = core.LibXCFunctional('XC_GGA_C_PBE', restricted)
    pbe_c3c.set_tweak([0.03])
    sup.add_c_functional(pbe_c3c)

    # set gks up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.42)
    sup.set_c_alpha(0.0)

    # Call this last
    sup.allocate()
    return (sup, ("PBEH3C", "-d3bj"))


hyb_superfunc_list = {
          "pbeh3c"   : build_pbeh3c_superfunctional,
          "pbe0"     : build_pbe0_superfunctional,
          "wpbe"     : build_wpbe_superfunctional,
          "wpbe0"    : build_wpbe0_superfunctional,
          "b5050lyp" : build_b5050lyp_superfunctional,
          "wb97x-d"  : build_wb97xd_superfunctional,
          "hf-d"     : build_hfd_superfunctional,
          "hf"       : build_hf_superfunctional,
          "hf3c"     : build_hf3c_superfunctional,

}





# def build_wpbesol_superfunctional(name, npoints, deriv, restricted):

#     # Call this first
#     sup = core.SuperFunctional.blank()
#     sup.set_max_points(npoints)
#     sup.set_deriv(deriv)

#     # => User-Customization <= #

#     # No spaces, keep it short and according to convention
#     sup.set_name('wPBEsol')
#     # Tab in, trailing newlines
#     sup.set_description('    PBEsol SR-XC Functional (HJS Model)\n')
#     # Tab in, trailing newlines
#     sup.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n')

#     # Add member functionals
#     sup.add_x_functional(build_functional('wPBEsol_X'))
#     sup.add_c_functional(build_functional('PBE_C'))

#     # Set GKS up after adding functionals
#     sup.set_x_omega(0.4)
#     sup.set_c_omega(0.0)
#     sup.set_x_alpha(0.0)
#     sup.set_c_alpha(0.0)

#     # => End User-Customization <= #

#     # Call this last
#     sup.allocate()
#     return (sup, False)


# def build_wpbesol0_superfunctional(name, npoints, deriv, restricted):

#     sup = build_wpbesol_superfunctional(name, npoints, deriv)[0]
#     sup.set_name('wPBEsol0')
#     sup.set_description('    PBEsol0 SR-XC Functional (HJS Model)\n')
#     sup.set_x_omega(0.3)
#     sup.set_x_alpha(0.25)
#     return (sup, False)
