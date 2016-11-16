  #
  # @BEGIN LICENSE
  #
  # Psi4: an open-source quantum chemistry software package
  #
  # Copyright (c) 2007-2016 The Psi4 Developers.
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
List of GGA SuperFunctionals built from LibXC primitives.
"""

from psi4 import core

def build_svwn_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('SVWN')
    sup.set_description('    SVWN3 (RPA) LSDA Functional\n')
    sup.set_citation('    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n')

    # Add member functionals
    sup.add_x_functional(core.Functional.build_base('XC_LDA_X'))
    sup.add_c_functional(core.Functional.build_base('XC_LDA_C_VWN_RPA'))

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
    sup.set_description('    BLYP GGA Exchange-Correlation Functional\n')
    sup.set_citation('    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n    B. Miehlich et. al., Chem. Phys. Lett., 157(3), 200-206 1989\n')

    # Add member functionals
    sup.add_x_functional(core.Functional.build_base('XC_GGA_X_B88'))
    sup.add_c_functional(core.Functional.build_base('XC_GGA_C_LYP'))

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
    sup.set_description('    B86BPBE GGA Exchange-Correlation Functional\n')
    sup.set_citation('    A. D. Becke, J. Chem. Phys. 85:7184, 1986.\n')

    # Add member functionals
    sup.add_x_functional(core.Functional.build_base('XC_GGA_X_B86_MGC'))
    sup.add_c_functional(core.Functional.build_base('XC_GGA_C_PBE'))

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
    sup.set_description('    PW86PBE GGA Exchange-Correlation Functional\n')
    sup.set_citation('    J. P. Perdew and W. Yue, Phys. Rev. B 33:8800(R), 1986.\n')

    # Add member functionals
    sup.add_x_functional(core.Functional.build_base('XC_GGA_X_PW86'))
    sup.add_c_functional(core.Functional.build_base('XC_GGA_C_PBE'))

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
    sup.set_description('    PBE GGA Exchange-Correlation Functional\n')
    sup.set_citation('    J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996\n')

    # Add member functionals
    sup.add_x_functional(core.Functional.build_base('XC_GGA_X_PBE'))
    sup.add_c_functional(core.Functional.build_base('XC_GGA_C_PBE'))

    # Call this last
    sup.allocate()
    return (sup, False)

# def build_wsvwn_superfunctional(name, npoints, deriv):

#     # Call this first
#     sup = core.SuperFunctional.blank()
#     sup.set_max_points(npoints)
#     sup.set_deriv(deriv)

#     # => User-Customization <= #

#     # No spaces, keep it short and according to convention
#     sup.set_name('wSVWN')
#     sup.set_description('    LSDA SR-XC Functional\n')
#     sup.set_citation('    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n')

#     # Add member functionals
#     wS_X = core.Functional.build_base('wS_X')
#     sup.add_x_functional(wS_X)
#     sup.add_c_functional(core.Functional.build_base('VWN3RPA_C'))

#     # Set GKS up after adding functionals
#     sup.set_x_omega(0.3)
#     sup.set_c_omega(0.0)
#     sup.set_x_alpha(0.0)
#     sup.set_c_alpha(0.0)

#     # => End User-Customization <= #

#     # Call this last
#     sup.allocate()
#     return (sup, False)

def build_pw91_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('PW91')
    sup.set_description('    PW91 GGA Exchange-Correlation Functional\n')
    sup.set_citation('    J.P. Perdew et. al., Phys. Rev. B., 46(11), 6671-6687, 1992\n')

    # Add member functionals
    sup.add_x_functional(core.Functional.build_base('XC_GGA_X_PW91'))
    sup.add_c_functional(core.Functional.build_base('XC_GGA_C_PW91'))

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
    sup.set_description('    BP86 GGA Exchange-Correlation Functional\n')
    sup.set_citation('    Null\n')

    # Add member functionals
    sup.add_x_functional(core.Functional.build_base('XC_GGA_X_B88'))
    sup.add_c_functional(core.Functional.build_base('XC_GGA_C_P86'))

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
    sup.set_description('   FT97 GGA Exchange-Correlation Functional\n')
    sup.set_citation('    M. Filatov and W. Theil, Int. J. Quant. Chem., 62, 603-616, 1997\n')

    # Add member functionals
    sup.add_x_functional(core.Functional.build_base('XC_GGA_X_FT97_B'))
    sup.add_c_functional(core.Functional.build_base('XC_GGA_C_FT97'))


    # Call this last
    sup.allocate()
    return (sup, False)


gga_superfunc_list = {
          "b86bpbe" : build_b86bpbe_superfunctional,
          "blyp"    : build_blyp_superfunctional,
          "svwn"    : build_svwn_superfunctional,
          "pw86pbe" : build_pw86pbe_superfunctional,
          "pbe"     : build_pbe_superfunctional,
          "bp86"    : build_bp86_superfunctional,
          "pw91"    : build_pw91_superfunctional,
          "ft97"    : build_ft97_superfunctional,
}
