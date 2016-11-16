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

def build_pbe0_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('PBE0')
    # Tab in, trailing newlines
    sup.set_description('    PBE GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('    J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996\n')

    # Add member functionals
    pbe_x = core.Functional.build_base('XC_GGA_X_PBE')
    pbe_x.set_alpha(0.75)
    sup.add_x_functional(pbe_x)
    sup.add_c_functional(core.Functional.build_base('XC_GGA_C_PBE'))

    sup.set_x_alpha(0.25)

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
    pbe_x = core.Functional.build_base('wPBE_X')
    pbe_x.set_omega(0.4)
    sup.add_x_functional(pbe_x)
    sup.add_c_functional(core.Functional.build_base('XC_GGA_C_PBE'))

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
    pbe_x = core.Functional.build_base('wPBE_X')
    pbe_x.set_omega(0.3)
    pbe_x.set_alpha(0.75)
    sup.add_x_functional(pbe_x)
    sup.add_c_functional(core.Functional.build_base('XC_GGA_C_PBE'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.25)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)

def build_wb97xd_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.XC_build("XC_HYB_GGA_XC_WB97X_D")
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


hyb_gga_superfunc_list = {
          "pbe0"    : build_pbe0_superfunctional,
          "wpbe"    : build_wpbe_superfunctional,
          "wpbe0"   : build_wpbe0_superfunctional,
          "wb97x-d" : build_wb97xd_superfunctional,

}
