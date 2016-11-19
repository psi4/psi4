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
List of MGGA SuperFunctionals built from LibXC primitives.
"""

from psi4 import core

def build_dldf_superfunctional(name, npoints, deriv):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('dlDF')
    sup.set_description('    Dispersionless Hybrid Meta-GGA XC Functional\n')
    sup.set_citation('    Pernal et. al., Phys. Rev. Lett., 103, 263201, 2009\n')

    # Add member functionals
    x_coef = 0.6144129
    dldf_x = core.Functional.build_base('XC_HYB_MGGA_X_DLDF')
    # dldf_x.set_alpha(1.0 - x_coef)
    sup.add_x_functional(dldf_x)
    sup.add_c_functional(core.Functional.build_base('XC_MGGA_C_DLDF'))

    sup.set_x_alpha(x_coef)

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


mgga_superfunc_list = {
          "dldf" : build_dldf_superfunctional,
          "dldf+d09" : build_dldfd09_superfunctional,
          "dldf+d" : build_dldfd10_superfunctional,
}
