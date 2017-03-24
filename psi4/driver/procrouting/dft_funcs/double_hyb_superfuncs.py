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


def build_b2plyp_superfunctional(name, npoints, deriv, restricted):

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
    becke = core.LibXCFunctional('XC_GGA_X_B88', restricted)
    becke.set_alpha(0.47)
    sup.add_x_functional(becke)
    lyp = core.LibXCFunctional('XC_GGA_C_LYP', restricted)
    lyp.set_alpha(0.73)
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



def build_dsd_blyp_superfunctional(name, npoints, deriv, restricted):

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
    X = core.LibXCFunctional('XC_GGA_X_B88', restricted)
    X.set_alpha(0.29)
    sup.add_x_functional(X)
    C = core.LibXCFunctional('XC_GGA_C_LYP', restricted)
    C.set_alpha(0.55) #  Fix this!
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


def build_pbe0_2_superfunctional(name, npoints, deriv, restricted):

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
    X = core.LibXCFunctional('XC_GGA_X_PBE', restricted)
    X.set_alpha(0.206299)
    sup.add_x_functional(X)
    C = core.LibXCFunctional('XC_GGA_C_PBE', restricted)
    C.set_alpha(0.5)
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

def build_dsd_pbep86_superfunctional(name, npoints, deriv, restricted):

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
    X = core.LibXCFunctional('XC_GGA_X_PBE', restricted)
    X.set_alpha(0.32)
    sup.add_x_functional(X)
    C = core.LibXCFunctional('XC_GGA_C_P86', restricted)
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


def build_dsd_pbepbe_superfunctional(name, npoints, deriv, restricted):

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
    X = core.LibXCFunctional('XC_GGA_X_PBE', restricted)
    X.set_alpha(0.34)
    sup.add_x_functional(X)
    C = core.LibXCFunctional('XC_GGA_C_PBE', restricted)
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


double_hyb_superfunc_list = {
          "b2plyp"     : build_b2plyp_superfunctional,
          "pbe0-2"     : build_pbe0_2_superfunctional,
          # "dsd-blyp"   : build_dsd_blyp_superfunctional,
          "dsd-pbep86" : build_dsd_pbep86_superfunctional,
          "dsd-pbepbe" : build_dsd_pbepbe_superfunctional,
}


# def build_wb97x_2tqz_superfunctional(name, npoints, deriv, restricted):

#     # Call this first
#     sup = core.SuperFunctional.blank()
#     sup.set_max_points(npoints)
#     sup.set_deriv(deriv)

#     # => User-Customization <= #

#     # No spaces, keep it short and according to convention
#     sup.set_name('wB97X-2(TQZ)')
#     # Tab in, trailing newlines
#     sup.set_description('    Double Hybrid LRC B97 GGA XC Functional (TQZ parametrization)\n')
#     # Tab in, trailing newlines
#     sup.set_citation('    J.-D. Chai and M. Head-Gordon, J. Chem. Phys., 131, 174105, 2009\n')

#     # Add member functionals
#     X = build_functional('wB97_X')
#     X.set_name('wB97X_X')
#     X.set_alpha(1.0 / (1.0 - 0.636158))

#     X.set_parameter('B97_gamma', 0.004)
#     X.set_parameter('B97_a0', 3.15503E-1)
#     X.set_parameter('B97_a1', 1.04772E0)
#     X.set_parameter('B97_a2', -2.33506E0)
#     X.set_parameter('B97_a3', 3.19909E0)

#     C = build_functional('B_C')
#     C.set_name('wB97X_C')

#     C.set_parameter('B97_os_gamma', 0.006)
#     C.set_parameter('B97_os_a0', 5.18198E-1)
#     C.set_parameter('B97_os_a1', -5.85956E-1)
#     C.set_parameter('B97_os_a2', 4.27080E0)
#     C.set_parameter('B97_os_a3', -6.48897E0)

#     C.set_parameter('B97_ss_gamma', 0.2)
#     C.set_parameter('B97_ss_a0', 9.08460E-1)
#     C.set_parameter('B97_ss_a1', -2.80936E0)
#     C.set_parameter('B97_ss_a2', 6.02676E0)
#     C.set_parameter('B97_ss_a3', -4.56981E0)

#     sup.add_x_functional(X)
#     sup.add_c_functional(C)

#     # Set GKS up after adding functionals
#     sup.set_x_omega(0.3)
#     sup.set_c_omega(0.0)
#     sup.set_x_alpha(0.636158)
#     sup.set_c_alpha(1.0)
#     sup.set_c_os_alpha(0.447105)
#     sup.set_c_ss_alpha(0.529319)

#     # => End User-Customization <= #

#     # Call this last
#     sup.allocate()
#     return (sup, False)


# def build_wb97x_2lp_superfunctional(name, npoints, deriv, restricted):

#     # Call this first
#     sup = core.SuperFunctional.blank()
#     sup.set_max_points(npoints)
#     sup.set_deriv(deriv)

#     # => User-Customization <= #

#     # No spaces, keep it short and according to convention
#     sup.set_name('wB97X-2(LP)')
#     # Tab in, trailing newlines
#     sup.set_description('    Double Hybrid LRC B97 GGA XC Functional (Large Pople parametrization)\n')
#     # Tab in, trailing newlines
#     sup.set_citation('    J.-D. Chai and M. Head-Gordon, J. Chem. Phys., 131, 174105, 2009\n')

#     # Add member functionals
#     X = build_functional('wB97_X')
#     X.set_name('wB97X_X')
#     X.set_alpha(1.0 / (1.0 - 0.678792))

#     X.set_parameter('B97_gamma', 0.004)
#     X.set_parameter('B97_a0', 2.51767E-1)
#     X.set_parameter('B97_a1', 1.57375E0)
#     X.set_parameter('B97_a2', -5.26624E0)
#     X.set_parameter('B97_a3', 6.74313E0)

#     C = build_functional('B_C')
#     C.set_name('wB97X_C')

#     C.set_parameter('B97_os_gamma', 0.006)
#     C.set_parameter('B97_os_a0', 5.53261E-1)
#     C.set_parameter('B97_os_a1', -1.16626E0)
#     C.set_parameter('B97_os_a2', 6.84409E0)
#     C.set_parameter('B97_os_a3', -8.90640E0)

#     C.set_parameter('B97_ss_gamma', 0.2)
#     C.set_parameter('B97_ss_a0', 1.15698E0)
#     C.set_parameter('B97_ss_a1', -3.31669E0)
#     C.set_parameter('B97_ss_a2', 6.27265E0)
#     C.set_parameter('B97_ss_a3', -4.51464E0)

#     sup.add_x_functional(X)
#     sup.add_c_functional(C)

#     # Set GKS up after adding functionals
#     sup.set_x_omega(0.3)
#     sup.set_c_omega(0.0)
#     sup.set_x_alpha(0.678792)
#     sup.set_c_alpha(1.0)
#     sup.set_c_os_alpha(0.477992)
#     sup.set_c_ss_alpha(0.581569)

#     # => End User-Customization <= #

#     # Call this last
#     sup.allocate()
#     return (sup, False)