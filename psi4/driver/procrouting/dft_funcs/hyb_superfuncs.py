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
    sup.set_citation(
        '    J.P. Perdew et. al., J. Chem. Phys., 105(22), 9982-9985, 1996\n    C. Adamo et. a., J. Chem Phys., 110(13), 6158-6170, 1999\n'
    )

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
    sup.set_citation(
        '    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754, 2009\n'
    )

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
    sup.set_citation(
        '    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754, 2009\n'
    )

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
    sup.set_name('SCF')
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


def build_sogga11_x_superfunctional(name, npoints, deriv, restricted):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('SOGGA11-X')
    sup.set_description('   SOGGA11-X Hybrid Exchange-Correlation Functional\n')
    sup.set_citation('    R. Peverati and D. G. Truhlar, J. Chem. Phys. 135, 191102, 2011\n')

    # Add member functionals
    sup.add_x_functional(core.LibXCFunctional('XC_HYB_GGA_X_SOGGA11_X', restricted))
    sup.add_c_functional(core.LibXCFunctional('XC_GGA_C_SOGGA11_X', restricted))
    sup.set_x_alpha(0.4015)

    # Call this last
    sup.allocate()
    return (sup, False)


def build_scan0_superfunctional(name, npoints, deriv, restricted):
    # Disabled below, no SCAN correlation in LibXC 3.0.0
    # SCAN correlation results unreliable.

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('SCAN0')
    sup.set_description('   SCAN0 Hybrid Exchange-Correlation Functional\n')
    sup.set_citation('    K. Hui and J.-D. Chai, J. Chem. Phys. 144, 044114, 2016\n')

    # Add member functionals
    sup.add_x_functional(core.LibXCFunctional('XC_HYB_MGGA_X_SCAN0', restricted))
    sup.add_c_functional(core.LibXCFunctional('XC_MGGA_C_SCAN', restricted))
    sup.set_x_alpha(0.25)

    # Call this last
    sup.allocate()
    return (sup, False)


def build_n12_sx_superfunctional(name, npoints, deriv, restricted):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('N12-SX')
    sup.set_description('   N12-SX Hybrid Screened Exchange-Correlation Functional\n')
    sup.set_citation('    R. Peverati, D. G. Truhlar, Phys. Chem. Chem. Phys 14, 16187, 2012\n')

    # Add member functionals
    sup.add_x_functional(core.LibXCFunctional('XC_HYB_GGA_X_N12_SX', restricted))
    sup.add_c_functional(core.LibXCFunctional('XC_GGA_C_N12_SX', restricted))
    sup.set_x_alpha(0.25)
    sup.set_x_beta(-0.25)
    sup.set_x_omega(0.11)
    

    # Call this last
    sup.allocate()
    return (sup, False)


def build_mn12_sx_superfunctional(name, npoints, deriv, restricted):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('MN12-SX')
    sup.set_description('   MN12-SX Meta-GGA Hybrid Screened Exchange-Correlation Functional\n')
    sup.set_citation('    R. Peverati, D. G. Truhlar, Phys. Chem. Chem. Phys 14, 16187, 2012\n')

    # Add member functionals
    sup.add_x_functional(core.LibXCFunctional('XC_HYB_MGGA_X_MN12_SX', restricted))
    sup.add_c_functional(core.LibXCFunctional('XC_MGGA_C_MN12_SX', restricted))
    sup.set_x_alpha(0.25)
    sup.set_x_beta(-0.25)
    sup.set_x_omega(0.11)
    

    # Call this last
    sup.allocate()
    return (sup, False)    


def build_mn15_superfunctional(name, npoints, deriv, restricted):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('MN15')
    sup.set_description('   MN15 Hybrid Exchange-Correlation Functional\n')
    sup.set_citation('    H. S. Yu, X. He, S. L. Li, and D. G. Truhlar, Chem. Sci. 7, 5032-5051, 2016\n')

    # Add member functionals
    sup.add_x_functional(core.LibXCFunctional('XC_HYB_MGGA_X_MN15', restricted))
    sup.add_c_functional(core.LibXCFunctional('XC_MGGA_C_MN15', restricted))
    sup.set_x_alpha(0.44)

    # Call this last
    sup.allocate()
    return (sup, False)


def build_pw6b95_superfunctional(name, npoints, deriv, restricted):

    # Call this first
    sup = core.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('pw6b95')
    # Tab in, trailing newlines
    sup.set_description('    PW6B95 Hybrid-meta XC Functional\n')
    # Tab in, trailing newlines
    sup.set_citation('  Y. Zhao and D. Truhlar, J. Phys. Chem. A., 109,5656-5667, 2005\n')
    #PW6B95(hybrid) 0.00538      1.7382  3.8901  0.00262 0.03668 0.28
    # Add member functionals
    pw6 = core.LibXCFunctional('XC_GGA_X_PW91', restricted)
    # modify PW91 suitable for libxc b=1/X2S is unchanged
    #  a    =  6.0*bt/X2S;
    #  b    =  1.0/X2S;
    #  c    =  bt/(X_FACTOR_C*X2S*X2S);
    #  d    = -(bt - beta)/(X_FACTOR_C*X2S*X2S);
    #  f    = 1.0e-6/(X_FACTOR_C*POW(X2S, expo));
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
    pw6.set_tweak([a_pw6, b_pw6, c_pw6, d_pw6, f_pw6, alpha_pw6, expo_pw6])
    pw6.set_alpha(0.72)
    sup.add_x_functional(pw6)

    mb95 = core.LibXCFunctional('XC_MGGA_C_BC95', restricted)
    copp = 0.00262
    css = 0.03668
    mb95.set_tweak([css, copp])
    sup.add_c_functional(mb95)

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.28)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return (sup, False)


hyb_superfunc_list = {
    "pbeh3c": build_pbeh3c_superfunctional,
    "pbe0": build_pbe0_superfunctional,
    "wpbe": build_wpbe_superfunctional,
    "wpbe0": build_wpbe0_superfunctional,
    "b5050lyp": build_b5050lyp_superfunctional,
    "wb97x-d": build_wb97xd_superfunctional,
    #          "wb97x-d3" : build_wb97xd3_superfunctional,
    "hf-d": build_hfd_superfunctional,
    "hf": build_hf_superfunctional,
    "scf": build_hf_superfunctional,
    "hf3c": build_hf3c_superfunctional,
    #    "scan0": build_scan0_superfunctional, # XC_MGGA_C_SCAN not present in LibXC 3.0.0
    "sogga11-x": build_sogga11_x_superfunctional,
    "n12-sx": build_n12_sx_superfunctional,
    "mn12-sx": build_mn12_sx_superfunctional,
    "mn15": build_mn15_superfunctional,
    "pw6b95": build_pw6b95_superfunctional,
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

#def build_wb97xd3_superfunctional(name, npoints, deriv, restricted):
#   UNFINISHED/IMPOSSIBLE
#   # this needs custom parameters since B97 is re-parametrized again.
#   # But there is no interface in LibXC
#   # Call this first
#   sup = core.SuperFunctional.blank()
#   sup.set_max_points(npoints)
#   sup.set_deriv(deriv)
#
#   # => User-Customization <= #
#
#   # No spaces, keep it short and according to convention
#   sup.set_name('wB97X-D3')
#   # Tab in, trailing newlines
#   sup.set_description('    Parameterized Hybrid LRC B97 GGA XC Functional with D3(0) Dispersion\n')
#   # Tab in, trailing newlines
#   sup.set_citation(' Y.-S. Lin, G.-D. Li, S.-P. M., and J.-D. Chai J. Chem. Theory and Comput., 9,  263-272, 2013 \n')
#
#   # Add member functionals
#   wb97_x = core.LibXCFunctional('XC_GGA_X_', restricted)
#   wb97_x.set_omega(0.3)
#   wb97_x.set_alpha(0.804272)
#   sup.add_x_functional(wb97_x)
#   sup.add_c_functional(core.LibXCFunctional('XC_GGA_', restricted))
#
#   # Set GKS up after adding functionals
#   sup.set_x_omega(0.25)
#   sup.set_c_omega(0.0)
#   sup.set_x_alpha(0.195728)
#   sup.set_c_alpha(0.0)
#
#   # Call this last
#   sup.allocate()
#   return (sup, ('wB97', '-d3zero'))
