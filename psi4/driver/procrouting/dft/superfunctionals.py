#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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
Module to provide lightweight definitions of functionals and
SuperFunctionals
"""
import re
import os

from psi4 import core
from psi4.driver.p4util.exceptions import ValidationError
from . import dft_builder


def build_superfunctional(name, restricted, npoints=None, deriv=1):
    if npoints is None:
        npoints = core.get_option("SCF", "DFT_BLOCK_MAX_POINTS")

    # We are a XC generating function

    if hasattr(name, '__call__'):
        custom_error = "SCF: Custom functional type must either be a SuperFunctional or a tuple of (SuperFunctional, (base_name, dashparam))."
        sfunc = name("name", npoints, deriv, restricted)

        # Without Dispersion
        if isinstance(sfunc, core.SuperFunctional):
            sup = (sfunc, False)
        # With Dispersion
        elif isinstance(sfunc[0], core.SuperFunctional):
            sup = sfunc
            # Can we validate dispersion?
        else:
            raise ValidationError(custom_error)

        # Double check that the SuperFunctional is correctly sized (why dont we always do this?)
        sup[0].set_max_points(npoints)
        sup[0].set_deriv(deriv)
        sup[0].allocate()

    # Check for supplied dict_func functionals
    elif isinstance(name, dict):
        sup = dft_builder.build_superfunctional_from_dictionary(name, npoints, deriv, restricted)
    # Check for pre-defined dict-based functionals
    elif name.lower() in dft_builder.functionals:
        sup = dft_builder.build_superfunctional_from_dictionary(dft_builder.functionals[name.lower()], npoints, deriv,
                                                                restricted)
    else:
        raise ValidationError("SCF: Functional (%s) not found!" % name)

    if (core.get_global_option('INTEGRAL_PACKAGE') == 'ERD') and (sup[0].is_x_lrc() or sup[0].is_c_lrc()):
        raise ValidationError("INTEGRAL_PACKAGE ERD does not play nicely with omega ERI's, so stopping.")

    # Lock and unlock the functional
    sup[0].set_lock(False)

    # set LibXC density screening
    dens_tol = core.get_option("SCF", "DFT_DENSITY_TOLERANCE")
    if (dens_tol > 0.0):
        sup[0].set_density_tolerance(dens_tol)

    # Set options
    if core.has_option_changed("SCF", "DFT_OMEGA") and sup[0].is_x_lrc():
        omega = core.get_option("SCF", "DFT_OMEGA")
        sup[0].set_x_omega(omega)

        # We also need to loop through all of the exchange functionals
        if sup[0].is_libxc_func():
            # Full libxc funcs are dropped in c_functionals (smooth move!)
            sup[0].c_functionals()[0].set_omega(omega)
        else:
            for x_func in sup[0].x_functionals():
                x_func.set_omega(omega)
    if core.has_option_changed("SCF", "DFT_OMEGA_C") and sup[0].is_c_lrc():
        sup[0].set_c_omega(core.get_option("SCF", "DFT_OMEGA_C"))

    if core.has_option_changed("SCF", "DFT_ALPHA"):
        sup[0].set_x_alpha(core.get_option("SCF", "DFT_ALPHA"))
    if core.has_option_changed("SCF", "DFT_ALPHA_C"):
        sup[0].set_c_alpha(core.get_option("SCF", "DFT_ALPHA_C"))

    # add VV10 correlation to any functional or modify existing
    # custom procedures using name 'scf' without any quadrature grid like HF will fail and are not detected
    if (core.has_option_changed("SCF", "NL_DISPERSION_PARAMETERS") and sup[0].vv10_b() > 0.0):
        if not isinstance(name, dict):
            if (name.lower() == 'hf'):
                raise ValidationError("SCF: HF with -NL not implemented")
        nl_tuple = core.get_option("SCF", "NL_DISPERSION_PARAMETERS")
        sup[0].set_vv10_b(nl_tuple[0])
        if len(nl_tuple) > 1:
            sup[0].set_vv10_c(nl_tuple[1])
        if len(nl_tuple) > 2:
            raise ValidationError("too many entries in NL_DISPERSION_PARAMETERS for DFT-NL")
    elif core.has_option_changed("SCF", "DFT_VV10_B"):
        if not isinstance(name, dict):
            if (name.lower() == 'hf'):
                raise ValidationError("SCF: HF with -NL not implemented")
        vv10_b = core.get_option("SCF", "DFT_VV10_B")
        sup[0].set_vv10_b(vv10_b)
        if core.has_option_changed("SCF", "DFT_VV10_C"):
            vv10_c = core.get_option("SCF", "DFT_VV10_C")
            sup[0].set_vv10_c(vv10_c)
        if (abs(sup[0].vv10_c() - 0.0) <= 1e-8):
            core.print_out("SCF: VV10_C not specified. Using default (C=0.0093)!")
            sup[0].set_vv10_c(0.0093)

    if (core.has_option_changed("SCF", "NL_DISPERSION_PARAMETERS") and core.has_option_changed("SCF", "DFT_VV10_B")):
        raise ValidationError("SCF: Decide between NL_DISPERSION_PARAMETERS and DFT_VV10_B !!")

    # Check SCF_TYPE
    if sup[0].is_x_lrc() and (core.get_global_option("SCF_TYPE")
                              not in ["DISK_DF", "MEM_DF", "DIRECT", "DF", "OUT_OF_CORE", "PK"]):
        raise ValidationError(
            "SCF: SCF_TYPE (%s) not supported for range-separated functionals, plese use SCF_TYPE = 'DF' to automatically select the correct JK build."
            % core.get_global_option("SCF_TYPE"))

    if (core.get_global_option('INTEGRAL_PACKAGE') == 'ERD') and (sup[0].is_x_lrc()):
        raise ValidationError('INTEGRAL_PACKAGE ERD does not play nicely with LRC DFT functionals, so stopping.')

    sup[0].set_lock(True)

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

    rho_line = re.compile(
        r'^\s*rhoa=\s*(-?\d+\.\d+E[+-]\d+)\s*rhob=\s*(-?\d+\.\d+E[+-]\d+)\s*sigmaaa=\s*(-?\d+\.\d+E[+-]\d+)\s*sigmaab=\s*(-?\d+\.\d+E[+-]\d+)\s*sigmabb=\s*(-?\d+\.\d+E[+-]\d+)\s*'
    )
    val_line = re.compile(r'^\s*(\w*)\s*=\s*(-?\d+\.\d+E[+-]\d+)')

    aliases = {
        'zk': 'v',
        'vrhoa': 'v_rho_a',
        'vrhob': 'v_rho_b',
        'vsigmaaa': 'v_gamma_aa',
        'vsigmaab': 'v_gamma_ab',
        'vsigmabb': 'v_gamma_bb',
        'v2rhoa2': 'v_rho_a_rho_a',
        'v2rhoab': 'v_rho_a_rho_b',
        'v2rhob2': 'v_rho_b_rho_b',
        'v2rhoasigmaaa': 'v_rho_a_gamma_aa',
        'v2rhoasigmaab': 'v_rho_a_gamma_ab',
        'v2rhoasigmabb': 'v_rho_a_gamma_bb',
        'v2rhobsigmaaa': 'v_rho_b_gamma_aa',
        'v2rhobsigmaab': 'v_rho_b_gamma_ab',
        'v2rhobsigmabb': 'v_rho_b_gamma_bb',
        'v2sigmaaa2': 'v_gamma_aa_gamma_aa',
        'v2sigmaaaab': 'v_gamma_aa_gamma_ab',
        'v2sigmaaabb': 'v_gamma_aa_gamma_bb',
        'v2sigmaab2': 'v_gamma_ab_gamma_ab',
        'v2sigmaabbb': 'v_gamma_ab_gamma_bb',
        'v2sigmabb2': 'v_gamma_bb_gamma_bb',
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

    super = build_superfunctional(functional, True, npoints=N, deriv=1)
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
        core.print_out('rho_a= %11.3E, rho_b= %11.3E, gamma_aa= %11.3E, gamma_ab= %11.3E, gamma_bb= %11.3E\n' %
                       (rho_a[index], rho_b[index], gamma_aa[index], gamma_ab[index], gamma_bb[index]))

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

            core.print_out('\t%-15s %24.16E %24.16E %24.16E %24.16E %6s\n' %
                           (task, v_ref, v_obs, delta, epsilon, passed))

        index = index + 1

    core.print_out('\n')
    return check
