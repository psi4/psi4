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
from . import libxc_xc_funcs
from . import gga_superfuncs
from . import hyb_superfuncs
from . import mgga_superfuncs
from . import double_hyb_superfuncs

## ==> SuperFunctionals <== ##

superfunctionals = {}
superfunctionals.update(libxc_xc_funcs.libxc_xc_functional_list)
superfunctionals.update(gga_superfuncs.gga_superfunc_list)
superfunctionals.update(hyb_superfuncs.hyb_superfunc_list)
superfunctionals.update(mgga_superfuncs.mgga_superfunc_list)
superfunctionals.update(double_hyb_superfuncs.double_hyb_superfunc_list)

## ==> SuperFunctional List <== ##

superfunctional_list = []
for key in superfunctionals.keys():
    sup = superfunctionals[key](key, 1, 1, True)[0]
    superfunctional_list.append(sup)

## ==> Dispersion SuperFunctional List <== ##

p4_funcs = set([x for x in list(superfunctionals)])
p4_funcs -= set(['b97-d'])
for dashlvl, dashparam_dict in dftd3.dashcoeff.items():
    func_list = (set(dashparam_dict) & p4_funcs)
    for func in func_list:
        sup = superfunctionals[func](func, 1, 1, True)[0]
        sup.set_name(sup.name() + '-' + dashlvl.upper())
        superfunctional_list.append(sup)

        if dashlvl == 'd2p4':
            # -D2 overide
            sup = superfunctionals[func](func, 1, 1, True)[0]
            sup.set_name(sup.name() + '-D2')
            superfunctional_list.append(sup)

            # -D overide
            sup = superfunctionals[func](func, 1, 1, True)[0]
            sup.set_name(sup.name() + '-D')
            superfunctional_list.append(sup)

        if dashlvl == 'd3zero':
            sup = superfunctionals[func](func, 1, 1, True)[0]
            sup.set_name(sup.name() + '-D3')
            superfunctional_list.append(sup)

        if dashlvl == 'd3mzero':
            sup = superfunctionals[func](func, 1, 1, True)[0]
            sup.set_name(sup.name() + '-D3M')
            superfunctional_list.append(sup)

# # B97D is an odd one
for dashlvl in dftd3.full_dash_keys:
    if dashlvl == 'd2p4': continue

    sup = superfunctionals['b97-d']('b97-d', 1, 1, True)[0]
    sup.set_name('B97-' + dashlvl.upper())
    superfunctional_list.append(sup)

# wPBE, grr need a new scheme
for dashlvl in ['d3', 'd3m', 'd3zero', 'd3mzero', 'd3bj', 'd3mbj']:
    sup = superfunctionals['wpbe']('wpbe', 1, 1, True)[0]
    sup.set_name(sup.name() + '-' + dashlvl.upper())
    superfunctional_list.append(sup)


## ==> SuperFunctional Builder <== ##

def build_superfunctional(alias, restricted):
    name = alias.lower()

    npoints = core.get_option("SCF", "DFT_BLOCK_MAX_POINTS");
    deriv = 1 # Default depth for now

    # Grab out superfunctional
    if name in ["gen", ""]:
        sup = (core.get_option("DFT_CUSTOM_FUNCTIONAL"), False)
        if not isinstance(sup[0], core.SuperFunctional):
            raise KeyError("SCF: Custom Functional requested, but nothing provided in DFT_CUSTOM_FUNCTIONAL")

    elif name in superfunctionals.keys():
        sup = superfunctionals[name](name, npoints, deriv, restricted)

    elif name.upper() in superfunctionals.keys():
        sup = superfunctionals[name.upper()](name, npoints, deriv, restricted)



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

        func = superfunctionals[base_name](base_name, npoints, deriv, restricted)[0]

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
