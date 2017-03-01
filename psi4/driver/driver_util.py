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

from __future__ import print_function
import math
import re
from psi4 import core
from psi4.driver import qcdb
from psi4.driver import p4util
from psi4.driver.p4util.exceptions import *
from psi4.driver.procrouting import *


def _method_exists(ptype, method_name):
    r"""
    Quick check to see if this method exists, if it does not exist we raise a convenient flag.
    """
    if method_name not in procedures[ptype].keys():
        alternatives = ""
        alt_method_name = p4util.text.find_approximate_string_matches(method_name,
                                                                procedures[ptype].keys(), 2)
        if len(alt_method_name) > 0:
            alternatives = " Did you mean? %s" % (" ".join(alt_method_name))
        Cptype = ptype[0].upper() + ptype[1:]
        raise ValidationError('%s method "%s" is not available.%s' % (Cptype, method_name, alternatives))


def _set_convergence_criterion(ptype, method_name, scf_Ec, pscf_Ec, scf_Dc, pscf_Dc, gen_Ec, verbose=1):
    r"""
    This function will set local SCF and global energy convergence criterion
    to the defaults listed at:
    http://www.psicode.org/psi4manual/master/scf.html#convergence-and-
    algorithm-defaults. SCF will be converged more tightly if a post-SCF
    method is select (pscf_Ec, and pscf_Dc) else the looser (scf_Ec, and
    scf_Dc convergence criterion will be used).

    ptype -         Procedure type (energy, gradient, etc). Nearly always test on
                    procedures['energy'] since that's guaranteed to exist for a method.
    method_name -   Name of the method
    scf_Ec -        E convergence criterion for scf target method
    pscf_Ec -       E convergence criterion for scf of post-scf target method
    scf_Dc -        D convergence criterion for scf target method
    pscf_Dc -       D convergence criterion for scf of post-scf target method
    gen_Ec -        E convergence criterion for post-scf target method

    """
    optstash = p4util.OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['E_CONVERGENCE'])

    # Kind of want to move this out of here
    _method_exists(ptype, method_name)

    if verbose >= 2:
        print('      Setting convergence', end=' ')
    # Set method-dependent scf convergence criteria, check against energy routines
    if not core.has_option_changed('SCF', 'E_CONVERGENCE'):
        if procedures['energy'][method_name] in [proc.run_scf, proc.run_dft]:
            core.set_local_option('SCF', 'E_CONVERGENCE', scf_Ec)
            if verbose >= 2:
                print(scf_Ec, end=' ')
        else:
            core.set_local_option('SCF', 'E_CONVERGENCE', pscf_Ec)
            if verbose >= 2:
                print(pscf_Ec, end=' ')
    else:
        if verbose >= 2:
            print('CUSTOM', core.get_option('SCF', 'E_CONVERGENCE'), end=' ')

    if not core.has_option_changed('SCF', 'D_CONVERGENCE'):
        if procedures['energy'][method_name] in [proc.run_scf, proc.run_dft]:
            core.set_local_option('SCF', 'D_CONVERGENCE', scf_Dc)
            if verbose >= 2:
                print(scf_Dc, end=' ')
        else:
            core.set_local_option('SCF', 'D_CONVERGENCE', pscf_Dc)
            if verbose >= 2:
                print(pscf_Dc, end=' ')
    else:
        if verbose >= 2:
            print('CUSTOM', core.get_option('SCF', 'D_CONVERGENCE'), end=' ')

    # Set post-scf convergence criteria (global will cover all correlated modules)
    if not core.has_global_option_changed('E_CONVERGENCE'):
        if procedures['energy'][method_name] not in [proc.run_scf, proc.run_dft]:
            core.set_global_option('E_CONVERGENCE', gen_Ec)
            if verbose >= 2:
                print(gen_Ec, end=' ')
    else:
        if procedures['energy'][method_name] not in [proc.run_scf, proc.run_dft]:
            if verbose >= 2:
                print('CUSTOM', core.get_global_option('E_CONVERGENCE'), end=' ')

    if verbose >= 2:
        print('')
    return optstash


def parse_arbitrary_order(name):
    r"""Function to parse name string into a method family like CI or MRCC and specific
    level information like 4 for CISDTQ or MRCCSDTQ.

    """

    name = name.lower()

    # matches 'mrccsdt(q)'
    if name.startswith('mrcc'):

        # avoid undoing fn's good work when called twice
        if name == 'mrcc':
            return name, None

        # grabs 'sdt(q)'
        ccfullname = name[4:]

        # A negative order indicates perturbative method
        methods = {
            'sd'          : { 'method': 1, 'order':  2, 'fullname': 'CCSD'         },
            'sdt'         : { 'method': 1, 'order':  3, 'fullname': 'CCSDT'        },
            'sdtq'        : { 'method': 1, 'order':  4, 'fullname': 'CCSDTQ'       },
            'sdtqp'       : { 'method': 1, 'order':  5, 'fullname': 'CCSDTQP'      },
            'sdtqph'      : { 'method': 1, 'order':  6, 'fullname': 'CCSDTQPH'     },
            'sd(t)'       : { 'method': 3, 'order': -3, 'fullname': 'CCSD(T)'      },
            'sdt(q)'      : { 'method': 3, 'order': -4, 'fullname': 'CCSDT(Q)'     },
            'sdtq(p)'     : { 'method': 3, 'order': -5, 'fullname': 'CCSDTQ(P)'    },
            'sdtqp(h)'    : { 'method': 3, 'order': -6, 'fullname': 'CCSDTQP(H)'   },
            'sd(t)_l'     : { 'method': 4, 'order': -3, 'fullname': 'CCSD(T)_L'    },
            'sdt(q)_l'    : { 'method': 4, 'order': -4, 'fullname': 'CCSDT(Q)_L'   },
            'sdtq(p)_l'   : { 'method': 4, 'order': -5, 'fullname': 'CCSDTQ(P)_L'  },
            'sdtqp(h)_l'  : { 'method': 4, 'order': -6, 'fullname': 'CCSDTQP(H)_L' },
            'sdt-1a'      : { 'method': 5, 'order':  3, 'fullname': 'CCSDT-1a'     },
            'sdtq-1a'     : { 'method': 5, 'order':  4, 'fullname': 'CCSDTQ-1a'    },
            'sdtqp-1a'    : { 'method': 5, 'order':  5, 'fullname': 'CCSDTQP-1a'   },
            'sdtqph-1a'   : { 'method': 5, 'order':  6, 'fullname': 'CCSDTQPH-1a'  },
            'sdt-1b'      : { 'method': 6, 'order':  3, 'fullname': 'CCSDT-1b'     },
            'sdtq-1b'     : { 'method': 6, 'order':  4, 'fullname': 'CCSDTQ-1b'    },
            'sdtqp-1b'    : { 'method': 6, 'order':  5, 'fullname': 'CCSDTQP-1b'   },
            'sdtqph-1b'   : { 'method': 6, 'order':  6, 'fullname': 'CCSDTQPH-1b'  },
            '2'           : { 'method': 7, 'order':  2, 'fullname': 'CC2'          },
            '3'           : { 'method': 7, 'order':  3, 'fullname': 'CC3'          },
            '4'           : { 'method': 7, 'order':  4, 'fullname': 'CC4'          },
            '5'           : { 'method': 7, 'order':  5, 'fullname': 'CC5'          },
            '6'           : { 'method': 7, 'order':  6, 'fullname': 'CC6'          },
            'sdt-3'       : { 'method': 8, 'order':  3, 'fullname': 'CCSDT-3'      },
            'sdtq-3'      : { 'method': 8, 'order':  4, 'fullname': 'CCSDTQ-3'     },
            'sdtqp-3'     : { 'method': 8, 'order':  5, 'fullname': 'CCSDTQP-3'    },
            'sdtqph-3'    : { 'method': 8, 'order':  6, 'fullname': 'CCSDTQPH-3'   }
        }

        # looks for 'sdt(q)' in dictionary
        if ccfullname in methods:
            return 'mrcc', methods[ccfullname]
        else:
            raise ValidationError('MRCC method \'%s\' invalid.' % (name))

    elif re.match(r'^[a-z]+\d+$', name):
        decompose = re.compile(r'^([a-z]+)(\d+)$').match(name)
        namestump = decompose.group(1)
        namelevel = int(decompose.group(2))

        if namestump in ['mp', 'zapt', 'ci']:
            # Let mp2, mp3, mp4 pass through to select functions
            if namestump == 'mp' and namelevel in [2, 3, 4]:
                return name, None
            # Otherwise return method and order
            else:
                return namestump, namelevel
        else:
            return name, None
    else:
        return name, None


def parse_cotton_irreps(irrep, point_group):
    r"""Function to return validated Cotton ordering index for molecular
    *point_group* from string or integer irreducible representation *irrep*.

    """
    cotton = {
        'c1': {
            'a': 1,
            '1': 1
        },
        'ci': {
            'ag': 1,
            'au': 2,
            '1': 1,
            '2': 2
        },
        'c2': {
            'a': 1,
            'b': 2,
            '1': 1,
            '2': 2
        },
        'cs': {
            'ap': 1,
            'app': 2,
            '1': 1,
            '2': 2
        },
        'd2': {
            'a': 1,
            'b1': 2,
            'b2': 3,
            'b3': 4,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4
        },
        'c2v': {
            'a1': 1,
            'a2': 2,
            'b1': 3,
            'b2': 4,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4
        },
        'c2h': {
            'ag': 1,
            'bg': 2,
            'au': 3,
            'bu': 4,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4,
        },
        'd2h': {
            'ag': 1,
            'b1g': 2,
            'b2g': 3,
            'b3g': 4,
            'au': 5,
            'b1u': 6,
            'b2u': 7,
            'b3u': 8,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4,
            '5': 5,
            '6': 6,
            '7': 7,
            '8': 8
        }
    }

    try:
        return cotton[point_group.lower()][str(irrep).lower()]
    except KeyError:
        raise ValidationError("""Irrep '%s' not valid for point group '%s'.""" % (str(irrep), point_group))
