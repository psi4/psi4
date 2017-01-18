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

"""Module with functions that call upon those in modules
:py:mod:`proc`, :py:mod:`driver`, and :py:mod:`wrappers`.

Place in this file quickly defined procedures such as
   - aliases for complex methods
   - simple modifications to existing methods

"""
from __future__ import print_function
from __future__ import absolute_import
import re
import os
import math
import warnings
from psi4.driver import driver_cbs
#from wrappers import *
#from gaussian_n import *  # CU
#from wrappers_cfour import *  # CU
#from qmmm import *  # CU

# Python procedures like these can be run directly from the input file or integrated
# with the energy(), etc. routines by means of lines like those at the end
# of this file.


def fake_file11(wfn, filename='fake_file11.dat', **kwargs):
    r"""Function to print a file *filename* of the old file11 format
    from molecule and gradient information in *wfn*.

    .. versionadded:: 0.6
       *wfn* parameter passed explicitly

    :returns: None

    :type filename: string
    :param filename: destination file name for file11 file

    :type wfn: :py:class:`~psi4.core.Wavefunction`
    :param wfn: set of molecule, gradient from which to generate file11

    :examples:

    >>> # [1] file11 for CISD calculation
    >>> G, wfn = gradient('cisd', return_wfn=True)
    >>> fake_file11(wfn, 'mycalc.11')

    """
    molecule = wfn.molecule()
    molecule.update_geometry()
    gradient = wfn.gradient()

    with open(filename, 'w') as handle:
        handle.write('%d\n' % (molecule.natom()))

        for at in range(molecule.natom()):
            handle.write('%6s %16.8f %16.8f %16.8f\n' % (molecule.symbol(
                at), molecule.x(at), molecule.y(at), molecule.z(at)))

        for at in range(molecule.natom()):
            handle.write('%6s %16.8f %16.8f %16.8f\n' % (
                '', gradient.get(at, 0), gradient.get(at, 1), gradient.get(at, 2)))


def sherrill_gold_standard(func, label, **kwargs):
    r"""Function to call the quantum chemical method known as 'Gold Standard'
    in the Sherrill group. Uses :py:func:`~driver_cbs.complete_basis_set` to evaluate
    the following expression. Two-point extrapolation of the correlation energy
    performed according to :py:func:`~driver_cbs.corl_xtpl_helgaker_2`.

    .. math:: E_{total}^{\text{Au\_std}} = E_{total,\; \text{SCF}}^{\text{aug-cc-pVQZ}} \; + E_{corl,\; \text{MP2}}^{\text{aug-cc-pV[TQ]Z}} \; + \delta_{\text{MP2}}^{\text{CCSD(T)}}\big\vert_{\text{aug-cc-pVTZ}}

    >>> # [1] single-point energy by this composite method
    >>> energy('sherrill_gold_standard')

    >>> # [2] finite-difference geometry optimization
    >>> optimize('sherrill_gold_standard')

    >>> # [3] finite-difference geometry optimization, overwriting some pre-defined sherrill_gold_standard options
    >>> optimize('sherrill_gold_standard', corl_basis='cc-pV[DT]Z', delta_basis='3-21g')

    """
    kwargs['scf_basis'] = kwargs.get('scf_basis', 'aug-cc-pVQZ')
    kwargs['scf_scheme'] = kwargs.get('scf_scheme', driver_cbs.xtpl_highest_1)

    kwargs['corl_wfn'] = kwargs.get('corl_wfn', 'mp2')
    kwargs['corl_basis'] = kwargs.get('corl_basis', 'aug-cc-pV[TQ]Z')
    kwargs['corl_scheme'] = kwargs.get('corl_scheme', driver_cbs.corl_xtpl_helgaker_2)

    kwargs['delta_wfn'] = kwargs.get('delta_wfn', 'ccsd(t)')
    kwargs['delta_wfn_lesser'] = kwargs.get('delta_wfn_lesser', 'mp2')
    kwargs['delta_basis'] = kwargs.get('delta_basis', 'aug-cc-pVTZ')
    kwargs['delta_scheme'] = kwargs.get('delta_scheme', driver_cbs.xtpl_highest_1)

    if label == 'custom_function':
        label = 'Sherrill Group Gold Standard'
    return driver_cbs.cbs(func, label, **kwargs)


def allen_focal_point(func, label, **kwargs):
    r"""Function to call Wes Allen-style Focal
    Point Analysis. JCP 127 014306.  Uses
    :py:func:`~driver_cbs.complete_basis_set` to evaluate the following
    expression. SCF employs a three-point extrapolation according
    to :py:func:`~driver_cbs.scf_xtpl_helgaker_3`. MP2, CCSD, and
    CCSD(T) employ two-point extrapolation performed according to
    :py:func:`~driver_cbs.corl_xtpl_helgaker_2`.  CCSDT and CCSDT(Q)
    are plain deltas. This wrapper requires :ref:`Kallay's MRCC code <sec:mrcc>`.

    .. math:: E_{total}^{\text{FPA}} = E_{total,\; \text{SCF}}^{\text{cc-pV[Q56]Z}} \; + E_{corl,\; \text{MP2}}^{\text{cc-pV[56]Z}} \; + \delta_{\text{MP2}}^{\text{CCSD}}\big\vert_{\text{cc-pV[56]Z}} \; + \delta_{\text{CCSD}}^{\text{CCSD(T)}}\big\vert_{\text{cc-pV[56]Z}} \; + \delta_{\text{CCSD(T)}}^{\text{CCSDT}}\big\vert_{\text{cc-pVTZ}} \; + \delta_{\text{CCSDT}}^{\text{CCSDT(Q)}}\big\vert_{\text{cc-pVDZ}}

    >>> # [1] single-point energy by this composite method
    >>> energy('allen_focal_point')

    >>> # [2] finite-difference geometry optimization embarrasingly parallel
    >>> optimize('allen_focal_point', mode='sow')

    """

    # SCF
    kwargs['scf_basis'] = kwargs.get('scf_basis', 'cc-pV[Q56]Z')
    kwargs['scf_scheme'] = kwargs.get('scf_scheme', driver_cbs.scf_xtpl_helgaker_3)

    # delta MP2 - SCF
    kwargs['corl_wfn'] = kwargs.get('corl_wfn', 'mp2')
    kwargs['corl_basis'] = kwargs.get('corl_basis', 'cc-pV[56]Z')
    kwargs['corl_scheme'] = kwargs.get('corl_scheme', driver_cbs.corl_xtpl_helgaker_2)

    # delta CCSD - MP2
    kwargs['delta_wfn'] = kwargs.get('delta_wfn', 'mrccsd')
    kwargs['delta_wfn_lesser'] = kwargs.get('delta_wfn_lesser', 'mp2')
    kwargs['delta_basis'] = kwargs.get('delta_basis', 'cc-pV[56]Z')
    kwargs['delta_scheme'] = kwargs.get('delta_scheme', driver_cbs.corl_xtpl_helgaker_2)

    # delta CCSD(T) - CCSD
    kwargs['delta2_wfn'] = kwargs.get('delta2_wfn', 'mrccsd(t)')
    kwargs['delta2_wfn_lesser'] = kwargs.get('delta2_wfn_lesser', 'mrccsd')
    kwargs['delta2_basis'] = kwargs.get('delta2_basis', 'cc-pV[56]Z')
    kwargs['delta2_scheme'] = kwargs.get('delta2_scheme', driver_cbs.corl_xtpl_helgaker_2)

    # delta CCSDT - CCSD(T)
    kwargs['delta3_wfn'] = kwargs.get('delta3_wfn', 'mrccsdt')
    kwargs['delta3_wfn_lesser'] = kwargs.get('delta3_wfn_lesser', 'mrccsd(t)')
    kwargs['delta3_basis'] = kwargs.get('delta3_basis', 'cc-pVTZ')
    kwargs['delta3_scheme'] = kwargs.get('delta3_scheme', driver_cbs.xtpl_highest_1)

    # delta CCSDT(Q) - CCSDT
    kwargs['delta4_wfn'] = kwargs.get('delta4_wfn', 'mrccsdt(q)')
    kwargs['delta4_wfn_lesser'] = kwargs.get('delta4_wfn_lesser', 'mrccsdt')
    kwargs['delta4_basis'] = kwargs.get('delta4_basis', 'cc-pVDZ')
    kwargs['delta4_scheme'] = kwargs.get('delta4_scheme', driver_cbs.xtpl_highest_1)

    if label == 'custom_function':
        label = 'Allen Focal Point'
    return driver_cbs.cbs(func, label, **kwargs)

