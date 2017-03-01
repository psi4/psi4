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
from __future__ import absolute_import
from __future__ import division
import re
import math
import sys
import numpy as np

from psi4 import core

from psi4.driver import qcdb
from psi4.driver import p4util
from psi4.driver import driver_util
from psi4.driver import constants

from psi4.driver.p4util.exceptions import *
from psi4.driver.procrouting.interface_cfour import cfour_psivar_list

zeta_values = ['d', 't', 'q', '5', '6', '7', '8']
zeta_val2sym = {k + 2: v for k, v in zip(range(7), zeta_values)}
zeta_sym2val = {v: k for k, v in zeta_val2sym.items()}


def _expand_bracketed_basis(basisstring, molecule=None):
    r"""Function to transform and validate basis series specification
    *basisstring* for cbs(). A basis set with no paired square brackets is
    passed through with zeta level 0 (e.g., '6-31+G(d,p)' is returned as
    [6-31+G(d,p)] and [0]). A basis set with square brackets is checked
    for sensible sequence and Dunning-ness and returned as separate basis
    sets (e.g., 'cc-pV[Q5]Z' is returned as [cc-pVQZ, cc-pV5Z] and [4,
    5]). This function checks that the basis is valid by trying to build
    the qcdb.BasisSet object for *molecule* or for H2 if None. Allows
    out-of-order zeta specification (e.g., [qtd]) and numeral for number
    (e.g., [23]) but not skipped zetas (e.g., [dq]) or zetas outside [2,
    8] or non-Dunning sets or non-findable .gbs sets.

    """
    BSET = []
    ZSET = []
    legit_compound_basis = re.compile(r'^(?P<pre>.*cc-.*)\[(?P<zeta>[dtq2345678,]*)\](?P<post>.*z)$', re.IGNORECASE)

    if legit_compound_basis.match(basisstring):
        basisname = legit_compound_basis.match(basisstring)
        # filter out commas and be forgiving of e.g., t5q or 3q
        bn_gz = basisname.group('zeta')
        zetas = [z for z in zeta_values if (z in bn_gz or str(zeta_values.index(z) + 2) in bn_gz)]
        for b in zetas:
            if ZSET and (int(ZSET[len(ZSET) - 1]) - zeta_values.index(b)) != 1:
                    raise ValidationError("""Basis set '%s' has skipped zeta level '%s'.""" % (basisstring, b))
            BSET.append(basisname.group('pre') + b + basisname.group('post'))
            ZSET.append(zeta_values.index(b) + 2)
    elif re.match(r'.*\[.*\].*$', basisstring, flags=re.IGNORECASE):
        raise ValidationError("""Basis series '%s' invalid. Specify a basis series matching"""
                              """ '*cc-*[dtq2345678,]*z'.""" % (basisstring))
    else:
        BSET.append(basisstring)
        ZSET.append(0)

    if molecule is None:
        molecule = """\nH\nH 1 1.00\n"""

    for basis in BSET:
        try:
            qcdb.BasisSet.pyconstruct(molecule, "BASIS", basis)
        except qcdb.BasisSetNotFound:
            e=sys.exc_info()[1]
            raise ValidationError("""Basis set '%s' not available for molecule.""" % (basis))

    return (BSET, ZSET)


def _contract_bracketed_basis(basisarray):
    r"""Function to reform a bracketed basis set string from a sequential series
    of basis sets *basisarray* (e.g, form 'cc-pv[q5]z' from array [cc-pvqz, cc-pv5z]).
    Used to print a nicely formatted basis set string in the results table.

    """
    if len(basisarray) == 1:
        return basisarray[0]

    else:
        zetaindx = [i for i in range(len(basisarray[0])) if basisarray[0][i] != basisarray[1][i]][0]
        ZSET = [bas[zetaindx] for bas in basisarray]

        pre = basisarray[0][:zetaindx]
        post = basisarray[0][zetaindx + 1:]
        basisstring = pre + '[' + ''.join(ZSET) + ']' + post
        return basisstring


def xtpl_highest_1(functionname, zHI, valueHI, verbose=True):
    r"""Scheme for total or correlation energies with a single basis or the highest
    zeta-level among an array of bases. Used by :py:func:`~psi4.cbs`.

    .. math:: E_{total}^X = E_{total}^X

    """
    if isinstance(valueHI, float):

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> %s <==\n\n""" % (functionname.upper())
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)

            core.print_out(cbsscheme)

        return valueHI

    elif isinstance(valueHI, (core.Matrix, core.Vector)):

        if verbose > 2:
            core.print_out("""   HI-zeta (%s) Total Energy:\n""" % (str(zHI)))
            valueHI.print_out()

        return valueHI


def scf_xtpl_helgaker_2(functionname, zLO, valueLO, zHI, valueHI, verbose=True, alpha=1.63):
    r"""Extrapolation scheme for reference energies with two adjacent zeta-level bases.
    Used by :py:func:`~psi4.cbs`.

    .. math:: E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha X}, \alpha = 1.63

    """

    if type(valueLO) != type(valueHI):
        raise ValidationError("scf_xtpl_helgaker_2: Inputs must be of the same datatype! (%s, %s)"
                              % (type(valueLO), type(valueHI)))

    beta_division = 1 / (math.exp(-1 * alpha * zLO) * (math.exp(-1 * alpha) - 1))
    beta_mult = math.exp(-1 * alpha * zHI)

    if isinstance(valueLO, float):
        beta = (valueHI - valueLO) / (math.exp(-1 * alpha * zLO) * (math.exp(-1 * alpha) - 1))
        value = valueHI - beta * math.exp(-1 * alpha * zHI)

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Helgaker 2-point SCF extrapolation for method: %s <==\n\n""" % (functionname.upper())
            cbsscheme += """   LO-zeta (%s) Energy:               % 16.12f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Alpha (exponent) Value:           % 16.12f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n\n""" % (beta)

            name_str = "%s/(%s,%s)" % (functionname.upper(), zeta_val2sym[zLO].upper(), zeta_val2sym[zHI].upper())
            cbsscheme += """   @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (18 - len(name_str))
            cbsscheme += """% 16.12f\n\n""" % value
            core.print_out(cbsscheme)

        return value

    elif isinstance(valueLO, (core.Matrix, core.Vector)):
        beta = valueHI.clone()
        beta.name = 'Helgaker SCF (%s, %s) beta' % (zLO, zHI)
        beta.subtract(valueLO)
        beta.scale(beta_division)
        beta.scale(beta_mult)

        value = valueHI.clone()
        value.subtract(beta)
        value.name = 'Helgaker SCF (%s, %s) data' % (zLO, zHI)

        if verbose > 2:
            core.print_out("""\n   ==> Helgaker 2-point SCF extrapolation for method: %s <==\n\n""" % (functionname.upper()))
            core.print_out("""   LO-zeta (%s)""" % str(zLO))
            core.print_out("""   LO-zeta Data""")
            valueLO.print_out()
            core.print_out("""   HI-zeta (%s)""" % str(zHI))
            core.print_out("""   HI-zeta Data""")
            valueHI.print_out()
            core.print_out("""   Extrapolated Data:\n""")
            value.print_out()
            core.print_out("""   Alpha (exponent) Value:          %16.8f\n""" % (alpha))
            core.print_out("""   Beta Data:\n""")
            beta.print_out()

        return value

    else:
        raise ValidationError("scf_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))


def scf_xtpl_helgaker_3(functionname, zLO, valueLO, zMD, valueMD, zHI, valueHI, verbose=True):
    r"""Extrapolation scheme for reference energies with three adjacent zeta-level bases.
    Used by :py:func:`~psi4.cbs`.

    .. math:: E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha X}
    """

    if (type(valueLO) != type(valueMD)) or (type(valueMD) != type(valueHI)):
        raise ValidationError("scf_xtpl_helgaker_3: Inputs must be of the same datatype! (%s, %s, %s)"
                              % (type(valueLO), type(valueMD), type(valueHI)))

    if isinstance(valueLO, float):

        ratio = (valueHI - valueMD) / (valueMD - valueLO)
        alpha = -1 * math.log(ratio)
        beta = (valueHI - valueMD) / (math.exp(-1 * alpha * zMD) * (ratio - 1))
        value = valueHI - beta * math.exp(-1 * alpha * zHI)

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Helgaker 3-point SCF extrapolation for method: %s <==\n\n""" % (functionname.upper())
            cbsscheme += """   LO-zeta (%s) Energy:               % 16.12f\n""" % (str(zLO), valueLO)
            cbsscheme += """   MD-zeta (%s) Energy:               % 16.12f\n""" % (str(zMD), valueMD)
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Alpha (exponent) Value:           % 16.12f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n\n""" % (beta)

            name_str = "%s/(%s,%s,%s)" % (functionname.upper(), zeta_val2sym[zLO].upper(), zeta_val2sym[zMD].upper(),
                                                             zeta_val2sym[zHI].upper())
            cbsscheme += """   @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (18 - len(name_str))
            cbsscheme += """% 16.12f\n\n""" % value
            core.print_out(cbsscheme)

        return value

    elif isinstance(valueLO, (core.Matrix, core.Vector)):
        valueLO = np.array(valueLO)
        valueMD = np.array(valueMD)
        valueHI = np.array(valueHI)

        nonzero_mask = np.abs(valueHI) > 1.e-14
        top = (valueHI - valueMD)[nonzero_mask]
        bot = (valueMD - valueLO)[nonzero_mask]

        ratio = top / bot
        alpha = -1 * np.log(np.abs(ratio))
        beta = top / (np.exp(-1 * alpha * zMD) * (ratio - 1))
        np_value = valueHI.copy()
        np_value[nonzero_mask] -= beta * np.exp(-1 * alpha * zHI)
        np_value[~nonzero_mask] = 0.0

        # Build and set from numpy routines
        value = core.Matrix(*valueHI.shape)
        value_view = np.asarray(value)
        value_view[:] = np_value
        return value

    else:
        raise ValidationError("scf_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))


#def corl_xtpl_helgaker_2(functionname, valueSCF, zLO, valueLO, zHI, valueHI, verbose=True):
def corl_xtpl_helgaker_2(functionname, zLO, valueLO, zHI, valueHI, verbose=True):
    r"""Extrapolation scheme for correlation energies with two adjacent zeta-level bases.
    Used by :py:func:`~psi4.cbs`.

    .. math:: E_{corl}^X = E_{corl}^{\infty} + \beta X^{-3}

    """
    if type(valueLO) != type(valueHI):
        raise ValidationError("corl_xtpl_helgaker_2: Inputs must be of the same datatype! (%s, %s)"
                              % (type(valueLO), type(valueHI)))

    if isinstance(valueLO, float):
        value = (valueHI * zHI ** 3 - valueLO * zLO ** 3) / (zHI ** 3 - zLO ** 3)
        beta = (valueHI - valueLO) / (zHI ** (-3) - zLO ** (-3))

#        final = valueSCF + value
        final = value
        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = """\n\n   ==> Helgaker 2-point correlated extrapolation for method: %s <==\n\n""" % (functionname.upper())
#            cbsscheme += """   HI-zeta (%1s) SCF Energy:           % 16.12f\n""" % (str(zHI), valueSCF)
            cbsscheme += """   LO-zeta (%s) Energy:               % 16.12f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)
#            cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n""" % beta
            cbsscheme += """   Extrapolated Energy:              % 16.12f\n\n""" % value
            #cbsscheme += """   LO-zeta (%s) Correlation Energy:   % 16.12f\n""" % (str(zLO), valueLO)
            #cbsscheme += """   HI-zeta (%s) Correlation Energy:   % 16.12f\n""" % (str(zHI), valueHI)
            #cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n""" % beta
            #cbsscheme += """   Extrapolated Correlation Energy:  % 16.12f\n\n""" % value

            name_str = "%s/(%s,%s)" % (functionname.upper(), zeta_val2sym[zLO].upper(), zeta_val2sym[zHI].upper())
            cbsscheme += """   @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (19 - len(name_str))
            cbsscheme += """% 16.12f\n\n""" % final
            core.print_out(cbsscheme)

        return final

    elif isinstance(valueLO, (core.Matrix, core.Vector)):

        beta = valueHI.clone()
        beta.subtract(valueLO)
        beta.scale(1 / (zHI ** (-3) - zLO ** (-3)))
        beta.name = 'Helgaker SCF (%s, %s) beta' % (zLO, zHI)

        value = valueHI.clone()
        value.scale(zHI ** 3)

        tmp = valueLO.clone()
        tmp.scale(zLO ** 3)
        value.subtract(tmp)

        value.scale(1 / (zHI ** 3 - zLO ** 3))
        value.name = 'Helgaker Corr (%s, %s) data' % (zLO, zHI)

        if verbose > 2:
            core.print_out("""\n   ==> Helgaker 2-point correlated extrapolation for """
                           """method: %s <==\n\n""" % (functionname.upper()))
            core.print_out("""   LO-zeta (%s) Data\n""" % (str(zLO)))
            valueLO.print_out()
            core.print_out("""   HI-zeta (%s) Data\n""" % (str(zHI)))
            valueHI.print_out()
            core.print_out("""   Extrapolated Data:\n""")
            value.print_out()
            core.print_out("""   Beta Data:\n""")
            beta.print_out()

#        value.add(valueSCF)
        return value

    else:
        raise ValidationError("scf_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))


def return_energy_components():
    VARH = {}
    VARH['scf'] = {
                            'scf': 'SCF TOTAL ENERGY'}
    VARH['hf'] = {
                             'hf': 'HF TOTAL ENERGY'}
    VARH['mp2'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY'}
    VARH['mp2.5'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY'}
    VARH['mp3'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY'}
    VARH['mp4(sdq)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY',
                       'mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY'}
    VARH['mp4'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY',
                       'mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY',
                            'mp4': 'MP4(SDTQ) TOTAL ENERGY'}
    VARH['omp2'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'omp2': 'OMP2 TOTAL ENERGY'}
    VARH['omp2.5'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                         'omp2.5': 'OMP2.5 TOTAL ENERGY'}
    VARH['omp3'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY',
                           'omp3': 'OMP3 TOTAL ENERGY'}
    VARH['olccd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'olccd': 'OLCCD TOTAL ENERGY'}
    VARH['lccd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'lccd': 'LCCD TOTAL ENERGY'}
    VARH['lccsd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'lccsd': 'LCCSD TOTAL ENERGY'}
    VARH['cepa(0)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                        'cepa(0)': 'CEPA(0) TOTAL ENERGY'}
    VARH['cepa(1)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                        'cepa(1)': 'CEPA(1) TOTAL ENERGY'}
    VARH['cepa(3)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                        'cepa(3)': 'CEPA(3) TOTAL ENERGY'}
    VARH['acpf'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'acpf': 'ACPF TOTAL ENERGY'}
    VARH['aqcc'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'aqcc': 'AQCC TOTAL ENERGY'}
    VARH['qcisd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY',
                       'mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY',
                          'qcisd': 'QCISD TOTAL ENERGY'}
    VARH['cc2'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                            'cc2': 'CC2 TOTAL ENERGY'}
    VARH['ccsd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'ccsd': 'CCSD TOTAL ENERGY'}
    VARH['bccd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'bccd': 'CCSD TOTAL ENERGY'}
    VARH['cc3'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                            'cc3': 'CC3 TOTAL ENERGY'}
    VARH['fno-ccsd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                       'fno-ccsd': 'CCSD TOTAL ENERGY'}
    VARH['fno-ccsd(t)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'ccsd': 'CCSD TOTAL ENERGY',
                    'fno-ccsd(t)': 'CCSD(T) TOTAL ENERGY'}
    VARH['qcisd(t)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY',
                       'mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY',
                          'qcisd': 'QCISD TOTAL ENERGY',
                       'qcisd(t)': 'QCISD(T) TOTAL ENERGY'}
    VARH['ccsd(t)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'ccsd': 'CCSD TOTAL ENERGY',
                        'ccsd(t)': 'CCSD(T) TOTAL ENERGY'}
    VARH['bccd(t)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'ccsd': 'CCSD TOTAL ENERGY',
                        'bccd(t)': 'CCSD(T) TOTAL ENERGY'}
    VARH['cisd'] = {
                             'hf': 'HF TOTAL ENERGY',
                           'cisd': 'CISD TOTAL ENERGY'}
    VARH['cisdt'] = {
                             'hf': 'HF TOTAL ENERGY',
                          'cisdt': 'CISDT TOTAL ENERGY'}
    VARH['cisdtq'] = {
                             'hf': 'HF TOTAL ENERGY',
                         'cisdtq': 'CISDTQ TOTAL ENERGY'}
    VARH['fci'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'fci': 'FCI TOTAL ENERGY'}
    VARH['mrccsd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                         'mrccsd': 'CCSD TOTAL ENERGY'}
    VARH['mrccsd(t)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                         'mrccsd': 'CCSD TOTAL ENERGY',
                      'mrccsd(t)': 'CCSD(T) TOTAL ENERGY'}
    VARH['mrccsdt'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                        'mrccsdt': 'CCSDT TOTAL ENERGY'}
    VARH['mrccsdt(q)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                        'mrccsdt': 'CCSDT TOTAL ENERGY',
                     'mrccsdt(q)': 'CCSDT(Q) TOTAL ENERGY'}

    for cilevel in range(2, 99):
        VARH['ci%s' % (str(cilevel))] = {
                             'hf': 'HF TOTAL ENERGY',
          'ci%s' % (str(cilevel)): 'CI TOTAL ENERGY'}

    for mplevel in range(5, 99):
        VARH['mp%s' % (str(mplevel))] = {
                             'hf': 'HF TOTAL ENERGY',
          'mp%s' % (str(mplevel)): 'MP%s TOTAL ENERGY' % (str(mplevel))}
        for mplevel2 in range(2, mplevel):
            VARH['mp%s' % (str(mplevel))]['mp%s' % (str(mplevel2))] = \
                                      'MP%s TOTAL ENERGY' % (str(mplevel2))

    # Integrate CFOUR methods
    VARH.update(cfour_psivar_list())
    return VARH

VARH = return_energy_components()


###################################
##  Start of Complete Basis Set  ##
###################################

def cbs(func, label, **kwargs):
    r"""Function to define a multistage energy method from combinations of
    basis set extrapolations and delta corrections and condense the
    components into a minimum number of calculations.

    :aliases: complete_basis_set()

    :returns: (*float*) -- Total electronic energy in Hartrees

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`CBS TOTAL ENERGY <CBSTOTALENERGY>`
       * :psivar:`CBS REFERENCE ENERGY <CBSREFERENCEENERGY>`
       * :psivar:`CBS CORRELATION ENERGY <CBSCORRELATIONENERGY>`
       * :psivar:`CURRENT ENERGY <CURRENTENERGY>`
       * :psivar:`CURRENT REFERENCE ENERGY <CURRENTREFERENCEENERGY>`
       * :psivar:`CURRENT CORRELATION ENERGY <CURRENTCORRELATIONENERGY>`

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - No way to tell function to boost fitting basis size for all calculations.

       - No way to extrapolate def2 family basis sets

       - Need to add more extrapolation schemes

    As represented in the equation below, a CBS energy method is defined in several
    sequential stages (scf, corl, delta, delta2, delta3, delta4, delta5) covering treatment
    of the reference total energy, the correlation energy, a delta correction to the
    correlation energy, and a second delta correction, etc.. Each is activated by its
    stage_wfn keyword and is only allowed if all preceding stages are active.

    .. include:: ../cbs_eqn.rst

    * Energy Methods
        The presence of a stage_wfn keyword is the indicator to incorporate
        (and check for stage_basis and stage_scheme keywords) and compute
        that stage in defining the CBS energy.

        The cbs() function requires, at a minimum, ``name='scf'`` and ``scf_basis``
        keywords to be specified for reference-step only jobs and ``name`` and
        ``corl_basis`` keywords for correlated jobs.

        The following energy methods have been set up for cbs().

        .. hlist::
           :columns: 5

           * scf
           * hf
           * mp2
           * mp2.5
           * mp3
           * mp4(sdq)
           * mp4
           * mp\ *n*
           * omp2
           * omp2.5
           * omp3
           * olccd
           * lccd
           * lccsd
           * cepa(0)
           * cepa(1)
           * cepa(3)
           * acpf
           * aqcc
           * qcisd
           * cc2
           * ccsd
           * fno-ccsd
           * bccd
           * cc3
           * qcisd(t)
           * ccsd(t)
           * fno-ccsd(t)
           * bccd(t)
           * cisd
           * cisdt
           * cisdtq
           * ci\ *n*
           * fci
           * mrccsd
           * mrccsd(t)
           * mrccsdt
           * mrccsdt(q)

    :type name: string
    :param name: ``'scf'`` || ``'ccsd'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        for the correlation energy, unless only reference step to be performed,
        in which case should be ``'scf'``. Overruled if stage_wfn keywords supplied.

    :type scf_wfn: string
    :param scf_wfn: |dl| ``'scf'`` |dr| || ``'c4-scf'`` || etc.

        Indicates the energy method for which the reference energy is to be
        obtained. Generally unnecessary, as 'scf' is *the* scf in |PSIfour| but
        can be used to direct lone scf components to run in |PSIfour| or Cfour
        in a mixed-program composite method.

    :type corl_wfn: string
    :param corl_wfn: ``'mp2'`` || ``'ccsd(t)'`` || etc.

        Indicates the energy method for which the correlation energy is to be
        obtained. Can also be specified with ``name`` or as the unlabeled
        first argument to the function.

    :type delta_wfn: string
    :param delta_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a delta correction
        to the correlation energy is to be obtained.

    :type delta_wfn_lesser: string
    :param delta_wfn_lesser: |dl| ``corl_wfn`` |dr| || ``'mp2'`` || etc.

        Indicates the inferior energy method for which a delta correction
        to the correlation energy is to be obtained.

    :type delta2_wfn: string
    :param delta2_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a second delta correction
        to the correlation energy is to be obtained.

    :type delta2_wfn_lesser: string
    :param delta2_wfn_lesser: |dl| ``delta_wfn`` |dr| || ``'ccsd(t)'`` || etc.

        Indicates the inferior energy method for which a second delta correction
        to the correlation energy is to be obtained.

    :type delta3_wfn: string
    :param delta3_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a third delta correction
        to the correlation energy is to be obtained.

    :type delta3_wfn_lesser: string
    :param delta3_wfn_lesser: |dl| ``delta2_wfn`` |dr| || ``'ccsd(t)'`` || etc.

        Indicates the inferior energy method for which a third delta correction
        to the correlation energy is to be obtained.

    :type delta4_wfn: string
    :param delta4_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a fourth delta correction
        to the correlation energy is to be obtained.

    :type delta4_wfn_lesser: string
    :param delta4_wfn_lesser: |dl| ``delta3_wfn`` |dr| || ``'ccsd(t)'`` || etc.

        Indicates the inferior energy method for which a fourth delta correction
        to the correlation energy is to be obtained.

    :type delta5_wfn: string
    :param delta5_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a fifth delta correction
        to the correlation energy is to be obtained.

    :type delta5_wfn_lesser: string
    :param delta5_wfn_lesser: |dl| ``delta4_wfn`` |dr| || ``'ccsd(t)'`` || etc.

        Indicates the inferior energy method for which a fifth delta correction
        to the correlation energy is to be obtained.

    * Basis Sets
        Currently, the basis set set through ``set`` commands have no influence
        on a cbs calculation.

    :type scf_basis: :ref:`basis string <apdx:basisElement>`
    :param scf_basis: |dl| ``corl_basis`` |dr| || ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the reference energy.
        If any correlation method is specified, ``scf_basis`` can default
        to ``corl_basis``.

    :type corl_basis: :ref:`basis string <apdx:basisElement>`
    :param corl_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the correlation energy.

    :type delta_basis: :ref:`basis string <apdx:basisElement>`
    :param delta_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the delta correction
        to the correlation energy.

    :type delta2_basis: :ref:`basis string <apdx:basisElement>`
    :param delta2_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the second delta correction
        to the correlation energy.

    :type delta3_basis: :ref:`basis string <apdx:basisElement>`
    :param delta3_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the third delta correction
        to the correlation energy.

    :type delta4_basis: :ref:`basis string <apdx:basisElement>`
    :param delta4_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the fourth delta correction
        to the correlation energy.

    :type delta5_basis: :ref:`basis string <apdx:basisElement>`
    :param delta5_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the fifth delta correction
        to the correlation energy.

    * Schemes
        Transformations of the energy through basis set extrapolation for each
        stage of the CBS definition. A complaint is generated if number of basis
        sets in stage_basis does not exactly satisfy requirements of stage_scheme.
        An exception is the default, ``'xtpl_highest_1'``, which uses the best basis
        set available. See :ref:`sec:cbs_xtpl` for all available schemes.

    :type scf_scheme: function
    :param scf_scheme: |dl| ``xtpl_highest_1`` |dr| || ``scf_xtpl_helgaker_3`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the reference energy.
        Defaults to :py:func:`~scf_xtpl_helgaker_3` if three valid basis sets
        present in ``scf_basis``, :py:func:`~scf_xtpl_helgaker_2` if two valid basis
        sets present in ``scf_basis``, and :py:func:`~xtpl_highest_1` otherwise.

    :type corl_scheme: function
    :param corl_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``corl_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type delta_scheme: function
    :param delta_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the delta correction
        to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type delta2_scheme: function
    :param delta2_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the second delta correction
        to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta2_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type delta3_scheme: function
    :param delta3_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the third delta correction
        to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta3_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type delta4_scheme: function
    :param delta4_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the fourth delta correction
        to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta4_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type delta5_scheme: function
    :param delta5_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the fifth delta correction
        to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta5_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :examples:


    >>> # [1] replicates with cbs() the simple model chemistry scf/cc-pVDZ: set basis cc-pVDZ energy('scf')
    >>> cbs(name='scf', scf_basis='cc-pVDZ')

    >>> # [2] replicates with cbs() the simple model chemistry mp2/jun-cc-pVDZ: set basis jun-cc-pVDZ energy('mp2')
    >>> cbs(name='mp2', corl_basis='jun-cc-pVDZ')

    >>> # [3] DTQ-zeta extrapolated scf reference energy
    >>> cbs(name='scf', scf_basis='cc-pV[DTQ]Z', scf_scheme=scf_xtpl_helgaker_3)

    >>> # [4] DT-zeta extrapolated mp2 correlation energy atop a T-zeta reference
    >>> cbs(corl_wfn='mp2', corl_basis='cc-pv[dt]z', corl_scheme=corl_xtpl_helgaker_2)

    >>> # [5] a DT-zeta extrapolated coupled-cluster correction atop a TQ-zeta extrapolated mp2 correlation energy atop a Q-zeta reference (both equivalent)
    >>> cbs(corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z')
    >>> cbs(energy, wfn='mp2', corl_basis='aug-cc-pv[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z', delta_scheme=corl_xtpl_helgaker_2)

    >>> # [6] a D-zeta ccsd(t) correction atop a DT-zeta extrapolated ccsd cluster correction atop a TQ-zeta extrapolated mp2 correlation energy atop a Q-zeta reference
    >>> cbs(name='mp2', corl_basis='aug-cc-pv[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd', delta_basis='aug-cc-pv[dt]z', delta_scheme=corl_xtpl_helgaker_2, delta2_wfn='ccsd(t)', delta2_wfn_lesser='ccsd', delta2_basis='aug-cc-pvdz')

    >>> # [7] cbs() coupled with database()
    >>> TODO database('mp2', 'BASIC', subset=['h2o','nh3'], symm='on', func=cbs, corl_basis='cc-pV[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd(t)', delta_basis='sto-3g')

    >>> # [8] cbs() coupled with optimize()
    >>> TODO optimize('mp2', corl_basis='cc-pV[DT]Z', corl_scheme=corl_xtpl_helgaker_2, func=cbs)

    """
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    verbose = kwargs.pop('verbose', 0)
    ptype = kwargs.pop('ptype')

    # Establish function to call (only energy makes sense for cbs)
    if ptype not in ['energy', 'gradient', 'hessian']:
        raise ValidationError("""Wrapper complete_basis_set is unhappy to be calling function '%s' instead of 'energy'.""" % ptype)

    optstash = p4util.OptionsState(
        ['BASIS'],
        ['WFN'],
        ['WRITER_FILE_LABEL'])

    # Define some quantum chemical knowledge, namely what methods are subsumed in others

    do_scf = True
    do_corl = False
    do_delta = False
    do_delta2 = False
    do_delta3 = False
    do_delta4 = False
    do_delta5 = False

    user_writer_file_label = core.get_global_option('WRITER_FILE_LABEL')

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()
    molstr = molecule.create_psi4_string_from_molecule()
    natom = molecule.natom()

    # Establish method for reference energy
    cbs_scf_wfn = kwargs.pop('scf_wfn', 'hf').lower()

    if do_scf:
        if cbs_scf_wfn not in VARH.keys():
            raise ValidationError("""Requested SCF method '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" % (cbs_scf_wfn))

    # Establish method for correlation energy
    cbs_corl_wfn = kwargs.pop('corl_wfn', '').lower()
    if cbs_corl_wfn:
        do_corl = True

    if do_corl:
        if cbs_corl_wfn not in VARH.keys():
            raise ValidationError("""Requested CORL method '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" % (cbs_corl_wfn))

        cbs_corl_wfn_lesser = kwargs.get('corl_wfn_lesser', cbs_scf_wfn).lower()
        if cbs_corl_wfn_lesser not in VARH.keys():
            raise ValidationError("""Requested CORL method lesser '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" % (cbs_delta_wfn_lesser))

    # Establish method for delta correction energy
    if 'delta_wfn' in kwargs:
        do_delta = True
        cbs_delta_wfn = kwargs['delta_wfn'].lower()
        if cbs_delta_wfn not in VARH.keys():
            raise ValidationError("""Requested DELTA method '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" % (cbs_delta_wfn))

        cbs_delta_wfn_lesser = kwargs.get('delta_wfn_lesser', cbs_corl_wfn).lower()
        if cbs_delta_wfn_lesser not in VARH.keys():
            raise ValidationError("""Requested DELTA method lesser '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" % (cbs_delta_wfn_lesser))

    # Establish method for second delta correction energy
    if 'delta2_wfn' in kwargs:
        do_delta2 = True
        cbs_delta2_wfn = kwargs['delta2_wfn'].lower()
        if cbs_delta2_wfn not in VARH.keys():
            raise ValidationError("""Requested DELTA2 method '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" % (cbs_delta2_wfn))

        cbs_delta2_wfn_lesser = kwargs.get('delta2_wfn_lesser', cbs_delta_wfn).lower()
        if cbs_delta2_wfn_lesser not in VARH.keys():
            raise ValidationError("""Requested DELTA2 method lesser '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" % (cbs_delta2_wfn_lesser))

#    # Establish method for third delta correction energy
#    if 'delta3_wfn' in kwargs:
#        do_delta3 = True
#        cbs_delta3_wfn = kwargs['delta3_wfn'].lower()
#        if cbs_delta3_wfn not in VARH.keys():
#            raise ValidationError("""Requested DELTA3 method '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" % (cbs_delta3_wfn))
#
#        cbs_delta3_wfn_lesser = kwargs.get('delta3_wfn_lesser', cbs_delta2_wfn).lower()
#        if not (cbs_delta3_wfn_lesser in VARH.keys()):
#            raise ValidationError("""Requested DELTA3 method lesser '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" % (cbs_delta3_wfn_lesser))
#
#    # Establish method for fourth delta correction energy
#    if 'delta4_wfn' in kwargs:
#        do_delta4 = True
#        cbs_delta4_wfn = kwargs['delta4_wfn'].lower()
#        if not (cbs_delta4_wfn in VARH.keys()):
#            raise ValidationError('Requested DELTA4 method \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta4_wfn))
#
#        if 'delta4_wfn_lesser' in kwargs:
#            cbs_delta4_wfn_lesser = kwargs['delta4_wfn_lesser'].lower()
#        else:
#            cbs_delta4_wfn_lesser = cbs_delta3_wfn
#        if not (cbs_delta4_wfn_lesser in VARH.keys()):
#            raise ValidationError('Requested DELTA4 method lesser \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta4_wfn_lesser))
#
#    # Establish method for fifth delta correction energy
#    if 'delta5_wfn' in kwargs:
#        do_delta5 = True
#        cbs_delta5_wfn = kwargs['delta5_wfn'].lower()
#        if not (cbs_delta5_wfn in VARH.keys()):
#            raise ValidationError('Requested DELTA5 method \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta5_wfn))
#
#        if 'delta5_wfn_lesser' in kwargs:
#            cbs_delta5_wfn_lesser = kwargs['delta5_wfn_lesser'].lower()
#        else:
#            cbs_delta5_wfn_lesser = cbs_delta4_wfn
#        if not (cbs_delta5_wfn_lesser in VARH.keys()):
#            raise ValidationError('Requested DELTA5 method lesser \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta5_wfn_lesser))

    # Check that user isn't skipping steps in scf + corl + delta + delta2 sequence
    if   do_scf and not do_corl and not do_delta and not do_delta2 and not do_delta3 and not do_delta4 and not do_delta5:
        pass
    elif do_scf and do_corl and not do_delta and not do_delta2 and not do_delta3 and not do_delta4 and not do_delta5:
        pass
    elif do_scf and do_corl and do_delta and not do_delta2 and not do_delta3 and not do_delta4 and not do_delta5:
        pass
    elif do_scf and do_corl and do_delta and do_delta2 and not do_delta3 and not do_delta4 and not do_delta5:
        pass
    #elif do_scf and do_corl and do_delta and do_delta2 and do_delta3 and not do_delta4 and not do_delta5:
    #    pass
    #elif do_scf and do_corl and do_delta and do_delta2 and do_delta3 and do_delta4 and not do_delta5:
    #    pass
    #elif do_scf and do_corl and do_delta and do_delta2 and do_delta3 and do_delta4 and do_delta5:
    #    pass
    else:
        raise ValidationError('Requested scf (%s) + corl (%s) + delta (%s) + delta2 (%s) + delta3 (%s) + delta4 (%s) + delta5 (%s) not valid. These steps are cummulative.' %
            (do_scf, do_corl, do_delta, do_delta2, do_delta3, do_delta4, do_delta5))

    # Establish list of valid basis sets for correlation energy
    if do_corl:
        if 'corl_basis' in kwargs:
            BSTC, ZETC = _expand_bracketed_basis(kwargs['corl_basis'].lower(), molecule=molstr)
        else:
            raise ValidationError("""CORL basis sets through keyword '%s' are required.""" % ('corl_basis'))

    # Establish list of valid basis sets for scf energy
    if 'scf_basis' in kwargs:
        BSTR, ZETR = _expand_bracketed_basis(kwargs['scf_basis'].lower(), molecule=molstr)
    elif do_corl:
        BSTR = BSTC[:]
        ZETR = ZETC[:]
    else:
        raise ValidationError("""SCF basis sets through keyword '%s' are required. Or perhaps you forgot the '%s'.""" % ('scf_basis', 'corl_wfn'))

    # Establish list of valid basis sets for delta correction energy
    if do_delta:
        if 'delta_basis' in kwargs:
            BSTD, ZETD = _expand_bracketed_basis(kwargs['delta_basis'].lower(), molecule=molstr)
        else:
            raise ValidationError("""DELTA basis sets through keyword '%s' are required.""" % ('delta_basis'))

    # Establish list of valid basis sets for second delta correction energy
    if do_delta2:
        if 'delta2_basis' in kwargs:
            BSTD2, ZETD2 = _expand_bracketed_basis(kwargs['delta2_basis'].lower(), molecule=molstr)
        else:
            raise ValidationError("""DELTA2 basis sets through keyword '%s' are required.""" % ('delta2_basis'))

#    # Establish list of valid basis sets for third delta correction energy
#    if do_delta3:
#        if 'delta3_basis' in kwargs:
#            BSTD3, ZETD3 = validate_bracketed_basis(kwargs['delta3_basis'].lower())
#        else:
#            raise ValidationError('DELTA3 basis sets through keyword \'%s\' are required.' % ('delta3_basis'))
#
#    # Establish list of valid basis sets for fourth delta correction energy
#    if do_delta4:
#        if 'delta4_basis' in kwargs:
#            BSTD4, ZETD4 = validate_bracketed_basis(kwargs['delta4_basis'].lower())
#        else:
#            raise ValidationError('DELTA4 basis sets through keyword \'%s\' are required.' % ('delta4_basis'))
#
#    # Establish list of valid basis sets for fifth delta correction energy
#    if do_delta5:
#        if 'delta5_basis' in kwargs:
#            BSTD5, ZETD5 = validate_bracketed_basis(kwargs['delta5_basis'].lower())
#        else:
#            raise ValidationError('DELTA5 basis sets through keyword \'%s\' are required.' % ('delta5_basis'))

    # Establish treatment for scf energy (validity check useless since python will catch it long before here)
    if (len(BSTR) == 3) and ('scf_basis' in kwargs):
        cbs_scf_scheme = scf_xtpl_helgaker_3
    elif (len(BSTR) == 2) and ('scf_basis' in kwargs):
        cbs_scf_scheme = scf_xtpl_helgaker_2
    elif (len(BSTR) == 1) and ('scf_basis' in kwargs):
        cbs_scf_scheme = xtpl_highest_1
    elif 'scf_basis' in kwargs:
        raise ValidationError("""SCF basis sets of number %d cannot be handled.""" % (len(BSTR)))
    elif do_corl:
        cbs_scf_scheme = xtpl_highest_1
        BSTR = [BSTC[-1]]
        ZETR = [ZETC[-1]]
    if 'scf_scheme' in kwargs:
        cbs_scf_scheme = kwargs['scf_scheme']

    # Establish treatment for correlation energy
    if do_corl:
        if len(BSTC) == 2:
            cbs_corl_scheme = corl_xtpl_helgaker_2
        else:
            cbs_corl_scheme = xtpl_highest_1
        if 'corl_scheme' in kwargs:
            cbs_corl_scheme = kwargs['corl_scheme']

    # Establish treatment for delta correction energy
    if do_delta:
        if len(BSTD) == 2:
            cbs_delta_scheme = corl_xtpl_helgaker_2
        else:
            cbs_delta_scheme = xtpl_highest_1
        if 'delta_scheme' in kwargs:
            cbs_delta_scheme = kwargs['delta_scheme']

    # Establish treatment for delta2 correction energy
    if do_delta2:
        if len(BSTD2) == 2:
            cbs_delta2_scheme = corl_xtpl_helgaker_2
        else:
            cbs_delta2_scheme = xtpl_highest_1
        if 'delta2_scheme' in kwargs:
            cbs_delta2_scheme = kwargs['delta2_scheme']

#    # Establish treatment for delta3 correction energy
#    if do_delta3:
#        if len(BSTD3) == 2:
#            cbs_delta3_scheme = corl_xtpl_helgaker_2
#        else:
#            cbs_delta3_scheme = xtpl_highest_1
#        if 'delta3_scheme' in kwargs:
#            cbs_delta3_scheme = kwargs['delta3_scheme']
#
#    # Establish treatment for delta4 correction energy
#    if do_delta4:
#        if len(BSTD4) == 2:
#            cbs_delta4_scheme = corl_xtpl_helgaker_2
#        else:
#            cbs_delta4_scheme = xtpl_highest_1
#        if 'delta4_scheme' in kwargs:
#            cbs_delta4_scheme = kwargs['delta4_scheme']
#
#    # Establish treatment for delta5 correction energy
#    if do_delta5:
#        if len(BSTD5) == 2:
#            cbs_delta5_scheme = corl_xtpl_helgaker_2
#        else:
#            cbs_delta5_scheme = xtpl_highest_1
#        if 'delta5_scheme' in kwargs:
#            cbs_delta5_scheme = kwargs['delta5_scheme']

    # Build string of title banner
    cbsbanners = ''
    cbsbanners += """core.print_out('\\n')\n"""
    cbsbanners += """p4util.banner(' CBS Setup: %s ' % label)\n"""
    cbsbanners += """core.print_out('\\n')\n\n"""
    exec(cbsbanners)

    # Call schemes for each portion of total energy to 'place orders' for calculations needed
    d_fields = ['d_stage', 'd_scheme', 'd_basis', 'd_wfn', 'd_need', 'd_coef', 'd_energy', 'd_gradient', 'd_hessian']
    f_fields = ['f_wfn', 'f_basis', 'f_zeta', 'f_energy', 'f_gradient', 'f_hessian']
    GRAND_NEED = []
    MODELCHEM = []
    if do_scf:
        NEED = _expand_scheme_orders(cbs_scf_scheme, BSTR, ZETR, cbs_scf_wfn, natom)
        GRAND_NEED.append(dict(zip(d_fields, ['scf', cbs_scf_scheme,
            _contract_bracketed_basis(BSTR), cbs_scf_wfn, NEED, +1, 0.0, None, None])))

    if do_corl:
        NEED = _expand_scheme_orders(cbs_corl_scheme, BSTC, ZETC, cbs_corl_wfn, natom)
        GRAND_NEED.append(dict(zip(d_fields, ['corl', cbs_corl_scheme,
            _contract_bracketed_basis(BSTC), cbs_corl_wfn, NEED, +1, 0.0, None, None])))

        NEED = _expand_scheme_orders(cbs_corl_scheme, BSTC, ZETC, cbs_corl_wfn_lesser, natom)
        GRAND_NEED.append(dict(zip(d_fields, ['corl', cbs_corl_scheme,
            _contract_bracketed_basis(BSTC), cbs_corl_wfn_lesser, NEED, -1, 0.0, None, None])))

    if do_delta:
        NEED = _expand_scheme_orders(cbs_delta_scheme, BSTD, ZETD, cbs_delta_wfn, natom)
        GRAND_NEED.append(dict(zip(d_fields, ['delta', cbs_delta_scheme,
            _contract_bracketed_basis(BSTD), cbs_delta_wfn, NEED, +1, 0.0, None, None])))

        NEED = _expand_scheme_orders(cbs_delta_scheme, BSTD, ZETD, cbs_delta_wfn_lesser, natom)
        GRAND_NEED.append(dict(zip(d_fields, ['delta', cbs_delta_scheme,
            _contract_bracketed_basis(BSTD), cbs_delta_wfn_lesser, NEED, -1, 0.0, None, None])))

    if do_delta2:
        NEED = _expand_scheme_orders(cbs_delta2_scheme, BSTD2, ZETD2, cbs_delta2_wfn, natom)
        GRAND_NEED.append(dict(zip(d_fields, ['delta2', cbs_delta2_scheme,
            _contract_bracketed_basis(BSTD2), cbs_delta2_wfn, NEED, +1, 0.0, None, None])))

        NEED = _expand_scheme_orders(cbs_delta2_scheme, BSTD2, ZETD2, cbs_delta2_wfn_lesser, natom)
        GRAND_NEED.append(dict(zip(d_fields, ['delta2', cbs_delta2_scheme,
            _contract_bracketed_basis(BSTD2), cbs_delta2_wfn_lesser, NEED, -1, 0.0, None, None])))

#    if do_delta3:
#        NEED = call_function_in_1st_argument(cbs_delta3_scheme,
#            mode='requisition', basisname=BSTD3, basiszeta=ZETD3, wfnname=cbs_delta3_wfn)
#        GRAND_NEED.append(dict(zip(d_fields, ['delta3', cbs_delta3_scheme,
#            reconstitute_bracketed_basis(NEED), cbs_delta3_wfn, NEED, +1, 0.0])))
#
#        NEED = call_function_in_1st_argument(cbs_delta3_scheme,
#            mode='requisition', basisname=BSTD3, basiszeta=ZETD3, wfnname=cbs_delta3_wfn_lesser)
#        GRAND_NEED.append(dict(zip(d_fields, ['delta3', cbs_delta3_scheme,
#            reconstitute_bracketed_basis(NEED), cbs_delta3_wfn_lesser, NEED, -1, 0.0])))
#
#    if do_delta4:
#        NEED = call_function_in_1st_argument(cbs_delta4_scheme,
#            mode='requisition', basisname=BSTD4, basiszeta=ZETD4, wfnname=cbs_delta4_wfn)
#        GRAND_NEED.append(dict(zip(d_fields, ['delta4', cbs_delta4_scheme,
#            reconstitute_bracketed_basis(NEED), cbs_delta4_wfn, NEED, +1, 0.0])))
#
#        NEED = call_function_in_1st_argument(cbs_delta4_scheme,
#            mode='requisition', basisname=BSTD4, basiszeta=ZETD4, wfnname=cbs_delta4_wfn_lesser)
#        GRAND_NEED.append(dict(zip(d_fields, ['delta4', cbs_delta4_scheme,
#            reconstitute_bracketed_basis(NEED), cbs_delta4_wfn_lesser, NEED, -1, 0.0])))
#
#    if do_delta5:
#        NEED = call_function_in_1st_argument(cbs_delta5_scheme,
#            mode='requisition', basisname=BSTD5, basiszeta=ZETD5, wfnname=cbs_delta5_wfn)
#        GRAND_NEED.append(dict(zip(d_fields, ['delta5', cbs_delta5_scheme,
#            reconstitute_bracketed_basis(NEED), cbs_delta5_wfn, NEED, +1, 0.0])))
#
#        NEED = call_function_in_1st_argument(cbs_delta5_scheme,
#            mode='requisition', basisname=BSTD5, basiszeta=ZETD5, wfnname=cbs_delta5_wfn_lesser)
#        GRAND_NEED.append(dict(zip(d_fields, ['delta5', cbs_delta5_scheme,
#            reconstitute_bracketed_basis(NEED), cbs_delta5_wfn_lesser, NEED, -1, 0.0])))

    for stage in GRAND_NEED:
        for lvl in stage['d_need'].items():
            MODELCHEM.append(lvl[1])

    # Apply chemical reasoning to choose the minimum computations to run
    JOBS = MODELCHEM[:]

    addlremark = {'energy': '', 'gradient': ', GRADIENT', 'hessian': ', HESSIAN'}
    instructions = ''
    instructions += """    Naive listing of computations required.\n"""
    for mc in JOBS:
        instructions += """   %12s / %-24s for  %s%s\n""" % \
            (mc['f_wfn'], mc['f_basis'], VARH[mc['f_wfn']][mc['f_wfn']], addlremark[ptype])

    #     Remove duplicate modelchem portion listings
    for mc in MODELCHEM:
        dups = -1
        for indx_job, job in enumerate(JOBS):
            if (job['f_wfn'] == mc['f_wfn']) and (job['f_basis'] == mc['f_basis']):
                dups += 1
                if dups >= 1:
                    del JOBS[indx_job]

    #     Remove chemically subsumed modelchem portion listings
    if ptype == 'energy':
        for mc in MODELCHEM:
            for wfn in VARH[mc['f_wfn']]:
                for indx_job, job in enumerate(JOBS):
                    if (VARH[mc['f_wfn']][wfn] == VARH[job['f_wfn']][job['f_wfn']]) and \
                       (mc['f_basis'] == job['f_basis']) and not \
                       (mc['f_wfn'] == job['f_wfn']):
                        del JOBS[indx_job]

    instructions += """\n    Enlightened listing of computations required.\n"""
    for mc in JOBS:
        instructions += """   %12s / %-24s for  %s%s\n""" % \
            (mc['f_wfn'], mc['f_basis'], VARH[mc['f_wfn']][mc['f_wfn']], addlremark[ptype])

    #     Expand listings to all that will be obtained
    JOBS_EXT = []
    for job in JOBS:
        for wfn in VARH[job['f_wfn']]:
            JOBS_EXT.append(dict(zip(f_fields, [wfn, job['f_basis'], job['f_zeta'],
                                                0.0,
                                                core.Matrix(natom, 3),
                                                core.Matrix(3 * natom, 3 * natom)])))

    instructions += """\n    Full listing of computations to be obtained (required and bonus).\n"""
    for mc in JOBS_EXT:
        instructions += """   %12s / %-24s for  %s%s\n""" % \
            (mc['f_wfn'], mc['f_basis'], VARH[mc['f_wfn']][mc['f_wfn']], addlremark[ptype])
    core.print_out(instructions)

    psioh = core.IOManager.shared_object()
    psioh.set_specific_retention(constants.PSIF_SCF_MOS, True)
    # projection across point groups not allowed and cbs() usually a mix of symm-enabled and symm-tol calls
    #   needs to be communicated to optimize() so reset by that optstash
    core.set_local_option('SCF', 'GUESS_PERSIST', True)

    Njobs = 0
    # Run necessary computations
    for mc in JOBS:
        kwargs['name'] = mc['f_wfn']

        # Build string of title banner
        cbsbanners = ''
        cbsbanners += """core.print_out('\\n')\n"""
        cbsbanners += """p4util.banner(' CBS Computation: %s / %s%s ')\n""" % \
            (mc['f_wfn'].upper(), mc['f_basis'].upper(), addlremark[ptype])
        cbsbanners += """core.print_out('\\n')\n\n"""
        exec(cbsbanners)

        # Build string of molecule and commands that are dependent on the database
        commands = '\n'
        commands += """\ncore.set_global_option('BASIS', '%s')\n""" % (mc['f_basis'])
        commands += """core.set_global_option('WRITER_FILE_LABEL', '%s')\n""" % \
            (user_writer_file_label + ('' if user_writer_file_label == '' else '-') + mc['f_wfn'].lower() + '-' + mc['f_basis'].lower())
        exec(commands)

        # Make energy(), etc. call
        response = func(molecule=molecule, **kwargs)
        if ptype == 'energy':
            mc['f_energy'] = response
        elif ptype == 'gradient':
            mc['f_gradient'] = response
            mc['f_energy'] = core.get_variable('CURRENT ENERGY')
            if verbose > 1:
                mc['f_gradient'].print_out()
        elif ptype == 'hessian':
            mc['f_hessian'] = response
            mc['f_energy'] = core.get_variable('CURRENT ENERGY')
            if verbose > 1:
                mc['f_hessian'].print_out()
        Njobs += 1
        if verbose > 1:
            core.print_out("\nCURRENT ENERGY: %14.16f\n" % mc['f_energy'])

        # Fill in energies for subsumed methods
        if ptype == 'energy':
            for wfn in VARH[mc['f_wfn']]:
                for job in JOBS_EXT:
                    if (wfn == job['f_wfn']) and (mc['f_basis'] == job['f_basis']):
                        job['f_energy'] = core.get_variable(VARH[wfn][wfn])

        if verbose > 1:
            core.print_variables()
        core.clean_variables()
        core.clean()

        # Copy data from 'run' to 'obtained' table
        for mce in JOBS_EXT:
            if (mc['f_wfn'] == mce['f_wfn']) and (mc['f_basis'] == mce['f_basis']):
                mce['f_energy'] = mc['f_energy']
                mce['f_gradient'] = mc['f_gradient']
                mce['f_hessian'] = mc['f_hessian']

    psioh.set_specific_retention(constants.PSIF_SCF_MOS, False)

    # Build string of title banner
    cbsbanners = ''
    cbsbanners += """core.print_out('\\n')\n"""
    cbsbanners += """p4util.banner(' CBS Results: %s ' % label)\n"""
    cbsbanners += """core.print_out('\\n')\n\n"""
    exec(cbsbanners)

    # Insert obtained energies into the array that stores the cbs stages
    for stage in GRAND_NEED:
        for lvl in stage['d_need'].items():
            MODELCHEM.append(lvl[1])

            for job in JOBS_EXT:
                # Dont ask
                if (((lvl[1]['f_wfn'] == job['f_wfn']) or
                     ((lvl[1]['f_wfn'][3:] == job['f_wfn']) and lvl[1]['f_wfn'].startswith('c4-')) or
                     ((lvl[1]['f_wfn'] == job['f_wfn'][3:]) and job['f_wfn'].startswith('c4-')) or
                     (('c4-' + lvl[1]['f_wfn']) == job['f_wfn']) or
                     (lvl[1]['f_wfn'] == ('c4-' + job['f_wfn']))) and
                     (lvl[1]['f_basis'] == job['f_basis'])):
                    lvl[1]['f_energy'] = job['f_energy']
                    lvl[1]['f_gradient'] = job['f_gradient']
                    lvl[1]['f_hessian'] = job['f_hessian']

    # Make xtpl() call
    finalenergy = 0.0
    finalgradient = core.Matrix(natom, 3)
    finalhessian = core.Matrix(3 * natom, 3 * natom)
    for stage in GRAND_NEED:
        hiloargs = _contract_scheme_orders(stage['d_need'], 'f_energy')
        stage['d_energy'] = stage['d_scheme'](**hiloargs)
        finalenergy += stage['d_energy'] * stage['d_coef']

        if ptype == 'gradient':
            hiloargs = _contract_scheme_orders(stage['d_need'], 'f_gradient')
            stage['d_gradient'] = stage['d_scheme'](**hiloargs)
            work = stage['d_gradient'].clone()
            work.scale(stage['d_coef'])
            finalgradient.add(work)

        elif ptype == 'hessian':
            hiloargs = _contract_scheme_orders(stage['d_need'], 'f_hessian')
            stage['d_hessian'] = stage['d_scheme'](**hiloargs)
            work = stage['d_hessian'].clone()
            work.scale(stage['d_coef'])
            finalhessian.add(work)

    # Build string of results table
    table_delimit = '  ' + '-' * 105 + '\n'
    tables = ''
    tables += """\n   ==> %s <==\n\n""" % ('Components')
    tables += table_delimit
    tables += """     %6s %20s %1s %-26s %3s %16s   %-s\n""" % ('', 'Method', '/', 'Basis', 'Rqd', 'Energy [Eh]', 'Variable')
    tables += table_delimit
    for job in JOBS_EXT:
        star = ''
        for mc in MODELCHEM:
            if (job['f_wfn'] == mc['f_wfn']) and (job['f_basis'] == mc['f_basis']):
                star = '*'
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % ('', job['f_wfn'],
                  '/', job['f_basis'], star, job['f_energy'], VARH[job['f_wfn']][job['f_wfn']])
    tables += table_delimit

    tables += """\n   ==> %s <==\n\n""" % ('Stages')
    tables += table_delimit
    tables += """     %6s %20s %1s %-27s %2s %16s   %-s\n""" % ('Stage', 'Method', '/', 'Basis', 'Wt', 'Energy [Eh]', 'Scheme')
    tables += table_delimit
    for stage in GRAND_NEED:
        tables += """     %6s %20s %1s %-27s %2d %16.8f   %-s\n""" % (stage['d_stage'], stage['d_wfn'],
                  '/', stage['d_basis'], stage['d_coef'], stage['d_energy'], stage['d_scheme'].__name__)
    tables += table_delimit

    tables += """\n   ==> %s <==\n\n""" % ('CBS')
    tables += table_delimit
    tables += """     %6s %20s %1s %-27s %2s %16s   %-s\n""" % ('Stage', 'Method', '/', 'Basis', '', 'Energy [Eh]', 'Scheme')
    tables += table_delimit
    if do_scf:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[0]['d_stage'],
                                                                      GRAND_NEED[0]['d_wfn'], '/', GRAND_NEED[0]['d_basis'], '',
                                                                      GRAND_NEED[0]['d_energy'],
                                                                      GRAND_NEED[0]['d_scheme'].__name__)
    if do_corl:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[1]['d_stage'],
                                                                      GRAND_NEED[1]['d_wfn'], '/', GRAND_NEED[1]['d_basis'], '',
                                                                      GRAND_NEED[1]['d_energy'] - GRAND_NEED[2]['d_energy'],
                                                                      GRAND_NEED[1]['d_scheme'].__name__)
    if do_delta:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[3]['d_stage'],
                                                                      GRAND_NEED[3]['d_wfn'] + ' - ' + GRAND_NEED[4]['d_wfn'], '/', GRAND_NEED[3]['d_basis'], '',
                                                                      GRAND_NEED[3]['d_energy'] - GRAND_NEED[4]['d_energy'],
                                                                      GRAND_NEED[3]['d_scheme'].__name__)
    if do_delta2:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[5]['d_stage'],
                                                                      GRAND_NEED[5]['d_wfn'] + ' - ' + GRAND_NEED[6]['d_wfn'], '/', GRAND_NEED[5]['d_basis'], '',
                                                                      GRAND_NEED[5]['d_energy'] - GRAND_NEED[6]['d_energy'],
                                                                      GRAND_NEED[5]['d_scheme'].__name__)
#    if do_delta3:
#        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[6]['d_stage'], GRAND_NEED[6]['d_wfn'] + ' - ' + GRAND_NEED[7]['d_wfn'],
#                  '/', GRAND_NEED[6]['d_basis'], '', GRAND_NEED[6]['d_energy'] - GRAND_NEED[7]['d_energy'], GRAND_NEED[6]['d_scheme'].__name__)
#    if do_delta4:
#        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[8]['d_stage'], GRAND_NEED[8]['d_wfn'] + ' - ' + GRAND_NEED[9]['d_wfn'],
#                  '/', GRAND_NEED[8]['d_basis'], '', GRAND_NEED[8]['d_energy'] - GRAND_NEED[9]['d_energy'], GRAND_NEED[8]['d_scheme'].__name__)
#    if do_delta5:
#        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[10]['d_stage'], GRAND_NEED[10]['d_wfn'] + ' - ' + GRAND_NEED[11]['d_wfn'],
#                  '/', GRAND_NEED[10]['d_basis'], '', GRAND_NEED[10]['d_energy'] - GRAND_NEED[11]['d_energy'], GRAND_NEED[10]['d_scheme'].__name__)
    tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % ('total', 'CBS', '', '', '', finalenergy, '')
    tables += table_delimit

    core.print_out(tables)

    core.set_variable('CBS REFERENCE ENERGY', GRAND_NEED[0]['d_energy'])
    core.set_variable('CBS CORRELATION ENERGY', finalenergy - GRAND_NEED[0]['d_energy'])
    core.set_variable('CBS TOTAL ENERGY', finalenergy)
    core.set_variable('CURRENT REFERENCE ENERGY', GRAND_NEED[0]['d_energy'])
    core.set_variable('CURRENT CORRELATION ENERGY', finalenergy - GRAND_NEED[0]['d_energy'])
    core.set_variable('CURRENT ENERGY', finalenergy)
    core.set_variable('CBS NUMBER', Njobs)

    # new skeleton wavefunction w/mol, highest-SCF basis (just to choose one), & not energy
    basis = core.BasisSet.build(molecule, "ORBITAL", 'sto-3g')
    wfn = core.Wavefunction(molecule, basis)

    optstash.restore()

    if ptype == 'energy':
        finalquantity = finalenergy
    elif ptype == 'gradient':
        finalquantity = finalgradient
        wfn.set_gradient(finalquantity)
        if finalquantity.rows(0) < 20:
            core.print_out('CURRENT GRADIENT')
            finalquantity.print_out()
    elif ptype == 'hessian':
        finalquantity = finalhessian
        wfn.set_hessian(finalquantity)
        if finalquantity.rows(0) < 20:
            core.print_out('CURRENT HESSIAN')
            finalquantity.print_out()

    if return_wfn:
        return (finalquantity, wfn)
    else:
        return finalquantity


_lmh_labels = {1: ['HI'],
               2: ['LO', 'HI'],
               3: ['LO', 'MD', 'HI'],
               4: ['LO', 'MD', 'M2', 'HI'],
               5: ['LO', 'MD', 'M2', 'M3', 'HI']}


def _expand_scheme_orders(scheme, basisname, basiszeta, wfnname, natom):
    """Check that the length of *basiszeta* array matches the implied degree of
    extrapolation in *scheme* name. Return a dictionary of same length as
    basiszeta, with *basisname* and *basiszeta* distributed therein.

    """
    Nxtpl = len(basiszeta)

    if int(scheme.__name__.split('_')[-1]) != Nxtpl:
        raise ValidationError("""Call to '%s' not valid with '%s' basis sets.""" % (scheme.__name__, len(basiszeta)))

    f_fields = ['f_wfn', 'f_basis', 'f_zeta', 'f_energy', 'f_gradient', 'f_hessian']
    NEED = {}
    for idx in range(Nxtpl):
        NEED[_lmh_labels[Nxtpl][idx]] = dict(zip(f_fields, [wfnname, basisname[idx], basiszeta[idx],
                                                            0.0,
                                                            core.Matrix(natom, 3),
                                                            core.Matrix(3 * natom, 3 * natom)]))
    return NEED


def _contract_scheme_orders(needdict, datakey='f_energy'):
    """Prepared named arguments for extrapolation functions by
    extracting zetas and values (which one determined by *datakey*) out
    of *needdict* and returning a dictionary whose keys are contructed
    from _lmh_labels.

    """
    largs = {}
    largs['functionname'] = needdict['HI']['f_wfn']
    Nxtpl = len(needdict)
    zlabels = _lmh_labels[Nxtpl]  # e.g., ['LO', 'HI']

    for zeta in range(Nxtpl):
        zlab = zlabels[zeta]  # e.g., LO
        largs['z' + zlab] = needdict[zlab]['f_zeta']
        largs['value' + zlab] = needdict[zlab][datakey]

    return largs

##  Aliases  ##
complete_basis_set = cbs


def _cbs_wrapper_methods(**kwargs):
    cbs_method_kwargs = ['scf_wfn', 'corl_wfn', 'delta_wfn']
    cbs_method_kwargs += ['delta%d_wfn' % x for x in range(2, 6)]

    cbs_methods = []
    for method in cbs_method_kwargs:
        if method in kwargs:
            cbs_methods.append(kwargs[method])
    return cbs_methods


def _parse_cbs_gufunc_string(method_name):
    method_name_list = re.split( """\+(?![^\[\]]*\]|[^\(\)]*\))""", method_name)
    if len(method_name_list) > 2:
        raise ValidationError("CBS gufunc: Text parsing is only valid for a single delta, please use the CBS wrapper directly")

    method_list = []
    basis_list = []
    for num, method_str in enumerate(method_name_list):
        if (method_str.count("[") > 1) or (method_str.count("]") > 1):
            raise ValidationError("""CBS gufunc: Too many brackets given! %s """ % method_str)

        if method_str.count('/') != 1:
            raise ValidationError("""CBS gufunc: All methods must specify a basis with '/'. %s""" % method_str)

        if num > 0:
            method_str = method_str.strip()
            if method_str[:2].lower() != 'd:':
                raise ValidationError("""CBS gufunc: Delta method must start with 'D:'.""")
            else:
                method_str = method_str[2:]
        method, basis = method_str.split('/')
        method_list.append(method)
        basis_list.append(basis)
    return method_list, basis_list


def _cbs_gufunc(func, total_method_name, **kwargs):
    """
    Text based wrapper of the CBS function.
    """

    # Catch kwarg issues
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    core.clean_variables()
    user_dertype = kwargs.pop('dertype', None)
    cbs_verbose = kwargs.pop('cbs_verbose', False)
    ptype = kwargs.pop('ptype', None)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    # Sanitize total_method_name
    label = total_method_name
    total_method_name = total_method_name.lower()
    total_method_name = total_method_name.replace(' ', '')

    # Split into components
    method_list, basis_list = _parse_cbs_gufunc_string(total_method_name)

    # Single energy call?
    single_call = len(method_list) == 1
    single_call &= '[' not in basis_list[0]
    single_call &= ']' not in basis_list[0]

    if single_call:
        method_name = method_list[0]
        basis = basis_list[0]

        # Save some global variables so we can reset them later
        optstash = p4util.OptionsState(['BASIS'])
        core.set_global_option('BASIS', basis)

        ptype_value, wfn = func(method_name, return_wfn=True, molecule=molecule, **kwargs)
        core.clean()

        optstash.restore()

        if return_wfn:
            return (ptype_value, wfn)
        else:
            return ptype_value

    # If we are not a single call, let CBS wrapper handle it!
    cbs_kwargs = {}
    cbs_kwargs['ptype'] = ptype
    cbs_kwargs['return_wfn'] = True
    cbs_kwargs['molecule'] = molecule
    cbs_kwargs['verbose'] = cbs_verbose

    # Find method and basis
    if method_list[0] in ['scf', 'hf']:
        cbs_kwargs['scf_wfn'] = method_list[0]
        cbs_kwargs['scf_basis'] = basis_list[0]
    else:
        cbs_kwargs['corl_wfn'] = method_list[0]
        cbs_kwargs['corl_basis'] = basis_list[0]

    if len(method_list) > 1:
        cbs_kwargs['delta_wfn'] = method_list[1]
        cbs_kwargs['delta_basis'] = basis_list[1]

    ptype_value, wfn = cbs(func, label, **cbs_kwargs)

    if return_wfn:
        return (ptype_value, wfn)
    else:
        return ptype_value
