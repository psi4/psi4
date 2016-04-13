#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
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
#@END LICENSE
#

import math
import re
from p4util.exceptions import *
import psi4
import qcdb
import numpy as np

def validate_bracketed_basis(basisstring, molecule=None):
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
    ZETA = ['d', 't', 'q', '5', '6', '7', '8']
    BSET = []
    ZSET = []
    legit_compound_basis = re.compile(r'^(?P<pre>.*cc-.*)\[(?P<zeta>[dtq2345678,]*)\](?P<post>.*z)$', re.IGNORECASE)

    if legit_compound_basis.match(basisstring):
        basisname = legit_compound_basis.match(basisstring)
        # filter out commas and be forgiving of e.g., t5q or 3q
        zetas = [z for z in ZETA if (z in basisname.group('zeta') or str(ZETA.index(z) + 2) in basisname.group('zeta'))]
        for b in zetas:
            if ZSET:
                if (int(ZSET[len(ZSET) - 1]) - ZETA.index(b)) != 1:
                    raise ValidationError("""Basis set '%s' has skipped zeta level '%s'.""" % (basisstring, b))
            BSET.append(basisname.group('pre') + b + basisname.group('post'))
            ZSET.append(ZETA.index(b) + 2)
    elif re.match(r'.*\[.*\].*$', basisstring, flags=re.IGNORECASE):
        raise ValidationError("""Basis series '%s' invalid. Specify a basis series matching"""
                              """ '*cc-*[dtq2345678,]*z'.""" % (basisstring))
    else:
        BSET.append(basisstring)
        ZSET.append(0)

    if molecule is None:
        molecule = """\nH\nH 1 1.00\n"""
    for b in BSET:
        try:
            bsdict = qcdb.BasisSet.pyconstruct(molecule, "BASIS", b)
        except qcdb.BasisSetNotFound, e:
            raise ValidationError("""Basis set '%s' not available for molecule.""" % (b))

    return [BSET, ZSET]


def scf_xtpl_helgaker_2(functionname, zLO, valueLO, zHI, valueHI, verbose=True, alpha=1.63):
    r"""Extrapolation scheme for reference energies with two adjacent zeta-level bases.
    Used by :py:func:`~wrappers.complete_basis_set`.

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
        print(valueLO, valueHI)

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Helgaker 2-point SCF extrapolation for method: %s <==\n\n""" % (functionname)
            cbsscheme += """   LO-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Extrapolated Correlation Energy: %16.8f\n""" % (value)
            cbsscheme += """   Alpha (exponent) Value:          %16.8f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:        %16.8f\n""" % (beta)
            psi4.print_out(cbsscheme)

        return value

    elif isinstance(valueLO, (psi4.Matrix, psi4.Vector)):
        beta = valueHI.clone()
        beta.set_name('Helgaker SCF (%s, %s) beta' % (zLO, zHI))
        beta.subtract(valueLO)
        beta.scale(beta_division)
        beta.scale(beta_mult)

        value = valueHI.clone()
        value.subtract(beta)        
        value.set_name('Helgaker SCF (%s, %s) data' % (zLO, zHI))

        if verbose > 2:
            psi4.print_out( """\n   ==> Helgaker 2-point SCF extrapolation for method: %s <==\n\n""" % (functionname))
            psi4.print_out( """   LO-zeta (%s)""" % str(zLO))
            psi4.print_out( """   LO-zeta Data""")
            valueLO.print_out()
            psi4.print_out( """   HI-zeta (%s)""" % str(zHI))
            psi4.print_out( """   HI-zeta Data""")
            valueHI.print_out()
            psi4.print_out( """   Extrapolated Data:\n""")
            value.print_out()
            psi4.print_out( """   Alpha (exponent) Value:          %16.8f\n""" % (alpha))
            psi4.print_out( """   Beta Data:\n""")
            beta.print_out()
        
        return value

    else:
        raise ValidationError("scf_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))


def scf_xtpl_helgaker_3(functionname, zLO, valueLO, zMD, valueMD, zHI, valueHI, verbose=True):
    r"""Extrapolation scheme for reference energies with three adjacent zeta-level bases.
    Used by :py:func:`~wrappers.complete_basis_set`.

    .. math:: E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha X}
    """

    if (type(valueLO) != type(valueMD)) or (type(valueMD) != type(valueHI)):
        raise ValidationError("scf_xtpl_helgaker_3: Inputs must be of the same datatype! (%s, %s, %s)"
                              % (type(valueLO), type(valueMD), type(valueHI)))

    if isinstance(valueLO, float):

        ratio = (valueHI - valueMD) / (valueMD - valueLO)
        print(valueLO, valueMD, valueHI)
        print(ratio)
        print((valueHI - valueMD))
        print((valueMD - valueLO))
        alpha = -1 * math.log(ratio)
        beta = (valueHI - valueMD) / (math.exp(-1 * alpha * zMD) * (ratio - 1))
        value = valueHI - beta * math.exp(-1 * alpha * zHI)

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Helgaker 3-point SCF extrapolation for method: %s <==\n\n""" % (functionname)
            cbsscheme += """   LO-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zLO), valueLO)
            cbsscheme += """   MD-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zMD), valueMD)
            cbsscheme += """   HI-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Extrapolated Correlation Energy: %16.8f\n""" % (value)
            cbsscheme += """   Alpha (exponent) Value:          %16.8f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:        %16.8f\n""" % (beta)
            psi4.print_out(cbsscheme)

        return value

    elif isinstance(valueLO, (psi4.Matrix, psi4.Vector)):

        ratio = (np.array(valueHI) - np.array(valueMD) / (np.array(valueMD) - np.array(valueLO)))
        np_alpha = -1 * np.log(ratio)
        np_beta = (np.array(valueHI) - valueMD) / (math.exp(-1 * alpha * zMD) * (ratio - 1))
        np_value = np.array(valueHI) - np_beta * np.exp(-1 * np_alpha *zHI)

        # Build and set from numpy routines
        alpha = valueHI.clone()
        beta = valueHI.clone()
        value = valueHI.clone()

        alpha_view = np.asarray(alpha)
        beta_view = np.asarray(beta)
        value_view = np.asarray(value)

        alpha_view[:] = np_alpha
        beta_view[:] = np_beta
        value_view[:] = np_value

        if verbose > 2:
            psi4.print_out( """\n   ==> Helgaker 3-point SCF extrapolation for method: %s <==\n\n""" % (functionname))

            psi4.print_out( """   LO-zeta (%s)""" % str(zLO))
            psi4.print_out( """   LO-zeta Data""")
            valueLO.print_out()

            psi4.print_out( """   MD-zeta (%s)""" % str(zMD))
            psi4.print_out( """   MD-zeta Data""")
            valueMD.print_out()

            psi4.print_out( """   HI-zeta (%s)""" % str(zHI))
            psi4.print_out( """   HI-zeta Data""")
            valueHI.print_out()

            psi4.print_out( """   Extrapolated Data:\n""")
            value.print_out()
            psi4.print_out( """   Alpha (exponent) Value:          %16.8f\n""" % (alpha))
            psi4.print_out( """   Beta Data:\n""")
            beta.print_out()
        
        return value

    else:
        raise ValidationError("scf_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))


def corl_xtpl_helgaker_2(functionname, zLO, valueLO, zHI, valueHI, verbose=True):
    r"""Extrapolation scheme for correlation energies with two adjacent zeta-level bases.
    Used by :py:func:`~wrappers.complete_basis_set`.

    .. math:: E_{corl}^X = E_{corl}^{\infty} + \beta X^{-3}

    """

    if type(valueLO) != type(valueHI):
        raise ValidationError("scf_xtpl_helgaker_2: Inputs must be of the same datatype! (%s, %s)"
                              % (type(valueLO), type(valueHI)))

    if isinstance(valueLO, float):
        value = (valueHI * zHI ** 3 - valueLO * zLO ** 3) / (zHI ** 3 - zLO ** 3)
        beta = (valueHI - valueLO) / (zHI ** (-3) - zLO ** (-3))

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> %s <==\n\n""" % (functionname)
            cbsscheme += """   LO-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Extrapolated Correlation Energy: %16.8f\n""" % (value)
            cbsscheme += """   Beta (coefficient) Value:        %16.8f\n""" % (beta)
            psi4.print_out(cbsscheme)

        return value

    elif isinstance(valueLO, (psi4.Matrix, psi4.Vector)):

        beta = valueHI.clone()
        beta.subtract(valueLO)
        beta.scale(1 / (zHI ** (-3) - zLO ** (-3)))
        beta.set_name('Helgaker SCF (%s, %s) beta' % (zLO, zHI))

        value = valueHI.clone()
        value.scale(zHI ** 3)
        
        tmp = valueLO.clone()
        tmp.scale(zLO ** 3)
        value.subtract(tmp)

        value.scale(1 / (zHI ** 3 - zLO ** 3))
        value.set_name('Helgaker Corr (%s, %s) data' % (zLO, zHI))

        if verbose > 2:
            psi4.print_out( """\n   ==> %s <==\n\n""" % (functionname))
            psi4.print_out( """   LO-zeta (%s)\n""" % str(zLO))
            psi4.print_out( """   LO-zeta Data\n""" % str(zLO))
            valueLO.print_out()
            psi4.print_out( """   HI-zeta (%s)\n""" % str(zHI))
            valueHI.print_out()
            psi4.print_out( """   Extrapolated Data:\n""")
            value.print_out()
            psi4.print_out( """   Alpha (exponent) Value:          %16.8f\n""" % (alpha))
            psi4.print_out( """   Beta Data:\n""")
            beta.print_out()
        
        return value

    else:
        raise ValidationError("scf_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))

