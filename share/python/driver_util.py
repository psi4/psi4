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

zeta_values = ['d', 't', 'q', '5', '6', '7', '8']
zeta_val2sym = {k+2:v for k, v in zip(range(7), zeta_values)}
zeta_sym2val = {v:k for k, v in zeta_val2sym.items()}

def expand_bracketed_basis(basisstring, molecule=None):
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
        except qcdb.BasisSetNotFound, e:
            raise ValidationError("""Basis set '%s' not available for molecule.""" % (basis))

    return (BSET, ZSET)


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

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Helgaker 2-point SCF extrapolation for method: %s <==\n\n""" % (functionname.upper())
            cbsscheme += """   LO-zeta (%s) Energy:             % 16.14f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%s) Energy:             % 16.14f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Alpha (exponent) Value:          % 16.14f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:        % 16.14f\n\n""" % (beta)

            name_str = "%s/(%s,%s)" % (functionname.upper(), zeta_val2sym[zLO].upper(), zeta_val2sym[zHI].upper())
            cbsscheme += """  @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (18 - len(name_str))
            cbsscheme += """% 16.14f\n\n""" % value
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
            psi4.print_out( """\n   ==> Helgaker 2-point SCF extrapolation for method: %s <==\n\n""" % (functionname.upper()))
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
        alpha = -1 * math.log(ratio)
        beta = (valueHI - valueMD) / (math.exp(-1 * alpha * zMD) * (ratio - 1))
        value = valueHI - beta * math.exp(-1 * alpha * zHI)

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Helgaker 3-point SCF extrapolation for method: %s <==\n\n""" % (functionname.upper())
            cbsscheme += """   LO-zeta (%s) Energy:             % 16.14f\n""" % (str(zLO), valueLO)
            cbsscheme += """   MD-zeta (%s) Energy:             % 16.14f\n""" % (str(zMD), valueMD)
            cbsscheme += """   HI-zeta (%s) Energy:             % 16.14f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Alpha (exponent) Value:          % 16.14f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:        % 16.14f\n\n""" % (beta)

            name_str = "%s/(%s,%s,%s)" % (functionname.upper(), zeta_val2sym[zLO].upper(), zeta_val2sym[zMD].upper(),
                                                             zeta_val2sym[zHI].upper())
            cbsscheme += """  @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (18 - len(name_str))
            cbsscheme += """% 16.14f\n\n""" % value
            psi4.print_out(cbsscheme)

        return value

    elif isinstance(valueLO, (psi4.Matrix, psi4.Vector)):
        valueLO = np.array(valueLO)
        valueMD = np.array(valueMD)
        valueHI = np.array(valueHI)

        nonzero_mask = np.abs(valueHI) > 1.e-14
        top = (valueHI - valueMD)[nonzero_mask]
        bot = (valueMD - valueLO)[nonzero_mask]

        ratio = top/bot
        alpha = -1 * np.log(np.abs(ratio))
        beta = top / (np.exp(-1 * alpha * zMD) * (ratio - 1))
        np_value = valueHI.copy()
        np_value[nonzero_mask] -= beta * np.exp(-1 * alpha * zHI)
        np_value[~nonzero_mask] = 0.0

        # Build and set from numpy routines
        value = psi4.Matrix(*valueHI.shape)
        value_view = np.asarray(value)
        value_view[:] = np_value 
        return value

    else:
        raise ValidationError("scf_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))


def corl_xtpl_helgaker_2(functionname, valueSCF, zLO, valueLO, zHI, valueHI, verbose=True):
    r"""Extrapolation scheme for correlation energies with two adjacent zeta-level bases.
    Used by :py:func:`~wrappers.complete_basis_set`.

    .. math:: E_{corl}^X = E_{corl}^{\infty} + \beta X^{-3}

    """

    if type(valueLO) != type(valueHI):
        raise ValidationError("corl_xtpl_helgaker_2: Inputs must be of the same datatype! (%s, %s)"
                              % (type(valueLO), type(valueHI)))

    if isinstance(valueLO, float):
        value = (valueHI * zHI ** 3 - valueLO * zLO ** 3) / (zHI ** 3 - zLO ** 3)
        beta = (valueHI - valueLO) / (zHI ** (-3) - zLO ** (-3))

        final = valueSCF + value
        if verbose:
            # Output string with extrapolation parameters
            cbsscheme  = """\n\n   ==> Helgaker 2-point correlated extrapolation for method: %s <==\n\n""" % (functionname.upper())
            cbsscheme += """   HI-zeta (%1s) SCF Energy:           % 16.14f\n""" % (str(zHI), valueSCF)
            cbsscheme += """   LO-zeta (%1s) Correlation Energy:   % 16.14f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%1s) Correlation Energy:   % 16.14f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Beta (coefficient) Value:         % 16.14f\n""" % beta
            cbsscheme += """   Extrapolated Correlation Energy:  % 16.14f\n\n""" % value
        
            name_str = "%s/(%s,%s)" % (functionname.upper(), zeta_val2sym[zLO].upper(), zeta_val2sym[zHI].upper())
            cbsscheme += """  @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (19 - len(name_str))
            cbsscheme += """% 16.14f\n\n""" % final
            psi4.print_out(cbsscheme)

        return final

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
            psi4.print_out( """\n   ==> Helgaker 2-point correlated extrapolation for """
                            """method: %s <==\n\n""" % (functionname.upper()))
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
       

        value.add(valueSCF) 
        return value

    else:
        raise ValidationError("scf_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))

