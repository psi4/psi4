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

import math
import re
import sys
from typing import Callable, List, Union, Tuple

import numpy as np

from psi4 import core
from psi4.driver import qcdb
from psi4.driver import p4util
from psi4.driver import driver_util
from psi4.driver import psifiles as psif
from psi4.driver.p4util.exceptions import *
from psi4.driver.procrouting.interface_cfour import cfour_psivar_list
from psi4.driver.procrouting.dft import functionals, build_superfunctional_from_dictionary

zeta_values = 'dtq5678'
_zeta_val2sym = {k + 2: v for k, v in enumerate(zeta_values)}
_zeta_sym2val = {v: k for k, v in _zeta_val2sym.items()}
_addlremark = {'energy': '', 'gradient': ', GRADIENT', 'hessian': ', HESSIAN'}
_lmh_labels = {
    1: ['HI'],
    2: ['LO', 'HI'],
    3: ['LO', 'MD', 'HI'],
    4: ['LO', 'MD', 'M2', 'HI'],
    5: ['LO', 'MD', 'M2', 'M3', 'HI']
}


def _expand_bracketed_basis(
    basisstring: str, 
    molecule: Union[qcdb.molecule.Molecule, core.Molecule] = None
) -> Tuple[List[str], List[int]]:
    """Function to transform and validate basis series specification for cbs().

    Parameters
    ----------
    basisstring
        A string containing the basis sets to be expanded.
        A basis set with no paired square brackets is passed through
        with zeta level 0 (e.g., ``'6-31+G(d,p)'`` is returned as
        ``(["6-31+G(d,p)"], [0])``). A basis set with square brackets is checked
        for sensible sequence and returned as separate basis sets
        (e.g., ``'cc-pV[Q5]Z'` is returned as ``(["cc-pVQZ", "cc-pV5Z"], [4, 5])``).
        Allows out-of-order zeta specification (e.g., ``[qtd]``) and numeral for
        number (e.g., ``[23]``). Does not allow skipped zetas (e.g., ``[dq]``), 
        zetas outside the [2,8] range, non-Dunning, non-Ahlrichs, or non-Jensen 
        sets, or non-findable .gbs sets.
    molecule : qcdb.molecule or psi4.core.Molecule
        This function checks that the basis is valid by trying to build
        the qcdb.BasisSet object for *molecule* or for H2 if None.

    Returns
    -------
    (BSET, ZSET) : Tuple[List[str], List[int]]
        Tuple in the ``([basis set names], [basis set zetas])`` format.

    """
    BSET = []
    ZSET = []
    legit_compound_basis = re.compile(
        r'^(?P<pre>.*cc-.*|def2-|.*pcs+eg-|.*)\[(?P<zeta>[dtq2345678,s1]*)\](?P<post>.*z.*|)$', re.IGNORECASE)
    pc_basis = re.compile(r'.*pcs+eg-$', re.IGNORECASE)
    def2_basis = re.compile(r'def2-', re.IGNORECASE)
    zapa_basis = re.compile(r'.*zapa.*', re.IGNORECASE)

    if legit_compound_basis.match(basisstring):
        basisname = legit_compound_basis.match(basisstring)
        # handle def2-svp* basis sets as double-zeta
        if def2_basis.match(basisname.group('pre')):
            bn_gz = basisname.group('zeta').replace("s", "d")
        # handle pc-n basis set polarisation -> zeta conversion
        elif pc_basis.match(basisname.group('pre')):
            bn_gz = basisname.group('zeta').replace("4", "5").replace("3", "4").replace("2", "3").replace("1", "2")
        else:
            bn_gz = basisname.group('zeta')
        # filter out commas and be forgiving of e.g., t5q or 3q
        zetas = [z for z in zeta_values if (z in bn_gz or str(zeta_values.index(z) + 2) in bn_gz)]
        for b in zetas:
            if ZSET and (int(ZSET[len(ZSET) - 1]) - zeta_values.index(b)) != 1:
                raise ValidationError("""Basis set '%s' has skipped zeta level '%s'.""" %
                                      (basisstring, _zeta_val2sym[_zeta_sym2val[b] - 1]))
            # reassemble def2-svp* properly instead of def2-dzvp*
            if def2_basis.match(basisname.group('pre')) and b == "d":
                BSET.append(basisname.group('pre') + "s" + basisname.group('post')[1:])
            # reassemble pc-n basis sets properly
            elif pc_basis.match(basisname.group('pre')):
                BSET.append(basisname.group('pre') + "{0:d}".format(_zeta_sym2val[b] - 1))
            # assemble nZaPa basis sets
            elif zapa_basis.match(basisname.group('post')):
                bzapa = b.replace("d", "2").replace("t", "3").replace("q", "4")
                BSET.append(basisname.group('pre') + bzapa + basisname.group('post'))
            else:
                BSET.append(basisname.group('pre') + b + basisname.group('post'))
            ZSET.append(zeta_values.index(b) + 2)
    elif re.match(r'.*\[.*\].*$', basisstring, flags=re.IGNORECASE):
        raise ValidationError(
            """Basis series '%s' invalid. Specify a basis series matching"""
            """ '*cc-*[dtq2345678,]*z*'. or 'def2-[sdtq]zvp*' or '*pcs[s]eg-[1234]' or '[1234567]ZaPa' """ %
            (basisstring))
    else:
        BSET.append(basisstring)
        ZSET.append(0)

    if molecule is None:
        molecule = """\nH\nH 1 1.00\n"""
    elif isinstance(molecule, core.Molecule):
        molecule = qcdb.Molecule(molecule.to_dict())

    for basis in BSET:
        try:
            qcdb.BasisSet.pyconstruct(molecule, "BASIS", basis)
        except qcdb.BasisSetNotFound:
            e = sys.exc_info()[1]
            raise ValidationError(f"""Basis set '{basis}' not available for molecule.""")

    return (BSET, ZSET)


def _contract_bracketed_basis(basisarray: List[str]) -> str:
    """Function to reform a bracketed basis set string from a sequential series
    of basis sets. Essentially the inverse of _expand_bracketed_basis(). Used to
    print a nicely formatted basis set string in the results table.

    Parameters
    ----------
    basisarray
        Basis set names, differing by zeta level, e.g. ``["cc-pvqz", "cc-pv5z"]``.

    Returns
    -------
    basisstring : str
        A nicely formatted basis set string, e.g. ``"cc-pv[q5]z"`` for the above example.

    """

    if len(basisarray) == 1:
        return basisarray[0]

    else:
        zetaindx = [i for i in range(len(basisarray[0])) if basisarray[0][i] != basisarray[1][i]][0]
        ZSET = [bas[zetaindx] for bas in basisarray]

        pre = basisarray[1][:zetaindx]
        post = basisarray[1][zetaindx + 1:]
        basisstring = pre + '[' + ''.join(ZSET) + ']' + post
        return basisstring


### GENERIC FUNCTIONS
def xtpl_highest_1(
    functionname: str,
    zHI: int,
    valueHI: Union[float, core.Matrix, core.Vector],
    verbose: bool = True,
    **kwargs
) -> Union[float, core.Matrix, core.Vector]:
    r"""Scheme for energies with a single basis or the highest zeta-level among 
    an array of bases. Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname
        Name of the CBS component.
    zHI
        Zeta-level, only used for printing.
    valueHI
        Value of the CBS component. Type :class:`float` for energy, 
        :class:`core.Matrix` or :class:`core.Vector` for gradients.

    Returns
    -------
    valueHI : Union[float, core.Matrix, core.Vector]
        Returns :math:`E_{total}^{\infty}` which is equal to valueHI.

    Notes
    -----
    .. math:: E_{total}^X = E_{total}^{\infty}

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
            core.print_out("""   HI-zeta (%s) Total Data:\n""" % (str(zHI)))
            valueHI.print_out()

        return valueHI


def xtpl_exponential_2(
    functionname: str,  
    zLO: int, 
    valueLO: Union[float, core.Matrix, core.Vector],
    zHI: int,
    valueHI: Union[float, core.Matrix, core.Vector],
    alpha: float,
    verbose: bool = True
) -> Union[float, core.Matrix, core.Vector]:
    r"""Extrapolation scheme using exponential form for energies with two 
    adjacent zeta-level bases. Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname
        Name of the CBS component.
    zLO
        Lower zeta level.
    valueLO
        Lower value used for extrapolation. Type :class:`float` for energy, 
        :class:`core.Matrix` or :class:`core.Vector` for gradients.
    zHI
        Higher zeta level. Should be equal to zLO + 1.
    valueHI
        Higher value used for extrapolation. Type :class:`float` for energy, 
        :class:`core.Matrix` or :class:`core.Vector` for gradients.
    alpha
        Extrapolation parameter :math:`\alpha`.

    Returns
    -------
    value : Union[float, core.Matrix, core.Vector]
        Returns :math:`E_{total}^{\infty}`, or the analogue for derivatives,
        see below.

    Notes
    -----
    The extrapolation is calculated according to:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha X}`
    
    """
    if type(valueLO) != type(valueHI):
        raise ValidationError("xtpl_exponential_2: Inputs must be of the same datatype! (%s, %s)" %
                              (type(valueLO), type(valueHI)))

    if not isinstance(alpha, float):
        raise ValidationError(f"xtpl_exponential_2: Provided alpha must be a float, not {alpha}.")
    
    beta_division = 1 / (math.exp(-1 * alpha * zLO) * (math.exp(-1 * alpha) - 1))
    beta_mult = math.exp(-1 * alpha * zHI)

    if isinstance(valueLO, float):
        beta = (valueHI - valueLO) * beta_division
        value = valueHI - beta * beta_mult

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> 2-point exponential extrapolation for method: %s <==\n\n""" % (
                functionname.upper())
            cbsscheme += """   LO-zeta (%s) Energy:               % 16.12f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Alpha (exponent) Value:           % 16.12f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n\n""" % (beta)

            name_str = "%s/(%s,%s)" % (functionname.upper(), _zeta_val2sym[zLO].upper(), _zeta_val2sym[zHI].upper())
            cbsscheme += """   @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (18 - len(name_str))
            cbsscheme += """% 16.12f\n\n""" % value
            core.print_out(cbsscheme)

        return value

    elif isinstance(valueLO, (core.Matrix, core.Vector)):
        beta = valueHI.clone()
        beta.name = 'Exponential extrapolation (%s, %s) beta' % (zLO, zHI)
        beta.subtract(valueLO)
        beta.scale(beta_division)
        beta.scale(beta_mult)

        value = valueHI.clone()
        value.subtract(beta)
        value.name = 'Exponential extrapolation (%s, %s) data' % (zLO, zHI)

        if verbose > 2:
            core.print_out("""\n   ==> 2-point exponential extrapolation for method: %s <==\n\n""" %
                           (functionname.upper()))
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
        raise ValidationError("xtpl_exponential_2: datatype is not recognized '%s'." % type(valueLO))


def xtpl_power_2(
    functionname: str,
    zLO: int,
    valueLO: Union[float, core.Matrix, core.Vector],
    zHI: int,
    valueHI: Union[float, core.Matrix, core.Vector],
    alpha: float,
    verbose: bool = True
) -> Union[float, core.Matrix, core.Vector]:
    r"""Extrapolation scheme using power form for energies with two adjacent 
    zeta-level bases. Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname
        Name of the CBS component.
    zLO
        Lower zeta level.
    valueLO
        Lower value used for extrapolation. Type :class:`float` for energy, 
        :class:`core.Matrix` or :class:`core.Vector` for gradients.
    zHI
        Higher zeta level. Should be equal to zLO + 1.
    valueHI
        Higher value used for extrapolation. Type :class:`float` for energy, 
        :class:`core.Matrix` or :class:`core.Vector` for gradients.
    alpha
        Extrapolation parameter :math:`\alpha`.

    Returns
    -------
    value : Union[float, core.Matrix, core.Vector]
        Returns :math:`E_{total}^{\infty}`, or the analogue for derivatives,
        see below.

    Notes
    -----
    The extrapolation is calculated according to:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta X^{-\alpha}`

    """
    if type(valueLO) != type(valueHI):
        raise ValidationError("xtpl_power_2: Inputs must be of the same datatype! (%s, %s)" %
                              (type(valueLO), type(valueHI)))

    if alpha is None:
        raise ValidationError("xtpl_power_2: alpha must be provided.")

    beta_division = 1 / (zHI**(-1 * alpha) - zLO**(-1 * alpha))
    beta_mult = zHI**(-1 * alpha)

    if isinstance(valueLO, float):
        beta = (valueHI - valueLO) * beta_division
        value = valueHI - beta * beta_mult

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> 2-point power form extrapolation for method: %s <==\n\n""" % (
                functionname.upper())
            cbsscheme += """   LO-zeta (%s) Energy:               % 16.12f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Alpha (exponent) Value:           % 16.12f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n\n""" % (beta)

            name_str = "%s/(%s,%s)" % (functionname.upper(), _zeta_val2sym[zLO].upper(), _zeta_val2sym[zHI].upper())
            cbsscheme += """   @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (18 - len(name_str))
            cbsscheme += """% 16.12f\n\n""" % value
            core.print_out(cbsscheme)

        return value

    elif isinstance(valueLO, (core.Matrix, core.Vector)):
        beta = valueHI.clone()
        beta.name = 'Power extrapolation (%s, %s) beta' % (zLO, zHI)
        beta.subtract(valueLO)
        beta.scale(beta_division)
        beta.scale(beta_mult)

        value = valueHI.clone()
        value.subtract(beta)
        value.name = 'Power extrapolation (%s, %s) data' % (zLO, zHI)

        if verbose > 2:
            core.print_out("""\n   ==> 2-point power from SCF extrapolation for method: %s <==\n\n""" %
                           (functionname.upper()))
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
        raise ValidationError("xtpl_power_2: datatype is not recognized '%s'." % type(valueLO))


def xtpl_expsqrt_2(
    functionname: str,
    zLO: int,
    valueLO: Union[float, core.Matrix, core.Vector],
    zHI: int,
    valueHI: Union[float, core.Matrix, core.Vector],
    alpha: float = None,
    verbose: bool = True
) -> Union[float, core.Matrix, core.Vector]:
    r"""Extrapolation scheme using square root-exponential form for energies 
    with two adjacent zeta-level bases. Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname
        Name of the CBS component.
    zLO
        Lower zeta level.
    valueLO
        Lower value used for extrapolation. Type :class:`float` for energy, 
        :class:`core.Matrix` or :class:`core.Vector` for gradients.
    zHI
        Higher zeta level. Should be equal to zLO + 1.
    valueHI
        Higher value used for extrapolation. Type :class:`float` for energy, 
        :class:`core.Matrix` or :class:`core.Vector` for gradients.
    alpha
        Extrapolation parameter :math:`\alpha`.

    Returns
    -------
    value : Union[float, core.Matrix, core.Vector]
        Returns :math:`E_{total}^{\infty}`, or the analogue for derivatives,
        see below.

    Notes
    -----
    The extrapolation is calculated according to:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha\sqrt{X}}`

    """

    if type(valueLO) != type(valueHI):
        raise ValidationError("xtpl_expsqrt_2: Inputs must be of the same datatype! (%s, %s)" %
                              (type(valueLO), type(valueHI)))

    if alpha is None:
        raise ValidationError("xtpl_expsqrt_2: alpha must be provided.")

    beta_division = 1 / (math.exp(-1 * alpha) * (math.exp(math.sqrt(zHI)) - math.exp(math.sqrt(zLO))))
    beta_mult = math.exp(-1 * alpha * math.sqrt(zHI))

    if isinstance(valueLO, float):
        beta = (valueHI - valueLO) * beta_division
        value = valueHI - beta * beta_mult

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> 2-point exp.-sqrt. form SCF extrapolation for method: %s <==\n\n""" % (
                functionname.upper())
            cbsscheme += """   LO-zeta (%s) Energy:               % 16.12f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Alpha (exponent) Value:           % 16.12f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n\n""" % (beta)

            name_str = "%s/(%s,%s)" % (functionname.upper(), _zeta_val2sym[zLO].upper(), _zeta_val2sym[zHI].upper())
            cbsscheme += """   @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (18 - len(name_str))
            cbsscheme += """% 16.12f\n\n""" % value
            core.print_out(cbsscheme)

        return value

    elif isinstance(valueLO, (core.Matrix, core.Vector)):
        beta = valueHI.clone()
        beta.name = 'Expsqrt extrapolation (%s, %s) beta' % (zLO, zHI)
        beta.subtract(valueLO)
        beta.scale(beta_division)
        beta.scale(beta_mult)

        value = valueHI.clone()
        value.subtract(beta)
        value.name = 'Expsqrt extrapolation (%s, %s) data' % (zLO, zHI)

        if verbose > 2:
            core.print_out("""\n   ==> 2-point exp.-sqrt. from SCF extrapolation for method: %s <==\n\n""" %
                           (functionname.upper()))
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
        raise ValidationError("xtpl_expsqrt_2: datatype is not recognized '%s'." % type(valueLO))


def xtpl_exponential_3(
    functionname: str,
    zLO: int,
    valueLO: Union[float, core.Matrix, core.Vector],
    zMD: int,
    valueMD: Union[float, core.Matrix, core.Vector],
    zHI: int,
    valueHI: Union[float, core.Matrix, core.Vector],
    verbose: bool = True,
    **kwargs
) -> Union[float, core.Matrix, core.Vector]:
    r"""Extrapolation scheme using exponential form for energies with three 
    adjacent zeta-level bases. Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname
        Name of the CBS component.
    zLO
        Lower zeta level.
    valueLO
        Lower value used for extrapolation. Type :class:`float` for energy, 
        :class:`core.Matrix` or :class:`core.Vector` for gradients.
    zMD
        Intermediate zeta level. Should be equal to zLO + 1.
    valueMD
        Intermediate value used for extrapolation. Type :class:`float` for energy, 
        :class:`core.Matrix` or :class:`core.Vector` for gradients.
    zHI
        Higher zeta level. Should be equal to zLO + 2.
    valueHI
        Higher value used for extrapolation.

    Returns
    -------
    value : Union[float, core.Matrix, core.Vector]
        Returns :math:`E_{total}^{\infty}`, or the analogue for derivatives,
        see below.

    Notes
    -----
    The extrapolation is calculated according to [4]_:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha X}`

    """

    if (type(valueLO) != type(valueMD)) or (type(valueMD) != type(valueHI)):
        raise ValidationError("xtpl_exponential_3: Inputs must be of the same datatype! (%s, %s, %s)" %
                              (type(valueLO), type(valueMD), type(valueHI)))

    if isinstance(valueLO, float):

        ratio = (valueHI - valueMD) / (valueMD - valueLO)
        alpha = -1 * math.log(ratio)
        beta = (valueHI - valueMD) / (math.exp(-1 * alpha * zMD) * (ratio - 1))
        value = valueHI - beta * math.exp(-1 * alpha * zHI)

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> 3-point exponential extrapolation for method: %s <==\n\n""" % (
                functionname.upper())
            cbsscheme += """   LO-zeta (%s) Energy:               % 16.12f\n""" % (str(zLO), valueLO)
            cbsscheme += """   MD-zeta (%s) Energy:               % 16.12f\n""" % (str(zMD), valueMD)
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Alpha (exponent) Value:           % 16.12f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n\n""" % (beta)

            name_str = "%s/(%s,%s,%s)" % (functionname.upper(), _zeta_val2sym[zLO].upper(), _zeta_val2sym[zMD].upper(),
                                          _zeta_val2sym[zHI].upper())
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
        raise ValidationError("xtpl_exponential_3: datatype is not recognized '%s'." % type(valueLO))


### NAMED FUNCTION ALIASES
def scf_xtpl_helgaker_2(
    functionname: str,
    zLO: int,
    valueLO: Union[float, core.Matrix, core.Vector],
    zHI: int,
    valueHI: Union[float, core.Matrix, core.Vector],
    alpha: float = None,
    verbose: bool = True
) -> Union[float, core.Matrix, core.Vector]:
    r"""Alias for the exponential form of the extrapolation scheme for reference 
    energies with two adjacent zeta-level bases, with :math:`\alpha = 1.63`, see 
    [1]_. Based on the average of 24 :math:`\alpha`s obtained by fitting 2-6, 
    2-5, and 3-6 zeta results with the exponential form for H2, C2, N2, BH, HF, 
    BF, CO, and NO+, see [1]_. This is the standard extrapolation function for 
    reference energies in Psi4. Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname
        Name of the CBS component.
    zLO
        Lower zeta level.
    valueLO
        Lower value used for extrapolation.
    zHI
        Higher zeta level. Should be equal to zLO + 1.
    valueHI
        Higher value used for extrapolation.
    alpha
        Overrides the default :math:`\alpha = 1.63`

    Returns
    -------
    value : Union[float, core.Matrix, core.Vector]
        Returns :math:`E_{total}^{\infty}`, or the analogue for derivatives,
        see below.

    Notes
    -----
    The extrapolation is calculated according to [1]_:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha X}, \alpha = 1.63`

    References
    ----------

    .. [1] Halkier, Helgaker, Jorgensen, Klopper, & Olsen, 
       Chem. Phys. Lett. 302 (1999) 437-446,
       DOI: 10.1016/S0009-2614(99)00179-7

    """

    if alpha is None:
        alpha = 1.63

    return xtpl_exponential_2(functionname, zLO, valueLO, zHI, valueHI, alpha, verbose=verbose)


def scf_xtpl_truhlar_2(
    functionname: str,
    zLO: int,
    valueLO: Union[float, core.Matrix, core.Vector],
    zHI: int,
    valueHI: Union[float, core.Matrix, core.Vector],
    alpha: float = None,
    verbose: bool = True
) -> Union[float, core.Matrix, core.Vector]:
    r"""Alias for the power form of the extrapolation scheme for reference 
    energies with two adjacent zeta-level bases with :math:`\alpha = 3.4`, see 
    [2]_. Derived by minimizing the RMS error for Ne, HF, H2O; particularly 
    useful for extrapolation from [2,3]-zeta Dunning-style basis sets, see [2]_. 
    Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname : string
        Name of the CBS component.
    zLO
        Lower zeta level.
    valueLO
        Lower value used for extrapolation.
    zHI
        Higher zeta level. Should be equal to zLO + 1.
    valueHI
        Higher value used for extrapolation.
    alpha
        Overrides the default :math:`\alpha = 3.4`

    Returns
    -------
    value : Union[float, core.Matrix, core.Vector]
        Returns :math:`E_{total}^{\infty}`, or the analogue for derivatives,
        see below.

    Notes
    -----
    The extrapolation is calculated according to [2]_:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta X^{-\alpha}, \alpha = 3.4`

    References
    ----------

    .. [2] Truhlar, Chem. Phys. Lett. 294 (1998) 45-48,
       DOI: 10.1016/S0009-2614(98)00866-5

    """

    if alpha is None:
        alpha = 3.40
    
    return xtpl_power_2(functionname, zLO, valueLO, zHI, valueHI, alpha, verbose=verbose)


def scf_xtpl_karton_2(
    functionname: str,
    zLO: int,
    valueLO: Union[float, core.Matrix, core.Vector],
    zHI: int,
    valueHI: Union[float, core.Matrix, core.Vector],
    alpha: float = None,
    verbose: bool = True
) -> Union[float, core.Matrix, core.Vector]:
    r"""Alias for the exponential-square root of the extrapolation scheme for reference energies with two adjacent
    zeta-level bases, see e.g. [3]_. Origin of :math:`\alpha = 6.3` unknown. The exponential-square root function is 
    known to be the correct scaling for HF and DFT, see [4]_; this is due to the nuclear cusp, see [5]_ and [6]_. 
    Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname
        Name of the CBS component.
    zLO
        Lower zeta level.
    valueLO
        Lower value used for extrapolation.
    zHI
        Higher zeta level. Should be equal to zLO + 1.
    valueHI
        Higher value used for extrapolation.
    alpha
        Overrides the default :math:`\alpha = 6.3`

    Returns
    -------
    value : Union[float, core.Matrix, core.Vector]
        Returns :math:`E_{total}^{\infty}`, or the analogue for derivatives,
        see below.

    Notes
    -----
    The extrapolation is calculated according to [3]_:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha\sqrt{X}}, \alpha = 6.3`

    References
    ----------

    .. [3] Karton, Martin, Theor. Chem. Acc. 115 (2006) 330-333,
       DOI: 10.1007/s00214-005-0028-6
    
    .. [4] Shaw, Int. J. Quantum. Chem. 120 (2020) e26264,
       DOI: 10.1002/qua.26264

    .. [5] McKemmish, Gill, J. Chem. Theory Comput. 8 (2012) 267-273,
       DOI: 10.1021/ct300559t
    
    .. [6] Klopper, Kutzelnigg, J. Mol. Struct. 135 (1986) 339-356,
       DOI: 10.1016/0166-1280(86)80068-9

    """

    if alpha is None:
        alpha = 6.3

    return xtpl_expsqrt_2(functionname, zLO, valueLO, zHI, valueHI, alpha, verbose=verbose)


def scf_xtpl_helgaker_3(
    functionname: str,
    zLO: int,
    valueLO: Union[float, core.Matrix, core.Vector],
    zMD: int,
    valueMD: Union[float, core.Matrix, core.Vector],
    zHI: int,
    valueHI: Union[float, core.Matrix, core.Vector],
    verbose: bool = True,
    **kwargs
) -> Union[float, core.Matrix, core.Vector]:
    r"""Extrapolation scheme for reference energies with three adjacent 
    zeta-level bases. Here, :math:`\alpha` is calculated directly from the three 
    datapoints in addition to  the pre-factor :math:`\beta`, see [7]_. Unreliable, 
    as the smallest basis set is too small to provide much benefit to extrapolate 
    the largest basis set; even [7]_ suggests to use a two-point extrapolation 
    with :math:`\alpha = 1.63`, i.e. :py:func:`scf_xtpl_helgaker_2`, instead. 
    This  is the default extrapolation function for three point extrapolations 
    in Psi4. Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname
        Name of the CBS component.
    zLO
        Lower zeta level.
    valueLO
        Lower value used for extrapolation.
    zMD
        Intermediate zeta level. Should be equal to zLO + 1.
    valueMD
        Intermediate value used for extrapolation.
    zHI
        Higher zeta level. Should be equal to zLO + 2.
    valueHI
        Higher value used for extrapolation.
    
    Returns
    -------
    value : Union[float, core.Matrix, core.Vector]
        Returns :math:`E_{total}^{\infty}`, or the analogue for derivatives,
        see below.

    Notes
    -----
    The extrapolation is calculated according to [4]_:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha X}`
    with :math:`\alpha` fitted to data.

    References
    ----------

    .. [7] Halkier, Helgaker, Jorgensen, Klopper, & Olsen, Chem. Phys. Lett. 302 (1999) 437-446,
       DOI: 10.1016/S0009-2614(99)00179-7

    """

    return xtpl_exponential_3(functionname, zLO, valueLO, zMD, valueMD, zHI, valueHI, verbose=verbose)


def corl_xtpl_helgaker_2(
    functionname: str,
    zLO: int,
    valueLO: Union[float, core.Matrix, core.Vector],
    zHI: int,
    valueHI: Union[float, core.Matrix, core.Vector],
    alpha: float = None,
    verbose: bool = True
) -> Union[float, core.Matrix, core.Vector]:
    r"""Alias for the cubic power extrapolation scheme for correlation energies 
    with two adjacent zeta-level bases. The :math:`\alpha = 3.0` is from 
    correlated calculations of H2O, see [8]_; it has been confirmed to work well
    for HF, Ne, H2O with Dunning basis sets, see [9]_. This is the default 
    extrapolation method for correlated calculations in Psi4. Used by 
    :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname
        Name of the CBS component.
    zLO
        Lower zeta level.
    valueLO
        Lower value used for extrapolation.
    zHI
        Higher zeta level. Should be equal to zLO + 1.
    valueHI
        Higher value used for extrapolation.
    alpha
        Overrides the default :math:`\alpha = 3.0`

    Returns
    -------
    value : Union[float, core.Matrix, core.Vector]
        Returns :math:`E_{total}^{\infty}`, or the analogue for derivatives,
        see below.

    Notes
    -----
    The extrapolation is calculated according to [5]_:
    :math:`E_{corl}^X = E_{corl}^{\infty} + \beta X^{-\alpha}`

    References
    ----------
    
    .. [8] Helgaker, Klopper, Koch, Noga, J. Chem. Phys. 106 (1997) 9639-9646,
       DOI: 10.1063/1.473863

    .. [9] Halkier, Helgaker, Jorgensen, Klopper, Koch, Olsen, & Wilson,
       Chem. Phys. Lett. 286 (1998) 243-252,
       DOI: 10.1016/S0009-2614(99)00179-7

    """

    if alpha is None:
        alpha = 3.0

    return xtpl_power_2(functionname, zLO, valueLO, zHI, valueHI, alpha, verbose=verbose)


def return_energy_components() -> dict:
    """
    Define some quantum chemical knowledge, namely what methods are subsumed in others.
    The keys of the returned dict are callable methods, the values contain a set of
    component: Psi4 variable pairs. In the wider context of :func:`~psi4.cbs`,
    ``"wfn"`` refers to the top-level keys, ``"component"`` refers to the nested keys.
    """

    # yapf: disable
    VARH = {}
    VARH['scf'] = {
                            'scf': 'SCF TOTAL ENERGY'}
    VARH['hf'] = {
                             'hf': 'HF TOTAL ENERGY'}
    VARH['mp2'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY'}
    VARH['dlpno-mp2'] = {
                             'hf': 'HF TOTAL ENERGY',
                      'dlpno-mp2': 'MP2 TOTAL ENERGY'}
    VARH['mp2d'] = {
                             'hf': 'HF TOTAL ENERGY',
                           'mp2d': 'MP2D TOTAL ENERGY'}
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
                            'mp4': 'MP4 TOTAL ENERGY'}
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
                       'fno-ccsd': 'CCSD TOTAL ENERGY',
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
    VARH['ccsd(at)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'ccsd': 'CCSD TOTAL ENERGY',
                       'ccsd(at)': 'CCSD(AT) TOTAL ENERGY'}
    VARH["a-ccsd(t)"] = VARH["ccsd(at)"]
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
        VARH[f'ci{cilevel}'] = {
                             'hf': 'HF TOTAL ENERGY',
                   f'ci{cilevel}': 'CI TOTAL ENERGY'}

    for mplevel in range(5, 99):
        VARH[f'mp{mplevel}'] = {
                             'hf': 'HF TOTAL ENERGY',
                   f'mp{mplevel}': f'MP{mplevel} TOTAL ENERGY'}
        for mplevel2 in range(2, mplevel):
            VARH[f'mp{mplevel}'][f'mp{mplevel2}'] = f'MP{mplevel2} TOTAL ENERGY'

    # Incorporate CFOUR methods
    VARH.update(cfour_psivar_list())

    # Incorporate density functionals
    # Somewhat confusingly, we need both {funcname} and {funcname}-dft, even 
    # though they point at the same energy:
    #   - {funcname} matches the name that is extracted as a wfn, either by
    #     specifications in the metadata directly, or from a method/basis format
    #   - {funcname}-dft is the component that corresponds to the total energy
    #     of the functional recipe, useful for "naive" extrapolations
    # The other components -fctl, -disp, -dh, and -nl are self-explanatory,
    # however note that currently -fctl includes -nl but does not include -disp
    for funcname in functionals:
        if funcname == "hf" or funcname == "scf":
            continue
        VARH[funcname] = {
                         funcname: 'DFT TOTAL ENERGY',
                f'{funcname}-dft': 'DFT TOTAL ENERGY',
               f'{funcname}-fctl': 'DFT FUNCTIONAL TOTAL ENERGY'}
        if 'dispersion' in functionals[funcname].keys():
            if functionals[funcname]['dispersion']['type'] == 'nl':
                VARH[funcname][f'{funcname}-nl'] = 'DFT VV10 ENERGY'
            else:
                VARH[funcname][f'{funcname}-disp'] = 'DISPERSION CORRECTION ENERGY'
        if 'c_mp2' in functionals[funcname].keys():
            VARH[funcname][f'{funcname}-dh'] = 'DOUBLE-HYBRID CORRECTION ENERGY'
    return VARH
    # yapf: enable


VARH = return_energy_components()


def _get_default_xtpl(nbasis: int, xtpl_type: Union[str, None]) -> Callable:
    """ A helper function to determine default extrapolation type.

    Parameters
    ----------
    nbasis
        Number of basis sets
    xtpl_type 
        {'scf', 'corl', 'fctl', 'dh', None}
        Extrapolation type: 'scf' and 'fctl' for the total energy, 'corl' and 'dh'
        for just the correlation component, None if no extrapolation requested.

    Returns
    -------
    Callable
        Extrapolation function to be used.
    """

    if nbasis == 1:
        allowed = ["scf", "corl", "fctl", "dh", "dft", None]
        if xtpl_type in allowed:
            return xtpl_highest_1
        else:
            raise ValidationError(f"Stage treatment must be one of {allowed}, not '{xtpl_type}'")
    elif nbasis == 2:
        allowed = ["scf", "corl", "fctl", "dh", "dft"]
        if xtpl_type == "scf":
            return scf_xtpl_helgaker_2
        elif xtpl_type == "corl":
            return corl_xtpl_helgaker_2
        elif xtpl_type == "fctl":
            return xtpl_expsqrt_2
        elif xtpl_type == "dft":
            return xtpl_expsqrt_2
        elif xtpl_type == "dh":
            return xtpl_power_2
        else:
            raise ValidationError(f"Stage treatment must be one of {allowed}, not '{xtpl_type}'")
    elif nbasis == 3:
        allowed = ["scf"]
        if xtpl_type == "scf":
            return scf_xtpl_helgaker_3
        else:
            raise ValidationError(f"Stage treatment must be one of {allowed}, not '{xtpl_type}'")
    else:
        raise ValidationError(f"Wrong number of basis sets supplied to scf_xtpl: {nbasis}")


def _get_dfa_alpha(
    xtpl_type: str, 
    bdata: Tuple[List[str], List[int]], 
    funcname: str
) -> float:
    """ A helper function to determine default extrapolation alpha for DFT.

    The parameters for the 'fctl' component are from Kraus, [1]_ based on a fit
    using numerical results of 30 singlet diatomics with PBE-like functionals
    with varying amount of HF exchange and MP2 correlation.

    The 'dh' values are based on the standard correlation alpha of 3.0 for [3,4]
     extrapolation, [2]_ and Truhlar's MP2 alpha of 2.2 for [2,3] extrapolation. [3]_
     
    Parameters
    ----------
    xtpl_type
        {'fctl', 'dh'}
        Extrapolation type: 'fctl' for the functional energy,
        'dh' for the double-hybrid correlation component.
    bdata
        Names and zetas of the basis sets to extrapolate from.
    funcname
        Name of the functional.

    Returns
    -------
    alpha : float
        Extrapolation alpha to be used.

    References
    ----------
    .. [1] Kraus, J. Chem. Theor. Comput. 17 (2021) 5651-5660,
       DOI: 10.1021/acs.jctc.1c00542
    .. [2] Halkier, Helgaker, Jorgensen, Klopper, Koch, Olsen, & Wilson,
       Chem. Phys. Lett. 286 (1998) 243-252,
       DOI: 10.1016/S0009-2614(99)00179-7
    .. [3] Truhlar, Chem. Phys. Lett. 294 (1998) 45-48,
       DOI: 10.1016/S0009-2614(98)00866-5
    """
    bnames, bzetas = bdata
    if funcname not in functionals:
        raise ValidationError(f"Functional name {funcname} is undefined.")
    elif xtpl_type in ["fctl", "dft"]:
        if len({f"cc-pv{x}z" for x in "dtq56"}.intersection(bnames)) >= 2:
            alphadata = [3.622, 1.511, 0.005]
        elif len({f"cc-pwcv{x}z" for x in "dtq5"}.intersection(bnames)) >= 2:
            alphadata = [4.157, 1.192, -0.048]
        elif len({f"aug-cc-pv{x}z" for x in "dtq56"}.intersection(bnames)) >= 2:
            alphadata = [3.676, 1.887, 0.139]
        elif len({f"aug-cc-pwcv{x}z" for x in "dtq5"}.intersection(bnames)) >= 2:
            alphadata = [4.485, 1.445, 0.085]
        elif len({"def2-svp", "def2-tzvp", "def2-qzvp"}.intersection(bnames)) >= 2:
            alphadata = [7.406, 1.266, -0.046]
        elif len({"def2-svp", "def2-tzvpp", "def2-qzvpp"}.intersection(bnames)) >= 2:
            alphadata = [7.408, 1.267, -0.046]
        elif len({"def2-svpd", "def2-tzvpd", "def2-qzvpd"}.intersection(bnames)) >= 2:
            alphadata = [7.925, 1.370, 0.101]
        elif len({"def2-svpd", "def2-tzvppd", "def2-qzvppd"}.intersection(bnames)) >= 2:
            alphadata = [7.927, 1.371, 0.101]
        elif len({f"pc-{N}" for N in "01234"}.intersection(bnames)) >= 2:
            alphadata = [6.172, -1.623, 0.183]
        elif len({f"pcseg-{N}" for N in "01234"}.intersection(bnames)) >= 2:
            alphadata = [5.883, -1.825, 0.227]
        elif len({f"aug-pc-{N}" for N in "01234"}.intersection(bnames)) >= 2:
            alphadata = [6.390, -1.874, 0.260]
        elif len({f"aug-pcseg-{N}" for N in "01234"}.intersection(bnames)) >= 2:
            alphadata = [6.166, -2.137, 0.296]
        elif len({f"jorge-{x}zp" for x in "dtq56"}.intersection(bnames)) >= 2:
            alphadata = [3.531, 0.338, -0.011]
        elif len({f"jorge-a{x}zp" for x in "dtq5"}.intersection(bnames)) >= 2:
            alphadata = [3.386, 0.245, 0.033]
        elif len({f"{x}zapa-nr" for x in "23456"}.intersection(bnames)) >= 2:
            alphadata = [3.306, 2.525, -0.040]
        elif len({f"{x}zapa-nr-cv" for x in "23456"}.intersection(bnames)) >= 2:
            alphadata = [5.618, 0.490, -0.016]
        else:
            raise ValidationError(f"DFT fctl extrapolation alpha for basis sets {bnames} undefined.")
        sup = build_superfunctional_from_dictionary(functionals[funcname], 1, 1, True)[0]
        alpha = alphadata[0] + alphadata[1] * sup.x_alpha() + alphadata[2] * sup.c_alpha()
        return alpha
    elif xtpl_type == "dh":
        if bzetas == [2, 3]:
            return 2.2
        elif bzetas == [3, 4]:
            return 3.0
        else:
            raise ValidationError(f"DFT dh extrapolation alpha for basis sets {bnames} undefined.")


def _get_default_alpha(
    xtpl_type: Union[str, None], 
    bdata: Tuple[List[str], List[int]], 
    funcname: str = None
) -> Union[float, None]:
    """ 
    A helper function to determine default extrapolation alpha. Currently only
    used by :py:func:`_interpret_cbs_inputs` to get default alphas for DFA, as both
    "scf" and "corl" types have their default alpha set.


    Parameters
    ----------
    xtpl_type
        {'scf', 'corl', 'fctl', 'dh', None}
        Extrapolation type: 'scf' and 'fctl' for the total energy, 'corl' and 'dh'
        for just the correlation component, None if no extrapolation requested.
    bdata
        Names and zetas of the basis sets to extrapolate from.
    funcname
        Functional name for DFT extrapolation.

    Returns
    -------
    alpha : Union[float, None]
        Extrapolation alpha to be used.
    """
    nbasis = len(bdata[0])
    if nbasis == 1:
        return None
    elif nbasis == 2:
        if xtpl_type == "scf":
            return 1.63  # this is the scf_xtpl_helgaker_2 default
        elif xtpl_type == "corl":
            return 3.00  # this is the corl_xtpl_helgaker_2 default
        elif xtpl_type in ["fctl", "dh", "dft"]:
            return _get_dfa_alpha(xtpl_type, bdata, funcname)
        else:
            raise ValidationError(
                f"Stage treatment must be one of ['scf','corl','fctl','dh',"
                f"'dft'], not '{xtpl_type}'"
            )
    elif nbasis == 3:
        return None  # we don't need alpha for _3 methods


def _interpret_cbs_inputs(
    cbs_metadata: List[dict], 
    molecule: Union[qcdb.molecule.Molecule, core.Molecule],
) -> List[dict]:
    """ A helper function which validates the ``cbs_metadata`` format,
    expands basis sets, and provides sensible defaults for optional arguments.

    There are currently three ways to run a CBS calculation:

    1) By specifying arguments, such as ``scf_wfn``, ``corl_basis``, or ``delta_scheme``.
       This is the original interface, however option (3) is preferred nowadays.
    2) By using the "method/basis+delta" syntax. This option offeres limited features,
       as the default treatments and alphas are imposed.
    3) By explicitly passing ``cbs_metadata`` as one of the arguments. This is the
       full-featured user-facing interface.

    These three are interpreted into a single, internal specification here. For the
    description of the user-facing interfaces see :py:func:`~psi4.cbs`.

    The main role of this function is to interpret the user-provided arguments,
    infer any missing stages (if possible) such as the SCF stage if the user
    provides "MP2/cc-pv[dt]z" or the split stages if user requests an extrapolation
    of a DFA.

    The resulting ``metadata`` :class:`(list[dict])` contains the list of all stages,
    with the following entries set in each of the dicts:

     - ``"wfn"`` :class:`(str)`: method name for upper level of theory (l.o.t.)
     - ``"basis"`` :class:`(tuple)`: definition of the basis sets for the upper
       l.o.t., provided as an "expanded" tuple of names and zeta levels
     - ``"options"`` :class:`(dict)`: specification of additional Psi4 options to
       switch on in the upper l.o.t.
     - ``"component"`` :class:`(str)`: component of the total E/G/H to use for
       extrapolation as the upper l.o.t.; for WFT it's usually the same as ``"wfn"``,
       for DFT it allows splitting of the total DFT energy into its components. 
       See :func:`~psi4.driver.driver_cbs.return_energy_components` for an
       enumeration of all components of each wfn.
     - ``"isdelta"`` :class:`(bool)`: a toggle determining whether the current 
       stage is to be computed as a difference between upper and lower l.o.t.;
       if ``True``, the following keys will be determined:
         
         - ``"wfn_lo"``: ``"wfn"`` analogue for the lower l.o.t. Defaults to the
           value of ``"wfn"`` from the **previous** stage.
         - ``"basis_lo"``: ``"basis"`` analogue for the lower l.o.t. Defaults to 
           the value of ``"basis"`` from the **current** stage.
         - ``"options_lo"``: ``"options"`` analogue for the lower l.o.t. Defaults
           to an empty ``{}``, i.e. no special options compared to upper l.o.t.
         - ``"component_lo"``: lower l.o.t. analogue for component. Defaults to
           the value of ``"component"`` from the **current** stage
     
     - ``"stage"`` :class:`(str)`: a string tag of the stage, used only in printing
     - ``"scheme"`` :class:`(Callable)`: the function to be used for extrapolation,
       applies for both upper and lower l.o.t.
     - ``"alpha"`` :class:`(float)`: the value of :math:`\alpha` to be passed into 
       the requested ``"scheme"``, applies for both upper and lower l.o.t.

    Parameters
    ----------
    cbs_metadata
        List of dicts containing CBS stage keywords.
    molecule
        Molecule to be passed to _expand_bracketed_basis()

    Returns
    -------
    metadata : List[dict]
        Validated list of dictionaries, with each item consisting of an extrapolation
        stage. All validation takes place here.
    """

    metadata = []
    for iitem, item in enumerate(cbs_metadata):
        # 1) all items must have wfn and basis set
        if "wfn" not in item:
            raise ValidationError(f"Stage {iitem} doesn't have defined level of theory!")
        if "basis" not in item:
            raise ValidationError(f"Stage {iitem} doesn't have defined basis sets!")
        # 2) process required stage parameters and assign defaults
        stage = {}
        stage["wfn"] = item["wfn"].lower()
        stage["basis"] = _expand_bracketed_basis(item["basis"].lower(), molecule)
        # 2a) process DFT methods
        if stage["wfn"] in functionals and stage["wfn"] not in ["hf", "scf"]:
            # 2ai) first stage - split it into components, unless "component" == "dft"
            #     This allows users to request method/basis extrapolations, such
            # as "BLYP-D3(0)/pc-[12]", without specifying extra parameters. If 
            # "component" is not "dft", the "fctl", "disp", "dh", and "nl" stages 
            # are generated automatically as required, and the "DFT TOTAL ENERGY" 
            # variable is not used in the extrapolation.
            #     In cases where "component" is "dft", that setting gets processed,
            # and further generation of component stages is skipped.
            #     Note that none of these stages are delta-stages, as only one E/G/H
            # call is made and the components are picked from the returned wfn
            # automatically.
            if len(metadata) == 0:
                if "component" not in item:
                    fctl = {
                        "wfn": stage["wfn"],
                        "basis": stage["basis"],
                        "component": f'{stage["wfn"]}-fctl',
                        "stage": item.get("stage", "fctl"),
                        "isdelta": False,
                        "scheme": _get_default_xtpl(len(stage["basis"][1]), "fctl"),
                        "options": item.get("options", False),
                        "alpha": item.get("alpha", _get_default_alpha("fctl", stage["basis"], stage["wfn"]))
                    }
                    metadata.append(fctl)
                    if f'{stage["wfn"]}-dh' in VARH[stage["wfn"]]:
                        dh = {
                            "wfn": stage["wfn"],
                            "basis": stage["basis"],
                            "component": f'{stage["wfn"]}-dh',
                            "stage": 'dh',
                            "isdelta": False,
                            "scheme": _get_default_xtpl(len(stage["basis"][1]), "dh"),
                            "options": item.get("options", False),
                            "alpha": item.get("alpha", _get_default_alpha("dh", stage["basis"], stage["wfn"]))
                        }
                        metadata.append(dh)
                    if f'{stage["wfn"]}-nl' in VARH[stage["wfn"]]:
                        nl = {
                            "wfn": stage["wfn"],
                            "basis": ([stage["basis"][0][-1]], [stage["basis"][1][-1]]),
                            "component": f'{stage["wfn"]}-nl',
                            "stage": 'nl',
                            "isdelta": False,
                            "scheme": _get_default_xtpl(1, None),
                            "options": item.get("options", False),
                            "alpha": item.get("alpha", None)
                        }
                        metadata.append(nl)
                    if f'{stage["wfn"]}-disp' in VARH[stage["wfn"]]:
                        di = {
                            "wfn": stage["wfn"],
                            "basis": ([stage["basis"][0][-1]], [stage["basis"][1][-1]]),
                            "component": f'{stage["wfn"]}-disp',
                            "stage": 'disp',
                            "isdelta": False,
                            "scheme": _get_default_xtpl(1, None),
                            "options": item.get("options", False),
                            "alpha": item.get("alpha", None)
                        }
                        metadata.append(di)
                elif item["component"] in ["dft", "fctl"]:
                    fctl = {
                        "wfn": stage["wfn"],
                        "basis": stage["basis"],
                        "component": f'{stage["wfn"]}-{item["component"]}',
                        "stage": item.get("stage", item["component"]),
                        "isdelta": False,
                        "scheme": _get_default_xtpl(len(stage["basis"][1]), item["component"]),
                        "options": item.get("options", False),
                        "alpha": item.get("alpha", _get_default_alpha(item["component"], stage["basis"], stage["wfn"]))
                    }
                    metadata.append(fctl)
                else:
                    raise ValidationError(
                        f"Component '{item['component']}' is not a valid choice "
                        "for the first stage in extrapolation."
                    )
            # 2aii) further stages - treat as delta stages
            else:
                delta = {"wfn": stage["wfn"], "basis": stage["basis"]}
                delta["component"] = f'{delta["wfn"]}-{item.get("component", "dft")}'
                delta["wfn_lo"] = item.get("wfn_lo", metadata[-1]["wfn"]).lower()
                delta["basis_lo"] = _expand_bracketed_basis(item.get("basis_lo", item["basis"]).lower(), molecule)
                if len(delta["basis"][0]) != len(delta["basis_lo"][0]):
                    raise ValidationError("""Number of basis sets inconsistent
                                                between high ({}) and low ({}) levels.""".format(
                        len(stage["basis"][0]), len(stage["basis_lo"][0])))
                delta["stage"] = item.get("stage", f"delta_{len(metadata)}")
                # This allows for doing true delta energies (from two calculations, by passing no component or "dft")
                # as well as single-stage extrapolated "deltas" for ["dh", "disp", "nl"] components.
                delta["isdelta"] = True if item.get("component", "dft") in ["fctl", "dft"] else False
                delta["scheme"] = _get_default_xtpl(len(stage["basis"][1]), item.get("treatment", "scf"))
                delta["options"] = item.get("options", False)
                delta["options_lo"] = item.get("options_lo", False)
                delta["alpha"] = item.get("alpha", None)
                delta["component_lo"] = f'{delta["wfn_lo"]}-{item.get("component", "dft")}'
                metadata.append(delta)
        # 2b) process WFT methods
        else:
            # 2bi) generate first item ("scf" stage) if necessary
            if len(metadata) == 0 and stage["wfn"] not in ["hf", "c4-hf", "scf", "c4-scf"]:
                scf = {}
                if stage["wfn"].startswith("c4"):
                    scf["wfn"] = "c4-hf"
                else:
                    scf["wfn"] = "hf"
                scf["basis"] = ([stage["basis"][0][-1]], [stage["basis"][1][-1]])
                scf["isdelta"] = False
                scf["stage"] = "scf"
                scf["scheme"] = _get_default_xtpl(1, "scf")
                scf["alpha"] = None
                scf["options"] = False
                scf["options_lo"] = False
                scf["component"] = scf["wfn"]
                metadata.append(scf)
            # 2bii) keep processing current stage
            stage["component"] = stage["wfn"]
            stage["isdelta"] = False
            stage["stage"] = item.get("stage", False)
            if not stage["stage"]:
                if len(metadata) == 0:
                    stage["stage"] = "scf"
                elif len(metadata) == 1:
                    stage["stage"] = "corl"
                else:
                    stage["stage"] = f"delta{len(metadata) - 1}"
            if len(metadata) > 0:
                stage["isdelta"] = True
                stage["wfn_lo"] = item.get("wfn_lo", metadata[-1].get("wfn")).lower()
                stage["component_lo"] = stage["wfn_lo"]
                stage["basis_lo"] = _expand_bracketed_basis(item.get("basis_lo", item["basis"]).lower(), molecule)
                if len(stage["basis"][0]) != len(stage["basis_lo"][0]):
                    raise ValidationError("""Number of basis sets inconsistent
                                                between high ({}) and low ({}) levels.""".format(
                        len(stage["basis"][0]), len(stage["basis_lo"][0])))
            stage["scheme"] = item.get("scheme", None)
            stage["alpha"] = item.get("alpha", None)
            if stage["scheme"] is None:
                stage["scheme"] = _get_default_xtpl(len(stage["basis"][1]),
                                                    item.get("treatment", "scf" if len(metadata) == 0 else "corl"))
                if stage["alpha"] is None:
                    stage["alpha"] = _get_default_alpha(item.get("treatment", "scf" if len(metadata) == 0 else "corl"),
                                                        stage["basis"], stage["wfn"])
            stage["options"] = item.get("options", False)
            stage["options_lo"] = item.get("options_lo", False)
            metadata.append(stage)
    return (metadata)


def _process_cbs_kwargs(kwargs: dict) -> List[dict]:
    """ A helper function which translates supplied kwargs into the
    ``cbs_metadata`` format and passes it for further validation.

    Parameters
    ----------
    kwargs
        kwargs containing the CBS function specification.

    Returns
    -------
    metadata : List[dict]
        List of dictionaries, with each item consisting of an extrapolation stage.
    """

    molecule = kwargs.get('molecule', core.get_active_molecule())

    if "cbs_metadata" in kwargs:
        # if we passed in a dict, validate it right away
        cbs_metadata = kwargs["cbs_metadata"]
    else:
        # if we passed in options, check for consecutive correlations first
        if "delta_wfn" in kwargs and "corl_wfn" not in kwargs:
            raise ValidationError("Delta function supplied without corl_wfn defined.")
        if "delta2_wfn" in kwargs and "delta_wfn" not in kwargs:
            raise ValidationError("Second delta function supplied without delta_wfn defined.")
        cbs_metadata = []
        possible_stages = ["scf", "corl"]
        while len(possible_stages) > 0:
            sn = possible_stages.pop(0)
            if f"{sn}_wfn" in kwargs and f"{sn}_basis" in kwargs:
                # either both *_wfn and *_basis have to be specified
                stage = {"wfn": kwargs[f"{sn}_wfn"], "basis": kwargs[f"{sn}_basis"]}
            elif sn == "scf" and f"{sn}_basis" in kwargs:
                # or we're at a scf stage which can be implied with a provided scf_basis
                stage = {"wfn": "hf", "basis": kwargs[f"{sn}_basis"]}
            else:
                # otherwise go to the next possible stage
                continue
            # if we made it here, stage exists - parse other keywords
            if f"{sn}_scheme" in kwargs:
                stage["scheme"] = kwargs[f"{sn}_scheme"]
            if f"{sn}_wfn_lesser" in kwargs:
                stage["wfn_lo"] = kwargs[f"{sn}_wfn_lesser"]
            if f"cbs_{sn}_alpha" in kwargs:
                stage["alpha"] = kwargs[f"cbs_{sn}_alpha"]
            elif f"{sn}_alpha" in kwargs:
                stage["alpha"] = kwargs[f"{sn}_alpha"]
            cbs_metadata.append(stage)
            if sn == "corl":
                possible_stages.append("delta")
            elif sn == "delta":
                possible_stages.append("delta2")

    return _interpret_cbs_inputs(cbs_metadata, molecule)


###################################
##  Start of Complete Basis Set  ##
###################################


def cbs(func, label, **kwargs):
    r"""Function to define a multistage energy method from combinations of
    basis set extrapolations and delta corrections and condense the components 
    into a minimum number of calculations.

    :aliases: complete_basis_set()

    :returns: (*float*) -- Total electronic energy in Hartrees

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`CBS TOTAL ENERGY`
       * :psivar:`CBS REFERENCE ENERGY`
       * :psivar:`CBS CORRELATION ENERGY`
       * :psivar:`CURRENT ENERGY`
       * :psivar:`CURRENT REFERENCE ENERGY`
       * :psivar:`CURRENT CORRELATION ENERGY`

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - No way to tell function to boost fitting basis size for all calculations.

       - Need to add more extrapolation schemes

    As represented in the equation below, a CBS energy method is defined in 
    several sequential stages. For WFT they are ``scf``, ``corl``, ``delta1``, 
    ``delta2``, ... , covering the treatment of the reference total energy, 
    the correlation energy, and any further delta corrections to the correlation 
    or the reference energies. Each stage may be activated by supplying the
    stage_wfn argument, or it can be defined using a :class:`dict` in a field of 
    the ``cbs_metadata`` list. Delta stages are only allowed if all preceding 
    stages are active.

    .. include:: /cbs_eqn.rst

    .. note::
        In the following, ``<stage>`` corresponds to one of ``scf``, ``corl``,
        ``delta1``, or ``delta2``, as appropriate.

    The extrapolation of DFT calculations from basis sets of successive zeta 
    levels is also implemented. By default, the components of the DFA recipe
    (``fctl``, ``dh``, ``disp``, ``nl``) are treated separately: the ``fctl`` 
    component, which is analogous to the ``scf`` stage in WFT, is extrapolated 
    using the exponential -- square root function; the ``dh`` component, which 
    is analogous to the ``corl`` stage in WFT, is extrapolated using the power 
    function; and the ``disp`` and ``nl`` components are not extrapolated at all.
    For further details, see :ref:`below<cbs-dft>`.

    .. _cbs-energy-methods:

    Energy Methods
    --------------
    
    The presence of any ``<stage>_wfn`` argument is the indicator to incorporate 
    (and check for matching ``<stage>_basis`` and ``<stage>_scheme`` args) and 
    compute that stage as an element of the CBS energy.

    The cbs() function requires, at a minimum, ``name='scf'`` and ``scf_basis``
    arguments to be specified for reference-stage only jobs and ``name`` and
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
           * all built-in DFAs

    :type name: str
    :param name: ``'scf'`` || ``'ccsd'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        for the correlation energy, unless only reference step to be performed,
        in which case should be ``'scf'``. Overruled if any ``<stage>_wfn`` arg 
        is supplied.

    :type scf_wfn: str
    :param scf_wfn: |dl| ``'scf'`` |dr| || ``'c4-scf'`` || etc.

        Indicates the energy method for which the reference energy is to be
        obtained. Generally unnecessary, as 'scf' is *the* scf in |PSIfour| but
        can be used to direct lone scf components to run in |PSIfour| or Cfour
        in a mixed-program composite method.

    :type corl_wfn: str
    :param corl_wfn: ``'mp2'`` || ``'ccsd(t)'`` || etc.

        Indicates the energy method for which the correlation energy is to be
        obtained. Can also be specified with ``name`` or as the unlabeled
        first argument to the function.

    :type delta_wfn: str
    :param delta_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a delta correction
        to the correlation energy is to be obtained.

    :type delta_wfn_lesser: str
    :param delta_wfn_lesser: |dl| ``corl_wfn`` |dr| || ``'mp2'`` || etc.

        Indicates the inferior energy method for which a delta correction
        to the correlation energy is to be obtained.

    :type delta2_wfn: str
    :param delta2_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a second delta correction
        to the correlation energy is to be obtained.

    :type delta2_wfn_lesser: str
    :param delta2_wfn_lesser: |dl| ``delta_wfn`` |dr| || ``'ccsd(t)'`` || etc.

        Indicates the inferior energy method for which a second delta correction
        to the correlation energy is to be obtained.

    Basis Sets
    ----------

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
    
    .. note::
        Currently, setting the basis set using the ``set`` command has no influence
        on a cbs calculation.

    Schemes
    -------

    Schemes are functions that define the transformations of the energy (or gradient) 
    through basis set extrapolation for each stage of the CBS definition. A complaint 
    is generated if number of basis sets in ``<stage>_basis`` does not exactly 
    satisfy the requirements of ``<stage>_scheme``. An exception is the default, 
    ``'xtpl_highest_1'``, which uses the highest-zeta basis set from those supplied.
    See :ref:`sec:cbs_xtpl` for all available schemes.

    :type scf_scheme: Callable
    :param scf_scheme: |dl| ``xtpl_highest_1`` |dr| || ``scf_xtpl_helgaker_3`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the 
        reference energy. Defaults to :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_3` 
        if three valid basis sets present in ``psi4.driver.driver_cbs.scf_basis``, 
        :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_2` if two valid basis
        sets present in ``scf_basis``, and :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1` 
        otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs.xtpl_exponential_2` 
           * :py:func:`~psi4.driver.driver_cbs.xtpl_power_2` 
           * :py:func:`~psi4.driver.driver_cbs.xtpl_expsqrt_2` 
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_2`
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_truhlar_2`
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_karton_2`
           * :py:func:`~psi4.driver.driver_cbs.xtpl_exponential_3`
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_3`

    :type corl_scheme: Callable
    :param corl_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Te basis set extrapolation scheme to be applied to the 
        correlation energy. Defaults to :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2` 
        if two valid basis sets present in ``corl_basis`` and 
        :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1` otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs.xtpl_power_2`
           * :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2`

    :type delta_scheme: Callable
    :param delta_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        The basis set extrapolation scheme to be applied to the delta 
        correction to the correlation energy. Defaults to 
        :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2` if two valid basis 
        sets present in ``delta_basis`` and :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1` 
        otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs.xtpl_power_2` 
           * :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2`

    :type delta2_scheme: Callable
    :param delta2_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        The basis set extrapolation scheme to be applied to the second 
        delta correction to the correlation energy. Defaults to 
        :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2` if two valid 
        basis sets present in ``delta2_basis`` and 
        :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1` otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs.xtpl_power_2` 
           * :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2`

    :type scf_alpha: float
    :param scf_alpha: |dl| ``1.63`` |dr|

        Overrides the default :math:`\alpha` parameter used in the listed SCF 
        extrapolation procedures. Does not affect other procedures, including 
        :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1` and 
        :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_3`.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_2`
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_truhlar_2`
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_karton_2`

    :type corl_alpha: float
    :param corl_alpha: |dl| ``3.00`` |dr|

        Overrides the default :math:`\alpha` parameter used in the listed 
        :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2` correlation
        extrapolation to the corl stage. The supplied :math:`\alpha` does not 
        impact delta or any further stages.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2`

    :type delta_alpha: float
    :param delta_alpha: |dl| ``3.00`` |dr|

        Overrides the default :math:`\alpha` parameter used in the listed
        :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2` correlation 
        extrapolation for the delta correction. Useful when delta correction is 
        performed using smaller basis sets for which a different :math:`\alpha` 
        might be more appropriate.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2`

    Combined interface
    ------------------

    This user-facing interface allows for the specification of more detailed 
    extrapolation recipes, and is a superset of the keyword functionality above.
    Internally, all calls to cbs() get first translated into this user-facing form, 
    before being validated and interpreted into an internal specification discussed
    in :py:func:`~psi4.driver.driver_cbs._interpret_cbs_inputs`.

    :type cbs_metadata: list(dict)
    :param cbs_metadata: |dl| autogenerated from above keywords |dr| || ``[{"wfn": "hf", "basis": "cc-pv[TQ5]z"}]`` || etc.

        Each item in this list is a dict, corresponding to a stage in the 
        extrapolation recipe. For discussion of DFT, see :ref:`further below <cbs-dft>`.
        For WFT, the first item in the list is always defining the SCF contribution to 
        the total energy, and will be generated automatically if not supplied. 
        The latter list elements then define correlated or delta stages, with 
        keywords filled-in automatically from previous stages, unless supplied 
        by the user. 
        
        The required keywords in each stage are:

        * ``"wfn"``: method name, for WFT this is typically ``"HF"`` or one of the
          other energy methods in :ref:`cbs-energy-methods`. Note that HF is 
          subsumed in correlated methods anyway.
        * ``"basis"``: basis set, can be in a bracketed form (``"cc-pv[tq]z"``)

        Other supported keywords for all stages are:

        * ``"scheme"``: the extrapolation function (exponential, power, etc.). By 
          default it is assigned based on the number of basis sets (1 - 3) 
          supplied in ``"basis"`` and the ``"treatment"`` keyword for the stage.
        * ``"alpha"``: the :math:`\alpha\` for the above scheme, if the default 
          is to be overriden.
        * ``"options"``: if special Psi4 options are required for a stage, they 
          can be entered as a dict here. For corl or delta stages, if an option
          should be used for both parts of the stage, it should be entered in 
          both ``"options"`` and ``"options_lo"``. This functionality is helpful
          for calculating all electron corrections in otherwise frozen core 
          calculations, or relativistic (DKH) Hamiltionian corrections for 
          otherwise nonrelativistic calculations.
        * ``"treatment"``: treat the extrapolation stage as ``"scf"`` or ``"corl"``. 
          Used only for WFT calculations. By default, only the first stage is 
          treated as ``"scf"``, every latter one is treated as ``"corl"``. 
          Processed by the metadata validator to infer the ``"scheme"``.
        * ``"stage"``: a str tag for the stage written in the output tables.

        From second stage onwards, the following keywords can be specified:

        * ``"wfn_lo"``: the lower method from which a delta correction is to be 
          calculated. By default, it is set to ``"wfn"`` from the previous stage
          in the ``cbs_metadata`` list: therefore a ``"corl"`` stage would be 
          ``"wfn": "MP2" - "wfn_lo":"HF"``.
        * ```basis_lo```: basis set to be used for the delta correction. By 
          default, it is the same as the ``"basis"`` in the current stage 
          specified above.
        * ``"options_lo"``: special options for the lower method in a given stage. 
          This is useful when options for both the lower and upper method in a
          single stage should differ from the global options, e.g. to calculate 
          a stage using direct integrals in an otherwise density-fitted 
          calculation.
    
    .. _cbs-dft:

        For DFT calculations, the components of the total DFT energy obtained in 
        a single calculation may be extrapolated using different formulas. To 
        avoid recomputing, the ``"component"`` keyword can be specified for each
        stage. The allowed values are:

        * ``"dft"``, specifying the total energy, i.e. the ``DFT TOTAL ENERGY``,
        * ``"fctl"``, specifying the base functional SCF energy excluding non-
          local corrections, i.e. ``DFT FUNCTIONAL TOTAL ENERGY`` less any
          ``DFT VV10 ENERGY``
        * ``"dh"``, specifying the perturbative correlation in double hybrid
          functionals, i.e. the ``DOUBLE-HYBRID CORRECTION ENERGY``
        * ``"disp"``, specifying dispersion corrections such as -D3(0), i.e.
          the ``DISPERSION CORRECTION ENERGY``
        * ``"nl"``, specifying the non-local correction energy, i.e. the
          ``DFT VV10 ENERGY``
        
        Unless specified by the user, any ``cbs()`` calculation requesting a DFA
        as a method will generate a stage with ``"fctl"`` treatment automatically,
        using :py:func:`~psi4.driver.driver_cbs.xtpl_expsqrt_2` and 
        :py:func:`~psi4.driver.driver_cbs._get_dfa_alpha` to extrapolate. The
        ``"dh"``, ``"disp"``, or ``"nl"`` stages are generated if those components
        are present in the DFA. The ``"dh"`` stage is extrapolated using 
        :py:func:`~psi4.driver.driver_cbs.xtpl_power_2`, while the ``"disp"`` and 
        ``"nl"`` stages are not extrapolated by default. See example [10] below.

        DFA delta stages, i.e. stages where both ``"wfn"`` and ``"wfn_lo"`` are 
        supplied, are understood by the parser and can be computed; see example 
        [11] below.

    Others
    ------

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :examples:


    >>> # [1] replicates with cbs() the simple model chemistry scf/cc-pVDZ: set basis cc-pVDZ energy('scf')
    >>> energy(cbs, scf_wfn='scf', scf_basis='cc-pVDZ')

    >>> # [2] replicates with cbs() the simple model chemistry mp2/jun-cc-pVDZ: set basis jun-cc-pVDZ energy('mp2')
    >>> energy(cbs, corl_wfn='mp2', corl_basis='jun-cc-pVDZ')

    >>> # [3] DTQ-zeta extrapolated scf reference energy
    >>> energy(cbs, scf_wfn='scf', scf_basis='cc-pV[DTQ]Z', scf_scheme=scf_xtpl_helgaker_3)

    >>> # [4] DT-zeta extrapolated mp2 correlation energy atop a T-zeta reference
    >>> energy(cbs, corl_wfn='mp2', corl_basis='cc-pv[dt]z', corl_scheme=corl_xtpl_helgaker_2)

    >>> # [5] a DT-zeta extrapolated coupled-cluster correction atop a TQ-zeta extrapolated mp2 correlation energy atop a Q-zeta reference (both equivalent)
    >>> energy(cbs, corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z')
    >>> energy(cbs, corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z', delta_scheme=corl_xtpl_helgaker_2)

    >>> # [6] a D-zeta ccsd(t) correction atop a DT-zeta extrapolated ccsd cluster correction atop a TQ-zeta extrapolated mp2 correlation energy atop a Q-zeta reference
    >>> energy(cbs, corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd', delta_basis='aug-cc-pv[dt]z', delta_scheme=corl_xtpl_helgaker_2, delta2_wfn='ccsd(t)', delta2_wfn_lesser='ccsd', delta2_basis='aug-cc-pvdz')

    >>> # [7] a Q5-zeta MP2 calculation, corrected by CCSD(T) at the TQ-zeta extrapolated level, and all-electron CCSD(T) correlation at T-zeta level
    >>> energy(cbs, cbs_metadata=[{"wfn": "hf", "basis": "cc-pv5z"}, {"wfn": "mp2", "basis": "cc-pv[q5]z"}, {"wfn": "ccsd(t)", "basis": "cc-pv[tq]z"}, {"wfn": "ccsd(t)", "basis": "cc-pvtz", "options": {"freeze_core": "False"}}])

    >>> # [8] cbs() coupled with database()
    >>> TODO database('mp2', 'BASIC', subset=['h2o','nh3'], symm='on', func=cbs, corl_basis='cc-pV[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd(t)', delta_basis='sto-3g')

    >>> # [9] cbs() coupled with optimize()
    >>> TODO optimize('mp2', corl_basis='cc-pV[DT]Z', corl_scheme=corl_xtpl_helgaker_2, func=cbs)

    >>> # [10] basis set extrapolation for a dispersion-corrected double-hybrid DFA with fctl and dh extrapolated separately
    >>> energy('B2PLYP-D3BJ/def2-[tq]zvpd')

    >>> # [11] cbs() of BLYP-D3/cc-pv[dt]z with a delta-correction from B2PLYP-D3
    >>> E = energy(cbs, cbs_metadata=[{'wfn': "blyp-d3", 'basis': "cc-pv[dt]z"},
    >>>                               {'wfn': "b2plyp-d3", 'basis': "cc-pvdz"}])


    """
    kwargs = p4util.kwargs_lower(kwargs)

    # All kwargs get processed using _process_cbs_kwargs(), which calls 
    # _interpret_cbs_inputs(), ensuring that the user-facing spec is always correctly
    # converted into the internal spec. Note that the "method/basis" type inputs
    # get processed into kwargs in _cbs_gufunc() before ever reaching cbs().
    metadata = _process_cbs_kwargs(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    verbose = kwargs.pop('verbose', 0)
    ptype = kwargs.pop('ptype')

    if ptype not in ['energy', 'gradient', 'hessian']:
        raise ValidationError("""Wrapper complete_basis_set is unhappy to be calling
                                 function '%s' instead of 'energy', 'gradient' or 'hessian'.""" % ptype)
    # Abort here for dft extrapolations, as analytic Hessians are not yet implemented.
    elif ptype in ['hessian'] and metadata[0]["wfn"] in functionals and metadata[0]["wfn"] not in ["hf", "scf"]:
        raise ValidationError("""Wrapper complete_basis_set is unhappy to be calling
                                 function '%s' instead of 'energy'. Try with 'dertype=0' """ % ptype)

    optstash = p4util.OptionsState(['BASIS'], ['WFN'], ['WRITER_FILE_LABEL'])

    # Define some quantum chemical knowledge, namely what methods are subsumed in others

    user_writer_file_label = core.get_global_option('WRITER_FILE_LABEL')

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()
    natom = molecule.natom()

    if metadata[0]["wfn"] not in VARH.keys():
        raise ValidationError(
            """Requested SCF method '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" %
            (metadata[0]["wfn"]))

    if len(metadata) > 1:
        for delta in metadata[1:]:
            if delta["wfn"] not in VARH.keys():
                raise ValidationError(
                    """Requested higher %s method '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" %
                    (delta["treatment"], delta["wfn"]))
            if delta["isdelta"] and delta["wfn_lo"] not in VARH.keys():
                raise ValidationError(
                    """Requested lesser %s method '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" %
                    (delta["treatment"], delta["wfn_lo"]))

    # Build string of title banner
    instructions = "\n" + p4util.banner(f" CBS Setup{':' + label if label else ''} ", strNotOutfile=True) + "\n"
    core.print_out(instructions)

    # Call schemes for each portion of total energy to 'place orders' for calculations needed
    d_fields = [
        'd_stage', 'd_isdelta', 'd_scheme', 'd_basis', 'd_wfn', 'd_component', 'd_need', 'd_coef', 'd_energy',
        'd_gradient', 'd_hessian', 'd_alpha'
    ]
    f_fields = ['f_wfn', 'f_component', 'f_basis', 'f_zeta', 'f_energy', 'f_gradient', 'f_hessian', 'f_options']
    GRAND_NEED = []
    MODELCHEM = []

    for stage in metadata:
        NEED = _expand_scheme_orders(
            stage["scheme"], 
            stage["basis"][0], 
            stage["basis"][1], 
            stage["wfn"],
            stage["options"],
            natom
        )
        GRAND_NEED.append(dict(zip(
            d_fields,
            [
                stage["stage"], 
                stage["isdelta"], 
                stage["scheme"],
                _contract_bracketed_basis(stage["basis"][0]), 
                stage["wfn"], 
                stage["component"], 
                NEED, 
                +1, 
                0.0,
                None, 
                None, 
                stage["alpha"]
            ]
        )))
        if stage["isdelta"]:
            NEED = _expand_scheme_orders(
                stage["scheme"], 
                stage["basis_lo"][0], 
                stage["basis_lo"][1], 
                stage["wfn_lo"],
                stage["options_lo"],
                natom
            )
            GRAND_NEED.append(dict(zip(
                d_fields, 
                [
                    stage["stage"], 
                    stage["isdelta"], 
                    stage["scheme"],
                    _contract_bracketed_basis(stage["basis_lo"][0]), 
                    stage["wfn_lo"], 
                    stage["component_lo"], 
                    NEED,
                    -1, 
                    0.0, 
                    None, 
                    None, 
                    stage["alpha"]
                ]
            )))

    for stage in GRAND_NEED:
        for lvl in stage['d_need'].items():
            lvl[1]['f_component'] = stage["d_component"]
            MODELCHEM.append(lvl[1])

    # Apply chemical reasoning to choose the minimum computations to run
    JOBS = MODELCHEM[:]

    addlremark = {'energy': '', 'gradient': ', GRADIENT', 'hessian': ', HESSIAN'}
    instructions = ''
    instructions += """    Naive listing of computations required.\n"""
    for mc in JOBS:
        instructions += """   %12s / %-24s for  %s%s\n""" % \
            (mc['f_wfn'], mc['f_basis'] + " + options"*bool(mc['f_options']),
             VARH[mc['f_wfn']][mc['f_component']], addlremark[ptype])

    #     Remove duplicate modelchem portion listings
    for mc in MODELCHEM:
        dups = []
        for indx_job, job in enumerate(JOBS):
            if (
                job['f_wfn'] == mc['f_wfn'] and 
                job['f_basis'] == mc['f_basis'] and 
                job['f_options'] == mc['f_options'] and 
                (job['f_options'] != False or ptype == 'energy')
            ):
                # the highest priority "duplicate" with matching f_wfn is prepended
                dups.insert(0, indx_job) 
    #     Remove chemically subsumed modelchem portion listings for energies and DFT
            elif (
                ptype == "energy" or 
                (mc['f_wfn'] in functionals and mc['f_wfn'] not in ["scf", "hf"])
            ):
                if ( 
                    job['f_component'] != job['f_wfn'] and 
                    job['f_wfn'] == mc['f_wfn'] and
                    job['f_basis'] == mc['f_basis'] and 
                    job['f_options'] == False and 
                    job['f_options'] == mc['f_options']
                ):
                    dups.append(indx_job)
                else:
                    for component in VARH[mc['f_wfn']]:
                        if (
                            job['f_basis'] == mc['f_basis'] and 
                            job['f_options'] == False and 
                            job['f_options'] == mc['f_options'] and
                            job['f_component'] == component and 
                            job['f_wfn'] != mc['f_wfn']
                        ):
                            dups.append(indx_job)
        if len(dups) > 1:
            # keep the first item and delete others in reverse order to avoid oob error
            for d in dups[1:][::-1]:
                del JOBS[d]
            
    instructions += """\n    Enlightened listing of computations required.\n"""
    for mc in JOBS:
        instructions += """   %12s / %-24s for  %s%s\n""" % \
            (mc['f_wfn'], mc['f_basis'] + " + options"*bool(mc['f_options']),
             VARH[mc['f_wfn']][mc['f_component']], addlremark[ptype])

    #     Expand listings to all that will be obtained
    JOBS_EXT = []
    for job in JOBS:
        for component in VARH[job['f_wfn']]:
            JOBS_EXT.append(
                dict(
                    zip(f_fields, [
                        job['f_wfn'], component, job['f_basis'], job['f_zeta'], 0.0,
                        core.Matrix(natom, 3),
                        core.Matrix(3 * natom, 3 * natom), job['f_options']
                    ])))

    instructions += """\n    Full listing of computations to be obtained (required and bonus).\n"""
    for mc in JOBS_EXT:
        instructions += """   %12s / %-24s for  %s%s\n""" % \
            (mc['f_wfn'], mc['f_basis'] + " + options"*bool(mc['f_options']),
             VARH[mc['f_wfn']][mc['f_component']], addlremark[ptype])
    core.print_out(instructions)

    psioh = core.IOManager.shared_object()
    psioh.set_specific_retention(psif.PSIF_SCF_MOS, True)
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
            (mc['f_wfn'].upper(), mc['f_basis'].upper() + " + opts."*bool(mc['f_options']), addlremark[ptype])
        cbsbanners += """core.print_out('\\n')\n\n"""
        exec(cbsbanners)

        # Build string of molecule and commands that are dependent on the database
        commands = '\n'
        commands += """\ncore.set_global_option('BASIS', '%s')\n""" % (mc['f_basis'])
        commands += """core.set_global_option('WRITER_FILE_LABEL', '%s')\n""" % \
            (user_writer_file_label + ('' if user_writer_file_label == '' else '-') + mc['f_wfn'].lower() + '-' + mc['f_basis'].lower().replace('*', 's'))
        exec(commands)

        # Stash and set options if any
        if mc["f_options"]:
            optionstash = p4util.OptionsState(*[[opt] for opt in list(mc["f_options"])])
            for k, v, in mc["f_options"].items():
                core.set_global_option(k.upper(), v)
        else:
            optionstash = False

        # Make energy(), etc. call
        response, wfn = func(molecule=molecule, return_wfn = True, **kwargs)
        if ptype == 'energy':
            mc['f_energy'] = wfn.variable(VARH[mc['f_wfn']][mc['f_component']])
        elif ptype == 'gradient':
            # For DFT gradients, we handle all functional components below.
            mc['f_gradient'] = response
            mc['f_energy'] = wfn.variable('CURRENT ENERGY')
            if verbose > 1:
                mc['f_gradient'].print_out()
        elif ptype == 'hessian':
            mc['f_hessian'] = response
            mc['f_energy'] = wfn.variable('CURRENT ENERGY')
            if verbose > 1:
                mc['f_hessian'].print_out()
        Njobs += 1
        if verbose > 1:
            core.print_out("\nCURRENT ENERGY: %14.16f\n" % mc['f_energy'])

        # Restore modified options
        if optionstash:
            optionstash.restore()

        # Fill in energies for subsumed methods
        # For DFT, we need to "correct" the -fctl energy by subtracting the -nl component.
        # The "DFT VV10 ENERGY" is by default = 0.0, even for DFAs without -nl
        if ptype == 'energy':
            for component in VARH[mc['f_wfn']]:
                for job in JOBS_EXT:
                    if (mc['f_wfn'] == job['f_wfn']) and (mc['f_basis'] == job['f_basis']) and \
                       (component == job['f_component']) and (mc['f_options'] == job['f_options']):
                        job['f_energy'] = wfn.variable(VARH[mc['f_wfn']][component])
                        if component.endswith("fctl"):
                            job['f_energy'] -= wfn.variable("DFT VV10 ENERGY")
        # For DFT, we have deleted duplicate gradient calls but the gradient components
        # are available as variables, with ENERGY-> GRADIENT. So if they are required, let us fill them.
        # At the moment, only "disp" gradients are possible ("nl" and "dh" not yet implemented):
        elif ptype == 'gradient':
            for component in VARH[mc['f_wfn']]:
                if component.endswith("disp") or component.endswith("fctl"):
                    for job in JOBS_EXT:
                        if (mc['f_wfn'] == job['f_wfn']) and (mc['f_basis'] == job['f_basis']) and \
                        (component == job['f_component']) and (mc['f_options'] == job['f_options']):
                            job['f_energy'] = wfn.variable(VARH[mc['f_wfn']][component])
                            job['f_gradient'] = wfn.variable(VARH[mc['f_wfn']][component].replace(
                                "ENERGY", "GRADIENT"))

        if verbose > 1:
            core.print_variables()
        core.clean_variables()
        core.clean()

        # Copy data from 'run' to 'obtained' table. Here we exclude DFT components, as they have been correctly
        # assigned to JOBS_EXT in the block above, and we'd overwrite them with the wrong variable.
        for mce in JOBS_EXT:
            if (mc['f_wfn'] == mce['f_wfn']) and (mc['f_basis'] == mce['f_basis']) and \
               (mc['f_component'] == mce['f_component']) and (mc['f_options'] == mce['f_options']) and \
               not (mc['f_component'].endswith("disp") or mc['f_component'].endswith("fctl")):
                mce['f_energy'] = mc['f_energy']
                mce['f_gradient'] = mc['f_gradient']
                mce['f_hessian'] = mc['f_hessian']

    psioh.set_specific_retention(psif.PSIF_SCF_MOS, False)

    # Build string of title banner
    instructions = "\n" + p4util.banner(f" CBS Results{':' + label if label else ''} ", strNotOutfile=True) + "\n"
    core.print_out(instructions)

    # Insert obtained energies into the array that stores the cbs stages
    for stage in GRAND_NEED:
        for lvl in stage['d_need'].items():
            MODELCHEM.append(lvl[1])

            for job in JOBS_EXT:
                # Dont ask
                if (((lvl[1]['f_component'] == job['f_component']) or
                     ((lvl[1]['f_component'][3:] == job['f_component']) and lvl[1]['f_component'].startswith('c4-')) or
                     ((lvl[1]['f_component'] == job['f_component'][3:]) and job['f_component'].startswith('c4-')) or
                     (('c4-' + lvl[1]['f_component']) == job['f_component']) or
                     (lvl[1]['f_component'] == ('c4-' + job['f_component']))) and (lvl[1]['f_basis'] == job['f_basis'])
                        and (lvl[1]['f_options'] == job['f_options'])):
                    lvl[1]['f_energy'] = job['f_energy']
                    lvl[1]['f_gradient'] = job['f_gradient']
                    lvl[1]['f_hessian'] = job['f_hessian']

    # Make xtpl() call
    finalenergy = 0.0
    finalgradient = core.Matrix(natom, 3)
    finalhessian = core.Matrix(3 * natom, 3 * natom)

    for stage in GRAND_NEED:
        hiloargs = {'alpha': stage['d_alpha']}

        hiloargs.update(_contract_scheme_orders(stage['d_need'], 'f_energy'))
        stage['d_energy'] = stage['d_scheme'](**hiloargs)
        finalenergy += stage['d_energy'] * stage['d_coef']

        if ptype == 'gradient':
            hiloargs.update(_contract_scheme_orders(stage['d_need'], 'f_gradient'))
            stage['d_gradient'] = stage['d_scheme'](**hiloargs)
            work = stage['d_gradient'].clone()
            work.scale(stage['d_coef'])
            finalgradient.add(work)

        elif ptype == 'hessian':
            hiloargs.update(_contract_scheme_orders(stage['d_need'], 'f_hessian'))
            stage['d_hessian'] = stage['d_scheme'](**hiloargs)
            work = stage['d_hessian'].clone()
            work.scale(stage['d_coef'])
            finalhessian.add(work)

    # Build string of results table
    table_delimit = '  ' + '-' * 105 + '\n'
    tables = ''
    tables += """\n   ==> %s <==\n\n""" % ('Components')
    tables += table_delimit
    tables += """     %6s %20s %1s %-26s %3s %16s   %-s\n""" % ('', 'Method', '/', 'Basis', 'Rqd', 'Energy [Eh]',
                                                                'Variable')
    tables += table_delimit
    for job in JOBS_EXT:
        star = ''
        for mc in MODELCHEM:
            if (job['f_component'] == mc['f_component']) and (job['f_basis'] == mc['f_basis']):
                star = '*'
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (
            '', job['f_wfn'], '/', job['f_basis'] + " + options" * bool(job['f_options']), star, job['f_energy'],
            VARH[job['f_wfn']][job['f_component']])
    tables += table_delimit

    tables += """\n   ==> %s <==\n\n""" % ('Stages')
    tables += table_delimit
    tables += """     %6s %20s %1s %-27s %2s %16s   %-s\n""" % ('Stage', 'Method', '/', 'Basis', 'Wt', 'Energy [Eh]',
                                                                'Scheme')
    tables += table_delimit
    for stage in GRAND_NEED:
        tables += """     %6s %20s %1s %-27s %2d %16.8f   %-s\n""" % (stage['d_stage'], stage['d_wfn'], '/',
                                                                      stage['d_basis'], stage['d_coef'],
                                                                      stage['d_energy'], stage['d_scheme'].__name__)
    tables += table_delimit

    tables += """\n   ==> %s <==\n\n""" % ('CBS')
    tables += table_delimit
    tables += """     %6s %20s %1s %-27s %2s %16s   %-s\n""" % ('Stage', 'Method', '/', 'Basis', '', 'Energy [Eh]',
                                                                'Scheme')
    tables += table_delimit
    tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (
        GRAND_NEED[0]['d_stage'], GRAND_NEED[0]['d_wfn'], '/', GRAND_NEED[0]['d_basis'], '', GRAND_NEED[0]['d_energy'],
        GRAND_NEED[0]['d_scheme'].__name__)
    if len(metadata) > 1:
        if GRAND_NEED[1]['d_isdelta']:
            tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (
                GRAND_NEED[1]['d_stage'], GRAND_NEED[1]['d_wfn'], '/', GRAND_NEED[1]['d_basis'], '',
                GRAND_NEED[1]['d_energy'] - GRAND_NEED[2]['d_energy'], GRAND_NEED[1]['d_scheme'].__name__)
            dc = 3
        else:
            tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (
                GRAND_NEED[1]['d_stage'], GRAND_NEED[1]['d_wfn'], '/', GRAND_NEED[1]['d_basis'], '',
                GRAND_NEED[1]['d_energy'], GRAND_NEED[1]['d_scheme'].__name__)
            dc = 2
    if len(metadata) > 2:
        while dc < len(GRAND_NEED):
            if GRAND_NEED[dc]['d_isdelta']:
                deltaE_total = GRAND_NEED[dc]['d_energy'] - GRAND_NEED[dc + 1]['d_energy']
                tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (
                    GRAND_NEED[dc]['d_stage'], GRAND_NEED[dc]['d_wfn'] + ' - ' + GRAND_NEED[dc + 1]['d_wfn'], '/',
                    GRAND_NEED[dc]['d_basis'], '', deltaE_total, GRAND_NEED[dc]['d_scheme'].__name__)
                core.set_variable(f"CBS {GRAND_NEED[dc]['d_stage'].upper()} TOTAL ENERGY", deltaE_total)
                dc += 2
            else:
                tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (
                    GRAND_NEED[dc]['d_stage'], GRAND_NEED[dc]['d_wfn'], '/', GRAND_NEED[dc]['d_basis'], '',
                    GRAND_NEED[dc]['d_energy'], GRAND_NEED[dc]['d_scheme'].__name__)
                core.set_variable(f"CBS {GRAND_NEED[dc]['d_stage'].upper()} TOTAL ENERGY", GRAND_NEED[dc]['d_energy'])
                dc += 1

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
    basis = core.BasisSet.build(molecule, "ORBITAL", 'def2-svp')
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
        wfn.set_gradient(finalgradient)
        wfn.set_hessian(finalquantity)
        if finalquantity.rows(0) < 20:
            core.print_out('CURRENT HESSIAN')
            finalquantity.print_out()

    if return_wfn:
        return (finalquantity, wfn)
    else:
        return finalquantity


######### COMPUTE / ASSEMBLE
######### ASSEMBLE / REPORT


def _expand_scheme_orders(
    scheme: Callable, 
    basisname: List[str], 
    basiszeta: List[int], 
    wfnname: str, 
    options: dict,
    natom: int
) -> dict:
    """Check that the length of *basiszeta* array matches the implied degree of
    extrapolation in *scheme* name. Return a dictionary of same length as
    basiszeta, with *basisname* and *basiszeta* distributed therein.

    """
    Nxtpl = len(basiszeta)

    if int(scheme.__name__.split('_')[-1]) != Nxtpl:
        raise ValidationError("""Call to '%s' not valid with '%s' basis sets.""" % (scheme.__name__, len(basiszeta)))

    f_fields = ['f_wfn', 'f_component', 'f_basis', 'f_zeta', 'f_energy', 'f_gradient', 'f_hessian', 'f_options']
    NEED = {}
    for idx in range(Nxtpl):
        NEED[_lmh_labels[Nxtpl][idx]] = dict(
            zip(f_fields, [
                wfnname,
                wfnname,
                basisname[idx],
                basiszeta[idx],
                0.0,
                core.Matrix(natom, 3),
                core.Matrix(3 * natom, 3 * natom),
                options,
            ]))
    return NEED


def _contract_scheme_orders(needdict: dict, datakey: str = 'f_energy') -> dict:
    """Prepared named arguments for extrapolation functions by
    extracting zetas and values (which one determined by *datakey*) out
    of *needdict* and returning a dictionary whose keys are constructed
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


def _cbs_wrapper_methods(**kwargs: dict) -> List[str]:
    """ A helper function for the driver to enumerate methods used in the
    stages of a cbs calculation.

    Parameters
    ----------
    kwargs
        kwargs containing cbs specification either in the ``cbs_metadata``
        format, or in separate keywords (``scf_wfn``, ``corl_wfn`` etc.).

    Returns
    -------
    cbs_methods : List[str]
        List containing method name for each active stage.
    """

    cbs_methods = []
    if "cbs_metadata" in kwargs:
        for item in kwargs["cbs_metadata"]:
            cbs_methods.append(item.get("wfn"))
    else:
        cbs_method_kwargs = ['scf_wfn', 'corl_wfn', 'delta_wfn']
        cbs_method_kwargs += [f'delta{x}_wfn' for x in range(2, 6)]
        for method in cbs_method_kwargs:
            if method in kwargs:
                cbs_methods.append(kwargs[method])
    return cbs_methods


def _parse_cbs_gufunc_string(method_name: str) -> Tuple[List[str], List[int]]:
    """ A helper function that parses a ``"method/basis"`` input string
    into separate method and basis components. Also handles delta corrections.

    Parameters
    ----------
    method_name
        A ``"method/basis"`` style string defining the calculation.

    Returns
    -------
    (methods_list, basis_list) : Tuple[List[str], List[int]]
        Tuple in the ``(method_list, basis_list)`` format, where ``method_list``
        is the list of the component methods, and ``basis_list`` is the list of
        basis sets forming the extrapolation for each specified method.
        E.g. ``"mp2/cc-pv[tq]z+D:ccsd(t)/cc-pvtz"`` would return:
        ``(["mp2", "ccsd(t)"], ["cc-pv[tq]z", "cc-pvtz"])``.
    """

    method_name_list = re.split(r"""\+(?=\s*[Dd]:)""", method_name)
    if len(method_name_list) > 2:
        raise ValidationError(
            "CBS gufunc: Text parsing is only valid for a single delta, please use the CBS wrapper directly")

    method_list = []
    basis_list = []
    for num, method_str in enumerate(method_name_list):
        if (method_str.count("[") > 1) or (method_str.count("]") > 1):
            raise ValidationError(f"""CBS gufunc: Too many brackets given! {method_str}""")

        if method_str.count('/') != 1:
            raise ValidationError(f"""CBS gufunc: All methods must specify a basis with '/'. {method_str}""")

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


def _cbs_gufunc(
    func: Callable, 
    total_method_name: str, 
    **kwargs: dict
) -> Union[Tuple[float, core.Wavefunction], float]:
    """
    A text based parser of the CBS method string. Provided to handle "method/basis"
    specification of the requested calculations. Also handles "simple" (i.e.
    one-method and one-basis) calls.

    Parameters
    ----------
    func
        Function to be called (energy, gradient, frequency or cbs).
    total_method_name
        String in a ``"method/basis"`` syntax. Simple calls (e.g. ``"blyp/sto-3g"``) are
        bounced out of CBS. More complex calls (e.g. ``"mp2/cc-pv[tq]z"`` or
        ``"mp2/cc-pv[tq]z+D:ccsd(t)/cc-pvtz"``) are expanded by `_parse_cbs_gufunc_string()`
        and pushed through :py:func:`~psi4.cbs`.

    Returns
    -------
    ptype_value : Union[Tuple[float, core.Wavefunction], float]
        Float, or if ``return_wfn`` is specified, a tuple of ``(value, wavefunction)``.
    """

    # Catch kwarg issues for all methods
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    core.clean_variables()
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
        if core.get_option("SCF", "DF_INTS_IO") != "SAVE":
            core.clean()
        optstash.restore()

        if return_wfn:
            return (ptype_value, wfn)
        else:
            return ptype_value

    # Drop out for unsupported calls
    if ptype not in ["energy", "gradient", "hessian"]:
        raise ValidationError("%s: Cannot extrapolate or delta correct %s yet." % (ptype.title(), ptype))

    # Catch kwarg issues for CBS methods only
    user_dertype = kwargs.pop('dertype', None)
    cbs_verbose = kwargs.pop('cbs_verbose', False)

    # If we are not a single call, let CBS wrapper handle it!
    cbs_kwargs = {}
    cbs_kwargs['ptype'] = ptype
    cbs_kwargs['return_wfn'] = True
    cbs_kwargs['molecule'] = molecule
    cbs_kwargs['verbose'] = cbs_verbose

    if user_dertype != None:
        cbs_kwargs['dertype'] = user_dertype

    # Find method and basis
    metadata = []
    if method_list[0] in ['scf', 'hf', 'c4-scf', 'c4-hf']:
        stage = {}
        stage['wfn'] = method_list[0]
        stage['basis'] = basis_list[0]
        if 'scf_scheme' in kwargs:
            stage['scheme'] = kwargs.pop('scf_scheme')
    else:
        # _interpret_cbs_inputs will produce scf stage automatically
        stage = {}
        stage['wfn'] = method_list[0]
        stage['basis'] = basis_list[0]
        if 'corl_scheme' in kwargs:
            stage['scheme'] = kwargs.pop('corl_scheme')
    metadata.append(stage)

    # "method/basis" syntax only allows for one delta correction
    # via "method/basis+D:delta/basis". Maximum length of method_list is 2.
    if len(method_list) == 2:
        stage = {}
        stage['wfn'] = method_list[1]
        stage['basis'] = basis_list[1]
        if 'delta_scheme' in kwargs:
            stage['scheme'] = kwargs.pop('delta_scheme')
        metadata.append(stage)

    cbs_kwargs["cbs_metadata"] = metadata
    ptype_value, wfn = cbs(func, label, **cbs_kwargs)

    if return_wfn:
        return (ptype_value, wfn)
    else:
        return ptype_value
