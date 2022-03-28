#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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
import logging

import numpy as np

from psi4 import core
from psi4.driver import nppp
from psi4.driver.p4util.exceptions import ValidationError
from psi4.driver.aliases import sherrill_gold_standard, allen_focal_point

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

_zeta_val2sym = {k + 2: v for k, v in enumerate('dtq5678')}


def xtpl_highest_1(functionname, zHI, valueHI, verbose=1, **kwargs):
    r"""Scheme for total or correlation energies with a single basis or the highest
    zeta-level among an array of bases. Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname : string
        Name of the CBS component (e.g., 'mp2') used in summary printing.
    zHI : int
        Zeta-level, only used for printing.
    valueHI : float or ndarray
        Energy, gradient, or Hessian value at the basis set.
    verbose : int, optional
        Controls volume of printing.

    Returns
    -------
    float or ndarray
        Returns :math:`E_{total}^{\infty}` which is equal to valueHI.
        Eponymous function applied to input zetas and values; type from `valueHI`.

    Notes
    -----
    .. math:: E_{total}^X = E_{total}^{\infty}

    Examples
    --------
    >>> # [1] Fancy way to get HF/cc-pCVQZ
    >>> psi4.energy('cbs', scf_wfn='hf', scf_basis='cc-pcvqz', scf_scheme='xtpl_highest_1')

    """
    if isinstance(valueHI, float):

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> {} <==\n\n""".format(functionname.upper())
            cbsscheme += """   HI-zeta ({}) Energy:               {:16.12f}\n""".format(zHI, valueHI)

            core.print_out(cbsscheme)
            logger.debug(cbsscheme)

        return valueHI

    elif isinstance(valueHI, np.ndarray):

        if verbose > 2:
            cbsscheme = f"""\n   HI-zeta ({zHI}) Total Energy\n"""
            cbsscheme += nppp(valueHI)

        return valueHI


def scf_xtpl_helgaker_2(functionname, zLO, valueLO, zHI, valueHI, verbose=1, alpha=None):
    r"""Extrapolation scheme using exponential form for reference energies with two adjacent
    zeta-level bases. Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname : string
        Name of the CBS component (e.g., 'HF') used in summary printing.
    zLO : int
        Zeta number of the smaller basis set in 2-point extrapolation.
    valueLO : float or ndarray
        Energy, gradient, or Hessian value at the smaller basis set in 2-point.
        extrapolation.
    zHI : int
        Zeta number of the larger basis set in 2-point extrapolation.
        Must be `zLO + 1`.
    valueHI : float or numpy.ndarray
        Energy, gradient, or Hessian value at the larger basis set in 2-point
        extrapolation.
    verbose : int, optional
        Controls volume of printing.
    alpha : float, optional
        Fitted 2-point parameter. Overrides the default :math:`\alpha = 1.63`

    Returns
    -------
    float or ndarray
        Eponymous function applied to input zetas and values; type from `valueLO`.

    Notes
    -----
    The extrapolation is calculated according to [1]_:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha X}, \alpha = 1.63`

    References
    ----------

    .. [1] Halkier, Helgaker, Jorgensen, Klopper, & Olsen, Chem. Phys. Lett. 302 (1999) 437-446,
       DOI: 10.1016/S0009-2614(99)00179-7

    Examples
    --------
    >>> # [1] Hartree-Fock extrapolation
    >>> psi4.energy('cbs', scf_wfn='hf', scf_basis='cc-pV[DT]Z', scf_scheme='scf_xtpl_helgaker_2')

    """

    if type(valueLO) != type(valueHI):
        raise ValidationError(
            f"scf_xtpl_helgaker_2: Inputs must be of the same datatype! ({type(valueLO)}, {type(valueHI)})")

    if alpha is None:
        alpha = 1.63

    beta_division = 1 / (math.exp(-1 * alpha * zLO) * (math.exp(-1 * alpha) - 1))
    beta_mult = math.exp(-1 * alpha * zHI)

    if isinstance(valueLO, float):
        beta = (valueHI - valueLO) * beta_division
        value = valueHI - beta * beta_mult

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Helgaker 2-point exponential SCF extrapolation for method: %s <==\n\n""" % (
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
            logger.debug(cbsscheme)

        return value

    elif isinstance(valueLO, np.ndarray):
        beta = (valueHI - valueLO) * beta_division
        value = valueHI - beta * beta_mult
        #beta.name = 'Helgaker SCF (%s, %s) beta' % (zLO, zHI)
        #value.name = 'Helgaker SCF (%s, %s) data' % (zLO, zHI)

        if verbose > 2:
            cbsscheme = f"""\n   ==> Helgaker 2-point exponential SCF extrapolation for method: {functionname.upper()} <==\n"""
            cbsscheme += f"""\n   LO-zeta ({zLO}) Data\n"""
            cbsscheme += nppp(valueLO)
            cbsscheme += f"""\n   HI-zeta ({zHI}) Data\n"""
            cbsscheme += nppp(valueHI)

            cbsscheme += """   Alpha (exponent) Value:          %16.8f\n""" % (alpha)
            cbsscheme += f"""\n   Beta Data\n"""
            cbsscheme += nppp(beta)
            cbsscheme += f"""\n   Extrapolated Data\n"""
            cbsscheme += nppp(value)
            core.print_out(cbsscheme)
            logger.debug(cbsscheme)

        return value

    else:
        raise ValidationError("scf_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))


def scf_xtpl_truhlar_2(functionname, zLO, valueLO, zHI, valueHI, verbose=1, alpha=None):
    r"""Extrapolation scheme using power form for reference energies with two adjacent
    zeta-level bases. Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname : string
        Name of the CBS component.
    zLO : int
        Lower zeta level.
    valueLO : float
        Lower value used for extrapolation.
    zHI : int
        Higher zeta level. Should be equal to zLO + 1.
    valueHI : float
        Higher value used for extrapolation.
    alpha : float, optional
        Overrides the default :math:`\alpha = 3.4`

    Returns
    -------
    float
        Returns :math:`E_{total}^{\infty}`, see below.

    Notes
    -----
    The extrapolation is calculated according to [2]_:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta X^{-\alpha}, \alpha = 3.4`

    References
    ----------

    .. [2] Truhlar, Chem. Phys. Lett. 294 (1998) 45-48,
       DOI: 10.1016/S0009-2614(98)00866-5

    """

    if type(valueLO) != type(valueHI):
        raise ValidationError(
            f"scf_xtpl_truhlar_2: Inputs must be of the same datatype! ({type(valueLO)}, {type(valueHI)})")

    if alpha is None:
        alpha = 3.40

    beta_division = 1 / (zHI**(-1 * alpha) - zLO**(-1 * alpha))
    beta_mult = zHI**(-1 * alpha)

    if isinstance(valueLO, float):
        beta = (valueHI - valueLO) * beta_division
        value = valueHI - beta * beta_mult

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Truhlar 2-point power form SCF extrapolation for method: %s <==\n\n""" % (
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
            logger.debug(cbsscheme)

        return value

    elif isinstance(valueLO, (core.Matrix, core.Vector)):
        beta = valueHI.clone()
        beta.name = 'Truhlar SCF (%s, %s) beta' % (zLO, zHI)
        beta.subtract(valueLO)
        beta.scale(beta_division)
        beta.scale(beta_mult)

        value = valueHI.clone()
        value.subtract(beta)
        value.name = 'Truhlar SCF (%s, %s) data' % (zLO, zHI)

        if verbose > 2:
            cbsscheme = f"""\n   ==> Truhlar 2-point power SCF extrapolation for method: {functionname.upper()} <==\n"""
            cbsscheme += f"""\n   LO-zeta ({zLO}) Data\n"""
            cbsscheme += nppp(valueLO)
            cbsscheme += f"""\n   HI-zeta ({zHI}) Data\n"""
            cbsscheme += nppp(valueHI)

            cbsscheme += """   Alpha (exponent) Value:          %16.8f\n""" % (alpha)
            cbsscheme += f"""\n   Beta Data\n"""
            cbsscheme += nppp(beta)
            cbsscheme += f"""\n   Extrapolated Data\n"""
            cbsscheme += nppp(value)
            core.print_out(cbsscheme)
            logger.debug(cbsscheme)

        return value

    else:
        raise ValidationError("scf_xtpl_truhlar_2: datatype is not recognized '%s'." % type(valueLO))


def scf_xtpl_karton_2(functionname, zLO, valueLO, zHI, valueHI, verbose=1, alpha=None):
    r"""Extrapolation scheme using root-power form for reference energies with two adjacent
    zeta-level bases. Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname : str
        Name of the CBS component.
    zLO : int
        Lower zeta level.
    valueLO : float
        Lower value used for extrapolation.
    zHI : int
        Higher zeta level. Should be equal to zLO + 1.
    valueHI : float
        Higher value used for extrapolation.
    alpha : float, optional
        Overrides the default :math:`\alpha = 6.3`

    Returns
    -------
    float
        Returns :math:`E_{total}^{\infty}`, see below.

    Notes
    -----
    The extrapolation is calculated according to [3]_:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha\sqrt{X}}, \alpha = 6.3`

    References
    ----------

    .. [3] Karton, Martin, Theor. Chem. Acc. 115 (2006) 330-333,
       DOI: 10.1007/s00214-005-0028-6

    """

    if type(valueLO) != type(valueHI):
        raise ValidationError(
            f"scf_xtpl_karton_2: Inputs must be of the same datatype! ({type(valueLO)}, {type(valueHI)})")

    if alpha is None:
        alpha = 6.30

    beta_division = 1 / (math.exp(-1 * alpha) * (math.exp(math.sqrt(zHI)) - math.exp(math.sqrt(zLO))))
    beta_mult = math.exp(-1 * alpha * math.sqrt(zHI))

    if isinstance(valueLO, float):
        beta = (valueHI - valueLO) * beta_division
        value = valueHI - beta * beta_mult

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Karton 2-point power form SCF extrapolation for method: %s <==\n\n""" % (
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
            logger.debug(cbsscheme)

        return value

    elif isinstance(valueLO, (core.Matrix, core.Vector)):
        beta = valueHI.clone()
        beta.name = 'Karton SCF (%s, %s) beta' % (zLO, zHI)
        beta.subtract(valueLO)
        beta.scale(beta_division)
        beta.scale(beta_mult)

        value = valueHI.clone()
        value.subtract(beta)
        value.name = 'Karton SCF (%s, %s) data' % (zLO, zHI)

        if verbose > 2:
            cbsscheme = f"""\n   ==> Karton 2-point power SCF extrapolation for method: {functionname.upper()} <==\n"""
            cbsscheme += f"""\n   LO-zeta ({zLO}) Data\n"""
            cbsscheme += nppp(valueLO)
            cbsscheme += f"""\n   HI-zeta ({zHI}) Data\n"""
            cbsscheme += nppp(valueHI)

            cbsscheme += """   Alpha (exponent) Value:          %16.8f\n""" % (alpha)
            cbsscheme += f"""\n   Beta Data\n"""
            cbsscheme += nppp(beta)
            cbsscheme += f"""\n   Extrapolated Data\n"""
            cbsscheme += nppp(value)
            core.print_out(cbsscheme)
            logger.debug(cbsscheme)

        return value

    else:
        raise ValidationError("scf_xtpl_Karton_2: datatype is not recognized '%s'." % type(valueLO))


def scf_xtpl_helgaker_3(functionname, zLO, valueLO, zMD, valueMD, zHI, valueHI, verbose=1, alpha=None):
    r"""Extrapolation scheme for reference energies with three adjacent zeta-level bases.
    Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    functionname : str
        Name of the CBS component (e.g., 'HF') used in summary printing.
    zLO : int
        Zeta number of the smaller basis set in 3-point extrapolation.
    valueLO : float or ndarray
        Energy, gradient, or Hessian value at the smaller basis set in 3-point
        extrapolation.
    zMD : int
        Zeta number of the medium basis set in 3-point extrapolation.
        Must be `zLO + 1`.
    valueMD : float or ndarray
        Energy, gradient, or Hessian value at the medium basis set in 3-point
        extrapolation.
    zHI : int
        Zeta number of the larger basis set in 3-point extrapolation.
        Must be `zLO + 2`.
    valueHI : float or ndarray
        Energy, gradient, or Hessian value at the larger basis set in 3-point
        extrapolation.
    verbose : int, optional
        Controls volume of printing.
    alpha : float, optional
        Not used.

    Returns
    -------
    float or ndarray
        Eponymous function applied to input zetas and values; type from `valueLO`.

    Notes
    -----
    The extrapolation is calculated according to [4]_:
    :math:`E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha X}, \alpha = 3.0`

    References
    ----------

    .. [4] Halkier, Helgaker, Jorgensen, Klopper, & Olsen, Chem. Phys. Lett. 302 (1999) 437-446,
       DOI: 10.1016/S0009-2614(99)00179-7

    Examples
    --------
    >>> # [1] Hartree-Fock extrapolation
    >>> psi4.energy('cbs', scf_wfn='hf', scf_basis='cc-pV[DTQ]Z', scf_scheme='scf_xtpl_helgaker_3')

    """

    if (type(valueLO) != type(valueMD)) or (type(valueMD) != type(valueHI)):
        raise ValidationError(
            f"scf_xtpl_helgaker_3: Inputs must be of the same datatype! ({type(valueLO)}, {type(valueMD)}, {type(valueHI)})"
        )

    if isinstance(valueLO, float):

        ratio = (valueHI - valueMD) / (valueMD - valueLO)
        alpha = -1 * math.log(ratio)
        beta = (valueHI - valueMD) / (math.exp(-1 * alpha * zMD) * (ratio - 1))
        value = valueHI - beta * math.exp(-1 * alpha * zHI)

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Helgaker 3-point SCF extrapolation for method: %s <==\n\n""" % (
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
            logger.debug(cbsscheme)

        return value

    elif isinstance(valueLO, np.ndarray):
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

        if verbose > 2:
            cbsscheme = f"""\n   ==> Helgaker 3-point correlated extrapolation for method: {functionname.upper()} <==\n"""
            cbsscheme += f"""\n   LO-zeta ({zLO}) Data\n"""
            cbsscheme += nppp(valueLO)
            cbsscheme += f"""\n   MD-zeta ({zMD}) Data\n"""
            cbsscheme += nppp(valueMD)
            cbsscheme += f"""\n   HI-zeta ({zHI}) Data\n"""
            cbsscheme += nppp(valueHI)

            cbsscheme += f"""\n   Alpha Data\n"""
            cbsscheme += nppp(alpha)
            cbsscheme += f"""\n   Beta Data\n"""
            cbsscheme += nppp(beta)
            cbsscheme += f"""\n   Extrapolated Data\n"""
            cbsscheme += nppp(np_value)
            core.print_out(cbsscheme)
            logger.debug(cbsscheme)

        return np_value

    else:
        raise ValidationError(f"scf_xtpl_helgaker_3: datatype is not recognized '{type(valueLO)}'.")


def corl_xtpl_helgaker_2(functionname, zLO, valueLO, zHI, valueHI, verbose=1, alpha=None):
    r"""Extrapolation scheme for correlation energies with two adjacent zeta-level bases.
    Used by :py:func:`~psi4.cbs`.

    Parameters
    ----------
    mtdname : str
        Name of the CBS components (e.g., 'MP2') used in summary printing.
    zLO : int
        Zeta number of the smaller basis set in 2-point extrapolation.
    valueLO : float or ndarray
        Energy, gradient, or Hessian value at the smaller basis set in 2-point
        extrapolation.
    zHI : int
        Zeta number of the larger basis set in 2-point extrapolation.
        Must be `zLO + 1`.
    valueHI : float or numpy.ndarray
        Energy, gradient, or Hessian value at the larger basis set in 2-point
        extrapolation.
    verbose : int, optional
        Controls volume of printing.
    alpha : float, optional
        Overrides the default :math:`\alpha = 3.0`

    Returns
    -------
    float or numpy.ndarray
        Eponymous function applied to input zetas and values; type from `valueLO`.

    Notes
    -----
    The extrapolation is calculated according to [5]_:
    :math:`E_{corl}^X = E_{corl}^{\infty} + \beta X^{-alpha}`

    References
    ----------

    .. [5] Halkier, Helgaker, Jorgensen, Klopper, Koch, Olsen, & Wilson,
       Chem. Phys. Lett. 286 (1998) 243-252,
       DOI: 10.1016/S0009-2614(99)00179-7

    Examples
    --------
    >>> # [1] CISD extrapolation
    >>> energy('cbs', corl_wfn='cisd', corl_basis='cc-pV[DT]Z', corl_scheme='corl_xtpl_helgaker_2')

    """
    if type(valueLO) != type(valueHI):
        raise ValidationError(
            f"corl_xtpl_helgaker_2: Inputs must be of the same datatype! ({type(valueLO)}, {type(valueHI)})")

    if alpha is None:
        alpha = 3.0

    if isinstance(valueLO, float):
        value = (valueHI * zHI**alpha - valueLO * zLO**alpha) / (zHI**alpha - zLO**alpha)
        beta = (valueHI - valueLO) / (zHI**(-alpha) - zLO**(-alpha))  # lgtm[py/unused-local-variable]

        final = value
        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = f"""\n\n   ==> Helgaker 2-point correlated extrapolation for method: {functionname.upper()} <==\n\n"""
            cbsscheme += """   LO-zeta (%s) Energy:               % 16.12f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Alpha (exponent) Value:           % 16.12f\n""" % alpha
            cbsscheme += """   Extrapolated Energy:              % 16.12f\n\n""" % value

            name_str = "%s/(%s,%s)" % (functionname.upper(), _zeta_val2sym[zLO].upper(), _zeta_val2sym[zHI].upper())
            cbsscheme += """   @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (19 - len(name_str))
            cbsscheme += """% 16.12f\n\n""" % final
            core.print_out(cbsscheme)
            logger.info(cbsscheme)

        return final

    elif isinstance(valueLO, np.ndarray):
        value = (valueHI * zHI**alpha - valueLO * zLO**alpha) / (zHI**alpha - zLO**alpha)

        if verbose > 2:
            cbsscheme = f"""\n   ==> Helgaker 2-point correlated extrapolation for method: {functionname.upper()} <==\n"""
            cbsscheme += f"""\n   LO-zeta ({zLO}) Data\n"""
            cbsscheme += nppp(valueLO)
            cbsscheme += f"""\n   HI-zeta ({zHI}) Data\n"""
            cbsscheme += nppp(valueHI)
            cbsscheme += f"""\n   Extrapolated Data\n"""
            cbsscheme += nppp(value)
            core.print_out(cbsscheme)
            logger.debug(cbsscheme)

        return value

    else:
        raise ValidationError("corl_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))


xtpl_procedures = {
    'xtpl_highest_1': xtpl_highest_1,
    'scf_xtpl_helgaker_2': scf_xtpl_helgaker_2,
    'scf_xtpl_truhlar_2': scf_xtpl_truhlar_2,
    'scf_xtpl_karton_2': scf_xtpl_karton_2,
    'scf_xtpl_helgaker_3': scf_xtpl_helgaker_3,
    'corl_xtpl_helgaker_2': corl_xtpl_helgaker_2,
}

composite_procedures = {
    'sherrill_gold_standard': sherrill_gold_standard,
    'allen_focal_point': allen_focal_point,
}


def register_xtpl_function(func):
    """Enable user-defined extrapolation `func`."""

    if func.__name__.split('_')[-1].isdigit():
        xtpl_procedures[func.__name__] = func
    else:
        raise ValidationError("Extrapolation function names follow <scf|corl>_xtpl_<scientist>_<#basis>")


def register_composite_function(func):
    """Enable user-defined alias to composite method `func`."""

    composite_procedures[func.__name__] = func
