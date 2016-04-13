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

from p4util.exceptions import *
import psi4


def scf_xtpl_helgaker_2(funcitonname, zLO, valueLO, zHI, valueHI, verbose=True, alpha=1.63):
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
            cbsscheme += """\n   ==> %s <==\n\n""" % (functionname)
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

#def scf_xtpl_helgaker_3(**largs):
#    r"""Extrapolation scheme for reference energies with three adjacent zeta-level bases.
#    Used by :py:func:`~wrappers.complete_basis_set`.
#
#    .. math:: E_{total}^X = E_{total}^{\infty} + \beta e^{-\alpha X}
#
#    """
#    energypiece = 0.0
#    functionname = sys._getframe().f_code.co_name
#    f_fields = ['f_wfn', 'f_portion', 'f_basis', 'f_zeta', 'f_energy']
#    [mode, NEED, wfnname, BSET, ZSET] = validate_scheme_args(functionname, **largs)
#
#    if (mode == 'requisition'):
#
#        # Impose restrictions on zeta sequence
#        if (len(ZSET) != 3):
#            raise ValidationError('Call to \'%s\' not valid with \'%s\' basis sets.' % (functionname, len(ZSET)))
#
#        # Return array that logs the requisite jobs
#        NEED = {'HI': dict(zip(f_fields, [wfnname, 'tot', BSET[2], ZSET[2], 0.0])),
#                'MD': dict(zip(f_fields, [wfnname, 'tot', BSET[1], ZSET[1], 0.0])),
#                'LO': dict(zip(f_fields, [wfnname, 'tot', BSET[0], ZSET[0], 0.0]))}
#
#        return NEED
#
#    elif (mode == 'evaluate'):
#
#        # Extract required energies and zeta integers from array
#        eHI = NEED['HI']['f_energy']
#        eMD = NEED['MD']['f_energy']
#        eLO = NEED['LO']['f_energy']
#        zHI = NEED['HI']['f_zeta']
#        zMD = NEED['MD']['f_zeta']
#        zLO = NEED['LO']['f_zeta']
#
#        # Compute extrapolated energy
#        ratio = (eHI - eMD) / (eMD - eLO)
#        alpha = -1 * math.log(ratio)
#        beta = (eHI - eMD) / (math.exp(-1 * alpha * zMD) * (ratio - 1))
#        energypiece = eHI - beta * math.exp(-1 * alpha * zHI)
#
#        # Output string with extrapolation parameters
#        cbsscheme = ''
#        cbsscheme += """\n   ==> %s <==\n\n""" % (functionname)
#        cbsscheme += """   LO-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zLO), eLO)
#        cbsscheme += """   MD-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zMD), eMD)
#        cbsscheme += """   HI-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zHI), eHI)
#        cbsscheme += """   Extrapolated Correlation Energy: %16.8f\n""" % (energypiece)
#        cbsscheme += """   Alpha (exponent) Value:          %16.8f\n""" % (alpha)
#        cbsscheme += """   Beta (coefficient) Value:        %16.8f\n""" % (beta)
#        psi4.print_out(cbsscheme)
#
#

    # Compute extrapolated energy

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
        value = (valueHI * zHI ** 3 - valueLO * zLO ** 3) / (zHI ** 3 - zLO ** 3)


        beta = (valueHI - valueLO) / (zHI ** (-3) - zLO ** (-3))
        beta = valueHI.clone()
        beta.set_name('Helgaker SCF (%s, %s) beta' % (zLO, zHI))
        beta.subtract(valueLO)
        beta.scale(beta_division)
        beta.scale(beta_mult)

        value = valueHI.clone()
        value.subtract(beta)        
        value.set_name('Helgaker SCF (%s, %s) data' % (zLO, zHI))

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

