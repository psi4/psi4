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

import psi4
import p4const
from math import sqrt, pi
from diatomic_fits import *

def anharmonicity(rvals, energies, mol = None):
    """Generates spectroscopic constants for a diatomic molecules.
       Fits a diatomic potential energy curve using either a 5 or 9 point Legendre fit, locates the minimum
       energy point, and then applies second order vibrational perturbation theory to obtain spectroscopic
       constants.  The r values provided must bracket the minimum energy point, or an error will result.

       A dictionary with the following keys, which correspond to spectroscopic constants, is returned:

       :type rvals: list
       :param rvals: The bond lengths (in Angstrom) for which energies are
           provided of length either 5 or 9 but must be the same length as
           the energies array

       :type energies: list
       :param energies: The energies (Eh) computed at the bond lengths in the rvals list

       :returns: (*dict*) Keys: "re", "r0", "we", "wexe", "nu", "ZPVE(harmonic)", "ZPVE(anharmonic)", "Be", "B0", "ae", "De"
                 corresponding to the spectroscopic constants in cm-1
    """

    angstrom_to_bohr = 1.0 / p4const.psi_bohr2angstroms
    angstrom_to_meter = 10e-10;

    if len(rvals) != len(energies):
        raise Exception("The number of energies must match the number of distances")

    npoints = len(rvals)

    if npoints != 5 and npoints != 9:
        raise Exception("Only 5- or 9-point fits are implemented right now")

    psi4.print_out("\n\nPerforming a %d-point fit\n" % npoints)

    psi4.print_out("\nOptimizing geometry based on current surface:\n\n");
    if (npoints == 5):
        optx = rvals[2]
    elif (npoints == 9):
        optx = rvals[4]

    # Molecule can be passed in be user. Look at the function definition above.
    if mol == None:
        mol = psi4.get_active_molecule()
    natoms = mol.natom()
    if natoms != 2:
        raise Exception("The current molecule must be a diatomic for this code to work!")
    m1 = mol.mass(0)
    m2 = mol.mass(1)

    maxit = 30
    thres = 1.0e-9
    for i in range(maxit):
        if (npoints == 5):
            grad= first_deriv_5pt(rvals, energies, optx)
            secd = second_deriv_5pt(rvals, energies, optx)
            energy = function_5pt(rvals, energies, optx)
        elif (npoints == 9):
            grad = first_deriv_9pt(rvals, energies, optx)
            secd = second_deriv_9pt(rvals, energies, optx)
            energy = function_9pt(rvals, energies, optx)
        psi4.print_out("       E = %20.14f, x = %14.7f, grad = %20.14f\n" % (energy, optx, grad))
        if abs(grad) < thres:
            break
        optx -= grad / secd;
    psi4.print_out(" Final E = %20.14f, x = %14.7f, grad = %20.14f\n" % (function_5pt(rvals, energies, optx), optx, grad));

    if optx < min(rvals):
        raise Exception("Minimum energy point is outside range of points provided.  Use a lower range of r values.")
    if optx > max(rvals):
        raise Exception("Minimum energy point is outside range of points provided.  Use a higher range of r values.")

    if (npoints == 5):
        energy = function_5pt(rvals, energies, optx)
        first = first_deriv_5pt(rvals, energies, optx)
        secd = second_deriv_5pt(rvals, energies, optx) * p4const.psi_hartree2aJ
        third = third_deriv_5pt(rvals, energies, optx) * p4const.psi_hartree2aJ
        fourth = fourth_deriv_5pt(rvals, energies, optx) * p4const.psi_hartree2aJ
    elif (npoints == 9):
        energy = function_9pt(rvals, energies, optx)
        first = first_deriv_9pt(rvals, energies, optx)
        secd = second_deriv_9pt(rvals, energies, optx) * p4const.psi_hartree2aJ
        third = third_deriv_9pt(rvals, energies, optx) * p4const.psi_hartree2aJ
        fourth = fourth_deriv_9pt(rvals, energies, optx) * p4const.psi_hartree2aJ

    psi4.print_out("\nEquilibrium Energy %20.14f Hartrees\n" % energy)
    psi4.print_out("Gradient           %20.14f\n" % first)
    psi4.print_out("Quadratic Force Constant %14.7f MDYNE/A\n" % secd)
    psi4.print_out("Cubic Force Constant     %14.7f MDYNE/A**2\n" % third)
    psi4.print_out("Quartic Force Constant   %14.7f MDYNE/A**3\n" % fourth)

    hbar = p4const.psi_h / (2.0 * pi)
    mu = ((m1*m2)/(m1+m2))*p4const.psi_amu2kg
    we = 5.3088375e-11*sqrt(secd/mu)
    wexe = (1.2415491e-6)*(we/secd)**2 * ((5.0*third*third)/(3.0*secd)-fourth)

    # Rotational constant: Be
    I = ((m1*m2)/(m1+m2)) * p4const.psi_amu2kg * (optx * angstrom_to_meter)**2
    B = p4const.psi_h / (8.0 * pi**2 * p4const.psi_c * I)

    # alpha_e and quartic centrifugal distortion constant
    ae = -(6.0 * B**2 / we) * ((1.05052209e-3*we*third)/(sqrt(B * secd**3))+1.0)
    de = 4.0*B**3 / we**2

    # B0 and r0 (plus re check using Be)
    B0 = B - ae / 2.0
    r0 = sqrt(p4const.psi_h / (8.0 * pi**2 * mu * p4const.psi_c * B0))
    recheck = sqrt(p4const.psi_h / (8.0 * pi**2 * mu * p4const.psi_c * B))
    r0 /= angstrom_to_meter;
    recheck /= angstrom_to_meter;

    # Fundamental frequency nu
    nu = we - 2.0 * wexe;
    zpve_nu = 0.5 * we - 0.25 * wexe;

    psi4.print_out("\nre     = %10.6f A  check: %10.6f\n" % (optx, recheck))
    psi4.print_out("r0       = %10.6f A\n" % r0)
    psi4.print_out("we       = %10.4f cm-1\n" % we)
    psi4.print_out("wexe     = %10.4f cm-1\n" % wexe)
    psi4.print_out("nu       = %10.4f cm-1\n" % nu)
    psi4.print_out("ZPVE(nu) = %10.4f cm-1\n" % zpve_nu)
    psi4.print_out("Be       = %10.4f cm-1\n" % B)
    psi4.print_out("B0       = %10.4f cm-1\n" % B0)
    psi4.print_out("ae       = %10.4f cm-1\n" % ae)
    psi4.print_out("De       = %10.7f cm-1\n" % de)
    results = {
               "re"               :  optx,
               "r0"               :  r0,
               "we"               :  we,
               "wexe"             :  wexe,
               "nu"               :  nu,
               "ZPVE(harmonic)"   :  zpve_nu,
               "ZPVE(anharmonic)" :  zpve_nu,
               "Be"               :  B,
               "B0"               :  B0,
               "ae"               :  ae,
               "De"               :  de
              }
    return results

