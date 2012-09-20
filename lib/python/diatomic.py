import PsiMod
from physconst import psi_bohr2angstroms, psi_hartree2aJ, psi_amu2kg, psi_h, psi_c
from math import sqrt, pi
from diatomic_fits import *

def diatomic_anharmonicity(rvals, energies):
    """Generates spectroscopic constants for a diatomic molecules.  Input parameters are

       rvals: a list containing the bond lengths (in Angstrom) for which energies are provided
       energies: a list containing the energies (Eh) at the bond lengths in the rvals list

       The above lists may be either 5 or 9 elements long, and the r values chosen must bracket
       the minimum energy point, or an error will result. A dictionary with the following keys,
       which correspond to spectroscopic constants, is returned:
               "re"              
               "r0"              
               "we"              
               "wexe"            
               "nu"              
               "ZPVE(harmonic)"  
               "ZPVE(anharmonic)"
               "Be"              
               "B0"              
               "ae"              
               "De"              
    """

    angstrom_to_bohr = 1.0 / psi_bohr2angstroms
    angstrom_to_meter = 10e-10;
   
    if len(rvals) != len(energies):
        raise Exception("The number of energies must match the number of distances")
   
    npoints = len(rvals)

    if npoints != 5 and npoints != 9:
        raise Exception("Only 5- or 9-point fits are implemented right now")

    PsiMod.print_out("\n\nPerforming a %d-point fit\n" % npoints)
    
    PsiMod.print_out("\nOptimizing geometry based on current surface:\n\n");
    if (npoints == 5):
        optx = rvals[2]
    elif (npoints == 9):
        optx = rvals[4]

    mol = PsiMod.get_active_molecule()
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
        PsiMod.print_out("       E = %20.14f, x = %14.7f, grad = %20.14f\n" % (energy, optx, grad))
        if abs(grad) < thres:
            break
        optx -= grad / secd;
    PsiMod.print_out(" Final E = %20.14f, x = %14.7f, grad = %20.14f\n" % (function_5pt(rvals, energies, optx), optx, grad));

    if optx < min(rvals):
        raise Exception("Minimum energy point is outside range of points provided.  Use a lower range of r values.")
    if optx > max(rvals):
        raise Exception("Minimum energy point is outside range of points provided.  Use a higher range of r values.")
    
    if (npoints == 5):
        energy = function_5pt(rvals, energies, optx)
        first = first_deriv_5pt(rvals, energies, optx)
        secd = second_deriv_5pt(rvals, energies, optx) * psi_hartree2aJ
        third = third_deriv_5pt(rvals, energies, optx) * psi_hartree2aJ
        fourth = fourth_deriv_5pt(rvals, energies, optx) * psi_hartree2aJ
    elif (npoints == 9):
        energy = function_9pt(rvals, energies, optx)
        first = first_deriv_9pt(rvals, energies, optx)
        secd = second_deriv_9pt(rvals, energies, optx) * psi_hartree2aJ
        third = third_deriv_9pt(rvals, energies, optx) * psi_hartree2aJ
        fourth = fourth_deriv_9pt(rvals, energies, optx) * psi_hartree2aJ

    PsiMod.print_out("\nEquilibrium Energy %20.14f Hartrees\n" % energy)
    PsiMod.print_out("Gradient           %20.14f\n" % first)
    PsiMod.print_out("Quadratic Force Constant %14.7f MDYNE/A\n" % secd)
    PsiMod.print_out("Cubic Force Constant     %14.7f MDYNE/A**2\n" % third)
    PsiMod.print_out("Quartic Force Constant   %14.7f MDYNE/A**3\n" % fourth)
    
    hbar = psi_h / (2.0 * pi)
    mu = ((m1*m2)/(m1+m2))*psi_amu2kg
    we = 5.3088375e-11*sqrt(secd/mu)
    wexe = (1.2415491e-6)*(we/secd)**2 * ((5.0*third*third)/(3.0*secd)-fourth)
    
    # Rotational constant: Be
    I = ((m1*m2)/(m1+m2)) * psi_amu2kg * (optx * angstrom_to_meter)**2
    B = psi_h / (8.0 * pi**2 * psi_c * I)
    
    # alpha_e and quartic centrifugal distortion constant
    ae = -(6.0 * B**2 / we) * ((1.05052209e-3*we*third)/(sqrt(B * secd**3))+1.0)
    de = 4.0*B**3 / we**2

    # B0 and r0 (plus re check using Be)
    B0 = B - ae / 2.0
    r0 = sqrt(psi_h / (8.0 * pi**2 * mu * psi_c * B0))
    recheck = sqrt(psi_h / (8.0 * pi**2 * mu * psi_c * B))
    r0 /= angstrom_to_meter; 
    recheck /= angstrom_to_meter;
    
    # Fundamental frequency nu
    nu = we - 2.0 * wexe;
    zpve_nu = 0.5 * we - 0.25 * wexe;

    PsiMod.print_out("\nre     = %10.6f A  check: %10.6f\n" % (optx, recheck))
    PsiMod.print_out("r0       = %10.6f A\n" % r0)
    PsiMod.print_out("we       = %10.4f cm-1\n" % we)
    PsiMod.print_out("wexe     = %10.4f cm-1\n" % wexe)
    PsiMod.print_out("nu       = %10.4f cm-1\n" % nu)
    PsiMod.print_out("ZPVE(nu) = %10.4f cm-1\n" % zpve_nu)
    PsiMod.print_out("Be       = %10.4f cm-1\n" % B)
    PsiMod.print_out("B0       = %10.4f cm-1\n" % B0)
    PsiMod.print_out("ae       = %10.4f cm-1\n" % ae)
    PsiMod.print_out("De       = %10.7f cm-1\n" % de)
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

