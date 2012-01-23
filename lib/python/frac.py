import PsiMod
import os
import input
import math
from molecule import * 
from driver import * 

def frac_traverse(mol, **kwargs):

    # The molecule is required, and should be the neutral species
    mol.update_geometry()
    activate(mol)
    charge0 = mol.molecular_charge()
    mult0   = mol.multiplicity() 

    chargep = charge0 + 1
    chargem = charge0 - 1
    
    # By default, the multiplicity of the cation/anion are mult0 + 1
    # These are overridden with the cation_mult and anioin_mult kwargs
    multp = mult0 + 1
    multm = mult0 + 1
    if kwargs.has_key('cation_mult'):
        multp = kwargs['cation_mult'] 
    if kwargs.has_key('anion_mult'):
        multm = kwargs['anion_mult'] 
    
    # By default, we start the frac procedure on the 25th iteration 
    # when not reading a previous guess
    frac_start = 25
    if kwargs.has_key('frac_start'):
        frac_start = kwargs['frac_start']

    # By default, we occupy by tenths of electrons
    LUMO_occs = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
    HOMO_occs  = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
    if kwargs.has_key('HOMO_occs'):
        HOMO_occs = kwargs['HOMO_occs']    
    if kwargs.has_key('LUMO_occs'):
        LUMO_occs = kwargs['LUMO_occs']    

    # By default, HOMO and LUMO are both in alpha
    Z = 0;
    for A in range(mol.natom()):
        Z += mol.Z(A)
    Z -= charge0
    if (Z%2):
        HOMO = Z/2+1
    else:
        HOMO = Z/2
    LUMO = HOMO+1
    if kwargs.has_key('HOMO'):
        HOMO = kwargs['HOMO']
    if kwargs.has_key('LUMO'):
        LUMO = kwargs['LUMO']

    # By default, DIIS in FRAC (1.0 occupation is always DIIS'd)
    frac_diis = True
    if kwargs.has_key('frac_diis'):
        frac_diis = kwargs['frac_diis']

    # By default, use the neutral orbitals as a guess for the anion
    neutral_guess = True
    if kwargs.has_key('neutral_guess'):
        neutral_guess = kwargs['neutral_guess']

    # => Traverse <= #
    occs = []
    energies = []

    # => Run the neutral for its orbitals, if requested <= #

    if (neutral_guess):
        energy('scf')
        old_guess = PsiMod.get_global_option("GUESS")
        PsiMod.set_global_option("GUESS", "READ")

    # => Run the anion first <= #

    mol.set_molecular_charge(chargem)
    mol.set_multiplicity(multm)

    PsiMod.set_global_option("FRAC_START", frac_start)
    PsiMod.set_global_option("FRAC_RENORMALIZE", True)
    PsiMod.set_global_option("FRAC_LOAD", False)

    for occ in LUMO_occs:
       
        PsiMod.set_global_option("FRAC_OCC", [LUMO])
        PsiMod.set_global_option("FRAC_VAL", [occ]) 

        E = energy('scf')

        occs.append(occ)
        energies.append(E)

        PsiMod.set_global_option("FRAC_START", 2)
        PsiMod.set_global_option("FRAC_LOAD", True)
        PsiMod.set_global_option("FRAC_DIIS", frac_diis)

    # Reset the old guess    
    if (neutral_guess):
        PsiMod.set_global_option("GUESS", old_guess)

    # => Run the neutral next <= #

    mol.set_molecular_charge(charge0)
    mol.set_multiplicity(mult0)

    PsiMod.set_global_option("FRAC_START", frac_start)
    PsiMod.set_global_option("FRAC_RENORMALIZE", True)
    PsiMod.set_global_option("FRAC_LOAD", False)

    for occ in LUMO_occs:
       
        PsiMod.set_global_option("FRAC_OCC", [HOMO])
        PsiMod.set_global_option("FRAC_VAL", [occ]) 

        E = energy('scf')

        occs.append(occ - 1.0)
        energies.append(E)

        PsiMod.set_global_option("FRAC_START", 2)
        PsiMod.set_global_option("FRAC_LOAD", True)
        PsiMod.set_global_option("FRAC_DIIS", frac_diis)

    # => Print the results out <= #
    E = {}
    PsiMod.print_out('\n    ==> Fractional Occupation Traverse Results <==\n\n')
    PsiMod.print_out('\t%-11s %-24s\n' %('Nelec', 'Energy'))
    for k in range(len(occs)):
        PsiMod.print_out('\t%11.3E %24.16E\n' % (occs[k], energies[k]))
        E[occs[k]] = energies[k]

    return E    
