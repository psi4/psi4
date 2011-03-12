import PsiMod
import re
import math
import warnings
from driver import *
from molecule import *
from text import *

def cp(name, **kwargs):

    check_bsse = False
    if (kwargs.has_key('check_bsse')):
        check_bsse = kwargs['check_bsse']

    molecule = PsiMod.get_active_molecule() 
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')    
 
    if not molecule:
        raise ValueNotSet("no molecule found")
    
    activate(molecule) 
    molecule.update_geometry()

    PsiMod.print_out("\n")
    banner("CP Computation: Complex.\nFull Basis Set.")
    PsiMod.print_out("\n")
    e_dimer = energy(name, **kwargs)


    #All monomers with ghosts
    monomers = extract_clusters(molecule, True, 1)
    e_monomer_full = []
    
    cluster_n = 0
    for cluster in monomers:
        activate(cluster)
        PsiMod.print_out("\n")
        banner(("CP Computation: Monomer %d.\n Full Basis Set." % (cluster_n + 1)))
        PsiMod.print_out("\n")
        e_monomer_full.append(energy(name,**kwargs))
        cluster_n = cluster_n + 1
        
    if (check_bsse): 
        #All monomers without ghosts
        monomers = extract_clusters(molecule, False, 1)
        e_monomer_bsse = []
   
        cluster_n = 0
        for cluster in monomers:
            activate(cluster)
            PsiMod.print_out("\n")
            #cluster.print_to_output()
            banner(("CP Computation: Monomer %d.\n Monomer Set." % (cluster_n + 1)))
            PsiMod.print_out("\n")
            e_monomer_bsse.append(energy(name,**kwargs))
            cluster_n = cluster_n + 1

    activate(molecule) 
        
    if (check_bsse != True):
        cp_table = Table(rows = ["System:"], cols = ["Energy (full):"])
        cp_table["Complex"] = [e_dimer]
        for cluster_n in range(0,len(monomers)):
            key = "Monomer %d" % (cluster_n+1)
            cp_table[key] = [e_monomer_full[cluster_n]]
        
        e_full = e_dimer
        for cluster_n in range(0,len(monomers)):
            e_full = e_full - e_monomer_full[cluster_n]
        cp_table["Interaction"] = [e_full]
    
    else:
        cp_table = Table(rows = ["System:"], cols = ["Energy (full):","Energy (monomer):","BSSE:"])
        cp_table["Complex"] = [e_dimer, 0.0, 0.0]
        for cluster_n in range(0,len(monomers)):
            key = "Monomer %d" % (cluster_n+1)
            cp_table[key] = [e_monomer_full[cluster_n], e_monomer_bsse[cluster_n], \
                e_monomer_full[cluster_n] - e_monomer_bsse[cluster_n]]
        
        e_full = e_dimer
        e_bsse = e_dimer
        for cluster_n in range(0,len(monomers)):
            e_full = e_full - e_monomer_full[cluster_n]
            e_bsse = e_bsse - e_monomer_bsse[cluster_n]
        cp_table["Totals:"] = [e_full, e_bsse, e_full-e_bsse]
    
    PsiMod.print_out("\n")
    banner("CP Computation: Results.")
    PsiMod.print_out("\n")
    
    banner("Hartree",2)
    PsiMod.print_out("\n")
    
    PsiMod.print_out(str(cp_table))

    PsiMod.print_out("\n")
    banner("kcal*mol^-1",2)
    PsiMod.print_out("\n")

    cp_table.scale()   
 
    PsiMod.print_out(str(cp_table))
    return e_full

#############
# End of cp #
#############

## Aliases ##
# Easily remembered longhand
counterpoise_correct = cp
counterpoise_correction = cp





# Basis Set Extrapolation
def basis_set_extrapolate(basis_name, **kwargs):
    """ Performs a basis set extrapolation substituting the tag [X] in 
        the basis name for letters given in an array by the optional 'labels' keyword.  The indices in
        labels must correspond to the number of zetas (or other extrapolation parameter).  The default is
             ['','', 'D', 'T', 'Q', '5', '6', '7', '8', '9'].  Other optional keywords include:
            * smallest:\t\tindicates the smallest number of zetas 
                        (or other extrapolation parameter) to include in the calculation.  The default is 2.
            * largest: \t\tindicates the largest number of zetas 
                        (or other extrapolation parameter) to include in the extrapolation.  The default is 4.
            * molecule:\t\toptionally specify a specific molecule to check.  If no molecule is specified,
                        the last molecule from the input file will be used.
            * wfn:     \t\tThe type of wavefunction to use for the extrapolation.  If this keyword is not
                        specified, the wfn from the globals block is used.
            * largest_correlated:\tSpecify a different largest number of zetas for the correlated portion of
                        the basis set extrapolation.
    """
    # TODO Implement more extrapolation schemes and put them in as an "extrapolation type" option   
    

    # Back up global options that may be changed
    backup_basis = PsiMod.get_global_option("BASIS")
    backup_wfn = PsiMod.get_global_option("WFN")


    labels = ['','','D','T','Q','5','6','7','8','9']
    if (kwargs.has_key('labels')):
        labels = kwargs['labels']    

    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')    
 
    if not molecule:
        raise ValueNotSet("no molecule found")
    
    activate(molecule) 
    molecule.update_geometry()

    PsiMod.print_out("\n")
    banner("Basis Set Extrapolation")
    PsiMod.print_out("\n")
    
    smallest = 2
    if(kwargs.has_key('smallest')):
        smallest = kwargs['smallest'];

    largest = 4
    if(kwargs.has_key('largest')):
        largest = kwargs['largest'];
    
    largest_correlated = largest
    if(kwargs.has_key('largest_correlated')):
        largest_correlated = kwargs['largest_correlated'];

    is_hf = False
    wavefunction = PsiMod.get_global_option("WFN");
    if(kwargs.has_key('wfn')):
        wavefunction = kwargs['wfn'];
    if not wavefunction:
        is_hf = True  # This is just a HF basis set extrapolation
    if(wavefunction.lower() == "mp2"):
        PsiMod.print_out("Warning: Plain mp2 energy not available.  Using dfmp2 instead.")
        wavefunction = "dfmp2"

    reference = PsiMod.get_global_option("REFERENCE");
    
    if(is_hf and kwargs.has_key('largest_correlated')):
        PsiMod.print_out("Warning: largest_correlated was specified in basis set extrapolation, but no WFN was specified,\n\
                and so a HF extrapolation was assumed.  \
                Are you sure this is what you meant to do?\n")

    if(largest < largest_correlated):
        PsiMod.print_out("Warning: largest_correlated (=" + str(largest_correlated) + \
                ") was set to be larger than largest (=" + str(largest) +") basis set parameter.  This doesn't really make sense.\n\
                Are you sure this is what you meant to do?  largest will be set to largest_correlated\n")
        largest = largest_correlated
    

    # Need at least 3 calculations to do HF extrapolation
    if(largest - smallest < 2):
        raise "Not enough basis sets to do SCF extrapolation"
    # Need at least 2 calculations to do correlated extrapolation
    if(largest_correlated - smallest < 1):
        raise "Not enough basis sets to do correlated extrapolation"
    

    def basis_for(zetas):
        return re.sub('\[X\]', labels[zetas], basis_name)

    def get_energy(n_zetas, hf_only=is_hf):
        basis = basis_for(n_zetas)
        PsiMod.print_out("\n")
        banner("Finding " + wavefunction + " energy with basis set " + basis)
        PsiMod.print_out("\n")
        PsiMod.set_local_option("SCF", "REFERENCE", reference)
        if (hf_only):
            PsiMod.set_local_option("SCF", "BASIS", basis)
            PsiMod.clean()
            PsiMod.scf()
            ret_val = PsiMod.get_variable("CURRENT ENERGY") 
            PsiMod.clean()
            return ret_val
        else:
            PsiMod.clean()
            PsiMod.set_global_option("BASIS", basis)
            PsiMod.set_global_option("WFN", wavefunction)  # Can't set local option since we don't know which module will be run necessarily
            energy(wavefunction.lower())
            ret_val = PsiMod.get_variable("SCF ENERGY")
            cor_val = PsiMod.get_variable("CURRENT ENERGY") 
            PsiMod.clean()
            return ret_val, cor_val
        

    def print_energies(e_array, extrapolated_vals, e_corr_array = []):
        #TODO print extrapolation coefficients
        PsiMod.print_out("\n")
        banner("Basis Set Extrapolation Results")
        PsiMod.print_out("\n")
        if(is_hf):
            width = 22
            PsiMod.print_out("-"*(width * 2 + 5) + "\n")
            PsiMod.print_out((("|%"+ str(width) + "s |") % "Basis Set") + (("%" + str(width) + "s |\n") % (reference + " Reference Energy"))) 
            PsiMod.print_out("-"*(width * 2 + 5) + "\n")
            for i in range(largest+1 - smallest):
                PsiMod.print_out((("|%"+ str(width) + "s |") % basis_for(smallest + i)) + (("%" + str(width)+"."+str(width-6)+"f |\n") % e_array[i]))
            PsiMod.print_out("-"*(width * 2 + 5) + "\n")
            PsiMod.print_out((("|%"+ str(width) + "s |") % "Extrapolated Value") + (("%" + str(width)+"."+str(width-6)+"f |\n") % extrapolated_vals[0]))
            PsiMod.print_out("-"*(width * 2 + 5) + "\n")
        else:
            width = 25
            str_form = "%" + str(width) + "s"
            float_form = "%" + str(width) + "." + str(width-6) + "f"
            PsiMod.print_out("-"*(width * 4 + 11) + "\n")
            PsiMod.print_out((("|" + str_form + " |") % "Basis Set") + ((str_form + " | ") % (reference + " Reference Energy"))) 
            PsiMod.print_out(((str_form + " | ") % (wavefunction + " Correlation Energy")) + ((str_form + " |\n") % ("Total Energy"))) 
            PsiMod.print_out("-"*(width * 4 + 11) + "\n")
            for i in range(largest+1 - smallest):
                PsiMod.print_out((("|"+ str_form + " |") % basis_for(smallest + i)) + ((float_form + " | ") % e_array[i]))
                if(smallest + i <= largest_correlated):
                    PsiMod.print_out(((float_form + " | ") % e_corr_array[i]) + ((float_form + " |\n") % (e_array[i] + e_corr_array[i])))
                else:
                    PsiMod.print_out(((str_form + " | ") % "[Not Calculated]") + ((str_form + " |\n") % "[Not Calculated]"))
            PsiMod.print_out("-"*(width * 4 + 11) + "\n")
            PsiMod.print_out((("|"+ str_form + " |") % "Extrapolated Values") + ((float_form + " | ") % extrapolated_vals[0]))
            PsiMod.print_out(((float_form + " | ") % extrapolated_vals[1]) + ((float_form + " |\n") % (extrapolated_vals[0] + extrapolated_vals[1])))
            PsiMod.print_out("-"*(width * 4 + 11) + "\n")


    e = []
    ehf = []
    
    a = 0.0
    ahf = 0.0

    if(is_hf):  # Then all we need to do is the HF extrapolation
        for i in range(smallest, largest + 1):
            ehf.append(get_energy(i))
        x = range(smallest, largest + 1)
        # For now, just a 3-point extrapolation.
        delta_e = (ehf[-2] - ehf[-3])/(ehf[-1] - ehf[-2])
        if(delta_e < 0.0):
            PsiMod.print_out("\nSCF extrapolation problem: (ehf[n-1] - ehf[n-2])/(ehf[n] - ehf[n-1]) is negative:\n\
                    \tcannot take the log of a negative number.  Printing what I can...\n\n")
            print_energies(ehf, [float("nan"), float("nan")], e)
            raise Exception("SCF extrapolation problem: (ehf[n-1] - ehf[n-2])/(ehf[n] - ehf[n-1]) is negative; \
                    \ncannot take the log of a negative number: " \
                    + str(delta_e))  
        c = math.log(delta_e)
        b = (ehf[-1] - ehf[-2])/(math.exp(-c * x[-1]) - math.exp(-c * x[-2]))
        a = ehf[-1] - (b * math.exp(-c * x[-1]))
        print_energies(ehf, [a])
    else:
        for i in range(smallest, largest_correlated + 1):
            tmp_hf, tmp = get_energy(i)
            e.append(tmp - tmp_hf)
            ehf.append(tmp_hf)
        
        for i in range(largest_correlated + 1, largest + 1):
            tmp_hf = get_energy(i, True)
            ehf.append(tmp_hf)
        
        # Three-point HF extrapolation
        x = range(smallest, largest + 1)
        delta_e = (ehf[-2] - ehf[-3])/(ehf[-1] - ehf[-2])
        if(delta_e < 0.0):
            PsiMod.print_out("\nSCF extrapolation problem: (ehf[n-1] - ehf[n-2])/(ehf[n] - ehf[n-1]) is negative:\n\
                    \tcannot take the log of a negative number.  Printing what I can...\n\n")
            print_energies(ehf, [float("nan"), float("nan")], e)
            raise Exception("SCF extrapolation problem: (ehf[n-1] - ehf[n-2])/(ehf[n] - ehf[n-1]) is negative; \
                    \ncannot take the log of a negative number: " \
                    + str(delta_e))

        
        c = math.log(delta_e)
        bhf = (ehf[-1] - ehf[-2])/(math.exp(-c * x[-1]) - math.exp(-c * x[-2]))
        ahf = ehf[-1] - (bhf * math.exp(-c * x[-1]))

        
        # Two-point Correlated extrapolation
        x = range(smallest, largest_correlated + 1)
        b = (e[-1] - e[-2])/(x[-1]**(-3) - x[-2]**(-3))
        a = e[-1] - (b * x[-1]**(-3))


        print_energies(ehf, [ahf, a], e)
        

        
    # Restore global options that may be changed
    PsiMod.set_global_option("BASIS", backup_basis)
    PsiMod.set_global_option("WFN", backup_wfn)
    
    # Return the correlated extrapolation value if we have it
    if not is_hf:
        return a + ahf
    else:
        return ahf

################################
# End of basis_set_extrapolate #
################################

## Aliases ##
# Helpful shorthand
bse = basis_set_extrapolate
# Easily mistyped long-hand
basis_set_extrapolation = basis_set_extrapolate
