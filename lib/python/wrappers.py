import PsiMod
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

