import PsiMod
import re
import os
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
    
    ri_ints_io = PsiMod.get_option('RI_INTS_IO')
    PsiMod.set_global_option('RI_INTS_IO','SAVE')
    activate(molecule) 
    molecule.update_geometry()

    PsiMod.print_out("\n")
    banner("CP Computation: Complex.\nFull Basis Set.")
    PsiMod.print_out("\n")
    e_dimer = energy(name, **kwargs)
    PsiMod.set_global_option('RI_INTS_IO','LOAD')

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
        

    PsiMod.set_global_option('RI_INTS_IO','NONE')
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

    PsiMod.set_global_option('RI_INTS_IO',ri_ints_io)
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
            * smallest: indicates the smallest number of zetas 
                        (or other extrapolation parameter) to include in the calculation.  The default is 2.
            * largest: indicates the largest number of zetas 
                        (or other extrapolation parameter) to include in the extrapolation.  The default is 4.
            * molecule: optionally specify a specific molecule to check.  If no molecule is specified,
                        the last molecule from the input file will be used.
            * wfn: The type of wavefunction to use for the extrapolation.  If this keyword is not
                        specified, the wfn from the globals block is used.
            * largest_correlated: Specify a different largest number of zetas for the correlated portion of
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
        #  Blow away the user's file32 if they like to keep such things.  They shouldn't need it for this type of thing anyway
        #  TODO devise a more elegant solution (if necessary?) for the file32 keeping problem
        path32 = PsiMod.IOManager.shared_object().get_file_path(32);
        filename = path32 + "psi." + str(os.getpid()) + ".32"
        if(os.path.isfile(filename)):
            os.remove(filename)
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
        ahf = a
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
        return (a + ahf)
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


#########################
##  Start of Database  ##
#########################

def database_execute(e_name, db_name, **kwargs):

    import input
    #hartree2kcalmol = 627.508924  # consistent with constants in physconst.h
    hartree2kcalmol = 627.509469  # consistent with perl SETS scripts 

    # Define path and load module for requested database
    sys.path.append('%sdatabases' % (PsiMod.Process.environment["PSIDATADIR"]))
    sys.path.append('./../../../lib/databases')  # for the test suite (better soln needed)
    sys.path.append('./../../../source/lib/databases')  # for the test suite with other directory structures
    try: 
        database = __import__(db_name)
    except ImportError: 
        PsiMod.print_out('\nPython module for database %s failed to load\n\n' % (db_name))
        raise Exception("Python module loading problem for database " + str(db_name))
    else:
        dbse = database.dbse
        HRXN = database.HRXN
        ACTV = database.ACTV
        RXNM = database.RXNM
        BIND = database.BIND
        TAGL = database.TAGL
        GEOS = database.GEOS

    # Must collect (here) and set (below) basis sets after every new molecule activation
    user_basis = PsiMod.get_option('BASIS')
    user_ri_basis_scf = PsiMod.get_option('RI_BASIS_SCF')
    user_ri_basis_mp2 = PsiMod.get_option('RI_BASIS_MP2')
    user_ri_basis_cc = PsiMod.get_option('RI_BASIS_CC')
    user_ri_basis_sapt = PsiMod.get_option('RI_BASIS_SAPT')
    user_reference = PsiMod.get_option('REFERENCE')

    # Temporary alarms until more things are working
    if re.match(r'^dfmp2$', e_name, re.IGNORECASE):
        raise Exception('Databases not yet compatible with DF-MP2 calculations (who knows why).')
    if re.match(r'^mp2c$', e_name, re.IGNORECASE):
        raise Exception('Databases not yet compatible with SAPT or MP2C calculations (supermolecular vs sapt structure).')
    if re.match(r'^ccsd', e_name, re.IGNORECASE):
        raise Exception('Databases not yet compatible with CC calculations (checkpoint file, perhaps).')

    # Configuration based upon e_name & db_name options
    #   Force non-supramolecular if needed
    sapt_override = 0
    if re.match(r'^sapt', e_name, re.IGNORECASE):
        try:
            database.ACTV_SA
        except AttributeError:
            raise Exception('Database %s not suitable for non-supramolecular calculation.' % (db_name))
        else:
            sapt_override = 1
            ACTV = database.ACTV_SA
    #   Force open-shell if needed
    openshell_override = 0
    if user_reference == 'RHF':
        try:
            database.isOS
        except AttributeError:
            pass
        else:
            openshell_override = 1
            PsiMod.print_out('\nSome reagents in database %s require an open-shell reference; will be reset to UHF as needed.\n' % (db_name))

    # Configuration based upon database keyword options
    #   Option mode of operation- whether db run in one job or files farmed out
    db_mode = 'continuous'
    if(kwargs.has_key('mode')):
        db_mode = kwargs['mode'];

    if re.match(r'^sow$', db_mode, re.IGNORECASE):
        raise Exception('Database execution mode \'sow\' not yet implemented.')
    elif re.match(r'^reap$', db_mode, re.IGNORECASE):
        raise Exception('Database execution mode \'reap\' not yet implemented.')
    elif re.match(r'^continuous$', db_mode, re.IGNORECASE):
        pass
    else:
        raise Exception('Database execution mode \'%s\' not valid.' % (db_mode))
      
    #   Option counterpoise- whether for interaction energy databases run in bsse-corrected or not
    db_cp = 'no'
    if(kwargs.has_key('cp')):
        db_cp = kwargs['cp'];

    if input.yes.match(str(db_cp)):
        try:
            database.ACTV_CP
        except AttributeError:
            raise Exception('Counterpoise correction mode \'yes\' invalid for database %s.' % (db_name))
        else:
            ACTV = database.ACTV_CP
    elif input.no.match(str(db_cp)):
        pass
    else:
        raise Exception('Counterpoise correction mode \'%s\' not valid.' % (db_cp))

    #   Option zero-point-correction- whether for thermochem databases jobs are corrected by zpe
    db_zpe = 'no'
    if(kwargs.has_key('zpe')):
        db_zpe = kwargs['zpe'];

    if input.yes.match(str(db_zpe)):
        raise Exception('Zero-point-correction mode \'yes\' not yet implemented.')
    elif input.no.match(str(db_zpe)):
        pass
    else:
        raise Exception('Zero-point-correction \'mode\' %s not valid.' % (db_zpe))

    #   Option subset- whether all of the database or just a portion is run
    db_subset = HRXN
    if(kwargs.has_key('subset')):
        db_subset = kwargs['subset'];

    if isinstance(db_subset, basestring):
        if re.match(r'^small$', db_subset, re.IGNORECASE):
            try:
                database.HRXN_SM
            except AttributeError:
                raise Exception('Special subset \'small\' not available for database %s.' % (db_name))
            else:
                HRXN = database.HRXN_SM
        elif re.match(r'^large$', db_subset, re.IGNORECASE):
            try:
                database.HRXN_LG
            except AttributeError:
                raise Exception('Special subset \'large\' not available for database %s.' % (db_name))
            else:
                HRXN = database.HRXN_LG
        elif re.match(r'^equilibrium$', db_subset, re.IGNORECASE):
            try:
                database.HRXN_EQ
            except AttributeError:
                raise Exception('Special subset \'equilibrium\' not available for database %s.' % (db_name))
            else:
                HRXN = database.HRXN_EQ
        else:
            try:
                getattr(database, db_subset)
            except AttributeError:
                raise Exception('Special subset \'%s\' not available for database %s.' % (db_subset, db_name))
            else:
                HRXN = getattr(database, db_subset)
    else:
        temp = []
        for rxn in db_subset:
            if rxn in HRXN:
                temp.append(rxn)
            else:
                raise Exception('Subset element \'%s\' not a member of database %s.' % (str(rxn), db_name))
        HRXN = temp

    temp = []
    for rxn in HRXN:
        temp.append(ACTV['%s-%s' % (dbse, rxn)])
    HSYS = drop_duplicates(sum(temp, []))

    # Sow all the necessary reagent computations
    PsiMod.print_out("\n")
    banner(("Database %s Computation" % (db_name)))
    PsiMod.print_out("\n")

    ERGT = {}
    ERXN = {}

    #   Loop through chemical systems
    for rgt in HSYS:

        PsiMod.print_out("\n")
        banner(("Database %s Computation:\nReagent %s" % (db_name, TAGL[rgt])))
        PsiMod.print_out("\n")

        exec GEOS[rgt]

        PsiMod.set_global_option('BASIS', user_basis)
        if not((user_ri_basis_scf == "") or (user_ri_basis_scf == 'NONE')):
            PsiMod.set_global_option('RI_BASIS_SCF', user_ri_basis_scf)
        if not((user_ri_basis_mp2 == "") or (user_ri_basis_mp2 == 'NONE')):
            PsiMod.set_global_option('RI_BASIS_MP2', user_ri_basis_mp2)
        if not((user_ri_basis_cc == "") or (user_ri_basis_cc == 'NONE')):
            PsiMod.set_global_option('RI_BASIS_SCF', user_ri_basis_scf)
        if not((user_ri_basis_sapt == "") or (user_ri_basis_sapt == 'NONE')):
            PsiMod.set_global_option('RI_BASIS_SAPT', user_ri_basis_sapt)

        molecule = PsiMod.get_active_molecule()
        molecule.update_geometry()

        if sapt_override:
            molecule.reset_point_group('c1')
            molecule.fix_orientation(1) 
            molecule.update_geometry()

        if (openshell_override) and (molecule.multiplicity() != 1):
            PsiMod.set_global_option('REFERENCE', 'UHF')

        if not molecule:
            raise ValueNotSet("No molecule found.")

        #print 'MOLECULE LIVES %23s %8s %4d %4d %4s' % (rgt, PsiMod.get_option('REFERENCE'), molecule.molecular_charge(), molecule.multiplicity(), molecule.schoenflies_symbol())
        elecenergy = energy(e_name, **kwargs)
        #print elecenergy
        ERGT[rgt] = elecenergy
        PsiMod.set_global_option("REFERENCE", user_reference)
        PsiMod.clean()

    # Reap all the necessary reaction computations
    PsiMod.print_out("\n")
    banner(("Database %s Results" % (db_name)))
    PsiMod.print_out("\n")

    maxactv = []
    for rxn in HRXN:
        maxactv.append(len(ACTV[dbse+'-'+str(rxn)]))
    maxrgt = max(maxactv) 
    table_delimit = '-' * (52+20*maxrgt)

    minDerror = 10000.0
    maxDerror = 0.0
    MSDerror  = 0.0
    MADerror  = 0.0
    RMSDerror = 0.0

    #print '\n   %s' % (table_delimit)
    #print '%23s %19s %8s' % ('Reaction', 'Reaction Energy', 'Error'),
    PsiMod.print_out('\n   %s\n' % (table_delimit))
    PsiMod.print_out('%23s %19s %8s' % ('Reaction', 'Reaction Energy', 'Error'),)
    for i in range(maxrgt):
        #print '%20s' % ('Reagent '+str(i+1)),
        PsiMod.print_out('%20s' % ('Reagent '+str(i+1)),)
    #print '\n%23s %10s %17s' % ('', 'Ref', '[kcal/mol]'),
    PsiMod.print_out('\n%23s %10s %17s' % ('', 'Ref', '[kcal/mol]'),)
    for i in range(maxrgt):
        #print '%20s' % ('[H] Wt'),
        PsiMod.print_out('%20s' % ('[H] Wt'),)
    #print '\n   %s' % (table_delimit),
    PsiMod.print_out('\n   %s' % (table_delimit),)

    count_rxn = 0
    for rxn in HRXN:

        db_rxn = dbse + '-' + str(rxn)
        fail_rxn = 0
        for i in range(len(ACTV[db_rxn])):
            if abs(ERGT[ACTV[db_rxn][i]]) < 1.0e-12:
                fail_rxn = 1

        if fail_rxn:
            #print '\n%23s   %8.4f %8s %8s' % (db_rxn, BIND[db_rxn], '****', '****'),
            PsiMod.print_out('\n%23s   %8.4f %8s %8s' % (db_rxn, BIND[db_rxn], '****', '****'),)
            for i in range(len(ACTV[db_rxn])):
                #print ' %16.8f %2.0f' % (ERGT[ACTV[db_rxn][i]], RXNM[db_rxn][ACTV[db_rxn][i]]),
                PsiMod.print_out(' %16.8f %2.0f' % (ERGT[ACTV[db_rxn][i]], RXNM[db_rxn][ACTV[db_rxn][i]]),)

        else:
            ERXN[db_rxn] = 0.0
            for i in range(len(ACTV[db_rxn])):
                ERXN[db_rxn] += ERGT[ACTV[db_rxn][i]] * RXNM[db_rxn][ACTV[db_rxn][i]]
            error = hartree2kcalmol * ERXN[db_rxn] - BIND[db_rxn]
        
            #print '\n%23s   %8.4f %8.4f %8.4f' % (db_rxn, BIND[db_rxn], hartree2kcalmol*ERXN[db_rxn], error),
            PsiMod.print_out('\n%23s   %8.4f %8.4f %8.4f' % (db_rxn, BIND[db_rxn], hartree2kcalmol*ERXN[db_rxn], error),)
            for i in range(len(ACTV[db_rxn])):
                #print ' %16.8f %2.0f' % (ERGT[ACTV[db_rxn][i]], RXNM[db_rxn][ACTV[db_rxn][i]]),
                PsiMod.print_out(' %16.8f %2.0f' % (ERGT[ACTV[db_rxn][i]], RXNM[db_rxn][ACTV[db_rxn][i]]),)

            if abs(error) < abs(minDerror): minDerror = error
            if abs(error) > abs(maxDerror): maxDerror = error
            MSDerror += error
            MADerror += abs(error)
            RMSDerror += error*error
            count_rxn += 1

    #print '\n   %s' % (table_delimit)
    #print '%23s   %17s %8.4f' % ('Minimal Dev', '', minDerror)
    #print '%23s   %17s %8.4f' % ('Maximal Dev', '', maxDerror)
    #print '%23s   %17s %8.4f' % ('Mean Signed Dev', '', MSDerror/float(count_rxn))
    #print '%23s   %17s %8.4f' % ('Mean Absolute Dev', '', MADerror/float(count_rxn))
    #print '%23s   %17s %8.4f' % ('RMS Dev', '', sqrt(RMSDerror/float(count_rxn)))
    #print '   %s\n' % (table_delimit)

    PsiMod.print_out('\n   %s' % (table_delimit))
    PsiMod.print_out('\n%23s   %17s %8.4f' % ('Minimal Dev', '', minDerror))
    PsiMod.print_out('\n%23s   %17s %8.4f' % ('Maximal Dev', '', maxDerror))
    PsiMod.print_out('\n%23s   %17s %8.4f' % ('Mean Signed Dev', '', MSDerror/float(count_rxn)))
    PsiMod.print_out('\n%23s   %17s %8.4f' % ('Mean Absolute Dev', '', MADerror/float(count_rxn)))
    PsiMod.print_out('\n%23s   %17s %8.4f' % ('RMS Dev', '', sqrt(RMSDerror/float(count_rxn))))
    PsiMod.print_out('\n   %s\n' % (table_delimit))

    return MADerror

def drop_duplicates(seq): 
    noDupes = []
    [noDupes.append(i) for i in seq if not noDupes.count(i)]
    return noDupes

## Aliases  ##
database = database_execute
db = database_execute

#######################
##  End of Database  ##
#######################


