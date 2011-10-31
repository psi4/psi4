import PsiMod
import re
import os
import math
import warnings
from driver import *
from molecule import *
from text import *

###################
##  Start of cp  ##
###################

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

##  Aliases  ##
counterpoise_correct = cp
counterpoise_correction = cp

#################
##  End of cp  ##
#################



######################################
##  Start of basis_set_extrapolate  ##
######################################

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

##  Aliases  ##
bse = basis_set_extrapolate
basis_set_extrapolation = basis_set_extrapolate

####################################
##  End of basis_set_extrapolate  ##
####################################



#########################
##  Start of Database  ##
#########################

def database(name, db_name, **kwargs):
    """Wrapper to access the molecule objects and reference energies of popular chemical databases.

    Required Arguments:
    -------------------
    * name (or unlabeled first argument) indicates the computational method to be applied to the database.
    * db_name (or unlabeled second argument) is a string of the requested database name. This name matches the
        python file in psi4/lib/databases , wherein literature citations for the database are also documented.

    Optional Arguments:  --> 'default_option' <-- 
    ------------------
    * mode = --> 'continuous' <-- | 'sow' | 'reap'
        Indicates whether the calculation required to complete the database are to be run in one
        file ('continuous') or are to be farmed out in an embarrassingly parallel fashion ('sow'/'reap').
        For the latter, run an initial job with 'sow' and follow instructions in its output file.
    * cp = 'on' | --> 'off' <--
        Indicates whether counterpoise correction is employed in computing interaction energies.
        Use this option and NOT the cp() wrapper for BSSE correction in the database() wrapper.
        Option valid only for databases consisting of bimolecular complexes.
    * symm = --> 'on' <-- | 'off'
        Indicates whether the native symmetry of the database molecules is employed ('on') or whether
        it is forced to c1 symmetry ('off'). Some computational methods (e.g., SAPT) require no symmetry,
        and this will be set by the database() wrapper.
    * zpe = 'on' | --> 'off' <--
        Indicates whether zero-point-energy corrections are appended to single-point energy values. Option
        valid only for certain thermochemical databases.
        Disabled until Hessians ready.
    * benchmark = --> 'default' <-- | etc.
        Indicates whether a non-default set of reference energies, if available, are employed for the 
        calculation of error statistics.
    * tabulate = --> [] <-- | etc.
        Indicates whether to form tables of variables other than the primary requested energy
    * subset
        Indicates a subset of the full database to run. This is a very flexible option and can used in
        three distinct ways, outlined below. Note that two take a string and the last takes a list.
        * subset = 'small' | 'large' | 'equilibrium'
            Calls predefined subsets of the requested database, either 'small', a few of the smallest
            database members, 'large', the largest of the database members, or 'equilibrium', the
            equilibrium geometries for a database composed of dissociation curves.
        * subset = 'BzBz_S' | 'FaOOFaON' | 'ArNe' | etc.
            For databases composed of dissociation curves, individual dissociation curves can be called
            by name. Consult the database python files for available molecular systems.
        * subset = [1,2,5] | ['1','2','5'] | ['BzMe-3.5', 'MeMe-5.0'] | etc.
            Specify a list of database members to run. Consult the database python files for available 
            molecular systems.

    Examples:
    ---------
    A database job requires, at a minimum, the basis set to be set in its input file. The following
    are valid example calls for the database() wrapper.
        database('scf','S22')
    """

    import input
    import pickle
    from collections import defaultdict
    #hartree2kcalmol = 627.508924  # consistent with constants in physconst.h
    hartree2kcalmol = 627.509469  # consistent with perl SETS scripts 

    # Wrap any positional arguments into kwargs (for intercalls among wrappers)
    if not('name' in kwargs) and name:
        kwargs['name'] = name
    if not('db_name' in kwargs) and db_name:
        kwargs['db_name'] = db_name
    if not('func' in kwargs):
        kwargs['func'] = energy

    # Define path and load module for requested database
    sys.path.append('%sdatabases' % (PsiMod.Process.environment["PSIDATADIR"]))
    sys.path.append('%s/lib/databases' % PsiMod.psi_top_srcdir())
    try: 
        database = __import__(db_name)
    except ImportError: 
        PsiMod.print_out('\nPython module for database %s failed to load\n\n' % (db_name))
        PsiMod.print_out('\nSearch path that was tried:\n')
        PsiMod.print_out(", ".join(map(str, sys.path)))
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
    user_puream = PsiMod.get_option('PUREAM')
    user_ri_basis_scf = PsiMod.get_option('RI_BASIS_SCF')
    user_ri_basis_mp2 = PsiMod.get_option('RI_BASIS_MP2')
    user_ri_basis_cc = PsiMod.get_option('RI_BASIS_CC')
    user_ri_basis_sapt = PsiMod.get_option('RI_BASIS_SAPT')
    user_ri_basis_elst = PsiMod.get_option('RI_BASIS_ELST')
    user_reference = PsiMod.get_option('REFERENCE')

    user_memory = PsiMod.get_memory()

    # Configuration based upon e_name & db_name options
    #   Force non-supramolecular if needed
    symmetry_override = 0
    if re.match(r'^sapt', name, re.IGNORECASE):
        try:
            database.ACTV_SA
        except AttributeError:
            raise Exception('Database %s not suitable for non-supramolecular calculation.' % (db_name))
        else:
            symmetry_override = 1
            ACTV = database.ACTV_SA
    #   Force open-shell if needed
    openshell_override = 0
    if (user_reference == 'RHF') or (user_reference == 'RKS'):
        try:
            database.isOS
        except AttributeError:
            pass
        else:
            if input.yes.match(str(database.isOS)):
                openshell_override = 1
                PsiMod.print_out('\nSome reagents in database %s require an open-shell reference; will be reset to UHF/UKS as needed.\n' % (db_name))

    # Configuration based upon database keyword options
    #   Option symmetry- whether symmetry treated normally or turned off (currently req'd for dfmp2 & dft)
    db_symm = 'yes'
    if(kwargs.has_key('symm')):
        db_symm = kwargs['symm']

    if input.no.match(str(db_symm)):
        symmetry_override = 1
    elif input.yes.match(str(db_symm)):
        pass
    else:
        raise Exception('Symmetry mode \'%s\' not valid.' % (db_symm))

    #   Option mode of operation- whether db run in one job or files farmed out
    db_mode = 'continuous'
    if(kwargs.has_key('mode')):
        db_mode = kwargs['mode']

    if re.match(r'^continuous$', db_mode.lower()):
        pass
    elif re.match(r'^sow$', db_mode.lower()):
        pass
    elif re.match(r'^reap$', db_mode.lower()):
        if(kwargs.has_key('linkage')):
            db_linkage = kwargs['linkage']
        else:
            raise Exception('Database execution mode \'reap\' requires a linkage option.')
    else:
        raise Exception('Database execution mode \'%s\' not valid.' % (db_mode))

    #   Option counterpoise- whether for interaction energy databases run in bsse-corrected or not
    db_cp = 'no'
    if(kwargs.has_key('cp')):
        db_cp = kwargs['cp']

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
        db_zpe = kwargs['zpe']

    if input.yes.match(str(db_zpe)):
        raise Exception('Zero-point-correction mode \'yes\' not yet implemented.')
    elif input.no.match(str(db_zpe)):
        pass
    else:
        raise Exception('Zero-point-correction \'mode\' %s not valid.' % (db_zpe))

    #   Option benchmark- whether error statistics computed wrt alternate reference energies
    db_benchmark = 'default'
    if(kwargs.has_key('benchmark')):
        db_benchmark = kwargs['benchmark']

        if re.match(r'^default$', db_benchmark, re.IGNORECASE):
            pass
        else:
            try:
                getattr(database, 'BIND_' + db_benchmark)
            except AttributeError:
                raise Exception('Special benchmark \'%s\' not available for database %s.' % (db_benchmark, db_name))
            else:
                BIND = getattr(database, 'BIND_' + db_benchmark)

    #   Option tabulate- whether tables of variables other than primary energy method are formed
    db_tabulate = []
    if(kwargs.has_key('tabulate')):
        db_tabulate = kwargs['tabulate']

    #   Option subset- whether all of the database or just a portion is run
    db_subset = HRXN
    if(kwargs.has_key('subset')):
        db_subset = kwargs['subset']

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
    PsiMod.print_out("\n\n")
    banner(("Database %s Computation" % (db_name)))
    PsiMod.print_out("\n")

    #   write index of calcs to output file
    if re.match('continuous', db_mode.lower()):
        instructions  = """\n    The database single-job procedure has been selected through mode='continuous'.\n"""
        instructions +=   """    Calculations for the reagents will proceed in the order below and will be followed\n"""
        instructions +=   """    by summary results for the database.\n\n"""
        for rgt in HSYS:
            instructions += """                    %-s\n""" % (rgt)
        instructions += """\n    Alternatively, a farming-out of the database calculations may be accessed through\n"""
        instructions +=   """    the database wrapper option mode='sow'/'reap'.\n\n"""
        PsiMod.print_out(instructions)

    #   write sow/reap instructions and index of calcs to output file and reap input file
    if re.match('sow', db_mode.lower()):
        instructions  = """\n    The database sow/reap procedure has been selected through mode='sow'. In addition\n"""
        instructions +=   """    to this output file (which contains no quantum chemical calculations), this job\n"""
        instructions +=   """    has produced a number of input files (%s-*.in) for individual database members\n""" % (dbse)
        instructions +=   """    and a single input file (%s-master.in) with a database(mode='reap') command.\n""" % (dbse)
        instructions +=   """    The former may look very peculiar since processed and pickled python rather than\n"""
        instructions +=   """    raw input is written. Follow the instructions below to continue.\n\n"""
        instructions +=   """    (1)  Run all of the %s-*.in input files on any variety of computer architecture.\n""" % (dbse)
        instructions +=   """       The output file names must be as given below.\n\n"""
        for rgt in HSYS:
            instructions += """             psi4 -i %-27s -o %-27s\n""" % (rgt + '.in', rgt + '.out')
        instructions += """\n    (2)  Gather all the resulting output files in a directory. Place input file\n"""
        instructions +=   """         %s-master.in into that directory and run it. The job will be trivial in\n""" % (dbse)
        instructions +=   """         length and give summary results for the database in its output file.\n\n"""
        instructions +=   """             psi4 -i %-27s -o %-27s\n\n""" % (dbse + '-master.in', dbse + '-master.out')
        instructions +=   """    Alternatively, a single-job execution of the database may be accessed through\n"""
        instructions +=   """    the database wrapper option mode='continuous'.\n\n"""
        PsiMod.print_out(instructions)

        fmaster = open('%s-master.in' % (dbse), 'w')
        fmaster.write('# This is a psi4 input file auto-generated from the database() wrapper.\n\n')
        fmaster.write("database('%s', '%s', mode='reap', cp='%s', zpe='%s', benchmark='%s', linkage=%d, subset=%s, tabulate=%s)\n\n" %
            (name, db_name, db_cp, db_zpe, db_benchmark, os.getpid(), HRXN, db_tabulate))
        fmaster.close()

    #   Loop through chemical systems
    ERGT = {}
    ERXN = {}
    VRGT = {}
    VRXN = {}
    for rgt in HSYS:
        VRGT[rgt] = {}

        # extra definition of molecule so that logic in building commands string has something to act on
        exec GEOS[rgt]
        molecule = PsiMod.get_active_molecule()
        if not molecule:
            raise ValueNotSet("No molecule found.")

        # build string of title banner
        banners = ''
        banners += """PsiMod.print_out('\\n')\n"""
        banners += """banner(' Database %s Computation: Reagent %s \\n   %s')\n""" % (db_name, rgt, TAGL[rgt])
        banners += """PsiMod.print_out('\\n')\n\n"""

        # build string of lines that defines contribution of rgt to each rxn
        actives = ''
        actives += """PsiMod.print_out('   Database Contributions Map:\\n   %s\\n')\n""" % ('-' * 75)
        for rxn in HRXN:
            db_rxn = dbse + '-' + str(rxn)
            if rgt in ACTV[db_rxn]:
                actives += """PsiMod.print_out('   reagent %s contributes by %.4f to reaction %s\\n')\n""" \
                   % (rgt, RXNM[db_rxn][rgt], db_rxn)
        actives += """PsiMod.print_out('\\n')\n\n"""

        # build string of commands for options from the input file  TODO: handle local options too
        commands  = ''
        commands += """\nPsiMod.set_memory(%s)\n\n""" % (user_memory)
        for chgdopt in PsiMod.get_global_option_list():
            if PsiMod.has_option_changed(chgdopt):
                if chgdopt != 'PUREAM':
                    commands += """PsiMod.set_global_option('%s', '%s')\n""" % (chgdopt, PsiMod.get_global_option(chgdopt))

        # build string of molecule and commands that are dependent on the database
        commands += '\n'
        commands += """PsiMod.set_global_option('BASIS', '%s')\n""" % (user_basis)
        if input.yes.match(str(user_puream)):
            commands += """PsiMod.set_global_option('PUREAM', 1)\n""" 
        elif input.no.match(str(user_puream)):
            commands += """PsiMod.set_global_option('PUREAM', 0)\n""" 
        if not((user_ri_basis_scf == "") or (user_ri_basis_scf == 'NONE')):
            commands += """PsiMod.set_global_option('RI_BASIS_SCF', '%s')\n""" % (user_ri_basis_scf)
        if not((user_ri_basis_mp2 == "") or (user_ri_basis_mp2 == 'NONE')):
            commands += """PsiMod.set_global_option('RI_BASIS_MP2', '%s')\n""" % (user_ri_basis_mp2)
        if not((user_ri_basis_cc == "") or (user_ri_basis_cc == 'NONE')):
            commands += """PsiMod.set_global_option('RI_BASIS_CC', '%s')\n""" % (user_ri_basis_cc)
        if not((user_ri_basis_sapt == "") or (user_ri_basis_sapt == 'NONE')):
            commands += """PsiMod.set_global_option('RI_BASIS_SAPT', '%s')\n""" % (user_ri_basis_sapt)
        if not((user_ri_basis_elst == "") or (user_ri_basis_elst == 'NONE')):
            commands += """PsiMod.set_global_option('RI_BASIS_ELST', '%s')\n""" % (user_ri_basis_elst)
        commands += """molecule = PsiMod.get_active_molecule()\n"""
        commands += """molecule.update_geometry()\n"""

        if symmetry_override:
            commands += """molecule.reset_point_group('c1')\n"""
            commands += """molecule.fix_orientation(1)\n"""
            commands += """molecule.update_geometry()\n"""

        if (openshell_override) and (molecule.multiplicity() != 1):
            if user_reference == 'RHF': 
                commands += """PsiMod.set_global_option('REFERENCE', 'UHF')\n"""
            elif user_reference == 'RKS': 
                commands += """PsiMod.set_global_option('REFERENCE', 'UKS')\n"""

        # all modes need to step through the reagents but all for different purposes
        # continuous: defines necessary commands, executes energy(method) call, and collects results into dictionary
        # sow: opens individual reagent input file, writes the necessary commands, and writes energy(method) call
        # reap: opens individual reagent output file, collects results into a dictionary
        if re.match('continuous', db_mode.lower()):
            exec banners
            exec GEOS[rgt]
            exec commands
            #print 'MOLECULE LIVES %23s %8s %4d %4d %4s' % (rgt, PsiMod.get_option('REFERENCE'),
            #    molecule.molecular_charge(), molecule.multiplicity(), molecule.schoenflies_symbol())
            lesserkwargs = kwargs.copy() 
            del lesserkwargs['func'] 
            ERGT[rgt] = assemble_function2call_dropfunc(kwargs['func'], **lesserkwargs)
            #print ERGT[rgt]
            PsiMod.print_variables()
            exec actives
            for envv in db_tabulate:
               VRGT[rgt][envv] = PsiMod.get_variable(envv)
            PsiMod.set_global_option("REFERENCE", user_reference)
            PsiMod.clean()
            
        elif re.match('sow', db_mode.lower()):
            freagent = open('%s.in' % (rgt), 'w')
            freagent.write('# This is a psi4 input file auto-generated from the database() wrapper.\n\n')
            freagent.write(banners)
            freagent.write(GEOS[rgt])
            freagent.write(commands)
            lesserkwargs = kwargs.copy() 
            del lesserkwargs['func'] 
            # non-pickle route
            #freagent.write("""\nkwargs = %s\n""" % (lesserkwargs)) 
            # pickle route (conflics w/Andy's parenthesis matching)
            freagent.write('''\npickle_kw = ("""''')
            pickle.dump(lesserkwargs, freagent)
            freagent.write('''""")\n''')
            freagent.write("""\nkwargs = pickle.loads(pickle_kw)\n""")
            # end routes
            freagent.write("""electronic_energy = %s(**kwargs)\n\n""" % (kwargs['func'].__name__))
            freagent.write("""PsiMod.print_variables()\n""")
            freagent.write("""PsiMod.print_out('\\nDATABASE RESULT: computation %d for reagent %s """
                % (os.getpid(), rgt))
            freagent.write("""yields electronic energy %20.12f\\n' % (electronic_energy))\n\n""")
            for envv in db_tabulate:
                freagent.write("""PsiMod.print_out('DATABASE RESULT: computation %d for reagent %s """
                    % (os.getpid(), rgt))
                freagent.write("""yields variable value    %20.12f for variable %s\\n' % (PsiMod.get_variable(""")
                freagent.write("""'%s'), '%s'))\n""" % (envv.upper(), envv.upper()))
            freagent.close()

        elif re.match('reap', db_mode.lower()):
            ERGT[rgt] = 0.0
            for envv in db_tabulate:
               VRGT[rgt][envv] = 0.0
            exec banners
            exec actives
            try:
                freagent = open('%s.out' % (rgt), 'r')
            except IOError:
                PsiMod.print_out('Warning: Output file \'%s.out\' not found.\n' % (rgt))
                PsiMod.print_out('         Database summary will have 0.0 and **** in its place.\n')
            else:
                while 1:
                    line = freagent.readline()
                    if not line:
                        if ERGT[rgt] == 0.0:
                           PsiMod.print_out('Warning: Output file \'%s.out\' has no DATABASE RESULT line.\n' % (rgt))
                           PsiMod.print_out('         Database summary will have 0.0 and **** in its place.\n')
                        break
                    s = line.split()
                    if (len(s) != 0) and (s[0:3] == ['DATABASE', 'RESULT:', 'computation']):
                        if int(s[3]) != db_linkage:
                           raise Exception('Output file \'%s.out\' has linkage %s incompatible with master.in linkage %s.' 
                               % (rgt, str(s[3]), str(db_linkage)))
                        if s[6] != rgt:
                           raise Exception('Output file \'%s.out\' has nominal affiliation %s incompatible with reagent %s.' 
                               % (rgt, s[6], rgt))
                        if (s[8:10] == ['electronic', 'energy']):
                            ERGT[rgt] = float(s[10])
                            PsiMod.print_out('DATABASE RESULT: electronic energy = %20.12f\n' % (ERGT[rgt]))
                        elif (s[8:10] == ['variable', 'value']):
                            for envv in db_tabulate:
                                if (s[13:] == envv.upper().split()):
                                    VRGT[rgt][envv] = float(s[10])
                                    PsiMod.print_out('DATABASE RESULT: variable value    = %20.12f\n' % (VRGT[rgt][envv]))
                freagent.close()

    #   end sow after writing files 
    if re.match('sow', db_mode.lower()):
        return 0.0

    # Reap all the necessary reaction computations
    PsiMod.print_out("\n")
    banner(("Database %s Results" % (db_name)))
    PsiMod.print_out("\n")

    maxactv = []
    for rxn in HRXN:
        maxactv.append(len(ACTV[dbse+'-'+str(rxn)]))
    maxrgt = max(maxactv) 
    table_delimit = '-' * (52+20*maxrgt)
    tables = ''

    #   find any reactions that are incomplete
    FAIL = defaultdict(int)
    for rxn in HRXN:
        db_rxn = dbse + '-' + str(rxn)
        for i in range(len(ACTV[db_rxn])):
            if abs(ERGT[ACTV[db_rxn][i]]) < 1.0e-12:
                FAIL[rxn] = 1

    #   tabulate requested process::environment variables
    tables += """   For each VARIABLE requested by tabulate, a 'Reaction Value' will be formed from\n"""
    tables += """   'Reagent' values according to weightings 'Wt', as for the REQUESTED ENERGY below.\n"""
    tables += """   Depending on the nature of the variable, this may or may not make any physical sense.\n"""
    for envv in db_tabulate:
        tables += """\n\n   %s""" % (envv.upper())
        tables += tblhead(maxrgt, table_delimit, 2)

        for rxn in HRXN:
            db_rxn = dbse + '-' + str(rxn)
            VRXN[db_rxn] = {}
    
            if FAIL[rxn]:
                tables += """\n%23s   %8.4f %8s %8s""" % (db_rxn, BIND[db_rxn], '****', '****')
                for i in range(len(ACTV[db_rxn])):
                    tables += """ %16.8f %2.0f""" % (VRGT[ACTV[db_rxn][i]][envv], RXNM[db_rxn][ACTV[db_rxn][i]])
    
            else:
                VRXN[db_rxn][envv] = 0.0
                for i in range(len(ACTV[db_rxn])):
                    VRXN[db_rxn][envv] += VRGT[ACTV[db_rxn][i]][envv] * RXNM[db_rxn][ACTV[db_rxn][i]]
            
                tables += """\n%23s        %16.8f     """ % (db_rxn, VRXN[db_rxn][envv])
                for i in range(len(ACTV[db_rxn])):
                    tables += """ %16.8f %2.0f""" % (VRGT[ACTV[db_rxn][i]][envv], RXNM[db_rxn][ACTV[db_rxn][i]])
        tables += """\n   %s\n""" % (table_delimit)

    #   tabulate primary requested energy variable with statistics
    count_rxn = 0
    minDerror = 100000.0
    maxDerror = 0.0
    MSDerror  = 0.0
    MADerror  = 0.0
    RMSDerror = 0.0

    tables += """\n\n   %s""" % ("REQUESTED ENERGY")
    tables += tblhead(maxrgt, table_delimit, 1)
    for rxn in HRXN:
        db_rxn = dbse + '-' + str(rxn)

        if FAIL[rxn]:
            tables += """\n%23s   %8.4f %8s %8s""" % (db_rxn, BIND[db_rxn], '****', '****')
            for i in range(len(ACTV[db_rxn])):
                tables += """ %16.8f %2.0f""" % (ERGT[ACTV[db_rxn][i]], RXNM[db_rxn][ACTV[db_rxn][i]])

        else:
            ERXN[db_rxn] = 0.0
            for i in range(len(ACTV[db_rxn])):
                ERXN[db_rxn] += ERGT[ACTV[db_rxn][i]] * RXNM[db_rxn][ACTV[db_rxn][i]]
            error = hartree2kcalmol * ERXN[db_rxn] - BIND[db_rxn]
        
            tables += """\n%23s   %8.4f %8.4f %8.4f""" % (db_rxn, BIND[db_rxn], hartree2kcalmol*ERXN[db_rxn], error)
            for i in range(len(ACTV[db_rxn])):
                tables += """ %16.8f %2.0f""" % (ERGT[ACTV[db_rxn][i]], RXNM[db_rxn][ACTV[db_rxn][i]])

            if abs(error) < abs(minDerror): minDerror = error
            if abs(error) > abs(maxDerror): maxDerror = error
            MSDerror += error
            MADerror += abs(error)
            RMSDerror += error*error
            count_rxn += 1
    tables += """\n   %s\n""" % (table_delimit)

    if count_rxn:
        tables += """%23s   %17s %8.4f\n""" % ('Minimal Dev', '', minDerror)
        tables += """%23s   %17s %8.4f\n""" % ('Maximal Dev', '', maxDerror)
        tables += """%23s   %17s %8.4f\n""" % ('Mean Signed Dev', '', MSDerror/float(count_rxn))
        tables += """%23s   %17s %8.4f\n""" % ('Mean Absolute Dev', '', MADerror/float(count_rxn))
        tables += """%23s   %17s %8.4f\n""" % ('RMS Dev', '', sqrt(RMSDerror/float(count_rxn)))
        tables += """   %s\n""" % (table_delimit)

        #print tables
        PsiMod.print_out(tables)
        return MADerror/float(count_rxn)

    else:
        return 0.0

def drop_duplicates(seq): 
    noDupes = []
    [noDupes.append(i) for i in seq if not noDupes.count(i)]
    return noDupes

def assemble_function2call(**kwargs):
    function2call = kwargs['func']
    return function2call(**kwargs)

def assemble_function2call_dropfunc(funcarg, **kwargs):
    function2call = funcarg
    return function2call(**kwargs)

def tblhead(tbl_maxrgt, tbl_delimit, ttype):
    tbl_str = ''
    tbl_str += """\n   %s""" % (tbl_delimit)
    if   ttype == 1: tbl_str += """\n%23s %19s %8s""" % ('Reaction', 'Reaction Energy', 'Error')
    elif ttype == 2: tbl_str += """\n%23s     %19s %4s""" % ('Reaction', 'Reaction Value', '')
    for i in range(tbl_maxrgt):
        tbl_str += """%20s""" % ('Reagent '+str(i+1))
    if   ttype == 1: tbl_str += """\n%23s %10s %17s""" % ('', 'Ref', '[kcal/mol]')
    elif ttype == 2: tbl_str += """\n%23s %10s %17s""" % ('', '', '')
    for i in range(tbl_maxrgt):
        if   ttype == 1: tbl_str += """%20s""" % ('[H] Wt')
        elif ttype == 2: tbl_str += """%20s""" % ('Value Wt')
    tbl_str += """\n   %s""" % (tbl_delimit)
    return tbl_str

##  Aliases  ##
db = database

#######################
##  End of Database  ##
#######################



###################################
##  Start of Complete Basis Set  ##
###################################

def complete_basis_set(name, **kwargs):

    # Wrap any positional arguments into kwargs (for intercalls among wrappers)
    if not('name' in kwargs) and name:
        kwargs['name'] = name
    if not('cbs_func' in kwargs):
        if ('func' in kwargs):
            kwargs['cbs_func'] = kwargs['func']
        else:
            kwargs['cbs_func'] = energy

    # Define some quantum chemical knowledge, namely what methods are subsumed in others
    VARH = {}
    VARH['scf']     = {'scftot'      : 'SCF ENERGY'                }
    VARH['mp2']     = {'scftot'      : 'SCF ENERGY',
                       'mp2corl'     : 'MP2 CORRELATION ENERGY'    }
    VARH['ccsd']    = {'scftot'      : 'SCF ENERGY',
                       'mp2corl'     : 'MP2 CORRELATION ENERGY',
                       'ccsdcorl'    : 'CCSD CORRELATION ENERGY'   }
    VARH['ccsd(t)'] = {'scftot'      : 'SCF ENERGY',
                       'mp2corl'     : 'MP2 CORRELATION ENERGY',
                       'ccsdcorl'    : 'CCSD CORRELATION ENERGY',
                       'ccsd(t)corl' : 'CCSD(T) CORRELATION ENERGY'}

    domax_scf = 1
    domax_xtpl = 0
    domax_delta = 0
    domax_delta2 = 0
    finalenergy = 0.0
    do_scf = 1
    do_xtpl = 0
    do_delta = 0
    do_delta2 = 0

    # Must collect (here) and set (below) basis sets after every new molecule activation
    user_basis = PsiMod.get_option('BASIS')
    #user_puream = PsiMod.get_option('PUREAM')
    #user_ri_basis_scf = PsiMod.get_option('RI_BASIS_SCF')
    #user_ri_basis_mp2 = PsiMod.get_option('RI_BASIS_MP2')
    #user_ri_basis_cc = PsiMod.get_option('RI_BASIS_CC')
    #user_ri_basis_sapt = PsiMod.get_option('RI_BASIS_SAPT')
    #user_ri_basis_elst = PsiMod.get_option('RI_BASIS_ELST')
    #user_reference = PsiMod.get_option('REFERENCE')

    # Establish method for correlation energy
    #cbs_corl_wfn = 'none'
    if(kwargs.has_key('name')):
        cbs_corl_wfn = kwargs['name'].lower()
        if re.match(r'^scf$', cbs_corl_wfn.lower()):
            pass
            #cbs_corl_wfn = 'none'
        else:
            domax_scf = 0
            domax_xtpl = 1
            do_xtpl = 1
    if(kwargs.has_key('corl_wfn')):
        cbs_corl_wfn = kwargs['corl_wfn'].lower()
        domax_scf = 0
        domax_xtpl = 1
        do_xtpl = 1
    if domax_xtpl:
        if not (cbs_corl_wfn in VARH.keys()):
            raise Exception('Requested CORL method \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_corl_wfn))

    # Establish method for delta correction energy
    #cbs_delta_wfn = 'none'
    if(kwargs.has_key('delta_wfn')):
        cbs_delta_wfn = kwargs['delta_wfn'].lower()
        domax_xtpl = 0
        domax_delta = 1
        do_delta = 1
    if domax_delta:
        if not (cbs_delta_wfn in VARH.keys()):
            raise Exception('Requested DELTA method \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta_wfn))

        if(kwargs.has_key('delta_wfn_lesser')):
            cbs_delta_wfn_lesser = kwargs['delta_wfn_lesser'].lower()
        else:
            cbs_delta_wfn_lesser = 'mp2'
        if not (cbs_delta_wfn_lesser in VARH.keys()):
            raise Exception('Requested DELTA method lesser \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta_wfn_lesser))

    # Establish method for second delta correction energy
    if(kwargs.has_key('delta2_wfn')):
        cbs_delta2_wfn = kwargs['delta2_wfn'].lower()
        domax_delta = 0
        domax_delta2 = 1
        do_delta2 = 1
    if domax_delta2:
        if not (cbs_delta2_wfn in VARH.keys()):
            raise Exception('Requested DELTA2 method \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta2_wfn))

        if(kwargs.has_key('delta2_wfn_lesser')):
            cbs_delta2_wfn_lesser = kwargs['delta2_wfn_lesser'].lower()
        else:
            cbs_delta2_wfn_lesser = 'mp2'
        if not (cbs_delta2_wfn_lesser in VARH.keys()):
            raise Exception('Requested DELTA2 method lesser \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta2_wfn_lesser))

    # Establish list of valid basis sets for correlation energy
    if(kwargs.has_key('corl_basis')):
        BSTC, ZETC = validate_bracketed_basis(kwargs['corl_basis'].lower())
    else:
        if not domax_scf:
            raise Exception('CORL basis sets through keyword \'%s\' are required.' % ('corl_basis'))

    # Establish list of valid basis sets for scf energy
    if(kwargs.has_key('scf_basis')):
        BSTR, ZETR = validate_bracketed_basis(kwargs['scf_basis'].lower())
    else:
        if domax_scf:
            raise Exception('SCF basis sets through keyword \'%s\' are required.' % ('scf_basis'))
        else:
            if(kwargs.has_key('corl_basis')):
                BSTR = BSTC[:]
                ZETR = ZETC[:]

    # Establish list of valid basis sets for delta correction energy
    if domax_delta or domax_delta2:
        if(kwargs.has_key('delta_basis')):
            BSTD, ZETD = validate_bracketed_basis(kwargs['delta_basis'].lower())
        else:
            raise Exception('DELTA basis sets through keyword \'%s\' are required.' % ('delta_basis'))

    # Establish list of valid basis sets for second delta correction energy
    if domax_delta2:
        if(kwargs.has_key('delta2_basis')):
            BSTD2, ZETD2 = validate_bracketed_basis(kwargs['delta2_basis'].lower())
        else:
            raise Exception('DELTA2 basis sets through keyword \'%s\' are required.' % ('delta2_basis'))

    # Establish treatment for scf energy
    cbs_scf_scheme = highest_1
    if(kwargs.has_key('scf_scheme')):
        cbs_scf_scheme = kwargs['scf_scheme']

    if cbs_scf_scheme:
        pass
    else:
        raise Exception('SCF extrapolation mode \'%s\' not valid.' % (cbs_scf_scheme))

    # Establish treatment for correlation energy
    cbs_corl_scheme = highest_1
    if(kwargs.has_key('corl_scheme')):
        cbs_corl_scheme = kwargs['corl_scheme']

    if cbs_corl_scheme:
        pass
    else:
        raise Exception('CORL extrapolation mode \'%s\' not valid.' % (cbs_corl_scheme))
  
    # Establish treatment for delta correction energy
    cbs_delta_scheme = highest_1
    if(kwargs.has_key('delta_scheme')):
        cbs_delta_scheme = kwargs['delta_scheme']

    if cbs_delta_scheme:
        pass
    else:
        raise Exception('DELTA extrapolation mode \'%s\' not valid.' % (cbs_delta_scheme))
  
    # Establish treatment for delta2 correction energy
    cbs_delta2_scheme = highest_1
    if(kwargs.has_key('delta2_scheme')):
        cbs_delta2_scheme = kwargs['delta2_scheme']

    if cbs_delta2_scheme:
        pass
    else:
        raise Exception('DELTA2 extrapolation mode \'%s\' not valid.' % (cbs_delta2_scheme))
  
    # build string of title banner
    cbsbanners = ''
    cbsbanners += """PsiMod.print_out('\\n')\n"""
    cbsbanners += """banner(' CBS Setup ')\n"""
    cbsbanners += """PsiMod.print_out('\\n')\n\n"""
    exec cbsbanners

    #print '\n'
    #if do_scf:    print 'SCF:    ', cbs_scf_scheme, ZETR, BSTR
    #if do_xtpl:   print 'CORL:   ', cbs_corl_scheme, ZETC, BSTC
    #if do_delta:  print 'DELTA:  ', cbs_delta_scheme, ZETD, BSTD
    #if do_delta2: print 'DELTA2: ', cbs_delta2_scheme, ZETD2, BSTD2
    #print '\n'

    d_fields = ['d_stage', 'd_scheme', 'd_need', 'd_coef', 'd_energy']
    f_fields = ['f_wfn', 'f_portion', 'f_basis', 'f_zeta', 'f_energy']
    GRAND_NEED = []
    MODELCHEM = []
    if do_scf:
        NEED = assemble_function2call_dropfunc(cbs_scf_scheme, 
            mode='requisition', basisname=BSTR, basiszeta=ZETR, wfnname='scf')
        GRAND_NEED.append(dict(zip(d_fields, ['scf', cbs_scf_scheme, NEED, +1, 0.0])))

    if do_xtpl:
        NEED = assemble_function2call_dropfunc(cbs_corl_scheme, 
            mode='requisition', basisname=BSTC, basiszeta=ZETC, wfnname=cbs_corl_wfn)
        GRAND_NEED.append(dict(zip(d_fields, ['corl', cbs_corl_scheme, NEED, +1, 0.0])))

    if do_delta:
        NEED = assemble_function2call_dropfunc(cbs_delta_scheme, 
            mode='requisition', basisname=BSTD, basiszeta=ZETD, wfnname=cbs_delta_wfn)
        GRAND_NEED.append(dict(zip(d_fields, ['delta', cbs_delta_scheme, NEED, +1, 0.0])))

        NEED = assemble_function2call_dropfunc(cbs_delta_scheme, 
            mode='requisition', basisname=BSTD, basiszeta=ZETD, wfnname=cbs_delta_wfn_lesser)
        GRAND_NEED.append(dict(zip(d_fields, ['delta', cbs_delta_scheme, NEED, -1, 0.0])))

    if do_delta2:
        NEED = assemble_function2call_dropfunc(cbs_delta2_scheme, 
            mode='requisition', basisname=BSTD2, basiszeta=ZETD2, wfnname=cbs_delta2_wfn)
        GRAND_NEED.append(dict(zip(d_fields, ['delta2', cbs_delta2_scheme, NEED, +1, 0.0])))

        NEED = assemble_function2call_dropfunc(cbs_delta2_scheme, 
            mode='requisition', basisname=BSTD2, basiszeta=ZETD2, wfnname=cbs_delta2_wfn_lesser)
        GRAND_NEED.append(dict(zip(d_fields, ['delta2', cbs_delta2_scheme, NEED, -1, 0.0])))

    for stage in GRAND_NEED:
        for lvl in stage['d_need'].iteritems():
            MODELCHEM.append(lvl[1])

    #print '\nGRAND_NEED', GRAND_NEED
    #print '\nMODELCHEM', MODELCHEM


    # Apply chemical reasoning to choose the minimum computations to run
    JOBS = MODELCHEM[:]

    instructions  = ''
    instructions += """    Naive listing of computations required.\n"""
    for mc in JOBS:
        instructions += """   %12s / %-24s for  %s\n""" % (mc['f_wfn'], mc['f_basis'], VARH[mc['f_wfn']][mc['f_wfn']+mc['f_portion']])

    #     Remove duplicate modelchem portion listings
    for indx_mc, mc in enumerate(MODELCHEM):
        dups = -1
        for indx_job, job in enumerate(JOBS):
            if (job['f_wfn'] == mc['f_wfn']) and (job['f_basis'] == mc['f_basis']):
                dups += 1
                if (dups >= 1):
                    del JOBS[indx_job]

    #     Remove chemically subsumed modelchem portion listings
    for indx_mc, mc in enumerate(MODELCHEM):
        for menial in VARH[mc['f_wfn']]:
            for indx_job, job in enumerate(JOBS):
                if (menial == job['f_wfn']+job['f_portion']) and (mc['f_basis'] == job['f_basis']) and not (mc['f_wfn'] == job['f_wfn']):
                    del JOBS[indx_job]

    instructions += """\n    Enlightened listing of computations required.\n"""
    for mc in JOBS:
        instructions += """   %12s / %-24s for  %s\n""" % (mc['f_wfn'], mc['f_basis'], VARH[mc['f_wfn']][mc['f_wfn']+mc['f_portion']])

    #     Expand listings to all that will be obtained
    JOBS_EXT = []
    for indx_job, job in enumerate(JOBS):
        for menial in VARH[job['f_wfn']]:
            temp_wfn, temp_portion = split_menial(menial)
            JOBS_EXT.append(dict(zip(f_fields, [temp_wfn, temp_portion, job['f_basis'], job['f_zeta'], 0.0])))

    instructions += """\n    Full listing of computations to be obtained (required and bonus).\n"""
    for mc in JOBS_EXT:
        instructions += """   %12s / %-24s for  %s\n""" % (mc['f_wfn'], mc['f_basis'], VARH[mc['f_wfn']][mc['f_wfn']+mc['f_portion']])
    #print instructions
    PsiMod.print_out(instructions)


    # Run necessary computations
    for mc in JOBS:
        kwargs['name'] = mc['f_wfn']

        # build string of title banner
        cbsbanners = ''
        cbsbanners += """PsiMod.print_out('\\n')\n"""
        cbsbanners += """banner(' CBS Computation: %s / %s ')\n""" % (mc['f_wfn'].upper(), mc['f_basis'].upper())
        cbsbanners += """PsiMod.print_out('\\n')\n\n"""
        exec cbsbanners

        # build string of molecule and commands that are dependent on the database
        commands  = '\n'
        commands += """\nPsiMod.set_global_option('BASIS', '%s')\n""" % (mc['f_basis'])
        exec commands

        # make energy() call
        lesserkwargs = kwargs.copy()
        del lesserkwargs['cbs_func']
        mc['f_energy'] = assemble_function2call_dropfunc(kwargs['cbs_func'], **lesserkwargs)

        # fill in energies for subsumed methods
        for menial in VARH[mc['f_wfn']]:
            temp_wfn, temp_portion = split_menial(menial)
            for job in JOBS_EXT:
                if (temp_wfn == job['f_wfn']) and (temp_portion == job['f_portion']) and (mc['f_basis'] == job['f_basis']):
                    job['f_energy'] = PsiMod.get_variable(VARH[temp_wfn][menial])

        PsiMod.clean()

    # Insert obtained energies into the array that stores the cbs stages
    for stage in GRAND_NEED:
        for lvl in stage['d_need'].iteritems():
            MODELCHEM.append(lvl[1])

            for job in JOBS_EXT:
                if (lvl[1]['f_wfn'] == job['f_wfn']) and (lvl[1]['f_portion'] == job['f_portion']) and \
                   (lvl[1]['f_basis'] == job['f_basis']) and (lvl[1]['f_zeta'] == job['f_zeta']):
                    lvl[1]['f_energy'] = job['f_energy']

    for stage in GRAND_NEED:
        stage['d_energy'] = assemble_function2call_dropfunc(stage['d_scheme'], needname = stage['d_need'], mode = 'evaluate')
        finalenergy += stage['d_energy'] * stage['d_coef']

    #print '\nJOBS', JOBS
    #print '\nJOBS_EXT', JOBS_EXT
    #print '\nGRAND_NEED', GRAND_NEED
    #print '\nMODELCHEM', MODELCHEM

    # build string of title banner
    cbsbanners = ''
    cbsbanners += """PsiMod.print_out('\\n')\n"""
    cbsbanners += """banner(' CBS Results ')\n"""
    cbsbanners += """PsiMod.print_out('\\n')\n\n"""
    exec cbsbanners

    # build string of results table
    table_delimit = '  ' + '-' * 90 + '\n'
    tables  = ''
    tables += table_delimit
    tables += """   %12s   %-18s %8s %16s  %-s\n""" % ('Method/Stage', 'Basis', 'Wt/Rqd', 'Energy [H]', 'Variable/Scheme')
    tables += table_delimit
    tables += """   Components\n"""
    for job in JOBS_EXT:
        star = ''
        for mc in MODELCHEM:
            if (job['f_wfn'] == mc['f_wfn']) and (job['f_basis'] == mc['f_basis']):
                star = '*'
        tables += """   %12s / %-24s %2s %16.8f  %s\n""" % (job['f_wfn'], job['f_basis'], star,
                  job['f_energy'], VARH[job['f_wfn']][job['f_wfn']+job['f_portion']])
    tables += table_delimit
    tables += """   Stages\n"""
    for stage in GRAND_NEED:
        tables += """   %12s   %24s %2d %16.8f  %-24s\n""" % (stage['d_stage'], '', stage['d_coef'], stage['d_energy'], stage['d_scheme'].__name__)
    tables += table_delimit
    tables += """   CBS\n"""
    tables += """   %42s %16.8f\n""" % ('', finalenergy)
    tables += table_delimit
     
    #print tables
    PsiMod.print_out(tables)

    # Restore global options that may be changed
    PsiMod.set_global_option('BASIS', user_basis)

    return finalenergy


# Transform and validate basis sets from 'cc-pV[Q5]Z' into [cc-pVQZ, cc-pV5Z] and [4, 5]
def validate_bracketed_basis(basisstring):

    ZETA = ['d', 't', 'q', '5', '6']
    BSET = []
    ZSET = []
    if re.match(r'.*cc-.*\[.*\]z$', basisstring, flags=re.IGNORECASE):
        basispattern = re.compile(r'^(.*)\[(.*)\](.*)$')
        basisname = basispattern.match(basisstring)
        for b in basisname.group(2):
            if b not in ZETA:
                raise Exception('Basis set \'%s\' has invalid zeta level \'%s\'.' % (basisstring, b))
            if len(ZSET) != 0:
                if (int(ZSET[len(ZSET)-1]) - ZETA.index(b)) != 1:
                    raise Exception('Basis set \'%s\' has out-of-order zeta level \'%s\'.' % (basisstring, b))
            BSET.append(basisname.group(1) + b + basisname.group(3))
            if b == 'd': b = '2'
            if b == 't': b = '3'
            if b == 'q': b = '4'
            ZSET.append(int(b))
    else:
        raise Exception('Basis set \'%s\' not valid.' % (basisstring))

    return [BSET, ZSET]


# Defining equation in LaTeX:  $E_{total}(\ell_{max}) =$ 
def highest_1(**largs):

    energypiece = 0.0
    functionname = sys._getframe().f_code.co_name
    f_fields = ['f_wfn', 'f_portion', 'f_basis', 'f_zeta', 'f_energy']
    [mode, NEED, wfnname, BSET, ZSET] = validate_scheme_args(functionname, **largs)

    if (mode == 'requisition'):

        # Impose restrictions on zeta sequence
        if (len(ZSET) == 0):
            raise Exception('Call to \'%s\' not valid with \'%s\' basis sets.' % (functionname, len(ZSET)))

        # Return array that logs the requisite jobs
        if (wfnname == 'scf'):
            portion = 'tot'
        else:
            portion = 'corl'
        NEED = {'HI' : dict(zip(f_fields, [wfnname, portion, BSET[len(ZSET)-1], ZSET[len(ZSET)-1], 0.0])) }

        return NEED 

    elif (mode == 'evaluate'):

        energypiece = NEED['HI']['f_energy']
        return energypiece

   
# Defining equation in LaTeX:  $E_{corl}^{X} = E_{corl}^{\infty} + \beta X^{-3}$
# Solution equation in LaTeX:  $E_{corl}^{\infty} = \frac{E_{corl}^{X} X^3 - E_{corl}^{X-1} (X-1)^3}{X^3 - (X-1)^3}$
# Solution equation in LaTeX:  $\beta = \frac{E_{corl}^{X} - E_{corl}^{X-1}}{X^{-3} - (X-1)^{-3}}$
def corl_xtpl_helgaker_2(**largs):

    energypiece = 0.0
    functionname = sys._getframe().f_code.co_name
    f_fields = ['f_wfn', 'f_portion', 'f_basis', 'f_zeta', 'f_energy']
    [mode, NEED, wfnname, BSET, ZSET] = validate_scheme_args(functionname, **largs)

    if (mode == 'requisition'):

        # Impose restrictions on zeta sequence
        if (len(ZSET) != 2):
            raise Exception('Call to \'%s\' not valid with \'%s\' basis sets.' % (functionname, len(ZSET)))

        # Return array that logs the requisite jobs
        NEED = {'HI' : dict(zip(f_fields, [wfnname, 'corl', BSET[1], ZSET[1], 0.0])),
                'LO' : dict(zip(f_fields, [wfnname, 'corl', BSET[0], ZSET[0], 0.0])) }

        return NEED

    elif (mode == 'evaluate'):

        # Extract required energies and zeta integers from array
        eHI = NEED['HI']['f_energy']
        zHI = NEED['HI']['f_zeta']
        eLO = NEED['LO']['f_energy']
        zLO = NEED['LO']['f_zeta']

        # Return extrapolated energy
        energypiece = (eHI * zHI**3 - eLO * zLO**3) / (zHI**3 - zLO**3) 
        beta = (eHI - eLO) / (zHI**(-3) - zLO**(-3))

        cbsscheme  = ''
        cbsscheme += """   %s\n""" % (functionname)
        cbsscheme += """   Energy Piece: %16.8f\n""" % (energypiece)
        cbsscheme += """   Beta:         %16.8f\n""" % (beta)
        cbsscheme += """\n"""
        PsiMod.print_out(cbsscheme)

        return energypiece


# Defining equation in LaTeX:  $E_{scf}(\ell_{max}) = E_{scf}^{\text{CBS}} + Ae^{b\ell_{max}}$
def scf_xtpl_helgaker_3(**largs):

    energypiece = 0.0
    functionname = sys._getframe().f_code.co_name
    f_fields = ['f_wfn', 'f_portion', 'f_basis', 'f_zeta', 'f_energy']
    [mode, NEED, wfnname, BSET, ZSET] = validate_scheme_args(functionname, **largs)

    if (mode == 'requisition'):

        # Impose restrictions on zeta sequence
        if (len(ZSET) != 3):
            raise Exception('Call to \'%s\' not valid with \'%s\' basis sets.' % (functionname, len(ZSET)))

        # Return array that logs the requisite jobs
        NEED = {'HI' : dict(zip(f_fields, [wfnname, 'tot', BSET[2], ZSET[2], 0.0])),
                'MD' : dict(zip(f_fields, [wfnname, 'tot', BSET[1], ZSET[1], 0.0])),
                'LO' : dict(zip(f_fields, [wfnname, 'tot', BSET[0], ZSET[0], 0.0])) }

        return NEED

    elif (mode == 'evaluate'):

        # Extract required energies and zeta integers from array
        eHI = NEED['HI']['f_energy']
        eMD = NEED['MD']['f_energy']
        eLO = NEED['LO']['f_energy']

        # Return extrapolated energy
        energypiece = (eHI * eLO - eMD * eMD) / (eHI + eLO - 2 * eMD)

        return energypiece


# Defining equation in LaTeX:  $E_{scf}(\ell_{max}) = E_{scf}^{\text{CBS}} + Ae^{b\ell_{max}}$
def scf_xtpl_helgaker_2(**largs):

    energypiece = 0.0
    functionname = sys._getframe().f_code.co_name
    f_fields = ['f_wfn', 'f_portion', 'f_basis', 'f_zeta', 'f_energy']
    [mode, NEED, wfnname, BSET, ZSET] = validate_scheme_args(functionname, **largs)

    if (mode == 'requisition'):

        # Impose restrictions on zeta sequence
        if (len(ZSET) != 2):
            raise Exception('Call to \'%s\' not valid with \'%s\' basis sets.' % (functionname, len(ZSET)))

        # Return array that logs the requisite jobs
        NEED = {'HI' : dict(zip(f_fields, [wfnname, 'tot', BSET[1], ZSET[1], 0.0])),
                'LO' : dict(zip(f_fields, [wfnname, 'tot', BSET[0], ZSET[0], 0.0])) }

        return NEED

    elif (mode == 'evaluate'):

        # Extract required energies and zeta integers from array
        eHI = NEED['HI']['f_energy']
        eLO = NEED['LO']['f_energy']

        alpha = 1.63
        if(largs.has_key('parameter')):
            alpha = largs['parameter']

        # Return extrapolated energy
        energypiece = (eHI - eLO * math.exp(-alpha)) / (1 - math.exp(-alpha))

        return energypiece


def validate_scheme_args(functionname, **largs):

    mode = ''
    NEED = []
    wfnname = ''
    BSET = []
    ZSET = []

    # Mode where function fills out a form NEED with the computations needed to fulfill its call
    if re.match(r'^requisition$', largs['mode'].lower()):
        mode = largs['mode'].lower()

        if(largs.has_key('wfnname')):
            wfnname = largs['wfnname']
        else:
            raise Exception('Call to \'%s\' has keyword \'wfnname\' missing.' % (functionname))

        if(largs.has_key('basisname')):
            BSET = largs['basisname']
        else:
            raise Exception('Call to \'%s\' has keyword \'basisname\' missing.' % (functionname))

        if(largs.has_key('basiszeta')):
            ZSET = largs['basiszeta']
        else:
            raise Exception('Call to \'%s\' has keyword \'basiszeta\' missing.' % (functionname))

    # Mode where function reads the now-filled-in energies from that same form and performs the sp, xtpl, delta, etc.
    elif re.match(r'^evaluate$', largs['mode'].lower()):
        mode = largs['mode'].lower()

        if(largs.has_key('needname')):
            NEED = largs['needname']
        else:
            raise Exception('Call to \'%s\' has keyword \'needname\' missing.' % (functionname))

    else:
        raise Exception('Call to \'%s\' has keyword \'mode\' missing or invalid.' % (functionname))

    return [mode, NEED, wfnname, BSET, ZSET]


def split_menial(menial):

    PTYP = ['tot', 'corl']
    for temp in PTYP:
        if menial.endswith(temp):
            temp_wfn = menial[:-len(temp)]
            temp_portion = temp

    return [temp_wfn, temp_portion]


##  Aliases  ##
cbs = complete_basis_set

#################################
##  End of Complete Basis Set  ##
#################################

