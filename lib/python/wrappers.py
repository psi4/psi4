"""Module with functions that call the four main :py:mod:`driver`
functions: :py:mod:`driver.energy`, :py:mod:`driver.optimize`,
:py:mod:`driver.response`, and :py:mod:`driver.frequency`.

"""


import PsiMod
import re
import os
import input
import math
import warnings
import pickle
import copy
import physconst
from driver import *
from molutil import *
from text import *
from collections import defaultdict
from procutil import *


# Function to make calls among wrappers(), energy(), optimize(), etc.
def call_function_in_1st_argument(funcarg, **largs):
    """Function to make primary function call to energy(), opt(), etc.
    with options dictionary *largs*.
    Useful when *funcarg* to call is stored in variable.

    """
    return funcarg(**largs)


#######################
##  Start of n_body  ##
#######################

def n_body(name, **kwargs):
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # Wrap any positional arguments into kwargs (for intercalls among wrappers)
    if not('name' in kwargs) and name:
        kwargs['name'] = name.lower()

    # Establish function to call
    if not('n_body_func' in kwargs):
        if ('func' in kwargs):
            kwargs['n_body_func'] = kwargs['func']
            del kwargs['func']
        else:
            kwargs['n_body_func'] = energy
    func = kwargs['n_body_func']
    if not func:
        raise ValidationError('Function \'%s\' does not exist to be called by wrapper n_body.' % (func.__name__))
    if (func is db):
        raise ValidationError('Wrapper n_body is unhappy to be calling function \'%s\'.' % (func.__name__))

    # Make sure the molecule the user provided is the active one
    if 'molecule' in kwargs:
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    PsiMod.set_global_option("BASIS", PsiMod.get_global_option("BASIS"))

    # N-body run configuration
    bsse = 'on'
    if 'bsse' in kwargs:
        bsse = kwargs['bsse']

    max_n_body = molecule.nfragments()
    if 'max_n_body' in kwargs:
        max_n_body = kwargs['max_n_body']

    do_total = False
    if 'do_total' in kwargs:
        do_total = kwargs['do_total']

    external = None
    external_indices = []
    if 'external' in kwargs:
        external = kwargs['external']
        external_indices = [molecule.nfragments()]
        if 'external_monomers' in kwargs:
            external_indices = kwargs['external_monomers']

    # Check input args
    if not bsse == 'off' and not bsse == 'on' and not bsse == 'both':
        raise ValidationError('n_body: bsse argument is one of on, off, or both')
    if max_n_body < 1:
        raise ValidationError('n_body: max_n_body must be at least 1')
    if max_n_body > molecule.nfragments():
        raise ValidationError('n_body: max_n_body must be <= to the number of fragments in the molecule')

    # Set to save RI integrals for repeated full-basis computations
    ri_ints_io = PsiMod.get_option('DF_INTS_IO')
    PsiMod.set_global_option('DF_INTS_IO', 'SAVE')
    psioh = PsiMod.IOManager.shared_object()
    psioh.set_specific_retention(97, True)

    # Tell 'em what you're gonna tell 'em
    has_external = 'No'
    if (external):
        has_external = 'Yes'
    PsiMod.print_out('\n')
    PsiMod.print_out('    ==> N-Body Interaction Energy Analysis <==\n\n')
    PsiMod.print_out('        BSSE Treatment:              %s\n' % (bsse))
    PsiMod.print_out('        Maximum N-Body Interactions: %d\n' % (max_n_body))
    PsiMod.print_out('        Compute Total Energy:        %s\n' % (do_total))
    PsiMod.print_out('        External Field:              %s\n' % (has_external))
    if (external):
        PsiMod.print_out('        External Field Monomers:     ')
        for k in external_indices:
            PsiMod.print_out('%-3d ' % (k))
        PsiMod.print_out('\n')
    PsiMod.print_out('\n')

    # Run the total molecule, if required
    energies_full = {}
    energies_mon = {}
    N = molecule.nfragments()
    Etotal = 0.0
    if do_total or max_n_body == molecule.nfragments():
        PsiMod.print_out('    => Total Cluster Energy <=\n')
        # Full cluster always gets the external field
        if (external):
            PsiMod.set_global_option_python("EXTERN", external)
        Etotal = call_function_in_1st_argument(func, **kwargs)
        if (external):
            PsiMod.set_global_option_python("EXTERN", None)
        energies_full[N] = []
        energies_full[N].append(Etotal)
        energies_mon[N] = []
        energies_mon[N].append(Etotal)
        PsiMod.set_global_option('DF_INTS_IO', 'LOAD')
        PsiMod.clean()

    max_effective = max_n_body
    if (max_effective == N):
        max_effective = N - 1

    # Build the combos for indexing purposes
    Ns = []
    if (max_n_body == N or do_total):
        Ns.append(N)
    for n in range(max_effective, 0, -1):
        Ns.append(n)

    combos = {}
    for n in Ns:

        combos[n] = []

        # Loop through combinations in lexical order #

        # initialize the reals list
        reals = []
        #setup first combination [3,2,1] lexical ordering
        #fragments indexing is 1's based, bloody hell
        for index in range(n, 0, -1):
            reals.append(index)
        #start loop through lexical promotion
        counter = 0
        while True:

            counter = counter + 1

            # Append the current combo
            combos[n].append(copy.deepcopy(reals))

            #reset rank
            rank = 0

            #look for lexical promotion opportunity
            #i.e.: [4 2 1] has a promotion opportunity at
            # index 1 to produce [4 3 1]
            for k in range(n - 2, -1, -1):
                if (reals[k] != reals[k + 1] + 1):
                    rank = k + 1
                    break

            #do the promotion
            reals[rank] = reals[rank] + 1

            #demote the right portion of the register
            val = 1
            for k in range(n - 1, rank, -1):
                reals[k] = val
                val = val + 1

            #boundary condition is promotion into
            #[nfrag+1 nfrag-1 ...]
            if (reals[0] > N):
                break

    # Hack for external
    externNone = PsiMod.ExternalPotential()

    # Run the clusters in the full basis
    if bsse == 'on' or bsse == 'both':
        for n in range(max_effective, 0, -1):
            energies_full[n] = []
            clusters = extract_clusters(molecule, True, n)
            for k in range(len(clusters)):
                activate(clusters[k])
                # Do the external field for this cluster or not?
                if (external):
                    do_extern = False
                    for mon in combos[n][k]:
                        if (mon in external_indices):
                            do_extern = True
                            break
                    if do_extern:
                        PsiMod.set_global_option_python("EXTERN", external)
                PsiMod.print_out('\n    => Cluster (N-Body %4d, Combination %4d) Energy (Full Basis) <=\n' % (n, k + 1))
                energies_full[n].append(call_function_in_1st_argument(func, **kwargs))
                # Turn the external field off
                if (external):
                    PsiMod.set_global_option_python("EXTERN", externNone)
                PsiMod.set_global_option('DF_INTS_IO', 'LOAD')
                PsiMod.clean()

    # Run the clusters in the minimal cluster bases
    PsiMod.set_global_option('DF_INTS_IO', 'NONE')
    if bsse == 'off' or bsse == 'both':
        for n in range(max_effective, 0, -1):
            energies_mon[n] = []
            clusters = extract_clusters(molecule, False, n)
            for k in range(len(clusters)):
                activate(clusters[k])
                # Do the external field for this cluster or not?
                if (external):
                    do_extern = False
                    for mon in combos[n][k]:
                        if (mon in external_indices):
                            do_extern = True
                            break
                    if do_extern:
                        PsiMod.set_global_option_python("EXTERN", external)
                PsiMod.print_out('\n    => Cluster (N-Body %4d, Combination %4d) Energy (Cluster Basis) <=\n' % (n, k + 1))
                energies_mon[n].append(call_function_in_1st_argument(func, **kwargs))
                # Turn the external field off
                if (external):
                    PsiMod.set_global_option_python("EXTERN", externNone)
                PsiMod.clean()

    # Report the energies
    PsiMod.print_out('\n    ==> N-Body Interaction Energy Analysis: Combination Definitions <==\n\n')

    PsiMod.print_out('     %6s %6s | %-24s\n' % ("N-Body", "Combo", "Monomers"))
    for n in Ns:
        for k in range(len(combos[n])):
            PsiMod.print_out('     %6d %6d | ' % (n, k + 1))
            for l in combos[n][k]:
                PsiMod.print_out('%-3d ' % (l))
            PsiMod.print_out('\n')
    PsiMod.print_out('\n')

    PsiMod.print_out('    ==> N-Body Interaction Energy Analysis: Total Energies <==\n\n')

    if bsse == 'on' or bsse == 'both':
        PsiMod.print_out('     => Full Basis Set Results <=\n\n')
        PsiMod.print_out('     %6s %6s %24s %24s\n' % ("N-Body", "Combo", "E [H]", "E [kcal mol^-1]"))
        for n in Ns:
            for k in range(len(energies_full[n])):
                PsiMod.print_out('     %6d %6d %24.16E %24.16E\n' % (n, k + 1, energies_full[n][k],
                   physconst.psi_hartree2kcalmol * energies_full[n][k]))
        PsiMod.print_out('\n')

    if bsse == 'off' or bsse == 'both':
        PsiMod.print_out('     => Cluster Basis Set Results <=\n\n')
        PsiMod.print_out('     %6s %6s %24s %24s\n' % ("N-Body", "Combo", "E [H]", "E [kcal mol^-1]"))
        for n in Ns:
            for k in range(len(energies_mon[n])):
                PsiMod.print_out('     %6d %6d %24.16E %24.16E\n' % (n, k + 1, energies_mon[n][k],
                   physconst.psi_hartree2kcalmol * energies_mon[n][k]))
        PsiMod.print_out('\n')

    if bsse == 'both':
        PsiMod.print_out('     => BSSE Results <=\n\n')
        PsiMod.print_out('     %6s %6s %24s %24s\n' % ("N-Body", "Combo", "Delta E [H]", "Delta E [kcal mol^-1]"))
        for n in Ns:
            for k in range(len(energies_mon[n])):
                PsiMod.print_out('     %6d %6d %24.16E %24.16E\n' % (n, k + 1, energies_full[n][k] - energies_mon[n][k],
                   physconst.psi_hartree2kcalmol * (energies_full[n][k] - energies_mon[n][k])))
        PsiMod.print_out('\n')

    PsiMod.print_out('    ==> N-Body Interaction Energy Analysis: N-Body Energies <==\n\n')

    if bsse == 'on' or bsse == 'both':
        PsiMod.print_out('     => Full Basis Set Results <=\n\n')
        PsiMod.print_out('     %6s %6s %24s %24s\n' % ("N-Body", "Combo", "E [H]", "E [kcal mol^-1]"))
        energies_n_full = {}
        for n in Ns:
            if n == 1:
                continue
            En = 0.0
            for k in range(len(energies_full[n])):
                E = energies_full[n][k]
                for l in range(len(combos[n][k])):
                    E -= energies_full[1][combos[n][k][l] - 1]
                PsiMod.print_out('     %6d %6d %24.16E %24.16E\n' % (n, k + 1, E, physconst.psi_hartree2kcalmol * E))
                En += E
            energies_n_full[n] = En
        for n in Ns:
            if n == 1:
                continue
            nn = molecule.nfragments() - 2
            kk = n - 2
            energies_n_full[n] /= (math.factorial(nn) / (math.factorial(kk) * math.factorial(nn - kk)))
            PsiMod.print_out('     %6d %6s %24.16E %24.16E\n' % (n, 'Total', energies_n_full[n],
               physconst.psi_hartree2kcalmol * energies_n_full[n]))
        PsiMod.print_out('\n')

    if bsse == 'off' or bsse == 'both':
        PsiMod.print_out('     => Cluster Basis Set Results <=\n\n')
        PsiMod.print_out('     %6s %6s %24s %24s\n' % ("N-Body", "Combo", "E [H]", "E [kcal mol^-1]"))
        energies_n_mon = {}
        for n in Ns:
            if n == 1:
                continue
            En = 0.0
            for k in range(len(energies_mon[n])):
                E = energies_mon[n][k]
                for l in range(len(combos[n][k])):
                    E -= energies_mon[1][combos[n][k][l] - 1]
                PsiMod.print_out('     %6d %6d %24.16E %24.16E\n' % (n, k + 1, E, physconst.psi_hartree2kcalmol * E))
                En += E
            energies_n_mon[n] = En
        for n in Ns:
            if n == 1:
                continue
            nn = molecule.nfragments() - 2
            kk = n - 2
            energies_n_mon[n] /= (math.factorial(nn) / (math.factorial(kk) * math.factorial(nn - kk)))
            PsiMod.print_out('     %6d %6s %24.16E %24.16E\n' % (n, 'Total', energies_n_mon[n],
               physconst.psi_hartree2kcalmol * energies_n_mon[n]))
        PsiMod.print_out('\n')

    if bsse == 'both':
        PsiMod.print_out('     => BSSE Results <=\n\n')
        PsiMod.print_out('     %6s %6s %24s %24s\n' % ("N-Body", "Combo", "Delta E [H]", "Delta E [kcal mol^-1]"))
        energies_n_bsse = {}
        for n in Ns:
            if n == 1:
                continue
            En = 0.0
            for k in range(len(energies_mon[n])):
                E = energies_full[n][k] - energies_mon[n][k]
                for l in range(len(combos[n][k])):
                    E -= energies_full[1][combos[n][k][l] - 1]
                    E += energies_mon[1][combos[n][k][l] - 1]
                PsiMod.print_out('     %6d %6d %24.16E %24.16E\n' % (n, k + 1, E, physconst.psi_hartree2kcalmol * E))
                En += E
            energies_n_bsse[n] = En
        for n in Ns:
            if n == 1:
                continue
            nn = molecule.nfragments() - 2
            kk = n - 2
            energies_n_bsse[n] /= (math.factorial(nn) / (math.factorial(kk) * math.factorial(nn - kk)))
            PsiMod.print_out('     %6d %6s %24.16E %24.16E\n' % (n, 'Total', energies_n_bsse[n],
               physconst.psi_hartree2kcalmol * energies_n_bsse[n]))
        PsiMod.print_out('\n')

    PsiMod.print_out('    ==> N-Body Interaction Energy Analysis: Non-Additivities <==\n\n')

    if bsse == 'on' or bsse == 'both':
        energies_n_full[1] = 0.0
        PsiMod.print_out('     => Full Basis Set Results <=\n\n')
        PsiMod.print_out('     %6s %24s %24s\n' % ("N-Body", "E [H]", "E [kcal mol^-1]"))
        for k in range(len(Ns)):
            n = Ns[k]
            if n == 1:
                continue
            E = energies_n_full[Ns[k]] - energies_n_full[Ns[k + 1]]
            PsiMod.print_out('     %6s %24.16E %24.16E\n' % (n, E, physconst.psi_hartree2kcalmol * E))
        PsiMod.print_out('\n')

    if bsse == 'off' or bsse == 'both':
        energies_n_mon[1] = 0.0
        PsiMod.print_out('     => Cluster Basis Set Results <=\n\n')
        PsiMod.print_out('     %6s %24s %24s\n' % ("N-Body", "E [H]", "E [kcal mol^-1]"))
        for k in range(len(Ns)):
            n = Ns[k]
            if n == 1:
                continue
            E = energies_n_mon[Ns[k]] - energies_n_mon[Ns[k + 1]]
            PsiMod.print_out('     %6s %24.16E %24.16E\n' % (n, E, physconst.psi_hartree2kcalmol * E))
        PsiMod.print_out('\n')

    if bsse == 'both':
        energies_n_bsse[1] = 0.0
        PsiMod.print_out('     => BSSE Results <=\n\n')
        PsiMod.print_out('     %6s %24s %24s\n' % ("N-Body", "Delta E [H]", "Delta E [kcal mol^-1]"))
        for k in range(len(Ns)):
            n = Ns[k]
            if n == 1:
                continue
            E = energies_n_bsse[Ns[k]] - energies_n_bsse[Ns[k + 1]]
            PsiMod.print_out('     %6s %24.16E %24.16E\n' % (n, E, physconst.psi_hartree2kcalmol * E))
        PsiMod.print_out('\n')

    # Put everything back the way it was
    PsiMod.set_global_option('DF_INTS_IO', ri_ints_io)
    psioh.set_specific_retention(97, False)
    PsiMod.clean()
    activate(molecule)

    if bsse == 'on' or bsse == 'both':
        return energies_n_full[Ns[0]]
    else:
        return energies_n_mon[Ns[0]]

##  Aliases  ##
nbody = n_body

#####################
##  End of n_body  ##
#####################


###################
##  Start of cp  ##
###################

def cp(name, **kwargs):
    """The cp function computes counterpoise-corrected two-body interaction energies
    for complexes composed of arbitrary numbers of monomers.

    :aliases: counterpoise_correct(), counterpoise_correction()

    :returns: (*float*) Counterpoise-corrected interaction energy in kcal/mol

    :PSI variables:

    .. envvar:: CP-CORRECTED 2-BODY INTERACTION ENERGY
        UNCP-CORRECTED 2-BODY INTERACTION ENERGY

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - No values of func besides energy have been tested.

       - Table print-out needs improving. Add some PSI variables.

    **Keywords**

    :type name: string
    :param name: ``'scf'`` || ``'ccsd(t)'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the molecule. May be any valid argument to
        :py:func:`driver.energy`; however, SAPT is not appropriate.

    :type func: function
    :param func: |dl| ``energy`` |dr| || ``optimize`` || ``cbs``

        Indicates the type of calculation to be performed on the molecule
        and each of its monomers. The default performs a single-point
        ``energy('name')``, while ``optimize`` perfoms a geometry optimization
        on each system, and ``cbs`` performs a compound single-point energy.
        If a nested series of python functions is intended
        (see `Function Intercalls`_), use keyword ``cp_func`` instead of ``func``.

    :type check_bsse: bool
    :param check_bsse: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether to additionally compute un-counterpoise corrected
        monomers and thus obtain an estimate for the basis set superposition error.

    **Examples**

    >>> # [1] counterpoise-corrected mp2 interaction energy
    >>> cp('dfmp2')

    """

    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # Wrap any positional arguments into kwargs (for intercalls among wrappers)
    if not('name' in kwargs) and name:
        kwargs['name'] = name.lower()

    # Establish function to call
    if not('cp_func' in kwargs):
        if ('func' in kwargs):
            kwargs['cp_func'] = kwargs['func']
            del kwargs['func']
        else:
            kwargs['cp_func'] = energy
    func = kwargs['cp_func']
    if not func:
        raise ValidationError('Function \'%s\' does not exist to be called by wrapper counterpoise_correct.' % (func.__name__))
    if (func is db):
        raise ValidationError('Wrapper counterpoise_correct is unhappy to be calling function \'%s\'.' % (func.__name__))

    if 'check_bsse' in kwargs and input.yes.match(str(kwargs['check_bsse'])):
        check_bsse = True
    else:
        check_bsse = False

    # Make sure the molecule the user provided is the active one
    if 'molecule' in kwargs:
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    PsiMod.set_global_option("BASIS", PsiMod.get_global_option("BASIS"))

    df_ints_io = PsiMod.get_option('DF_INTS_IO')
    PsiMod.set_global_option('DF_INTS_IO', 'SAVE')
    psioh = PsiMod.IOManager.shared_object()
    psioh.set_specific_retention(97, True)

    activate(molecule)
    molecule.update_geometry()

    PsiMod.print_out("\n")
    banner("CP Computation: Complex.\nFull Basis Set.")
    PsiMod.print_out("\n")
    e_dimer = call_function_in_1st_argument(func, **kwargs)
    #e_dimer = energy(name, **kwargs)

    PsiMod.clean()
    PsiMod.set_global_option('DF_INTS_IO', 'LOAD')

    # All monomers with ghosts
    monomers = extract_clusters(molecule, True, 1)
    e_monomer_full = []

    cluster_n = 0
    for cluster in monomers:
        activate(cluster)
        PsiMod.print_out("\n")
        banner(("CP Computation: Monomer %d.\n Full Basis Set." % (cluster_n + 1)))
        PsiMod.print_out("\n")
        e_monomer_full.append(call_function_in_1st_argument(func, **kwargs))
        #e_monomer_full.append(energy(name,**kwargs))
        cluster_n = cluster_n + 1
        PsiMod.clean()

    PsiMod.set_global_option('DF_INTS_IO', 'NONE')
    if (check_bsse):
        # All monomers without ghosts
        monomers = extract_clusters(molecule, False, 1)
        e_monomer_bsse = []

        cluster_n = 0
        for cluster in monomers:
            activate(cluster)
            PsiMod.print_out("\n")
            #cluster.print_to_output()
            banner(("CP Computation: Monomer %d.\n Monomer Set." % (cluster_n + 1)))
            PsiMod.print_out("\n")
            e_monomer_bsse.append(call_function_in_1st_argument(func, **kwargs))
            #e_monomer_bsse.append(energy(name,**kwargs))
            cluster_n = cluster_n + 1

    PsiMod.set_global_option('DF_INTS_IO', df_ints_io)
    psioh.set_specific_retention(97, False)

    activate(molecule)

    if (check_bsse == False):
        cp_table = Table(rows=["System:"], cols=["Energy (full):"])
        cp_table["Complex"] = [e_dimer]
        for cluster_n in range(0, len(monomers)):
            key = "Monomer %d" % (cluster_n + 1)
            cp_table[key] = [e_monomer_full[cluster_n]]

        e_full = e_dimer
        for cluster_n in range(0, len(monomers)):
            e_full = e_full - e_monomer_full[cluster_n]
        cp_table["Interaction"] = [e_full]

        PsiMod.set_variable('CP-CORRECTED 2-BODY INTERACTION ENERGY', e_full)

    else:
        cp_table = Table(rows=["System:"], cols=["Energy (full):", "Energy (monomer):", "BSSE:"])
        cp_table["Complex"] = [e_dimer, 0.0, 0.0]
        for cluster_n in range(0, len(monomers)):
            key = "Monomer %d" % (cluster_n + 1)
            cp_table[key] = [e_monomer_full[cluster_n], e_monomer_bsse[cluster_n], \
                e_monomer_full[cluster_n] - e_monomer_bsse[cluster_n]]

        e_full = e_dimer
        e_bsse = e_dimer
        for cluster_n in range(0, len(monomers)):
            e_full = e_full - e_monomer_full[cluster_n]
            e_bsse = e_bsse - e_monomer_bsse[cluster_n]
        cp_table["Totals:"] = [e_full, e_bsse, e_full - e_bsse]

        PsiMod.set_variable('UNCP-CORRECTED 2-BODY INTERACTION ENERGY', e_full)

    PsiMod.print_out("\n")
    banner("CP Computation: Results.")
    PsiMod.print_out("\n")

    banner("Hartree", 2)
    PsiMod.print_out("\n")

    PsiMod.print_out(str(cp_table))

    PsiMod.print_out("\n")
    banner("kcal*mol^-1", 2)
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


#########################
##  Start of Database  ##
#########################

def database(name, db_name, **kwargs):
    """Function to access the molecule objects and reference energies of
    popular chemical databases.

    :aliases: db()

    :returns: (*float*) Mean absolute deviation of the database in kcal/mol

    :PSI variables:

    .. envvar:: db_name DATABASE MEAN SIGNED DEVIATION
        db_name DATABASE MEAN ABSOLUTE DEVIATION
        db_name DATABASE ROOT-MEAN-SQUARE DEVIATION

    .. note:: It is very easy to make a database from a collection of xyz files
        using the script ``$PSIDATADIR/databases/ixyz2database.pl``.
        See `Creating a New Database`_ for details.

    .. caution:: Some features are not yet implemented. Buy a developer some coffee.

       - In sow/reap mode, use only global options (e.g., the local option set by ``set scf scf_type df`` will not be respected).

    **Keywords**

    :type name: string
    :param name: ``'scf'`` || ``'sapt0'`` || ``'ccsd(t)'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the database. May be any valid argument to
        :py:func:`driver.energy`.

    :type db_name: string
    :param db_name: ``'BASIC'`` || ``'S22'`` || ``'HTBH'`` || etc.

        Second argument, usually unlabeled. Indicates the requested database
        name, matching the name of a python file in ``psi4/lib/databases``.
        Consult that directory for available databases and literature citations.

    :type func: function
    :param func: |dl| ``energy`` |dr| || ``optimize`` || ``cbs``

        Indicates the type of calculation to be performed on each database
        member. The default performs a single-point ``energy('name')``, while
        ``optimize`` perfoms a geometry optimization on each reagent, and
        ``cbs`` performs a compound single-point energy. If a nested series
        of python functions is intended (see `Function Intercalls`_), use
        keyword ``db_func`` instead of ``func``.

    :type mode: string
    :param mode: |dl| ``'continuous'`` |dr| || ``'sow'`` || ``'reap'``

        Indicates whether the calculations required to complete the
        database are to be run in one file (``'continuous'``) or are to be
        farmed out in an embarrassingly parallel fashion
        (``'sow'``/``'reap'``).  For the latter, run an initial job with
        ``'sow'`` and follow instructions in its output file.

    :type cp: bool
    :param cp: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether counterpoise correction is employed in computing
        interaction energies. Use this option and NOT the :py:func:`wrappers.cp`
        function for BSSE correction in database().  Option
        valid only for databases consisting of bimolecular complexes.

    :type rlxd: bool
    :param rlxd: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether correction for deformation energy is
        employed in computing interaction energies.  Option valid only for
        databases consisting of bimolecular complexes with non-frozen
        monomers, e.g., HBC6.

    :type symm: bool
    :param symm: |dl| ``'on'`` |dr| || ``'off'``

        Indicates whether the native symmetry of the database reagents is
        employed (``'on'``) or whether it is forced to :math:`C_1` symmetry
        (``'off'``). Some computational methods (e.g., SAPT) require no
        symmetry, and this will be set by database().

    :type zpe: bool
    :param zpe: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether zero-point-energy corrections are appended to
        single-point energy values. Option valid only for certain
        thermochemical databases. Disabled until Hessians ready.

    :type benchmark: string
    :param benchmark: |dl| ``'default'`` |dr| || ``'S22A'`` || etc.

        Indicates whether a non-default set of reference energies, if
        available, are employed for the calculation of error statistics.

    :type tabulate: array of strings
    :param tabulate: |dl| ``[]`` |dr| || ``['scf total energy', 'natom']`` || etc.

        Indicates whether to form tables of variables other than the
        primary requested energy.  Available for any PSI variable.

    :type subset: string or array of strings
    :param subset:

        Indicates a subset of the full database to run. This is a very
        flexible option and can be used in three distinct ways, outlined
        below. Note that two take a string and the last takes an array.

        * ``'small'`` || ``'large'`` || ``'equilibrium'``
            Calls predefined subsets of the requested database, either
            ``'small'``, a few of the smallest database members,
            ``'large'``, the largest of the database members, or
            ``'equilibrium'``, the equilibrium geometries for a database
            composed of dissociation curves.
        * ``'BzBz_S'`` || ``'FaOOFaON'`` || ``'ArNe'`` || etc.
            For databases composed of dissociation curves, individual
            curves can be called by name. Consult the database python
            files for available molecular systems.  The choices for this
            keyword are case sensitive and must match the database python file
        * ``[1,2,5]`` || ``['1','2','5']`` || ``['BzMe-3.5', 'MeMe-5.0']`` || etc.
            Specify a list of database members to run. Consult the
            database python files for available molecular systems.  The
            choices for this keyword are case sensitive and must match the
            database python file

    **Examples**

    >>> # [1] Two-stage SCF calculation on short, equilibrium, and long helium dimer
    >>> db('scf','RGC10',cast_up='sto-3g',subset=['HeHe-0.85','HeHe-1.0','HeHe-1.5'], tabulate=['scf total energy','natom'])

    >>> # [2] Counterpoise-corrected interaction energies for three complexes in S22
    >>> #     Error statistics computed wrt an old benchmark, S22A
    >>> database('dfmp2','S22',cp=1,subset=[16,17,8],benchmark='S22A')

    >>> # [3] SAPT0 on the neon dimer dissociation curve
    >>> db('sapt0',subset='NeNe',cp=0,symm=0,db_name='RGC10')

    >>> # [4] Optimize system 1 in database S22, producing tables of scf and mp2 energy
    >>> db('mp2','S22',db_func=optimize,subset=[1], tabulate=['mp2 total energy','current energy'])

    >>> # [5] CCSD on the smallest systems of HTBH, a hydrogen-transfer database
    >>> database('ccsd','HTBH',subset='small', tabulate=['ccsd total energy', 'mp2 total energy'])

    """

    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)
    #hartree2kcalmol = 627.509469  # consistent with perl SETS scripts

    # Wrap any positional arguments into kwargs (for intercalls among wrappers)
    if not('name' in kwargs) and name:
        kwargs['name'] = name.lower()
    if not('db_name' in kwargs) and db_name:
        kwargs['db_name'] = db_name

    # Establish function to call
    if not('db_func' in kwargs):
        if ('func' in kwargs):
            kwargs['db_func'] = kwargs['func']
            del kwargs['func']
        else:
            kwargs['db_func'] = energy
    func = kwargs['db_func']
    if not func:
        raise ValidationError('Function \'%s\' does not exist to be called by wrapper database.' % (func.__name__))
    if (func is cp):
        raise ValidationError('Wrapper database is unhappy to be calling function \'%s\'. Use the cp keyword within database instead.' % (func.__name__))

    # Define path and load module for requested database
    sys.path.append('%sdatabases' % (PsiMod.Process.environment["PSIDATADIR"]))
    sys.path.append('%s/lib/databases' % PsiMod.psi_top_srcdir())
    try:
        database = __import__(db_name)
    except ImportError:
        PsiMod.print_out('\nPython module for database %s failed to load\n\n' % (db_name))
        PsiMod.print_out('\nSearch path that was tried:\n')
        PsiMod.print_out(", ".join(map(str, sys.path)))
        raise ValidationError("Python module loading problem for database " + str(db_name))
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
    user_df_basis_scf = PsiMod.get_option('DF_BASIS_SCF')
    user_df_basis_mp2 = PsiMod.get_option('DF_BASIS_MP2')
    user_df_basis_cc = PsiMod.get_option('DF_BASIS_CC')
    user_df_basis_sapt = PsiMod.get_option('DF_BASIS_SAPT')
    user_df_basis_elst = PsiMod.get_option('DF_BASIS_ELST')

    b_user_reference = PsiMod.has_global_option_changed('REFERENCE')
    user_reference = PsiMod.get_option('REFERENCE')
    user_memory = PsiMod.get_memory()

    user_molecule = PsiMod.get_active_molecule()

    # Configuration based upon e_name & db_name options
    #   Force non-supramolecular if needed
    symmetry_override = 0
    if re.match(r'^sapt', lowername):
        try:
            database.ACTV_SA
        except AttributeError:
            raise ValidationError('Database %s not suitable for non-supramolecular calculation.' % (db_name))
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
    if 'symm' in kwargs:
        db_symm = kwargs['symm']

    if input.no.match(str(db_symm)):
        symmetry_override = 1
    elif input.yes.match(str(db_symm)):
        pass
    else:
        raise ValidationError('Symmetry mode \'%s\' not valid.' % (db_symm))

    #   Option mode of operation- whether db run in one job or files farmed out
    #db_mode = 'continuous'
    #if(kwargs.has_key('mode')):
    #    db_mode = kwargs['mode']
    if not('db_mode' in kwargs):
        if ('mode' in kwargs):
            kwargs['db_mode'] = kwargs['mode']
            del kwargs['mode']
        else:
            kwargs['db_mode'] = 'continuous'
    db_mode = kwargs['db_mode']

    if (db_mode.lower() == 'continuous'):
        pass
    elif (db_mode.lower() == 'sow'):
        pass
    elif (db_mode.lower() == 'reap'):
        if 'linkage' in kwargs:
            db_linkage = kwargs['linkage']
        else:
            raise ValidationError('Database execution mode \'reap\' requires a linkage option.')
    else:
        raise ValidationError('Database execution mode \'%s\' not valid.' % (db_mode))

    #   Option counterpoise- whether for interaction energy databases run in bsse-corrected or not
    db_cp = 'no'
    if 'cp' in kwargs:
        db_cp = kwargs['cp']

    if input.yes.match(str(db_cp)):
        try:
            database.ACTV_CP
        except AttributeError:
            raise ValidationError('Counterpoise correction mode \'yes\' invalid for database %s.' % (db_name))
        else:
            ACTV = database.ACTV_CP
    elif input.no.match(str(db_cp)):
        pass
    else:
        raise ValidationError('Counterpoise correction mode \'%s\' not valid.' % (db_cp))

    #   Option relaxed- whether for non-frozen-monomer interaction energy databases include deformation correction or not?
    db_rlxd = 'no'
    if 'rlxd' in kwargs:
        db_rlxd = kwargs['rlxd']

    if input.yes.match(str(db_rlxd)):
        if input.yes.match(str(db_cp)):
            try:
                database.ACTV_CPRLX
                database.RXNM_CPRLX
            except AttributeError:
                raise ValidationError('Deformation and counterpoise correction mode \'yes\' invalid for database %s.' % (db_name))
            else:
                ACTV = database.ACTV_CPRLX
                RXNM = database.RXNM_CPRLX
        elif input.no.match(str(db_cp)):
            try:
                database.ACTV_RLX
            except AttributeError:
                raise ValidationError('Deformation correction mode \'yes\' invalid for database %s.' % (db_name))
            else:
                ACTV = database.ACTV_RLX
    elif input.no.match(str(db_rlxd)):
        pass
    else:
        raise ValidationError('Deformation correction mode \'%s\' not valid.' % (db_rlxd))

    #   Option zero-point-correction- whether for thermochem databases jobs are corrected by zpe
    db_zpe = 'no'
    if 'zpe' in kwargs:
        db_zpe = kwargs['zpe']

    if input.yes.match(str(db_zpe)):
        raise ValidationError('Zero-point-correction mode \'yes\' not yet implemented.')
    elif input.no.match(str(db_zpe)):
        pass
    else:
        raise ValidationError('Zero-point-correction \'mode\' %s not valid.' % (db_zpe))

    #   Option benchmark- whether error statistics computed wrt alternate reference energies
    db_benchmark = 'default'
    if 'benchmark' in kwargs:
        db_benchmark = kwargs['benchmark']

        if (db_benchmark.lower() == 'default'):
            pass
        else:
            try:
                getattr(database, 'BIND_' + db_benchmark)
            except AttributeError:
                raise ValidationError('Special benchmark \'%s\' not available for database %s.' % (db_benchmark, db_name))
            else:
                BIND = getattr(database, 'BIND_' + db_benchmark)

    #   Option tabulate- whether tables of variables other than primary energy method are formed
    db_tabulate = []
    if 'tabulate' in kwargs:
        db_tabulate = kwargs['tabulate']

    #   Option subset- whether all of the database or just a portion is run
    db_subset = HRXN
    if 'subset' in kwargs:
        db_subset = kwargs['subset']

    if isinstance(db_subset, basestring):
        if (db_subset.lower() == 'small'):
            try:
                database.HRXN_SM
            except AttributeError:
                raise ValidationError('Special subset \'small\' not available for database %s.' % (db_name))
            else:
                HRXN = database.HRXN_SM
        elif (db_subset.lower() == 'large'):
            try:
                database.HRXN_LG
            except AttributeError:
                raise ValidationError('Special subset \'large\' not available for database %s.' % (db_name))
            else:
                HRXN = database.HRXN_LG
        elif (db_subset.lower() == 'equilibrium'):
            try:
                database.HRXN_EQ
            except AttributeError:
                raise ValidationError('Special subset \'equilibrium\' not available for database %s.' % (db_name))
            else:
                HRXN = database.HRXN_EQ
        else:
            try:
                getattr(database, db_subset)
            except AttributeError:
                raise ValidationError('Special subset \'%s\' not available for database %s.' % (db_subset, db_name))
            else:
                HRXN = getattr(database, db_subset)
    else:
        temp = []
        for rxn in db_subset:
            if rxn in HRXN:
                temp.append(rxn)
            else:
                raise ValidationError('Subset element \'%s\' not a member of database %s.' % (str(rxn), db_name))
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
    if (db_mode.lower() == 'continuous'):
        instructions = """\n    The database single-job procedure has been selected through mode='continuous'.\n"""
        instructions += """    Calculations for the reagents will proceed in the order below and will be followed\n"""
        instructions += """    by summary results for the database.\n\n"""
        for rgt in HSYS:
            instructions += """                    %-s\n""" % (rgt)
        instructions += """\n    Alternatively, a farming-out of the database calculations may be accessed through\n"""
        instructions += """    the database wrapper option mode='sow'/'reap'.\n\n"""
        PsiMod.print_out(instructions)

    #   write sow/reap instructions and index of calcs to output file and reap input file
    if (db_mode.lower() == 'sow'):
        instructions = """\n    The database sow/reap procedure has been selected through mode='sow'. In addition\n"""
        instructions += """    to this output file (which contains no quantum chemical calculations), this job\n"""
        instructions += """    has produced a number of input files (%s-*.in) for individual database members\n""" % (dbse)
        instructions += """    and a single input file (%s-master.in) with a database(mode='reap') command.\n""" % (dbse)
        instructions += """    The former may look very peculiar since processed and pickled python rather than\n"""
        instructions += """    raw input is written. Follow the instructions below to continue.\n\n"""
        instructions += """    (1)  Run all of the %s-*.in input files on any variety of computer architecture.\n""" % (dbse)
        instructions += """       The output file names must be as given below.\n\n"""
        for rgt in HSYS:
            instructions += """             psi4 -i %-27s -o %-27s\n""" % (rgt + '.in', rgt + '.out')
        instructions += """\n    (2)  Gather all the resulting output files in a directory. Place input file\n"""
        instructions += """         %s-master.in into that directory and run it. The job will be trivial in\n""" % (dbse)
        instructions += """         length and give summary results for the database in its output file.\n\n"""
        instructions += """             psi4 -i %-27s -o %-27s\n\n""" % (dbse + '-master.in', dbse + '-master.out')
        instructions += """    Alternatively, a single-job execution of the database may be accessed through\n"""
        instructions += """    the database wrapper option mode='continuous'.\n\n"""
        PsiMod.print_out(instructions)

        fmaster = open('%s-master.in' % (dbse), 'w')
        fmaster.write('# This is a psi4 input file auto-generated from the database() wrapper.\n\n')
        fmaster.write("database('%s', '%s', mode='reap', cp='%s', rlxd='%s', zpe='%s', benchmark='%s', linkage=%d, subset=%s, tabulate=%s)\n\n" %
            (name, db_name, db_cp, db_rlxd, db_zpe, db_benchmark, os.getpid(), HRXN, db_tabulate))
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
        commands = ''
        commands += """\nPsiMod.set_memory(%s)\n\n""" % (user_memory)
        for chgdopt in PsiMod.get_global_option_list():
            if PsiMod.has_option_changed(chgdopt):
                chgdoptval = PsiMod.get_global_option(chgdopt)
                #chgdoptval = PsiMod.get_option(chgdopt)
                if isinstance(chgdoptval, basestring):
                    commands += """PsiMod.set_global_option('%s', '%s')\n""" % (chgdopt, chgdoptval)
                elif isinstance(chgdoptval, int) or isinstance(chgdoptval, float):
                    commands += """PsiMod.set_global_option('%s', %s)\n""" % (chgdopt, chgdoptval)
                else:
                    raise ValidationError('Option \'%s\' is not of a type (string, int, float, bool) that can be processed by database wrapper.' % (chgdopt))

        # build string of molecule and commands that are dependent on the database
        commands += '\n'
        commands += """PsiMod.set_global_option('BASIS', '%s')\n""" % (user_basis)
        if not((user_df_basis_scf == "") or (user_df_basis_scf == 'NONE')):
            commands += """PsiMod.set_global_option('DF_BASIS_SCF', '%s')\n""" % (user_df_basis_scf)
        if not((user_df_basis_mp2 == "") or (user_df_basis_mp2 == 'NONE')):
            commands += """PsiMod.set_global_option('DF_BASIS_MP2', '%s')\n""" % (user_df_basis_mp2)
        if not((user_df_basis_cc == "") or (user_df_basis_cc == 'NONE')):
            commands += """PsiMod.set_global_option('DF_BASIS_CC', '%s')\n""" % (user_df_basis_cc)
        if not((user_df_basis_sapt == "") or (user_df_basis_sapt == 'NONE')):
            commands += """PsiMod.set_global_option('DF_BASIS_SAPT', '%s')\n""" % (user_df_basis_sapt)
        if not((user_df_basis_elst == "") or (user_df_basis_elst == 'NONE')):
            commands += """PsiMod.set_global_option('DF_BASIS_ELST', '%s')\n""" % (user_df_basis_elst)
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
        if (db_mode.lower() == 'continuous'):
            exec banners
            exec GEOS[rgt]
            exec commands
            #print 'MOLECULE LIVES %23s %8s %4d %4d %4s' % (rgt, PsiMod.get_option('REFERENCE'),
            #    molecule.molecular_charge(), molecule.multiplicity(), molecule.schoenflies_symbol())
            PsiMod.set_variable('NATOM', molecule.natom())
            ERGT[rgt] = call_function_in_1st_argument(func, **kwargs)
            #print ERGT[rgt]
            PsiMod.print_variables()
            exec actives
            for envv in db_tabulate:
                VRGT[rgt][envv] = PsiMod.get_variable(envv)
            PsiMod.set_global_option("REFERENCE", user_reference)
            PsiMod.clean()

        elif (db_mode.lower() == 'sow'):
            freagent = open('%s.in' % (rgt), 'w')
            freagent.write('# This is a psi4 input file auto-generated from the database() wrapper.\n\n')
            freagent.write(banners)
            freagent.write(GEOS[rgt])
            freagent.write(commands)
            freagent.write('''\npickle_kw = ("""''')
            pickle.dump(kwargs, freagent)
            freagent.write('''""")\n''')
            freagent.write("""\nkwargs = pickle.loads(pickle_kw)\n""")
            freagent.write("""electronic_energy = %s(**kwargs)\n\n""" % (func.__name__))
            freagent.write("""PsiMod.print_variables()\n""")
            freagent.write("""PsiMod.print_out('\\nDATABASE RESULT: computation %d for reagent %s """
                % (os.getpid(), rgt))
            freagent.write("""yields electronic energy %20.12f\\n' % (electronic_energy))\n\n""")
            freagent.write("""PsiMod.set_variable('NATOM', molecule.natom())\n""")
            for envv in db_tabulate:
                freagent.write("""PsiMod.print_out('DATABASE RESULT: computation %d for reagent %s """
                    % (os.getpid(), rgt))
                freagent.write("""yields variable value    %20.12f for variable %s\\n' % (PsiMod.get_variable(""")
                freagent.write("""'%s'), '%s'))\n""" % (envv.upper(), envv.upper()))
            freagent.close()

        elif (db_mode.lower() == 'reap'):
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
                            raise ValidationError('Output file \'%s.out\' has linkage %s incompatible with master.in linkage %s.'
                                % (rgt, str(s[3]), str(db_linkage)))
                        if s[6] != rgt:
                            raise ValidationError('Output file \'%s.out\' has nominal affiliation %s incompatible with reagent %s.'
                                % (rgt, s[6], rgt))
                        if (s[8:10] == ['electronic', 'energy']):
                            ERGT[rgt] = float(s[10])
                            PsiMod.print_out('DATABASE RESULT: electronic energy = %20.12f\n' % (ERGT[rgt]))
                        elif (s[8:10] == ['variable', 'value']):
                            for envv in db_tabulate:
                                if (s[13:] == envv.upper().split()):
                                    VRGT[rgt][envv] = float(s[10])
                                    PsiMod.print_out('DATABASE RESULT: variable %s value    = %20.12f\n' % (envv.upper(), VRGT[rgt][envv]))
                freagent.close()

    #   end sow after writing files
    if (db_mode.lower() == 'sow'):
        return 0.0

    # Reap all the necessary reaction computations
    PsiMod.print_out("\n")
    banner(("Database %s Results" % (db_name)))
    PsiMod.print_out("\n")

    maxactv = []
    for rxn in HRXN:
        maxactv.append(len(ACTV[dbse + '-' + str(rxn)]))
    maxrgt = max(maxactv)
    table_delimit = '-' * (54 + 20 * maxrgt)
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
        tables += """\n   ==> %s <==\n\n""" % (envv.title())
        tables += tblhead(maxrgt, table_delimit, 2)

        for rxn in HRXN:
            db_rxn = dbse + '-' + str(rxn)
            VRXN[db_rxn] = {}

            if FAIL[rxn]:
                tables += """\n%23s   %8s %8s   %8s""" % (db_rxn, '', '****', '')
                for i in range(len(ACTV[db_rxn])):
                    tables += """ %16.8f %2.0f""" % (VRGT[ACTV[db_rxn][i]][envv], RXNM[db_rxn][ACTV[db_rxn][i]])

            else:
                VRXN[db_rxn][envv] = 0.0
                for i in range(len(ACTV[db_rxn])):
                    VRXN[db_rxn][envv] += VRGT[ACTV[db_rxn][i]][envv] * RXNM[db_rxn][ACTV[db_rxn][i]]

                tables += """\n%23s        %16.8f       """ % (db_rxn, VRXN[db_rxn][envv])
                for i in range(len(ACTV[db_rxn])):
                    tables += """ %16.8f %2.0f""" % (VRGT[ACTV[db_rxn][i]][envv], RXNM[db_rxn][ACTV[db_rxn][i]])
        tables += """\n   %s\n""" % (table_delimit)

    #   tabulate primary requested energy variable with statistics
    count_rxn = 0
    minDerror = 100000.0
    maxDerror = 0.0
    MSDerror = 0.0
    MADerror = 0.0
    RMSDerror = 0.0

    tables += """\n   ==> %s <==\n\n""" % ('Requested Energy')
    tables += tblhead(maxrgt, table_delimit, 1)
    for rxn in HRXN:
        db_rxn = dbse + '-' + str(rxn)

        if FAIL[rxn]:
            tables += """\n%23s   %8.4f %8s   %8s""" % (db_rxn, BIND[db_rxn], '****', '****')
            for i in range(len(ACTV[db_rxn])):
                tables += """ %16.8f %2.0f""" % (ERGT[ACTV[db_rxn][i]], RXNM[db_rxn][ACTV[db_rxn][i]])

        else:
            ERXN[db_rxn] = 0.0
            for i in range(len(ACTV[db_rxn])):
                ERXN[db_rxn] += ERGT[ACTV[db_rxn][i]] * RXNM[db_rxn][ACTV[db_rxn][i]]
            error = physconst.psi_hartree2kcalmol * ERXN[db_rxn] - BIND[db_rxn]

            tables += """\n%23s   %8.4f %8.4f   %8.4f""" % (db_rxn, BIND[db_rxn], physconst.psi_hartree2kcalmol * ERXN[db_rxn], error)
            for i in range(len(ACTV[db_rxn])):
                tables += """ %16.8f %2.0f""" % (ERGT[ACTV[db_rxn][i]], RXNM[db_rxn][ACTV[db_rxn][i]])

            if abs(error) < abs(minDerror):
                minDerror = error
            if abs(error) > abs(maxDerror):
                maxDerror = error
            MSDerror += error
            MADerror += abs(error)
            RMSDerror += error * error
            count_rxn += 1
    tables += """\n   %s\n""" % (table_delimit)

    if count_rxn:

        MSDerror /= float(count_rxn)
        MADerror /= float(count_rxn)
        RMSDerror = sqrt(RMSDerror / float(count_rxn))

        tables += """%23s   %19s %8.4f\n""" % ('Minimal Dev', '', minDerror)
        tables += """%23s   %19s %8.4f\n""" % ('Maximal Dev', '', maxDerror)
        tables += """%23s   %19s %8.4f\n""" % ('Mean Signed Dev', '', MSDerror)
        tables += """%23s   %19s %8.4f\n""" % ('Mean Absolute Dev', '', MADerror)
        tables += """%23s   %19s %8.4f\n""" % ('RMS Dev', '', RMSDerror)
        tables += """   %s\n""" % (table_delimit)

        PsiMod.set_variable('%s DATABASE MEAN SIGNED DEVIATION' % (db_name), MSDerror)
        PsiMod.set_variable('%s DATABASE MEAN ABSOLUTE DEVIATION' % (db_name), MADerror)
        PsiMod.set_variable('%s DATABASE ROOT-MEAN-SQUARE DEVIATION' % (db_name), RMSDerror)

        #print tables
        PsiMod.print_out(tables)
        finalenergy = MADerror

    else:
        finalenergy = 0.0

    # restore molecule and options
    activate(user_molecule)
    user_molecule.update_geometry()
    PsiMod.set_global_option("BASIS", user_basis)
    PsiMod.set_global_option("REFERENCE", user_reference)
    if not b_user_reference:
        PsiMod.revoke_global_option_changed('REFERENCE')

    return finalenergy


def drop_duplicates(seq):
    """Function that given an array *seq*, returns an array without any duplicate
    entries. There is no guarantee of which duplicate entry is dropped.

    """
    noDupes = []
    [noDupes.append(i) for i in seq if not noDupes.count(i)]
    return noDupes


def tblhead(tbl_maxrgt, tbl_delimit, ttype):
    """Function that prints the header for the changable-width results tables in db().
    *tbl_maxrgt* is the number of reagent columns the table must plan for. *tbl_delimit*
    is a string of dashes of the correct length to set off the table. *ttype* is 1 for
    tables comparing the computed values to the reference or 2 for simple tabulation
    and sum of the computed values.

    """
    tbl_str = ''
    tbl_str += """   %s""" % (tbl_delimit)
    if ttype == 1:
        tbl_str += """\n%23s %19s   %8s""" % ('Reaction', 'Reaction Energy', 'Error')
    elif ttype == 2:
        tbl_str += """\n%23s     %19s %6s""" % ('Reaction', 'Reaction Value', '')
    for i in range(tbl_maxrgt):
        tbl_str += """%20s""" % ('Reagent ' + str(i + 1))
    if ttype == 1:
        tbl_str += """\n%23s   %8s %8s %8s""" % ('', 'Ref', 'Calc', '[kcal/mol]')
    elif ttype == 2:
        tbl_str += """\n%54s""" % ('')
    for i in range(tbl_maxrgt):
        if ttype == 1:
            tbl_str += """%20s""" % ('[H] Wt')
        elif ttype == 2:
            tbl_str += """%20s""" % ('Value Wt')
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
    """Function to define a multistage energy method from combinations of
    basis set extrapolations and delta corrections and condense the
    components into a minimum number of calculations.

    :aliases: cbs()

    :returns: (*float*) -- Total electronic energy in Hartrees

    :PSI variables:

    .. envvar:: CBS TOTAL ENERGY
        CBS REFERENCE ENERGY
        CBS CORRELATION ENERGY
        CURRENT ENERGY
        CURRENT REFERENCE ENERGY
        CURRENT CORRELATION ENERGY

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Methods beyond basic scf, mp2, ccsd, ccsd(t) not yet hooked in through PSI variables, df-mp2 in particular.

       - No scheme defaults for given basis zeta number, so scheme must be specified explicitly.

       - No way to tell function to boost fitting basis size for all calculations.

       - No way to extrapolate def2 family basis sets

       - Need to add more extrapolation schemes

    As represented in the equation below, a CBS energy method is defined in four
    sequential stages (scf, corl, delta, delta2) covering treatment of the
    reference total energy, the correlation energy, a delta correction to the
    correlation energy, and a second delta correction. Each is activated by its
    stage_wfn keyword and is only allowed if all preceding stages are active.

    .. include:: cbs_eqn.rst

    **Keywords**

    * Energy Methods
        The presence of a stage_wfn keyword is the indicator to incorporate
        (and check for stage_basis and stage_scheme keywords) and compute
        that stage in defining the CBS energy.

    The cbs() function requires, at a minimum, ``name='scf'`` and ``scf_basis``
        keywords to be specified for reference-step only jobs and ``name`` and
        ``corl_basis`` keywords for correlated jobs.

    :type name: string
    :param name: ``'scf'`` || ``'ccsd'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        for the correlation energy, unless only reference step to be performed,
        in which case should be ``'scf'``. Overruled if stage_wfn keywords supplied.

    :type corl_wfn: string
    :param corl_wfn: ``'mp2'`` || ``'ccsd(t)'`` || etc.

        Indicates the energy method for which the correlation energy is to be
        obtained. Can also be specified with ``name`` or as the unlabeled
        first argument to the function.

    :type delta_wfn: string
    :param delta_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a delta correction
        to the correlation energy is to be obtained.

    :type delta_wfn_lesser: string
    :param delta_wfn_lesser: |dl| ``'mp2'`` |dr| || ``'ccsd'`` || etc.

        Indicates the inferior energy method for which a delta correction
        to the correlation energy is to be obtained.

    :type delta2_wfn: string
    :param delta2_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a second delta correction
        to the correlation energy is to be obtained.

    :type delta2_wfn_lesser: string
    :param delta2_wfn_lesser: |dl| ``'mp2'`` |dr| || ``'ccsd(t)'`` || etc.

        Indicates the inferior energy method for which a second delta correction
        to the correlation energy is to be obtained.

    * Basis Sets
        Currently, the basis set set through ``set`` commands have no influence
        on a cbs calculation.

    :type scf_basis: string
    :param scf_basis: |dl| ``corl_basis`` |dr| || ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the reference energy.
        If any correlation method is specified, ``scf_basis`` can default
        to ``corl_basis``.

    :type corl_basis: string
    :param corl_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the correlation energy.

    :type delta_basis: string
    :param delta_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the delta correction
        to the correlation energy.

    :type delta2_basis: string
    :param delta2_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the second delta correction
        to the correlation energy.

    * Schemes
        Transformations of the energy through basis set extrapolation for each
        stage of the CBS definition. A complaint is generated if number of basis
        sets in stage_basis does not exactly satisfy requirements of stage_scheme.
        An exception is the default, ``'highest_1'``, which uses the best basis
        set available. See `Extrapolation Schemes`_ for all available schemes.

    :type scf_scheme: function
    :param scf_scheme: |dl| ``highest_1`` |dr| || ``scf_xtpl_helgaker_3`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the reference energy.

    :type corl_scheme: function
    :param corl_scheme: |dl| ``highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the correlation energy.

    :type delta_scheme: function
    :param delta_scheme: |dl| ``highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the delta correction
        to the correlation energy.

    :type delta2_scheme: function
    :param delta2_scheme: |dl| ``highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the second delta correction
        to the correlation energy.

    **Examples**

    >>> # [1] replicates with cbs() the simple model chemistry scf/cc-pVDZ: set basis cc-pVDZ energy('scf')
    >>> cbs('scf', scf_basis='cc-pVDZ')

    >>> # [2] replicates with cbs() the simple model chemistry mp2/jun-cc-pVDZ: set basis jun-cc-pVDZ energy('mp2')
    >>> cbs('mp2', corl_basis='jun-cc-pVDZ')

    >>> # [3] DTQ-zeta extrapolated scf reference energy
    >>> cbs('scf', scf_basis='cc-pV[DTQ]Z', scf_scheme=scf_xtpl_helgaker_3)

    >>> # [4] DT-zeta extrapolated mp2 correlation energy atop a T-zeta reference
    >>> cbs('mp2', corl_basis='cc-pv[dt]z', corl_scheme=corl_xtpl_helgaker_2)

    >>> # [5] a DT-zeta extrapolated coupled-cluster correction atop a TQ-zeta extrapolated mp2 correlation energy atop a Q-zeta reference
    >>> cbs('mp2', corl_basis='aug-cc-pv[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z', delta_scheme=corl_xtpl_helgaker_2)

    >>> # [6] a D-zeta ccsd(t) correction atop a DT-zeta extrapolated ccsd cluster correction atop a TQ-zeta extrapolated mp2 correlation energy atop a Q-zeta reference
    >>> cbs('mp2', corl_basis='aug-cc-pv[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd', delta_basis='aug-cc-pv[dt]z', delta_scheme=corl_xtpl_helgaker_2, delta2_wfn='ccsd(t)', delta2_wfn_lesser='ccsd', delta2_basis='aug-cc-pvdz')

    >>> # [7] cbs() coupled with database()
    >>> database('mp2', 'BASIC', subset=['h2o','nh3'], symm='on', func=cbs, corl_basis='cc-pV[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd(t)', delta_basis='sto-3g')

    """

    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # Wrap any positional arguments into kwargs (for intercalls among wrappers)
    if not('name' in kwargs) and name:
        kwargs['name'] = name.lower()

    # Establish function to call (only energy makes sense for cbs)
    if not('cbs_func' in kwargs):
        if ('func' in kwargs):
            kwargs['cbs_func'] = kwargs['func']
            del kwargs['func']
        else:
            kwargs['cbs_func'] = energy
    func = kwargs['cbs_func']
    if not func:
        raise ValidationError('Function \'%s\' does not exist to be called by wrapper complete_basis_set.' % (func.__name__))
    if not(func is energy):
        raise ValidationError('Wrapper complete_basis_set is unhappy to be calling function \'%s\' instead of \'energy\'.' % (func.__name__))

    # Define some quantum chemical knowledge, namely what methods are subsumed in others
    VARH = {}
    VARH['scf'] = {         'scftot': 'SCF TOTAL ENERGY'}
    VARH['mp2'] = {         'scftot': 'SCF TOTAL ENERGY',
                           'mp2corl': 'MP2 CORRELATION ENERGY'}
    VARH['ccsd'] = {        'scftot': 'SCF TOTAL ENERGY',
                           'mp2corl': 'MP2 CORRELATION ENERGY',
                          'ccsdcorl': 'CCSD CORRELATION ENERGY'}
    VARH['ccsd(t)'] = {     'scftot': 'SCF TOTAL ENERGY',
                           'mp2corl': 'MP2 CORRELATION ENERGY',
                          'ccsdcorl': 'CCSD CORRELATION ENERGY',
                       'ccsd(t)corl': 'CCSD(T) CORRELATION ENERGY'}

    finalenergy = 0.0
    do_scf = 1
    do_corl = 0
    do_delta = 0
    do_delta2 = 0

    # Must collect (here) and set (below) basis sets after every new molecule activation
    b_user_basis = PsiMod.has_global_option_changed('BASIS')
    user_basis = PsiMod.get_option('BASIS')
    #user_df_basis_scf = PsiMod.get_option('DF_BASIS_SCF')
    #user_df_basis_mp2 = PsiMod.get_option('DF_BASIS_MP2')
    #user_df_basis_cc = PsiMod.get_option('DF_BASIS_CC')
    #user_df_basis_sapt = PsiMod.get_option('DF_BASIS_SAPT')
    #user_df_basis_elst = PsiMod.get_option('DF_BASIS_ELST')
    b_user_wfn = PsiMod.has_global_option_changed('WFN')
    user_wfn = PsiMod.get_option('WFN')

    # Make sure the molecule the user provided is the active one
    if 'molecule' in kwargs:
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    PsiMod.set_global_option("BASIS", PsiMod.get_global_option("BASIS"))

    # Establish method for correlation energy
    if 'name' in kwargs:
        if (lowername == 'scf'):
            pass
        else:
            do_corl = 1
            cbs_corl_wfn = kwargs['name'].lower()
    if 'corl_wfn' in kwargs:
        do_corl = 1
        cbs_corl_wfn = kwargs['corl_wfn'].lower()
    if do_corl:
        if not (cbs_corl_wfn in VARH.keys()):
            raise ValidationError('Requested CORL method \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_corl_wfn))

    # Establish method for delta correction energy
    if 'delta_wfn' in kwargs:
        do_delta = 1
        cbs_delta_wfn = kwargs['delta_wfn'].lower()
        if not (cbs_delta_wfn in VARH.keys()):
            raise ValidationError('Requested DELTA method \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta_wfn))

        if 'delta_wfn_lesser' in kwargs:
            cbs_delta_wfn_lesser = kwargs['delta_wfn_lesser'].lower()
        else:
            cbs_delta_wfn_lesser = 'mp2'
        if not (cbs_delta_wfn_lesser in VARH.keys()):
            raise ValidationError('Requested DELTA method lesser \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta_wfn_lesser))

    # Establish method for second delta correction energy
    if 'delta2_wfn' in kwargs:
        do_delta2 = 1
        cbs_delta2_wfn = kwargs['delta2_wfn'].lower()
        if not (cbs_delta2_wfn in VARH.keys()):
            raise ValidationError('Requested DELTA2 method \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta2_wfn))

        if 'delta2_wfn_lesser' in kwargs:
            cbs_delta2_wfn_lesser = kwargs['delta2_wfn_lesser'].lower()
        else:
            cbs_delta2_wfn_lesser = 'mp2'
        if not (cbs_delta2_wfn_lesser in VARH.keys()):
            raise ValidationError('Requested DELTA2 method lesser \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta2_wfn_lesser))

    # Check that user isn't skipping steps in scf + corl + delta + delta2 sequence
    if do_scf and not do_corl and not do_delta and not do_delta2:
        pass
    elif do_scf and do_corl and not do_delta and not do_delta2:
        pass
    elif do_scf and do_corl and do_delta and not do_delta2:
        pass
    elif do_scf and do_corl and do_delta and do_delta2:
        pass
    else:
        raise ValidationError('Requested scf (%s) + corl (%s) + delta (%s) + delta2 (%s) not valid. These steps are cummulative.' %
            (do_scf, do_corl, do_delta, do_delta2))

    # Establish list of valid basis sets for correlation energy
    if do_corl:
        if 'corl_basis' in kwargs:
            BSTC, ZETC = validate_bracketed_basis(kwargs['corl_basis'].lower())
        else:
            raise ValidationError('CORL basis sets through keyword \'%s\' are required.' % ('corl_basis'))

    # Establish list of valid basis sets for scf energy
    if 'scf_basis' in kwargs:
        BSTR, ZETR = validate_bracketed_basis(kwargs['scf_basis'].lower())
    else:
        if do_corl:
            BSTR = BSTC[:]
            ZETR = ZETC[:]
        else:
            raise ValidationError('SCF basis sets through keyword \'%s\' are required. Or perhaps you forgot the \'%s\'.' % ('scf_basis', 'corl_wfn'))

    # Establish list of valid basis sets for delta correction energy
    if do_delta:
        if 'delta_basis' in kwargs:
            BSTD, ZETD = validate_bracketed_basis(kwargs['delta_basis'].lower())
        else:
            raise ValidationError('DELTA basis sets through keyword \'%s\' are required.' % ('delta_basis'))

    # Establish list of valid basis sets for second delta correction energy
    if do_delta2:
        if 'delta2_basis' in kwargs:
            BSTD2, ZETD2 = validate_bracketed_basis(kwargs['delta2_basis'].lower())
        else:
            raise ValidationError('DELTA2 basis sets through keyword \'%s\' are required.' % ('delta2_basis'))

    # Establish treatment for scf energy (validity check useless since python will catch it long before here)
    cbs_scf_scheme = highest_1
    if 'scf_scheme' in kwargs:
        cbs_scf_scheme = kwargs['scf_scheme']

    # Establish treatment for correlation energy
    cbs_corl_scheme = highest_1
    if 'corl_scheme' in kwargs:
        cbs_corl_scheme = kwargs['corl_scheme']

    # Establish treatment for delta correction energy
    cbs_delta_scheme = highest_1
    if 'delta_scheme' in kwargs:
        cbs_delta_scheme = kwargs['delta_scheme']

    # Establish treatment for delta2 correction energy
    cbs_delta2_scheme = highest_1
    if 'delta2_scheme' in kwargs:
        cbs_delta2_scheme = kwargs['delta2_scheme']

    # Build string of title banner
    cbsbanners = ''
    cbsbanners += """PsiMod.print_out('\\n')\n"""
    cbsbanners += """banner(' CBS Setup ')\n"""
    cbsbanners += """PsiMod.print_out('\\n')\n\n"""
    exec cbsbanners

    # Call schemes for each portion of total energy to 'place orders' for calculations needed
    d_fields = ['d_stage', 'd_scheme', 'd_basis', 'd_wfn', 'd_need', 'd_coef', 'd_energy']
    f_fields = ['f_wfn', 'f_portion', 'f_basis', 'f_zeta', 'f_energy']
    GRAND_NEED = []
    MODELCHEM = []
    bstring = ''
    if do_scf:
        NEED = call_function_in_1st_argument(cbs_scf_scheme,
            mode='requisition', basisname=BSTR, basiszeta=ZETR, wfnname='scf')
        GRAND_NEED.append(dict(zip(d_fields, ['scf', cbs_scf_scheme, reconstitute_bracketed_basis(NEED), 'scf', NEED, +1, 0.0])))

    if do_corl:
        NEED = call_function_in_1st_argument(cbs_corl_scheme,
            mode='requisition', basisname=BSTC, basiszeta=ZETC, wfnname=cbs_corl_wfn)
        GRAND_NEED.append(dict(zip(d_fields, ['corl', cbs_corl_scheme, reconstitute_bracketed_basis(NEED), cbs_corl_wfn, NEED, +1, 0.0])))

    if do_delta:
        NEED = call_function_in_1st_argument(cbs_delta_scheme,
            mode='requisition', basisname=BSTD, basiszeta=ZETD, wfnname=cbs_delta_wfn)
        GRAND_NEED.append(dict(zip(d_fields, ['delta', cbs_delta_scheme, reconstitute_bracketed_basis(NEED), cbs_delta_wfn, NEED, +1, 0.0])))

        NEED = call_function_in_1st_argument(cbs_delta_scheme,
            mode='requisition', basisname=BSTD, basiszeta=ZETD, wfnname=cbs_delta_wfn_lesser)
        GRAND_NEED.append(dict(zip(d_fields, ['delta', cbs_delta_scheme, reconstitute_bracketed_basis(NEED), cbs_delta_wfn_lesser, NEED, -1, 0.0])))

    if do_delta2:
        NEED = call_function_in_1st_argument(cbs_delta2_scheme,
            mode='requisition', basisname=BSTD2, basiszeta=ZETD2, wfnname=cbs_delta2_wfn)
        GRAND_NEED.append(dict(zip(d_fields, ['delta2', cbs_delta2_scheme, reconstitute_bracketed_basis(NEED), cbs_delta2_wfn, NEED, +1, 0.0])))

        NEED = call_function_in_1st_argument(cbs_delta2_scheme,
            mode='requisition', basisname=BSTD2, basiszeta=ZETD2, wfnname=cbs_delta2_wfn_lesser)
        GRAND_NEED.append(dict(zip(d_fields, ['delta2', cbs_delta2_scheme, reconstitute_bracketed_basis(NEED), cbs_delta2_wfn_lesser, NEED, -1, 0.0])))

    for stage in GRAND_NEED:
        for lvl in stage['d_need'].iteritems():
            MODELCHEM.append(lvl[1])

    # Apply chemical reasoning to choose the minimum computations to run
    JOBS = MODELCHEM[:]

    instructions = ''
    instructions += """    Naive listing of computations required.\n"""
    for mc in JOBS:
        instructions += """   %12s / %-24s for  %s\n""" % (mc['f_wfn'], mc['f_basis'], VARH[mc['f_wfn']][mc['f_wfn'] + mc['f_portion']])

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
                if (menial == job['f_wfn'] + job['f_portion']) and (mc['f_basis'] == job['f_basis']) and not (mc['f_wfn'] == job['f_wfn']):
                    del JOBS[indx_job]

    instructions += """\n    Enlightened listing of computations required.\n"""
    for mc in JOBS:
        instructions += """   %12s / %-24s for  %s\n""" % (mc['f_wfn'], mc['f_basis'], VARH[mc['f_wfn']][mc['f_wfn'] + mc['f_portion']])

    #     Expand listings to all that will be obtained
    JOBS_EXT = []
    for indx_job, job in enumerate(JOBS):
        for menial in VARH[job['f_wfn']]:
            temp_wfn, temp_portion = split_menial(menial)
            JOBS_EXT.append(dict(zip(f_fields, [temp_wfn, temp_portion, job['f_basis'], job['f_zeta'], 0.0])))

    #instructions += """\n    Full listing of computations to be obtained (required and bonus).\n"""
    #for mc in JOBS_EXT:
    #    instructions += """   %12s / %-24s for  %s\n""" % (mc['f_wfn'], mc['f_basis'], VARH[mc['f_wfn']][mc['f_wfn']+mc['f_portion']])
    PsiMod.print_out(instructions)

    # Run necessary computations
    for mc in JOBS:
        kwargs['name'] = mc['f_wfn']

        # Build string of title banner
        cbsbanners = ''
        cbsbanners += """PsiMod.print_out('\\n')\n"""
        cbsbanners += """banner(' CBS Computation: %s / %s ')\n""" % (mc['f_wfn'].upper(), mc['f_basis'].upper())
        cbsbanners += """PsiMod.print_out('\\n')\n\n"""
        exec cbsbanners

        # Build string of molecule and commands that are dependent on the database
        commands = '\n'
        commands += """\nPsiMod.set_global_option('BASIS', '%s')\n""" % (mc['f_basis'])
        exec commands

        # Make energy() call
        mc['f_energy'] = call_function_in_1st_argument(func, **kwargs)

        # Fill in energies for subsumed methods
        for menial in VARH[mc['f_wfn']]:
            temp_wfn, temp_portion = split_menial(menial)
            for job in JOBS_EXT:
                if (temp_wfn == job['f_wfn']) and (temp_portion == job['f_portion']) and (mc['f_basis'] == job['f_basis']):
                    job['f_energy'] = PsiMod.get_variable(VARH[temp_wfn][menial])

        PsiMod.clean()

    # Build string of title banner
    cbsbanners = ''
    cbsbanners += """PsiMod.print_out('\\n')\n"""
    cbsbanners += """banner(' CBS Results ')\n"""
    cbsbanners += """PsiMod.print_out('\\n')\n\n"""
    exec cbsbanners

    # Insert obtained energies into the array that stores the cbs stages
    for stage in GRAND_NEED:
        for lvl in stage['d_need'].iteritems():
            MODELCHEM.append(lvl[1])

            for job in JOBS_EXT:
                if ((lvl[1]['f_wfn'] == job['f_wfn']) and (lvl[1]['f_portion'] == job['f_portion']) and
                   (lvl[1]['f_basis'] == job['f_basis'])):
                    lvl[1]['f_energy'] = job['f_energy']

    for stage in GRAND_NEED:
        stage['d_energy'] = call_function_in_1st_argument(stage['d_scheme'], needname=stage['d_need'], mode='evaluate')
        finalenergy += stage['d_energy'] * stage['d_coef']

    # Build string of results table
    table_delimit = '  ' + '-' * 105 + '\n'
    tables = ''
    tables += """\n   ==> %s <==\n\n""" % ('Components')
    tables += table_delimit
    tables += """     %6s %20s %1s %-26s %3s %16s   %-s\n""" % ('', 'Method', '/', 'Basis', 'Rqd', 'Energy [H]', 'Variable')
    tables += table_delimit
    for job in JOBS_EXT:
        star = ''
        for mc in MODELCHEM:
            if (job['f_wfn'] == mc['f_wfn']) and (job['f_basis'] == mc['f_basis']):
                star = '*'
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % ('', job['f_wfn'],
                  '/', job['f_basis'], star, job['f_energy'], VARH[job['f_wfn']][job['f_wfn'] + job['f_portion']])
    tables += table_delimit

    tables += """\n   ==> %s <==\n\n""" % ('Stages')
    tables += table_delimit
    tables += """     %6s %20s %1s %-27s %2s %16s   %-s\n""" % ('Stage', 'Method', '/', 'Basis', 'Wt', 'Energy [H]', 'Scheme')
    tables += table_delimit
    for stage in GRAND_NEED:
        tables += """     %6s %20s %1s %-27s %2d %16.8f   %-s\n""" % (stage['d_stage'], stage['d_wfn'],
                  '/', stage['d_basis'], stage['d_coef'], stage['d_energy'], stage['d_scheme'].__name__)
    tables += table_delimit

    tables += """\n   ==> %s <==\n\n""" % ('CBS')
    tables += table_delimit
    tables += """     %6s %20s %1s %-27s %2s %16s   %-s\n""" % ('Stage', 'Method', '/', 'Basis', '', 'Energy [H]', 'Scheme')
    tables += table_delimit
    if do_scf:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[0]['d_stage'], GRAND_NEED[0]['d_wfn'],
                  '/', GRAND_NEED[0]['d_basis'], '', GRAND_NEED[0]['d_energy'], GRAND_NEED[0]['d_scheme'].__name__)
    if do_corl:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[1]['d_stage'], GRAND_NEED[1]['d_wfn'],
                  '/', GRAND_NEED[1]['d_basis'], '', GRAND_NEED[1]['d_energy'], GRAND_NEED[1]['d_scheme'].__name__)
    if do_delta:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[2]['d_stage'], GRAND_NEED[2]['d_wfn'] + ' - ' + GRAND_NEED[3]['d_wfn'],
                  '/', GRAND_NEED[2]['d_basis'], '', GRAND_NEED[2]['d_energy'] - GRAND_NEED[3]['d_energy'], GRAND_NEED[2]['d_scheme'].__name__)
    if do_delta2:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[4]['d_stage'], GRAND_NEED[4]['d_wfn'] + ' - ' + GRAND_NEED[5]['d_wfn'],
                  '/', GRAND_NEED[4]['d_basis'], '', GRAND_NEED[4]['d_energy'] - GRAND_NEED[5]['d_energy'], GRAND_NEED[4]['d_scheme'].__name__)
    tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % ('total', 'CBS', '', '', '', finalenergy, '')
    tables += table_delimit

    #print tables
    PsiMod.print_out(tables)

    # Restore molecule and options
    #PsiMod.set_local_option('SCF', "WFN", user_wfn)    # TODO refuses to set global option WFN - rejects SCF as option
    PsiMod.set_global_option('BASIS', user_basis)

    PsiMod.set_global_option('WFN', user_wfn)
    if not b_user_wfn:
        PsiMod.revoke_global_option_changed('WFN')

    PsiMod.set_variable('CBS REFERENCE ENERGY', GRAND_NEED[0]['d_energy'])
    PsiMod.set_variable('CBS CORRELATION ENERGY', finalenergy - GRAND_NEED[0]['d_energy'])
    PsiMod.set_variable('CBS TOTAL ENERGY', finalenergy)
    PsiMod.set_variable('CURRENT REFERENCE ENERGY', GRAND_NEED[0]['d_energy'])
    PsiMod.set_variable('CURRENT CORRELATION ENERGY', finalenergy - GRAND_NEED[0]['d_energy'])
    PsiMod.set_variable('CURRENT ENERGY', finalenergy)
    return finalenergy


# Transform and validate basis sets from 'cc-pV[Q5]Z' into [cc-pVQZ, cc-pV5Z] and [4, 5]
def validate_bracketed_basis(basisstring):
    """Function to transform and validate basis sets for cbs(). A basis set with no
    paired square brackets is passed through with zeta level 0 (e.g., '6-31+G(d,p)'
    is returned as [6-31+G(d,p)] and [0]). A basis set with square brackets is
    checked for sensible sequence and Dunning-ness and returned as separate basis
    sets (e.g., 'cc-pV[Q5]Z' is returned as [cc-pVQZ, cc-pV5Z] and [4, 5]). Note
    that this function has no communication with the basis set library to check
    if the basis actually exists. Used by :py:func:`wrappers.complete_basis_set`.

    """
    ZETA = ['d', 't', 'q', '5', '6']
    BSET = []
    ZSET = []
    if re.match(r'.*cc-.*\[.*\].*z$', basisstring, flags=re.IGNORECASE):
        basispattern = re.compile(r'^(.*)\[(.*)\](.*)$')
        basisname = basispattern.match(basisstring)
        for b in basisname.group(2):
            if b not in ZETA:
                raise ValidationError('Basis set \'%s\' has invalid zeta level \'%s\'.' % (basisstring, b))
            if len(ZSET) != 0:
                if (int(ZSET[len(ZSET) - 1]) - ZETA.index(b)) != 1:
                    raise ValidationError('Basis set \'%s\' has out-of-order zeta level \'%s\'.' % (basisstring, b))
            BSET.append(basisname.group(1) + b + basisname.group(3))
            if b == 'd':
                b = '2'
            if b == 't':
                b = '3'
            if b == 'q':
                b = '4'
            ZSET.append(int(b))
    elif re.match(r'.*\[.*\].*$', basisstring, flags=re.IGNORECASE):
        raise ValidationError('Basis set surrounding series indicator [] in \'%s\' is invalid.' % (basisstring))
    else:
        BSET.append(basisstring)
        ZSET.append(0)

    return [BSET, ZSET]


# Reform string basis set descriptor from basis set strings, 'cc-pv[q5]z' from [cc-pvqz, cc-pv5z]
def reconstitute_bracketed_basis(needarray):
    """Function to reform a bracketed basis set string from a sequential series
    of basis sets (e.g, form 'cc-pv[q5]z' from array [cc-pvqz, cc-pv5z]). The
    basis set array is extracted from the *f_basis* field of a *NEED* dictionary in
    :py:func:`wrappers.complete_basis_set`. Result is used to print a nicely
    formatted basis set string in the results table.

    """
    ZETA = {'d': 2, 't': 3, 'q': 4, '5': 5, '6': 6}
    ZSET = [''] * len(ZETA)
    BSET = []

    for lvl in needarray.iteritems():
        BSET.append(lvl[1]['f_basis'])

    if (len(BSET) == 1):
        basisstring = BSET[0]
    else:
        indx = 0
        while indx < len(BSET[0]):
            if (BSET[0][indx] != BSET[1][indx]):
                zetaindx = indx
            indx += 1
        for basis in BSET:
            ZSET[ZETA[basis[zetaindx]] - 2] = basis[zetaindx]

        pre = BSET[0][:zetaindx]
        post = BSET[0][zetaindx + 1:]
        basisstring = pre + '[' + ''.join(ZSET) + ']' + post

    return basisstring


def highest_1(**largs):
    """Scheme for total or correlation energies with a single basis or the highest
    zeta-level among an array of bases. Used by :py:func:`wrappers.complete_basis_set`.

    .. math:: E_{total}^X = E_{total}^X

    """

    energypiece = 0.0
    functionname = sys._getframe().f_code.co_name
    f_fields = ['f_wfn', 'f_portion', 'f_basis', 'f_zeta', 'f_energy']
    [mode, NEED, wfnname, BSET, ZSET] = validate_scheme_args(functionname, **largs)

    if (mode == 'requisition'):

        # Impose restrictions on zeta sequence
        if (len(ZSET) == 0):
            raise ValidationError('Call to \'%s\' not valid with \'%s\' basis sets.' % (functionname, len(ZSET)))

        # Return array that logs the requisite jobs
        if (wfnname == 'scf'):
            portion = 'tot'
        else:
            portion = 'corl'
        NEED = {'HI': dict(zip(f_fields, [wfnname, portion, BSET[len(ZSET) - 1], ZSET[len(ZSET) - 1], 0.0]))}

        return NEED

    elif (mode == 'evaluate'):

        # Extract required energies and zeta integers from array
        # Compute extrapolated energy
        energypiece = NEED['HI']['f_energy']

        # Output string with extrapolation parameters
        cbsscheme = ''
        cbsscheme += """\n   ==> %s <==\n\n""" % (functionname)
        if (NEED['HI']['f_wfn'] == 'scf'):
            cbsscheme += """   HI-zeta (%s) Total Energy:        %16.8f\n""" % (str(NEED['HI']['f_zeta']), energypiece)
        else:
            cbsscheme += """   HI-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(NEED['HI']['f_zeta']), energypiece)
        PsiMod.print_out(cbsscheme)

        return energypiece


# Solution equation in LaTeX:  $E_{corl}^{\infty} = \frac{E_{corl}^{X} X^3 - E_{corl}^{X-1} (X-1)^3}{X^3 - (X-1)^3}$
# Solution equation in LaTeX:  $\beta = \frac{E_{corl}^{X} - E_{corl}^{X-1}}{X^{-3} - (X-1)^{-3}}$
def corl_xtpl_helgaker_2(**largs):
    """Extrapolation scheme for correlation energies with two adjacent zeta-level bases.
    Used by :py:func:`wrappers.complete_basis_set`.

    .. math:: E_{corl}^X = E_{corl}^{\infty} + \\beta X^{-3}

    """

    energypiece = 0.0
    functionname = sys._getframe().f_code.co_name
    f_fields = ['f_wfn', 'f_portion', 'f_basis', 'f_zeta', 'f_energy']
    [mode, NEED, wfnname, BSET, ZSET] = validate_scheme_args(functionname, **largs)

    if (mode == 'requisition'):

        # Impose restrictions on zeta sequence
        if (len(ZSET) != 2):
            raise ValidationError('Call to \'%s\' not valid with \'%s\' basis sets.' % (functionname, len(ZSET)))

        # Return array that logs the requisite jobs
        NEED = {'HI': dict(zip(f_fields, [wfnname, 'corl', BSET[1], ZSET[1], 0.0])),
                'LO': dict(zip(f_fields, [wfnname, 'corl', BSET[0], ZSET[0], 0.0]))}

        return NEED

    elif (mode == 'evaluate'):

        # Extract required energies and zeta integers from array
        eHI = NEED['HI']['f_energy']
        zHI = NEED['HI']['f_zeta']
        eLO = NEED['LO']['f_energy']
        zLO = NEED['LO']['f_zeta']

        # Compute extrapolated energy
        energypiece = (eHI * zHI ** 3 - eLO * zLO ** 3) / (zHI ** 3 - zLO ** 3)
        beta = (eHI - eLO) / (zHI ** (-3) - zLO ** (-3))

        # Output string with extrapolation parameters
        cbsscheme = ''
        cbsscheme += """\n   ==> %s <==\n\n""" % (functionname)
        cbsscheme += """   LO-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zLO), eLO)
        cbsscheme += """   HI-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zHI), eHI)
        cbsscheme += """   Extrapolated Correlation Energy: %16.8f\n""" % (energypiece)
        cbsscheme += """   Beta (coefficient) Value:        %16.8f\n""" % (beta)
        PsiMod.print_out(cbsscheme)

        return energypiece


def scf_xtpl_helgaker_3(**largs):
    """Extrapolation scheme for reference energies with three adjacent zeta-level bases.
    Used by :py:func:`wrappers.complete_basis_set`.

    .. math:: E_{total}^X = E_{total}^{\infty} + \\beta e^{-\\alpha X}

    """

    energypiece = 0.0
    functionname = sys._getframe().f_code.co_name
    f_fields = ['f_wfn', 'f_portion', 'f_basis', 'f_zeta', 'f_energy']
    [mode, NEED, wfnname, BSET, ZSET] = validate_scheme_args(functionname, **largs)

    if (mode == 'requisition'):

        # Impose restrictions on zeta sequence
        if (len(ZSET) != 3):
            raise ValidationError('Call to \'%s\' not valid with \'%s\' basis sets.' % (functionname, len(ZSET)))

        # Return array that logs the requisite jobs
        NEED = {'HI': dict(zip(f_fields, [wfnname, 'tot', BSET[2], ZSET[2], 0.0])),
                'MD': dict(zip(f_fields, [wfnname, 'tot', BSET[1], ZSET[1], 0.0])),
                'LO': dict(zip(f_fields, [wfnname, 'tot', BSET[0], ZSET[0], 0.0]))}

        return NEED

    elif (mode == 'evaluate'):

        # Extract required energies and zeta integers from array
        eHI = NEED['HI']['f_energy']
        eMD = NEED['MD']['f_energy']
        eLO = NEED['LO']['f_energy']
        zHI = NEED['HI']['f_zeta']
        zMD = NEED['MD']['f_zeta']
        zLO = NEED['LO']['f_zeta']

        # Compute extrapolated energy
        ratio = (eHI - eMD) / (eMD - eLO)
        alpha = -1 * math.log(ratio)
        beta = (eHI - eMD) / (math.exp(-1 * alpha * zMD) * (ratio - 1))
        energypiece = eHI - beta * math.exp(-1 * alpha * zHI)

        # Output string with extrapolation parameters
        cbsscheme = ''
        cbsscheme += """\n   ==> %s <==\n\n""" % (functionname)
        cbsscheme += """   LO-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zLO), eLO)
        cbsscheme += """   MD-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zMD), eMD)
        cbsscheme += """   HI-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zHI), eHI)
        cbsscheme += """   Extrapolated Correlation Energy: %16.8f\n""" % (energypiece)
        cbsscheme += """   Alpha (exponent) Value:          %16.8f\n""" % (alpha)
        cbsscheme += """   Beta (coefficient) Value:        %16.8f\n""" % (beta)
        PsiMod.print_out(cbsscheme)

        return energypiece


def scf_xtpl_helgaker_2(**largs):
    """Extrapolation scheme for reference energies with two adjacent zeta-level bases.
    Used by :py:func:`wrappers.complete_basis_set`.

    .. math:: E_{total}^X = E_{total}^{\infty} + \\beta e^{-\\alpha X}, \\alpha = 1.63

    """

    energypiece = 0.0
    functionname = sys._getframe().f_code.co_name
    f_fields = ['f_wfn', 'f_portion', 'f_basis', 'f_zeta', 'f_energy']
    [mode, NEED, wfnname, BSET, ZSET] = validate_scheme_args(functionname, **largs)

    if (mode == 'requisition'):

        # Impose restrictions on zeta sequence
        if (len(ZSET) != 2):
            raise ValidationError('Call to \'%s\' not valid with \'%s\' basis sets.' % (functionname, len(ZSET)))

        # Return array that logs the requisite jobs
        NEED = {'HI': dict(zip(f_fields, [wfnname, 'tot', BSET[1], ZSET[1], 0.0])),
                'LO': dict(zip(f_fields, [wfnname, 'tot', BSET[0], ZSET[0], 0.0]))}

        return NEED

    elif (mode == 'evaluate'):

        # Extract required energies and zeta integers from array
        eHI = NEED['HI']['f_energy']
        eLO = NEED['LO']['f_energy']
        zHI = NEED['HI']['f_zeta']
        zLO = NEED['LO']['f_zeta']

        # LAB TODO add ability to pass alternate parameter values in

        # Return extrapolated energy
        alpha = 1.63
        beta = (eHI - eLO) / (math.exp(-1 * alpha * zLO) * (math.exp(-1 * alpha) - 1))
        energypiece = eHI - beta * math.exp(-1 * alpha * zHI)

        # Output string with extrapolation parameters
        cbsscheme = ''
        cbsscheme += """\n   ==> %s <==\n\n""" % (functionname)
        cbsscheme += """   LO-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zLO), eLO)
        cbsscheme += """   HI-zeta (%s) Correlation Energy:  %16.8f\n""" % (str(zHI), eHI)
        cbsscheme += """   Extrapolated Correlation Energy: %16.8f\n""" % (energypiece)
        cbsscheme += """   Alpha (exponent) Value:          %16.8f\n""" % (alpha)
        cbsscheme += """   Beta (coefficient) Value:        %16.8f\n""" % (beta)
        PsiMod.print_out(cbsscheme)

        return energypiece


def validate_scheme_args(functionname, **largs):
    """Function called by each extrapolation scheme in :py:func:`wrappers.complete_basis_set`.
    Checks that all the input arguments are present and suitable so that
    the scheme function can focus on defining the extrapolation.

    """
    mode = ''
    NEED = []
    wfnname = ''
    BSET = []
    ZSET = []

    # Mode where function fills out a form NEED with the computations needed to fulfill its call
    if (largs['mode'].lower() == 'requisition'):
        mode = largs['mode'].lower()

        if 'wfnname' in largs:
            wfnname = largs['wfnname']
        else:
            raise ValidationError('Call to \'%s\' has keyword \'wfnname\' missing.' % (functionname))

        if re.match(r'scf_.*$', functionname) and (wfnname != 'scf'):
            raise ValidationError('Call to \'%s\' is intended for scf portion of calculation.' % (functionname))
        if re.match(r'corl_.*$', functionname) and (wfnname == 'scf'):
            raise ValidationError('Call to \'%s\' is not intended for scf portion of calculation.' % (functionname))

        if 'basisname' in largs:
            BSET = largs['basisname']
        else:
            raise ValidationError('Call to \'%s\' has keyword \'basisname\' missing.' % (functionname))

        if 'basiszeta' in largs:
            ZSET = largs['basiszeta']
        else:
            raise ValidationError('Call to \'%s\' has keyword \'basiszeta\' missing.' % (functionname))

    # Mode where function reads the now-filled-in energies from that same form and performs the sp, xtpl, delta, etc.
    elif (largs['mode'].lower() == 'evaluate'):
        mode = largs['mode'].lower()

        if 'needname' in largs:
            NEED = largs['needname']
        else:
            raise ValidationError('Call to \'%s\' has keyword \'needname\' missing.' % (functionname))

    else:
        raise ValidationError('Call to \'%s\' has keyword \'mode\' missing or invalid.' % (functionname))

    return [mode, NEED, wfnname, BSET, ZSET]


def split_menial(menial):
    """Function used by :py:func:`wrappers.complete_basis_set` to separate
    *menial* 'scftot' into [scf, tot] and 'mp2corl' into [mp2, corl].

    """
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
