from __future__ import print_function
"""Module with functions that encode the sequence of PSI module
calls for each of the *name* values of the energy(), optimize(),
response(), and frequency() function.

"""

import PsiMod
import shutil
import os
import subprocess
import re
import input
import physconst
from molutil import *
from text import *
from procutil import *
from basislist import *
from functional import *
# never import driver, wrappers, or aliases into this file


def run_dcft(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density cumulant functional theory calculation.

    """
    oldref = PsiMod.get_global_option('REFERENCE')
    PsiMod.set_global_option('REFERENCE', 'UHF')
    PsiMod.scf()
    return PsiMod.dcft()
    PsiMod.set_global_option('REFERENCE', oldref)


def run_dcft_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    DCFT gradient calculation.

    """
    PsiMod.set_global_option('DERTYPE', 'FIRST')

    run_dcft(name, **kwargs)
    PsiMod.deriv()

def run_omp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an orbital-optimized MP2 computation

    """
    oldref = PsiMod.get_global_option('REFERENCE')
    PsiMod.set_global_option('REFERENCE', 'UHF')
    PsiMod.scf()
    return PsiMod.omp2()
    PsiMod.set_global_option('REFERENCE', oldref)

def run_scf(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a self-consistent-field theory (HF & DFT) calculation.

    """
    lowername = name.lower()

    g_user_ref = PsiMod.get_global_option('REFERENCE')
    bg_user_ref = PsiMod.has_global_option_changed('REFERENCE')

    user_fctl_scf = PsiMod.get_option('SCF', 'DFT_FUNCTIONAL')
    g_user_fctl = PsiMod.get_global_option('DFT_FUNCTIONAL')
    l_user_fctl_scf = PsiMod.get_local_option('SCF', 'DFT_FUNCTIONAL')
    bg_user_fctl = PsiMod.has_global_option_changed('DFT_FUNCTIONAL')
    bl_user_fctl_scf = PsiMod.has_local_option_changed('SCF', 'DFT_FUNCTIONAL')

    user_scftype_scf = PsiMod.get_option('SCF', 'SCF_TYPE')
    g_user_scftype = PsiMod.get_global_option('SCF_TYPE')
    l_user_scftype_scf = PsiMod.get_local_option('SCF', 'SCF_TYPE')
    bg_user_scftype = PsiMod.has_global_option_changed('SCF_TYPE')
    bl_user_scftype_scf = PsiMod.has_local_option_changed('SCF', 'SCF_TYPE')

    if lowername == 'df-scf':
        PsiMod.set_local_option('SCF', 'SCF_TYPE', 'DF')
    elif lowername == 'hf':
        if PsiMod.get_option('SCF', 'REFERENCE') == 'RKS':
            PsiMod.set_global_option('REFERENCE', 'RHF')
        elif PsiMod.get_option('SCF', 'REFERENCE') == 'UKS':
            PsiMod.set_global_option('REFERENCE', 'UHF')
        else:
            pass
    elif lowername == 'rhf':
        PsiMod.set_global_option('REFERENCE', 'RHF')
    elif lowername == 'uhf':
        PsiMod.set_global_option('REFERENCE', 'UHF')
    elif lowername == 'rohf':
        PsiMod.set_global_option('REFERENCE', 'ROHF')
    elif lowername == 'rscf':
        if (len(PsiMod.get_option('SCF', 'DFT_FUNCTIONAL')) > 0) or PsiMod.get_option('SCF', 'DFT_CUSTOM_FUNCTIONAL') is not None:
            PsiMod.set_global_option('REFERENCE', 'RKS')
        else:
            PsiMod.set_global_option('REFERENCE', 'RHF')
    elif lowername == 'uscf':
        if (len(PsiMod.get_option('SCF', 'DFT_FUNCTIONAL')) > 0) or PsiMod.get_option('SCF', 'DFT_CUSTOM_FUNCTIONAL') is not None:
            PsiMod.set_global_option('REFERENCE', 'UKS')
        else:
            PsiMod.set_global_option('REFERENCE', 'UHF')
    elif lowername == 'roscf':
        if (len(PsiMod.get_option('SCF', 'DFT_FUNCTIONAL')) > 0) or PsiMod.get_option('SCF', 'DFT_CUSTOM_FUNCTIONAL') is not None:
            raise ValidationError('ROHF reference for DFT is not available.')
        else:
            PsiMod.set_global_option('REFERENCE', 'ROHF')

    returnvalue = scf_helper(name, **kwargs)

    PsiMod.set_global_option('REFERENCE', g_user_ref)
    if not bg_user_ref:
        PsiMod.revoke_global_option_changed('REFERENCE')

    PsiMod.set_global_option('DFT_FUNCTIONAL', g_user_fctl)
    if not bg_user_fctl:
        PsiMod.revoke_global_option_changed('DFT_FUNCTIONAL')
    PsiMod.set_local_option('SCF', 'DFT_FUNCTIONAL', l_user_fctl_scf)
    if not bl_user_fctl_scf:
        PsiMod.revoke_local_option_changed('SCF', 'DFT_FUNCTIONAL')

    PsiMod.set_global_option('SCF_TYPE', g_user_scftype)
    if not bg_user_scftype:
        PsiMod.revoke_global_option_changed('SCF_TYPE')
    PsiMod.set_local_option('SCF', 'SCF_TYPE', l_user_scftype_scf)
    if not bl_user_scftype_scf:
        PsiMod.revoke_local_option_changed('SCF', 'SCF_TYPE')

    return returnvalue


def run_scf_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a SCF gradient calculation.

    """

    run_scf(name, **kwargs)

    if (PsiMod.get_global_option('SCF_TYPE') == 'DF'):
        PsiMod.scfgrad()
    else:
        PsiMod.deriv()


def run_libfock(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a calculation through libfock, namely RCPHF,
    RCIS, RTDHF, RTDA, and RTDDFT.

    """
    if (name.lower() == 'cphf'):
        PsiMod.set_global_option('MODULE', 'RCPHF')
    if (name.lower() == 'cis'):
        PsiMod.set_global_option('MODULE', 'RCIS')
    if (name.lower() == 'tdhf'):
        PsiMod.set_global_option('MODULE', 'RTDHF')
    if (name.lower() == 'cpks'):
        PsiMod.set_global_option('MODULE', 'RCPKS')
    if (name.lower() == 'tda'):
        PsiMod.set_global_option('MODULE', 'RTDA')
    if (name.lower() == 'tddft'):
        PsiMod.set_global_option('MODULE', 'RTDDFT')

    PsiMod.libfock()


def run_mcscf(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a multiconfigurational self-consistent-field calculation.

    """
    return PsiMod.mcscf()


def scf_helper(name, **kwargs):
    """Function serving as helper to SCF, choosing whether to cast
    up or just run SCF with a standard guess. This preserves
    previous SCF options set by other procedures (e.g., SAPT
    output file types for SCF).

    """

    g_user_basis = PsiMod.get_global_option('BASIS')
    bg_user_basis = PsiMod.has_global_option_changed('BASIS')

    g_user_puream = PsiMod.get_global_option('PUREAM')
    bg_user_puream = PsiMod.has_global_option_changed('PUREAM')

    user_scftype_scf = PsiMod.get_option('SCF', 'SCF_TYPE')
    g_user_scftype = PsiMod.get_global_option('SCF_TYPE')
    l_user_scftype_scf = PsiMod.get_local_option('SCF', 'SCF_TYPE')
    bg_user_scftype = PsiMod.has_global_option_changed('SCF_TYPE')
    bl_user_scftype_scf = PsiMod.has_local_option_changed('SCF', 'SCF_TYPE')

    user_scf_dfbasisscf = PsiMod.get_option('SCF', 'DF_BASIS_SCF')
    b_user_scf_dfbasisscf = PsiMod.has_local_option_changed('SCF', 'DF_BASIS_SCF')
    user_scf_guess = PsiMod.get_option('SCF', 'GUESS')
    b_user_scf_guess = PsiMod.has_local_option_changed('SCF', 'GUESS')
    user_scf_dfintsio = PsiMod.get_option('SCF', 'DF_INTS_IO')
    b_user_scf_dfintsio = PsiMod.has_local_option_changed('SCF', 'DF_INTS_IO')

    cast = False
    if 'cast_up' in kwargs:
        cast = kwargs.pop('cast_up')
        if input.yes.match(str(cast)):
            cast = True
        elif input.no.match(str(cast)):
            cast = False
        
        if user_scftype_scf == 'DF':
            castdf = True
        else:
            castdf = False
    
        if 'cast_up_df' in kwargs:
            castdf = kwargs.pop('cast_up_df')
            if input.yes.match(str(castdf)):
                castdf = True
            elif input.no.match(str(castdf)):
                castdf = False

    precallback = None
    if 'precallback' in kwargs:
        precallback = kwargs.pop('precallback')

    postcallback = None
    if 'postcallback' in kwargs:
        postcallback = kwargs.pop('postcallback')

    # Hack to ensure cartesian or pure are used throughout
    # Note that can't query PUREAM option directly, as it only
    #   reflects user changes to value, so load basis and
    #   read effective PUREAM setting off of it
    # This if statement is only to handle generic, unnamed basis cases
    #   (like mints2) that should be departing soon
    if bg_user_basis:
        PsiMod.set_global_option('BASIS', g_user_basis)
        PsiMod.set_global_option('PUREAM', PsiMod.MintsHelper().basisset().has_puream())

    if (cast):

        if input.yes.match(str(cast)):
            guessbasis = '3-21G'
        else:
            guessbasis = cast

        if (castdf):
            if input.yes.match(str(castdf)):
                guessbasisdf = corresponding_jkfit(guessbasis)
            else:
                guessbasisdf = castdf

        # Switch to the guess namespace
        namespace = PsiMod.IO.get_default_namespace()
        PsiMod.IO.set_default_namespace((namespace + '.guess'))

        # Setup initial SCF
        PsiMod.set_global_option('BASIS', guessbasis)
        if (castdf):
            PsiMod.set_local_option('SCF', 'SCF_TYPE', 'DF')
            PsiMod.set_local_option('SCF', 'DF_INTS_IO', 'none')
            PsiMod.set_local_option('SCF', 'DF_BASIS_SCF', guessbasisdf)

        # Print some info about the guess
        PsiMod.print_out('\n')
        banner('Guess SCF, %s Basis' % (guessbasis))
        PsiMod.print_out('\n')

        # Perform the guess scf
        PsiMod.scf()

        # Move files to proper namespace
        PsiMod.IO.change_file_namespace(180, (namespace + '.guess'), namespace)
        PsiMod.IO.set_default_namespace(namespace)

        # DF-BASIS-SCF keyword could be seeded here from the BasisFamily
        #   defaults so that more calcs could run w/o fitting basis error.
        #   However, fitting basis defaults are currently handled in-module
        #   in several modules in addition to scf, so will wait for std soln.

        # Set to read and project, and reset bases to final ones
        PsiMod.set_global_option('BASIS', g_user_basis)
        PsiMod.set_local_option('SCF', 'SCF_TYPE', user_scftype_scf)
        PsiMod.set_local_option('SCF', 'GUESS', 'READ')
        if (user_scftype_scf == 'DF'):
            PsiMod.set_local_option('SCF', 'DF_INTS_IO', user_scf_dfintsio)
            PsiMod.set_local_option('SCF', 'DF_BASIS_SCF', user_scf_dfbasisscf)

        # Print the banner for the standard operation
        PsiMod.print_out('\n')
        banner(name.upper())
        PsiMod.print_out('\n')

        # Do the full scf
        e_scf = PsiMod.scf(precallback, postcallback)

        PsiMod.set_local_option('SCF', 'GUESS', user_scf_guess)

    else:

        e_scf = PsiMod.scf(precallback, postcallback)

    PsiMod.set_global_option('BASIS', g_user_basis)
    if not bg_user_basis:
        PsiMod.revoke_global_option_changed('BASIS')
    PsiMod.set_global_option('PUREAM', g_user_puream)
    if not bg_user_puream:
        PsiMod.revoke_global_option_changed('PUREAM')

    PsiMod.set_global_option('SCF_TYPE', g_user_scftype)
    if not bg_user_scftype:
        PsiMod.revoke_global_option_changed('SCF_TYPE')
    PsiMod.set_local_option('SCF', 'SCF_TYPE', l_user_scftype_scf)
    if not bl_user_scftype_scf:
        PsiMod.revoke_local_option_changed('SCF', 'SCF_TYPE')

    PsiMod.set_local_option('SCF', 'DF_BASIS_SCF', user_scf_dfbasisscf)
    if not b_user_scf_dfbasisscf:
        PsiMod.revoke_local_option_changed('SCF', 'DF_BASIS_SCF')
    PsiMod.set_local_option('SCF', 'GUESS', user_scf_guess)
    if not b_user_scf_guess:
        PsiMod.revoke_local_option_changed('SCF', 'GUESS')
    PsiMod.set_local_option('SCF', 'DF_INTS_IO', user_scf_dfintsio)
    if not b_user_scf_dfintsio:
        PsiMod.revoke_local_option_changed('SCF', 'DF_INTS_IO')

    return e_scf


def run_mp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a MP2 calculation.

    """
    PsiMod.set_global_option('WFN', 'MP2')

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and input.yes.match(str(kwargs['bypass_scf']))):
        run_scf('scf', **kwargs)

        # If the scf type is DF, then the AO integrals were never generated
        if PsiMod.get_option('SCF', 'SCF_TYPE') == 'DF':
            mints = PsiMod.MintsHelper()
            mints.integrals()

    PsiMod.transqt2()
    PsiMod.ccsort()
    returnvalue = PsiMod.mp2()

    PsiMod.set_global_option('WFN', 'SCF')
    PsiMod.revoke_global_option_changed('WFN')

    return returnvalue


def run_mp2_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a MP2 gradient calculation.

    """
    PsiMod.set_global_option('DERTYPE', 'FIRST')

    run_mp2(name, **kwargs)
    PsiMod.set_global_option('WFN', 'MP2')

    PsiMod.deriv()

    PsiMod.set_global_option('WFN', 'SCF')
    PsiMod.revoke_global_option_changed('WFN')
    PsiMod.set_global_option('DERTYPE', 'NONE')
    PsiMod.revoke_global_option_changed('DERTYPE')

def run_dfmp2_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a DFMP2 gradient calculation.

    """

    if 'restart_file' in kwargs:
        restartfile = kwargs.pop('restart_file')
        # Rename the checkpoint file to be consistent with psi4's file system
        psioh = PsiMod.IOManager.shared_object()
        psio = PsiMod.IO.shared_object()
        filepath = psioh.get_file_path(32)
        namespace = psio.get_default_namespace()
        pid = str(os.getpid())
        prefix = 'psi'
        targetfile = filepath + prefix + '.' + pid + '.' + namespace + '.32'
        if(PsiMod.me() == 0):
            shutil.copy(restartfile, targetfile)
    else:
        run_scf('RHF', **kwargs)

    PsiMod.print_out('\n')
    banner('DFMP2')
    PsiMod.print_out('\n')
    PsiMod.dfmp2grad()
    e_dfmp2 = PsiMod.get_variable('DF-MP2 ENERGY')
    e_scs_dfmp2 = PsiMod.get_variable('SCS-DF-MP2 ENERGY')
    if (name.upper() == 'SCS-DFMP2'):
        return e_scs_dfmp2
    elif (name.upper() == 'DF-MP2'):
        return e_dfmp2
    
def run_ccenergy(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a CCSD, CC2, and CC3 calculation.

    """
    if (name.lower() == 'ccsd'):
        PsiMod.set_global_option('WFN', 'CCSD')
    elif (name.lower() == 'ccsd(t)'):
        PsiMod.set_global_option('WFN', 'CCSD_T')
    elif (name.lower() == 'cc2'):
        PsiMod.set_global_option('WFN', 'CC2')
    elif (name.lower() == 'cc3'):
        PsiMod.set_global_option('WFN', 'CC3')
    elif (name.lower() == 'eom-cc2'):
        PsiMod.set_global_option('WFN', 'EOM_CC2')
    elif (name.lower() == 'eom-ccsd'):
        PsiMod.set_global_option('WFN', 'EOM_CCSD')
    # Call a plain energy('ccenergy') and have full control over options,
    # incl. wfn
    elif(name.lower() == 'ccenergy'):
        pass

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and input.yes.match(str(kwargs['bypass_scf']))):
        run_scf('scf', **kwargs)

        # If the scf type is DF, then the AO integrals were never generated
        if PsiMod.get_option('SCF', 'SCF_TYPE') == 'DF':
            mints = PsiMod.MintsHelper()
            mints.integrals()

    PsiMod.transqt2()
    PsiMod.ccsort()
    returnvalue = PsiMod.ccenergy()

    if (name.lower() != 'ccenergy'):
        PsiMod.set_global_option('WFN', 'SCF')
        PsiMod.revoke_global_option_changed('WFN')

    return returnvalue


def run_cc_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a CCSD and CCSD(T) gradient calculation.

    """
    PsiMod.set_global_option('DERTYPE', 'FIRST')

    run_ccenergy(name, **kwargs)
    if (name.lower() == 'ccsd'):
        PsiMod.set_global_option('WFN', 'CCSD')
    elif (name.lower() == 'ccsd(t)'):
        PsiMod.set_global_option('WFN', 'CCSD_T')

    PsiMod.cchbar()
    PsiMod.cclambda()
    PsiMod.ccdensity()
    PsiMod.deriv()

    if (name.lower() != 'ccenergy'):
        PsiMod.set_global_option('WFN', 'SCF')
        PsiMod.revoke_global_option_changed('WFN')
        PsiMod.set_global_option('DERTYPE', 'NONE')
        PsiMod.revoke_global_option_changed('DERTYPE')


def run_bccd(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a Brueckner CCD calculation.

    """
    if (name.lower() == 'bccd'):
        PsiMod.set_global_option('WFN', 'BCCD')

    # Bypass routine scf if user did something special to get it to
    # converge
    if not (('bypass_scf' in kwargs) and input.yes.match(str(kwargs['bypass_scf']))):
        run_scf('scf', **kwargs)

        # If the scf type is DF, then the AO integrals were never generated
        if PsiMod.get_option('scf', 'scf_type') == 'DF':
            mints = PsiMod.MintsHelper()
            mints.integrals()

    PsiMod.set_global_option('DELETE_TEI', 'false')

    while True:
        PsiMod.transqt2()
        PsiMod.ccsort()
        returnvalue = PsiMod.ccenergy()
        PsiMod.print_out('Brueckner convergence check: %d\n' % PsiMod.get_variable('BRUECKNER CONVERGED'))
        if (PsiMod.get_variable('BRUECKNER CONVERGED') == True):
            break

    return returnvalue


def run_bccd_t(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a Brueckner CCD(T) calculation.

    """
    PsiMod.set_global_option('WFN', 'BCCD_T')
    run_bccd(name, **kwargs)

    return PsiMod.cctriples()


def run_scf_property(name, **kwargs):

    run_scf(name, **kwargs)


def run_cc_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    all CC property calculations.

    """
    oneel_properties = ['dipole', 'quadrupole']
    twoel_properties = []
    response_properties = ['polarizability', 'rotation', 'roa']
    excited_properties = ['oscillator_strength', 'rotational_strength']

    one = []
    two = []
    response = []
    excited = []
    invalid = []

    if 'properties' in kwargs:
        properties = kwargs.pop('properties')
        properties = drop_duplicates(properties)

        for prop in properties:
            if prop in oneel_properties:
                one.append(prop)
            elif prop in twoel_properties:
                two.append(prop)
            elif prop in response_properties:
                response.append(prop)
            elif prop in excited_properties:
                excited.append(prop)
            else:
                invalid.append(prop)
    else:
        print("The \"properties\" keyword is required with the property() function.")
        exit(1)

    n_one = len(one)
    n_two = len(two)
    n_response = len(response)
    n_excited = len(excited)
    n_invalid = len(invalid)

    if (n_invalid > 0):
        print("The following properties are not currently supported: %s" % invalid)

    if (n_excited > 0 and (name.lower() != 'eom-ccsd' and name.lower() != 'eom-cc2')):
        print("Excited state CC properties require EOM-CC2 or EOM-CCSD.")
        exit(1)

    if ((name.lower() == 'eom-ccsd' or name.lower() == 'eom-cc2') and n_response > 0):
        print("Cannot (yet) compute response properties for excited states.")
        exit(1)

    if (n_one > 0 or n_two > 0) and (n_response > 0):
        print("Computing both density- and response-based properties.")

    if (name.lower() == 'ccsd'):
        PsiMod.set_global_option('WFN', 'CCSD')
        run_ccenergy('ccsd', **kwargs)
        PsiMod.set_global_option('WFN', 'CCSD')
    elif (name.lower() == 'cc2'):
        PsiMod.set_global_option('WFN', 'CC2')
        run_ccenergy('cc2', **kwargs)
        PsiMod.set_global_option('WFN', 'CC2')
    elif (name.lower() == 'eom-ccsd'):
        PsiMod.set_global_option('WFN', 'EOM_CCSD')
        run_ccenergy('eom-ccsd', **kwargs)
        PsiMod.set_global_option('WFN', 'EOM_CCSD')
    elif (name.lower() == 'eom-cc2'):
        PsiMod.set_global_option('WFN', 'EOM_CC2')
        run_ccenergy('eom-cc2', **kwargs)
        PsiMod.set_global_option('WFN', 'EOM_CC2')

    # Need cchbar for everything
    PsiMod.cchbar()

    # Need ccdensity at this point only for density-based props
    if (n_one > 0 or n_two > 0):
        if (name.lower() == 'eom-ccsd'):
            PsiMod.set_global_option('WFN', 'EOM_CCSD')
            PsiMod.set_global_option('DERTYPE', 'NONE')
            PsiMod.set_global_option('ONEPDM', 'TRUE')
            PsiMod.cceom()
        elif (name.lower() == 'eom-cc2'):
            PsiMod.set_global_option('WFN', 'EOM_CC2')
            PsiMod.set_global_option('DERTYPE', 'NONE')
            PsiMod.set_global_option('ONEPDM', 'TRUE')
            PsiMod.cceom()
        PsiMod.set_global_option('DERTYPE', 'NONE')
        PsiMod.set_global_option('ONEPDM', 'TRUE')
        PsiMod.cclambda()
        PsiMod.ccdensity()

    # Need ccresponse only for response-type props
    if (n_response > 0):
        PsiMod.set_global_option('DERTYPE', 'RESPONSE')
        PsiMod.cclambda()
        for prop in response:
            PsiMod.set_global_option('PROPERTY', prop)
            PsiMod.ccresponse()

    # Excited-state transition properties
    if (n_excited > 0):
        if (name.lower() == 'eom-ccsd'):
            PsiMod.set_global_option('WFN', 'EOM_CCSD')
        elif (name.lower() == 'eom-cc2'):
            PsiMod.set_global_option('WFN', 'EOM_CC2')
        else:
            print("Unknown excited-state CC wave function.")
            exit(1)
        PsiMod.set_global_option('DERTYPE', 'NONE')
        PsiMod.set_global_option('ONEPDM', 'TRUE')
        PsiMod.cceom()
        PsiMod.cclambda()
        PsiMod.ccdensity()

    PsiMod.set_global_option('WFN', 'SCF')
    PsiMod.revoke_global_option_changed('WFN')
    PsiMod.set_global_option('DERTYPE', 'NONE')
    PsiMod.revoke_global_option_changed('DERTYPE')


def run_eom_cc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an EOM-CC calculation, namely EOM-CC2, EOM-CCSD, and EOM-CC3.

    """
    if (name.lower() == 'eom-ccsd'):
        PsiMod.set_global_option('WFN', 'EOM_CCSD')
        run_ccenergy('ccsd', **kwargs)
        PsiMod.set_global_option('WFN', 'EOM_CCSD')
    elif (name.lower() == 'eom-cc2'):
        PsiMod.set_global_option('WFN', 'EOM_CC2')
        run_ccenergy('cc2', **kwargs)
        PsiMod.set_global_option('WFN', 'EOM_CC2')
    elif (name.lower() == 'eom-cc3'):
        PsiMod.set_global_option('WFN', 'EOM_CC3')
        run_ccenergy('cc3', **kwargs)
        PsiMod.set_global_option('WFN', 'EOM_CC3')

    PsiMod.cchbar()
    returnvalue = PsiMod.cceom()

    PsiMod.set_global_option('WFN', 'SCF')
    PsiMod.revoke_global_option_changed('WFN')

    return returnvalue


def run_eom_cc_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an EOM-CCSD gradient calculation.

    """
    PsiMod.set_global_option('DERTYPE', 'FIRST')

    if (name.lower() == 'eom-ccsd'):
        PsiMod.set_global_option('WFN', 'EOM_CCSD')
        energy = run_eom_cc(name, **kwargs)
        PsiMod.set_global_option('WFN', 'EOM_CCSD')

    PsiMod.set_global_option('WFN', 'EOM_CCSD')
    PsiMod.set_global_option('ZETA', 'FALSE')
    PsiMod.cclambda()
    PsiMod.set_global_option('XI', 'TRUE')
    PsiMod.ccdensity()
    PsiMod.set_global_option('ZETA', 'TRUE')
    PsiMod.cclambda()
    PsiMod.set_global_option('XI', 'FALSE')
    PsiMod.ccdensity()
    PsiMod.deriv()

    PsiMod.set_global_option('WFN', 'SCF')
    PsiMod.revoke_global_option_changed('WFN')
    PsiMod.set_global_option('DERTYPE', 'NONE')
    PsiMod.revoke_global_option_changed('DERTYPE')


def run_adc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an algebraic diagrammatic construction calculation.

    .. caution:: Get rid of active molecule lines- should be handled in energy.

    """
    molecule = PsiMod.get_active_molecule()
    if 'molecule' in kwargs:
        molecule = kwargs.pop('molecule')

    if not molecule:
        raise ValueNotSet('no molecule found')

    PsiMod.scf()

    return PsiMod.adc()


def run_dft(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-functional-theory calculation.

    """
    lowername = name.lower()

    user_fctl = PsiMod.get_option('SCF', 'DFT_FUNCTIONAL')
    b_user_fctl = PsiMod.has_option_changed('SCF', 'DFT_FUNCTIONAL')
    user_ref = PsiMod.get_option('SCF', 'REFERENCE')
    b_user_ref = PsiMod.has_option_changed('SCF', 'REFERENCE')


    PsiMod.set_global_option('DFT_FUNCTIONAL', lowername)

    if (user_ref == 'RHF'):
        PsiMod.set_global_option('REFERENCE', 'RKS')
    elif (user_ref == 'UHF'):
        PsiMod.set_global_option('REFERENCE', 'UKS')
    elif (user_ref == 'ROHF'):
        raise ValidationError('ROHF reference for DFT is not available.')
    elif (user_ref == 'CUHF'):
        raise ValidationError('CUHF reference for DFT is not available.')

    returnvalue = run_scf(name,**kwargs) 

    PsiMod.set_global_option('DFT_FUNCTIONAL', user_fctl)
    if not b_user_fctl:
        PsiMod.revoke_global_option_changed('DFT_FUNCTIONAL')
    PsiMod.set_global_option('REFERENCE', user_ref)
    if not b_user_ref:
        PsiMod.revoke_global_option_changed('REFERENCE')

    return returnvalue

def run_dft_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-functional-theory gradient calculation.

    """
    lowername = name.lower()

    user_fctl = PsiMod.get_global_option('DFT_FUNCTIONAL')
    user_ref = PsiMod.get_global_option('REFERENCE')

    PsiMod.set_global_option('DFT_FUNCTIONAL', lowername)

    if (user_ref == 'RHF'):
        PsiMod.set_global_option('REFERENCE', 'RKS')
    elif (user_ref == 'UHF'):
        PsiMod.set_global_option('REFERENCE', 'UKS')
    elif (user_ref == 'ROHF'):
        raise ValidationError('ROHF reference for DFT is not available.')
    elif (user_ref == 'CUHF'):
        raise ValidationError('CUHF reference for DFT is not available.')

    if (PsiMod.get_global_option('SCF_TYPE') != 'DF'):
        raise ValidationError('SCF_TYPE must be DF for DFT gradient (for now).')

    run_scf_gradient(name,**kwargs) 

    PsiMod.set_global_option('DFT_FUNCTIONAL', user_fctl)
    PsiMod.revoke_global_option_changed('DFT_FUNCTIONAL')
    PsiMod.set_global_option('REFERENCE', user_ref)
    PsiMod.revoke_global_option_changed('REFERENCE')

def run_detci(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a configuration interaction calculation, namely FCI,
    CIn, MPn, and ZAPTn.

    """
    if (name.lower() == 'zapt'):
        PsiMod.set_global_option('WFN', 'ZAPTN')
        level = kwargs['level']
        maxnvect = (level + 1) / 2 + (level + 1) % 2
        PsiMod.set_global_option('MAX_NUM_VECS', maxnvect)
        if ((level + 1) % 2):
            PsiMod.set_global_option('MPN_ORDER_SAVE', 2)
        else:
            PsiMod.set_global_option('MPN_ORDER_SAVE', 1)
    elif (name.lower() == 'mp'):
        PsiMod.set_global_option('WFN', 'DETCI')
        PsiMod.set_global_option('MPN', 'TRUE')

        level = kwargs['level']
        maxnvect = (level + 1) / 2 + (level + 1) % 2
        PsiMod.set_global_option('MAX_NUM_VECS', maxnvect)
        if ((level + 1) % 2):
            PsiMod.set_global_option('MPN_ORDER_SAVE', 2)
        else:
            PsiMod.set_global_option('MPN_ORDER_SAVE', 1)
    elif (name.lower() == 'fci'):
            PsiMod.set_global_option('WFN', 'DETCI')
            PsiMod.set_global_option('FCI', 'TRUE')
    elif (name.lower() == 'cisd'):
            PsiMod.set_global_option('WFN', 'DETCI')
            PsiMod.set_global_option('EX_LEVEL', 2)
    elif (name.lower() == 'cisdt'):
            PsiMod.set_global_option('WFN', 'DETCI')
            PsiMod.set_global_option('EX_LEVEL', 3)
    elif (name.lower() == 'cisdtq'):
            PsiMod.set_global_option('WFN', 'DETCI')
            PsiMod.set_global_option('EX_LEVEL', 4)
    elif (name.lower() == 'ci'):
        PsiMod.set_global_option('WFN', 'DETCI')
        level = kwargs['level']
        PsiMod.set_global_option('EX_LEVEL', level)
    # Call a plain energy('detci') and have full control over options
    elif(name.lower() == 'detci'):
        pass

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and input.yes.match(str(kwargs['bypass_scf']))):
        run_scf('scf', **kwargs)

        # If the scf type is DF, then the AO integrals were never generated
        if PsiMod.get_option('SCF', 'SCF_TYPE') == 'DF':
            mints = PsiMod.MintsHelper()
            mints.integrals()

    PsiMod.transqt2()
    returnvalue = PsiMod.detci()

    if (name.lower() != 'detci'):
        PsiMod.set_global_option('WFN', 'SCF')
        PsiMod.revoke_global_option_changed('WFN')
        PsiMod.set_global_option('MPN', 'FALSE')
        PsiMod.revoke_global_option_changed('MPN')
        PsiMod.set_global_option('MAX_NUM_VECS', 12)
        PsiMod.revoke_global_option_changed('MAX_NUM_VECS')
        PsiMod.set_global_option('MPN_ORDER_SAVE', 0)
        PsiMod.revoke_global_option_changed('MPN_ORDER_SAVE')
        PsiMod.set_global_option('FCI', 'FALSE')
        PsiMod.revoke_global_option_changed('FCI')
        PsiMod.set_global_option('EX_LEVEL', 2)
        PsiMod.revoke_global_option_changed('EX_LEVEL')

    return returnvalue


def run_dfmp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted MP2 calculation.

    .. caution:: Get rid of madness-era restart file

    """
    if 'restart_file' in kwargs:
        restartfile = kwargs.pop('restart_file')
        # Rename the checkpoint file to be consistent with psi4's file system
        psioh = PsiMod.IOManager.shared_object()
        psio = PsiMod.IO.shared_object()
        filepath = psioh.get_file_path(32)
        namespace = psio.get_default_namespace()
        pid = str(os.getpid())
        prefix = 'psi'
        targetfile = filepath + prefix + '.' + pid + '.' + namespace + '.32'
        if(PsiMod.me() == 0):
            shutil.copy(restartfile, targetfile)
    else:
        run_scf('RHF', **kwargs)

    PsiMod.print_out('\n')
    banner('DFMP2')
    PsiMod.print_out('\n')
    e_dfmp2 = PsiMod.dfmp2()
    e_scs_dfmp2 = PsiMod.get_variable('SCS-DF-MP2 ENERGY')
    if (name.upper() == 'SCS-DFMP2'):
        return e_scs_dfmp2
    elif (name.upper() == 'DF-MP2'):
        return e_dfmp2


def run_psimrcc(name, **kwargs):
    """Function encoding sequence of PSI module calls for a PSIMRCC computation
     using a reference from the MCSCF module

    """
    run_mcscf(name, **kwargs)
    return PsiMod.psimrcc()


def run_psimrcc_scf(name, **kwargs):
    """Function encoding sequence of PSI module calls for a PSIMRCC computation
     using a reference from the SCF module

    """

    run_scf(name, **kwargs)
    return PsiMod.psimrcc()


def run_mp2c(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a coupled MP2 calculation.

    """
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    monomerA = molecule.extract_subsets(1, 2)
    monomerA.set_name('monomerA')
    monomerB = molecule.extract_subsets(2, 1)
    monomerB.set_name('monomerB')

    ri = PsiMod.get_option('SCF', 'SCF_TYPE')
    df_ints_io = PsiMod.get_option('SCF', 'DF_INTS_IO')
    # inquire if above at all applies to dfmp2

    PsiMod.IO.set_default_namespace('dimer')
    PsiMod.set_local_option('SCF', 'SAPT', '2-dimer')
    PsiMod.print_out('\n')
    banner('Dimer HF')
    PsiMod.print_out('\n')
    PsiMod.set_global_option('DF_INTS_IO', 'SAVE')
    e_dimer = scf_helper('RHF', **kwargs)
    PsiMod.print_out('\n')
    banner('Dimer DFMP2')
    PsiMod.print_out('\n')
    e_dimer_mp2 = PsiMod.dfmp2()
    PsiMod.set_global_option('DF_INTS_IO', 'LOAD')

    activate(monomerA)
    if (ri == 'DF'):
        PsiMod.IO.change_file_namespace(97, 'dimer', 'monomerA')
    PsiMod.IO.set_default_namespace('monomerA')
    PsiMod.set_local_option('SCF', 'SAPT', '2-monomer_A')
    PsiMod.print_out('\n')
    banner('Monomer A HF')
    PsiMod.print_out('\n')
    e_monomerA = scf_helper('RHF', **kwargs)
    PsiMod.print_out('\n')
    banner('Monomer A DFMP2')
    PsiMod.print_out('\n')
    e_monomerA_mp2 = PsiMod.dfmp2()

    activate(monomerB)
    if (ri == 'DF'):
        PsiMod.IO.change_file_namespace(97, 'monomerA', 'monomerB')
    PsiMod.IO.set_default_namespace('monomerB')
    PsiMod.set_local_option('SCF', 'SAPT', '2-monomer_B')
    PsiMod.print_out('\n')
    banner('Monomer B HF')
    PsiMod.print_out('\n')
    e_monomerB = scf_helper('RHF', **kwargs)
    PsiMod.print_out('\n')
    banner('Monomer B DFMP2')
    PsiMod.print_out('\n')
    e_monomerB_mp2 = PsiMod.dfmp2()
    PsiMod.set_global_option('DF_INTS_IO', df_ints_io)

    PsiMod.IO.change_file_namespace(121, 'monomerA', 'dimer')
    PsiMod.IO.change_file_namespace(122, 'monomerB', 'dimer')

    activate(molecule)
    PsiMod.IO.set_default_namespace('dimer')
    PsiMod.set_local_option('SAPT', 'E_CONVERGENCE', 10e-10)
    PsiMod.set_local_option('SAPT', 'D_CONVERGENCE', 10e-10)
    PsiMod.set_local_option('SAPT', 'SAPT_LEVEL', 'MP2C')
    PsiMod.print_out('\n')
    banner('MP2C')
    PsiMod.print_out('\n')

    PsiMod.set_variable('MP2C DIMER MP2 ENERGY', e_dimer_mp2)
    PsiMod.set_variable('MP2C MONOMER A MP2 ENERGY', e_monomerA_mp2)
    PsiMod.set_variable('MP2C MONOMER B MP2 ENERGY', e_monomerB_mp2)

    e_sapt = PsiMod.sapt()
    return e_sapt


def run_sapt(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a SAPT calculation of any level.

    """

    molecule = PsiMod.get_active_molecule()
    user_pg = molecule.schoenflies_symbol()
    molecule.reset_point_group('c1')
    molecule.fix_orientation(True)
    molecule.update_geometry()

    nfrag = molecule.nfragments() 
    if nfrag != 2: 
        raise ValidationError('SAPT requires active molecule to have 2 fragments, not %s.' % (nfrag)) 

    sapt_basis = 'dimer'
    if 'sapt_basis' in kwargs:
        sapt_basis = kwargs.pop('sapt_basis')
    sapt_basis = sapt_basis.lower()

    if (sapt_basis == 'dimer'):
        molecule.update_geometry()
        monomerA = molecule.extract_subsets(1, 2)
        monomerA.set_name('monomerA')
        monomerB = molecule.extract_subsets(2, 1)
        monomerB.set_name('monomerB')
    elif (sapt_basis == 'monomer'):
        molecule.update_geometry()
        monomerA = molecule.extract_subsets(1)
        monomerA.set_name('monomerA')
        monomerB = molecule.extract_subsets(2)
        monomerB.set_name('monomerB')

    ri = PsiMod.get_option('SCF', 'SCF_TYPE')
    df_ints_io = PsiMod.get_option('SCF', 'DF_INTS_IO')
    # inquire if above at all applies to dfmp2

    PsiMod.IO.set_default_namespace('dimer')
    PsiMod.set_local_option('SCF', 'SAPT', '2-dimer')
    PsiMod.print_out('\n')
    banner('Dimer HF')
    PsiMod.print_out('\n')
    if (sapt_basis == 'dimer'):
        PsiMod.set_global_option('DF_INTS_IO', 'SAVE')
    e_dimer = scf_helper('RHF', **kwargs)
    if (sapt_basis == 'dimer'):
        PsiMod.set_global_option('DF_INTS_IO', 'LOAD')

    activate(monomerA)
    if (ri == 'DF' and sapt_basis == 'dimer'):
        PsiMod.IO.change_file_namespace(97, 'dimer', 'monomerA')
    PsiMod.IO.set_default_namespace('monomerA')
    PsiMod.set_local_option('SCF', 'SAPT', '2-monomer_A')
    PsiMod.print_out('\n')
    banner('Monomer A HF')
    PsiMod.print_out('\n')
    e_monomerA = scf_helper('RHF', **kwargs)

    activate(monomerB)
    if (ri == 'DF' and sapt_basis == 'dimer'):
        PsiMod.IO.change_file_namespace(97, 'monomerA', 'monomerB')
    PsiMod.IO.set_default_namespace('monomerB')
    PsiMod.set_local_option('SCF', 'SAPT', '2-monomer_B')
    PsiMod.print_out('\n')
    banner('Monomer B HF')
    PsiMod.print_out('\n')
    e_monomerB = scf_helper('RHF', **kwargs)
    PsiMod.set_global_option('DF_INTS_IO', df_ints_io)

    PsiMod.IO.change_file_namespace(121, 'monomerA', 'dimer')
    PsiMod.IO.change_file_namespace(122, 'monomerB', 'dimer')

    activate(molecule)
    PsiMod.IO.set_default_namespace('dimer')
    PsiMod.set_local_option('SAPT', 'E_CONVERGENCE', 10e-10)
    PsiMod.set_local_option('SAPT', 'D_CONVERGENCE', 10e-10)
    if (name.lower() == 'sapt0'):
        PsiMod.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT0')
    elif (name.lower() == 'sapt2'):
        PsiMod.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2')
    elif (name.lower() == 'sapt2+'):
        PsiMod.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+')
    elif (name.lower() == 'sapt2+(3)'):
        PsiMod.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        PsiMod.set_local_option('SAPT', 'DO_THIRD_ORDER', False)
    elif (name.lower() == 'sapt2+3'):
        PsiMod.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        PsiMod.set_local_option('SAPT', 'DO_THIRD_ORDER', True)
    PsiMod.print_out('\n')
    banner(name.upper())
    PsiMod.print_out('\n')
    e_sapt = PsiMod.sapt()

    molecule.reset_point_group(user_pg)
    molecule.update_geometry()

    return e_sapt


def run_sapt_ct(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a charge-transfer SAPT calcuation of any level.

    """
    molecule = PsiMod.get_active_molecule()
    user_pg = molecule.schoenflies_symbol()
    molecule.reset_point_group('c1')
    molecule.fix_orientation(True)
    molecule.update_geometry()

    nfrag = molecule.nfragments() 
    if nfrag != 2: 
        raise ValidationError('SAPT requires active molecule to have 2 fragments, not %s.' % (nfrag)) 

    monomerA = molecule.extract_subsets(1, 2)
    monomerA.set_name('monomerA')
    monomerB = molecule.extract_subsets(2, 1)
    monomerB.set_name('monomerB')
    molecule.update_geometry()
    monomerAm = molecule.extract_subsets(1)
    monomerAm.set_name('monomerAm')
    monomerBm = molecule.extract_subsets(2)
    monomerBm.set_name('monomerBm')

    ri = PsiMod.get_option('SCF', 'SCF_TYPE')
    df_ints_io = PsiMod.get_option('SCF', 'DF_INTS_IO')
    # inquire if above at all applies to dfmp2

    PsiMod.IO.set_default_namespace('dimer')
    PsiMod.set_local_option('SCF', 'SAPT', '2-dimer')
    PsiMod.print_out('\n')
    banner('Dimer HF')
    PsiMod.print_out('\n')
    PsiMod.set_global_option('DF_INTS_IO', 'SAVE')
    e_dimer = scf_helper('RHF', **kwargs)
    PsiMod.set_global_option('DF_INTS_IO', 'LOAD')

    activate(monomerA)
    if (ri == 'DF'):
        PsiMod.IO.change_file_namespace(97, 'dimer', 'monomerA')
    PsiMod.IO.set_default_namespace('monomerA')
    PsiMod.set_local_option('SCF', 'SAPT', '2-monomer_A')
    PsiMod.print_out('\n')
    banner('Monomer A HF (Dimer Basis)')
    PsiMod.print_out('\n')
    e_monomerA = scf_helper('RHF', **kwargs)

    activate(monomerB)
    if (ri == 'DF'):
        PsiMod.IO.change_file_namespace(97, 'monomerA', 'monomerB')
    PsiMod.IO.set_default_namespace('monomerB')
    PsiMod.set_local_option('SCF', 'SAPT', '2-monomer_B')
    PsiMod.print_out('\n')
    banner('Monomer B HF (Dimer Basis)')
    PsiMod.print_out('\n')
    e_monomerB = scf_helper('RHF', **kwargs)
    PsiMod.set_global_option('DF_INTS_IO', df_ints_io)

    activate(monomerAm)
    PsiMod.IO.set_default_namespace('monomerAm')
    PsiMod.set_local_option('SCF', 'SAPT', '2-monomer_A')
    PsiMod.print_out('\n')
    banner('Monomer A HF (Monomer Basis)')
    PsiMod.print_out('\n')
    e_monomerA = scf_helper('RHF', **kwargs)

    activate(monomerBm)
    PsiMod.IO.set_default_namespace('monomerBm')
    PsiMod.set_local_option('SCF', 'SAPT', '2-monomer_B')
    PsiMod.print_out('\n')
    banner('Monomer B HF (Monomer Basis)')
    PsiMod.print_out('\n')
    e_monomerB = scf_helper('RHF', **kwargs)

    activate(molecule)
    PsiMod.IO.set_default_namespace('dimer')
    PsiMod.set_local_option('SAPT', 'E_CONVERGENCE', 10e-10)
    PsiMod.set_local_option('SAPT', 'D_CONVERGENCE', 10e-10)
    if (name.lower() == 'sapt0-ct'):
        PsiMod.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT0')
    elif (name.lower() == 'sapt2-ct'):
        PsiMod.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2')
    elif (name.lower() == 'sapt2+-ct'):
        PsiMod.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+')
    elif (name.lower() == 'sapt2+(3)-ct'):
        PsiMod.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        PsiMod.set_local_option('SAPT', 'DO_THIRD_ORDER', False)
    elif (name.lower() == 'sapt2+3-ct'):
        PsiMod.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        PsiMod.set_local_option('SAPT', 'DO_THIRD_ORDER', True)
    PsiMod.print_out('\n')
    banner('SAPT Charge Transfer')
    PsiMod.print_out('\n')

    PsiMod.print_out('\n')
    banner('Dimer Basis SAPT')
    PsiMod.print_out('\n')
    PsiMod.IO.change_file_namespace(121, 'monomerA', 'dimer')
    PsiMod.IO.change_file_namespace(122, 'monomerB', 'dimer')
    e_sapt = PsiMod.sapt()
    CTd = PsiMod.get_variable('SAPT CT ENERGY')

    PsiMod.print_out('\n')
    banner('Monomer Basis SAPT')
    PsiMod.print_out('\n')
    PsiMod.IO.change_file_namespace(121, 'monomerAm', 'dimer')
    PsiMod.IO.change_file_namespace(122, 'monomerBm', 'dimer')
    e_sapt = PsiMod.sapt()
    CTm = PsiMod.get_variable('SAPT CT ENERGY')
    CT = CTd - CTm

    PsiMod.print_out('\n\n')
    PsiMod.print_out('    SAPT Charge Transfer Analysis\n')
    PsiMod.print_out('  -----------------------------------------------------------------------------\n')
    line1 = '    SAPT Induction (Dimer Basis)      %10.4lf mH    %10.4lf kcal mol^-1\n' % (CTd * 1000.0, CTd * physconst.psi_hartree2kcalmol)
    line2 = '    SAPT Induction (Monomer Basis)    %10.4lf mH    %10.4lf kcal mol^-1\n' % (CTm * 1000.0, CTm * physconst.psi_hartree2kcalmol)
    line3 = '    SAPT Charge Transfer              %10.4lf mH    %10.4lf kcal mol^-1\n\n' % (CT * 1000.0, CT * physconst.psi_hartree2kcalmol)
    PsiMod.print_out(line1)
    PsiMod.print_out(line2)
    PsiMod.print_out(line3)
    PsiMod.set_variable('SAPT CT ENERGY', CT)

    molecule.reset_point_group(user_pg)
    molecule.update_geometry()

    return e_sapt


def run_mrcc(name, **kwargs):
    """Function that prepares environment and input files
    for a calculation calling Kallay's MRCC code.

    """
    # TODO: Check to see if we really need to run the SCF code.
    run_scf(name, **kwargs)

    # The parse_arbitrary_order method provides us the following information
    # We require that level be provided. level is a dictionary
    # of settings to be passed to PsiMod.mrcc
    if not('level' in kwargs):
        raise ValidationError('level parameter was not provided.')

    level = kwargs['level']

    # Fullname is the string we need to search for in iface
    fullname = level['fullname']

    # User can provide 'keep' to the method.
    # When provided, do not delete the MRCC scratch directory.
    keep = False
    if 'keep' in kwargs:
        keep = kwargs['keep']

    # Save current directory location
    current_directory = os.getcwd()

    # Need to move to the scratch directory, perferrably into a separate directory in that location
    psi_io = PsiMod.IOManager.shared_object()
    os.chdir(psi_io.get_default_path())

    # Make new directory specifically for mrcc
    mrcc_tmpdir = 'mrcc_' + str(os.getpid())
    if 'path' in kwargs:
        mrcc_tmpdir = kwargs['path']

    # Check to see if directory already exists, if not, create.
    if os.path.exists(mrcc_tmpdir) == False:
        os.mkdir(mrcc_tmpdir)

    # Move into the new directory
    os.chdir(mrcc_tmpdir)

    # Generate integrals and input file (dumps files to the current directory)
    PsiMod.mrcc_generate_input(level)

    # Load the fort.56 file
    # and dump a copy into the outfile
    PsiMod.print_out('\n===== Begin fort.56 input for MRCC ======\n')
    PsiMod.print_out(open('fort.56', 'r').read())
    PsiMod.print_out('===== End   fort.56 input for MRCC ======\n')

    # Close output file
    PsiMod.close_outfile()

    # Modify the environment:
    #    PGI Fortan prints warning to screen if STOP is used
    os.environ['NO_STOP_MESSAGE'] = '1'

    # Obtain user's OMP_NUM_THREADS so that we don't blow it away.
    omp_num_threads_found = 'OMP_NUM_THREADS' in os.environ
    if omp_num_threads_found == True:
        omp_num_threads_user = os.environ['OMP_NUM_THREADS']

    # If the user provided MRCC_OMP_NUM_THREADS set the environ to it
    if PsiMod.has_option_changed('MRCC', 'MRCC_OMP_NUM_THREADS') == True:
        os.environ['OMP_NUM_THREADS'] = str(PsiMod.get_option('MRCC', 'MRCC_OMP_NUM_THREADS'))

    # Call dmrcc, directing all screen output to the output file
    try:
        if PsiMod.outfile_name() == 'stdout':
            retcode = subprocess.call('dmrcc', shell=True)
        else:
            retcode = subprocess.call('dmrcc >> ' + current_directory + '/' + PsiMod.outfile_name(), shell=True)

        if retcode < 0:
            print('MRCC was terminated by signal %d' % -retcode, file=sys.stderr)
            exit(1)
        elif retcode > 0:
            print('MRCC errored %d' % retcode, file=sys.stderr)
            exit(1)

    except OSError as e:
        print('Execution failed: %s' % e, file=sys.stderr)
        exit(1)

    # Restore the OMP_NUM_THREADS that the user set.
    if omp_num_threads_found == True:
        if PsiMod.has_option_changed('MRCC', 'MRCC_OMP_NUM_THREADS') == True:
            os.environ['OMP_NUM_THREADS'] = omp_num_threads_user

    # Scan iface file and grab the file energy.
    e = 0.0
    for line in file('iface'):
        if fullname in line:
            fields = line.split()
            e = float(fields[5])

    PsiMod.set_variable('CURRENT ENERGY', e)
    PsiMod.set_variable(fullname + ' ENERGY', e)

    # Load the iface file
    iface = open('iface', 'r')
    iface_contents = iface.read()

    # Delete mrcc tempdir
    os.chdir('..')
    try:
        # Delete unless we're told not to
        if (keep == False and not('path' in kwargs)):
            shutil.rmtree(mrcc_tmpdir)
    except OSerror as e:
        print('Unable to remove MRCC temporary directory %s' % e, file=sys.stderr)
        exit(1)

    # Revert to previous current directory location
    os.chdir(current_directory)

    # Reopen output file
    PsiMod.reopen_outfile()

    # If we're told to keep the files or the user provided a path, do nothing.
    if (keep != False or ('path' in kwargs)):
        PsiMod.print_out('\nMRCC scratch files have been kept.\n')
        PsiMod.print_out('They can be found in ' + mrcc_tmpdir)

    # Dump iface contents to output
    PsiMod.print_out('\n')
    banner('Full results from MRCC')
    PsiMod.print_out('\n')
    PsiMod.print_out(iface_contents)

    return e

def run_cepa(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a cepa-like calculation.

    >>> energy('cepa(1)')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # override symmetry if integral direct
    if PsiMod.get_global_option('CEPA_VABCD_DIRECT'):
       molecule = PsiMod.get_active_molecule()
       molecule.update_geometry()
       molecule.reset_point_group('c1')
       #molecule.fix_orientation(1)
       #molecule.update_geometry()

    # throw an exception for open-shells
    if (PsiMod.get_global_option('reference') != 'RHF' ):
       PsiMod.print_out("\n")
       PsiMod.print_out("Error: %s requires \"reference rhf\".\n" % lowername )
       PsiMod.print_out("\n")
       sys.exit(1)

    # what type of cepa?
    if (lowername == 'cepa(0)'):
        PsiMod.set_global_option('cepa_level', 'cepa0')
    if (lowername == 'cepa(1)'):
        PsiMod.set_global_option('cepa_level', 'cepa1')
    if (lowername == 'cepa(2)'):
        #PsiMod.set_global_option('cepa_level', 'cepa2')
        # throw an exception for cepa(2)
        PsiMod.print_out("\n")
        PsiMod.print_out("Error: %s not implemented\n" % lowername )
        PsiMod.print_out("\n")
    if (lowername == 'cepa(3)'):
        PsiMod.set_global_option('cepa_level', 'cepa3')
    if (lowername == 'sdci'):
        PsiMod.set_global_option('cepa_level', 'cisd')
    if (lowername == 'dci'):
        PsiMod.set_global_option('cepa_level', 'cisd')
        user_no_singles = PsiMod.get_global_option('cepa_no_singles')
        PsiMod.set_global_option('cepa_no_singles', 1)
    if (lowername == 'acpf'):
        PsiMod.set_global_option('cepa_level', 'acpf')
    if (lowername == 'aqcc'):
        PsiMod.set_global_option('cepa_level', 'aqcc')

    PsiMod.set_global_option('WFN', 'CCSD')
    run_scf('scf', **kwargs)

    # If the scf type is DF, then the AO integrals were never generated
    if (PsiMod.get_global_option('scf_type') == 'DF' or PsiMod.get_option('scf','scf_type') == 'DF'):
       mints = PsiMod.MintsHelper()
       mints.integrals()
   
    # only call transqt2() if (ac|bd) is not integral direct
    #if (PsiMod.get_global_option('cepa_vabcd_direct') == False):
    # never call transqt2 since the switch to libtrans
    #   PsiMod.transqt2()

    PsiMod.cepa()

    # if dci, make sure we didn't overwrite the cepa_no_singles options
    if (lowername == 'dci'):
       PsiMod.set_global_option('cepa_no_singles', user_no_singles)
       PsiMod.revoke_global_option_changed('cepa_no_singles')

    return PsiMod.get_variable("CURRENT ENERGY")

# General wrapper for property computations
def run_property(name, **kwargs):

    junk = 1
    return junk

def run_b2plyp(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a B2PLYP double-hybrid density-functional-theory calculation.

    """

    # use with c_alpha = 0.0 in functional.py
    #     with Rob's current normalization in superfunctional.cc

    user_fctl = PsiMod.get_local_option('SCF', 'DFT_FUNCTIONAL')
    b_user_fctl = PsiMod.has_option_changed('DFT_FUNCTIONAL')
    user_ref = PsiMod.get_local_option('SCF', 'REFERENCE')
    b_user_ref = PsiMod.has_option_changed('REFERENCE')

    PsiMod.set_global_option('DFT_FUNCTIONAL', 'b2plyp_xc')

    if (user_ref == 'RHF'):
        PsiMod.set_global_option('REFERENCE', 'RKS')
    elif (user_ref == 'UHF'):
        PsiMod.set_global_option('REFERENCE', 'UKS')
    elif (user_ref == 'ROHF'):
        raise ValidationError('ROHF reference for DFT is not available.')
    elif (user_ref == 'CUHF'):
        raise ValidationError('CUHF reference for DFT is not available.')

    e_dft    = run_scf(name, **kwargs) 
    PsiMod.dfmp2()
    e_dhdft  = e_dft + 0.27 * PsiMod.get_variable("DF-MP2 CORRELATION ENERGY")

    PsiMod.set_global_option('DFT_FUNCTIONAL', user_fctl)
    if not b_user_fctl:
        PsiMod.revoke_global_option_changed('DFT_FUNCTIONAL')
    PsiMod.set_global_option('REFERENCE', user_ref)
    if not b_user_ref:
        PsiMod.revoke_global_option_changed('REFERENCE')

    PsiMod.set_variable('Double Hybrid Energy', e_dhdft)

    return e_dhdft


def run2_b2plyp(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a B2PLYP double-hybrid density-functional-theory calculation.

    """
    lowername = name.lower()

    # use with c_alpha = 0.27 in functional.py
    #     and with normalization suppressed @268 in superfunctional.cc

    PsiMod.set_global_option('REFERENCE', 'RKS')
    fun = build_superfunctional('b2plyp_xc',5000,1)
    PsiMod.set_global_option_python('dft_custom_functional',fun)
    e_dft        = PsiMod.scf()
    e_dfmp2      = PsiMod.dfmp2()
    e_dfmp2_corr = PsiMod.get_variable("DF-MP2 CORRELATION ENERGY")
    e_dhdft      = e_dft + fun.c_alpha() * e_dfmp2_corr
    PsiMod.set_variable('Double Hybrid Energy', e_dhdft)
    return e_dhdft


def run3_b2plyp(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a B2PLYP double-hybrid density-functional-theory calculation.

    """
    lowername = name.lower()

    # use with c_alpha = 0.0 in functional.py
    #     with Rob's current normalization in superfunctional.cc

    PsiMod.set_global_option('REFERENCE', 'RKS')
    fun = build_superfunctional('b2plyp_xc',5000,1)
    PsiMod.set_global_option_python('dft_custom_functional',fun)
    e_dft        = PsiMod.scf()
    e_dfmp2      = PsiMod.dfmp2()
    e_dfmp2_corr = PsiMod.get_variable("DF-MP2 CORRELATION ENERGY")
    e_dhdft      = e_dft + 0.27 * e_dfmp2_corr
    PsiMod.set_variable('Double Hybrid Energy', e_dhdft)
    return e_dhdft


def run_pbe0_2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a B2PLYP double-hybrid density-functional-theory calculation.

    """
    lowername = name.lower()

    # use with c_alpha = 0.0 in functional.py
    #     with Rob's current normalization in superfunctional.cc

    PsiMod.set_global_option('REFERENCE', 'RKS')
    PsiMod.set_global_option('DFT_FUNCTIONAL', 'pbe0-2_xc')
    e_dft        = PsiMod.scf()
    e_dfmp2      = PsiMod.dfmp2()
    e_dfmp2_corr = PsiMod.get_variable("DF-MP2 CORRELATION ENERGY")
    e_dhdft      = e_dft + 0.5 * e_dfmp2_corr
    PsiMod.set_variable('Double Hybrid Energy', e_dhdft)
    return e_dhdft
