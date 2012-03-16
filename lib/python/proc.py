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


def run_dcft(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density cumulant functional theory calculation.

    """
    PsiMod.scf()
    return PsiMod.dcft()


def run_scf(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a self-consistent-field theory (HF & DFT) calculation.

    """

    return scf_helper(name, **kwargs)


def run_scf_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a SCF gradient calculation.

    """

    run_scf(name, **kwargs)
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
    cast = False
    if 'cast_up' in kwargs:
        cast = kwargs.pop('cast_up')

    precallback = None
    if 'precallback' in kwargs:
        precallback = kwargs.pop('precallback')

    postcallback = None
    if 'postcallback' in kwargs:
        postcallback = kwargs.pop('postcallback')

    if (cast):

        if input.yes.match(str(cast)):
            custom = 'Default'
            guessbasis = '3-21G'
        else:
            custom = 'Custom'
            guessbasis = cast

        # Hack to ensure cartesian or pure are used throughout
        # This touches the option, as if the user set it
        puream = PsiMod.get_global_option('PUREAM')
        PsiMod.set_global_option('PUREAM', puream)

        # Switch to the guess namespace
        namespace = PsiMod.IO.get_default_namespace()
        PsiMod.IO.set_default_namespace((namespace + '.guess'))

        # Are we in a DF algorithm here?
        scf_type = PsiMod.get_option('SCF_TYPE')
        guess_type = PsiMod.get_option('GUESS')
        df_basis_scf = PsiMod.get_option('DF_BASIS_SCF')
        df_ints = PsiMod.get_option('DF_INTS_IO')

        # Which basis is the final one
        basis = PsiMod.get_option('BASIS')

        # Setup initial SCF
        PsiMod.set_local_option('SCF', 'BASIS', guessbasis)
        if (scf_type == 'DF'):
            PsiMod.set_local_option('SCF', 'DF_BASIS_SCF', 'cc-pvdz-ri')
            PsiMod.set_global_option('DF_INTS_IO', 'none')

        # Print some info about the guess
        PsiMod.print_out('\n')
        banner('Guess SCF, %s Basis' % (guessbasis))
        PsiMod.print_out('\n')

        # Perform the guess scf
        PsiMod.scf()

        # Move files to proper namespace
        PsiMod.IO.change_file_namespace(180, (namespace + '.guess'), namespace)
        PsiMod.IO.set_default_namespace(namespace)

        # Set to read and project, and reset bases
        PsiMod.set_local_option('SCF', 'GUESS', 'READ')
        PsiMod.set_local_option('SCF', 'BASIS', basis)
        if (scf_type == 'DF'):
            PsiMod.set_local_option('SCF', 'DF_BASIS_SCF', df_basis_scf)
            PsiMod.set_global_option('DF_INTS_IO', df_ints)

        # Print the banner for the standard operation
        PsiMod.print_out('\n')
        banner(name.upper())
        PsiMod.print_out('\n')

        # Do the full scf
        e_scf = PsiMod.scf(precallback, postcallback)

        PsiMod.set_local_option('SCF', 'GUESS', guess_type)

    else:

        e_scf = PsiMod.scf(precallback, postcallback)

    return e_scf


def run_mp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a MP2 calculation.

    """
    PsiMod.set_global_option('WFN', 'MP2')

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and input.yes.match(str(kwargs['bypass_scf']))):
        run_scf('scf', **kwargs)

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


def run_cc_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    all CC property calculation.

    """


    PsiMod.set_global_option('DERTYPE', 'RESPONSE')

    if (name.lower() == 'ccsd'):
        PsiMod.set_global_option('WFN', 'CCSD')
        run_ccenergy('ccsd', **kwargs)
        PsiMod.set_global_option('WFN', 'CCSD')
    elif (name.lower() == 'cc2'):
        PsiMod.set_global_option('WFN', 'CC2')
        run_ccenergy('cc2', **kwargs)
        PsiMod.set_global_option('WFN', 'CC2')

    PsiMod.cchbar()
    PsiMod.cclambda()
    PsiMod.ccresponse()

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
    elif (name.upper() == 'DFMP2'):
        return e_dfmp2


def run_mp2drpa(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a MP2-DRPA calculation.

    """
    e_scf = run_scf('RHF', **kwargs)

    PsiMod.print_out('\n')
    banner('DFMP2')
    PsiMod.print_out('\n')
    e_dfmp2 = PsiMod.dfmp2()

    PsiMod.print_out('\n')
    banner('dRPA')
    PsiMod.print_out('\n')
    PsiMod.dfcc()
    e_delta = PsiMod.get_variable('RPA SCALED DELTA ENERGY')

    PsiMod.print_out('\n')
    banner('MP2/dRPA Analysis')
    PsiMod.print_out('\n')

    PsiMod.print_out('  Reference Energy               %20.14f\n' % (e_scf))
    PsiMod.print_out('  MP2 Correlation Energy         %20.14f\n' % (e_dfmp2 - e_scf))
    PsiMod.print_out('  MP2 Total Energy               %20.14f\n' % (e_dfmp2))
    PsiMod.print_out('  Scaled dRPA/MP2J Delta Energy  %20.14f\n' % (e_delta))
    PsiMod.print_out('  Total MP2/dRPA Energy          %20.14f\n\n' % (e_dfmp2 + e_delta))

    return e_dfmp2 + e_delta


def run_dfcc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted coupled-cluster calculation.

    """
    run_scf('RHF', **kwargs)

    PsiMod.print_out('\n')
    banner('DF-CC')
    PsiMod.print_out('\n')
    e_dfcc = PsiMod.dfcc()
    return e_dfcc

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

    ri = PsiMod.get_option('SCF_TYPE')
    df_ints_io = PsiMod.get_option('DF_INTS_IO')

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

    ri = PsiMod.get_option('SCF_TYPE')
    df_ints_io = PsiMod.get_option('DF_INTS_IO')

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
    return e_sapt


def run_sapt_ct(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a charge-transfer SAPT calcuation of any level.

    """
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    monomerA = molecule.extract_subsets(1, 2)
    monomerA.set_name('monomerA')
    monomerB = molecule.extract_subsets(2, 1)
    monomerB.set_name('monomerB')
    molecule.update_geometry()
    monomerAm = molecule.extract_subsets(1)
    monomerAm.set_name('monomerAm')
    monomerBm = molecule.extract_subsets(2)
    monomerBm.set_name('monomerBm')

    ri = PsiMod.get_option('SCF_TYPE')
    df_ints_io = PsiMod.get_option('DF_INTS_IO')

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
    if PsiMod.has_option_changed('MRCC_OMP_NUM_THREADS') == True:
        os.environ['OMP_NUM_THREADS'] = str(PsiMod.get_option('MRCC_OMP_NUM_THREADS'))

    # Call dmrcc, directing all screen output to the output file
    try:
        if PsiMod.outfile_name() == 'stdout':
            retcode = subprocess.call('dmrcc', shell=True)
        else:
            retcode = subprocess.call('dmrcc >> ' + current_directory + '/' + PsiMod.outfile_name(), shell=True)

        if retcode < 0:
            print >>sys.stderr, 'MRCC was terminated by signal', -retcode
            exit(1)
        elif retcode > 0:
            print >>sys.stderr, 'MRCC errored', retcode
            exit(1)

    except OSError, e:
        print >>sys.stderr, 'Execution failed:', e
        exit(1)

    # Restore the OMP_NUM_THREADS that the user set.
    if omp_num_threads_found == True:
        if PsiMod.has_option_changed('MRCC_OMP_NUM_THREADS') == True:
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
    except OSerror, e:
        print >>sys.stderr, 'Unable to remove MRCC temporary directory', e
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

# General wrapper for property computations
def run_property(name, **kwargs):

    junk = 1
    return junk
