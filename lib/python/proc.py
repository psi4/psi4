import PsiMod
import shutil
import os
from driver import *
from molecule import *
from text import *

def run_dcft(name, **kwargs):
    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
      molecule = kwargs.pop('molecule')

    if not molecule:
      raise ValueNotSet("no molecule found")

    PsiMod.scf()
    return PsiMod.dcft()

def run_scf(name, **kwargs):
    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')

    if not molecule:
        raise ValueNotSet("no molecule found")

    molecule.update_geometry()
    return scf_helper(name, **kwargs)

def run_scf_gradient(name, **kwargs):
    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')

    if not molecule:
        raise ValueNotSet("no molecule found")

    molecule.update_geometry()

    scf_helper(name, **kwargs)
    PsiMod.deriv()


# SCF helper chooses whether to cast up or just run SCF
# with a standard guess. This preserves previous SCF options
# set by other procedures (eg. SAPT output file types for SCF)

def run_mcscf(**kwargs):
    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
      molecule = kwargs.pop('molecule')

    if not molecule:
      raise ValueNotSet("no molecule found")

    return PsiMod.mcscf()

# SCF helper chooses whether to cast up or just run SCF
# with a standard guess. This preserves previous SCF options
# set by other procedures (eg. SAPT output file types for SCF)


def scf_helper(name, **kwargs):

    cast = False
    if (kwargs.has_key('cast_up')):
        cast = kwargs.pop('cast_up')

    if (kwargs.has_key('cast_up')):
        cast = kwargs.pop('cast_up')


    precallback = None
    if (kwargs.has_key('precallback')):
        precallback = kwargs.pop('precallback')

    postcallback = None
    if (kwargs.has_key('postcallback')):
        postcallback = kwargs.pop('postcallback')

    if (cast):

        if (cast == True):
            custom = 'Default'
            guessbasis = '3-21G'
        else:
            custom = 'Custom'
            guessbasis = cast

        # Hack to ensure cartesian or pure are used throughout
        # This touches the option, as if the user set it
        puream = PsiMod.get_global_option('PUREAM')
        PsiMod.set_global_option("PUREAM",puream)

        # Switch to the guess namespace
        namespace = PsiMod.IO.get_default_namespace()
        PsiMod.IO.set_default_namespace((namespace + ".guess"))

        # Are we in a DF algorithm here?
        scf_type = PsiMod.get_option('SCF_TYPE')
        guess_type = PsiMod.get_option('GUESS')
        ri_basis_scf = PsiMod.get_option('RI_BASIS_SCF')
        ri_ints = PsiMod.get_option("RI_INTS_IO")

        # Which basis is the final one
        basis = PsiMod.get_option('BASIS')

        # Setup initial SCF
        PsiMod.set_local_option('SCF','BASIS',guessbasis)
        if (scf_type == 'DF'):
            PsiMod.set_local_option('SCF','RI_BASIS_SCF','cc-pvdz-ri')
            PsiMod.set_global_option('RI_INTS_IO','none')

        # Print some info about the guess
        PsiMod.print_out('\n')
        banner('Guess SCF, %s Basis' %(guessbasis))
        PsiMod.print_out('\n')

        # Perform the guess scf
        PsiMod.scf()

        # Move files to proper namespace
        PsiMod.IO.change_file_namespace(100,(namespace + ".guess"),namespace)
        PsiMod.IO.set_default_namespace(namespace)

        # Set to read and project, and reset bases
        PsiMod.set_local_option('SCF','GUESS','READ')
        PsiMod.set_local_option('SCF','BASIS',basis)
        if (scf_type == 'DF'):
            PsiMod.set_local_option('SCF','RI_BASIS_SCF',ri_basis_scf)
            PsiMod.set_global_option('RI_INTS_IO',ri_ints)

        # Print the banner for the standard operation
        PsiMod.print_out('\n')
        banner(name.upper())
        PsiMod.print_out('\n')

        # Do the full scf
        e_scf = PsiMod.scf(precallback, postcallback)

        PsiMod.set_local_option('SCF','GUESS',guess_type)

    else:

        e_scf = PsiMod.scf(precallback, postcallback)

    return e_scf

def run_mp2(name, **kwargs):

    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')

    if not molecule:
        raise ValueNotSet("no molecule found")

    molecule.update_geometry()
    PsiMod.set_active_molecule(molecule)

    run_scf("scf", **kwargs);

    PsiMod.transqt2()
    PsiMod.ccsort()
    return PsiMod.mp2()

def run_ccsd(name, **kwargs):

    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')

    if not molecule:
        raise ValueNotSet("no molecule found")

    molecule.update_geometry()
    PsiMod.set_active_molecule(molecule)

    # For a CCSD energy, we need SCF to be run.
    # Could we somehow do a check to see if SCF was run?
    # This would be useful of the user had to do something special with SCF to get
    # it to converge.
    run_scf("scf", **kwargs);

    PsiMod.transqt2()
    PsiMod.ccsort()
    return PsiMod.ccenergy()

def run_ccsd_gradient(name, **kwargs):
    run_ccsd(name, **kwargs)
    PsiMod.cchbar()
    PsiMod.cclambda()
    PsiMod.ccdensity()
    PsiMod.deriv()

def run_ccsd_t(name, **kwargs):

    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')

    if not molecule:
        raise ValueNotSet("no molecule found")

    molecule.update_geometry()
    PsiMod.set_active_molecule(molecule)

    # Should probably do a check on the user's options to ensure wfn = ccsd_t

    # The new CCEnergyWavefunction object that is used to wrap ccenergy
    # automatically handles cctriples.
    return run_ccsd(name, **kwargs)

def run_ccsd_response(name, **kwargs):

    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')

    if not molecule:
        raise ValueNotSet("no molecule found")

    molecule.update_geometry()
    PsiMod.set_active_molecule(molecule)

    run_scf("scf", **kwargs);

    PsiMod.transqt2()
    PsiMod.ccsort()
    PsiMod.ccenergy()
    PsiMod.cchbar()
    PsiMod.cclambda()
    # Need oeprop integrals here -- code will fail at present
    PsiMod.ccresponse()

def run_detci(name, **kwargs):

    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')

    if not molecule:
        raise ValueNotSet("no molecule found")

    molecule.update_geometry()
    PsiMod.set_active_molecule(molecule)

    # For a CCSD energy, we need SCF to be run.
    # Could we somehow do a check to see if SCF was run?
    # This would be useful of the user had to do something special with SCF to get
    # it to converge.
    run_scf("scf", **kwargs);
    PsiMod.transqt()
    PsiMod.detci()

def run_dfmp2(name, **kwargs):

    if kwargs.has_key('restart_file'):
        restartfile = kwargs.pop('restart_file')
        # Rename the checkpoint file to be consistent with psi4's file system
        psioh      = PsiMod.IOManager.shared_object()
        psio       = PsiMod.IO.shared_object()
        filepath   = psioh.get_file_path(32)
        namespace  = psio.get_default_namespace()
        pid        = str(os.getpid())
        prefix     = 'psi'
        targetfile =  filepath + prefix + '.' + pid + '.' + namespace + ".32"
        if(PsiMod.me() == 0):
            shutil.copy(restartfile, targetfile)
    else:
        run_scf('RHF',**kwargs)

    PsiMod.print_out("\n")
    banner("DFMP2")
    PsiMod.print_out("\n")
    e_dfmp2 = PsiMod.dfmp2()
    e_scs_dfmp2 = PsiMod.get_variable('SCS-DF-MP2 ENERGY')
    if (name.upper() == 'SCS-DFMP2'):
        return e_scs_dfmp2
    elif (name.upper() == 'DFMP2'):
        return e_dfmp2

def run_mp2drpa(name, **kwargs):

    e_scf = run_scf('RHF',**kwargs)

    PsiMod.print_out("\n")
    banner("DFMP2")
    PsiMod.print_out("\n")
    e_dfmp2 = PsiMod.dfmp2()

    PsiMod.print_out("\n")
    banner("dRPA")
    PsiMod.print_out("\n")
    PsiMod.dfcc()
    e_delta = PsiMod.get_variable('RPA SCALED DELTA ENERGY')

    PsiMod.print_out("\n")
    banner("MP2/dRPA Analysis")
    PsiMod.print_out("\n")

    PsiMod.print_out('  Reference Energy               %20.14f\n' % (e_scf))
    PsiMod.print_out('  MP2 Correlation Energy         %20.14f\n' % (e_dfmp2 - e_scf))
    PsiMod.print_out('  MP2 Total Energy               %20.14f\n' % (e_dfmp2))
    PsiMod.print_out('  Scaled dRPA/MP2J Delta Energy  %20.14f\n' % (e_delta))
    PsiMod.print_out('  Total MP2/dRPA Energy          %20.14f\n\n' % (e_dfmp2 + e_delta))

    return e_dfmp2 + e_delta

def run_dfcc(name, **kwargs):

    run_scf('RHF',**kwargs)

    PsiMod.print_out("\n")
    banner("DF-CC")
    PsiMod.print_out("\n")
    e_dfcc = PsiMod.dfcc()
    return e_dfcc

def run_mp2c(name, **kwargs):

    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')

    if not molecule:
        raise ValueNotSet("no molecule found")

    molecule.update_geometry()
    monomerA = molecule.extract_subsets(1,2)
    monomerA.set_name("monomerA")
    monomerB = molecule.extract_subsets(2,1)
    monomerB.set_name("monomerB")

    ri = PsiMod.get_option('SCF_TYPE')
    ri_ints_io = PsiMod.get_option('RI_INTS_IO')

    PsiMod.IO.set_default_namespace("dimer")
    PsiMod.set_local_option("SCF","SAPT","2-dimer")
    PsiMod.print_out("\n")
    banner('Dimer HF')
    PsiMod.print_out("\n")
    PsiMod.set_global_option('RI_INTS_IO','SAVE')
    e_dimer = scf_helper('RHF',**kwargs)
    PsiMod.print_out("\n")
    banner('Dimer DFMP2')
    PsiMod.print_out("\n")
    e_dimer_mp2 = PsiMod.dfmp2()
    PsiMod.set_global_option('RI_INTS_IO','LOAD')

    activate(monomerA)
    if (ri == "DF"):
        PsiMod.IO.change_file_namespace(97,"dimer","monomerA")
    PsiMod.IO.set_default_namespace("monomerA")
    PsiMod.set_local_option("SCF","SAPT","2-monomer_A")
    PsiMod.print_out("\n")
    banner('Monomer A HF')
    PsiMod.print_out("\n")
    e_monomerA = scf_helper('RHF',**kwargs)
    PsiMod.print_out("\n")
    banner('Monomer A DFMP2')
    PsiMod.print_out("\n")
    e_monomerA_mp2 = PsiMod.dfmp2()

    activate(monomerB)
    if (ri == "DF"):
        PsiMod.IO.change_file_namespace(97,"monomerA","monomerB")
    PsiMod.IO.set_default_namespace("monomerB")
    PsiMod.set_local_option("SCF","SAPT","2-monomer_B")
    PsiMod.print_out("\n")
    banner('Monomer B HF')
    PsiMod.print_out("\n")
    e_monomerB = scf_helper('RHF',**kwargs)
    PsiMod.print_out("\n")
    banner('Monomer B DFMP2')
    PsiMod.print_out("\n")
    e_monomerB_mp2 = PsiMod.dfmp2()
    PsiMod.set_global_option('RI_INTS_IO',ri_ints_io)

    PsiMod.IO.change_file_namespace(121,"monomerA","dimer")
    PsiMod.IO.change_file_namespace(122,"monomerB","dimer")

    activate(molecule)
    PsiMod.IO.set_default_namespace("dimer")
    PsiMod.set_local_option("SAPT","E_CONVERGE",10)
    PsiMod.set_local_option("SAPT","D_CONVERGE",10)
    PsiMod.set_local_option("SAPT","SAPT_LEVEL","MP2C")
    PsiMod.print_out("\n")
    banner("MP2C")
    PsiMod.print_out("\n")

    PsiMod.set_variable("MP2C DIMER MP2 ENERGY",e_dimer_mp2)
    PsiMod.set_variable("MP2C MONOMER A MP2 ENERGY",e_monomerA_mp2)
    PsiMod.set_variable("MP2C MONOMER B MP2 ENERGY",e_monomerB_mp2)

    e_sapt = PsiMod.sapt()
    return e_sapt

def run_sapt(name, **kwargs):

    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')

    if not molecule:
        raise ValueNotSet("no molecule found")

    molecule.update_geometry()
    monomerA = molecule.extract_subsets(1,2)
    monomerA.set_name("monomerA")
    monomerB = molecule.extract_subsets(2,1)
    monomerB.set_name("monomerB")

    ri = PsiMod.get_option('SCF_TYPE')
    ri_ints_io = PsiMod.get_option('RI_INTS_IO')

    PsiMod.IO.set_default_namespace("dimer")
    PsiMod.set_local_option("SCF","SAPT","2-dimer")
    PsiMod.print_out("\n")
    banner('Dimer HF')
    PsiMod.print_out("\n")
    PsiMod.set_global_option('RI_INTS_IO','SAVE')
    e_dimer = scf_helper('RHF',**kwargs)
    PsiMod.set_global_option('RI_INTS_IO','LOAD')

    activate(monomerA)
    if (ri == "DF"):
        PsiMod.IO.change_file_namespace(97,"dimer","monomerA")
    PsiMod.IO.set_default_namespace("monomerA")
    PsiMod.set_local_option("SCF","SAPT","2-monomer_A")
    PsiMod.print_out("\n")
    banner('Monomer A HF')
    PsiMod.print_out("\n")
    e_monomerA = scf_helper('RHF',**kwargs)

    activate(monomerB)
    if (ri == "DF"):
        PsiMod.IO.change_file_namespace(97,"monomerA","monomerB")
    PsiMod.IO.set_default_namespace("monomerB")
    PsiMod.set_local_option("SCF","SAPT","2-monomer_B")
    PsiMod.print_out("\n")
    banner('Monomer B HF')
    PsiMod.print_out("\n")
    e_monomerB = scf_helper('RHF',**kwargs)
    PsiMod.set_global_option('RI_INTS_IO',ri_ints_io)

    PsiMod.IO.change_file_namespace(121,"monomerA","dimer")
    PsiMod.IO.change_file_namespace(122,"monomerB","dimer")

    activate(molecule)
    PsiMod.IO.set_default_namespace("dimer")
    PsiMod.set_local_option("SAPT","E_CONVERGE",10)
    PsiMod.set_local_option("SAPT","D_CONVERGE",10)
    if (name.lower() == 'sapt0'):
        PsiMod.set_local_option("SAPT","SAPT_LEVEL","SAPT0")
    elif (name.lower() == 'sapt2'):
        PsiMod.set_local_option("SAPT","SAPT_LEVEL","SAPT2")
    elif (name.lower() == 'sapt2+'):
        PsiMod.set_local_option("SAPT","SAPT_LEVEL","SAPT2+")
    elif (name.lower() == 'sapt2+3'):
        PsiMod.set_local_option("SAPT","SAPT_LEVEL","SAPT2+3")
    PsiMod.print_out("\n")
    banner(name.upper())
    PsiMod.print_out("\n")
    e_sapt = PsiMod.sapt()
    return e_sapt

