import PsiMod
from driver import *
from molecule import *
from text import *

def run_scf(name, **kwargs):

    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')

    if not molecule:
        raise ValueNotSet("no molecule found")

    molecule.update_geometry()
    PsiMod.set_default_options_for_module("SCF")
    return scf_helper(name, **kwargs)

# SCF helper chooses whether to cast up or just run SCF
# with a standard guess. This preserves previous SCF options
# set by other procedures (eg. SAPT output file types for SCF)

def run_mcscf(**kwargs):
    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
      molecule = kwargs.pop('molecule')

    if not molecule:
      raise ValueNotSet("no molecule found")

    molecule.update_geometry()
    PsiMod.set_default_options_for_module("MCSCF")
    return PsiMod.mcscf()

# SCF helper chooses whether to cast up or just run SCF
# with a standard guess. This preserves previous SCF options
# set by other procedures (eg. SAPT output file types for SCF)


def scf_helper(name, **kwargs):

    cast = False
    if (kwargs.has_key('cast_up')):
        cast = kwargs.pop('cast_up')

    precallback = None
    if (kwargs.has_key('precallback')):
        precallback = kwargs.pop('precallback')

    postcallback = None
    if (kwargs.has_key('postcallback')):
        postcallback = kwargs.pop('postcallback')

    if (cast):

        namespace = PsiMod.IO.get_default_namespace()
        basis = PsiMod.get_option('BASIS')
        ri_basis_scf = PsiMod.get_option('RI_BASIS_SCF')
        scf_type = PsiMod.get_option('SCF_TYPE')
        # Hack to ensure cartesian or pure are used throughout
        puream = PsiMod.get_global_option('PUREAM')
        PsiMod.set_global_option("PUREAM",puream)

        PsiMod.IO.set_default_namespace((namespace + ".guess"))
        PsiMod.set_option("NO_INPUT",True)

        PsiMod.set_option('GUESS','SAD')
        PsiMod.set_option('BASIS','3-21G')
        PsiMod.set_option('RI_BASIS_SCF','cc-pvdz-ri')
        PsiMod.set_option('DUAL_BASIS_SCF',basis)
        PsiMod.set_option('DUAL_BASIS',True)
        PsiMod.set_option('SCF_TYPE','DF')

        PsiMod.print_out('\n')
        banner('Cast-up SCF Computation')
        PsiMod.print_out('\n')

        PsiMod.print_out('\n')
        banner('3-21G/SAD Guess')
        PsiMod.print_out('\n')

        PsiMod.scf()

        PsiMod.IO.change_file_namespace(100,(namespace + ".guess"),namespace)
        PsiMod.IO.set_default_namespace(namespace)

        PsiMod.set_option("NO_INPUT",True)

        PsiMod.set_option('GUESS','DUAL_BASIS')
        PsiMod.set_option('BASIS',basis)
        PsiMod.set_option('RI_BASIS_SCF',ri_basis_scf)
        PsiMod.set_option('DUAL_BASIS_SCF','')
        PsiMod.set_option('DUAL_BASIS',False)
        PsiMod.set_option('SCF_TYPE',scf_type)

        PsiMod.print_out('\n')
        banner(name.upper())
        PsiMod.print_out('\n')

        e_scf = PsiMod.scf(precallback, postcallback)
    else:

        PsiMod.set_option("NO_INPUT",True)
        e_scf = PsiMod.scf(precallback, postcallback)

    return e_scf

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
    PsiMod.ccenergy()

def run_ccsd_t(name, **kwargs):

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
    PsiMod.ccenergy()
    PsiMod.cctriples()

def run_dfmp2(name, **kwargs):

    run_scf('RHF',**kwargs)

    PsiMod.set_default_options_for_module("DFMP2")
    PsiMod.set_option("NO_INPUT",True)

    PsiMod.print_out("\n")
    banner("DFMP2")
    PsiMod.print_out("\n")
    e_dfmp2 = PsiMod.dfmp2()
    e_scs_dfmp2 = PsiMod.get_variable('SCS-DF-MP2 ENERGY')
    if (name.upper() == 'SCS-DFMP2'):
        return e_scs_dfmp2
    elif (name.upper() == 'DFMP2'):
        return e_dfmp2

def run_dfcc(name, **kwargs):

    run_scf('RHF',**kwargs)

    PsiMod.set_default_options_for_module("DFCC")
    PsiMod.set_option("NO_INPUT",True)
    
    PsiMod.print_out("\n")
    banner("DF-CC")
    PsiMod.print_out("\n")
    e_dfcc = PsiMod.dfcc()
    return e_dfcc

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

    PsiMod.IO.set_default_namespace("dimer")
    PsiMod.set_default_options_for_module("SCF")
    PsiMod.set_option("NO_INPUT",True)
    PsiMod.set_option("SAPT","2-dimer")
    PsiMod.print_out("\n")
    banner('Dimer HF')
    PsiMod.print_out("\n")
    e_dimer = scf_helper('RHF',**kwargs)

    activate(monomerA)
    PsiMod.IO.set_default_namespace("monomerA")
    PsiMod.set_default_options_for_module("SCF")
    PsiMod.set_option("NO_INPUT",True)
    PsiMod.set_option("SAPT","2-monomer_A")
    PsiMod.print_out("\n")
    banner('Monomer A HF')
    PsiMod.print_out("\n")
    e_monomerA = scf_helper('RHF',**kwargs)

    activate(monomerB)
    PsiMod.IO.set_default_namespace("monomerB")
    PsiMod.set_default_options_for_module("SCF")
    PsiMod.set_option("NO_INPUT",True)
    PsiMod.set_option("SAPT","2-monomer_B")
    PsiMod.print_out("\n")
    banner('Monomer B HF')
    PsiMod.print_out("\n")
    e_monomerB = scf_helper('RHF',**kwargs)

    PsiMod.IO.change_file_namespace(121,"monomerA","dimer")
    PsiMod.IO.change_file_namespace(122,"monomerB","dimer")

    activate(molecule)
    PsiMod.IO.set_default_namespace("dimer")
    PsiMod.set_default_options_for_module("SAPT")
    PsiMod.set_option("NO_INPUT",True)
    PsiMod.set_option("E_CONVERGE",10)
    PsiMod.set_option("D_CONVERGE",10)
    if (name.lower() == 'sapt0'):
        PsiMod.set_option("SAPT_LEVEL","SAPT0")
    elif (name.lower() == 'scs-sapt'):
        PsiMod.set_option("SAPT_LEVEL","SCS_SAPT")
    elif (name.lower() == 'sapt2'):
        PsiMod.set_option("SAPT_LEVEL","SAPT2")
    elif (name.lower() == 'sapt2+'):
        PsiMod.set_option("SAPT_LEVEL","SAPT2+")
    elif (name.lower() == 'sapt2+3'):
        PsiMod.set_option("SAPT_LEVEL","SAPT2+3")
    elif (name.lower() == 'sapt_dft'):
        PsiMod.set_option("SAPT_LEVEL","SAPT_DFT")
    PsiMod.print_out("\n")
    banner(name.upper())
    PsiMod.print_out("\n")
    e_sapt = PsiMod.sapt()
    return e_sapt

