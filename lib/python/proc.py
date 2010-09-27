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
    PsiMod.set_option("NO_INPUT",True)

    PsiMod.print_out("\n")
    banner(name.upper())
    PsiMod.print_out("\n")
    e_scf = PsiMod.scf()
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

    PsiMod.transqt()
    PsiMod.ccsort()
    PsiMod.ccenergy()

def run_dfmp2(name, **kwargs):

    molecule = PsiMod.get_active_molecule()
    if (kwargs.has_key('molecule')):
        molecule = kwargs.pop('molecule')
    if not molecule:
        raise ValueNotSet("no molecule found")

    molecule.update_geometry()
    PsiMod.set_default_options_for_module("SCF")
    PsiMod.set_option("NO_INPUT",True)

    PsiMod.print_out("\n")
    banner("SCF")
    PsiMod.print_out("\n")
    e_scf = PsiMod.scf()

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
    e_dimer = PsiMod.scf()

    activate(monomerA)
    PsiMod.IO.set_default_namespace("monomerA")
    PsiMod.set_default_options_for_module("SCF")
    PsiMod.set_option("NO_INPUT",True)
    PsiMod.set_option("SAPT","2-monomer_A")
    PsiMod.print_out("\n")
    banner('Monomer A HF')
    PsiMod.print_out("\n")
    e_monomerA = PsiMod.scf()

    activate(monomerB)
    PsiMod.IO.set_default_namespace("monomerB")
    PsiMod.set_default_options_for_module("SCF")
    PsiMod.set_option("NO_INPUT",True)
    PsiMod.set_option("SAPT","2-monomer_B")
    PsiMod.print_out("\n")
    banner('Monomer B HF')
    PsiMod.print_out("\n")
    e_monomerB = PsiMod.scf()

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
    elif (name.lower() == 'sapt2'):
        PsiMod.set_option("SAPT_LEVEL","SAPT2")
    elif (name.lower() == 'sapt2+'):
        PsiMod.set_option("SAPT_LEVEL","SAPT2+")
    elif (name.lower() == 'sapt2+(3)'):
        PsiMod.set_option("SAPT_LEVEL","SAPT2+(3)")
    PsiMod.print_out("\n")
    banner(name.upper())
    PsiMod.print_out("\n")
    e_sapt = PsiMod.sapt()
    return e_sapt

