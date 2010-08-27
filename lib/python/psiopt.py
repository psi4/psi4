from PsiMod import *

def ccsd_energy():
    # Do checks to see if scf, transqt, and ccsort executed?

    transqt()
    ccsort()
    ccenergy()

def ccsd_gradient():
    # Do checks to see if ccsd executed?

    cchbar()
    cclambda()
    ccdensity()
    oeprop()
    transqt_to_ao()

    # Using one- and two-particle densities and lagrangian compute gradient
    #   NOTE: This code needs to be modified, currently coded to SCF gradient
    # was cints --deriv1
    deriv()

def scf_energy():
    # Would be nice to do a check that required modules were executed.

    # Need to check what the user asked for:
    #    direct and df (scf will compute integrals needed)
    #    pk and out of core [cints (or mints with modification) needs
    #                        to run first]

    scf()

def scf_gradient():
    # First, we'll need to determine if we're doing findif of energies
    # or analytic gradients.

    #if dertype == Energy:
    #    do steps needed for findif of energies
    #    Does this involve calling optking to generate displacements?
    #    or, will "opt(...)", below, tell us to use energies and provide
    #    displacements

    #elif dertype == Gradient:
    # At this point we will assume that the active molecule is the one we
    # should be using.

    # Compute energy
    scf_energy()       # or should we call scf() directly

    # Compute gradient:
    deriv()

def opt(method = scf_gradient):
    # Assume the active molecule is the one we want to use.

    # Tell the user what we're going to do:
    print_out("Optimizer:")

    # Compute the gradient
    method()

    # Tell the user that we don't know what to do with the gradient
    print_out("Optimizer: Don't know what to do with the gradient!")
