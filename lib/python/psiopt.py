from PsiMod import *

def scf_energy():
    # Would be nice to do a check that required modules were executed.
    scf()

def scf_gradient():
    # First, we'll need to determine if we're doing findif of energies
    # or analytic gradients.

    #if dertype == Energy:
    #    do steps needed for findif of energies
    #    Does this involve calling optking to generate displacements?
    #    or, will "opt(...)" tell us to use energies and provide
    #    displacements

    #elif dertype == Gradient:
    # At this point we will assume that the active molecule is the one we
    # should be using.

    # Compute energy
    scf_energy()

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
