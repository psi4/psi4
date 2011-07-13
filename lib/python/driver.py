import PsiMod
from proc import *
from text import *

#Procedure lookup tables
procedures = {
        'energy' : {
            'scf'           : run_scf,
            'mcscf'         : run_mcscf,
            'dfmp2'         : run_dfmp2,
            'dfcc'          : run_dfcc,
            'mp2-drpa'      : run_mp2drpa,
            'sapt0'         : run_sapt,
            'sapt2'         : run_sapt,
            'sapt2+'        : run_sapt,
            'sapt2+3'       : run_sapt,
            'mp2c'          : run_mp2c,
            'ccsd'          : run_ccsd,
            'ccsd(t)'       : run_ccsd_t,
            'detci'         : run_detci
        },
        'gradient' : {
            'scf'           : run_scf_gradient
        },
        'hessian' : {
        },
        'response' : {
            'scf' : run_scf,
            'ccsd' : run_ccsd_response
        }}

def energy(name, **kwargs):

    if (kwargs.has_key('molecule')):
        activate(kwargs['molecule'])
    if (kwargs.has_key('bases')):
        pass
    if (kwargs.has_key('functional')):
        pass

    try:
        return procedures['energy'][name](name,**kwargs)
    except KeyError:
        raise SystemExit('Energy Method %s Not Defined' %(name))

def gradient(name, **kwargs):
    if (kwargs.has_key('molecule')):
        activate(kwargs['molecule'])
    if (kwargs.has_key('bases')):
        pass
    if (kwargs.has_key('functional')):
        pass

    dertype = 1
    if (kwargs.has_key('dertype')):
        dertype = kwargs['dertype']

    # By default, set func to the energy function
    func = energy
    func_existed = False
    if (kwargs.has_key('func')):
        func = kwargs['func']
        func_existed = True

    # Does an analytic procedure exist for the requested method?
    if (procedures['gradient'].has_key(name) and dertype == 1 and func_existed == False):
        # Nothing to it but to do it. Gradient information is saved
        # into the current reference wavefunction
        procedures['gradient'][name](name, **kwargs)

        return PsiMod.reference_wavefunction().energy()
    else:
        # If not, perform finite difference of energies
        info = "Performing finite difference calculations"
        print info

        # Obtain the active molecule and update it.
        molecule = PsiMod.get_active_molecule()
        if not molecule:
            raise ValueNotSet("no molecule found")
        molecule.update_geometry()

        # Obtain list of displacements
        displacements = PsiMod.fd_geoms_1_0()

        ndisp = len(displacements)

        # This version is pretty dependent on the reference geometry being last (as it is now)
        print " %d displacments needed." % ndisp
        energies = []
        for n, displacment in enumerate(displacements):
            # Print information to output.dat
            PsiMod.print_out("\n")
            banner("Loading displacement %d of %d" % (n+1, ndisp))

            # Print information to the screen
            print "    displacment %d" % (n+1)

            # Load in displacement into the active molecule
            PsiMod.get_active_molecule().set_geometry(displacment)

            # Perform the energy calculation
            E = func(name, **kwargs)

            # Save the energy
            energies.append(E)

        # Obtain the gradient. This function stores the gradient into the reference wavefunction.
        PsiMod.fd_grad_1_0(energies)

        # The last item in the list is the reference energy, return it
        return energies[-1]

def hessian(name, **kwargs):
    try:
        return procedures['hessian'][name](name, **kwargs)
    except KeyError:
        raise SystemExit('Hessian Method %s Not Defined' %(name))

def response(name, **kwargs):
    try:
        return procedures['response'][name](name, **kwargs)
    except KeyError:
        raise SystemExit('Response Method %s Not Defined' %(name))

def optimize(name, **kwargs):
    for n in range(PsiMod.get_option("GEOM_MAXITER")):
        # Compute the gradient
        thisenergy = gradient(name, **kwargs)

        # Take step
        if PsiMod.optking() == PsiMod.PsiReturnType.EndLoop:
            print "Optimizer: Optimization complete!"
            PsiMod.opt_clean()
            return thisenergy

    PsiMod.print_out("\tOptimizer: Did not converge!")
    return 0.0
