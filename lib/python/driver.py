import PsiMod
from proc import *
from text import *

#Procedure lookup tables
procedures = {
        'energy' : {
            'scf'           : run_scf,
            'mcscf'         : run_mcscf,
            'dcft'          : run_dcft,
            'dfmp2'         : run_dfmp2,
            'dfcc'          : run_dfcc,
            'mp2'           : run_mp2,
            'mp2-drpa'      : run_mp2drpa,
            'sapt0'         : run_sapt,
            'sapt2'         : run_sapt,
            'sapt2+'        : run_sapt,
            'sapt2+3'       : run_sapt,
            'sapt0-ct'      : run_sapt_ct,
            'sapt2-ct'      : run_sapt_ct,
            'sapt2+-ct'     : run_sapt_ct,
            'sapt2+3-ct'    : run_sapt_ct,
            'mp2c'          : run_mp2c,
            'ccsd'          : run_ccsd,
            'ccsd(t)'       : run_ccsd_t,
            'detci'         : run_detci
        },
        'gradient' : {
            'scf'           : run_scf_gradient,
            'ccsd'          : run_ccsd_gradient,
            'mp2'           : run_mp2_gradient
        },
        'hessian' : {
        },
        'response' : {
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
        lowername = name.lower()
        return procedures['energy'][lowername](lowername,**kwargs)
    except KeyError:
        raise SystemExit('Energy Method %s Not Defined' %(name))

def gradient(name, **kwargs):
    # Set some defaults
    func = energy
    lowername = name.lower()
    # Order of precedence:
    #    1. Default for wavefunction
    #    2. Value obtained from liboptions, if user changed it
    #    3. If user provides a custom 'func' use that

    # 1. set the default to that of the provided name
    if (procedures['gradient'].has_key(lowername)):
        dertype = 1
    elif (procedures['energy'].has_key(lowername)):
        dertype = 0

    # 2. Check if the user set the global option
    if (PsiMod.has_option_changed('DERTYPE') and dertype != -1):
        option = PsiMod.get_option('DERTYPE')
        if (option == 'NONE'):
            dertype = 0
        elif (option == 'FIRST'):
            dertype = 1
        else:
            raise SystemExit("Value of DERTYPE option is not NONE or FIRST.")

    # 3. user passes dertype into this function
    if (kwargs.has_key('dertype')):
        dertype = kwargs['dertype']

    # 4. if the user provides a custom function THAT takes precendence
    if (kwargs.has_key('func')):
        dertype = 0
        func = kwargs['func']

    # Start handling the other options we support
    if (kwargs.has_key('molecule')):
        # Make sure the molecule the user provided is the active one
        activate(kwargs['molecule'])

    # dertype currently holds the type of calculation we need to run

    # First, check the values of dertype
    if (dertype < 0 and dertype > 1):
        raise SystemExit("The internal dertype is either less than 0 or greater than 1")

    # Does an analytic procedure exist for the requested method?
    if (dertype == 1):
        # Nothing to it but to do it. Gradient information is saved
        # into the current reference wavefunction
        procedures['gradient'][lowername](lowername, **kwargs)

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
        print " %d displacements needed." % ndisp
        energies = []
        for n, displacement in enumerate(displacements):
            # Print information to output.dat
            PsiMod.print_out("\n")
            banner("Loading displacement %d of %d" % (n+1, ndisp))

            # Print information to the screen
            print "    displacement %d" % (n+1)

            # Load in displacement into the active molecule
            PsiMod.get_active_molecule().set_geometry(displacement)

            # Perform the energy calculation
            E = func(lowername, **kwargs)

            # Save the energy
            energies.append(E)

        # Obtain the gradient. This function stores the gradient into the reference wavefunction.
        PsiMod.fd_1_0(energies)

        # The last item in the list is the reference energy, return it
        return energies[-1]

def hessian(name, **kwargs):
    lowername = name.lower()
    if (kwargs.has_key('molecule')):
        activate(kwargs['molecule'])

    dertype = 2
    if (kwargs.has_key('dertype')):
        dertype = kwargs['dertype']

    # By default, set func to the energy function
    func = energy
    func_existed = False
    if (kwargs.has_key('func')):
        func = kwargs['func']
        func_existed = True

    # Does an analytic procedure exist for the requested method?
    if (procedures['hessian'].has_key(lowername) and dertype == 2 and func_existed == False):
        # We have the desired method. Do it.
        procedures['hessian'][lowername](lowername, **kwargs)
        return PsiMod.reference_wavefunction().energy()
    elif (procedures['gradient'].has_key(lowername) and dertype == 1 and func_existed == False):
        # Ok, we're doing frequencies by gradients
        info = "Performing finite difference by gradient calculations"
        print info

        # Obtain the active molecule and update it
        molecule = PsiMod.get_active_molecule()
        if not molecule:
            raise ValueNotSet("no molecule found")
        molecule.update_geometry()

        func = procedures['gradient'][lowername]

        # Obtain list of displacements
        displacements = fd_geoms_freq_1()
        ndisp = len(displacements)

        print " %d displacements needed." % ndisp
        gradients = []
        for n, displacement in enumerate(displacements):
            # Print information to output.dat
            PsiMod.print_out("\n")
            banner("Loading displacement %d of %d" % (n+1, ndisp))

            # Print information to the screen
            print "    displacement %d" % (n+1)

            # Load in displacement into the active molecule
            PsiMod.get_active_molecule().set_geometry(displacement)

            # Perform the gradient calculation
            G = func(lowername, **kwargs)

            # Save the gradient
            gradients.append(G)

        # What I am supposed to do here?

    else: # Assume energy points
        # If not, perform finite difference of energies
        info = "Performing finite difference calculations"
        print info

        # Obtain the active molecule and update it.
        molecule = PsiMod.get_active_molecule()
        if not molecule:
            raise ValueNotSet("no molecule found")
        molecule.update_geometry()

        # Obtain list of displacements
        displacements = PsiMod.fd_geoms_freq_0()

        ndisp = len(displacements)

        # This version is pretty dependent on the reference geometry being last (as it is now)
        print " %d displacments needed." % ndisp
        energies = []
        for n, displacement in enumerate(displacements):
            # Print information to output.dat
            PsiMod.print_out("\n")
            banner("Loading displacement %d of %d" % (n+1, ndisp))

            # Print information to the screen
            print "    displacement %d" % (n+1)

            # Load in displacement into the active molecule
            PsiMod.get_active_molecule().set_geometry(displacement)

            # Perform the energy calculation
            E = func(lowername, **kwargs)

            # Save the energy
            energies.append(E)

        # Obtain the gradient. This function stores the gradient into the reference wavefunction.
        PsiMod.fd_freq_0(energies)

        print " Computation complete."

        # The last item in the list is the reference energy, return it
        return energies[-1]

def response(name, **kwargs):
    lowername = name.lower()
    try:
        return procedures['response'][lowername](lowername, **kwargs)
    except KeyError:
        raise SystemExit('Response Method %s Not Defined' %(lowername))

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

def frequencies(name, **kwargs):
    hessian(name, **kwargs)

