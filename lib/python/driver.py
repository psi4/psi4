import PsiMod
import input
from proc import *
from text import *
from procutil import *

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
            'ccenergy'      : run_ccenergy,  # full control over ccenergy
            'ccsd'          : run_ccenergy,
            'ccsd(t)'       : run_ccenergy,
            'cc2'           : run_ccenergy,
            'cc3'           : run_ccenergy,
            'mrcc'          : run_mrcc,      # interface to Kallay's MRCC program
            'bccd'          : run_bccd,
            'bccd(t)'       : run_bccd_t,
            'eom-ccsd'      : run_eom_cc,
            'eom-cc2'       : run_eom_cc,
            'eom-cc3'       : run_eom_cc,
            'detci'         : run_detci,  # full control over detci
            'mp'            : run_detci,  # arbitrary order mp(n)
            'zapt'          : run_detci,  # arbitrary order zapt(n)
            'cisd'          : run_detci,
            'cisdt'         : run_detci,
            'cisdtq'        : run_detci,
            'ci'            : run_detci,  # arbitrary order ci(n)
            'fci'           : run_detci,
            'adc'           : run_adc
        },
        'gradient' : {
            'scf'           : run_scf_gradient,
            'ccsd'          : run_cc_gradient,
            'ccsd(t)'       : run_cc_gradient,
            'mp2'           : run_mp2_gradient,
            'eom-ccsd'      : run_eom_cc_gradient
        },
        'hessian' : {
        },
        'response' : {
            'cc2'  : run_cc_response,
            'ccsd' : run_cc_response
        }}

def energy(name, **kwargs):
    lowername = name.lower()

    # Make sure the molecule the user provided is the active one
    if (kwargs.has_key('molecule')):
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()

    # Allow specification of methods to arbitrary order
    lowername, level, uppername = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level
    if uppername:
        kwargs['fullname'] = uppername   # For mrccsd(t) this will be CCSD(T)

    try:
        return procedures['energy'][lowername](lowername,**kwargs)
    except KeyError:
        raise ValidationError('Energy method %s not available.' % (lowername))

def gradient(name, **kwargs):
    lowername = name.lower()
    dertype = 1

    # Order of precedence:
    #    1. Default for wavefunction
    #    2. Value obtained from kwargs, if user changed it
    #    3. If user provides a custom 'func' use that

    # Allow specification of methods to arbitrary order
    lowername, level, uppername = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level
    if uppername:
        kwargs['fullname'] = uppername   # For mrccsd(t) this will be CCSD(T)

    # 1. set the default to that of the provided name
    if (procedures['gradient'].has_key(lowername)):
        dertype = 1
    elif (procedures['energy'].has_key(lowername)):
        dertype = 0
        func = energy

    # 2. Check if the user passes dertype into this function
    if (kwargs.has_key('dertype')):
        opt_dertype = kwargs['dertype']

        if input.der0th.match(str(opt_dertype)):
            dertype = 0
            func = energy
        elif input.der1st.match(str(opt_dertype)):
            dertype = 1
        else:
            raise ValidationError('Requested derivative level \'dertype\' %s not valid for helper function optimize.' % (opt_dertype))

    # 3. if the user provides a custom function THAT takes precendence
    if (kwargs.has_key('opt_func')) or (kwargs.has_key('func')):
        if (kwargs.has_key('func')):
            kwargs['opt_func'] = kwargs['func']
            del kwargs['func']
        dertype = 0
        func = kwargs['opt_func']

    # Summary validation
    if (dertype == 1) and (procedures['gradient'].has_key(lowername)):
        pass
    elif (dertype == 0) and (func is energy) and (procedures['energy'].has_key(lowername)):
        pass
    elif (dertype == 0) and not(func is energy):
        pass
    else:
        raise ValidationError('Requested method \'name\' %s and derivative level \'dertype\' %s are not available.'
            % (lowername, dertype))

    # Make sure the molecule the user provided is the active one
    if (kwargs.has_key('molecule')):
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    PsiMod.set_global_option("BASIS", PsiMod.get_global_option("BASIS"))

    # S/R: Mode of operation- whether finite difference opt run in one job or files farmed out
    opt_mode = 'continuous'
    if (kwargs.has_key('mode')) and (dertype == 0):
        opt_mode = kwargs['mode']

    if re.match(r'^continuous$', opt_mode.lower()):
        pass
    elif re.match(r'^sow$', opt_mode.lower()):
        pass
    elif re.match(r'^reap$', opt_mode.lower()):
        if(kwargs.has_key('linkage')):
            opt_linkage = kwargs['linkage']
        else:
            raise ValidationError('Optimize execution mode \'reap\' requires a linkage option.')
    else:
        raise ValidationError('Optimize execution mode \'%s\' not valid.' % (opt_mode))

    # Does dertype indicate an analytic procedure both exists and is wanted?
    if (dertype == 1):
        # Nothing to it but to do it. Gradient information is saved
        # into the current reference wavefunction
        procedures['gradient'][lowername](lowername, **kwargs)

        return PsiMod.reference_wavefunction().energy()
    else:
        # If not, perform finite difference of energies
        info = "Performing finite difference calculations"
        print info

        # Obtain list of displacements
        displacements = PsiMod.fd_geoms_1_0()
        ndisp = len(displacements)

        # This version is pretty dependent on the reference geometry being last (as it is now)
        print " %d displacements needed." % ndisp
        energies = []

        opt_iter = 1
        if (kwargs.has_key('opt_iter')):
            opt_iter = kwargs['opt_iter'] + 1

        # S/R: Write instructions for sow/reap procedure to output file and reap input file
        if re.match('sow', opt_mode.lower()):
            instructionsO  =   """\n    The optimization sow/reap procedure has been selected through mode='sow'. In addition\n"""
            instructionsO +=     """    to this output file (which contains no quantum chemical calculations), this job\n"""
            instructionsO +=     """    has produced a number of input files (OPT-%s-*.in) for individual components\n""" % (str(opt_iter))
            instructionsO +=     """    and a single input file (OPT-master.in) with an optimize(mode='reap') command.\n"""
            instructionsO +=     """    These files may look very peculiar since they contain processed and pickled python\n"""
            instructionsO +=     """    rather than normal input. Follow the instructions in OPT-master.in to continue.\n\n"""
            instructionsO +=     """    Alternatively, a single-job execution of the gradient may be accessed through\n"""
            instructionsO +=     """    the optimization wrapper option mode='continuous'.\n\n"""
            PsiMod.print_out(instructionsO)

            instructionsM  =   """\n#    Follow the instructions below to carry out this optimization cycle.\n#\n"""
            instructionsM +=     """#    (1)  Run all of the OPT-%s-*.in input files on any variety of computer architecture.\n""" % (str(opt_iter))
            instructionsM +=     """#       The output file names must be as given below.\n#\n"""
            for rgt in range(ndisp):
                pre = 'OPT-' + str(opt_iter) + '-' + str(rgt+1)
                instructionsM += """#             psi4 -i %-27s -o %-27s\n""" % (pre + '.in', pre + '.out')
            instructionsM +=  """#\n#    (2)  Gather all the resulting output files in a directory. Place input file\n"""
            instructionsM +=     """#         OPT-master.in into that directory and run it. The job will be minimal in\n"""
            instructionsM +=     """#         length and give summary results for the gradient step in its output file.\n#\n"""
            if opt_iter == 1:
                instructionsM += """#             psi4 -i %-27s -o %-27s\n#\n""" % ('OPT-master.in', 'OPT-master.out')
            else:
                instructionsM += """#             psi4 -a -i %-27s -o %-27s\n#\n""" % ('OPT-master.in', 'OPT-master.out')
            instructionsM +=     """#    After each optimization iteration, the OPT-master.in file is overwritten so return here\n"""
            instructionsM +=     """#    for new instructions. With the use of the psi4 -a flag, OPT-master.out is not\n"""
            instructionsM +=     """#    overwritten and so maintains a history of the job. To use the (binary) optimizer\n"""
            instructionsM +=     """#    data file to accelerate convergence, the OPT-master jobs must run on the same computer.\n\n"""

            fmaster = open('OPT-master.in', 'w')
            fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n')
            fmaster.write(format_molecule_for_input(molecule))
            fmaster.write(format_options_for_input())
            format_kwargs_for_input(fmaster, 2, **kwargs)
            fmaster.write("""%s('%s', **kwargs)\n\n""" % (optimize.__name__, lowername))
            fmaster.write(instructionsM)
            fmaster.close()

        for n, displacement in enumerate(displacements):
            rfile = 'OPT-%s-%s' % (opt_iter, n+1)

            # Build string of title banner
            banners = ''
            banners += """PsiMod.print_out('\\n')\n"""
            banners += """banner(' Gradient %d Computation: Displacement %d')\n""" % (opt_iter, n+1)
            banners += """PsiMod.print_out('\\n')\n\n"""

            if re.match('continuous', opt_mode.lower()):
                # Print information to output.dat
                PsiMod.print_out("\n")
                banner("Loading displacement %d of %d" % (n+1, ndisp))

                # Print information to the screen
                print "    displacement %d" % (n+1)

                # Load in displacement into the active molecule
                PsiMod.get_active_molecule().set_geometry(displacement)

                ## Wrap any positional arguments into kwargs (for intercalls among wrappers)
                #if not('name' in kwargs) and name:
                #    kwargs['name'] = lowername

                # Perform the energy calculation
                #E = func(lowername, **kwargs)
                func(lowername, **kwargs)
                E = PsiMod.get_variable("CURRENT ENERGY")
                #E = func(**kwargs)

                # Save the energy
                energies.append(E)

            # S/R: Write each displaced geometry to an input file
            elif re.match('sow', opt_mode.lower()):
                PsiMod.get_active_molecule().set_geometry(displacement)

                # S/R: Prepare molecule, options, and kwargs
                freagent = open('%s.in' % (rfile), 'w')
                freagent.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n')
                freagent.write(format_molecule_for_input(molecule))
                freagent.write(format_options_for_input())
                format_kwargs_for_input(freagent, **kwargs)

                # S/R: Prepare function call and energy save
                freagent.write("""electronic_energy = %s('%s', **kwargs)\n\n""" % (func.__name__, lowername))
                freagent.write("""PsiMod.print_out('\\nGRADIENT RESULT: computation %d for item %d """ % (os.getpid(), n+1))
                freagent.write("""yields electronic energy %20.12f\\n' % (electronic_energy))\n\n""")
                freagent.close()

            # S/R: Read energy from each displaced geometry output file and save in energies array
            elif re.match('reap', opt_mode.lower()):
                E = 0.0
                exec banners

                try:
                    freagent = open('%s.out' % (rfile), 'r')
                except IOError:
                    ValidationError('Aborting upon output file \'%s.out\' not found.\n' % (rfile))
                    return 0.0
                else:
                    while 1:
                        line = freagent.readline()
                        if not line:
                            if E == 0.0:
                               ValidationError('Aborting upon output file \'%s.out\' has no %s RESULT line.\n' % (rfile, 'GRADIENT'))
                            break
                        s = line.split()
                        if (len(s) != 0) and (s[0:3] == ['GRADIENT', 'RESULT:', 'computation']):
                            if int(s[3]) != opt_linkage:
                               raise ValidationError('Output file \'%s.out\' has linkage %s incompatible with master.in linkage %s.'
                                   % (rfile, str(s[3]), str(opt_linkage)))
                            if s[6] != str(n+1):
                               raise ValidationError('Output file \'%s.out\' has nominal affiliation %s incompatible with item %s.'
                                   % (rfile, s[6], str(n+1)))
                            if (s[8:10] == ['electronic', 'energy']):
                                E = float(s[10])
                                PsiMod.print_out('%s RESULT: electronic energy = %20.12f\n' % ('GRADIENT', E))
                    freagent.close()
                energies.append(E)

        # S/R: Quit sow after writing files
        if re.match('sow', opt_mode.lower()):
            return 0.0

        if re.match('reap', opt_mode.lower()):
            PsiMod.set_variable('CURRENT ENERGY', energies[-1])

        # Obtain the gradient
        PsiMod.fd_1_0(energies)

        # The last item in the list is the reference energy, return it
        return energies[-1]

def response(name, **kwargs):
    lowername = name.lower()

    # Make sure the molecule the user provided is the active one
    if (kwargs.has_key('molecule')):
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    PsiMod.set_global_option("BASIS", PsiMod.get_global_option("BASIS"))

    try:
        return procedures['response'][lowername](lowername, **kwargs)
    except KeyError:
        raise ValidationError('Response method %s not available.' %(lowername))

def optimize(name, **kwargs):
    n = 0
    if (kwargs.has_key('opt_iter')):
        n = kwargs['opt_iter']

    while n < PsiMod.get_option('GEOM_MAXITER'):
        kwargs['opt_iter'] = n

        # Compute the gradient
        thisenergy = gradient(name, **kwargs)

        # S/R: Quit after getting new displacements or if forming gradient fails
        if (kwargs.has_key('mode')) and (kwargs['mode'].lower() == 'sow'):
            return 0.0
        if (kwargs.has_key('mode')) and (kwargs['mode'].lower() == 'reap') and (thisenergy == 0.0):
            return 0.0

        # S/R: Move opt data file from last pass into namespace for this pass
        if (kwargs.has_key('mode')) and (kwargs['mode'].lower() == 'reap') and (n != 0):
            PsiMod.IOManager.shared_object().set_specific_retention(1, True)
            PsiMod.IOManager.shared_object().set_specific_path(1, './')
            if kwargs.has_key('opt_datafile'):
                restartfile = kwargs.pop('opt_datafile')
                if(PsiMod.me() == 0):
                    shutil.copy(restartfile, get_psifile(1))

        # Take step
        if PsiMod.optking() == PsiMod.PsiReturnType.EndLoop:
            print "Optimizer: Optimization complete!"
            PsiMod.opt_clean()
            PsiMod.clean()

            # S/R: Clean up opt input file
            if (kwargs.has_key('mode')) and (kwargs['mode'].lower() == 'reap'):
                fmaster = open('OPT-master.in', 'w')
                fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n')
                fmaster.write('# Optimization complete!\n\n')
                fmaster.close()
            return thisenergy

        # S/R: Preserve opt data file for next pass and switch modes to get new displacements
        if (kwargs.has_key('mode')) and (kwargs['mode'].lower() == 'reap'):
            kwargs['opt_datafile'] = get_psifile(1)
            kwargs['mode'] = 'sow'

        n += 1

    PsiMod.print_out("\tOptimizer: Did not converge!")
    return 0.0

##  Aliases  ##
opt = optimize

def parse_arbitrary_order(name):
    namelower = name.lower()

    if re.match(r'^mrcc(\w+)(\(\w\))?', namelower):
        ccorder = re.match(r'^mrcc(\w+)(\(\w\))?', namelower)

        orders = [ "sd", "sdt", "sdtq", "sdtqp", "stdqph" ]
        fullcc = ccorder.group(1)

        fullname = "CC" + fullcc
        parenfound = False
        if ccorder.group(2):
            parenfound = True
            parencc = ccorder.group(2)
            parenorder = parencc[1]
            fullcc = fullcc + parenorder
            fullname = fullname + parencc

        try:
            ind = orders.index(fullcc) + 2 # to account for SD being first
            if parenfound == True:
                ind = -ind

            return "mrcc", ind, fullname.upper()

        except ValueError:
            raise ValidationError('MRCC method \'%s\' invalid.' % (namelower))

    elif re.match(r'^[a-z]+\d+$', namelower):
        decompose = re.compile(r'^([a-z]+)(\d+)$').match(namelower)
        namestump = decompose.group(1)
        namelevel = int(decompose.group(2))

        if (namestump == 'mp') or (namestump == 'zapt') or (namestump == 'ci'):
            # Let 'mp2' pass through as itself
            if (namestump == 'mp') and (namelevel == 2):
                return namelower, None, None
            # Otherwise return method and order
            else:
                return namestump, namelevel, None
        else:
            return namelower, None, None
    else:
        return namelower, None, None

def frequencies(name, **kwargs):
    lowername = name.lower()

    # Make sure the molecule the user provided is the active one
    if (kwargs.has_key('molecule')):
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    PsiMod.set_global_option("BASIS", PsiMod.get_global_option("BASIS"))

    dertype = 2
    if (kwargs.has_key('dertype')):
        dertype = kwargs['dertype']

    if (kwargs.has_key('irrep')):
        irrep = kwargs['irrep'] - 1 # externally, A1 irrep is 1; internally 0
    else:
      irrep = -1; # -1 implies do all irreps

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

        func = procedures['gradient'][lowername]

        # Obtain list of displacements
        displacements = PsiMod.fd_geoms_freq_1(irrep)

        PsiMod.get_active_molecule().fix_orientation(True)
        PsiMod.get_active_molecule().fix_com(True)
        PsiMod.get_active_molecule().reinterpret_coordentry(False)

        ndisp = len(displacements)

        #print displacements to output.dat
        #for n, displacement in enumerate(displacements):
        #  displacement.print_out();

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
            func(lowername, **kwargs)

            # Save the gradient
            G = PsiMod.get_gradient()
            gradients.append(G)

            # clean may be necessary when changing irreps of displacements
            PsiMod.clean()

        PsiMod.fd_freq_1(gradients, irrep)

        print " Computation complete."

        # TODO: These need to be restored to the user specified setting
        PsiMod.get_active_molecule().fix_orientation(False)
        PsiMod.get_active_molecule().fix_com(False)
        # But not this one, it always goes back to True
        PsiMod.get_active_molecule().reinterpret_coordentry(True)

    else: # Assume energy points
        # If not, perform finite difference of energies
        info = "Performing finite difference calculations"
        print info

        # Obtain list of displacements
        displacements = PsiMod.fd_geoms_freq_0(irrep)
        PsiMod.get_active_molecule().fix_orientation(True)
        PsiMod.get_active_molecule().fix_com(True)
        PsiMod.get_active_molecule().reinterpret_coordentry(False)

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

            # clean may be necessary when changing irreps of displacements
            PsiMod.clean()

        # Obtain the gradient. This function stores the gradient into the reference wavefunction.
        PsiMod.fd_freq_0(energies, irrep)

        print " Computation complete."

        # TODO: These need to be restored to the user specified setting
        PsiMod.get_active_molecule().fix_orientation(False)
        PsiMod.get_active_molecule().fix_com(False)
        # But not this one, it always goes back to True
        PsiMod.get_active_molecule().reinterpret_coordentry(True)

        # The last item in the list is the reference energy, return it
        return energies[-1]

# hessian to be changed later to compute force constants
def hessian(name, **kwargs):
    frequencies(name, **kwargs)
