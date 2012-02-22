"""Module with a *procedures* dictionary specifying available quantum
chemical methods and functions driving the main quantum chemical
functionality, namely single-point energies, geometry optimizations,
response properties, and vibrational frequency calculations.

"""

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
            'sapt2+(3)'     : run_sapt,
            'sapt2+3'       : run_sapt,
            'sapt0-ct'      : run_sapt_ct,
            'sapt2-ct'      : run_sapt_ct,
            'sapt2+-ct'     : run_sapt_ct,
            'sapt2+(3)-ct'  : run_sapt_ct,
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
            'adc'           : run_adc,
            'cphf'          : run_libfock,
            'cis'           : run_libfock,
            'tdhf'          : run_libfock,
            'cpks'          : run_libfock,
            'tda'           : run_libfock,
            'tddft'         : run_libfock
            # Upon adding a method to this list, add it to the docstring in energy() below
        },
        'gradient' : {
            'scf'           : run_scf_gradient,
            'ccsd'          : run_cc_gradient,
            'ccsd(t)'       : run_cc_gradient,
            'mp2'           : run_mp2_gradient,
            'eom-ccsd'      : run_eom_cc_gradient
            # Upon adding a method to this list, add it to the docstring in optimize() below
        },
        'hessian' : {
            # Upon adding a method to this list, add it to the docstring in frequency() below
        },
        'response' : {
            'cc2'  : run_cc_response,
            'ccsd' : run_cc_response
            # Upon adding a method to this list, add it to the docstring in response() below
        }}

def energy(name, **kwargs):
    """Function to compute the single-point electronic energy.

    :returns: (*float*) Total electronic energy in Hartrees. SAPT returns interaction energy.

    :PSI variables:

    .. envvar:: CURRENT ENERGY
        CURRENT REFERENCE ENERGY
        CURRENT CORRELATION ENERGY

    +-------------------------+---------------------------------------------------------------------------------------+
    | name                    | calls method                                                                          | 
    +=========================+=======================================================================================+
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT)                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp2                     | 2nd-order Moller-Plesset perturbation theory (MP2)                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-mp2                  | MP2 with density fitting                                                              |
    +-------------------------+---------------------------------------------------------------------------------------+
    | dcft                    | density cumulant functional theory                                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mcscf                   | multiconfigurational self consistent field (SCF)                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | dfcc                    | coupled cluster with density fitting                                                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp2c                    | coupled MP2 (MP2C)                                                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp2-drpa                | random phase approximation?                                                           |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt0                   | 0th-order symmetry adapted perturbation theory (SAPT)                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2                   | 2nd-order SAPT, traditional definition                                                |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+                  | SAPT including all 2nd-order terms                                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+(3)               | SAPT including perturbative triples                                                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+3                 |                                                                                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt0-ct                | 0th-order SAPT plus charge transfer (CT) calculation                                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2-ct                | SAPT2 plus CT                                                                         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+-ct               | SAPT2+ plus CT                                                                        |   
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+(3)-ct            | SAPT2+(3) plus CT                                                                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+3-ct              | SAPT2+3 plus CT                                                                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cc2                     | approximate coupled cluster singles and doubles (CC2)                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD)                                            |
    +-------------------------+---------------------------------------------------------------------------------------+
    | bccd                    | Brueckner coupled cluster doubles (BCCD)                                              |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cc3                     | approximate coupled cluster singles, doubles, and triples (CC3)                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd(t)                 | CCSD with perturbative triples                                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | bccd(t)                 | BCCD with perturbative triples                                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccenergy                | **expert** full control over ccenergy module                                          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp *n*                  | *n* th-order Moller--Plesset perturbation theory                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | zapt *n*                | *n* th-order z-averaged perturbation theory (ZAPT)                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cisd                    | configuration interaction (CI) singles and doubles (CISD)                             |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cisdt                   | CI singles, doubles, and triples (CISDT)                                              |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cisdtq                  | CI singles, doubles, triples, and quadruples (CISDTQ)                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ci *n*                  | *n* th-order CI                                                                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | fci                     | full configuration interaction (FCI)                                                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | detci                   | **expert** full control over detci module                                             |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cphf                    | coupled-perturbed Hartree-Fock?                                                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cpks                    | coupled-perturbed Kohn-Sham?                                                          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cis                     | CI singles (CIS)                                                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | tda                     | Tamm-Dankoff approximation (TDA)                                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | tdhf                    | time-dependent HF (TDHF)                                                              |
    +-------------------------+---------------------------------------------------------------------------------------+
    | tddft                   | time-dependent DFT (TDDFT)                                                            |
    +-------------------------+---------------------------------------------------------------------------------------+
    | adc                     | 2nd-order algebraic diagrammatic construction (ADC)                                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-cc2                 | EOM-CC2                                                                               | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-ccsd                | equation of motion (EOM) CCSD                                                         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-cc3                 | EOM-CC3                                                                               |
    +-------------------------+---------------------------------------------------------------------------------------+


    +-------------------------+---------------------------------------------------------------------------------------+
    | name                    | calls method in Kallay's MRCC program                                                 | 
    +=========================+=======================================================================================+
    | mrccsd                  | CC through doubles                                                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt                 | CC through triples                                                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq                | CC through quadruples                                                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp               | CC through quintuples                                                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqph              | CC through sextuples                                                                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsd(t)               | CC through doubles with perturbative triples                                          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt(q)              | CC through triples with perturbative quadruples                                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq(p)             | CC through quadruples with pertubative quintuples                                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp(h)            | CC through quintuples with pertubative sextuples                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsd(t)_l             |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt(q)_l            |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq(p)_l           |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp(h)_l          |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt-1a              |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq-1a             |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp-1a            |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqph-1a           |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt-1b              |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq-1b             |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp-1b            |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqph-1b           |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc2                   |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc3                   |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc4                   |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc5                   |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc6                   |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt-3               |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq-3              |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp-3             |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqph-3            |                                                                                       | 
    +-------------------------+---------------------------------------------------------------------------------------+

    **Keywords**

    :type name: string
    :param name: ``'scf'`` || ``'df-mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method 
        to be applied to the system.

    :type bypass_scf: bool
    :param bypass_scf: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether, for *name* values built atop of scf calculations,
        the scf step is skipped. Suitable when special steps are taken to get
        the scf to converge in an explicit preceeding scf step.

    **Examples**

    >>> # [1] Coupled-cluster singles and doubles calculation with psi code
    >>> energy('ccsd')

    >>> # [2] Charge-transfer SAPT calculation with scf projection from small into requested basis
    >>> energy('sapt0-ct',cast_up=True)

    >>> # [3] Arbitrary-order MPn calculation
    >>> energy('mp4')

    """

    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # Make sure the molecule the user provided is the active one
    if (kwargs.has_key('molecule')):
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()

    # Allow specification of methods to arbitrary order
    lowername, level = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    try:
        return procedures['energy'][lowername](lowername,**kwargs)
    except KeyError:
        raise ValidationError('Energy method %s not available.' % (lowername))

def gradient(name, **kwargs):
    """Function complementary to optimize(). Carries out one gradient pass,
    deciding analytic or finite difference.

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)
    dertype = 1

    # Order of precedence:
    #    1. Default for wavefunction
    #    2. Value obtained from kwargs, if user changed it
    #    3. If user provides a custom 'func' use that

    # Allow specification of methods to arbitrary order
    lowername, level = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

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

    if (opt_mode.lower() == 'continuous'):
        pass
    elif (opt_mode.lower() == 'sow'):
        pass
    elif (opt_mode.lower() == 'reap'):
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

        opt_iter = 1
        if (kwargs.has_key('opt_iter')):
            opt_iter = kwargs['opt_iter'] + 1

        if opt_iter == 1:
            print 'Performing finite difference calculations'

        # Obtain list of displacements
        displacements = PsiMod.fd_geoms_1_0()
        ndisp = len(displacements)

        # This version is pretty dependent on the reference geometry being last (as it is now)
        print " %d displacements needed ..." % (ndisp),
        energies = []

        # S/R: Write instructions for sow/reap procedure to output file and reap input file
        if (opt_mode.lower() == 'sow'):
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

            if (opt_mode.lower() == 'continuous'):
                # Print information to output.dat
                PsiMod.print_out("\n")
                banner("Loading displacement %d of %d" % (n+1, ndisp))

                # Print information to the screen
                print " %d" % (n + 1),
                if (n + 1) == ndisp:
                    print "\n",

                # Load in displacement into the active molecule
                PsiMod.get_active_molecule().set_geometry(displacement)

                # Perform the energy calculation
                #E = func(lowername, **kwargs)
                func(lowername, **kwargs)
                E = PsiMod.get_variable("CURRENT ENERGY")
                #E = func(**kwargs)

                # Save the energy
                energies.append(E)

            # S/R: Write each displaced geometry to an input file
            elif (opt_mode.lower() == 'sow'):
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
            elif (opt_mode.lower() == 'reap'):
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
        if (opt_mode.lower() == 'sow'):
            return 0.0

        if (opt_mode.lower() == 'reap'):
            PsiMod.set_variable('CURRENT ENERGY', energies[-1])

        # Obtain the gradient
        PsiMod.fd_1_0(energies)

        # The last item in the list is the reference energy, return it
        return energies[-1]

def response(name, **kwargs):
    """Function to compute linear response properties.

    :returns: (*float*) Total electronic energy in Hartrees.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Check that energy is actually being returned.

       - Check if ther're some PSI variables that ought to be set.

    +-------------------------+---------------------------------------------------------------------------------------+
    | name                    | calls method                                                                          | 
    +=========================+=======================================================================================+
    | cc2                     | 2nd-order approximate CCSD                                                            |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD)                                            |
    +-------------------------+---------------------------------------------------------------------------------------+

    **Keywords**

    :type name: string
    :param name: ``'ccsd'`` || etc.

        First argument, usually unlabeled. Indicates the computational method 
        to be applied to the system.

    **Examples**

    >>> # [1] CCSD-LR properties calculation
    >>> response('ccsd')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

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
    """Function to perform a geometry optimization.

    :aliases: opt()

    :returns: (*float*) Total electronic energy of optimized structure in Hartrees.

    :PSI variables:
    .. envvar:: CURRENT ENERGY

    .. note:: Analytic gradients area available for all methods in the table
        below. Optimizations with other methods in the energy table proceed
        by finite differences.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Need to check that all methods do return electronic energy. I think gradient got changed at one point.

    +-------------------------+---------------------------------------------------------------------------------------+
    | name                    | calls method                                                                          | 
    +=========================+=======================================================================================+
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT)                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp2                     | 2nd-order Moller-Plesset perturbation theory (MP2)                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD)                                            |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd(t)                 | CCSD with perturbative triples                                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-ccsd                | equation of motion (EOM) CCSD                                                         |
    +-------------------------+---------------------------------------------------------------------------------------+

    **Keywords**

    :type name: string
    :param name: ``'scf'`` || ``'df-mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method 
        to be applied to the database. May be any valid argument to 
        :py:func:`driver.energy`.

    :type func: function
    :param func: |dl| ``gradient`` |dr| || ``energy`` || ``cbs``
    
        Indicates the type of calculation to be performed on the molecule.
        The default dertype accesses``'gradient'`` or ``'energy'``, while
        ``'cbs'`` performs a multistage finite difference calculation.
        If a nested series of python functions is intended (see `Function Intercalls`_),
        use keyword ``opt_func`` instead of ``func``.

    :type mode: string
    :param mode: |dl| ``'continuous'`` |dr| || ``'sow'`` || ``'reap'``

        Indicates whether the calculation required to complete the
        optimization are to be run in one file (``'continuous'``) or are to be
        farmed out in an embarrassingly parallel fashion
        (``'sow'``/``'reap'``).  For the latter, run an initial job with
        ``'sow'`` and follow instructions in its output file.

    :type dertype: dertype
    :param dertype: ``'gradient'`` || ``'energy'``

        Indicates whether analytic (if available) or finite difference
        optimization is to be performed.

    **Examples**

    >>> # [1] Analytic scf optimization
    >>> optimize('scf')

    >>> # [2] Finite difference mp3 optimization
    >>> opt('mp3')

    >>> # [3] Forced finite difference ccsd optimization
    >>> optimize('ccsd', dertype=1)

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    full_hess_every = PsiMod.get_global_option("FULL_HESS_EVERY")
    steps_since_last_hessian = 0

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

        # compute Hessian as requested; frequency wipes out gradient so stash it
        if ((full_hess_every > -1) and (n == 0)) or (steps_since_last_hessian == full_hess_every):
          G = PsiMod.get_gradient()
          PsiMod.IOManager.shared_object().set_specific_retention(1, True)
          PsiMod.IOManager.shared_object().set_specific_path(1, './')
          frequencies(name, **kwargs)
          steps_since_last_hessian = 0
          PsiMod.set_gradient(G)
          PsiMod.set_global_option("CART_HESS_READ", True)
        else:
          PsiMod.set_global_option("CART_HESS_READ", False)

        steps_since_last_hessian += 1

        # Take step
        if PsiMod.optking() == PsiMod.PsiReturnType.EndLoop:
            print "Optimizer: Optimization complete!"
            PsiMod.get_active_molecule().print_in_input_format()
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
    """Function to parse name string into a method family like CI or MRCC and specific
    level information like 4 for CISDTQ or MRCCSDTQ.

    """
    namelower = name.lower()

    # matches 'mrccsdt(q)'
    if namelower.startswith('mrcc'):
        # grabs 'sdt(q)'
        ccfullname = namelower[4:]

        # A negative order indicates perturbative method
        methods = {
            "sd"          : { "method" : 1, "order" :  2, "fullname" : "CCSD"         },
            "sdt"         : { "method" : 1, "order" :  3, "fullname" : "CCSDT"        },
            "sdtq"        : { "method" : 1, "order" :  4, "fullname" : "CCSDTQ"       },
            "sdtqp"       : { "method" : 1, "order" :  5, "fullname" : "CCSDTQP"      },
            "sdtqph"      : { "method" : 1, "order" :  6, "fullname" : "CCSDTQPH"     },
            "sd(t)"       : { "method" : 3, "order" : -3, "fullname" : "CCSD(T)"      },
            "sdt(q)"      : { "method" : 3, "order" : -4, "fullname" : "CCSDT(Q)"     },
            "sdtq(p)"     : { "method" : 3, "order" : -5, "fullname" : "CCSDTQ(P)"    },
            "sdtqp(h)"    : { "method" : 3, "order" : -6, "fullname" : "CCSDTQP(H)"   },
            "sd(t)_l"     : { "method" : 4, "order" : -3, "fullname" : "CCSD(T)_L"    },
            "sdt(q)_l"    : { "method" : 4, "order" : -4, "fullname" : "CCSDT(Q)_L"   },
            "sdtq(p)_l"   : { "method" : 4, "order" : -5, "fullname" : "CCSDTQ(P)_L"  },
            "sdtqp(h)_l"  : { "method" : 4, "order" : -6, "fullname" : "CCSDTQP(H)_L" },
            "sdt-1a"      : { "method" : 5, "order" :  3, "fullname" : "CCSDT-1a"     },
            "sdtq-1a"     : { "method" : 5, "order" :  4, "fullname" : "CCSDTQ-1a"    },
            "sdtqp-1a"    : { "method" : 5, "order" :  5, "fullname" : "CCSDTQP-1a"   },
            "sdtqph-1a"   : { "method" : 5, "order" :  6, "fullname" : "CCSDTQPH-1a"  },
            "sdt-1b"      : { "method" : 6, "order" :  3, "fullname" : "CCSDT-1b"     },
            "sdtq-1b"     : { "method" : 6, "order" :  4, "fullname" : "CCSDTQ-1b"    },
            "sdtqp-1b"    : { "method" : 6, "order" :  5, "fullname" : "CCSDTQP-1b"   },
            "sdtqph-1b"   : { "method" : 6, "order" :  6, "fullname" : "CCSDTQPH-1b"  },
            "2"           : { "method" : 7, "order" :  2, "fullname" : "CC2"          },
            "3"           : { "method" : 7, "order" :  3, "fullname" : "CC3"          },
            "4"           : { "method" : 7, "order" :  4, "fullname" : "CC4"          },
            "5"           : { "method" : 7, "order" :  5, "fullname" : "CC5"          },
            "6"           : { "method" : 7, "order" :  6, "fullname" : "CC6"          },
            "sdt-3"       : { "method" : 8, "order" :  3, "fullname" : "CCSDT-3"      },
            "sdtq-3"      : { "method" : 8, "order" :  4, "fullname" : "CCSDTQ-3"     },
            "sdtqp-3"     : { "method" : 8, "order" :  5, "fullname" : "CCSDTQP-3"    },
            "sdtqph-3"    : { "method" : 8, "order" :  6, "fullname" : "CCSDTQPH-3"   }
        }

        # looks for 'sdt(q)' in dictionary
        if methods.has_key(ccfullname):
            return "mrcc", methods[ccfullname]
        else:
            raise ValidationError('MRCC method \'%s\' invalid.' % (namelower))

    elif re.match(r'^[a-z]+\d+$', namelower):
        decompose = re.compile(r'^([a-z]+)(\d+)$').match(namelower)
        namestump = decompose.group(1)
        namelevel = int(decompose.group(2))

        if (namestump == 'mp') or (namestump == 'zapt') or (namestump == 'ci'):
            # Let 'mp2' pass through as itself
            if (namestump == 'mp') and (namelevel == 2):
                return namelower, None
            # Otherwise return method and order
            else:
                return namestump, namelevel
        else:
            return namelower, None
    else:
        return namelower, None

def frequency(name, **kwargs):
    """Function to compute harmonic vibrational frequencies.

    :aliases: frequencies(), freq()

    :returns: (*float*) Total electronic energy in Hartrees.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - RAK, why are you adding OPTKING options as GLOBALS? And shouldn't they be Py-side not C-side options?

       - Put in a dictionary, so IRREPS can be called by symmetry element or 'all'

       - Make frequency look analogous to gradient, especially in matching derivative levels. Make dertype actually a dertype type.

    **Keywords**

    :type name: string
    :param name: ``'scf'`` || ``'df-mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method 
        to be applied to the system.

    :type dertype: dertype
    :param dertype: |dl| ``'hessian'`` |dr| || ``'gradient'`` || ``'energy'``

        Indicates whether analytic (if available- they're not), finite
        difference of gradients (if available) or finite difference of 
        energies is to be performed.

    :type irrep: int
    :param irrep: |dl| ``-1`` |dr| || ``1`` || etc.

        Indicates which symmetry block of vibrational freqiencies to be
        computed. 1 represents :math:`a_1`, requesting only the totally symmetric modes.
        ``-1`` indicates a full frequency calculation.

    **Examples**

    >>> # [1] <example description>
    >>> <example python command>

    >>> # [2] Frequency calculation for b2 modes through finite difference of gradients
    >>> frequencies('scf', dertype=1, irrep=4)

    """

    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # Make sure the molecule the user provided is the active one
    if (kwargs.has_key('molecule')):
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    PsiMod.set_global_option("BASIS", PsiMod.get_global_option("BASIS"))

    types = [ "energy", "gradient", "hessian" ]

    dertype = 2
    if (kwargs.has_key('dertype')):
        dertype = kwargs['dertype']
        if procedures[types[dertype]].has_key(lowername) == False:
            print "Frequencies: dertype = %d for frequencies is not available, switching to automatic determination." % dertype
            dertype = -1

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

    if (kwargs.has_key('dertype') == False or dertype == -1):
        if procedures['hessian'].has_key(lowername):
            dertype = 2
        elif procedures['gradient'].has_key(lowername):
            dertype = 1
        else:
            dertype = 0

    # Does an analytic procedure exist for the requested method?
    if (dertype == 2 and func_existed == False):
        # We have the desired method. Do it.
        procedures['hessian'][lowername](lowername, **kwargs)
        return PsiMod.reference_wavefunction().energy()
    elif (dertype == 1 and func_existed == False):
        # Ok, we're doing frequencies by gradients
        info = "Performing finite difference by gradient calculations"
        print info

        func = procedures['gradient'][lowername]

        # Obtain list of displacements
        displacements = PsiMod.fd_geoms_freq_1(irrep)

        molecule.reinterpret_coordentry(False)
        molecule.fix_orientation(True)

        ndisp = len(displacements)
        print " %d displacements needed." % ndisp

        #print displacements to output.dat
        #for n, displacement in enumerate(displacements):
        #  displacement.print_out();

        gradients = []
        for n, displacement in enumerate(displacements):
            # Print information to output.dat
            PsiMod.print_out("\n")
            banner("Loading displacement %d of %d" % (n+1, ndisp))

            # Print information to the screen
            print "    displacement %d" % (n+1)

            # Load in displacement into the active molecule (xyz coordinates only)
            molecule.set_geometry(displacement)

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
        # But not this one, it always goes back to True
        PsiMod.get_active_molecule().reinterpret_coordentry(True)

    else: # Assume energy points
        # If not, perform finite difference of energies
        info = "Performing finite difference calculations by energies"
        print info

        # Obtain list of displacements
        displacements = PsiMod.fd_geoms_freq_0(irrep)
        PsiMod.get_active_molecule().fix_orientation(True)
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
        # But not this one, it always goes back to True
        PsiMod.get_active_molecule().reinterpret_coordentry(True)

        # The last item in the list is the reference energy, return it
        return energies[-1]

## Aliases ##
frequencies = frequency
freq = frequency

# hessian to be changed later to compute force constants
def hessian(name, **kwargs):
    """Function to compute force constants. Presently identical to frequency()."""
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)
    frequencies(name, **kwargs)

