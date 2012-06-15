from __future__ import print_function
"""Module with a *procedures* dictionary specifying available quantum
chemical methods and functions driving the main quantum chemical
functionality, namely single-point energies, geometry optimizations,
properties, and vibrational frequency calculations.

"""
import PsiMod
import input
from proc import *
from text import *
from procutil import *
from functional import *
# never import wrappers or aliases into this file


# Procedure lookup tables
procedures = {
        'energy': {
            'scf'           : run_scf,
            'mcscf'         : run_mcscf,
            'dcft'          : run_dcft,
            'dfmp2'         : run_dfmp2,
            'df-mp2'        : run_dfmp2,
            'mp2'           : run_mp2,
            'omp2'          : run_omp2,
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
            'eom_ccsd'      : run_eom_cc,
            'eom_cc2'       : run_eom_cc,
            'eom_cc3'       : run_eom_cc,
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
            'tddft'         : run_libfock,
            'psimrcc'       : run_psimrcc,
            'psimrcc_scf'   : run_psimrcc_scf,
            'hf'            : run_scf,
            'rhf'           : run_scf,
            'uhf'           : run_scf,
            'rohf'          : run_scf,
            'rscf'          : run_scf,
            'uscf'          : run_scf,
            'roscf'         : run_scf,
            'df-scf'        : run_scf,
            'cepa(0)'       : run_cepa,
            'cepa(1)'       : run_cepa,
            'cepa(3)'       : run_cepa,
            'acpf'          : run_cepa,
            'aqcc'          : run_cepa,
            'sdci'          : run_cepa,
            'dci'           : run_cepa,
            'b2plyp'        : run_b2plyp,
            'pbe0-2'        : run_pbe0_2,
            # Upon adding a method to this list, add it to the docstring in energy() below
        },
        'gradient' : {
            'scf'           : run_scf_gradient,
            'ccsd'          : run_cc_gradient,
            'ccsd(t)'       : run_cc_gradient,
            'mp2'           : run_mp2_gradient,
            'df-mp2'        : run_dfmp2_gradient,
            'dfmp2'         : run_dfmp2_gradient,
            'eom-ccsd'      : run_eom_cc_gradient,
            'dcft'          : run_dcft_gradient
            # Upon adding a method to this list, add it to the docstring in optimize() below
        },
        'hessian' : {
            # Upon adding a method to this list, add it to the docstring in frequency() below
        },
        'property' : {
            'scf'  : run_scf_property,
            'cc2'  : run_cc_property,
            'ccsd' : run_cc_property,
            'eom-cc2'  : run_cc_property,
            'eom-ccsd' : run_cc_property,
            'eom_cc2'  : run_cc_property,
            'eom_ccsd' : run_cc_property
            # Upon adding a method to this list, add it to the docstring in property() below
        }}

# Integrate DFT with driver routines
for ssuper in superfunctional_list():
    procedures['energy'][ssuper.name().lower()] = run_dft

for ssuper in superfunctional_list():
    procedures['gradient'][ssuper.name().lower()] = run_dft_gradient


def energy(name, **kwargs):
    r"""Function to compute the single-point electronic energy.

    :returns: (*float*) Total electronic energy in Hartrees. SAPT returns interaction energy.

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`CURRENT ENERGY <CURRENTENERGY>`
       * :psivar:`CURRENT REFERENCE ENERGY <CURRENTREFERENCEENERGY>`
       * :psivar:`CURRENT CORRELATION ENERGY <CURRENTCORRELATIONENERGY>`

    .. comment In this table immediately below, place methods that should only be called by 
    .. comment developers at present. This table won't show up in the manual.
    .. comment
    .. comment    .. _`table:energy_devel`:
    .. comment 
    .. comment    +-------------------------+---------------------------------------------------------------------------------------+
    .. comment    | name                    | calls method                                                                          |
    .. comment    +=========================+=======================================================================================+
    .. comment    | df-cc                   | coupled cluster with density fitting                                                  |
    .. comment    +-------------------------+---------------------------------------------------------------------------------------+
    .. comment    | mp2c                    | coupled MP2 (MP2C)                                                                    |
    .. comment    +-------------------------+---------------------------------------------------------------------------------------+
    .. comment    | mp2-drpa                | random phase approximation?                                                           |
    .. comment    +-------------------------+---------------------------------------------------------------------------------------+
    .. comment    | cphf                    | coupled-perturbed Hartree-Fock?                                                       |
    .. comment    +-------------------------+---------------------------------------------------------------------------------------+
    .. comment    | cpks                    | coupled-perturbed Kohn-Sham?                                                          |
    .. comment    +-------------------------+---------------------------------------------------------------------------------------+
    .. comment    | cis                     | CI singles (CIS)                                                                      |
    .. comment    +-------------------------+---------------------------------------------------------------------------------------+
    .. comment    | tda                     | Tamm-Dankoff approximation (TDA)                                                      |
    .. comment    +-------------------------+---------------------------------------------------------------------------------------+
    .. comment    | tdhf                    | time-dependent HF (TDHF)                                                              |
    .. comment    +-------------------------+---------------------------------------------------------------------------------------+
    .. comment    | tddft                   | time-dependent DFT (TDDFT)                                                            |
    .. comment    +-------------------------+---------------------------------------------------------------------------------------+

    .. _`table:energy_gen`:

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
    | mp\ *n*                 | *n*\ th-order Moller--Plesset perturbation theory                                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | zapt\ *n*               | *n*\ th-order z-averaged perturbation theory (ZAPT)                                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cisd                    | configuration interaction (CI) singles and doubles (CISD)                             |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cisdt                   | CI singles, doubles, and triples (CISDT)                                              |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cisdtq                  | CI singles, doubles, triples, and quadruples (CISDTQ)                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ci\ *n*                 | *n*\ th-order CI                                                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | fci                     | full configuration interaction (FCI)                                                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | detci                   | **expert** full control over detci module                                             |
    +-------------------------+---------------------------------------------------------------------------------------+
    | adc                     | 2nd-order algebraic diagrammatic construction (ADC)                                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-cc2                 | EOM-CC2                                                                               |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-ccsd                | equation of motion (EOM) CCSD                                                         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-cc3                 | EOM-CC3                                                                               |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cepa(n)                 | coupled electron pair approximation, variants 0, 1, and 3                             |
    +-------------------------+---------------------------------------------------------------------------------------+
    | acpf                    | averaged coupled-pair functional                                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | aqcc                    | averaged quadratic coupled cluster                                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | omp2                    | orbital-optimized second order Moller--Plesset perturbation theory                    |
    +-------------------------+---------------------------------------------------------------------------------------+

    .. _`table:energy_scf`:

    +-------------------------+---------------------------------------------------------------------------------------+
    | name                    | calls method (aliases to *name* = 'scf')                                              |
    +=========================+=======================================================================================+
    | hf                      | HF                                                                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | rhf                     | HF with restricted reference                                                          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | uhf                     | HF with unrestricted reference                                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | rohf                    | HF with restricted open-shell reference                                               |
    +-------------------------+---------------------------------------------------------------------------------------+
    | rscf                    | HF or DFT with restricted reference                                                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | uscf                    | HF or DFT with unrestricted reference                                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | roscf                   | HF or DFT with restricted open-shell reference                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-scf                  | HF or DFT with density fitting                                                        |
    +-------------------------+---------------------------------------------------------------------------------------+

    .. include:: autodoc_dft_energy.rst

    .. _`table:energy_mrcc`:

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
    | mrccsdt-1a              | CC through doubles with iterative triples (cheapest terms)                            |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq-1a             | CC through triples with iterative quadruples (cheapest terms)                         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp-1a            | CC through quadruples with iterative quintuples (cheapest terms)                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqph-1a           | CC through quintuples with iterative sextuples (cheapest terms)                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt-1b              | CC through doubles with iterative triples (cheaper terms)                             |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq-1b             | CC through triples with iterative quadruples (cheaper terms)                          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp-1b            | CC through quadruples with iterative quintuples (cheaper terms)                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqph-1b           | CC through quintuples with iterative sextuples (cheaper terms)                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc2                   | approximate CC through doubles                                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc3                   | approximate CC through triples                                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc4                   | approximate CC through quadruples                                                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc5                   | approximate CC through quintuples                                                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc6                   | approximate CC through sextuples                                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt-3               | CC through doubles with iterative triples (all but the most expensive terms)          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq-3              | CC through triples with iterative quadruples (all but the most expensive terms)       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp-3             | CC through quadruples with iterative quintuples (all but the most expensive terms)    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqph-3            | CC through quintuples with iterative sextuples (all but the most expensive terms)     |
    +-------------------------+---------------------------------------------------------------------------------------+

    :type name: string
    :param name: ``'scf'`` || ``'df-mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type bypass_scf: :ref:`boolean <op_py_boolean>`
    :param bypass_scf: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether, for *name* values built atop of scf calculations,
        the scf step is skipped. Suitable when special steps are taken to get
        the scf to converge in an explicit preceeding scf step.

    :examples:

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
    if 'molecule' in kwargs:
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()

    # Allow specification of methods to arbitrary order
    lowername, level = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    try:
        return procedures['energy'][lowername](lowername, **kwargs)
    except KeyError:
        raise ValidationError('Energy method %s not available.' % (lowername))


def gradient(name, **kwargs):
    r"""Function complementary to optimize(). Carries out one gradient pass,
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
    if lowername in procedures['gradient']:
        dertype = 1
    elif lowername in procedures['energy']:
        dertype = 0
        func = energy

    if (PsiMod.get_global_option('REFERENCE').lower() == 'rks') or (PsiMod.get_global_option('REFERENCE').lower() == 'uks'):
        dertype = 0
        func = energy

    # 2. Check if the user passes dertype into this function
    if 'dertype' in kwargs:
        opt_dertype = kwargs['dertype']

        if input.der0th.match(str(opt_dertype)):
            dertype = 0
            func = energy
        elif input.der1st.match(str(opt_dertype)):
            dertype = 1
        else:
            raise ValidationError('Requested derivative level \'dertype\' %s not valid for helper function optimize.' % (opt_dertype))

    # 3. if the user provides a custom function THAT takes precendence
    if ('opt_func' in kwargs) or ('func' in kwargs):
        if ('func' in kwargs):
            kwargs['opt_func'] = kwargs['func']
            del kwargs['func']
        dertype = 0
        func = kwargs['opt_func']

    # Summary validation
    if (dertype == 1) and (lowername in procedures['gradient']):
        pass
    elif (dertype == 0) and (func is energy) and (lowername in procedures['energy']):
        pass
    elif (dertype == 0) and not(func is energy):
        pass
    else:
        raise ValidationError('Requested method \'name\' %s and derivative level \'dertype\' %s are not available.'
            % (lowername, dertype))

    # Make sure the molecule the user provided is the active one
    if ('molecule' in kwargs):
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    PsiMod.set_global_option('BASIS', PsiMod.get_global_option('BASIS'))

    # S/R: Mode of operation- whether finite difference opt run in one job or files farmed out
    opt_mode = 'continuous'
    if ('mode' in kwargs) and (dertype == 0):
        opt_mode = kwargs['mode']

    if (opt_mode.lower() == 'continuous'):
        pass
    elif (opt_mode.lower() == 'sow'):
        pass
    elif (opt_mode.lower() == 'reap'):
        if('linkage' in kwargs):
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

        if 'mode' in kwargs and kwargs['mode'].lower() == 'sow':
            raise ValidationError('Optimize execution mode \'sow\' not valid for analytic gradient calculation.')
        PsiMod.reference_wavefunction().energy()
        return PsiMod.get_variable('CURRENT ENERGY')
    else:
        # If not, perform finite difference of energies

        opt_iter = 1
        if ('opt_iter' in kwargs):
            opt_iter = kwargs['opt_iter']

        if opt_iter == 1:
            print('Performing finite difference calculations')

        # Obtain list of displacements
        displacements = PsiMod.fd_geoms_1_0()
        ndisp = len(displacements)

        # This version is pretty dependent on the reference geometry being last (as it is now)
        print(' %d displacements needed ...' % (ndisp), end="")
        energies = []

        # S/R: Write instructions for sow/reap procedure to output file and reap input file
        if (opt_mode.lower() == 'sow'):
            instructionsO = """\n    The optimization sow/reap procedure has been selected through mode='sow'. In addition\n"""
            instructionsO += """    to this output file (which contains no quantum chemical calculations), this job\n"""
            instructionsO += """    has produced a number of input files (OPT-%s-*.in) for individual components\n""" % (str(opt_iter))
            instructionsO += """    and a single input file (OPT-master.in) with an optimize(mode='reap') command.\n"""
            instructionsO += """    These files may look very peculiar since they contain processed and pickled python\n"""
            instructionsO += """    rather than normal input. Follow the instructions in OPT-master.in to continue.\n\n"""
            instructionsO += """    Alternatively, a single-job execution of the gradient may be accessed through\n"""
            instructionsO += """    the optimization wrapper option mode='continuous'.\n\n"""
            PsiMod.print_out(instructionsO)

            instructionsM = """\n#    Follow the instructions below to carry out this optimization cycle.\n#\n"""
            instructionsM += """#    (1)  Run all of the OPT-%s-*.in input files on any variety of computer architecture.\n""" % (str(opt_iter))
            instructionsM += """#       The output file names must be as given below.\n#\n"""
            for rgt in range(ndisp):
                pre = 'OPT-' + str(opt_iter) + '-' + str(rgt + 1)
                instructionsM += """#             psi4 -i %-27s -o %-27s\n""" % (pre + '.in', pre + '.out')
            instructionsM += """#\n#    (2)  Gather all the resulting output files in a directory. Place input file\n"""
            instructionsM += """#         OPT-master.in into that directory and run it. The job will be minimal in\n"""
            instructionsM += """#         length and give summary results for the gradient step in its output file.\n#\n"""
            if opt_iter == 1:
                instructionsM += """#             psi4 -i %-27s -o %-27s\n#\n""" % ('OPT-master.in', 'OPT-master.out')
            else:
                instructionsM += """#             psi4 -a -i %-27s -o %-27s\n#\n""" % ('OPT-master.in', 'OPT-master.out')
            instructionsM += """#    After each optimization iteration, the OPT-master.in file is overwritten so return here\n"""
            instructionsM += """#    for new instructions. With the use of the psi4 -a flag, OPT-master.out is not\n"""
            instructionsM += """#    overwritten and so maintains a history of the job. To use the (binary) optimizer\n"""
            instructionsM += """#    data file to accelerate convergence, the OPT-master jobs must run on the same computer.\n\n"""

            fmaster = open('OPT-master.in', 'w')
            fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n')
            fmaster.write(format_molecule_for_input(molecule))
            fmaster.write(format_options_for_input())
            format_kwargs_for_input(fmaster, 2, **kwargs)
            fmaster.write("""%s('%s', **kwargs)\n\n""" % (optimize.__name__, lowername))
            fmaster.write(instructionsM)
            fmaster.close()

        for n, displacement in enumerate(displacements):
            rfile = 'OPT-%s-%s' % (opt_iter, n + 1)
            #rfile = 'OPT-fd-%s' % (n + 1)

            # Build string of title banner
            banners = ''
            banners += """PsiMod.print_out('\\n')\n"""
            banners += """banner(' Gradient %d Computation: Displacement %d')\n""" % (opt_iter, n + 1)
            banners += """PsiMod.print_out('\\n')\n\n"""

            if (opt_mode.lower() == 'continuous'):
                # Print information to output.dat
                PsiMod.print_out('\n')
                banner('Loading displacement %d of %d' % (n + 1, ndisp))

                # Print information to the screen
                print(' %d' % (n + 1), end="")
                if (n + 1) == ndisp:
                    print('\n', end="")

                # Load in displacement into the active molecule
                PsiMod.get_active_molecule().set_geometry(displacement)

                # Perform the energy calculation
                #E = func(lowername, **kwargs)
                func(lowername, **kwargs)
                E = PsiMod.get_variable('CURRENT ENERGY')
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
                freagent.write("""PsiMod.print_out('\\nGRADIENT RESULT: computation %d for item %d """ % (os.getpid(), n + 1))
                freagent.write("""yields electronic energy %20.12f\\n' % (electronic_energy))\n\n""")
                freagent.close()

            # S/R: Read energy from each displaced geometry output file and save in energies array
            elif (opt_mode.lower() == 'reap'):
                E = 0.0
                exec(banners)

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
                            if s[6] != str(n + 1):
                                raise ValidationError('Output file \'%s.out\' has nominal affiliation %s incompatible with item %s.'
                                    % (rfile, s[6], str(n + 1)))
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


def property(name, **kwargs):
    r"""Function to compute various properties.

    :aliases: prop()

    :returns: none.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - This function at present handles property functions only for CC methods.
         Consult the keywords sections for other modules for further property capabilities.

    +-------------------------+---------------------------------------------------------------------------------------+
    | name                    | calls method                                                                          |
    +=========================+=======================================================================================+
    | scf                     | Self-consistent field method(s)                                                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cc2                     | 2nd-order approximate CCSD                                                            |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD)                                            |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-cc2                 | 2nd-order approximate EOM-CCSD                                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-ccsd                | equation-of-motion coupled cluster singles and doubles (EOM-CCSD)                     |
    +-------------------------+---------------------------------------------------------------------------------------+

    :type name: string
    :param name: ``'ccsd'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type properties: array of strings
    :param properties: |dl| ``[]`` |dr| || ``['rotation', 'polarizability', 'oscillator_strength', 'roa']`` || etc.

        Indicates which properties should be computed.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :examples:

    >>> # [1] Optical rotation calculation
    >>> property('cc2', properties=['rotation'])

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # Make sure the molecule the user provided is the active one
    if ('molecule' in kwargs):
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    #PsiMod.set_global_option('BASIS', PsiMod.get_global_option('BASIS'))

    # Allow specification of methods to arbitrary order
    lowername, level = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    try:
        return procedures['property'][lowername](lowername, **kwargs)
    except KeyError:
        raise ValidationError('Property method %s not available.' % (lowername))

##  Aliases  ##
prop = property


def optimize(name, **kwargs):
    r"""Function to perform a geometry optimization.

    :aliases: opt()

    :returns: (*float*) Total electronic energy of optimized structure in Hartrees.

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`CURRENT ENERGY <CURRENTENERGY>`

    .. note:: Analytic gradients area available for all methods in the table
        below. Optimizations with other methods in the energy table proceed
        by finite differences.

    .. _`table:grad_gen`:

    +-------------------------+---------------------------------------------------------------------------------------+
    | name                    | calls method                                                                          |
    +=========================+=======================================================================================+
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT)                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | dcft                    | density cumulant functional theory                                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp2                     | 2nd-order Moller-Plesset perturbation theory (MP2)                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD)                                            |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd(t)                 | CCSD with perturbative triples                                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-ccsd                | equation of motion (EOM) CCSD                                                         |
    +-------------------------+---------------------------------------------------------------------------------------+

    :type name: string
    :param name: ``'scf'`` || ``'df-mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the database. May be any valid argument to
        :py:func:`driver.energy`.

    :type func: :ref:`function <op_py_function>`
    :param func: |dl| ``gradient`` |dr| || ``energy`` || ``cbs``

        Indicates the type of calculation to be performed on the molecule.
        The default dertype accesses``'gradient'`` or ``'energy'``, while
        ``'cbs'`` performs a multistage finite difference calculation.
        If a nested series of python functions is intended (see :ref:`sec:intercalls`),
        use keyword ``opt_func`` instead of ``func``.

    :type mode: string
    :param mode: |dl| ``'continuous'`` |dr| || ``'sow'`` || ``'reap'``

        For a finite difference of energies optimization, indicates whether
        the calculations required to complete the
        optimization are to be run in one file (``'continuous'``) or are to be
        farmed out in an embarrassingly parallel fashion
        (``'sow'``/``'reap'``).  For the latter, run an initial job with
        ``'sow'`` and follow instructions in its output file.

    :type dertype: :ref:`dertype <op_py_dertype>`
    :param dertype: ``'gradient'`` || ``'energy'``

        Indicates whether analytic (if available) or finite difference
        optimization is to be performed.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :examples:

    >>> # [1] Analytic scf optimization
    >>> optimize('scf')

    >>> # [2] Finite difference mp3 optimization
    >>> opt('mp3')

    >>> # [3] Forced finite difference ccsd optimization
    >>> optimize('ccsd', dertype=1)

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    full_hess_every = PsiMod.get_local_option('OPTKING', 'FULL_HESS_EVERY')
    steps_since_last_hessian = 0

    n = 1
    if ('opt_iter' in kwargs):
        n = kwargs['opt_iter']

    while n <= PsiMod.get_option('GEOM_MAXITER'):
        kwargs['opt_iter'] = n

        # Compute the gradient
        thisenergy = gradient(name, **kwargs)

        # S/R: Quit after getting new displacements or if forming gradient fails
        if ('mode' in kwargs) and (kwargs['mode'].lower() == 'sow'):
            return 0.0
        if ('mode' in kwargs) and (kwargs['mode'].lower() == 'reap') and (thisenergy == 0.0):
            return 0.0

        # S/R: Move opt data file from last pass into namespace for this pass
        if ('mode' in kwargs) and (kwargs['mode'].lower() == 'reap') and (n != 0):
            PsiMod.IOManager.shared_object().set_specific_retention(1, True)
            PsiMod.IOManager.shared_object().set_specific_path(1, './')
            if 'opt_datafile' in kwargs:
                restartfile = kwargs.pop('opt_datafile')
                if(PsiMod.me() == 0):
                    shutil.copy(restartfile, get_psifile(1))

        # print 'full_hess_every', full_hess_every
        # print 'steps_since_last_hessian', steps_since_last_hessian
        # compute Hessian as requested; frequency wipes out gradient so stash it
        if ((full_hess_every > -1) and (n == 1)) or (steps_since_last_hessian + 1 == full_hess_every):
            G = PsiMod.get_gradient()
            PsiMod.IOManager.shared_object().set_specific_retention(1, True)
            PsiMod.IOManager.shared_object().set_specific_path(1, './')
            frequencies(name, **kwargs)
            steps_since_last_hessian = 0
            PsiMod.set_gradient(G)
            PsiMod.set_global_option('CART_HESS_READ', True)
        elif ((full_hess_every == -1) and (PsiMod.get_global_option('CART_HESS_READ')) and (n == 1)):
            pass;
            # Do nothing; user said to read existing hessian once
        else:
            PsiMod.set_global_option('CART_HESS_READ', False)
            steps_since_last_hessian += 1

        # print 'cart_hess_read', PsiMod.get_global_option('CART_HESS_READ')
        # Take step
        if PsiMod.optking() == PsiMod.PsiReturnType.EndLoop:
            print('Optimizer: Optimization complete!')
            PsiMod.get_active_molecule().print_in_input_format()
            # Check if user wants to see the intcos; if so, don't delete them.
            if (PsiMod.get_option('INTCOS_GENERATE_EXIT') == False):
                PsiMod.opt_clean()
            PsiMod.clean()

            # S/R: Clean up opt input file
            if ('mode' in kwargs) and (kwargs['mode'].lower() == 'reap'):
                fmaster = open('OPT-master.in', 'w')
                fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n')
                fmaster.write('# Optimization complete!\n\n')
                fmaster.close()
            return thisenergy

        # S/R: Preserve opt data file for next pass and switch modes to get new displacements
        if ('mode' in kwargs) and (kwargs['mode'].lower() == 'reap'):
            kwargs['opt_datafile'] = get_psifile(1)
            kwargs['mode'] = 'sow'

        n += 1

    PsiMod.print_out('\tOptimizer: Did not converge!')
    return 0.0

##  Aliases  ##
opt = optimize


def parse_arbitrary_order(name):
    r"""Function to parse name string into a method family like CI or MRCC and specific
    level information like 4 for CISDTQ or MRCCSDTQ.

    """
    namelower = name.lower()

    # matches 'mrccsdt(q)'
    if namelower.startswith('mrcc'):
        # grabs 'sdt(q)'
        ccfullname = namelower[4:]

        # A negative order indicates perturbative method
        methods = {
            'sd'          : { 'method' : 1, 'order' :  2, 'fullname' : 'CCSD'         },
            'sdt'         : { 'method' : 1, 'order' :  3, 'fullname' : 'CCSDT'        },
            'sdtq'        : { 'method' : 1, 'order' :  4, 'fullname' : 'CCSDTQ'       },
            'sdtqp'       : { 'method' : 1, 'order' :  5, 'fullname' : 'CCSDTQP'      },
            'sdtqph'      : { 'method' : 1, 'order' :  6, 'fullname' : 'CCSDTQPH'     },
            'sd(t)'       : { 'method' : 3, 'order' : -3, 'fullname' : 'CCSD(T)'      },
            'sdt(q)'      : { 'method' : 3, 'order' : -4, 'fullname' : 'CCSDT(Q)'     },
            'sdtq(p)'     : { 'method' : 3, 'order' : -5, 'fullname' : 'CCSDTQ(P)'    },
            'sdtqp(h)'    : { 'method' : 3, 'order' : -6, 'fullname' : 'CCSDTQP(H)'   },
            'sd(t)_l'     : { 'method' : 4, 'order' : -3, 'fullname' : 'CCSD(T)_L'    },
            'sdt(q)_l'    : { 'method' : 4, 'order' : -4, 'fullname' : 'CCSDT(Q)_L'   },
            'sdtq(p)_l'   : { 'method' : 4, 'order' : -5, 'fullname' : 'CCSDTQ(P)_L'  },
            'sdtqp(h)_l'  : { 'method' : 4, 'order' : -6, 'fullname' : 'CCSDTQP(H)_L' },
            'sdt-1a'      : { 'method' : 5, 'order' :  3, 'fullname' : 'CCSDT-1a'     },
            'sdtq-1a'     : { 'method' : 5, 'order' :  4, 'fullname' : 'CCSDTQ-1a'    },
            'sdtqp-1a'    : { 'method' : 5, 'order' :  5, 'fullname' : 'CCSDTQP-1a'   },
            'sdtqph-1a'   : { 'method' : 5, 'order' :  6, 'fullname' : 'CCSDTQPH-1a'  },
            'sdt-1b'      : { 'method' : 6, 'order' :  3, 'fullname' : 'CCSDT-1b'     },
            'sdtq-1b'     : { 'method' : 6, 'order' :  4, 'fullname' : 'CCSDTQ-1b'    },
            'sdtqp-1b'    : { 'method' : 6, 'order' :  5, 'fullname' : 'CCSDTQP-1b'   },
            'sdtqph-1b'   : { 'method' : 6, 'order' :  6, 'fullname' : 'CCSDTQPH-1b'  },
            '2'           : { 'method' : 7, 'order' :  2, 'fullname' : 'CC2'          },
            '3'           : { 'method' : 7, 'order' :  3, 'fullname' : 'CC3'          },
            '4'           : { 'method' : 7, 'order' :  4, 'fullname' : 'CC4'          },
            '5'           : { 'method' : 7, 'order' :  5, 'fullname' : 'CC5'          },
            '6'           : { 'method' : 7, 'order' :  6, 'fullname' : 'CC6'          },
            'sdt-3'       : { 'method' : 8, 'order' :  3, 'fullname' : 'CCSDT-3'      },
            'sdtq-3'      : { 'method' : 8, 'order' :  4, 'fullname' : 'CCSDTQ-3'     },
            'sdtqp-3'     : { 'method' : 8, 'order' :  5, 'fullname' : 'CCSDTQP-3'    },
            'sdtqph-3'    : { 'method' : 8, 'order' :  6, 'fullname' : 'CCSDTQPH-3'   }
        }

        # looks for 'sdt(q)' in dictionary
        if ccfullname in methods:
            return 'mrcc', methods[ccfullname]
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
    r"""Function to compute harmonic vibrational frequencies.

    :aliases: frequencies(), freq()

    :returns: (*float*) Total electronic energy in Hartrees.

    .. note:: Analytic hessians are not available. Frequencies will proceed through
        finite differences according to availability of gradients or energies.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Make frequency look analogous to gradient, especially in matching derivative levels. Make dertype actually a dertype type.

    .. _`table:freq_gen`:

    :type name: string
    :param name: ``'scf'`` || ``'df-mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type dertype: :ref:`dertype <op_py_dertype>`
    :param dertype: |dl| ``'hessian'`` |dr| || ``'gradient'`` || ``'energy'``

        Indicates whether analytic (if available- they're not), finite
        difference of gradients (if available) or finite difference of
        energies is to be performed.

    :type irrep: int or string
    :param irrep: |dl| ``-1`` |dr| || ``1`` || ``'b2'`` || ``'App'`` || etc.

        Indicates which symmetry block (:ref:`Cotton <table:irrepOrdering>` ordering) of vibrational
        frequencies to be computed. ``1``, ``'1'``, or ``'a1'`` represents
        :math:`a_1`, requesting only the totally symmetric modes.
        ``-1`` indicates a full frequency calculation.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :examples:

    >>> # [1] <example description>
    >>> <example python command>

    >>> # [2] Frequency calculation for b2 modes through finite difference of gradients
    >>> frequencies('scf', dertype=1, irrep=4)

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # Make sure the molecule the user provided is the active one
    if ('molecule' in kwargs):
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    PsiMod.set_global_option('BASIS', PsiMod.get_global_option('BASIS'))

    types = ['energy', 'gradient', 'hessian']

    dertype = 2
    if ('dertype' in kwargs):
        dertype = kwargs['dertype']
        if not (lowername in procedures[types[dertype]]):
            print('Frequencies: dertype = %d for frequencies is not available, switching to automatic determination.' % dertype)
            dertype = -1

    if 'irrep' in kwargs:
        irrep = parse_cotton_irreps(kwargs['irrep']) - 1  # externally, A1 irrep is 1, internally 0
    else:
        irrep = -1  # -1 implies do all irreps

    # By default, set func to the energy function
    func = energy
    func_existed = False
    if 'func' in kwargs:
        func = kwargs['func']
        func_existed = True

    if (not('dertype' in kwargs) or dertype == -1):
        if lowername in procedures['hessian']:
            dertype = 2
        elif lowername in procedures['gradient']:
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
        info = 'Performing finite difference by gradient calculations'
        print(info)

        func = procedures['gradient'][lowername]

        # Obtain list of displacements
        displacements = PsiMod.fd_geoms_freq_1(irrep)

        molecule.reinterpret_coordentry(False)
        molecule.fix_orientation(True)
        # Make a note of the undisplaced molecule's symmetry
        PsiMod.set_parent_symmetry(molecule.schoenflies_symbol())

        ndisp = len(displacements)
        print(' %d displacements needed.' % ndisp)

        #print displacements to output.dat
        #for n, displacement in enumerate(displacements):
        #  displacement.print_out();

        gradients = []
        for n, displacement in enumerate(displacements):
            # Print information to output.dat
            PsiMod.print_out('\n')
            banner('Loading displacement %d of %d' % (n + 1, ndisp))

            # Print information to the screen
            print(' %d' % (n + 1), end="")
            if (n + 1) == ndisp:
                print('\n', end="")

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

        print(' Computation complete.')
        
        # Clear the "parent" symmetry now
        PsiMod.set_parent_symmetry("")

        # TODO: These need to be restored to the user specified setting
        PsiMod.get_active_molecule().fix_orientation(False)
        # But not this one, it always goes back to True
        PsiMod.get_active_molecule().reinterpret_coordentry(True)

    else:  # Assume energy points
        # If not, perform finite difference of energies
        info = 'Performing finite difference calculations by energies'
        print(info)

        # Obtain list of displacements
        displacements = PsiMod.fd_geoms_freq_0(irrep)
        molecule.fix_orientation(True)
        molecule.reinterpret_coordentry(False)
        # Make a note of the undisplaced molecule's symmetry
        PsiMod.set_parent_symmetry(molecule.schoenflies_symbol())

        ndisp = len(displacements)

        # This version is pretty dependent on the reference geometry being last (as it is now)
        print(' %d displacements needed.' % ndisp)
        energies = []
        for n, displacement in enumerate(displacements):
            # Print information to output.dat
            PsiMod.print_out('\n')
            banner('Loading displacement %d of %d' % (n + 1, ndisp))

            # Print information to the screen
            print(' %d' % (n + 1), end="")
            if (n + 1) == ndisp:
                print('\n', end='')

            # Load in displacement into the active molecule
            molecule.set_geometry(displacement)
   
            # Perform the energy calculation
            E = func(lowername, **kwargs)

            # Save the energy
            energies.append(E)

            # clean may be necessary when changing irreps of displacements
            PsiMod.clean()

        # Obtain the gradient. This function stores the gradient into the reference wavefunction.
        PsiMod.fd_freq_0(energies, irrep)

        print(' Computation complete.')
        
        # Clear the "parent" symmetry now
        PsiMod.set_parent_symmetry("")

        # TODO: These need to be restored to the user specified setting
        PsiMod.get_active_molecule().fix_orientation(False)
        # But not this one, it always goes back to True
        PsiMod.get_active_molecule().reinterpret_coordentry(True)

        # The last item in the list is the reference energy, return it
        return energies[-1]

##  Aliases  ##
frequencies = frequency
freq = frequency


# hessian to be changed later to compute force constants
def hessian(name, **kwargs):
    r"""Function to compute force constants. Presently identical to frequency()."""
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)
    frequencies(name, **kwargs)

def molden(filename):
    m = PsiMod.MoldenWriter(PsiMod.reference_wavefunction())
    m.write(filename)

def parse_cotton_irreps(irrep):
    r"""Function to return validated Cotton ordering index from string or integer
    irreducible representation *irrep*.

    """
    cotton = {
        'c1': {
            'a': 1,
            '1': 1
        },
        'ci': {
            'ag': 1,
            'au': 2,
            '1':  1,
            '2':  2
        },
        'c2': {
            'a': 1,
            'b': 2,
            '1': 1,
            '2': 2
        },
        'cs': {
            'ap':  1,
            'app': 2,
            '1':   1,
            '2':   2
        },
        'd2': {
            'a':  1,
            'b1': 2,
            'b2': 3,
            'b3': 4,
            '1':  1,
            '2':  2,
            '3':  3,
            '4':  4
        },
        'c2v': {
            'a1': 1,
            'a2': 2,
            'b1': 3,
            'b2': 4,
            '1':  1,
            '2':  2,
            '3':  3,
            '4':  4
        },
        'c2h': {
            'ag': 1,
            'bg': 2,
            'au': 3,
            'bu': 4,
            '1':  1,
            '2':  2,
            '3':  3,
            '4':  4,
        },
        'd2h': {
            'ag':  1,
            'b1g': 2,
            'b2g': 3,
            'b3g': 4,
            'au':  5,
            'b1u': 6,
            'b2u': 7,
            'b3u': 8,
            '1':   1,
            '2':   2,
            '3':   3,
            '4':   4,
            '5':   5,
            '6':   6,
            '7':   7,
            '8':   8
        }
    }

    point_group = PsiMod.get_active_molecule().schoenflies_symbol().lower()
    irreducible_representation = str(irrep).lower()

    try:
        return cotton[point_group][irreducible_representation]
    except KeyError:
        raise ValidationError("Irrep \'%s\' not valid for point group \'%s\'." % (str(irrep), point_group))
