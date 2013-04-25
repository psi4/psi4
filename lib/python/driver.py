from __future__ import print_function
"""Module with a *procedures* dictionary specifying available quantum
chemical methods and functions driving the main quantum chemical
functionality, namely single-point energies, geometry optimizations,
properties, and vibrational frequency calculations.

"""
import sys
import PsiMod
from proc import *
from text import *
from procutil import *
from functional import *
from psifiles import *
# never import wrappers or aliases into this file


# Procedure lookup tables
procedures = {
        'energy': {
            'scf'           : run_scf,
            'mcscf'         : run_mcscf,
            'dcft'          : run_dcft,
            'oldmp2'        : run_oldmp2,
            'dfmp2'         : run_dfmp2,
            'df-mp2'        : run_dfmp2,
            'conv-mp2'      : run_mp2,
            'mp3'           : run_mp3,
            'mp2.5'         : run_mp2_5,
            'mp2'           : run_mp2_select,
            'omp2'          : run_omp2,
            'scs-omp2'      : run_scs_omp2,
            'scsn-omp2'     : run_scs_omp2,
            'scs-mi-omp2'   : run_scs_omp2,
            'scs-omp2-vdw'  : run_scs_omp2,
            'sos-omp2'      : run_sos_omp2,
            'sos-pi-omp2'   : run_sos_omp2,
            'omp3'          : run_omp3,
            'scs-omp3'      : run_scs_omp3,
            'scsn-omp3'     : run_scs_omp3,
            'scs-mi-omp3'   : run_scs_omp3,
            'scs-omp3-vdw'  : run_scs_omp3,
            'sos-omp3'      : run_sos_omp3,
            'sos-pi-omp3'   : run_sos_omp3,
            'ocepa'         : run_ocepa,
            'cepa0'         : run_cepa0,
            'omp2.5'        : run_omp2_5,
            'dftsapt'       : run_dftsapt,
            'saptdft'       : run_dftsapt,
            'infsapt'       : run_infsapt,
            'sapt0'         : run_sapt,
            'sapt2'         : run_sapt,
            'sapt2+'        : run_sapt,
            'sapt2+(3)'     : run_sapt,
            'sapt2+3'       : run_sapt,
            'sapt2+(ccd)'   : run_sapt,
            'sapt2+(3)(ccd)': run_sapt,
            'sapt2+3(ccd)'  : run_sapt,
            'sapt0-ct'      : run_sapt_ct,
            'sapt2-ct'      : run_sapt_ct,
            'sapt2+-ct'     : run_sapt_ct,
            'sapt2+(3)-ct'  : run_sapt_ct,
            'sapt2+3-ct'    : run_sapt_ct,
            'sapt2+(ccd)-ct'     : run_sapt_ct,
            'sapt2+(3)(ccd)-ct'  : run_sapt_ct,
            'sapt2+3(ccd)-ct'    : run_sapt_ct,
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
            'detci-mp'      : run_detci,  # arbitrary order mp(n)
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
            'qcisd'         : run_fnocc,
            'qcisd(t)'      : run_fnocc,
            'mp4(sdq)'      : run_fnocc,
            'fno-ccsd'      : run_fnocc,
            'fno-ccsd(t)'   : run_fnocc,
            'fno-qcisd'     : run_fnocc,
            'fno-qcisd(t)'  : run_fnocc,
            'fno-mp3'       : run_fnocc,
            'fno-mp4(sdq)'  : run_fnocc,
            'fno-mp4'       : run_fnocc,
            'fnocc-mp'      : run_fnocc,
            'df-ccsd'       : run_fnodfcc,
            'df-ccsd(t)'    : run_fnodfcc,
            'fno-df-ccsd'   : run_fnodfcc,
            'fno-df-ccsd(t)': run_fnodfcc,
            'fno-cepa(0)'   : run_cepa,
            'fno-cepa(1)'   : run_cepa,
            'fno-cepa(3)'   : run_cepa,
            'fno-acpf'      : run_cepa,
            'fno-aqcc'      : run_cepa,
            'fno-sdci'      : run_cepa,
            'fno-dci'       : run_cepa,
            'cepa(0)'       : run_cepa,
            'cepa(1)'       : run_cepa,
            'cepa(3)'       : run_cepa,
            'acpf'          : run_cepa,
            'aqcc'          : run_cepa,
            'sdci'          : run_cepa,
            'dci'           : run_cepa,
            # Upon adding a method to this list, add it to the docstring in energy() below
            # If you must add an alias to this list (e.g., dfmp2/df-mp2), please search the
            #    whole driver to find uses of name in return values and psi variables and
            #    extend the logic to encompass the new alias.
        },
        'gradient' : {
            'scf'           : run_scf_gradient,
            'ccsd'          : run_cc_gradient,
            'ccsd(t)'       : run_cc_gradient,
            'mp2'           : run_mp2_select_gradient,
            'conv-mp2'      : run_mp2_gradient,
            'df-mp2'        : run_dfmp2_gradient,
            'dfmp2'         : run_dfmp2_gradient,
            'eom-ccsd'      : run_eom_cc_gradient,
            'dcft'          : run_dcft_gradient,
            'omp2'          : run_omp2_gradient,
            'omp3'          : run_omp3_gradient,
            'mp3'           : run_mp3_gradient,
            'mp2.5'         : run_mp2_5_gradient,
            'omp2.5'        : run_omp2_5_gradient,
            'cepa0'         : run_cepa0_gradient,
            'ocepa'         : run_ocepa_gradient
            # Upon adding a method to this list, add it to the docstring in optimize() below
        },
        'hessian' : {
            # Upon adding a method to this list, add it to the docstring in frequency() below
        },
        'property' : {
            'scf'  : run_scf_property,
            'cc2'  : run_cc_property,
            'ccsd' : run_cc_property,
            'df-mp2' : run_dfmp2_property,
            'dfmp2'  : run_dfmp2_property,
            'eom-cc2'  : run_cc_property,
            'eom-ccsd' : run_cc_property,
            # Upon adding a method to this list, add it to the docstring in property() below
        }}

# Integrate DFT with driver routines
for ssuper in superfunctional_list():
    procedures['energy'][ssuper.name().lower()] = run_dft

for ssuper in superfunctional_list():
    if ((not ssuper.is_c_hybrid()) and (not ssuper.is_c_lrc()) and (not ssuper.is_x_lrc())):
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
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT) :ref:`[manual] <sec:scf>`       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | dcft                    | density cumulant functional theory :ref:`[manual] <sec:dcft>`                         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mcscf                   | multiconfigurational self consistent field (SCF)                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp2                     | 2nd-order Moller-Plesset perturbation theory (MP2) :ref:`[manual] <sec:dfmp2>`        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-mp2                  | MP2 with density fitting :ref:`[manual] <sec:dfmp2>`                                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | conv-mp2                | conventional MP2 (non-density-fitting) :ref:`[manual] <sec:convocc>`                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp3                     | 3rd-order Moller-Plesset perturbation theory (MP3) :ref:`[manual] <sec:convocc>`      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp2.5                   | average of MP2 and MP3 :ref:`[manual] <sec:convocc>`                                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp4(sdq)                | 4th-order MP perturbation theory (MP4) less triples :ref:`[manual] <sec:fnompn>`      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp4                     | full MP4 :ref:`[manual] <sec:fnompn>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp\ *n*                 | *n*\ th-order Moller--Plesset (MP) perturbation theory :ref:`[manual] <sec:arbpt>`    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | zapt\ *n*               | *n*\ th-order z-averaged perturbation theory (ZAPT) :ref:`[manual] <sec:arbpt>`       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | omp2                    | orbital-optimized second-order MP perturbation theory :ref:`[manual] <sec:occ>`       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | omp3                    | orbital-optimized third-order MP perturbation theory :ref:`[manual] <sec:occ>`        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | omp2.5                  | orbital-optimized MP2.5 :ref:`[manual] <sec:occ>`                                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ocepa                   | orbital-optimized coupled electron pair approximation :ref:`[manual] <sec:occ>`       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cepa0                   | coupled electron pair approximation, equiv. linear. CCD :ref:`[manual] <sec:convocc>` |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cepa(0)                 | coupled electron pair approximation variant 0 :ref:`[manual] <sec:fnocepa>`           |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cepa(1)                 | coupled electron pair approximation variant 1 :ref:`[manual] <sec:fnocepa>`           |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cepa(3)                 | coupled electron pair approximation variant 3 :ref:`[manual] <sec:fnocepa>`           |
    +-------------------------+---------------------------------------------------------------------------------------+
    | acpf                    | averaged coupled-pair functional :ref:`[manual] <sec:fnocepa>`                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | aqcc                    | averaged quadratic coupled cluster :ref:`[manual] <sec:fnocepa>`                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | qcisd                   | quadratic CI singles doubles (QCISD) :ref:`[manual] <sec:fnocc>`                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cc2                     | approximate coupled cluster singles and doubles (CC2) :ref:`[manual] <sec:cc>`        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD) :ref:`[manual] <sec:cc>`                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | bccd                    | Brueckner coupled cluster doubles (BCCD) :ref:`[manual] <sec:cc>`                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | qcisd(t)                | QCISD with perturbative triples :ref:`[manual] <sec:fnocc>`                           |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd(t)                 | CCSD with perturbative triples (CCSD(T)) :ref:`[manual] <sec:cc>`                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | fno-df-ccsd(t)          | CCSD(T) with density fitting and frozen natural orbitals :ref:`[manual] <sec:fnocc>`  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | bccd(t)                 | BCCD with perturbative triples :ref:`[manual] <sec:cc>`                               |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cc3                     | approximate CC singles, doubles, and triples (CC3) :ref:`[manual] <sec:cc>`           |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccenergy                | **expert** full control over ccenergy module                                          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cisd                    | configuration interaction (CI) singles and doubles (CISD) :ref:`[manual] <sec:ci>`    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cisdt                   | CI singles, doubles, and triples (CISDT) :ref:`[manual] <sec:ci>`                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cisdtq                  | CI singles, doubles, triples, and quadruples (CISDTQ) :ref:`[manual] <sec:ci>`        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ci\ *n*                 | *n*\ th-order CI :ref:`[manual] <sec:ci>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------+
    | fci                     | full configuration interaction (FCI) :ref:`[manual] <sec:ci>`                         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | detci                   | **expert** full control over detci module                                             |
    +-------------------------+---------------------------------------------------------------------------------------+
    | gaussian-2 (g2)         | gaussian-2 composite method :ref:`[manual] <sec:fnogn>`                               |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt0                   | 0th-order symmetry adapted perturbation theory (SAPT) :ref:`[manual] <sec:sapt>`      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2                   | 2nd-order SAPT, traditional definition :ref:`[manual] <sec:sapt>`                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+                  | SAPT including all 2nd-order terms :ref:`[manual] <sec:sapt>`                         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+(3)               | SAPT including perturbative triples :ref:`[manual] <sec:sapt>`                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+3                 | SAPT including all 3rd-order terms :ref:`[manual] <sec:sapt>`                         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+(ccd)             | SAPT2+ with CC-based dispersion :ref:`[manual] <sec:sapt>`                            |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)          | SAPT2+(3) with CC-based dispersion :ref:`[manual] <sec:sapt>`                         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+3(ccd)            | SAPT2+3 with CC-based dispersion :ref:`[manual] <sec:sapt>`                           |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt0-ct                | 0th-order SAPT plus charge transfer (CT) calculation :ref:`[manual] <sec:saptct>`     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2-ct                | SAPT2 plus CT :ref:`[manual] <sec:saptct>`                                            |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+-ct               | SAPT2+ plus CT :ref:`[manual] <sec:saptct>`                                           |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+(3)-ct            | SAPT2+(3) plus CT :ref:`[manual] <sec:saptct>`                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+3-ct              | SAPT2+3 plus CT :ref:`[manual] <sec:saptct>`                                          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+(ccd)-ct          | SAPT2+(CCD) plus CT :ref:`[manual] <sec:saptct>`                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)-ct       | SAPT2+(3)(CCD) plus CT :ref:`[manual] <sec:saptct>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+3(ccd)-ct         | SAPT2+3(CCD) plus CT :ref:`[manual] <sec:saptct>`                                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | adc                     | 2nd-order algebraic diagrammatic construction (ADC) :ref:`[manual] <sec:adc>`         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-cc2                 | EOM-CC2 :ref:`[manual] <sec:eomcc>`                                                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-ccsd                | equation of motion (EOM) CCSD :ref:`[manual] <sec:eomcc>`                             |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-cc3                 | EOM-CC3 :ref:`[manual] <sec:eomcc>`                                                   |
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

    .. include:: autodoc_dft_energy.rst

    .. _`table:energy_mrcc`:

    +-------------------------+---------------------------------------------------------------------------------------+
    | name                    | calls method in Kallay's MRCC program :ref:`[manual] <sec:mrcc>`                      |
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

    .. comment :type cast_up: :ref:`boolean <op_py_boolean>` or string
    .. comment :param cast_up: ``'on'`` || |dl| ``'off'`` |dr| || ``'3-21g'`` || ``'cc-pVDZ'`` || etc.

    .. comment     Indicates whether, to accelerate convergence for the scf portion of
    .. comment     the *name* calculation, a preliminary scf should be performed with a
    .. comment     small basis set (3-21G if a basis name is not supplied as keyword
    .. comment     value) followed by projection into the full target basis.

    .. comment .. deprecated:: Sept-2012
    .. comment    Use option |scf__basis_guess| instead.

    .. comment :type cast_up_df: :ref:`boolean <op_py_boolean>` or string
    .. comment :param cast_up_df: ``'on'`` || |dl| ``'off'`` |dr| || ``'cc-pVDZ-RI'`` || ``'aug-cc-pVDZ-JKFIT'`` || etc.

    .. comment     Indicates whether, when *cast_up* is active, to run the preliminary
    .. comment     scf in density-fitted mode or what fitting basis to employ (when
    .. comment     available for all elements, cc-pVDZ-RI is the default).

    .. comment .. deprecated:: Sept-2012
    .. comment    Use option |scf__df_basis_guess| instead.

    :type bypass_scf: :ref:`boolean <op_py_boolean>`
    :param bypass_scf: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether, for *name* values built atop of scf calculations,
        the scf step is skipped. Suitable when special steps are taken to get
        the scf to converge in an explicit preceeding scf step.

    :examples:

    >>> # [1] Coupled-cluster singles and doubles calculation with psi code
    >>> energy('ccsd')

    >>> # [2] Charge-transfer SAPT calculation with scf projection from small into
    >>> #     requested basis, with specified projection fitting basis
    >>> set basis_guess true
    >>> set df_basis_guess jun-cc-pVDZ-JKFIT
    >>> energy('sapt0-ct')

    >>> # [3] Arbitrary-order MPn calculation
    >>> energy('mp4')

    >>> # [4] Converge scf as singlet, then run detci as triplet upon singlet reference
    >>> molecule H2 {\\n0 1\\nH\\nH 1 0.74\\n}
    >>> energy('scf')
    >>> H2.set_multiplicity(3)
    >>> energy('detci', bypass_scf=True)

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    optstash = OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['E_CONVERGENCE'])

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
        # Set method-dependent scf convergence criteria
        if not PsiMod.has_option_changed('SCF', 'E_CONVERGENCE'):
            if procedures['energy'][lowername] == run_scf or procedures['energy'][lowername] == run_dft:
                PsiMod.set_local_option('SCF', 'E_CONVERGENCE', 6)
            else:
                PsiMod.set_local_option('SCF', 'E_CONVERGENCE', 8)
        if not PsiMod.has_option_changed('SCF', 'D_CONVERGENCE'):
            if procedures['energy'][lowername] == run_scf or procedures['energy'][lowername] == run_dft:
                PsiMod.set_local_option('SCF', 'D_CONVERGENCE', 6)
            else:
                PsiMod.set_local_option('SCF', 'D_CONVERGENCE', 8)

        # Set post-scf convergence criteria (global will cover all correlated modules)
        if not PsiMod.has_global_option_changed('E_CONVERGENCE'):
            if not procedures['energy'][lowername] == run_scf and not procedures['energy'][lowername] == run_dft:
                PsiMod.set_global_option('E_CONVERGENCE', 6)

        procedures['energy'][lowername](lowername, **kwargs)

    except KeyError:
        raise ValidationError('Energy method %s not available.' % (lowername))

    optstash.restore()
    return PsiMod.get_variable('CURRENT ENERGY')


def gradient(name, **kwargs):
    r"""Function complementary to optimize(). Carries out one gradient pass,
    deciding analytic or finite difference.

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)
    dertype = 1

    optstash = OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['E_CONVERGENCE'])

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

    # 2. Check if the user passes dertype into this function
    if 'dertype' in kwargs:
        opt_dertype = kwargs['dertype']

        if der0th.match(str(opt_dertype)):
            dertype = 0
            func = energy
        elif der1st.match(str(opt_dertype)):
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

    # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
    if not PsiMod.has_option_changed('SCF', 'E_CONVERGENCE'):
        if procedures['energy'][lowername] == run_scf or procedures['energy'][lowername] == run_dft:
            PsiMod.set_local_option('SCF', 'E_CONVERGENCE', 8)
        else:
            PsiMod.set_local_option('SCF', 'E_CONVERGENCE', 10)
    if not PsiMod.has_option_changed('SCF', 'D_CONVERGENCE'):
        if procedures['energy'][lowername] == run_scf or procedures['energy'][lowername] == run_dft:
            PsiMod.set_local_option('SCF', 'D_CONVERGENCE', 8)
        else:
            PsiMod.set_local_option('SCF', 'D_CONVERGENCE', 10)

    # Set post-scf convergence criteria (global will cover all correlated modules)
    if not PsiMod.has_global_option_changed('E_CONVERGENCE'):
        if not procedures['energy'][lowername] == run_scf and not procedures['energy'][lowername] == run_dft:
            PsiMod.set_global_option('E_CONVERGENCE', 8)

    # Does dertype indicate an analytic procedure both exists and is wanted?
    if (dertype == 1):
        # Nothing to it but to do it. Gradient information is saved
        # into the current reference wavefunction
        procedures['gradient'][lowername](lowername, **kwargs)

        if 'mode' in kwargs and kwargs['mode'].lower() == 'sow':
            raise ValidationError('Optimize execution mode \'sow\' not valid for analytic gradient calculation.')
        PsiMod.wavefunction().energy()

        optstash.restore()
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
            optstash.restore()
            return 0.0

        if (opt_mode.lower() == 'reap'):
            PsiMod.set_variable('CURRENT ENERGY', energies[-1])

        # Obtain the gradient
        PsiMod.fd_1_0(energies)
        
        # The last item in the list is the reference energy, return it
        optstash.restore()
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
    | df-mp2                  | MP2 with density fitting                                                              |
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

    optstash = OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['E_CONVERGENCE'])

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
        # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
        #   SCF properties have been set as 6/5 so as to match those 
        #       run normally through OEProp so subject to change
        if not PsiMod.has_option_changed('SCF', 'E_CONVERGENCE'):
            if procedures['energy'][lowername] == run_scf or procedures['energy'][lowername] == run_dft:
                PsiMod.set_local_option('SCF', 'E_CONVERGENCE', 6)
            else:
                PsiMod.set_local_option('SCF', 'E_CONVERGENCE', 10)
        if not PsiMod.has_option_changed('SCF', 'D_CONVERGENCE'):
            if procedures['energy'][lowername] == run_scf or procedures['energy'][lowername] == run_dft:
                PsiMod.set_local_option('SCF', 'D_CONVERGENCE', 6)
            else:
                PsiMod.set_local_option('SCF', 'D_CONVERGENCE', 10)

        # Set post-scf convergence criteria (global will cover all correlated modules)
        if not PsiMod.has_global_option_changed('E_CONVERGENCE'):
            if not procedures['energy'][lowername] == run_scf and not procedures['energy'][lowername] == run_dft:
                PsiMod.set_global_option('E_CONVERGENCE', 8)

        returnvalue = procedures['property'][lowername](lowername, **kwargs)

    except KeyError:
        raise ValidationError('Property method %s not available.' % (lowername))

    optstash.restore()
    return returnvalue


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
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT) :ref:`[manual] <sec:scf>`       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | dcft                    | density cumulant functional theory :ref:`[manual] <sec:dcft>`                         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp2                     | 2nd-order Moller-Plesset perturbation theory (MP2) :ref:`[manual] <sec:dfmp2>`        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-mp2                  | MP2 with density fitting :ref:`[manual] <sec:dfmp2>`                                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | conv-mp2                | conventional MP2 (non-density-fitting) :ref:`[manual] <sec:convocc>`                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp2.5                   | MP2.5 :ref:`[manual] <sec:convocc>`                                                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mp3                     | third-order MP perturbation theory :ref:`[manual] <sec:convocc>`                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | omp2                    | orbital-optimized second-order MP perturbation theory :ref:`[manual] <sec:occ>`       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | omp2.5                  | orbital-optimized MP2.5 :ref:`[manual] <sec:occ>`                                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | omp3                    | orbital-optimized third-order MP perturbation theory :ref:`[manual] <sec:occ>`        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ocepa                   | orbital-optimized coupled electron pair approximation :ref:`[manual] <sec:occ>`       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cepa0                   | coupled electron pair approximation(0) :ref:`[manual] <sec:convocc>`                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD) :ref:`[manual] <sec:cc>`                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ccsd(t)                 | CCSD with perturbative triples (CCSD(T)) :ref:`[manual] <sec:cc>`                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | eom-ccsd                | equation of motion (EOM) CCSD :ref:`[manual] <sec:eomcc>`                             |
    +-------------------------+---------------------------------------------------------------------------------------+

    .. include:: autodoc_dft_opt.rst

    .. warning:: Optimizations where the molecule is specified in Z-matrix format 
       with dummy atoms will result in the geometry being converted to a Cartesian representation.

    :type name: string
    :param name: ``'scf'`` || ``'df-mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the database. May be any valid argument to
        :py:func:`~driver.energy`.

    :type func: :ref:`function <op_py_function>`
    :param func: |dl| ``gradient`` |dr| || ``energy`` || ``cbs``

        Indicates the type of calculation to be performed on the molecule.
        The default dertype accesses ``'gradient'`` or ``'energy'``, while
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

    >>> # [2] Finite difference mp5 optimization
    >>> opt('mp5')

    >>> # [3] Forced finite difference ccsd optimization
    >>> optimize('ccsd', dertype=1)

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    full_hess_every = PsiMod.get_local_option('OPTKING', 'FULL_HESS_EVERY')
    steps_since_last_hessian = 0

    # are we in sow/reap mode?
    isSowReap = False
    if ('mode' in kwargs) and (kwargs['mode'].lower() == 'sow'):
        isSowReap = True
    if ('mode' in kwargs) and (kwargs['mode'].lower() == 'reap'):
        isSowReap = True
    optstash = OptionsState(
        ['SCF', 'GUESS'])

    n = 1
    if ('opt_iter' in kwargs):
        n = kwargs['opt_iter']

    PsiMod.get_active_molecule().update_geometry()
    mol = PsiMod.get_active_molecule()
    mol.update_geometry()
    initial_sym = mol.schoenflies_symbol()
    while n <= PsiMod.get_global_option('GEOM_MAXITER'):
        mol = PsiMod.get_active_molecule()
        mol.update_geometry()
        current_sym = mol.schoenflies_symbol()
        if initial_sym != current_sym:
            raise Exception("Point group changed!  You should restart using " +\
                            "the last geometry in the output, after carefully "+\
                            "making sure all symmetry-dependent information in "+\
                            "the input, such as DOCC, is correct.")
        kwargs['opt_iter'] = n

        # Use orbitals from previous iteration as a guess
        if (n > 1) and (not isSowReap):
            PsiMod.set_local_option('SCF', 'GUESS', 'READ')

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
            pass
            # Do nothing; user said to read existing hessian once
        else:
            PsiMod.set_global_option('CART_HESS_READ', False)
            steps_since_last_hessian += 1

        # print 'cart_hess_read', PsiMod.get_global_option('CART_HESS_READ')
        # Take step
        if PsiMod.optking() == PsiMod.PsiReturnType.EndLoop:
            print('Optimizer: Optimization complete!')
            PsiMod.print_out('\n    Final optimized geometry and variables:\n')
            PsiMod.get_active_molecule().print_in_input_format()
            # Check if user wants to see the intcos; if so, don't delete them.
            if (PsiMod.get_option('OPTKING', 'INTCOS_GENERATE_EXIT') == False):
                PsiMod.opt_clean()
            PsiMod.clean()

            # S/R: Clean up opt input file
            if ('mode' in kwargs) and (kwargs['mode'].lower() == 'reap'):
                fmaster = open('OPT-master.in', 'w')
                fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n')
                fmaster.write('# Optimization complete!\n\n')
                fmaster.close()

            optstash.restore()
            return thisenergy
        PsiMod.print_out('\n    Structure for next step:\n')
        PsiMod.get_active_molecule().print_in_input_format()

        # S/R: Preserve opt data file for next pass and switch modes to get new displacements
        if ('mode' in kwargs) and (kwargs['mode'].lower() == 'reap'):
            kwargs['opt_datafile'] = get_psifile(1)
            kwargs['mode'] = 'sow'

        n += 1

    PsiMod.print_out('\tOptimizer: Did not converge!')

    optstash.restore()
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
            # Let 'mp2' and 'mp3' pass through as themselves to occ module
            if (namestump == 'mp') and ((namelevel == 2) or (namelevel == 3)):
                return namelower, None
            # Let 'mp4' be redirected to fnocc module if rhf
            elif (namestump == 'mp') and (namelevel == 4):
                if PsiMod.get_option('SCF','REFERENCE') == 'RHF':
                    return 'fnocc-mp', 4
                else:
                    return 'detci-mp', 4
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
    dertype = 2

    optstash = OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['E_CONVERGENCE'])

    # Order of precedence:
    #    1. Default for wavefunction
    #    2. Value obtained from kwargs, if user changed it
    #    3. If user provides a custom 'func' use that

    # Allow specification of methods to arbitrary order
    lowername, level = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    # 1. set the default to that of the provided name
    if lowername in procedures['hessian']:
        dertype = 2
    elif lowername in procedures['gradient']:
        dertype = 1
        func = gradient
    elif lowername in procedures['energy']:
        dertype = 0
        func = energy

    # 2. Check if the user passes dertype into this function
    if 'dertype' in kwargs:
        freq_dertype = kwargs['dertype']

        if der0th.match(str(freq_dertype)):
            dertype = 0
            func = energy
        elif der1st.match(str(freq_dertype)):
            dertype = 1
            func = gradient
        elif der2nd.match(str(freq_dertype)):
            dertype = 2
        else:
            raise ValidationError('Requested derivative level \'dertype\' %s not valid for helper function frequency.' % (freq_dertype))

    # 3. if the user provides a custom function THAT takes precedence
    if ('freq_func' in kwargs) or ('func' in kwargs):
        if ('func' in kwargs):
            kwargs['freq_func'] = kwargs['func']
            del kwargs['func']
        dertype = 0
        func = kwargs['freq_func']

    # Summary validation
    if (dertype == 2) and (lowername in procedures['hessian']):
        pass
    elif (dertype == 1) and (func is gradient) and (lowername in procedures['gradient']):
        pass
    elif (dertype == 1) and not(func is gradient):
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

    # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
    if not PsiMod.has_option_changed('SCF', 'E_CONVERGENCE'):
        if procedures['energy'][lowername] == run_scf or procedures['energy'][lowername] == run_dft:
            PsiMod.set_local_option('SCF', 'E_CONVERGENCE', 8)
        else:
            PsiMod.set_local_option('SCF', 'E_CONVERGENCE', 10)
    if not PsiMod.has_option_changed('SCF', 'D_CONVERGENCE'):
        if procedures['energy'][lowername] == run_scf or procedures['energy'][lowername] == run_dft:
            PsiMod.set_local_option('SCF', 'D_CONVERGENCE', 8)
        else:
            PsiMod.set_local_option('SCF', 'D_CONVERGENCE', 10)

    # Set post-scf convergence criteria (global will cover all correlated modules)
    if not PsiMod.has_global_option_changed('E_CONVERGENCE'):
        if not procedures['energy'][lowername] == run_scf and not procedures['energy'][lowername] == run_dft:
            PsiMod.set_global_option('E_CONVERGENCE', 8)

    # Select certain irreps
    if 'irrep' in kwargs:
        irrep = parse_cotton_irreps(kwargs['irrep']) - 1  # externally, A1 irrep is 1, internally 0
    else:
        irrep = -1  # -1 implies do all irreps

    # Does an analytic procedure exist for the requested method?
    if (dertype == 2):
        # We have the desired method. Do it.
        procedures['hessian'][lowername](lowername, **kwargs)
        optstash.restore()

        # call thermo module
        PsiMod.thermo()

        return PsiMod.wavefunction().energy()

    elif (dertype == 1):
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
            sys.stdout.flush()

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

        # call thermo module
        PsiMod.thermo()

        optstash.restore()
        # TODO: add return statement

    else:  # Assume energy points
        # If not, perform finite difference of energies
        info = 'Performing finite difference calculations by energies'
        print(info)

        # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
        optstash.restore()
        if not PsiMod.has_option_changed('SCF', 'E_CONVERGENCE'):
            if procedures['energy'][lowername] == run_scf or procedures['energy'][lowername] == run_dft:
                PsiMod.set_local_option('SCF', 'E_CONVERGENCE', 10)
            else:
                PsiMod.set_local_option('SCF', 'E_CONVERGENCE', 11)
        if not PsiMod.has_option_changed('SCF', 'D_CONVERGENCE'):
            if procedures['energy'][lowername] == run_scf or procedures['energy'][lowername] == run_dft:
                PsiMod.set_local_option('SCF', 'D_CONVERGENCE', 10)
            else:
                PsiMod.set_local_option('SCF', 'D_CONVERGENCE', 11)

        # Set post-scf convergence criteria (global will cover all correlated modules)
        if not PsiMod.has_global_option_changed('E_CONVERGENCE'):
            if not procedures['energy'][lowername] == run_scf and not procedures['energy'][lowername] == run_dft:
                PsiMod.set_global_option('E_CONVERGENCE', 10)

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
            sys.stdout.flush()

            # Load in displacement into the active molecule
            molecule.set_geometry(displacement)

            # Perform the energy calculation
            E = func(lowername, **kwargs)

            # Save the energy
            energies.append(E)

            # clean may be necessary when changing irreps of displacements
            PsiMod.clean()

        # Obtain the gradient. This function stores the gradient in the wavefunction.
        PsiMod.fd_freq_0(energies, irrep)

        print(' Computation complete.')

        # Clear the "parent" symmetry now
        PsiMod.set_parent_symmetry("")

        # TODO: These need to be restored to the user specified setting
        PsiMod.get_active_molecule().fix_orientation(False)
        # But not this one, it always goes back to True
        PsiMod.get_active_molecule().reinterpret_coordentry(True)

        # The last item in the list is the reference energy, return it
        optstash.restore()

        # Clear the "parent" symmetry now
        PsiMod.set_parent_symmetry("")

        # call thermo module
        PsiMod.thermo()

        optstash.restore()
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
    """Function to write wavefunction information in molden
    format to *filename*

    """
    m = PsiMod.MoldenWriter(PsiMod.wavefunction())
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
