
#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

from __future__ import print_function
"""Module with a *procedures* dictionary specifying available quantum
chemical methods and functions driving the main quantum chemical
functionality, namely single-point energies, geometry optimizations,
properties, and vibrational frequency calculations.

"""
from __future__ import absolute_import
import sys
import re
#CUimport psi4
#CUimport p4util
#CUimport p4const
from proc import *
from interface_cfour import *
#CUfrom functional import *
#CUfrom p4regex import *
# never import wrappers or aliases into this file


# Procedure lookup tables
procedures = {
        'energy': {
            'scf'           : run_scf,
            'mcscf'         : run_mcscf,
            'dcft'          : run_dcft,
            'mp3'           : select_mp3,
            'mp2.5'         : select_mp2p5,
            'mp2'           : select_mp2,
            'omp2'          : select_omp2,
            'scs-omp2'      : run_occ,
            'scs(n)-omp2'   : run_occ,
            'scs-omp2-vdw'  : run_occ,
            'sos-omp2'      : run_occ,
            'sos-pi-omp2'   : run_occ,
            'omp3'          : select_omp3,
            'scs-omp3'      : run_occ,
            'scs(n)-omp3'   : run_occ,
            'scs-omp3-vdw'  : run_occ,
            'sos-omp3'      : run_occ,
            'sos-pi-omp3'   : run_occ,
            'olccd'         : select_olccd,
            'omp2.5'        : select_omp2p5,
            'dfocc'         : run_dfocc,
            'qchf'          : run_qchf,
            'ccd'           : run_dfocc,
            'sapt0'         : run_sapt,
            'sapt2'         : run_sapt,
            'sapt2+'        : run_sapt,
            'sapt2+(3)'     : run_sapt,
            'sapt2+3'       : run_sapt,
            'sapt2+(ccd)'   : run_sapt,
            'sapt2+(3)(ccd)': run_sapt,
            'sapt2+3(ccd)'  : run_sapt,
            'sapt2+dmp2'    : run_sapt,
            'sapt2+(3)dmp2' : run_sapt,
            'sapt2+3dmp2'   : run_sapt,
            'sapt2+(ccd)dmp2' : run_sapt,
            'sapt2+(3)(ccd)dmp2' : run_sapt,
            'sapt2+3(ccd)dmp2' : run_sapt,
            'sapt0-ct'      : run_sapt_ct,
            'sapt2-ct'      : run_sapt_ct,
            'sapt2+-ct'     : run_sapt_ct,
            'sapt2+(3)-ct'  : run_sapt_ct,
            'sapt2+3-ct'    : run_sapt_ct,
            'sapt2+(ccd)-ct'     : run_sapt_ct,
            'sapt2+(3)(ccd)-ct'  : run_sapt_ct,
            'sapt2+3(ccd)-ct'    : run_sapt_ct,
            'fisapt0'       : run_fisapt,
            'ccenergy'      : run_ccenergy,  # full control over ccenergy
            'ccsd'          : select_ccsd,
            'ccsd(t)'       : select_ccsd_t_,
            'ccsd(at)'      : select_ccsd_at_,
            'cc2'           : run_ccenergy,
            'cc3'           : run_ccenergy,
            'mrcc'          : run_mrcc,  # interface to Kallay's MRCC program
            'bccd'          : run_bccd,
            'bccd(t)'       : run_bccd,
            'eom-ccsd'      : run_eom_cc,
            'eom-cc2'       : run_eom_cc,
            'eom-cc3'       : run_eom_cc,
            'detci'         : run_detci,  # full control over detci
            'mp'            : run_detci,  # arbitrary order mp(n)
            'zapt'          : run_detci,  # arbitrary order zapt(n)
            'cisd'          : select_cisd,
            'cisdt'         : run_detci,
            'cisdtq'        : run_detci,
            'ci'            : run_detci,  # arbitrary order ci(n)
            'fci'           : run_detci,
            'casscf'        : run_detcas,
            'rasscf'        : run_detcas,
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
            'qcisd'         : run_fnocc,
            'qcisd(t)'      : run_fnocc,
            'mp4'           : select_mp4,
            'mp4(sdq)'      : run_fnocc,
            'fno-ccsd'      : select_fnoccsd,
            'fno-ccsd(t)'   : select_fnoccsd_t_,
            'fno-qcisd'     : run_fnocc,
            'fno-qcisd(t)'  : run_fnocc,
            'fno-mp3'       : run_fnocc,
            'fno-mp4(sdq)'  : run_fnocc,
            'fno-mp4'       : run_fnocc,
            'fno-lccd'      : run_cepa,
            'fno-lccsd'     : run_cepa,
            'fno-cepa(0)'   : run_cepa,
            'fno-cepa(1)'   : run_cepa,
            'fno-cepa(3)'   : run_cepa,
            'fno-acpf'      : run_cepa,
            'fno-aqcc'      : run_cepa,
            'fno-cisd'      : run_cepa,
            'lccd'          : select_lccd,
            'lccsd'         : run_cepa,
            'cepa(0)'       : run_cepa,
            'cepa(1)'       : run_cepa,
            'cepa(3)'       : run_cepa,
            'acpf'          : run_cepa,
            'aqcc'          : run_cepa,
            'efp'           : run_efp,
            'dmrgscf'       : run_dmrgscf,
            'dmrgci'        : run_dmrgci,
            # Upon adding a method to this list, add it to the docstring in energy() below
            # Aliases are discouraged. If you must add an alias to this list (e.g.,
            #    lccsd/cepa(0)), please search the whole driver to find uses of
            #    name in return values and psi variables and extend the logic to
            #    encompass the new alias.
        },
        'gradient' : {
            'scf'           : run_scf_gradient,
            'ccsd'          : select_ccsd_gradient,
            'ccsd(t)'       : select_ccsd_t__gradient,
            'mp2'           : select_mp2_gradient,
            'eom-ccsd'      : run_eom_cc_gradient,
            'dcft'          : run_dcft_gradient,
            'omp2'          : select_omp2_gradient,
            'omp3'          : select_omp3_gradient,
            'mp3'           : select_mp3_gradient,
            'mp2.5'         : select_mp2p5_gradient,
            'omp2.5'        : select_omp2p5_gradient,
            'lccd'          : select_lccd_gradient,
            'olccd'         : select_olccd_gradient,
            'ccd'           : run_dfocc_gradient,
            'hf'            : run_scf_gradient,
            # Upon adding a method to this list, add it to the docstring in optimize() below
        },
        'hessian' : {
            # Upon adding a method to this list, add it to the docstring in frequency() below
        },
        'property' : {
            'scf'      : run_scf_property,
            'cc2'      : run_cc_property,
            'ccsd'     : run_cc_property,
            'df-mp2'   : run_dfmp2_property,
            'dfmp2'    : run_dfmp2_property,
            'ri-mp2'   : run_dfocc_property,
            'df-omp2'  : run_dfocc_property,
            'eom-cc2'  : run_cc_property,
            'eom-ccsd' : run_cc_property,
            'detci'    : run_detci_property,  # full control over detci
            'cisd'     : run_detci_property,
            'cisdt'    : run_detci_property,
            'cisdtq'   : run_detci_property,
            'ci'       : run_detci_property,  # arbitrary order ci(n)
            'fci'      : run_detci_property,
            'hf'       : run_scf_property,
            # Upon adding a method to this list, add it to the docstring in property() below
        }}

# Will only allow energy to be run for the following methods
energy_only_methods = [x for x in procedures['energy'].keys() if 'sapt' in x]
energy_only_methods += ['adc', 'efp', 'cphf', 'tdhf', 'cis'] 

# dictionary to register pre- and post-compute hooks for driver routines
hooks = dict((k1, dict((k2, []) for k2 in ['pre', 'post'])) for k1 in ['energy', 'optimize', 'frequency'])

# Integrate DFT with driver routines
for ssuper in superfunctional_list():
    procedures['energy'][ssuper.name().lower()] = run_dft

for ssuper in superfunctional_list():
    if ((not ssuper.is_c_hybrid()) and (not ssuper.is_c_lrc()) and (not ssuper.is_x_lrc())):
        procedures['gradient'][ssuper.name().lower()] = run_dft_gradient

# Integrate CFOUR with driver routines
for ssuper in cfour_list():
    procedures['energy'][ssuper] = run_cfour

for ssuper in cfour_gradient_list():
    procedures['gradient'][ssuper] = run_cfour


def energy(name, **kwargs):
    r"""Function to compute the single-point electronic energy.

    :returns: *float* |w--w| Total electronic energy in Hartrees. SAPT & EFP return interaction energy.

    :returns: (*float*, :ref:`Wavefunction<sec:psimod_Wavefunction>`) |w--w| energy and wavefunction when **return_wfn** specified.

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`CURRENT ENERGY <CURRENTENERGY>`
       * :psivar:`CURRENT REFERENCE ENERGY <CURRENTREFERENCEENERGY>`
       * :psivar:`CURRENT CORRELATION ENERGY <CURRENTCORRELATIONENERGY>`

    :type name: string
    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type return_wfn: :ref:`boolean <op_py_boolean>`
    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|

        Indicate to additionally return the :ref:`Wavefunction<sec:psimod_Wavefunction>`
        calculation result as the second element (after *float* energy) of a tuple.

    :type restart_file: string
    :param restart_file: ``['file.1, file.32]`` || ``./file`` || etc.

        Binary data files to be renamed for calculation restart.


    .. _`table:energy_gen`:

    +-------------------------+---------------------------------------------------------------------------------------+
    | name                    | calls method                                                                          |
    +=========================+=======================================================================================+
    | efp                     | effective fragment potential (EFP) :ref:`[manual] <sec:efp>`                          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT) :ref:`[manual] <sec:scf>`       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | hf                      | HF self consistent field (SCF) :ref:`[manual] <sec:scf>`                              |
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
    | df-omp2                 | density-fitted orbital-optimized MP2 :ref:`[manual] <sec:dfocc>`                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-omp3                 | density-fitted orbital-optimized MP3 :ref:`[manual] <sec:dfocc>`                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-omp2.5               | density-fitted orbital-optimized MP2.5 :ref:`[manual] <sec:dfocc>`                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cd-omp2                 | cholesky decomposed orbital-optimized MP2 :ref:`[manual] <sec:dfocc>`                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cd-omp3                 | cholesky decomposed orbital-optimized MP3 :ref:`[manual] <sec:dfocc>`                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cd-omp2.5               | cholesky decomposed orbital-optimized MP2.5 :ref:`[manual] <sec:dfocc>`               |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cd-mp2                  | cholesky decomposed MP2 :ref:`[manual] <sec:dfocc>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cd-mp3                  | cholesky decomposed MP3 :ref:`[manual] <sec:dfocc>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------+
    | cd-mp2.5                | cholesky decomposed MP2.5 :ref:`[manual] <sec:dfocc>`                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-ccsd2                | density-fitted CCSD from DFOCC module :ref:`[manual] <sec:dfocc>`                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | ri-ccsd(t)              | density-fitted CCSD(T) from DFOCC module :ref:`[manual] <sec:dfocc>`                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-ccsd(at)             | density-fitted Lambda-CCSD(T) from DFOCC module :ref:`[manual] <sec:dfocc>`           |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-ccd                  | density-fitted CCD from DFOCC module :ref:`[manual] <sec:dfocc>`                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-mp3                  | density-fitted MP3 from DFOCC module :ref:`[manual] <sec:dfocc>`                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-mp2.5                | density-fitted MP2.5 from DFOCC module :ref:`[manual] <sec:dfocc>`                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | qchf                    | density-fitted QC-HF from DFOCC module :ref:`[manual] <sec:dfocc>`                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-ccsdl                | density-fitted CCSDL from DFOCC module :ref:`[manual] <sec:dfocc>`                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-ccdl                 | density-fitted CCDL from DFOCC module :ref:`[manual] <sec:dfocc>`                     |
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
    | casscf                  | complete active space self consistent field (CASSCF)  :ref:`[manual] <sec:cas>`       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | rasscf                  | restricted active space self consistent field (RASSCF)  :ref:`[manual] <sec:cas>`     |
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
    | sapt2+dmp2              | SAPT including all 2nd-order terms and MP2 correction :ref:`[manual] <sec:sapt>`      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+(3)dmp2           | SAPT including perturbative triples and MP2 correction :ref:`[manual] <sec:sapt>`     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+3dmp2             | SAPT including all 3rd-order terms and MP2 correction :ref:`[manual] <sec:sapt>`      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+(ccd)dmp2         | SAPT2+ with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)dmp2      | SAPT2+(3) with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | sapt2+3(ccd)dmp2        | SAPT2+3 with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`        |
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

    .. include:: autodoc_dft_energy.rst

    .. include:: mrcc_table_energy.rst

    .. include:: cfour_table_energy.rst

    :examples:

    >>> # [1] Coupled-cluster singles and doubles calculation with psi code
    >>> energy('ccsd')

    >>> # [2] Charge-transfer SAPT calculation with scf projection from small into
    >>> #     requested basis, with specified projection fitting basis
    >>> set basis_guess true
    >>> set df_basis_guess jun-cc-pVDZ-JKFIT
    >>> energy('sapt0-ct')

    >>> # [3] Arbitrary-order MPn calculation
    >>> energy('mp7')

    >>> # [4] Converge scf as singlet, then run detci as triplet upon singlet reference
    >>> # Note that the integral transformation is not done automatically when detci is run in a separate step.
    >>> molecule H2 {\n0 1\nH\nH 1 0.74\n}
    >>> set global basis cc-pVDZ
    >>> set global reference rohf
    >>> energy('scf')
    >>> H2.set_multiplicity(3)
    >>> psi4.MintsHelper().integrals()
    >>> energy('detci', bypass_scf=True)

    >>> # [5] Run two CI calculations, keeping the integrals generated in the first one.
    >>> molecule ne {\\nNe\\n}
    >>> set globals  basis cc-pVDZ
    >>> set transqt2 delete_tei false
    >>> energy('cisd')
    >>> energy('fci', bypass_scf='True')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    psi4.clean_variables()

    optstash = p4util.OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['E_CONVERGENCE'])

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # Allow specification of methods to arbitrary order
    lowername, level = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    for precallback in hooks['energy']['pre']:
        precallback(lowername, **kwargs)

    try:
        # Set method-dependent scf convergence criteria
        if not psi4.has_option_changed('SCF', 'E_CONVERGENCE'):
            if procedures['energy'][lowername] in [run_scf, run_dft]:
                psi4.set_local_option('SCF', 'E_CONVERGENCE', 6)
            else:
                psi4.set_local_option('SCF', 'E_CONVERGENCE', 8)
        if not psi4.has_option_changed('SCF', 'D_CONVERGENCE'):
            if procedures['energy'][lowername] in [run_scf, run_dft]:
                psi4.set_local_option('SCF', 'D_CONVERGENCE', 6)
            else:
                psi4.set_local_option('SCF', 'D_CONVERGENCE', 8)

        # Set post-scf convergence criteria (global will cover all correlated modules)
        if not psi4.has_global_option_changed('E_CONVERGENCE'):
            if procedures['energy'][lowername] not in [run_scf, run_dft]:
                psi4.set_global_option('E_CONVERGENCE', 6)

# Before invoking the procedure, we rename any file that should be read.
# This is a workaround to do restarts with the current PSI4 capabilities
# before actual, clean restarts are put in there
# Restartfile is always converted to a single-element list if
# it contains a single string
        if 'restart_file' in kwargs:
            restartfile = kwargs['restart_file']  # Option still available for procedure-specific action
            if restartfile != list(restartfile):
                restartfile = [restartfile]
            # Rename the files to be read to be consistent with psi4's file system
            for item in restartfile:
                name_split = re.split(r'\.', item)
                filenum = name_split[len(name_split) - 1]
                try:
                    filenum = int(filenum)
                except ValueError:
                    filenum = 32  # Default file number is the checkpoint one
                psioh = psi4.IOManager.shared_object()
                psio = psi4.IO.shared_object()
                filepath = psioh.get_file_path(filenum)
                namespace = psio.get_default_namespace()
                pid = str(os.getpid())
                prefix = 'psi'
                targetfile = filepath + prefix + '.' + pid + '.' + namespace + '.' + str(filenum)
                shutil.copy(item, targetfile)

        wfn = procedures['energy'][lowername](lowername, molecule=molecule, **kwargs)

    except KeyError:
        alternatives = ""
        alt_lowername = p4util.text.find_approximate_string_matches(lowername, procedures['energy'].keys(), 2)
        if len(alt_lowername) > 0:
            alternatives = " Did you mean? %s" % (" ".join(alt_lowername))
        raise ValidationError('Energy method %s not available.%s' % (lowername, alternatives))

    for postcallback in hooks['energy']['post']:
        postcallback(lowername, **kwargs)

    optstash.restore()
    if return_wfn:  # TODO current energy safer than wfn.energy() for now, but should be revisited

        # TODO place this with the associated call, very awkward to call this in other areas at the moment
        if name.lower() in ['EFP', 'MRCC', 'DMRG', 'PSIMRCC']:
            psi4.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
            psi4.print_out("The returned wavefunction is the incoming reference wavefunction.\n\n")
        elif 'sapt' in name.lower():
            psi4.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
            psi4.print_out("The returned wavefunction is the dimer SCF wavefunction.\n\n")

        return (psi4.get_variable('CURRENT ENERGY'), wfn)
    else:
        return psi4.get_variable('CURRENT ENERGY')


def gradient(name, **kwargs):
    r"""Function complementary to optimize(). Carries out one gradient pass,
    deciding analytic or finite difference.

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    psi4.clean_variables()
    dertype = 1

    # Prevent methods that do not have associated energies 
    if lowername in energy_only_methods:
	raise ValidationError("gradient('%s') does not have an associated gradient" % name)

    optstash = p4util.OptionsState(
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
        if procedures['gradient'][lowername].__name__.startswith('select_'):
            try:
                procedures['gradient'][lowername](lowername, probe=True)
            except ManagedMethodError:
                dertype = 0
                func = energy
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
            raise ValidationError("""Derivative level 'dertype' %s not valid for helper function optimize.""" % (opt_dertype))

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
        alternatives = ''
        alt_lowername = p4util.text.find_approximate_string_matches(lowername, procedures['gradient'].keys(), 2)
        if len(alt_lowername) > 0:
            alternatives = " Did you mean? %s" % (" ".join(alt_lowername))
        raise ValidationError("""Derivative method 'name' %s and derivative level 'dertype' %s are not available.%s"""
            % (lowername, dertype, alternatives))

    # no analytic derivatives for scf_type cd
    if psi4.get_option('SCF', 'SCF_TYPE') == 'CD':
        if (dertype == 1):
            raise ValidationError("""No analytic derivatives for SCF_TYPE CD.""")

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # S/R: Mode of operation- whether finite difference opt run in one job or files farmed out
    opt_mode = kwargs.get('mode', 'continuous').lower()
    if opt_mode == 'continuous':
        pass
    elif opt_mode == 'sow':
        if dertype == 1:
            raise ValidationError("""Optimize execution mode 'sow' not valid for analytic gradient calculation.""")
    elif opt_mode == 'reap':
        opt_linkage = kwargs.get('linkage', None)
        if opt_linkage is None:
            raise ValidationError("""Optimize execution mode 'reap' requires a linkage option.""")
    else:
        raise ValidationError("""Optimize execution mode '%s' not valid.""" % (opt_mode))

    # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
    if not psi4.has_option_changed('SCF', 'E_CONVERGENCE'):
        if procedures['energy'][lowername] in [run_scf, run_dft]:
            psi4.set_local_option('SCF', 'E_CONVERGENCE', 8)
        else:
            psi4.set_local_option('SCF', 'E_CONVERGENCE', 10)
    if not psi4.has_option_changed('SCF', 'D_CONVERGENCE'):
        if procedures['energy'][lowername] in [run_scf, run_dft]:
            psi4.set_local_option('SCF', 'D_CONVERGENCE', 8)
        else:
            psi4.set_local_option('SCF', 'D_CONVERGENCE', 10)

    # Set post-scf convergence criteria (global will cover all correlated modules)
    if not psi4.has_global_option_changed('E_CONVERGENCE'):
        if procedures['energy'][lowername] not in [run_scf, run_dft]:
            psi4.set_global_option('E_CONVERGENCE', 8)

    # Does dertype indicate an analytic procedure both exists and is wanted?
    if dertype == 1:
        psi4.print_out("""gradient() will perform analytic gradient computation.\n""")

        # Perform the gradient calculation
        wfn = procedures['gradient'][lowername](lowername, molecule=molecule, **kwargs)

        optstash.restore()
        if return_wfn:
            return (wfn.gradient(), wfn)
        else:
            return wfn.gradient()

    else:
        psi4.print_out("""gradient() will perform gradient computation by finite difference of analytic energies.\n""")

        opt_iter = kwargs.get('opt_iter', 1)

        if opt_iter == 1:
            print('Performing finite difference calculations')

        # Shifting the geometry so need to copy the active molecule
        moleculeclone = molecule.clone()

        # Obtain list of displacements
        displacements = psi4.fd_geoms_1_0(moleculeclone)
        ndisp = len(displacements)

        # This version is pretty dependent on the reference geometry being last (as it is now)
        print(""" %d displacements needed ...""" % (ndisp), end='')
        energies = []

        # S/R: Write instructions for sow/reap procedure to output file and reap input file
        if opt_mode == 'sow':
            instructionsO = """\n    The optimization sow/reap procedure has been selected through mode='sow'. In addition\n"""
            instructionsO += """    to this output file (which contains no quantum chemical calculations), this job\n"""
            instructionsO += """    has produced a number of input files (OPT-%s-*.in) for individual components\n""" % (str(opt_iter))
            instructionsO += """    and a single input file (OPT-master.in) with an optimize(mode='reap') command.\n"""
            instructionsO += """    These files may look very peculiar since they contain processed and pickled python\n"""
            instructionsO += """    rather than normal input. Follow the instructions in OPT-master.in to continue.\n\n"""
            instructionsO += """    Alternatively, a single-job execution of the gradient may be accessed through\n"""
            instructionsO += """    the optimization wrapper option mode='continuous'.\n\n"""
            psi4.print_out(instructionsO)

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

            with open('OPT-master.in', 'wb') as fmaster:
                fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n'.encode('utf-8'))
                fmaster.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
                fmaster.write(p4util.format_options_for_input().encode('utf-8'))
                p4util.format_kwargs_for_input(fmaster, lmode=2, return_wfn=True, **kwargs)
                fmaster.write(("""retE, retwfn = %s('%s', **kwargs)\n\n""" % (optimize.__name__, lowername)).encode('utf-8'))
                fmaster.write(instructionsM.encode('utf-8'))

        for n, displacement in enumerate(displacements):
            rfile = 'OPT-%s-%s' % (opt_iter, n + 1)

            # Build string of title banner
            banners = ''
            banners += """psi4.print_out('\\n')\n"""
            banners += """p4util.banner(' Gradient %d Computation: Displacement %d ')\n""" % (opt_iter, n + 1)
            banners += """psi4.print_out('\\n')\n\n"""

            if opt_mode == 'continuous':

                # print progress to file and screen
                psi4.print_out('\n')
                p4util.banner('Loading displacement %d of %d' % (n + 1, ndisp))
                print(""" %d""" % (n + 1), end=('\n' if (n + 1 == ndisp) else ''))
                sys.stdout.flush()

                # Load in displacement into the active molecule
                moleculeclone.set_geometry(displacement)

                # Perform the energy calculation
                E, wfn = func(lowername, return_wfn=True, molecule=moleculeclone, **kwargs)
                energies.append(psi4.get_variable('CURRENT ENERGY'))

            # S/R: Write each displaced geometry to an input file
            elif opt_mode == 'sow':
                moleculeclone.set_geometry(displacement)

                # S/R: Prepare molecule, options, and kwargs
                with open('%s.in' % (rfile), 'wb') as freagent:
                    freagent.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n'.encode('utf-8'))
                    freagent.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
                    freagent.write(p4util.format_options_for_input().encode('utf-8'))
                    p4util.format_kwargs_for_input(freagent, **kwargs)

                    # S/R: Prepare function call and energy save
                    freagent.write(("""electronic_energy = %s('%s', **kwargs)\n\n""" % (func.__name__, lowername)).encode('utf-8'))
                    freagent.write(("""psi4.print_out('\\nGRADIENT RESULT: computation %d for item %d """ % (os.getpid(), n + 1)).encode('utf-8'))
                    freagent.write("""yields electronic energy %20.12f\\n' % (electronic_energy))\n\n""".encode('utf-8'))

            # S/R: Read energy from each displaced geometry output file and save in energies array
            elif opt_mode == 'reap':
                exec(banners)
                psi4.set_variable('NUCLEAR REPULSION ENERGY', moleculeclone.nuclear_repulsion_energy())
                energies.append(p4util.extract_sowreap_from_output(rfile, 'GRADIENT', n, opt_linkage, True))

        # S/R: Quit sow after writing files. Initialize skeleton wfn to receive grad for reap
        if opt_mode == 'sow':
            optstash.restore()
            if return_wfn:
                return (None, None)  # any point to building a dummy wfn here?
            else:
                return None
        elif opt_mode == 'reap':
            psi4.set_variable('CURRENT ENERGY', energies[-1])
            wfn = psi4.new_wavefunction(molecule, psi4.get_global_option('BASIS'))

        # Compute the gradient; last item in 'energies' is undisplaced
        psi4.set_local_option('FINDIF', 'GRADIENT_WRITE', True)
        G = psi4.fd_1_0(molecule, energies)
        G.print_out()
        wfn.set_gradient(G)

        optstash.restore()

        if return_wfn:
            return (wfn.gradient(), wfn)
        else:
            return wfn.gradient()


def property(name, **kwargs):
    r"""Function to compute various properties.

    :aliases: prop()

    :returns: none.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - This function at present has a limited functionality.
         Consult the keywords sections of other modules for further property capabilities.

    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | Name               | Calls Method                                  | Reference      | Supported Properties                                          |
    +====================+===============================================+================+===============================================================+
    | scf                | Self-consistent field method(s)               | RHF/ROHF/UHF   | Listed :ref:`here <sec:oeprop>`                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | hf                 | HF Self-consistent field method(s)            | RHF/ROHF/UHF   | Listed :ref:`here <sec:oeprop>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | cc2                | 2nd-order approximate CCSD                    | RHF            | dipole, quadrupole, polarizability, rotation, roa             |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | ccsd               | Coupled cluster singles and doubles (CCSD)    | RHF            | dipole, quadrupole, polarizability, rotation, roa             |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | df-mp2             | MP2 with density fitting                      | RHF            | dipole, quadrupole, mulliken_charges, no_occupations          |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | eom-cc2            | 2nd-order approximate EOM-CCSD                | RHF            | oscillator_strength, rotational_strength                      |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | eom-ccsd           | Equation-of-motion CCSD (EOM-CCSD)            | RHF            | oscillator_strength, rotational_strength                      |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | 'cisd', 'cisdt',   | Configuration interaction                     | RHF/ROHF       | dipole, quadrupole, transition_dipole, transition_quadrupole  |
    | 'cisdt', 'cisdtq', |                                               |                |                                                               |
    | 'ci5', etc...      |                                               |                |                                                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | 'fci'              | Full configuration interaction                | RHF/ROHF       | dipole, quadrupole, transition_dipole, transition_quadrupole  |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+

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
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)

    optstash = p4util.OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['E_CONVERGENCE'])

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # Allow specification of methods to arbitrary order
    lowername, level = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    try:
        # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
        #   SCF properties have been set as 6/5 so as to match those
        #       run normally through OEProp so subject to change
        if not psi4.has_option_changed('SCF', 'E_CONVERGENCE'):
            if procedures['energy'][lowername] in [run_scf, run_dft]:
                psi4.set_local_option('SCF', 'E_CONVERGENCE', 6)
            else:
                psi4.set_local_option('SCF', 'E_CONVERGENCE', 10)
        if not psi4.has_option_changed('SCF', 'D_CONVERGENCE'):
            if procedures['energy'][lowername] in [run_scf, run_dft]:
                psi4.set_local_option('SCF', 'D_CONVERGENCE', 6)
            else:
                psi4.set_local_option('SCF', 'D_CONVERGENCE', 10)

        # Set post-scf convergence criteria (global will cover all correlated modules)
        if not psi4.has_global_option_changed('E_CONVERGENCE'):
            if procedures['energy'][lowername] not in [run_scf, run_dft]:
                psi4.set_global_option('E_CONVERGENCE', 8)

        wfn = procedures['property'][lowername](lowername, **kwargs)

    except KeyError:
        alternatives = ''
        alt_lowername = p4util.text.find_approximate_string_matches(lowername, procedures['property'].keys(), 2)
        if len(alt_lowername) > 0:
            alternatives = """ Did you mean? %s""" % (' '.join(alt_lowername))
        raise ValidationError("""Property method %s not available.%s""" % (lowername, alternatives))

    optstash.restore()

    if return_wfn:
        return (psi4.get_variable('CURRENT ENERGY'), wfn)
    else:
        return psi4.get_variable('CURRENT ENERGY')


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
    | hf                      | Hartree--Fock (HF)  :ref:`[manual] <sec:scf>`                                         |
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
    | df-omp2                 | density-fitted orbital-optimized MP2 :ref:`[manual] <sec:dfocc>`                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-omp3                 | density-fitted orbital-optimized MP3 :ref:`[manual] <sec:dfocc>`                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-omp2.5               | density-fitted orbital-optimized MP2.5 :ref:`[manual] <sec:dfocc>`                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-ccsd                 | density-fitted CCSD (DF-CCSD) :ref:`[manual] <sec:dfocc>`                             |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-ccd                  | density-fitted CCD (DF-CCD) :ref:`[manual] <sec:dfocc>`                               |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-mp3                  | density-fitted MP3 from DFOCC module :ref:`[manual] <sec:dfocc>`                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | df-mp2.5                | density-fitted MP2.5 from DFOCC module :ref:`[manual] <sec:dfocc>`                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | efp                     | efp-only optimizations                                                                |
    +-------------------------+---------------------------------------------------------------------------------------+

    .. _`table:grad_scf`:


    .. include:: autodoc_dft_opt.rst

    .. include:: cfour_table_grad.rst

    .. warning:: Optimizations where the molecule is specified in Z-matrix format
       with dummy atoms will result in the geometry being converted to a Cartesian representation.

    :type name: string
    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

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
        (``'sow'``/``'reap'``). For the latter, run an initial job with
        ``'sow'`` and follow instructions in its output file. For maximum
        flexibility, ``return_wfn`` is always on in ``'reap'`` mode.

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
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)

    full_hess_every = psi4.get_option('OPTKING', 'FULL_HESS_EVERY')
    steps_since_last_hessian = 0
    hessian_with_method = kwargs.get('hessian_with', name)

    # are we in sow/reap mode?
    opt_mode = kwargs.get('mode', 'continuous').lower()
    if opt_mode not in ['continuous', 'sow', 'reap']:
        raise ValidationError("""Optimize execution mode '%s' not valid.""" % (opt_mode))

    optstash = p4util.OptionsState(
        ['OPTKING', 'INTRAFRAG_STEP_LIMIT'],
        ['FINDIF', 'HESSIAN_WRITE'],
        ['OPTKING', 'CART_HESS_READ'],
        ['SCF', 'GUESS_PERSIST'],  # handle on behalf of cbs()
        ['SCF', 'GUESS'])

    n = kwargs.get('opt_iter', 1)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # Shifting the geometry so need to copy the active molecule
    moleculeclone = molecule.clone()

    initial_sym = moleculeclone.schoenflies_symbol()
    while n <= psi4.get_option('OPTKING', 'GEOM_MAXITER'):
        current_sym = moleculeclone.schoenflies_symbol()
        if initial_sym != current_sym:
            raise ValidationError("""Point group changed! You should restart """
                                  """using the last geometry in the output, after """
                                  """carefully making sure all symmetry-dependent """
                                  """input, such as DOCC, is correct.""")
        kwargs['opt_iter'] = n

        # Use orbitals from previous iteration as a guess
        #   set within loop so that can be influenced by fns to optimize (e.g., cbs)
        if (n > 1) and (opt_mode == 'continuous') and (not psi4.get_option('SCF', 'GUESS_PERSIST')):
            psi4.set_local_option('SCF', 'GUESS', 'READ')

        # Before computing gradient, save previous molecule and wavefunction if this is an IRC optimization
        if (n > 1) and (psi4.get_option('OPTKING', 'OPT_TYPE') == 'IRC'):
            old_thisenergy = psi4.get_variable('CURRENT ENERGY')

        # Compute the gradient
        G, wfn = gradient(name, return_wfn=True, molecule=moleculeclone, **kwargs)
        thisenergy = psi4.get_variable('CURRENT ENERGY')

        # above, used to be getting energy as last of energy list from gradient()
        # thisenergy below should ultimately be testing on wfn.energy()

        # S/R: Quit after getting new displacements or if forming gradient fails
        if opt_mode  == 'sow':
            return (0.0, None)
        elif opt_mode == 'reap' and thisenergy == 0.0:
            return (0.0, None)

        psi4.set_gradient(G)

        # S/R: Move opt data file from last pass into namespace for this pass
        if opt_mode == 'reap' and n != 0:
            psi4.IOManager.shared_object().set_specific_retention(1, True)
            psi4.IOManager.shared_object().set_specific_path(1, './')
            if 'opt_datafile' in kwargs:
                restartfile = kwargs.pop('opt_datafile')
                #if psi4.me() == 0:  TODO ask Ryan
                shutil.copy(restartfile, p4util.get_psifile(1))

        opt_func = kwargs.get('opt_func', kwargs.get('func', energy))
        if opt_func.__name__ == 'complete_basis_set':
            psi4.IOManager.shared_object().set_specific_retention(1, True)

        if full_hess_every > -1:
            psi4.set_global_option('HESSIAN_WRITE', True)

        # compute Hessian as requested; frequency wipes out gradient so stash it
        if ((full_hess_every > -1) and (n == 1)) or (steps_since_last_hessian + 1 == full_hess_every):
            G = psi4.get_gradient()  # TODO
            psi4.IOManager.shared_object().set_specific_retention(1, True)
            psi4.IOManager.shared_object().set_specific_path(1, './')
            frequencies(hessian_with_method, **kwargs)
            steps_since_last_hessian = 0
            psi4.set_gradient(G)
            psi4.set_global_option('CART_HESS_READ', True)
        elif (full_hess_every == -1) and psi4.get_global_option('CART_HESS_READ') and (n == 1):
            pass
            # Do nothing; user said to read existing hessian once
        else:
            psi4.set_global_option('CART_HESS_READ', False)
            steps_since_last_hessian += 1

        # Take step. communicate to/from/within optking through legacy_molecule
        psi4.set_legacy_molecule(moleculeclone)
        optking_rval = psi4.optking()
        moleculeclone = psi4.get_legacy_molecule()
        moleculeclone.update_geometry()
        if optking_rval == psi4.PsiReturnType.EndLoop:
            # if this is the end of an IRC run, set wfn, energy, and molecule to that
            # of the last optimized IRC point
            if psi4.get_option('OPTKING', 'OPT_TYPE') == 'IRC':
                thisenergy = old_thisenergy
            print('Optimizer: Optimization complete!')
            psi4.print_out('\n    Final optimized geometry and variables:\n')
            moleculeclone.print_in_input_format()
            # Check if user wants to see the intcos; if so, don't delete them.
            if psi4.get_option('OPTKING', 'INTCOS_GENERATE_EXIT') == False:
                if psi4.get_option('OPTKING', 'KEEP_INTCOS') == False:
                    psi4.opt_clean()
            # Changing environment to optimized geometry as expected by user
            molecule.set_geometry(moleculeclone.geometry())
            for postcallback in hooks['optimize']['post']:
                postcallback(lowername, **kwargs)
            psi4.clean()

            # S/R: Clean up opt input file
            if opt_mode == 'reap':
                with open('OPT-master.in', 'wb') as fmaster:
                    fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n'.encode('utf-8'))
                    fmaster.write('# Optimization complete!\n\n'.encode('utf-8'))

            if opt_func.__name__ == 'complete_basis_set':
                psi4.IOManager.shared_object().set_specific_retention(1, False)

            optstash.restore()

            if return_wfn:
                return (thisenergy, wfn)
            else:
                return thisenergy

        elif optking_rval == psi4.PsiReturnType.Failure:
            print('Optimizer: Optimization failed!')
            if (psi4.get_option('OPTKING', 'KEEP_INTCOS') == False):
                psi4.opt_clean()
            molecule.set_geometry(moleculeclone.geometry())
            psi4.clean()
            optstash.restore()
            return thisenergy

        psi4.print_out('\n    Structure for next step:\n')
        moleculeclone.print_in_input_format()

        # S/R: Preserve opt data file for next pass and switch modes to get new displacements
        if opt_mode == 'reap':
            kwargs['opt_datafile'] = p4util.get_psifile(1)
            kwargs['mode'] = 'sow'

        n += 1

    psi4.print_out('\tOptimizer: Did not converge!')
    if psi4.get_option('OPTKING', 'INTCOS_GENERATE_EXIT') == False:
        if psi4.get_option('OPTKING', 'KEEP_INTCOS') == False:
            psi4.opt_clean()

    optstash.restore()


def parse_arbitrary_order(name):
    r"""Function to parse name string into a method family like CI or MRCC and specific
    level information like 4 for CISDTQ or MRCCSDTQ.

    """
    namelower = name.lower()

    # matches 'mrccsdt(q)'
    if namelower.startswith('mrcc'):

        # avoid undoing fn's good work when called twice
        if namelower == 'mrcc':
            return namelower, None

        # grabs 'sdt(q)'
        ccfullname = namelower[4:]

        # A negative order indicates perturbative method
        methods = {
            'sd'          : { 'method': 1, 'order':  2, 'fullname': 'CCSD'         },
            'sdt'         : { 'method': 1, 'order':  3, 'fullname': 'CCSDT'        },
            'sdtq'        : { 'method': 1, 'order':  4, 'fullname': 'CCSDTQ'       },
            'sdtqp'       : { 'method': 1, 'order':  5, 'fullname': 'CCSDTQP'      },
            'sdtqph'      : { 'method': 1, 'order':  6, 'fullname': 'CCSDTQPH'     },
            'sd(t)'       : { 'method': 3, 'order': -3, 'fullname': 'CCSD(T)'      },
            'sdt(q)'      : { 'method': 3, 'order': -4, 'fullname': 'CCSDT(Q)'     },
            'sdtq(p)'     : { 'method': 3, 'order': -5, 'fullname': 'CCSDTQ(P)'    },
            'sdtqp(h)'    : { 'method': 3, 'order': -6, 'fullname': 'CCSDTQP(H)'   },
            'sd(t)_l'     : { 'method': 4, 'order': -3, 'fullname': 'CCSD(T)_L'    },
            'sdt(q)_l'    : { 'method': 4, 'order': -4, 'fullname': 'CCSDT(Q)_L'   },
            'sdtq(p)_l'   : { 'method': 4, 'order': -5, 'fullname': 'CCSDTQ(P)_L'  },
            'sdtqp(h)_l'  : { 'method': 4, 'order': -6, 'fullname': 'CCSDTQP(H)_L' },
            'sdt-1a'      : { 'method': 5, 'order':  3, 'fullname': 'CCSDT-1a'     },
            'sdtq-1a'     : { 'method': 5, 'order':  4, 'fullname': 'CCSDTQ-1a'    },
            'sdtqp-1a'    : { 'method': 5, 'order':  5, 'fullname': 'CCSDTQP-1a'   },
            'sdtqph-1a'   : { 'method': 5, 'order':  6, 'fullname': 'CCSDTQPH-1a'  },
            'sdt-1b'      : { 'method': 6, 'order':  3, 'fullname': 'CCSDT-1b'     },
            'sdtq-1b'     : { 'method': 6, 'order':  4, 'fullname': 'CCSDTQ-1b'    },
            'sdtqp-1b'    : { 'method': 6, 'order':  5, 'fullname': 'CCSDTQP-1b'   },
            'sdtqph-1b'   : { 'method': 6, 'order':  6, 'fullname': 'CCSDTQPH-1b'  },
            '2'           : { 'method': 7, 'order':  2, 'fullname': 'CC2'          },
            '3'           : { 'method': 7, 'order':  3, 'fullname': 'CC3'          },
            '4'           : { 'method': 7, 'order':  4, 'fullname': 'CC4'          },
            '5'           : { 'method': 7, 'order':  5, 'fullname': 'CC5'          },
            '6'           : { 'method': 7, 'order':  6, 'fullname': 'CC6'          },
            'sdt-3'       : { 'method': 8, 'order':  3, 'fullname': 'CCSDT-3'      },
            'sdtq-3'      : { 'method': 8, 'order':  4, 'fullname': 'CCSDTQ-3'     },
            'sdtqp-3'     : { 'method': 8, 'order':  5, 'fullname': 'CCSDTQP-3'    },
            'sdtqph-3'    : { 'method': 8, 'order':  6, 'fullname': 'CCSDTQPH-3'   }
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

        if namestump in ['mp', 'zapt', 'ci']:
            # Let mp2, mp3, mp4 pass through to select functions
            if namestump == 'mp' and namelevel in [2, 3, 4]:
                return namelower, None
            # Otherwise return method and order
            else:
                return namestump, namelevel
        else:
            return namelower, None
    else:
        return namelower, None


def hessian(name, **kwargs):
    r"""Function complementary to :py:func:`~frequency`. Computes force
    constants, deciding analytic, finite difference of gradients, or
    finite difference of energies.

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    psi4.clean_variables()
    dertype = 2

    # Prevent methods that do not have associated energies 
    if lowername in energy_only_methods:
	raise ValidationError("hessian('%s') does not have an associated hessian" % name)

    optstash = p4util.OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['FINDIF', 'HESSIAN_WRITE'],
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
        if procedures['gradient'][lowername].__name__.startswith('select_'):
            try:
                procedures['gradient'][lowername](lowername, probe=True)
            except ManagedMethodError:
                dertype = 0
                func = energy
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
            raise ValidationError("""Derivative level 'dertype' %s not valid for helper function frequency.""" % (freq_dertype))

    # 3. if the user provides a custom function THAT takes precedence
    if ('freq_func' in kwargs) or ('func' in kwargs):
        if 'func' in kwargs:
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
        alternatives = ''
        alt_lowername = p4util.text.find_approximate_string_matches(lowername, procedures['energy'].keys(), 2)
        if len(alt_lowername) > 0:
            alternatives = """ Did you mean? %s""" % (' '.join(alt_lowername))

        raise ValidationError("""Derivative method 'name' %s and derivative level 'dertype' %s are not available.%s"""
            % (lowername, dertype, alternatives))

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # S/R: Mode of operation- whether finite difference freq run in one job or files farmed out
    freq_mode = kwargs.get('mode', 'continuous').lower()
    if freq_mode == 'continuous':
        pass
    elif freq_mode == 'sow':
        if dertype == 2:
            raise ValidationError("""Frequency execution mode 'sow' not valid for analytic Hessian calculation.""")
    elif freq_mode == 'reap':
        freq_linkage = kwargs.get('linkage', None)
        if freq_linkage is None:
            raise ValidationError("""Frequency execution mode 'reap' requires a linkage option.""")
    else:
        raise ValidationError("""Frequency execution mode '%s' not valid.""" % (freq_mode))

    # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
    if not psi4.has_option_changed('SCF', 'E_CONVERGENCE'):
        if procedures['energy'][lowername] in [run_scf, run_dft]:
            psi4.set_local_option('SCF', 'E_CONVERGENCE', 8)
        else:
            psi4.set_local_option('SCF', 'E_CONVERGENCE', 10)
    if not psi4.has_option_changed('SCF', 'D_CONVERGENCE'):
        if procedures['energy'][lowername] in [run_scf, run_dft]:
            psi4.set_local_option('SCF', 'D_CONVERGENCE', 8)
        else:
            psi4.set_local_option('SCF', 'D_CONVERGENCE', 10)

    # Set post-scf convergence criteria (global will cover all correlated modules)
    if not psi4.has_global_option_changed('E_CONVERGENCE'):
        if procedures['energy'][lowername] not in [run_scf, run_dft]:
            psi4.set_global_option('E_CONVERGENCE', 8)

    # Select certain irreps
    irrep = kwargs.get('irrep', -1)
    if irrep == -1:
        pass  # do all irreps
    else:
        irrep = parse_cotton_irreps(irrep, molecule.schoenflies_symbol())
        irrep -= 1  # A1 irrep is externally 1, internally 0

    # Does an analytic procedure exist for the requested method?
    if dertype == 2:
        psi4.print_out("""hessian() will perform analytic frequency computation.\n""")

        # We have the desired method. Do it.
        wfn = procedures['hessian'][lowername](lowername, molecule=molecule, **kwargs)
        optstash.restore()

        # TODO: check that current energy's being set to the right figure when this code is actually used
        psi4.set_variable('CURRENT ENERGY', wfn.energy())

        if return_wfn:
            return (wfn.hessian(), wfn)
        else:
            return wfn.hessian()

    elif dertype == 1:
        psi4.print_out("""hessian() will perform frequency computation by finite difference of analytic gradients.\n""")

        func = procedures['gradient'][lowername]

        # Shifting the geometry so need to copy the active molecule
        moleculeclone = molecule.clone()

        # Obtain list of displacements
        displacements = psi4.fd_geoms_freq_1(moleculeclone, irrep)
        moleculeclone.reinterpret_coordentry(False)
        moleculeclone.fix_orientation(True)

        # Record undisplaced symmetry for projection of diplaced point groups
        psi4.set_parent_symmetry(molecule.schoenflies_symbol())

        ndisp = len(displacements)
        print(""" %d displacements needed.""" % ndisp)
        gradients = []
        energies = []

        # S/R: Write instructions for sow/reap procedure to output file and reap input file
        if freq_mode == 'sow':
            instructionsO = """\n#    The frequency sow/reap procedure has been selected through mode='sow'. In addition\n"""
            instructionsO += """#    to this output file (which contains no quantum chemical calculations), this job\n"""
            instructionsO += """#    has produced a number of input files (FREQ-*.in) for individual components\n"""
            instructionsO += """#    and a single input file (FREQ-master.in) with a frequency(mode='reap') command.\n"""
            instructionsO += """#    These files may look very peculiar since they contain processed and pickled python\n"""
            instructionsO += """#    rather than normal input. Follow the instructions below (repeated in FREQ-master.in)\n"""
            instructionsO += """#    to continue.\n#\n"""
            instructionsO += """#    Alternatively, a single-job execution of the hessian may be accessed through\n"""
            instructionsO += """#    the frequency wrapper option mode='continuous'.\n#\n"""
            psi4.print_out(instructionsO)

            instructionsM = """\n#    Follow the instructions below to carry out this frequency computation.\n#\n"""
            instructionsM += """#    (1)  Run all of the FREQ-*.in input files on any variety of computer architecture.\n"""
            instructionsM += """#       The output file names must be as given below (these are the defaults when executed\n"""
            instructionsM += """#       as `psi4 FREQ-1.in`, etc.).\n#\n"""
            for rgt in range(ndisp):
                pre = 'FREQ-' + str(rgt + 1)
                instructionsM += """#             psi4 -i %-27s -o %-27s\n""" % (pre + '.in', pre + '.out')
            instructionsM += """#\n#    (2)  Gather all the resulting output files in a directory. Place input file\n"""
            instructionsM += """#         FREQ-master.in into that directory and run it. The job will be minimal in\n"""
            instructionsM += """#         length and give summary results for the frequency computation in its output file.\n#\n"""
            instructionsM += """#             psi4 -i %-27s -o %-27s\n#\n\n""" % ('FREQ-master.in', 'FREQ-master.out')

            with open('FREQ-master.in', 'wb') as fmaster:
                fmaster.write('# This is a psi4 input file auto-generated from the hessian() wrapper.\n\n'.encode('utf-8'))
                fmaster.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
                fmaster.write(p4util.format_options_for_input(moleculeclone, **kwargs))
                p4util.format_kwargs_for_input(fmaster, lmode=2, return_wfn=True, **kwargs)
                fmaster.write(("""retE, retwfn = %s('%s', **kwargs)\n\n""" % (frequency.__name__, lowername)).encode('utf-8'))
                fmaster.write(instructionsM.encode('utf-8'))
            psi4.print_out(instructionsM)

        for n, displacement in enumerate(displacements):
            rfile = 'FREQ-%s' % (n + 1)

            # Build string of title banner
            banners = ''
            banners += """psi4.print_out('\\n')\n"""
            banners += """p4util.banner(' Hessian Computation: Gradient Displacement %d ')\n""" % (n + 1)
            banners += """psi4.print_out('\\n')\n\n"""

            if freq_mode == 'continuous':

                # print progress to file and screen
                psi4.print_out('\n')
                p4util.banner('Loading displacement %d of %d' % (n + 1, ndisp))
                print(""" %d""" % (n + 1), end=('\n' if (n + 1 == ndisp) else ''))
                sys.stdout.flush()

                # Load in displacement into the active molecule (xyz coordinates only)
                moleculeclone.set_geometry(displacement)

                # Perform the gradient calculation
                wfn = func(lowername, molecule=moleculeclone, **kwargs)
                gradients.append(wfn.gradient())
                energies.append(psi4.get_variable('CURRENT ENERGY'))

                # clean may be necessary when changing irreps of displacements
                psi4.clean()

            # S/R: Write each displaced geometry to an input file
            elif freq_mode == 'sow':
                moleculeclone.set_geometry(displacement)

                # S/R: Prepare molecule, options, kwargs, function call and energy save
                #      forcexyz in molecule writer S/R enforcement of !reinterpret_coordentry above
                with open('%s.in' % (rfile), 'wb') as freagent:
                    freagent.write('# This is a psi4 input file auto-generated from the hessian() wrapper.\n\n')
                    freagent.write(p4util.format_molecule_for_input(moleculeclone, forcexyz=True).encode('utf-8'))
                    freagent.write(p4util.format_options_for_input(moleculeclone, **kwargs).encode('utf-8'))
                    p4util.format_kwargs_for_input(freagent, **kwargs)
                    freagent.write("""wfn = %s('%s', **kwargs)\n\n""" % (func.__name__, lowername))
                    freagent.write("""psi4.print_out('\\nHESSIAN RESULT: computation %d for item %d """ % (os.getpid(), n + 1))
                    freagent.write("""yields electronic gradient %r\\n' % (p4util.mat2arr(wfn.gradient())))\n\n""")
                    freagent.write("""psi4.print_out('\\nHESSIAN RESULT: computation %d for item %d """ % (os.getpid(), n + 1))
                    freagent.write("""yields electronic energy %20.12f\\n' % (get_variable('CURRENT ENERGY')))\n\n""")

            # S/R: Read energy from each displaced geometry output file and save in energies array
            elif freq_mode == 'reap':
                exec(banners)
                psi4.set_variable('NUCLEAR REPULSION ENERGY', moleculeclone.nuclear_repulsion_energy())
                pygrad = p4util.extract_sowreap_from_output(rfile, 'HESSIAN', n, freq_linkage, True, label='electronic gradient')
                p4mat = psi4.Matrix(moleculeclone.natom(), 3)
                p4mat.set(pygrad)
                p4mat.print_out()
                gradients.append(p4mat)
                energies.append(p4util.extract_sowreap_from_output(rfile, 'HESSIAN', n, freq_linkage, True))

        # S/R: Quit sow after writing files. Initialize skeleton wfn to receive grad for reap
        if freq_mode == 'sow':
            optstash.restore()
            if return_wfn:
                return (None, None)
            else:
                return None
        elif freq_mode == 'reap':
            wfn = psi4.new_wavefunction(molecule, psi4.get_global_option('BASIS'))

        # Assemble Hessian from gradients
        #   Final disp is undisp, so wfn has mol, G, H general to freq calc
        H = psi4.fd_freq_1(molecule, gradients, irrep)  # TODO or moleculeclone?
        wfn.set_hessian(H)
        wfn.set_frequencies(psi4.get_frequencies())

        # The last item in the list is the reference energy, return it
        psi4.set_variable('CURRENT ENERGY', energies[-1])

        psi4.set_parent_symmetry('')
        optstash.restore()

        if return_wfn:
            return (wfn.hessian(), wfn)
        else:
            return wfn.hessian()

    else:
        psi4.print_out("""hessian() will perform frequency computation by finite difference of analytic energies.\n""")

        # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
        optstash.restore()
        if not psi4.has_option_changed('SCF', 'E_CONVERGENCE'):
            if procedures['energy'][lowername] in [run_scf, run_dft]:
                psi4.set_local_option('SCF', 'E_CONVERGENCE', 10)
            else:
                psi4.set_local_option('SCF', 'E_CONVERGENCE', 11)
        if not psi4.has_option_changed('SCF', 'D_CONVERGENCE'):
            if procedures['energy'][lowername] in [run_scf, run_dft]:
                psi4.set_local_option('SCF', 'D_CONVERGENCE', 10)
            else:
                psi4.set_local_option('SCF', 'D_CONVERGENCE', 11)

        # Set post-scf convergence criteria (global will cover all correlated modules)
        if not psi4.has_global_option_changed('E_CONVERGENCE'):
            if procedures['energy'][lowername] not in [run_scf, run_dft]:
                psi4.set_global_option('E_CONVERGENCE', 10)

        # Shifting the geometry so need to copy the active molecule
        moleculeclone = molecule.clone()

        # Obtain list of displacements
        displacements = psi4.fd_geoms_freq_0(moleculeclone, irrep)
        moleculeclone.fix_orientation(True)
        moleculeclone.reinterpret_coordentry(False)

        # Record undisplaced symmetry for projection of diplaced point groups
        psi4.set_parent_symmetry(molecule.schoenflies_symbol())

        ndisp = len(displacements)

        # This version is pretty dependent on the reference geometry being last (as it is now)
        print(' %d displacements needed.' % ndisp)
        energies = []

        # S/R: Write instructions for sow/reap procedure to output file and reap input file
        if freq_mode == 'sow':
            instructionsO = """\n#    The frequency sow/reap procedure has been selected through mode='sow'. In addition\n"""
            instructionsO += """#    to this output file (which contains no quantum chemical calculations), this job\n"""
            instructionsO += """#    has produced a number of input files (FREQ-*.in) for individual components\n"""
            instructionsO += """#    and a single input file (FREQ-master.in) with a frequency(mode='reap') command.\n"""
            instructionsO += """#    These files may look very peculiar since they contain processed and pickled python\n"""
            instructionsO += """#    rather than normal input. Follow the instructions below (repeated in FREQ-master.in)\n"""
            instructionsO += """#    to continue.\n#\n"""
            instructionsO += """#    Alternatively, a single-job execution of the hessian may be accessed through\n"""
            instructionsO += """#    the frequency wrapper option mode='continuous'.\n#\n"""
            psi4.print_out(instructionsO)

            instructionsM = """\n#    Follow the instructions below to carry out this frequency computation.\n#\n"""
            instructionsM += """#    (1)  Run all of the FREQ-*.in input files on any variety of computer architecture.\n"""
            instructionsM += """#       The output file names must be as given below (these are the defaults when executed\n"""
            instructionsM += """#       as `psi4 FREQ-1.in`, etc.).\n#\n"""
            for rgt in range(ndisp):
                pre = 'FREQ-' + str(rgt + 1)
                instructionsM += """#             psi4 -i %-27s -o %-27s\n""" % (pre + '.in', pre + '.out')
            instructionsM += """#\n#    (2)  Gather all the resulting output files in a directory. Place input file\n"""
            instructionsM += """#         FREQ-master.in into that directory and run it. The job will be minimal in\n"""
            instructionsM += """#         length and give summary results for the frequency computation in its output file.\n#\n"""
            instructionsM += """#             psi4 -i %-27s -o %-27s\n#\n\n""" % ('FREQ-master.in', 'FREQ-master.out')

            with open('FREQ-master.in', 'wb') as fmaster:
                fmaster.write('# This is a psi4 input file auto-generated from the hessian() wrapper.\n\n'.encode('utf-8'))
                fmaster.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
                fmaster.write(p4util.format_options_for_input(moleculeclone, **kwargs))
                p4util.format_kwargs_for_input(fmaster, lmode=2, return_wfn=True, **kwargs)
                fmaster.write(("""retE, retwfn = %s('%s', **kwargs)\n\n""" % (frequency.__name__, lowername)).encode('utf-8'))
                fmaster.write(instructionsM.encode('utf-8'))
            psi4.print_out(instructionsM)

        for n, displacement in enumerate(displacements):
            rfile = 'FREQ-%s' % (n + 1)

            # Build string of title banner
            banners = ''
            banners += """psi4.print_out('\\n')\n"""
            banners += """p4util.banner(' Hessian Computation: Energy Displacement %d ')\n""" % (n + 1)
            banners += """psi4.print_out('\\n')\n\n"""

            if freq_mode == 'continuous':

                # print progress to file and screen
                psi4.print_out('\n')
                p4util.banner('Loading displacement %d of %d' % (n + 1, ndisp))
                print(""" %d""" % (n + 1), end=('\n' if (n + 1 == ndisp) else ''))
                sys.stdout.flush()

                # Load in displacement into the active molecule
                moleculeclone.set_geometry(displacement)

                # Perform the energy calculation
                E, wfn = func(lowername, return_wfn=True, molecule=moleculeclone, **kwargs)
                energies.append(psi4.get_variable('CURRENT ENERGY'))

                # clean may be necessary when changing irreps of displacements
                psi4.clean()

            # S/R: Write each displaced geometry to an input file
            elif freq_mode == 'sow':
                moleculeclone.set_geometry(displacement)

                # S/R: Prepare molecule, options, kwargs, function call and energy save
                with open('%s.in' % (rfile), 'wb') as freagent:
                    freagent.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n')
                    freagent.write(p4util.format_molecule_for_input(moleculeclone, forcexyz=True).encode('utf-8'))
                    freagent.write(p4util.format_options_for_input(moleculeclone, **kwargs).encode('utf-8'))
                    p4util.format_kwargs_for_input(freagent, **kwargs)
                    freagent.write("""electronic_energy = %s('%s', **kwargs)\n\n""" % (func.__name__, lowername))
                    freagent.write("""psi4.print_out('\\nHESSIAN RESULT: computation %d for item %d """ % (os.getpid(), n + 1))
                    freagent.write("""yields electronic energy %20.12f\\n' % (electronic_energy))\n\n""")

            # S/R: Read energy from each displaced geometry output file and save in energies array
            elif freq_mode == 'reap':
                exec(banners)
                psi4.set_variable('NUCLEAR REPULSION ENERGY', moleculeclone.nuclear_repulsion_energy())
                energies.append(p4util.extract_sowreap_from_output(rfile, 'HESSIAN', n, freq_linkage, True))

        # S/R: Quit sow after writing files. Initialize skeleton wfn to receive grad for reap
        if freq_mode == 'sow':
            optstash.restore()
            if return_wfn:
                return (None, None)
            else:
                return None
        elif freq_mode == 'reap':
        #    psi4.set_variable('CURRENT ENERGY', energies[-1])
            wfn = psi4.new_wavefunction(molecule, psi4.get_global_option('BASIS'))

        # Assemble Hessian from energies
        H = psi4.fd_freq_0(molecule, energies, irrep)
        wfn.set_hessian(H)
        wfn.set_frequencies(psi4.get_frequencies())

        # The last item in the list is the reference energy, return it
        psi4.set_variable('CURRENT ENERGY', energies[-1])

        psi4.set_parent_symmetry('')
        optstash.restore()

        if return_wfn:
            return (wfn.hessian(), wfn)
        else:
            return wfn.hessian()


def frequency(name, **kwargs):
    r"""Function to compute harmonic vibrational frequencies.

    :aliases: frequencies(), freq()

    :returns: *float* |w--w| Total electronic energy in Hartrees.

    :returns: (*float*, :ref:`Wavefunction<sec:psimod_Wavefunction>`) |w--w| energy and wavefunction when **return_wfn** specified.

    .. note:: Analytic hessians are not available. Frequencies will proceed through
        finite differences according to availability of gradients or energies.

    .. _`table:freq_gen`:

    :type name: string
    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type return_wfn: :ref:`boolean <op_py_boolean>`
    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|

        Indicate to additionally return the :ref:`Wavefunction<sec:psimod_Wavefunction>`
        calculation result as the second element (after *float* energy) of a tuple.
        Arrays of frequencies and the Hessian can be accessed through the wavefunction.

    :type dertype: :ref:`dertype <op_py_dertype>`
    :param dertype: |dl| ``'hessian'`` |dr| || ``'gradient'`` || ``'energy'``

        Indicates whether analytic (if available- they're not), finite
        difference of gradients (if available) or finite difference of
        energies is to be performed.

    :type mode: string
    :param mode: |dl| ``'continuous'`` |dr| || ``'sow'`` || ``'reap'``

        For a finite difference of energies or gradients frequency, indicates
        whether the calculations required to complete the frequency are to be run
        in one file (``'continuous'``) or are to be farmed out in an
        embarrassingly parallel fashion (``'sow'``/``'reap'``)/ For the latter,
        run an initial job with ``'sow'`` and follow instructions in its output file.
        For maximum flexibility, ``return_wfn`` is always on in ``'reap'`` mode.

    :type irrep: int or string
    :param irrep: |dl| ``-1`` |dr| || ``1`` || ``'b2'`` || ``'App'`` || etc.

        Indicates which symmetry block (:ref:`Cotton <table:irrepOrdering>` ordering) of vibrational
        frequencies to be computed. ``1``, ``'1'``, or ``'a1'`` represents
        :math:`a_1`, requesting only the totally symmetric modes.
        ``-1`` indicates a full frequency calculation.

    :examples:

    >>> # [1] Frequency calculation for all modes through highest available derivatives
    >>> frequency('ccsd')

    >>> # [2] Frequency calculation for b2 modes through finite difference of gradients
    >>> #     printing lowest mode frequency to screen and Hessian to output
    >>> E, wfn = frequencies('scf', dertype=1, irrep=4, return_wfn=True)
    >>> print wfn.frequencies().get(0, 0)
    >>> wfn.hessian().print_out()

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)

    # are we in sow/reap mode?
    freq_mode = kwargs.get('mode', 'continuous').lower()
    if freq_mode not in ['continuous', 'sow', 'reap']:
        raise ValidationError("""Frequency execution mode '%s' not valid.""" % (freq_mode))

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # Compute the hessian
    H, wfn = hessian(name, return_wfn=True, molecule=molecule, **kwargs)

    # S/R: Quit after getting new displacements
    if freq_mode == 'sow':
        return 0.0

    wfn.frequencies().print_out()
    psi4.thermo(wfn, wfn.frequencies())

    for postcallback in hooks['frequency']['post']:
        postcallback(lowername, **kwargs)

    if return_wfn:
        return (psi4.get_variable('CURRENT ENERGY'), wfn)
    else:
        return psi4.get_variable('CURRENT ENERGY')


def gdma(wfn, datafile=""):
    """Function to write wavefunction information in *wfn* to *filename* in
    molden format.

    .. versionadded:: 0.5
       *wfn* parameter passed explicitly

    :returns: None

    :type datafile: string
    :param datafile: optional control file (see GDMA manual) to peform more complicated DMA
                     analyses.  If this option is used, the File keyword must be set to read
                     a filename.fchk, where filename is provided by the WRITER_FILE_LABEL keyword.

    :type wfn: :ref:`Wavefunction<sec:psimod_Wavefunction>`
    :param wfn: set of molecule, basis, orbitals from which to generate DMA analysis

    :examples:

    >>> # [1] DMA analysis from MP2 wavefunction.  N.B. gradient must be requested to generate MP2 density.
    >>> grad, wfn = gradient('mp2', return_wfn=True)
    >>> gdma(wfn)

    """

    # Start by writing a G* checkpoint file, for the GDMA code to read in
    fw = psi4.FCHKWriter(wfn)
    molname = wfn.molecule().name()
    prefix = psi4.get_writer_file_prefix(molname)
    fchkfile = prefix + '.fchk'
    fw.write(fchkfile)

    if datafile:
        commands = datafile
    else:
        densname = wfn.name()
        if densname == "DFT":
            densname = "SCF"
        commands = 'psi4_dma_datafile.dma'
        radii = psi4.get_option('GDMA', 'GDMA_RADIUS')
        origin = psi4.get_option('GDMA', 'GDMA_ORIGIN')
        with open(commands, 'w') as f:
            f.write("File %s Density %s\n" % (fchkfile, densname))
            f.write("Angstrom\n")
            f.write("%s\n" % psi4.get_option('GDMA', 'GDMA_MULTIPOLE_UNITS'))
            f.write("Multipoles\n")
            if origin:
                try:
                    f.write("Origin %f %f %f\n" % (float(origin[0]), float(origin[1]), float(origin[2])))
                except:
                    raise ValidationError("The GDMA origin array should contain three entries: x, y, and z.")
            f.write("Switch %f\n" % psi4.get_option('GDMA', 'GDMA_SWITCH'))
            if radii:
                f.write("Radius %s\n" % " ".join([str(r) for r in radii]))
            f.write("Limit %d\n" % psi4.get_option('GDMA', 'GDMA_LIMIT') )
            f.write("Start\n")
            f.write("Finish\n")
    psi4.run_gdma(wfn, commands)

    os.remove(fchkfile)
    # If we generated the DMA control file, we should clean up here
    if not datafile:
        os.remove(commands)


def molden(wfn, filename):
    """Function to write wavefunction information in *wfn* to *filename* in
    molden format.

    .. versionadded:: 0.5
       *wfn* parameter passed explicitly

    :returns: None

    :type filename: string
    :param filename: destination file name for MOLDEN file

    :type wfn: :ref:`Wavefunction<sec:psimod_Wavefunction>`
    :param wfn: set of molecule, basis, orbitals from which to generate cube files

    :examples:

    >>> # [1] Molden file for DFT calculation
    >>> E, wfn = energy('b3lyp', return_wfn=True)
    >>> molden(wfn, 'mycalc.molden')

    """
    try:
        occa = wfn.occupation_a()
        occb = wfn.occupation_a()
    except AttributeError:
        psi4.print_out("\n!Molden warning: This wavefunction does not have occupation numbers.\n"
                       "Writing zero's for occupation numbers\n\n")
        occa = psi4.Vector(wfn.nmopi())
        occb = psi4.Vector(wfn.nmopi())

    # At this point occupation number will be difficult to build, lets set them to zero
    mw = psi4.MoldenWriter(wfn)
    mw.write(filename, wfn.Ca(), wfn.Cb(), wfn.epsilon_a(), wfn.epsilon_b(), occa, occb)

def parse_cotton_irreps(irrep, point_group):
    r"""Function to return validated Cotton ordering index for molecular
    *point_group* from string or integer irreducible representation *irrep*.

    """
    cotton = {
        'c1': {
            'a': 1,
            '1': 1
        },
        'ci': {
            'ag': 1,
            'au': 2,
            '1': 1,
            '2': 2
        },
        'c2': {
            'a': 1,
            'b': 2,
            '1': 1,
            '2': 2
        },
        'cs': {
            'ap': 1,
            'app': 2,
            '1': 1,
            '2': 2
        },
        'd2': {
            'a': 1,
            'b1': 2,
            'b2': 3,
            'b3': 4,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4
        },
        'c2v': {
            'a1': 1,
            'a2': 2,
            'b1': 3,
            'b2': 4,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4
        },
        'c2h': {
            'ag': 1,
            'bg': 2,
            'au': 3,
            'bu': 4,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4,
        },
        'd2h': {
            'ag': 1,
            'b1g': 2,
            'b2g': 3,
            'b3g': 4,
            'au': 5,
            'b1u': 6,
            'b2u': 7,
            'b3u': 8,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4,
            '5': 5,
            '6': 6,
            '7': 7,
            '8': 8
        }
    }

    try:
        return cotton[point_group.lower()][str(irrep).lower()]
    except KeyError:
        raise ValidationError("""Irrep '%s' not valid for point group '%s'.""" % (str(irrep), point_group))


# Aliases
opt = optimize
freq = frequency
frequencies = frequency
prop = property
