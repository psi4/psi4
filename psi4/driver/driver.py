#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""
Nexus of psi4.driver module with primary user-facing functions, including
single-point energies, geometry optimizations, properties, and vibrational
frequency calculations.

"""
import json
import os
import re
import copy
import shutil
import sys
import logging
from typing import Dict, Optional, Union
import logging

import numpy as np

from psi4 import core  # for typing
from psi4.driver import driver_util
from psi4.driver import driver_cbs
from psi4.driver import driver_nbody
from psi4.driver import driver_findif
from psi4.driver import task_planner
from psi4.driver import p4util
from psi4.driver import qcdb
from psi4.driver import pp, nppp, nppp10
from psi4.driver.p4util.exceptions import *
from psi4.driver.procrouting import *
from psi4.driver.mdi_engine import mdi_run
from psi4.driver.task_base import AtomicComputer

# never import wrappers or aliases into this file

logger = logging.getLogger(__name__)


def _energy_is_invariant(gradient_rms, stationary_criterion=1.e-2):
    """Polls options and probes `gradient` to return whether current method
    and system expected to be invariant to translations and rotations of
    the coordinate system.

    """
    stationary_point = gradient_rms < stationary_criterion  # 1.e-2 pulled out of a hat

    mol = core.get_active_molecule()
    efp_present = hasattr(mol, 'EFP')

    translations_projection_sound = (not core.get_option('SCF', 'EXTERN') and not core.get_option('SCF', 'PERTURB_H')
                                     and not efp_present)
    rotations_projection_sound = (translations_projection_sound and stationary_point)

    return translations_projection_sound, rotations_projection_sound


def _filter_renamed_methods(compute, method):
    r"""Raises UpgradeHelper when a method has been renamed."""
    if method == "dcft":
        raise UpgradeHelper(compute + "('dcft')", compute + "('dct')", 1.4,
                            " All instances of 'dcft' should be replaced with 'dct'.")


def energy(name, **kwargs):
    r"""Function to compute the single-point electronic energy.

    :returns: *float* |w--w| Total electronic energy in Hartrees. SAPT & EFP return interaction energy.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| energy and wavefunction when **return_wfn** specified.

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`CURRENT ENERGY`
       * :psivar:`CURRENT REFERENCE ENERGY`
       * :psivar:`CURRENT CORRELATION ENERGY`

    :type name: str
    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type return_wfn: :ref:`boolean <op_py_boolean>`
    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|

        Indicate to additionally return the :py:class:`~psi4.core.Wavefunction`
        calculation result as the second element (after *float* energy) of a tuple.

    :type write_orbitals: str, :ref:`boolean <op_py_boolean>`
    :param write_orbitals: ``filename`` || |dl| ``'on'`` |dr| || ``'off'`` 

        (str) Save wfn containing current orbitals to the given file name after each SCF iteration
        and retain after |PSIfour| finishes.

        (:ref:`boolean <op_py_boolean>`) Turns writing the orbitals after the converged SCF on/off.
        Orbital file will be deleted unless |PSIfour| is called with `-m` flag.

    :type restart_file: str
    :param restart_file: ``['file.1, file.32]`` || ``./file`` || etc.

        Existing files to be renamed and copied for calculation restart, e.g. a serialized wfn or module-specific binary data.

    .. _`table:energy_gen`:

    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | name                    | calls method                                                                                                                          |
    +=========================+=======================================================================================================================================+
    | efp                     | (with LibEFP) effective fragment potential (EFP) :ref:`[manual] <sec:libefp>`                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT) :ref:`[manual] <sec:scf>` :ref:`[details] <dd_b3lyp>`                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | hf                      | HF self consistent field (SCF) :ref:`[manual] <sec:scf>` :ref:`[details] <dd_hf>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | qchf                    | quadratically-convergent HF                                                                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | hf3c                    | HF with dispersion, BSSE, and basis set corrections :ref:`[manual] <sec:gcp>`                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | pbeh3c                  | PBEh with dispersion, BSSE, and basis set corrections :ref:`[manual] <sec:gcp>`                                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | dct                     | density cumulant (functional) theory :ref:`[manual] <sec:dct>`                                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | mp2                     | 2nd-order |MollerPlesset| perturbation theory (MP2) :ref:`[manual] <sec:dfmp2>` :ref:`[details] <dd_mp2>`                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | scs-mp2                 | spin-component scaled MP2 :ref:`[manual] <sec:occ_nonoo>`                                                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | scs(n)-mp2              | a special version of SCS-MP2 for nucleobase interactions :ref:`[manual] <sec:occ_nonoo>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | scs-mp2-vdw             | a special version of SCS-MP2 (from ethene dimers) :ref:`[manual] <sec:occ_nonoo>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sos-mp2                 | spin-opposite scaled MP2 :ref:`[manual] <sec:occ_nonoo>`                                                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | dlpno-mp2               | local MP2 with pair natural orbital domains (DLPNO) :ref:`[manual] <sec:dlpnomp2>`                                                    |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | scs-dlpno-mp2           | spin-component-scaled DLPNO MP2 :ref:`[manual] <sec:dlpnomp2>`                                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | mp3                     | 3rd-order |MollerPlesset| perturbation theory (MP3) :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <dd_mp3>`                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-mp3                 | MP3 with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | scs-mp3                 | spin-component scaled MP3 :ref:`[manual] <sec:occ_nonoo>`                                                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sos-pi-mp2              | A special version of SOS-MP2 for pi systems :ref:`[manual] <sec:occ_nonoo>`                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | mp2.5                   | average of MP2 and MP3 :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <dd_mp2p5>`                                                    |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | mp4(sdq)                | 4th-order MP perturbation theory (MP4) less triples :ref:`[manual] <sec:fnompn>` :ref:`[details] <dd_mp4_prsdq_pr>`                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-mp4(sdq)            | MP4 (less triples) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | mp4                     | full MP4 :ref:`[manual] <sec:fnompn>` :ref:`[details] <dd_mp4>`                                                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-mp4                 | full MP4 with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | mp\ *n*                 | *n*\ th-order |MollerPlesset| (MP) perturbation theory :ref:`[manual] <sec:arbpt>` :ref:`[details] <dd_mp4>`                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | zapt\ *n*               | *n*\ th-order z-averaged perturbation theory (ZAPT) :ref:`[manual] <sec:arbpt>` :ref:`[details] <dd_zapt2>`                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | omp2                    | orbital-optimized second-order MP perturbation theory :ref:`[manual] <sec:occ_oo>` :ref:`[details] <dd_omp2>`                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | scs-omp2                | spin-component scaled OMP2 :ref:`[manual] <sec:occ_oo>`                                                                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sos-omp2                | spin-opposite scaled OMP2 :ref:`[manual] <sec:occ_oo>`                                                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | omp3                    | orbital-optimized third-order MP perturbation theory :ref:`[manual] <sec:occ_oo>` :ref:`[details] <dd_omp3>`                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | scs-omp3                | spin-component scaled OMP3 :ref:`[manual] <sec:occ_oo>`                                                                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sos-omp3                | spin-opposite scaled OMP3 :ref:`[manual] <sec:occ_oo>`                                                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | omp2.5                  | orbital-optimized MP2.5 :ref:`[manual] <sec:occ_oo>` :ref:`[details] <dd_omp2p5>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | lccsd, cepa(0)          | coupled electron pair approximation variant 0 :ref:`[manual] <sec:fnocepa>` :ref:`[details] <dd_lccsd>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-lccsd, fno-cepa(0)  | CEPA(0) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | cepa(1)                 | coupled electron pair approximation variant 1 :ref:`[manual] <sec:fnocepa>` :ref:`[details] <dd_cepa_pr1_pr>`                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-cepa(1)             | CEPA(1) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | cepa(3)                 | coupled electron pair approximation variant 3 :ref:`[manual] <sec:fnocepa>` :ref:`[details] <dd_cepa_pr3_pr>`                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-cepa(3)             | CEPA(3) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | acpf                    | averaged coupled-pair functional :ref:`[manual] <sec:fnocepa>` :ref:`[details] <dd_acpf>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-acpf                | ACPF with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | aqcc                    | averaged quadratic coupled cluster :ref:`[manual] <sec:fnocepa>` :ref:`[details] <dd_aqcc>`                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-aqcc                | AQCC with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | qcisd                   | quadratic CI singles doubles (QCISD) :ref:`[manual] <sec:fnocc>` :ref:`[details] <dd_qcisd>`                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-qcisd               | QCISD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | remp2                   | 2nd-order retaining-the-excitation-degree MP hybrid perturbation theory :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <dd_remp2>`   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | oremp2                  | orbital-optimized REMP2 :ref:`[manual] <sec:occ_oo>` :ref:`[details] <dd_oremp2>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | lccd                    | Linear CCD :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <dd_lccd>`                                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-lccd                | LCCD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | olccd                   | orbital optimized LCCD :ref:`[manual] <sec:occ_oo>` :ref:`[details] <dd_olccd>`                                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | cc2                     | approximate coupled cluster singles and doubles (CC2) :ref:`[manual] <sec:cc>` :ref:`[details] <dd_cc2>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | ccd                     | coupled cluster doubles (CCD) :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <dd_ccd>`                                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD) :ref:`[manual] <sec:cc>` :ref:`[details] <dd_ccsd>`                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | bccd                    | Brueckner coupled cluster doubles (BCCD) :ref:`[manual] <sec:cc>` :ref:`[details] <dd_bccd>`                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-ccsd                | CCSD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | qcisd(t)                | QCISD with perturbative triples :ref:`[manual] <sec:fnocc>` :ref:`[details] <dd_qcisd_prt_pr>`                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-qcisd(t)            | QCISD(T) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | ccsd(t)                 | CCSD with perturbative triples (CCSD(T)) :ref:`[manual] <sec:cc>` :ref:`[details] <dd_ccsd_prt_pr>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | a-ccsd(t)               | CCSD with asymmetric perturbative triples (A-CCSD(T)) :ref:`[manual] <sec:cc>` :ref:`[details] <dd_accsd_prt_pr>`                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | bccd(t)                 | BCCD with perturbative triples :ref:`[manual] <sec:cc>` :ref:`[details] <dd_bccd_prt_pr>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-ccsd(t)             | CCSD(T) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | cc3                     | approximate CC singles, doubles, and triples (CC3) :ref:`[manual] <sec:cc>` :ref:`[details] <dd_cc3>`                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | ccenergy                | **expert** full control over ccenergy module                                                                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | cisd                    | configuration interaction (CI) singles and doubles (CISD) :ref:`[manual] <sec:ci>` :ref:`[details] <dd_cisd>`                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fno-cisd                | CISD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | cisdt                   | CI singles, doubles, and triples (CISDT) :ref:`[manual] <sec:ci>`                                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | cisdtq                  | CI singles, doubles, triples, and quadruples (CISDTQ) :ref:`[manual] <sec:ci>`                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | ci\ *n*                 | *n*\ th-order CI :ref:`[manual] <sec:ci>` :ref:`[details] <dd_cisd>`                                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fci                     | full configuration interaction (FCI) :ref:`[manual] <sec:ci>` :ref:`[details] <dd_fci>`                                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | detci                   | **expert** full control over detci module                                                                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | gaussian-2, g2          | Gaussian-2 composite method :ref:`[manual] <sec:fnogn>`                                                                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | casscf                  | complete active space self consistent field (CASSCF) :ref:`[manual] <sec:ci>`                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | rasscf                  | restricted active space self consistent field (RASSCF) :ref:`[manual] <sec:ci>`                                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | mcscf                   | multiconfigurational self consistent field (SCF) :ref:`[manual] <sec:psimrcc>`                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | psimrcc                 | Mukherjee multireference coupled cluster (Mk-MRCC) :ref:`[manual] <sec:psimrcc>`                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | psimrcc_scf             | Mk-MRCC with regular SCF module (convenience function) :ref:`[manual] <sec:psimrcc>`                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | dmrg-scf                | (with CheMPS2) density matrix renormalization group SCF :ref:`[manual] <sec:chemps2>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | dmrg-caspt2             | (with CheMPS2) density matrix renormalization group CASPT2 :ref:`[manual] <sec:chemps2>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | dmrg-ci                 | (with CheMPS2) density matrix renormalization group CI :ref:`[manual] <sec:chemps2>`                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt0                   | 0th-order symmetry adapted perturbation theory (SAPT) :ref:`[manual] <sec:sapt>`                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | ssapt0                  | 0th-order SAPT with special exchange scaling :ref:`[manual] <sec:sapt>`                                                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | fisapt0                 | 0th-order functional and/or intramolecular SAPT :ref:`[manual] <sec:fisapt>`                                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sf-sapt                 | 0th-order spin-flip SAPT :ref:`[manual] <sec:sfsapt>`                                                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt(dft)               | 0th-order SAPT upon KS reference :ref:`[manual] <sec:saptdft>`                                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2                   | 2nd-order SAPT, traditional definition :ref:`[manual] <sec:sapt>`                                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+                  | SAPT including all 2nd-order terms :ref:`[manual] <sec:sapt>`                                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)               | SAPT including perturbative triples :ref:`[manual] <sec:sapt>`                                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+3                 | SAPT including all 3rd-order terms :ref:`[manual] <sec:sapt>`                                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+(ccd)             | SAPT2+ with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                                            |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)          | SAPT2+(3) with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+3(ccd)            | SAPT2+3 with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+dmp2              | SAPT including all 2nd-order terms and MP2 correction :ref:`[manual] <sec:sapt>`                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)dmp2           | SAPT including perturbative triples and MP2 correction :ref:`[manual] <sec:sapt>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+3dmp2             | SAPT including all 3rd-order terms and MP2 correction :ref:`[manual] <sec:sapt>`                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+(ccd)dmp2         | SAPT2+ with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)dmp2      | SAPT2+(3) with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+3(ccd)dmp2        | SAPT2+3 with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt0-ct                | 0th-order SAPT plus charge transfer (CT) calculation :ref:`[manual] <sec:saptct>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2-ct                | SAPT2 plus CT :ref:`[manual] <sec:saptct>`                                                                                            |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+-ct               | SAPT2+ plus CT :ref:`[manual] <sec:saptct>`                                                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)-ct            | SAPT2+(3) plus CT :ref:`[manual] <sec:saptct>`                                                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+3-ct              | SAPT2+3 plus CT :ref:`[manual] <sec:saptct>`                                                                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+(ccd)-ct          | SAPT2+(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)-ct       | SAPT2+(3)(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | sapt2+3(ccd)-ct         | SAPT2+3(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | adc                     | 2nd-order algebraic diagrammatic construction (ADC), deprecated :ref:`[manual] <sec:adc>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | adc(1)                  | (with ADCC) 1st-order algebraic diagrammatic construction (ADC) :ref:`[manual] <sec:adc>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | adc(2)                  | (with ADCC) 2nd-order ADC :ref:`[manual] <sec:adc>`                                                                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | adc(2)-x                | (with ADCC) extended 2nd-order ADC :ref:`[manual] <sec:adc>`                                                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | adc(3)                  | (with ADCC) 3rd-order ADC :ref:`[manual] <sec:adc>`                                                                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | cvs-adc(1)              | (with ADCC) core-valence separation (CVS) 1st-order ADC :ref:`[manual] <sec:adc>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | cvs-adc(2)              | (with ADCC) CVS 2nd-order ADC :ref:`[manual] <sec:adc>`                                                                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | cvs-adc(2)-x            | (with ADCC) CVS extended 2nd-order ADC :ref:`[manual] <sec:adc>`                                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | cvs-adc(3)              | (with ADCC) CVS 3rd-order ADC :ref:`[manual] <sec:adc>`                                                                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | ep2                     | 2nd-order electron propagator theory                                                                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | eom-cc2                 | equation of motion (EOM) CC2 :ref:`[manual] <sec:eomcc>`                                                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | eom-ccsd                | EOM-CCSD :ref:`[manual] <sec:eomcc>`                                                                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
    | eom-cc3                 | EOM-CC3 :ref:`[manual] <sec:eomcc>`                                                                                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------------------------------+

    .. comment missing and why
    .. comment custom-*mp for occ --- not yet tested or mirrored in dfocc

    .. include:: /autodoc_dft_energy.rst

    .. include:: /mrcc_table_energy.rst

    .. include:: /cfour_table_energy.rst

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
    >>> set basis cc-pVDZ
    >>> set reference rohf
    >>> scf_e, scf_wfn = energy('scf', return_wfn=True)
    >>> H2.set_multiplicity(3)
    >>> core.MintsHelper(scf_wfn.basisset()).integrals()
    >>> energy('detci', ref_wfn=scf_wfn)

    >>> # [5] Run two CI calculations, keeping the integrals generated in the first one.
    >>> molecule ne {\nNe\n}
    >>> set basis cc-pVDZ
    >>> cisd_e, cisd_wfn = energy('cisd', return_wfn=True)
    >>> energy('fci', ref_wfn=cisd_wfn)

    >>> # [6] Can automatically perform complete basis set extrapolations
    >>> energy("CCSD/cc-pV[DT]Z")

    >>> # [7] Can automatically perform delta corrections that include extrapolations
    >>> # even with a user-defined extrapolation formula. See sample inputs named
    >>> # cbs-xtpl* for more examples of this input style
    >>> energy("MP2/aug-cc-pv([d,t]+d)z + d:ccsd(t)/cc-pvdz", corl_scheme=myxtplfn_2)

    """
    kwargs = p4util.kwargs_lower(kwargs)

    # Bounce to MDI (MolSSI driver interface) if mdi kwarg
    use_mdi = kwargs.pop('mdi', False)
    if use_mdi:
        return mdi_run(name, **kwargs)

    core.print_out("\nScratch directory: %s\n" % core.IOManager.shared_object().get_default_path())

    basisstash = p4util.OptionsState(['BASIS'])
    return_wfn = kwargs.pop('return_wfn', False)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    ## Pre-planning interventions

    # * Trip on function or alias as name
    lowername = driver_util.upgrade_interventions(name)

    _filter_renamed_methods("energy", lowername)

    # * Avert pydantic anger at incomplete modelchem spec
    userbas = core.get_global_option('BASIS') or kwargs.get('basis')
    if lowername in integrated_basis_methods and userbas is None:
        kwargs['basis'] = '(auto)'

    # Are we planning?
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
    logger.debug('ENERGY PLAN')
    logger.debug(pp.pformat(plan.dict()))

    if kwargs.get("return_plan", False):
        # Plan-only requested
        return plan

    elif not isinstance(plan, AtomicComputer):
        # Advanced "Computer" active
        plan.compute()
        return plan.get_psi_results(return_wfn=return_wfn)

    else:
        # We have unpacked to an AtomicInput
        lowername = plan.method
        basis = plan.basis
        core.set_global_option("BASIS", basis)

    ## Second half of this fn -- entry means program running exactly analytic 0th derivative

    # Commit to procedures['energy'] call hereafter
    core.clean_variables()

    #for precallback in hooks['energy']['pre']:
    #    precallback(lowername, **kwargs)

    # needed (+restore below) so long as AtomicComputer-s aren't run through json (where convcrit also lives)
    optstash = driver_util.negotiate_convergence_criterion((0, 0), lowername, return_optstash=True)
    optstash2 = p4util.OptionsState(['SCF', 'GUESS'])

    # Before invoking the procedure, we rename any file that should be read.
    # This is a workaround to do restarts with the current PSI4 capabilities
    # before actual, clean restarts are put in there
    # Restartfile is always converted to a single-element list if
    # it contains a single string
    # DGAS Note: This is hacked together at this point and should be revamped.
    if 'restart_file' in kwargs:
        restartfile = kwargs['restart_file']  # Option still available for procedure-specific action
        if not isinstance(restartfile, (list, tuple)):
            restartfile = (restartfile, )
        # Rename the files to be read to be consistent with psi4's file system
        for item in restartfile:
            is_numpy_file = (os.path.isfile(item) and item.endswith(".npy")) or os.path.isfile(item + ".npy")
            name_split = re.split(r'\.', item)
            if is_numpy_file:
                core.set_local_option('SCF', 'GUESS' ,'READ')
                core.print_out(" Found user provided orbital data. Setting orbital guess to READ")
                fname = os.path.split(os.path.abspath(core.get_writer_file_prefix(molecule.name())))[1]
                psi_scratch = core.IOManager.shared_object().get_default_path()
                file_num = item.split('.')[-2] if "180" in item else "180"
                targetfile = os.path.join(psi_scratch, fname + "." + file_num + ".npy")
                if not item.endswith(".npy"):
                    item = item + ".npy"
            else:
                filenum = name_split[-1]
                try:
                    filenum = int(filenum)
                except ValueError:
                    filenum = 32  # Default file number is the checkpoint one
                psioh = core.IOManager.shared_object()
                psio = core.IO.shared_object()
                filepath = psioh.get_file_path(filenum)
                namespace = psio.get_default_namespace()
                pid = str(os.getpid())
                prefix = 'psi'
                targetfile = filepath + prefix + '.' + pid + '.' + namespace + '.' + str(filenum)
            core.print_out(f" \n Copying restart file <{item}> to <{targetfile}> for internal processing\n")
            shutil.copy(item, targetfile)

    logger.info(f"Compute energy(): method={lowername}, basis={core.get_global_option('BASIS').lower()}, molecule={molecule.name()}, nre={'w/EFP' if hasattr(molecule, 'EFP') else molecule.nuclear_repulsion_energy()}")
    logger.debug("w/EFP" if hasattr(molecule, "EFP") else pp.pformat(molecule.to_dict()))
    wfn = procedures['energy'][lowername](lowername, molecule=molecule, **kwargs)
    logger.info(f"Return energy(): {core.variable('CURRENT ENERGY')}")

    for postcallback in hooks['energy']['post']:
        postcallback(lowername, wfn=wfn, **kwargs)

    basisstash.restore()
    optstash.restore()
    optstash2.restore()

    if return_wfn:  # TODO current energy safer than wfn.energy() for now, but should be revisited

        # TODO place this with the associated call, very awkward to call this in other areas at the moment
        if lowername in ['efp', 'mrcc', 'dmrg']:
            core.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
            core.print_out("The returned wavefunction is the incoming reference wavefunction.\n\n")
        elif 'sapt' in lowername:
            core.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
            core.print_out("The returned wavefunction is the dimer SCF wavefunction.\n\n")

        return (core.variable('CURRENT ENERGY'), wfn)
    else:
        return core.variable('CURRENT ENERGY')


def gradient(name, **kwargs):
    r"""Function complementary to :py:func:`~psi4.driver.optimize()`. Carries out one gradient pass,
    deciding analytic or finite difference.

    :returns: :py:class:`~psi4.core.Matrix` |w--w| Total electronic gradient in Hartrees/Bohr.

    :returns: (:py:class:`~psi4.core.Matrix`, :py:class:`~psi4.core.Wavefunction`) |w--w| gradient and wavefunction when **return_wfn** specified.

    :examples:

    >>> # [1] Single-point dft gradient getting the gradient
    >>> #     in file, core.Matrix, and np.array forms
    >>> set gradient_write on
    >>> G, wfn = gradient('b3lyp-d', return_wfn=True)
    >>> wfn.gradient().print_out()
    >>> np.array(G)

    """
    ## First half of this fn -- entry means user wants a 1st derivative by any means

    kwargs = p4util.kwargs_lower(kwargs)

    core.print_out("\nScratch directory: %s\n" % core.IOManager.shared_object().get_default_path())

    basisstash = p4util.OptionsState(['BASIS'])
    return_wfn = kwargs.pop('return_wfn', False)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    # Convert wrapper directives from options (where ppl know to find them) to kwargs (suitable for non-globals transmitting)
    kwargs['findif_verbose'] = core.get_option("FINDIF", "PRINT")
    kwargs['findif_stencil_size'] = core.get_option("FINDIF", "POINTS")
    kwargs['findif_step_size'] = core.get_option("FINDIF", "DISP_SIZE")

    ## Pre-planning interventions

    # * Trip on function or alias as name
    lowername = driver_util.upgrade_interventions(name)
    _filter_renamed_methods("gradient", lowername)

    # * Prevent methods that do not have associated derivatives
    if lowername in energy_only_methods:
        raise ValidationError(f"`gradient('{name}')` does not have an associated gradient.")

    # * Avert pydantic anger at incomplete modelchem spec
    userbas = core.get_global_option('BASIS') or kwargs.get('basis')
    if lowername in integrated_basis_methods and userbas is None:
        kwargs['basis'] = '(auto)'

    # Are we planning?
    plan = task_planner.task_planner("gradient", lowername, molecule, **kwargs)
    logger.debug('GRADIENT PLAN')
    logger.debug(pp.pformat(plan.dict()))

    if kwargs.get("return_plan", False):
        # Plan-only requested
        return plan

    elif not isinstance(plan, AtomicComputer):
        # Advanced "Computer" active
        plan.compute()
        return plan.get_psi_results(return_wfn=return_wfn)

    else:
        # We have unpacked to an AtomicInput
        lowername = plan.method
        basis = plan.basis
        core.set_global_option("BASIS", basis)

    ## Second half of this fn -- entry means program running exactly analytic 1st derivative

    # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
    optstash = driver_util.negotiate_convergence_criterion((1, 1), lowername, return_optstash=True)

    # Commit to procedures[] call hereafter
    core.clean_variables()

    # no analytic derivatives for scf_type cd
    if core.get_global_option('SCF_TYPE') == 'CD':
        raise ValidationError("""No analytic derivatives for SCF_TYPE CD.""")

    core.print_out("""gradient() will perform analytic gradient computation.\n""")

    # Perform the gradient calculation
    logger.info(f"Compute gradient(): method={lowername}, basis={core.get_global_option('BASIS').lower()}, molecule={molecule.name()}, nre={'w/EFP' if hasattr(molecule, 'EFP') else molecule.nuclear_repulsion_energy()}")
    logger.debug("w/EFP" if hasattr(molecule, "EFP") else pp.pformat(molecule.to_dict()))
    wfn = procedures['gradient'][lowername](lowername, molecule=molecule, **kwargs)
    logger.info(f"Return gradient(): {core.variable('CURRENT ENERGY')}")
    logger.info(nppp(wfn.gradient().np))

    basisstash.restore()
    optstash.restore()

    driver_findif.gradient_write(wfn)

    if return_wfn:
        return (wfn.gradient(), wfn)
    else:
        return wfn.gradient()


def properties(*args, **kwargs):
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
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | mp2                | MP2 with density fitting only (mp2_type df)   | RHF            | Listed :ref:`here <sec:oeprop>`                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | cc2                | 2nd-order approximate CCSD                    | RHF            | dipole, quadrupole, polarizability, rotation, roa_tensor      |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | ccsd               | Coupled cluster singles and doubles (CCSD)    | RHF            | dipole, quadrupole, polarizability, rotation, roa_tensor      |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | dct                | density cumulant (functional) theory          | RHF/UHF        | Listed :ref:`here <sec:oeprop>`                               |
    |                    | :ref:`[manual] <sec:dct>`                     |                |                                                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | omp2               | orbital-optimized second-order                | RHF/UHF        | Listed :ref:`here <sec:oeprop>`                               |
    |                    | MP perturbation theory                        |                | Density fitted only                                           |
    |                    | :ref:`[manual] <sec:occ_oo>`                  |                |                                                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | omp3               | orbital-optimized third-order                 | RHF/UHF        | Listed :ref:`here <sec:oeprop>`                               |
    |                    | MP perturbation theory                        |                | Density fitted only                                           |
    |                    | :ref:`[manual] <sec:occ_oo>`                  |                |                                                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | omp2.5             | orbital-optimized MP2.5                       | RHF/UHF        | Listed :ref:`here <sec:oeprop>`                               |
    |                    | :ref:`[manual] <sec:occ_oo>`                  |                | Density fitted only                                           |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | olccd              | orbital optimized LCCD                        | RHF/UHF        | Listed :ref:`here <sec:oeprop>`                               |
    |                    | :ref:`[manual] <sec:occ_oo>`                  |                | Density fitted only                                           |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | eom-cc2            | 2nd-order approximate EOM-CCSD                | RHF            | oscillator_strength, rotational_strength                      |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | eom-ccsd           | Equation-of-motion CCSD (EOM-CCSD)            | RHF            | oscillator_strength, rotational_strength                      |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | cisd, cisdt,       | Configuration interaction                     | RHF/ROHF       | Listed :ref:`here <sec:oeprop>`, transition_dipole,           |
    | cisdt, cisdtq,     |                                               |                | transition_quadrupole                                         |
    | ci5, ..., fci      |                                               |                |                                                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | casscf, rasscf     | Multi-configurational SCF                     | RHF/ROHF       | Listed :ref:`here <sec:oeprop>`, transition_dipole,           |
    |                    |                                               |                | transition_quadrupole                                         |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | adc(0), adc(1),    | Algebraic-diagrammatic construction methods   | RHF/UHF        | dipole, transition_dipole, oscillator_strength,               |
    | ..., adc(3),       | :ref:`[manual] <sec:adc>`                     |                | rotational_strength                                           |
    | cvs-adc(0), ...    |                                               |                |                                                               |
    | cvs-adc(3)         |                                               |                |                                                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+

    :type name: str
    :param name: ``'ccsd'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type properties: List[str]
    :param properties: |dl| ``[]`` |dr| || ``['rotation', 'polarizability', 'oscillator_strength', 'roa']`` || etc.

        Indicates which properties should be computed. Defaults to dipole and quadrupole.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :examples:

    >>> # [1] Optical rotation calculation
    >>> properties('cc2', properties=['rotation'])

    """
    kwargs = p4util.kwargs_lower(kwargs)

    basisstash = p4util.OptionsState(['BASIS'])
    return_wfn = kwargs.pop('return_wfn', False)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    ## Pre-planning interventions

    # * Trip on function or alias as name
    lowername = driver_util.upgrade_interventions(args[0])

    _filter_renamed_methods("properties", lowername)

    props = kwargs.get('properties', ['dipole', 'quadrupole'])
    if len(args) > 1:
        props += args[1:]
    kwargs['properties'] = p4util.drop_duplicates(props)

    # * Avert pydantic anger at incomplete modelchem spec
    userbas = core.get_global_option('BASIS') or kwargs.get('basis')
    if lowername in integrated_basis_methods and userbas is None:
        kwargs['basis'] = '(auto)'

    # Are we planning?
    plan = task_planner.task_planner("properties", lowername, molecule, **kwargs)
    logger.debug('PROPERTIES PLAN')
    logger.debug(pp.pformat(plan.dict()))

    if kwargs.get("return_plan", False):
        # Plan-only requested
        return plan

    elif not isinstance(plan, AtomicComputer):
        # Advanced "Computer" active
        plan.compute()
        return plan.get_psi_results(return_wfn=return_wfn)

    else:
        # We have unpacked to an AtomicInput
        lowername = plan.method
        basis = plan.basis
        core.set_global_option("BASIS", basis)

    ## Second half of this fn -- entry means program running exactly analytic property

    # Commit to procedures['properties'] call hereafter
    core.clean_variables()

    # needed (+restore below) so long as AtomicComputer-s aren't run through json (where convcrit also lives)
    optstash = driver_util.negotiate_convergence_criterion(("prop", "prop"), lowername, return_optstash=True)

    logger.info(f"Compute properties(): method={lowername}, basis={core.get_global_option('BASIS').lower()}, molecule={molecule.name()}, nre={'w/EFP' if hasattr(molecule, 'EFP') else molecule.nuclear_repulsion_energy()}")
    logger.debug("w/EFP" if hasattr(molecule, "EFP") else pp.pformat(molecule.to_dict()))
    wfn = procedures["properties"][lowername](lowername, molecule=molecule, **kwargs)
    logger.info(f"Return properties(): {core.variable('CURRENT ENERGY')}")

    basisstash.restore()
    optstash.restore()

    if return_wfn:
        return (core.variable('CURRENT ENERGY'), wfn)
    else:
        return core.variable('CURRENT ENERGY')


def optimize_geometric(name, **kwargs):

    import qcelemental as qcel
    from qcelemental.util import which_import

    if not which_import('geometric', return_bool=True):
        raise ModuleNotFoundError('Python module geometric not found. Solve by installing it: `conda install -c conda-forge geometric` or `pip install geometric`')
    import geometric

    class Psi4NativeEngine(geometric.engine.Engine):
        """
        Internally run an energy and gradient calculation for geometric 
        """
        def __init__(self, p4_name, p4_mol, p4_return_wfn, **p4_kwargs):
    
            self.p4_name = p4_name
            self.p4_mol = p4_mol
            self.p4_return_wfn = p4_return_wfn
            self.p4_kwargs = p4_kwargs
    
            molecule = geometric.molecule.Molecule()
            molecule.elem = [p4_mol.symbol(i).capitalize() for i in range(p4_mol.natom())]
            molecule.xyzs = [p4_mol.geometry().np * qcel.constants.bohr2angstroms] 
            molecule.build_bonds()
                                 
            super(Psi4NativeEngine, self).__init__(molecule)
    
        def calc(self, coords, dirname, read_data=False):
            self.p4_mol.set_geometry(core.Matrix.from_array(coords.reshape(-1,3)))
            self.p4_mol.update_geometry()
            if self.p4_return_wfn:
                g, wfn = gradient(self.p4_name, return_wfn=True, molecule=self.p4_mol, **self.p4_kwargs)
                self.p4_wfn = wfn
            else:
                g = gradient(self.p4_name, return_wfn=False, molecule=self.p4_mol, **self.p4_kwargs)
            e = core.variable('CURRENT ENERGY')
            return {'energy': e, 'gradient': g.np.ravel()}

    return_wfn = kwargs.pop('return_wfn', False)
    return_history = kwargs.pop('return_history', False)

    if return_history:
        step_energies = []
        step_gradients = []
        step_coordinates = []

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())

    # Do not change orientation or COM
    molecule.fix_orientation(True)
    molecule.fix_com(True)
    molecule.update_geometry()

    # Get geometric-specific options
    optimizer_keywords = {k.lower(): v for k, v in kwargs.get("optimizer_keywords", {}).items()}

    core.print_out('\n')
    core.print_out("\n  ==> GeomeTRIC Optimizer <==                                                                   ~\n")
                                 
    # Default to Psi4 maxiter unless overridden
    if 'maxiter' not in optimizer_keywords:
        optimizer_keywords['maxiter'] = core.get_global_option('GEOM_MAXITER')

    # Default to Psi4 geometry convergence criteria unless overridden 
    if 'convergence_set' not in optimizer_keywords:
        optimizer_keywords['convergence_set'] = core.get_global_option('G_CONVERGENCE')

        # GeomeTRIC doesn't know these convergence criterion
        if optimizer_keywords['convergence_set'] in ['CFOUR', 'QCHEM', 'MOLPRO']:
            core.print_out(f"\n  Psi4 convergence criteria {optimizer_keywords['convergence_set']:6s} not recognized by GeomeTRIC, switching to GAU_TIGHT          ~")
            optimizer_keywords['convergence_set'] = 'GAU_TIGHT'

    engine = Psi4NativeEngine(name, molecule, return_wfn, **kwargs)
    M = engine.M
    
    # Handle constraints
    constraints_dict = {k.lower(): v for k, v in optimizer_keywords.get("constraints", {}).items()}
    constraints_string = geometric.run_json.make_constraints_string(constraints_dict)
    Cons, CVals = None, None
    if constraints_string:
        if 'scan' in constraints_dict:
            raise ValueError("Coordinate scans are not yet available through the Psi4-GeomeTRIC interface")
        Cons, CVals = geometric.prepare.parse_constraints(M, constraints_string)
    
    # Set up the internal coordinate system
    coordsys = optimizer_keywords.get('coordsys', 'tric')
    CoordSysDict = {
        'cart': (geometric.internal.CartesianCoordinates, False, False),
        'prim': (geometric.internal.PrimitiveInternalCoordinates, True, False),
        'dlc': (geometric.internal.DelocalizedInternalCoordinates, True, False),
        'hdlc': (geometric.internal.DelocalizedInternalCoordinates, False, True),
        'tric': (geometric.internal.DelocalizedInternalCoordinates, False, False)
    }
    
    # Build internal coordinates
    CoordClass, connect, addcart = CoordSysDict[coordsys.lower()]
    IC = CoordClass(
        M,
        build=True,
        connect=connect,
        addcart=addcart,
        constraints=Cons,
        cvals=CVals[0] if CVals is not None else None)
    
    # Get initial coordinates in bohr
    coords = M.xyzs[0].flatten() / qcel.constants.bohr2angstroms

    # Setup an optimizer object
    params = geometric.optimize.OptParams(**optimizer_keywords)
    optimizer = geometric.optimize.Optimizer(coords, M, IC, engine, None, params)
    
    # TODO: print constraints
    # IC.printConstraints(coords, thre=-1)
    optimizer.calcEnergyForce()
    optimizer.prepareFirstStep()
    grms, gmax = optimizer.calcGradNorm()
    conv_gmax = '*' if gmax < params.Convergence_gmax else ' '
    conv_grms = '*' if grms < params.Convergence_grms else ' '
    core.print_out("\n  Measures of convergence in internal coordinates in au.                                        ~")
    core.print_out("\n  Criteria marked as inactive (o), active & met (*), and active & unmet ( ).                    ~")
    core.print_out("\n  --------------------------------------------------------------------------------------------- ~")
    core.print_out("\n   Step     Total Energy     Delta E     MAX Force     RMS Force      MAX Disp      RMS Disp    ~")
    core.print_out("\n  --------------------------------------------------------------------------------------------- ~")
    core.print_out((f"\n    Convergence Criteria  {params.Convergence_energy:10.2e}    "
                    f"{params.Convergence_gmax:10.2e}    {params.Convergence_grms:10.2e}    "
                    f"{params.Convergence_dmax:10.2e}    {params.Convergence_drms:10.2e}    ~"))
    core.print_out("\n  --------------------------------------------------------------------------------------------- ~")

    core.print_out((f"\n   {optimizer.Iteration:4d} {optimizer.E:16.8e}    --------    "
                    f"{gmax:10.2e} {conv_gmax}  {grms:10.2e} {conv_grms}    --------      --------    ~"))
    while True:
        if optimizer.state == geometric.optimize.OPT_STATE.CONVERGED:
            core.print_out("\n\n  Optimization converged!                                                                       ~\n")
            break
        elif optimizer.state == geometric.optimize.OPT_STATE.FAILED:
            core.print_out("\n\n  Optimization failed to converge!                                                              ~\n")
            break
        optimizer.step()
        optimizer.calcEnergyForce()
        optimizer.evaluateStep()
        grms, gmax = optimizer.calcGradNorm()
        drms, dmax = geometric.optimize.calc_drms_dmax(optimizer.X, optimizer.Xprev)
        conv_energy = '*' if np.abs(optimizer.E - optimizer.Eprev) < params.Convergence_energy else ' '
        conv_gmax = '*' if gmax < params.Convergence_gmax else ' '
        conv_grms = '*' if grms < params.Convergence_grms else ' '
        conv_dmax = '*' if dmax < params.Convergence_dmax else ' '
        conv_drms = '*' if drms < params.Convergence_drms else ' '
        core.print_out((f'\n   {optimizer.Iteration:4d} {optimizer.E:16.8e}  '
                        f'{optimizer.E-optimizer.Eprev:10.2e} {conv_energy}  {gmax:10.2e} {conv_gmax}  '
                        f'{grms:10.2e} {conv_grms}  {dmax:10.2e} {conv_dmax}  {drms:10.2e} {conv_drms}  ~'))

        if return_history:
            step_energies.append(optimizer.E)
            step_coordinates.append(core.Matrix.from_array(optimizer.X.reshape(-1,3)))
            step_gradients.append(core.Matrix.from_array(optimizer.gradx.reshape(-1,3)))

    return_energy = optimizer.E
    opt_geometry = core.Matrix.from_array(optimizer.X.reshape(-1,3))
    molecule.set_geometry(opt_geometry)
    molecule.update_geometry()
    core.print_out(f'\n  Final Energy : {return_energy} \n')
    core.print_out('\n  Final Geometry : \n')
    molecule.print_in_input_format()

    if return_history:
        history = {
            'energy': step_energies,
            'gradient': step_gradients,
            'coordinates': step_coordinates,
        }

    if return_wfn:
        wfn = engine.p4_wfn

    if return_wfn and return_history:
        return (return_energy, wfn, history)
    elif return_wfn and not return_history:
        return (return_energy, wfn)
    elif return_history and not return_wfn:
        return (return_energy, history)
    else:
        return return_energy


def optimize(name, **kwargs):
    r"""Function to perform a geometry optimization.

    :aliases: opt()

    :returns: *float* |w--w| Total electronic energy of optimized structure in Hartrees.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| energy and wavefunction when **return_wfn** specified.

    :raises: :py:class:`psi4.driver.OptimizationConvergenceError` if :term:`GEOM_MAXITER <GEOM_MAXITER (OPTKING)>` exceeded without reaching geometry convergence.

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`CURRENT ENERGY`

    :type name: str
    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the database. May be any valid argument to
        :py:func:`psi4.driver.energy`.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type return_wfn: :ref:`boolean <op_py_boolean>`
    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|

        Indicate to additionally return the :py:class:`~psi4.core.Wavefunction`
        calculation result as the second element (after *float* energy) of a tuple.

    :type return_history: :ref:`boolean <op_py_boolean>`
    :param return_history: ``'on'`` || |dl| ``'off'`` |dr|

        Indicate to additionally return dictionary of lists of geometries,
        energies, and gradients at each step in the optimization.

    :type engine: str
    :param engine: |dl| ``'optking'`` |dr| || ``'geometric'``

        Indicates the optimization engine to use, which can be either Psi4's
        native Optking optimizer or the GeomeTRIC program.

    :type optimizer_keywords: dict
    :param optimizer_keywords: Extra options passed to the GeomeTRIC or optking optimizers

        Indicates additional options to be passed to the GeomeTRIC optimizer if
        chosen as the optimization engine. Alternatively, can be used to set optking options
        that are not currently recognized by Psi4.

    :type func: :ref:`function <op_py_function>`
    :param func: |dl| ``gradient`` |dr| || ``energy`` || ``cbs``

        Indicates the type of calculation to be performed on the molecule.
        The default dertype accesses ``'gradient'`` or ``'energy'``, while
        ``'cbs'`` performs a multistage finite difference calculation.
        If a nested series of python functions is intended (see :ref:`sec:intercalls`),
        use keyword ``opt_func`` instead of ``func``.

    :type dertype: :ref:`dertype <op_py_dertype>`
    :param dertype: ``'gradient'`` || ``'energy'``

        Indicates whether analytic (if available) or finite difference
        optimization is to be performed.

    :type hessian_with: str
    :param hessian_with: ``'scf'`` || ``'mp2'`` || etc.

        Indicates the computational method with which to perform a hessian
        analysis to guide the geometry optimization.

    .. warning:: Optimizations where the molecule is specified in Z-matrix format
       with dummy atoms will result in the geometry being converted to a Cartesian representation.

    .. note:: Analytic gradients area available for all methods in the table
        below. Optimizations with other methods in the energy table proceed
        by finite differences.

    .. _`table:grad_gen`:

    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | name                    | calls method                                                                                                  |
    +=========================+===============================================================================================================+
    | efp                     | efp-only optimizations                                                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT) :ref:`[manual] <sec:scf>` :ref:`[details] <dd_b3lyp>`   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | hf                      | HF self consistent field (SCF) :ref:`[manual] <sec:scf>` :ref:`[details] <dd_hf>`                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dct                     | density cumulant (functional) theory :ref:`[manual] <sec:dct>`                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp2                     | 2nd-order |MollerPlesset| perturbation theory (MP2) :ref:`[manual] <sec:dfmp2>` :ref:`[details] <dd_mp2>`     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp3                     | 3rd-order |MollerPlesset| perturbation theory (MP3) :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <dd_mp3>` |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp2.5                   | average of MP2 and MP3 :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <dd_mp2p5>`                            |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp2                    | orbital-optimized second-order MP perturbation theory :ref:`[manual] <sec:occ_oo>` :ref:`[details] <dd_omp2>` |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp3                    | orbital-optimized third-order MP perturbation theory :ref:`[manual] <sec:occ_oo>` :ref:`[details] <dd_omp3>`  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp2.5                  | orbital-optimized MP2.5 :ref:`[manual] <sec:occ_oo>` :ref:`[details] <dd_omp2p5>`                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | oremp2                  | orbital-optimized REMP2 :ref:`[manual] <sec:occ_oo>` :ref:`[details] <dd_oremp2>`                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | lccd                    | Linear CCD :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <dd_lccd>`                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | olccd                   | orbital optimized LCCD :ref:`[manual] <sec:occ_oo>` :ref:`[details] <dd_olccd>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cc2                     | approximate coupled cluster singles and doubles (CC2) :ref:`[manual] <sec:cc>` :ref:`[details] <dd_cc2>`      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccd                     | coupled cluster doubles (CCD) :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <dd_ccd>`                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD) :ref:`[manual] <sec:cc>` :ref:`[details] <dd_ccsd>`                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd(t)                 | CCSD with perturbative triples (CCSD(T)) :ref:`[manual] <sec:cc>` :ref:`[details] <dd_ccsd_prt_pr>`           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | eom-ccsd                | equation of motion (EOM) CCSD :ref:`[manual] <sec:eomcc>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+

    .. _`table:grad_scf`:

    .. include:: /autodoc_dft_opt.rst

    .. include:: /cfour_table_grad.rst


    :examples:

    >>> # [1] Analytic hf optimization
    >>> optimize('hf')

    >>> # [2] Finite difference mp5 optimization with gradient
    >>> #     printed to output file
    >>> e, wfn = opt('mp5', return_wfn='yes')
    >>> wfn.gradient().print_out()

    >>> # [3] Can automatically perform complete basis set extrapolations
    >>> optimize('MP2/cc-pV([D,T]+d)Z')

    >>> # [4] Can automatically perform delta corrections that include extrapolations
    >>> # even with a user-defined extrapolation formula. See sample inputs named
    >>> # cbs-xtpl* for more examples of this input style
    >>> optimize("MP2/aug-cc-pv([d,t]+d)z + d:ccsd(t)/cc-pvdz", corl_scheme=myxtplfn_2)

    >>> # [5] Get info like geometry, gradient, energy back after an
    >>> #     optimization fails. Note that the energy and gradient
    >>> #     correspond to the last optimization cycle, whereas the
    >>> #     geometry (by default) is the anticipated *next* optimization step.
    >>> try:
    >>>     optimize('hf/cc-pvtz')
    >>> except psi4.OptimizationConvergenceError as ex:
    >>>     next_geom_coords_as_numpy_array = np.asarray(ex.wfn.molecule().geometry())

    """
    kwargs = p4util.kwargs_lower(kwargs)

    engine = kwargs.pop('engine', 'optking')
    if engine == 'geometric':
        return optimize_geometric(name, **kwargs)
    elif engine != 'optking':
        raise ValidationError(f"Optimizer {engine} is not supported.")

    import optking

    name = driver_util.upgrade_interventions(name)
    if hasattr(name, '__call__'):
        lowername = name
        custom_gradient = True
    else:
        lowername = name.lower()
        custom_gradient = False

    return_wfn = kwargs.pop('return_wfn', False)

    return_history = kwargs.pop('return_history', False)

    # This might be incorrect now?
    if custom_gradient and core.has_option_changed('OPTKING', 'FULL_HESS_EVERY'):
        raise ValidationError("Optimize: Does not support custom Hessian's yet.")
    else:
        hessian_with_method = kwargs.get('hessian_with', lowername)

    _filter_renamed_methods("optimize", lowername)

    optstash = p4util.OptionsState(
        ['FINDIF', 'HESSIAN_WRITE'],
        ['OPTKING', 'CART_HESS_READ'],
        ['SCF', 'GUESS_PERSIST'],  # handle on behalf of cbs()
        ['SCF', 'GUESS'])

    n = kwargs.get('opt_iter', 1)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())

    # If we are freezing cartesian, do not orient or COM
    if any([core.get_option('OPTKING', 'FROZEN_CARTESIAN'), core.get_option('OPTKING', 'EXT_FORCE_CARTESIAN')]):
        if molecule.has_zmatrix():
            raise ValidationError("Job includes cartesian coordinate constraints. This cannot be fully "
                                  "obeyed due to zmatrix in input. Please convert your zmatrix to cartesian "
                                  "coordinates if cartesian constraints are needed ")
        molecule.fix_orientation(True)
        molecule.fix_com(True)
    molecule.update_geometry()

    if core.get_option('OPTKING', 'OPT_RESTART'):
        # Recreate all of optking's internal classes to restart an optimization
        # This has not been well tested - Experimental
        opt_object = optking.opt_helper.CustomHelper(molecule)
        with open(f"{core.get_writer_file_prefix(molecule.name())}.1.dat", 'r') as f:
            stashed_opt = json.load(f)
        opt_object.from_dict(stashed_opt)
    else:
        # create an OptHelper to run an optimization through
        # Optking will ignore any keywords it doesn't recognize.
        params = p4util.prepare_options_for_modules()
        optimizer_params = {k: v.get('value') for k, v in params.pop("OPTKING").items() if v.get('has_changed')}
        optimizer_params.update(kwargs.get("optimizer_keywords", {}))
        opt_object = optking.opt_helper.CustomHelper(molecule, params=optimizer_params)

    initial_sym = molecule.schoenflies_symbol()
    while n <= core.get_option('OPTKING', 'GEOM_MAXITER'):
        current_sym = molecule.schoenflies_symbol()
        if initial_sym != current_sym:

            if any([core.get_option('OPTKING', 'FROZEN_CARTESIAN'), core.get_option('OPTKING', 'EXT_FORCE_CARTESIAN')]):
                raise ValidationError("Symmetrize cannot be called while cartesian constraints are active "
                                      "symmetrize was about to be called. Please check symmetry dependent input "
                                      ", such as DOCC, is correct or turn off symmetry")

            # Try to resymmetrize molecule if slightly broken.
            molecule.symmetrize(core.get_option("OPTKING", "CARTESIAN_SYM_TOLERANCE"))

            if molecule.schoenflies_symbol() != initial_sym:
                raise ValidationError("Point group changed! (%s <-- %s) You should restart "
                                      "using the last geometry in the output, after "
                                      "carefully making sure all symmetry-dependent "
                                      "input, such as DOCC, is correct." % (current_sym, initial_sym))

        kwargs['opt_iter'] = n
        core.set_variable('GEOMETRY ITERATIONS', n)

        # Use orbitals from previous iteration as a guess
        #   set within loop so that can be influenced by fns to optimize (e.g., cbs)
        if (n > 1) and (not core.get_option('SCF', 'GUESS_PERSIST')):
            core.set_local_option('SCF', 'GUESS', 'READ')

        # We'll currently ignore the possibility that the gradient isn't needed
        opt_calcs = opt_object.calculations_needed() # tuple of strings ('energy', 'gradient', etc)

        # Compute the gradient - no longer need to worry about opt_data being wiped
        G, wfn = gradient(lowername, return_wfn=True, molecule=molecule, **kwargs)
        thisenergy = core.variable('CURRENT ENERGY')
        opt_object.E = thisenergy
        opt_object.gX = G.np

        if core.get_option('OPTKING', 'CART_HESS_READ') and (n == 1):
            opt_object.params.cart_hess_read = True
            opt_object.params.hessian_file = f"{core.get_writer_file_prefix(molecule.name())}.hess"
                # compute Hessian as requested; frequency wipes out gradient so stash it
        elif 'hessian' in opt_calcs:
            # compute hessian as requested.

            # procedures proctable analytic hessians
            _, hess_wfn = frequencies(hessian_with_method,
                                      molecule=molecule,
                                      ref_gradient=G,
                                      return_wfn=True,
                                      **kwargs)
            opt_object.HX = hess_wfn.hessian().np

        # force optking to update its molecule to psi4's.
        # This allows for psi4 to rotate as desired. If optimizing in cartesians. rotation is not allowed
        # Process gradient / hessian. Take step. Print summary to output for user
        opt_object.molsys.geom = molecule.geometry().np
        core.print_out(opt_object.pre_step_str())  # print optking's molecule
        opt_object.compute()  # process E, gX, H
        try:
            opt_object.take_step()
        except optking.exceptions.AlgError:
            # Optking encountered an algorithm error and reset.
            if not opt_object.HX:
                n += 1
                continue
            else:
                raise ConvergenceError(
                        "Psi4 caught an AlgError. This should only happen after optking resets the history"
                        "and needs another Hessian",
                        n,
                        wfn
                )

        core.print_out(opt_object.post_step_str())  # print convergence and step info

        # Update psi4's molecule with new step. (Psi4 can rotate this molecule)
        molecule.set_geometry(core.Matrix.from_array(opt_object.molsys.geom))
        molecule.update_geometry()

        opt_status = opt_object.status()  # Query optking for convergence, failure or continuing opt.
        if opt_status == 'CONVERGED':

            # Last geom is normally last in history. For IRC last geom is last in IRC trajectory
            # Not sure how to handle ensuring that wfn corresponds to last point.
            final_energy, final_geom = opt_object.summarize_result()

            # Changing environment to optimized geometry as expected by user
            molecule.set_geometry(core.Matrix.from_array(final_geom))
            molecule.update_geometry()

            print('Optimizer: Optimization complete!')
            core.print_out('\n    Final optimized geometry and variables:\n')
            molecule.print_in_input_format()

            for postcallback in hooks['optimize']['post']:
                postcallback(lowername, wfn=wfn, **kwargs)
            core.clean()

            optstash.restore()

            if return_history:
                history = {
                    'energy': [step.E for step in opt_object.history.steps],
                    'gradient': [step.cart_grad for step in opt_object.history.steps],
                    'coordinates': [step.geom for step in opt_object.history.steps],
                }

            # Create OptimizationResult like Schema. Not validated since optimize() does not pass AtomicResults.
            opt_data = opt_object.close()
            if core.get_option('OPTKING', 'WRITE_OPT_HISTORY'):
                with open(f"{core.get_writer_file_prefix(molecule.name())}.opt.json", 'w+') as f:
                    json.dump(opt_data, f, indent=2)

            if return_wfn and return_history:
                return (thisenergy, wfn, history)
            elif return_wfn and not return_history:
                return (thisenergy, wfn)
            elif return_history and not return_wfn:
                return (thisenergy, history)
            else:
                return thisenergy

        elif opt_status == 'FAILED':

            print('Optimizer: Optimization failed!')
            molecule.set_geometry(core.Matrix.from_array(opt_object.molsys.geom))
            molecule.update_geometry()
            core.clean()
            optstash.restore()

            opt_data = opt_object.to_dict()
            if core.get_option('OPTKING', 'WRITE_OPT_HISTORY'):
                with open(f"{core.get_writer_file_prefix(molecule.name())}.opt.json", 'w+') as f:
                    json.dump(opt_data, f, indent=2)

            raise OptimizationConvergenceError("""geometry optimization""", n - 1, wfn)

        core.print_out('\n    Structure for next step:\n')
        molecule.print_in_input_format()

        n += 1

    optstash.restore()
    raise OptimizationConvergenceError("""geometry optimization""", n - 1, wfn)


def hessian(name, **kwargs):
    r"""Function complementary to :py:func:`~psi4.driver.frequency`. Computes force
    constants, deciding analytic, finite difference of gradients, or
    finite difference of energies.

    :returns: :py:class:`~psi4.core.Matrix` |w--w| Total non-mass-weighted electronic Hessian in Hartrees/Bohr/Bohr.

    :returns: (:py:class:`~psi4.core.Matrix`, :py:class:`~psi4.core.Wavefunction`) |w--w| Hessian and wavefunction when **return_wfn** specified.

    :examples:

    >>> # [1] Frequency calculation without thermochemical analysis
    >>> hessian('mp3')

    >>> # [2] Frequency calc w/o thermo analysis getting the Hessian
    >>> #     in file, core.Matrix, and np.array forms
    >>> set hessian_write on
    >>> H, wfn = hessian('ccsd', return_wfn=True)
    >>> wfn.hessian().print_out()
    >>> np.array(H)

    """
    ## First half of this fn -- entry means user wants a 2nd derivative by any means

    kwargs = p4util.kwargs_lower(kwargs)
    basisstash = p4util.OptionsState(['BASIS'])
    return_wfn = kwargs.pop('return_wfn', False)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    # Convert wrapper directives from options (where ppl know to find them) to kwargs (suitable for non-globals transmitting)
    kwargs['findif_verbose'] = core.get_option("FINDIF", "PRINT")
    kwargs['findif_stencil_size'] = core.get_option("FINDIF", "POINTS")
    kwargs['findif_step_size'] = core.get_option("FINDIF", "DISP_SIZE")

    # Select certain irreps
    irrep = kwargs.pop('irrep', -1)
    if irrep == -1:
        pass  # do all irreps
    else:
        irrep = driver_util.parse_cotton_irreps(irrep, molecule.schoenflies_symbol())
        irrep -= 1  # A1 irrep is externally 1, internally 0
    kwargs['findif_irrep'] = irrep

    ## Pre-planning interventions

    # * Trip on function or alias as name
    lowername = driver_util.upgrade_interventions(name)
    _filter_renamed_methods("hessian", lowername)

    # * Prevent methods that do not have associated derivatives
    if lowername in energy_only_methods:
        raise ValidationError(f"`hessian('{name}')` does not have an associated Hessian.")

    # * Avert pydantic anger at incomplete modelchem spec
    userbas = core.get_global_option('BASIS') or kwargs.get('basis')
    if lowername in integrated_basis_methods and userbas is None:
        kwargs['basis'] = '(auto)'

    # Are we planning?
    plan = task_planner.task_planner("hessian", lowername, molecule, **kwargs)
    logger.debug('HESSIAN PLAN')
    logger.debug(pp.pformat(plan.dict()))

    if kwargs.get("return_plan", False):
        # Plan-only requested
        return plan

    elif not isinstance(plan, AtomicComputer):
        # Advanced "Computer" active
        plan.compute()
        return plan.get_psi_results(return_wfn=return_wfn)

    else:
        # We have unpacked to an AtomicInput
        lowername = plan.method
        basis = plan.basis
        core.set_global_option("BASIS", basis)

    ## Second half of this fn -- entry means program running exactly analytic 2nd derivative

    _filter_renamed_methods("frequency", lowername)
    core.clean_variables()

    optstash = p4util.OptionsState(
        ['FINDIF', 'HESSIAN_WRITE'],
        ['FINDIF', 'FD_PROJECT'],
    )

    # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
    optstash_conv = driver_util.negotiate_convergence_criterion((2, 2), lowername, return_optstash=True)

    # At stationary point?
    if 'ref_gradient' in kwargs:
        core.print_out("""hessian() using ref_gradient to assess stationary point.\n""")
        G0 = kwargs['ref_gradient']
    else:
        tmpkwargs = copy.deepcopy(kwargs)
        tmpkwargs.pop('dertype', None)
        G0 = gradient(lowername, molecule=molecule, **tmpkwargs)
    translations_projection_sound, rotations_projection_sound = _energy_is_invariant(G0.rms())
    core.print_out(
        '\n  Based on options and gradient (rms={:.2E}), recommend {}projecting translations and {}projecting rotations.\n'
        .format(G0.rms(), '' if translations_projection_sound else 'not ',
                '' if rotations_projection_sound else 'not '))
    if not core.has_option_changed('FINDIF', 'FD_PROJECT'):
        core.set_local_option('FINDIF', 'FD_PROJECT', rotations_projection_sound)

    # We have the desired method. Do it.
    logger.info(f"Compute hessian(): method={lowername}, basis={core.get_global_option('BASIS').lower()}, molecule={molecule.name()}, nre={'w/EFP' if hasattr(molecule, 'EFP') else molecule.nuclear_repulsion_energy()}")
    logger.debug("w/EFP" if hasattr(molecule, "EFP") else pp.pformat(molecule.to_dict()))
    core.print_out("""hessian() will perform analytic frequency computation.\n""")
    wfn = procedures['hessian'][lowername](lowername, molecule=molecule, **kwargs)
    logger.info(f"Return hessian(): {wfn.energy()}")
    logger.info(nppp(wfn.hessian().np))

    wfn.set_gradient(G0)
    basisstash.restore()
    optstash.restore()
    optstash_conv.restore()

    #if isinstance(lowername, str) and lowername in procedures['energy']:
    #    # this correctly filters out cbs fn and "hf/cc-pvtz"
    #    # it probably incorrectly filters out mp5, but reconsider in DDD
    #    core.set_variable(f"CURRENT HESSIAN", H)
    #    core.set_variable(f"{lowername.upper()} TOTAL HESSIAN", H)
    #    core.set_variable(f"{lowername.upper()} TOTAL GRADIENT", G0)
    #    wfn.set_variable(f"{lowername.upper()} TOTAL HESSIAN", H)
    #    wfn.set_variable(f"{lowername.upper()} TOTAL GRADIENT", G0)
    # TODO: check that current energy's being set to the right figure when this code is actually used
    core.set_variable('CURRENT ENERGY', wfn.energy())
    core.set_variable("CURRENT GRADIENT", G0)
    driver_findif.hessian_write(wfn)

    if return_wfn:
        return (wfn.hessian(), wfn)
    else:
        return wfn.hessian()


def frequency(name, **kwargs):
    r"""Function to compute harmonic vibrational frequencies.

    :aliases: frequencies(), freq()

    :returns: *float* |w--w| Total electronic energy in Hartrees.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| energy and wavefunction when **return_wfn** specified.

    :type name: str
    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type return_wfn: :ref:`boolean <op_py_boolean>`
    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|

        Indicate to additionally return the :py:class:`~psi4.core.Wavefunction`
        calculation result as the second element (after *float* energy) of a tuple.
        Arrays of frequencies and the Hessian can be accessed through the wavefunction.

    :type func: :ref:`function <op_py_function>`
    :param func: |dl| ``gradient`` |dr| || ``energy`` || ``cbs``

        Indicates the type of calculation to be performed on the molecule.
        The default dertype accesses ``'gradient'`` or ``'energy'``, while
        ``'cbs'`` performs a multistage finite difference calculation.
        If a nested series of python functions is intended (see :ref:`sec:intercalls`),
        use keyword ``freq_func`` instead of ``func``.

    :type dertype: :ref:`dertype <op_py_dertype>`
    :param dertype: |dl| ``'hessian'`` |dr| || ``'gradient'`` || ``'energy'``

        Indicates whether analytic (if available- they're not), finite
        difference of gradients (if available) or finite difference of
        energies is to be performed.

    :type irrep: int or str
    :param irrep: |dl| ``-1`` |dr| || ``1`` || ``'b2'`` || ``'App'`` || etc.

        Indicates which symmetry block (:ref:`Cotton <table:irrepOrdering>` ordering) of vibrational
        frequencies to be computed. ``1``, ``'1'``, or ``'a1'`` represents
        :math:`a_1`, requesting only the totally symmetric modes.
        ``-1`` indicates a full frequency calculation.

    .. note:: Analytic hessians are only available for RHF and UHF. For all other methods, Frequencies will
        proceed through finite differences according to availability of gradients or energies.

    .. _`table:freq_gen`:

    +-------------------------+-----------------------------------------------------------------------------------------------------------------+
    | name                    | calls method                                                                                                    |
    +=========================+=================================================================================================================+
    | scf                     | Hartree--Fock (HF) or LSDA density functional theory (DFT) :ref:`[manual] <sec:scf>` :ref:`[details] <dd_svwn>` |
    +-------------------------+-----------------------------------------------------------------------------------------------------------------+

    :examples:

    >>> # [1] Frequency calculation for all modes through highest available derivatives
    >>> frequency('ccsd')

    >>> # [2] Frequency calculation for b2 modes through finite difference of gradients
    >>> #     printing lowest mode frequency to screen and Hessian to output
    >>> E, wfn = frequencies('scf', dertype=1, irrep=4, return_wfn=True)
    >>> print wfn.frequencies().get(0, 0)
    >>> wfn.hessian().print_out()

    >>> # [3] Frequency calculation at default conditions and Hessian reuse at STP
    >>> E, wfn = freq('mp2', return_wfn=True)
    >>> set t 273.15
    >>> set p 100000
    >>> thermo(wfn, wfn.frequencies())

    >>> # [4] Opt+Freq, skipping the gradient recalc at the start of the Hessian
    >>> e, wfn = optimize('hf', return_wfn=True)
    >>> frequencies('hf', ref_gradient=wfn.gradient())

    """
    kwargs = p4util.kwargs_lower(kwargs)

    return_wfn = kwargs.pop('return_wfn', False)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    # Compute the hessian
    H, wfn = hessian(name, return_wfn=True, molecule=molecule, **kwargs)

    # Project final frequencies?
    if wfn.gradient():  # available for analytic and any findif including totally symmetric space
        gradient_rms = wfn.gradient().rms()
    else:
        gradient_rms = 1  # choose to force non-projection of rotations
    translations_projection_sound, rotations_projection_sound = _energy_is_invariant(gradient_rms)

    project_trans = kwargs.get('project_trans', translations_projection_sound)
    project_rot = kwargs.get('project_rot', rotations_projection_sound)

    irrep = kwargs.get('irrep', None)
    vibinfo = vibanal_wfn(wfn, irrep=irrep, project_trans=project_trans, project_rot=project_rot)
    wfn.frequency_analysis = vibinfo

    for postcallback in hooks['frequency']['post']:
        postcallback(lowername, wfn=wfn, **kwargs)

    if return_wfn:
        return (core.variable('CURRENT ENERGY'), wfn)
    else:
        return core.variable('CURRENT ENERGY')


def vibanal_wfn(
    wfn: core.Wavefunction,
    hess: Optional[np.ndarray] = None,
    irrep: Optional[Union[int, str]] = None,
    molecule=None,
    project_trans: bool = True,
    project_rot: bool = True,
) -> Dict[str, np.ndarray]:
    """Function to perform analysis of a hessian or hessian block, specifically...
    calling for and printing vibrational and thermochemical analysis, setting thermochemical variables,
    and writing the vibrec and normal mode files.

    Parameters
    ----------
    wfn
        The wavefunction which had its Hessian computed.
    hess
        Hessian to analyze, if not the hessian in wfn.
        (3*nat, 3*nat) non-mass-weighted Hessian in atomic units, [Eh/a0/a0].
    irrep
        The irrep for which frequencies are calculated. Thermochemical analysis
        is skipped if this is specified (non-None),
        as only one symmetry block of the hessian has been computed.
    molecule : :py:class:`~psi4.core.Molecule` or qcdb.Molecule, optional
        The molecule to pull information from, if not the molecule in wfn. Must at least have similar
        geometry to the molecule in wfn.
    project_trans
        Should translations be projected in the harmonic analysis?
    project_rot
        Should rotations be projected in the harmonic analysis?

    Returns
    -------
    vibinfo : ~typing.Dict[str, ~numpy.ndarray]
        A dictionary of vibrational information. See :py:func:`~psi4.driver.qcdb.vib.harmonic_analysis`

    """

    if hess is None:
        nmwhess = np.asarray(wfn.hessian())
    else:
        nmwhess = hess

    dipder = wfn.variables().get("CURRENT DIPOLE GRADIENT", None)
    if dipder is not None:
        dipder = np.asarray(dipder).T

    mol = wfn.molecule()
    geom = np.asarray(mol.geometry())
    symbols = [mol.symbol(at) for at in range(mol.natom())]

    vibrec = {'molecule': mol.to_dict(np_out=False), 'hessian': nmwhess.tolist()}

    if molecule is not None:
        molecule.update_geometry()
        if mol.natom() != molecule.natom():
            raise ValidationError('Impostor molecule trying to be analyzed! natom {} != {}'.format(
                mol.natom(), molecule.natom()))
        if abs(mol.nuclear_repulsion_energy() - molecule.nuclear_repulsion_energy()) > 1.e-6:
            raise ValidationError('Impostor molecule trying to be analyzed! NRE {} != {}'.format(
                mol.nuclear_repulsion_energy(), molecule.nuclear_repulsion_energy()))
        if not np.allclose(np.asarray(mol.geometry()), np.asarray(molecule.geometry()), atol=1.e-6):
            core.print_out(
                'Warning: geometry center/orientation mismatch. Normal modes may not be in expected coordinate system.'
            )
        #    raise ValidationError('Impostor molecule trying to be analyzed! geometry\n{}\n   !=\n{}'.format(
        #        np.asarray(mol.geometry()), np.asarray(molecule.geometry())))
        mol = molecule

    m = np.asarray([mol.mass(at) for at in range(mol.natom())])
    irrep_labels = mol.irrep_labels()

    vibinfo, vibtext = qcdb.vib.harmonic_analysis(nmwhess,
                                                  geom,
                                                  m,
                                                  wfn.basisset(),
                                                  irrep_labels,
                                                  dipder=dipder,
                                                  project_trans=project_trans,
                                                  project_rot=project_rot)
    vibrec.update({k: qca.json() for k, qca in vibinfo.items()})

    core.print_out(vibtext)
    core.print_out(qcdb.vib.print_vibs(vibinfo, shortlong=True, normco='x', atom_lbl=symbols))

    if core.has_option_changed('THERMO', 'ROTATIONAL_SYMMETRY_NUMBER'):
        rsn = core.get_option('THERMO', 'ROTATIONAL_SYMMETRY_NUMBER')
    else:
        rsn = mol.rotational_symmetry_number()

    if irrep is None:
        therminfo, thermtext = qcdb.vib.thermo(
            vibinfo,
            T=core.get_option("THERMO", "T"),  # 298.15 [K]
            P=core.get_option("THERMO", "P"),  # 101325. [Pa]
            multiplicity=mol.multiplicity(),
            molecular_mass=np.sum(m),
            sigma=rsn,
            rotor_type=mol.rotor_type(),
            rot_const=np.asarray(mol.rotational_constants()),
            E0=core.variable('CURRENT ENERGY'))  # someday, wfn.energy()
        vibrec.update({k: qca.json() for k, qca in therminfo.items()})

        core.set_variable("ZPVE", therminfo['ZPE_corr'].data)  # P::e THERMO
        core.set_variable("THERMAL ENERGY CORRECTION", therminfo['E_corr'].data)  # P::e THERMO
        core.set_variable("ENTHALPY CORRECTION", therminfo['H_corr'].data)  # P::e THERMO
        core.set_variable("GIBBS FREE ENERGY CORRECTION", therminfo['G_corr'].data)  # P::e THERMO

        core.set_variable("ZERO K ENTHALPY", therminfo['ZPE_tot'].data)  # P::e THERMO
        core.set_variable("THERMAL ENERGY", therminfo['E_tot'].data)  # P::e THERMO
        core.set_variable("ENTHALPY", therminfo['H_tot'].data)  # P::e THERMO
        core.set_variable("GIBBS FREE ENERGY", therminfo['G_tot'].data)  # P::e THERMO

        core.print_out(thermtext)
    else:
        core.print_out('  Thermochemical analysis skipped for partial frequency calculation.\n')

    if core.get_option('FINDIF', 'HESSIAN_WRITE'):
        filename = core.get_writer_file_prefix(mol.name()) + ".vibrec"
        with open(filename, 'w') as handle:
            json.dump(vibrec, handle, sort_keys=True, indent=4)

    if core.get_option('FINDIF', 'NORMAL_MODES_WRITE'):
        filename = core.get_writer_file_prefix(mol.name()) + ".molden_normal_modes"
        with open(filename, 'w') as handle:
            handle.write(qcdb.vib.print_molden_vibs(vibinfo, symbols, geom, standalone=True))

    return vibinfo


def gdma(wfn, datafile=""):
    """Function to use wavefunction information in *wfn* and, if specified,
    additional commands in *filename* to run GDMA analysis.

    .. versionadded:: 0.6

    :returns: None

    :type wfn: :py:class:`~psi4.core.Wavefunction`
    :param wfn: set of molecule, basis, orbitals from which to generate DMA analysis

    :type datafile: str
    :param datafile: optional control file (see GDMA manual) to peform more complicated DMA
                     analyses.  If this option is used, the File keyword must be set to read
                     a filename.fchk, where filename is provided by :term:`WRITER_FILE_LABEL <WRITER_FILE_LABEL (GLOBALS)>` .

    :examples:

    >>> # [1] DMA analysis from MP2 wavefunction.  N.B. gradient must be requested to generate MP2 density.
    >>> grad, wfn = gradient('mp2', return_wfn=True)
    >>> gdma(wfn)

    """
    # Start by writing a G* checkpoint file, for the GDMA code to read in
    fw = core.FCHKWriter(wfn)
    molname = wfn.molecule().name()
    prefix = core.get_writer_file_prefix(molname)
    fchkfile = prefix + '.fchk'
    fw.write(fchkfile)

    if datafile:
        commands = datafile
    else:
        if wfn.reference_wavefunction():
            densname = "CC"
        else:
            densname = "SCF"
        commands = 'psi4_dma_datafile.dma'
        radii = core.get_option('GDMA', 'GDMA_RADIUS')
        origin = core.get_option('GDMA', 'GDMA_ORIGIN')
        with open(commands, 'w') as f:
            f.write("File %s Density %s\n" % (fchkfile, densname))
            f.write("Angstrom\n")
            f.write("%s\n" % core.get_option('GDMA', 'GDMA_MULTIPOLE_UNITS'))
            f.write("Multipoles\n")
            if origin:
                try:
                    f.write("Origin %f %f %f\n" % (float(origin[0]), float(origin[1]), float(origin[2])))
                except IndexError:
                    raise ValidationError("The GDMA origin array should contain three entries: x, y, and z.")
            f.write("Switch %f\n" % core.get_option('GDMA', 'GDMA_SWITCH'))
            if radii:
                f.write("Radius %s\n" % " ".join([str(r) for r in radii]))
            f.write("Limit %d\n" % core.get_option('GDMA', 'GDMA_LIMIT'))
            f.write("Start\n")
            f.write("Finish\n")
    core.run_gdma(wfn, commands)

    os.remove(fchkfile)
    # If we generated the DMA control file, we should clean up here
    if not datafile:
        os.remove(commands)


def fchk(wfn: core.Wavefunction, filename: str, *, debug: bool = False, strict_label: bool = True):
    """Function to write wavefunction information in *wfn* to *filename* in
    Gaussian FCHK format.

    .. versionadded:: 0.6

    :returns: None

    :param wfn: set of molecule, basis, orbitals from which to generate fchk file

    :param filename: destination file name for FCHK file

    :param debug: returns a dictionary to aid with debugging

    :param strict_label: If true set a density label compliant with what Gaussian would write. A warning will be printed if this is not possible.
                         Otherwise set the density label according to the method name.

    Notes
    -----
    * A description of the FCHK format is http://wild.life.nctu.edu.tw/~jsyu/compchem/g09/g09ur/f_formchk.htm
    * The allowed headers for methods are general and limited, i.e., "Total SCF|MP2|CI|CC Density",
      PSI4 will try to find the right one for the current calculation. If `strict_label=False` the PSI4 method name will be used as label.
    * Not all theory modules in PSI4 are compatible with the FCHK writer.
      A warning will be printed if a theory module is not supported.
    * Caution! For orbital-optimized correlated methods (e.g. DCT, OMP2) the 'Orbital Energy' field contains ambiguous data.

    :examples:

    >>> # [1] FCHK file for DFT calculation
    >>> E, wfn = energy('b3lyp', return_wfn=True)
    >>> fchk(wfn, 'mycalc.fchk')

    >>> # [2] FCHK file for correlated densities
    >>> E, wfn = gradient('ccsd', return_wfn=True)
    >>> fchk(wfn, 'mycalc.fchk')

    >>> # [2] Write FCHK file with non-standard label.
    >>> E, wfn = gradient('mp2.5', return_wfn=True)
    >>> fchk(wfn, 'mycalc.fchk', strict_label=False)

    """
    # * Known limitations and notes *
    #
    # OCC: (occ theory module only, not dfocc) is turned off as densities are not correctly set.
    # DFMP2: Contains natural orbitals in wfn.C() and wfn.epsilon() data. This is fixed to contain respective HF data.

    allowed = ['DFMP2', 'SCF', 'CCENERGY', 'DCT', 'DFOCC']
    module_ = wfn.module().upper()
    if module_ not in allowed:
        core.print_out(f"FCHKWriter: Theory module {module_} is currently not supported by the FCHK writer.")
        return None

    if (wfn.basisset().has_ECP()):
        core.print_out(f"FCHKWriter: Limited ECP support! No ECP data will be written to the FCHK file.")

    # fix orbital coefficients and energies for DFMP2
    if module_ in ['DFMP2']:
        wfn_ = core.Wavefunction.build(wfn.molecule(), core.get_global_option('BASIS'))
        wfn_.deep_copy(wfn)
        refwfn = wfn.reference_wavefunction()
        wfn_.set_reference_wavefunction(refwfn)  # refwfn not deep_copied
        wfn_.Ca().copy(refwfn.Ca())
        wfn_.Cb().copy(refwfn.Cb())
        wfn_.epsilon_a().copy(refwfn.epsilon_a())
        wfn_.epsilon_b().copy(refwfn.epsilon_b())
        fw = core.FCHKWriter(wfn_)
    else:
        fw = core.FCHKWriter(wfn)

    if module_ in ['DCT', 'DFOCC']:
        core.print_out("""FCHKWriter: Caution! For orbital-optimized correlated methods
            the 'Orbital Energy' field contains ambiguous data. \n""")

    # At this point we don't know the method name, so we try to search for it.
    # idea: get the method from the variable matching closely the 'current energy'
    # for varlist, wfn is long-term and to allow from-file wfns. core is b/c some modules not storing in wfn yet
    varlist = {**wfn.scalar_variables(), **core.scalar_variables()}
    current = varlist['CURRENT ENERGY']

    # delete problematic entries
    for key in ['CURRENT ENERGY', 'CURRENT REFERENCE ENERGY']:
        varlist.pop(key, None)

    # find closest matching energy
    for (key, val) in varlist.items():
        if (np.isclose(val, current, 1e-12)):
            method = key.split()[0]
            break

    # The 'official' list of labels for compatibility.
    # OMP2,MP2.5,OCCD, etc get reduced to MP2,CC.
    allowed_labels = {
        "HF": " SCF Density",
        "SCF": " SCF Density",
        "DFT": " SCF Density",
        "MP2": " MP2 Density",
        "MP3": " MP3 Density",
        "MP4": " MP4 Density",
        "CI": " CI Density",
        "CC": " CC Density",
    }
    # assign label from method name
    fchk_label = f" {method} Density"
    if strict_label:
        in_list = False
        for key in allowed_labels:
            if key in method:
                if key is not method:
                    core.print_out(f"FCHKWriter: !WARNING! method '{method}'' renamed to label '{key}'.\n")
                fchk_label = allowed_labels[key]
                in_list = True
        if not in_list:
            core.print_out(f"FCHKWriter: !WARNING! {method} is not recognized. Using non-standard label.\n")
    core.print_out(f"FCHKWriter: Writing {filename} with label '{fchk_label}'.\n")
    fw.set_postscf_density_label(fchk_label)

    fw.write(filename)
    # needed for the pytest. The SCF density below follows PSI4 ordering not FCHK ordering.
    if debug:
        ret = {
            "filename": filename,
            "detected energy": method,
            "selected label": fchk_label,
            "Total SCF Density": fw.SCF_Dtot().np,
        }
        return ret
    return None

def molden(wfn, filename=None, density_a=None, density_b=None, dovirtual=None):
    """Function to write wavefunction information in *wfn* to *filename* in
    molden format. Will write natural orbitals from *density* (MO basis) if supplied.
    Warning! Most post-SCF Wavefunctions do not build the density as this is often
    much more costly than the energy. In addition, the Wavefunction density attributes
    (Da and Db) return the SO density and must be transformed to the MO basis
    to use with this function.

    .. versionadded:: 0.5
       *wfn* parameter passed explicitly

    :returns: None

    :type wfn: :py:class:`~psi4.core.Wavefunction`
    :param wfn: set of molecule, basis, orbitals from which to generate cube files

    :type filename: str
    :param filename: destination file name for MOLDEN file (optional)

    :type density_a: :py:class:`~psi4.core.Matrix`
    :param density_a: density in the MO basis to build alpha NO's from (optional)

    :type density_b: :py:class:`~psi4.core.Matrix`
    :param density_b: density in the MO basis to build beta NO's from, assumes restricted if not supplied (optional)

    :type dovirtual: bool
    :param dovirtual: do write all the MOs to the MOLDEN file (true) or discard the unoccupied MOs, not valid for NO's (false) (optional)

    :examples:

    1. Molden file with the Kohn-Sham orbitals of a DFT calculation.

       >>> E, wfn = energy('b3lyp', return_wfn=True)
       >>> molden(wfn, 'mycalc.molden')

    2. Molden file for CI/MCSCF computation using NO roots.
       Any method returning a ``CIWavefunction`` object will work: ``detci``,
       ``fci``, ``casscf``, etc. The first two arguments of ``get_opdm`` can be
       set to ``n, n`` where n => 0 selects the root to write out, provided
       these roots were computed, see :term:`NUM_ROOTS <NUM_ROOTS (DETCI)>`. The
       third argument controls the spin (``"A"``, ``"B"`` or ``"SUM"``) and the final
       boolean option determines whether inactive orbitals are included.

       >>> E, wfn = energy('detci', return_wfn=True)
       >>> molden(wfn, 'no_root1.molden', density_a=wfn.get_opdm(0, 0, "A", True))

    3. The following produces **an INCORRECT Molden file**, because the
       ``molden`` function needs orbitals in the MO basis (which are internally
       converted and written to the Molden file in the AO basis). The correct
       usage is given in the next point.

       >>> E, wfn = energy('ccsd', return_wfn=True)
       >>> molden(wfn, 'ccsd_no.molden', density_a=wfn.Da())

    4. Molden file with the natural orbitals of the ground-state 1RDM of a
       Post-HF calculation. Note the required transformation of Da (SO->MO).

       >>> E, wfn = properties('ccsd', return_wfn=True)
       >>> Da_so = wfn.Da()
       >>> SCa = core.doublet(wfn.S(), wfn.Ca(), False, False)
       >>> Da_mo = core.triplet(SCa, Da_so, SCa, True, False, False)
       >>> molden(wfn, 'ccsd_no.molden', density_a=Da_mo)

    """

    if filename is None:
        filename = core.get_writer_file_prefix(wfn.molecule().name()) + ".molden"

    if dovirtual is None:
        dovirt = bool(core.get_option("SCF", "MOLDEN_WITH_VIRTUAL"))

    else:
        dovirt = dovirtual

    if density_a:
        nmopi = wfn.nmopi()
        nsopi = wfn.nsopi()

        NO_Ra = core.Matrix("NO Alpha Rotation Matrix", nmopi, nmopi)
        NO_occa = core.Vector(nmopi)
        density_a.diagonalize(NO_Ra, NO_occa, core.DiagonalizeOrder.Descending)
        NO_Ca = core.Matrix("Ca Natural Orbitals", nsopi, nmopi)
        NO_Ca.gemm(False, False, 1.0, wfn.Ca(), NO_Ra, 0)

        if density_b:
            NO_Rb = core.Matrix("NO Beta Rotation Matrix", nmopi, nmopi)
            NO_occb = core.Vector(nmopi)
            density_b.diagonalize(NO_Rb, NO_occb, core.DiagonalizeOrder.Descending)
            NO_Cb = core.Matrix("Cb Natural Orbitals", nsopi, nmopi)
            NO_Cb.gemm(False, False, 1.0, wfn.Cb(), NO_Rb, 0)

        else:
            NO_occb = NO_occa
            NO_Cb = NO_Ca

        mw = core.MoldenWriter(wfn)
        mw.write(filename, NO_Ca, NO_Cb, NO_occa, NO_occb, NO_occa, NO_occb, dovirt)

    else:
        try:
            occa = wfn.occupation_a()
            occb = wfn.occupation_b()
        except AttributeError:
            core.print_out("\n!Molden warning: This wavefunction does not have occupation numbers.\n"
                           "Writing zero's for occupation numbers\n\n")
            occa = core.Vector(wfn.nmopi())
            occb = core.Vector(wfn.nmopi())

        mw = core.MoldenWriter(wfn)
        mw.write(filename, wfn.Ca(), wfn.Cb(), wfn.epsilon_a(), wfn.epsilon_b(), occa, occb, dovirt)

def tdscf(wfn, **kwargs):
    return proc.run_tdscf_excitations(wfn,**kwargs)


# Aliases
opt = optimize
freq = frequency
frequencies = frequency
prop = properties
