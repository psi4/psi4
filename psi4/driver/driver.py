#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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
"""Module with a *procedures* dictionary specifying available quantum
chemical methods and functions driving the main quantum chemical
functionality, namely single-point energies, geometry optimizations,
properties, and vibrational frequency calculations.

"""
import os
import re
import sys
import json
import shutil

import numpy as np

from psi4.driver import driver_util
from psi4.driver import driver_cbs
from psi4.driver import driver_nbody
from psi4.driver import driver_findif
from psi4.driver import p4util
from psi4.driver import qcdb
from psi4.driver.procrouting import *
from psi4.driver.p4util.exceptions import *

# never import wrappers or aliases into this file


def _find_derivative_type(ptype, method_name, user_dertype):
    r"""
    Figures out the derivative type (0, 1, 2) for a given method_name. Will
    first use user default and then the highest available derivative type for
    a given method.
    """

    if ptype not in ['gradient', 'hessian']:
        raise ValidationError("_find_derivative_type: ptype must either be gradient or hessian.")

    dertype = "(auto)"

    # If user type is None, try to find the highest derivative
    if user_dertype is None:
        if (ptype == 'hessian') and (method_name in procedures['hessian']):
            dertype = 2
            # Will need special logic if we ever have managed Hessians
        elif method_name in procedures['gradient']:
            dertype = 1
            if procedures['gradient'][method_name].__name__.startswith('select_'):
                try:
                    procedures['gradient'][method_name](method_name, probe=True)
                except ManagedMethodError:
                    dertype = 0
        elif method_name in procedures['energy']:
            dertype = 0
    else:
        # Quick sanity check. Only *should* be able to be None or int, but hey, kids today...
        if not isinstance(user_dertype, int):
            raise ValidationError("_find_derivative_type: user_dertype should only be None or int!")
        dertype = user_dertype

    if (core.get_global_option('INTEGRAL_PACKAGE') == 'ERD') and (dertype != 0):
        raise ValidationError('INTEGRAL_PACKAGE ERD does not play nicely with derivatives, so stopping.')

    if (core.get_global_option('PCM')) and (dertype != 0):
        core.print_out('\nPCM analytic gradients are not implemented yet, re-routing to finite differences.\n')
        dertype = 0

    # Summary validation
    if (dertype == 2) and (method_name in procedures['hessian']):
        pass
    elif (dertype == 1) and (method_name in procedures['gradient']):
        pass
    elif (dertype == 0) and (method_name in procedures['energy']):
        pass
    else:
        alternatives = ''
        alt_method_name = p4util.text.find_approximate_string_matches(method_name, procedures['energy'].keys(), 2)
        if len(alt_method_name) > 0:
            alternatives = """ Did you mean? %s""" % (' '.join(alt_method_name))

        raise ValidationError("""Derivative method 'name' %s and derivative level 'dertype' %s are not available.%s"""
                              % (method_name, str(dertype), alternatives))

    return dertype


def _energy_is_invariant(gradient, stationary_criterion=1.e-2):
    """Polls options and probes `gradient` to return whether current method
    and system expected to be invariant to translations and rotations of
    the coordinate system.

    """
    stationary_point = gradient.rms() < stationary_criterion  # 1.e-2 pulled out of a hat

    mol = core.get_active_molecule()
    efp_present = hasattr(mol, 'EFP')

    translations_projection_sound = (not core.get_option('SCF', 'EXTERN') and not core.get_option('SCF', 'PERTURB_H')
                                     and not efp_present)
    rotations_projection_sound = (translations_projection_sound and stationary_point)

    return translations_projection_sound, rotations_projection_sound


def _process_displacement(derivfunc, method, molecule, displacement, n, ndisp, **kwargs):
    """A helper function to perform all processing for an individual finite
       difference computation.

       Parameters
       ----------
       derivfunc : func
           The function computing the target derivative.
       method : str
          A string specifying the method to be used for the computation.
       molecule: psi4.core.molecule or qcdb.molecule
          The molecule for the computation. All processing is handled internally.
          molecule must not be modified!
       displacement : dict
          A dictionary containing the necessary information for the displacement.
          See driver_findif/_geom_generator.py docstring for details.
       n : int
          The number of the displacement being computed, for print purposes.
       ndisp : int
           The total number of geometries, for print purposes.

       Returns
       -------
       wfn: :py:class:`~psi4.core.Wavefunction`
           The wavefunction computed.
    """

    # print progress to file and screen
    core.print_out('\n')
    p4util.banner('Loading displacement %d of %d' % (n, ndisp))
    print(""" %d""" % (n), end=('\n' if (n == ndisp) else ''))
    sys.stdout.flush()

    parent_group = molecule.point_group()
    clone = molecule.clone()
    clone.reinterpret_coordentry(False)
    clone.fix_orientation(True)

    # Load in displacement (flat list) into the active molecule
    geom_array = np.reshape(displacement["geometry"], (-1, 3))
    clone.set_geometry(core.Matrix.from_array(geom_array))

    # If the user insists on symmetry, weaken it if some is lost when displacing.
    if molecule.symmetry_from_input():
        disp_group = clone.find_highest_point_group()
        new_bits = parent_group.bits() & disp_group.bits()
        new_symm_string = qcdb.PointGroup.bits_to_full_name(new_bits)
        clone.reset_point_group(new_symm_string)

    # clean possibly necessary for n=1 if its irrep (unsorted in displacement list) different from initial G0 for freq
    core.clean()

    # Perform the derivative calculation
    derivative, wfn = derivfunc(method, return_wfn=True, molecule=clone, **kwargs)
    displacement["energy"] = core.variable('CURRENT ENERGY')

    # If we computed a first or higher order derivative, set it.
    if derivfunc == gradient:
        displacement["gradient"] = wfn.gradient().np.ravel().tolist()

    # clean may be necessary when changing irreps of displacements
    core.clean()

    return wfn


def energy(name, **kwargs):
    r"""Function to compute the single-point electronic energy.

    :returns: *float* |w--w| Total electronic energy in Hartrees. SAPT & EFP return interaction energy.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| energy and wavefunction when **return_wfn** specified.

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

        Indicate to additionally return the :py:class:`~psi4.core.Wavefunction`
        calculation result as the second element (after *float* energy) of a tuple.

    :type restart_file: string
    :param restart_file: ``['file.1, file.32]`` || ``./file`` || etc.

        Binary data files to be renamed for calculation restart.

    .. _`table:energy_gen`:

    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | name                    | calls method                                                                                                  |
    +=========================+===============================================================================================================+
    | efp                     | effective fragment potential (EFP) :ref:`[manual] <sec:libefp>`                                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT) :ref:`[manual] <sec:scf>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | hf                      | HF self consistent field (SCF) :ref:`[manual] <sec:scf>`                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | hf3c                    | HF with dispersion, BSSE, and basis set corrections :ref:`[manual] <sec:gcp>`                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | pbeh3c                  | PBEh with dispersion, BSSE, and basis set corrections :ref:`[manual] <sec:gcp>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dcft                    | density cumulant functional theory :ref:`[manual] <sec:dcft>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp2                     | 2nd-order |MollerPlesset| perturbation theory (MP2) :ref:`[manual] <sec:dfmp2>` :ref:`[details] <tlmp2>`      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp3                     | 3rd-order |MollerPlesset| perturbation theory (MP3) :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tlmp3>`  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-mp3                 | MP3 with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp2.5                   | average of MP2 and MP3 :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tlmp25>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp4(sdq)                | 4th-order MP perturbation theory (MP4) less triples :ref:`[manual] <sec:fnompn>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-mp4(sdq)            | MP4 (less triples) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp4                     | full MP4 :ref:`[manual] <sec:fnompn>` :ref:`[details] <tlmp4>`                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-mp4                 | full MP4 with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp\ *n*                 | *n*\ th-order |MollerPlesset| (MP) perturbation theory :ref:`[manual] <sec:arbpt>`                            |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | zapt\ *n*               | *n*\ th-order z-averaged perturbation theory (ZAPT) :ref:`[manual] <sec:arbpt>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp2                    | orbital-optimized second-order MP perturbation theory :ref:`[manual] <sec:occ_oo>`                            |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs-omp2                | spin-component scaled OMP2 :ref:`[manual] <sec:occ_oo>`                                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs(n)-omp2             | a special version of SCS-OMP2 for nucleobase interactions :ref:`[manual] <sec:occ_oo>`                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs-omp2-vdw            | a special version of SCS-OMP2 (from ethene dimers) :ref:`[manual] <sec:occ_oo>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sos-omp2                | spin-opposite scaled OMP2 :ref:`[manual] <sec:occ_oo>`                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sos-pi-omp2             | A special version of SOS-OMP2 for pi systems :ref:`[manual] <sec:occ_oo>`                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp3                    | orbital-optimized third-order MP perturbation theory :ref:`[manual] <sec:occ_oo>`                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs-omp3                | spin-component scaled OMP3 :ref:`[manual] <sec:occ_oo>`                                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs(n)-omp3             | a special version of SCS-OMP3 for nucleobase interactions :ref:`[manual] <sec:occ_oo>`                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs-omp3-vdw            | a special version of SCS-OMP3 (from ethene dimers) :ref:`[manual] <sec:occ_oo>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sos-omp3                | spin-opposite scaled OMP3 :ref:`[manual] <sec:occ_oo>`                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sos-pi-omp3             | A special version of SOS-OMP3 for pi systems :ref:`[manual] <sec:occ_oo>`                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp2.5                  | orbital-optimized MP2.5 :ref:`[manual] <sec:occ_oo>`                                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | lccsd, cepa(0)          | coupled electron pair approximation variant 0 :ref:`[manual] <sec:fnocepa>` :ref:`[details] <tllccsd>`        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-lccsd, fno-cepa(0)  | CEPA(0) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cepa(1)                 | coupled electron pair approximation variant 1 :ref:`[manual] <sec:fnocepa>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-cepa(1)             | CEPA(1) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cepa(3)                 | coupled electron pair approximation variant 3 :ref:`[manual] <sec:fnocepa>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-cepa(3)             | CEPA(3) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | acpf                    | averaged coupled-pair functional :ref:`[manual] <sec:fnocepa>`                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-acpf                | ACPF with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | aqcc                    | averaged quadratic coupled cluster :ref:`[manual] <sec:fnocepa>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-aqcc                | AQCC with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | qcisd                   | quadratic CI singles doubles (QCISD) :ref:`[manual] <sec:fnocc>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-qcisd               | QCISD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | lccd                    | Linear CCD :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tllccd>`                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-lccd                | LCCD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | olccd                   | orbital optimized LCCD :ref:`[manual] <sec:occ_oo>`                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cc2                     | approximate coupled cluster singles and doubles (CC2) :ref:`[manual] <sec:cc>`                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccd                     | coupled cluster doubles  (CCD) :ref:`[manual] <sec:occ_nonoo>`                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsd>`                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | bccd                    | Brueckner coupled cluster doubles (BCCD) :ref:`[manual] <sec:cc>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-ccsd                | CCSD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | qcisd(t)                | QCISD with perturbative triples :ref:`[manual] <sec:fnocc>`                                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-qcisd(t)            | QCISD(T) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd(t)                 | CCSD with perturbative triples (CCSD(T)) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsdt>`                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd(at)                | CCSD with asymmetric perturbative triples (CCSD(AT)) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsdat>`     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | bccd(t)                 | BCCD with perturbative triples :ref:`[manual] <sec:cc>`                                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-ccsd(t)             | CCSD(T) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cc3                     | approximate CC singles, doubles, and triples (CC3) :ref:`[manual] <sec:cc>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccenergy                | **expert** full control over ccenergy module                                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dfocc                   | **expert** full control over dfocc module                                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cisd                    | configuration interaction (CI) singles and doubles (CISD) :ref:`[manual] <sec:ci>` :ref:`[details] <tlcisd>`  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-cisd                | CISD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cisdt                   | CI singles, doubles, and triples (CISDT) :ref:`[manual] <sec:ci>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cisdtq                  | CI singles, doubles, triples, and quadruples (CISDTQ) :ref:`[manual] <sec:ci>`                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ci\ *n*                 | *n*\ th-order CI :ref:`[manual] <sec:ci>`                                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fci                     | full configuration interaction (FCI) :ref:`[manual] <sec:ci>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | detci                   | **expert** full control over detci module                                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | casscf                  | complete active space self consistent field (CASSCF)  :ref:`[manual] <sec:ci>`                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | rasscf                  | restricted active space self consistent field (RASSCF)  :ref:`[manual] <sec:ci>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mcscf                   | multiconfigurational self consistent field (SCF) :ref:`[manual] <sec:psimrcc>`                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | psimrcc                 | Mukherjee multireference coupled cluster (Mk-MRCC) :ref:`[manual] <sec:psimrcc>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dmrg-scf                | density matrix renormalization group SCF :ref:`[manual] <sec:chemps2>`                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dmrg-caspt2             | density matrix renormalization group CASPT2 :ref:`[manual] <sec:chemps2>`                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dmrg-ci                 | density matrix renormalization group CI :ref:`[manual] <sec:chemps2>`                                         |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt0                   | 0th-order symmetry adapted perturbation theory (SAPT) :ref:`[manual] <sec:sapt>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ssapt0                  | 0th-order SAPT with special exchange scaling :ref:`[manual] <sec:sapt>`                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fisapt0                 | 0th-order functional and/or intramolecular SAPT :ref:`[manual] <sec:fisapt>`                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2                   | 2nd-order SAPT, traditional definition :ref:`[manual] <sec:sapt>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+                  | SAPT including all 2nd-order terms :ref:`[manual] <sec:sapt>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)               | SAPT including perturbative triples :ref:`[manual] <sec:sapt>`                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3                 | SAPT including all 3rd-order terms :ref:`[manual] <sec:sapt>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(ccd)             | SAPT2+ with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                    |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)          | SAPT2+(3) with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3(ccd)            | SAPT2+3 with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+dmp2              | SAPT including all 2nd-order terms and MP2 correction :ref:`[manual] <sec:sapt>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)dmp2           | SAPT including perturbative triples and MP2 correction :ref:`[manual] <sec:sapt>`                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3dmp2             | SAPT including all 3rd-order terms and MP2 correction :ref:`[manual] <sec:sapt>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(ccd)dmp2         | SAPT2+ with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)dmp2      | SAPT2+(3) with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3(ccd)dmp2        | SAPT2+3 with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt0-ct                | 0th-order SAPT plus charge transfer (CT) calculation :ref:`[manual] <sec:saptct>`                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2-ct                | SAPT2 plus CT :ref:`[manual] <sec:saptct>`                                                                    |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+-ct               | SAPT2+ plus CT :ref:`[manual] <sec:saptct>`                                                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)-ct            | SAPT2+(3) plus CT :ref:`[manual] <sec:saptct>`                                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3-ct              | SAPT2+3 plus CT :ref:`[manual] <sec:saptct>`                                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(ccd)-ct          | SAPT2+(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)-ct       | SAPT2+(3)(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3(ccd)-ct         | SAPT2+3(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | adc                     | 2nd-order algebraic diagrammatic construction (ADC) :ref:`[manual] <sec:adc>`                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | eom-cc2                 | EOM-CC2 :ref:`[manual] <sec:eomcc>`                                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | eom-ccsd                | equation of motion (EOM) CCSD :ref:`[manual] <sec:eomcc>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | eom-cc3                 | EOM-CC3 :ref:`[manual] <sec:eomcc>`                                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+

    .. comment missing and why
    .. comment a certain isapt --- marginally released
    .. comment mrcc --- this is handled in its own table
    .. comment psimrcc_scf --- convenience fn

    .. include:: ../autodoc_dft_energy.rst

    .. include:: ../mrcc_table_energy.rst

    .. include:: ../cfour_table_energy.rst

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

    # Bounce to CP if bsse kwarg
    if kwargs.get('bsse_type', None) is not None:
        return driver_nbody.nbody_gufunc(energy, name, ptype='energy', **kwargs)

    # Bounce if name is function
    if hasattr(name, '__call__'):
        return name(energy, kwargs.pop('label', 'custom function'), ptype='energy', **kwargs)

    # Allow specification of methods to arbitrary order
    lowername = name.lower()
    lowername, level = driver_util.parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    # Bounce to CBS if "method/basis" name
    if "/" in lowername:
        return driver_cbs._cbs_gufunc(energy, name, ptype='energy', **kwargs)

    # Commit to procedures['energy'] call hereafter
    return_wfn = kwargs.pop('return_wfn', False)
    core.clean_variables()

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    #for precallback in hooks['energy']['pre']:
    #    precallback(lowername, **kwargs)

    optstash = driver_util._set_convergence_criterion('energy', lowername, 6, 8, 6, 8, 6)

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
            name_split = re.split(r'\.', item)
            if "npz" in item:
                fname = os.path.split(os.path.abspath(core.get_writer_file_prefix(molecule.name())))[1]
                psi_scratch = core.IOManager.shared_object().get_default_path()
                file_num = item.split('.')[-2]
                targetfile = os.path.join(psi_scratch, fname + "." + file_num + ".npz")
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
            shutil.copy(item, targetfile)

    wfn = procedures['energy'][lowername](lowername, molecule=molecule, **kwargs)

    for postcallback in hooks['energy']['post']:
        postcallback(lowername, wfn=wfn, **kwargs)

    optstash.restore()
    if return_wfn:  # TODO current energy safer than wfn.energy() for now, but should be revisited

        # TODO place this with the associated call, very awkward to call this in other areas at the moment
        if lowername in ['efp', 'mrcc', 'dmrg', 'psimrcc']:
            core.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
            core.print_out("The returned wavefunction is the incoming reference wavefunction.\n\n")
        elif 'sapt' in lowername:
            core.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
            core.print_out("The returned wavefunction is the dimer SCF wavefunction.\n\n")

        return (core.variable('CURRENT ENERGY'), wfn)
    else:
        return core.variable('CURRENT ENERGY')


def gradient(name, **kwargs):
    r"""Function complementary to :py:func:~driver.optimize(). Carries out one gradient pass,
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
    kwargs = p4util.kwargs_lower(kwargs)

    # Figure out what kind of gradient this is
    if hasattr(name, '__call__'):
        if name.__name__ in ['cbs', 'complete_basis_set']:
            gradient_type = 'cbs_wrapper'
        else:
            # Bounce to name if name is non-CBS function
            gradient_type = 'custom_function'
    elif kwargs.get('bsse_type', None) is not None:
        gradient_type = 'nbody_gufunc'
    elif '/' in name:
        gradient_type = 'cbs_gufunc'
    else:
        gradient_type = 'conventional'

    # Figure out lowername, dertype, and func
    # If we have analytical gradients we want to pass to our wrappers, otherwise we want to run
    # finite-diference energy or cbs energies
    # TODO MP5/cc-pv[DT]Z behavior unkown due to "levels"
    user_dertype = kwargs.pop('dertype', None)
    if gradient_type == 'custom_function':
        if user_dertype is None:
            dertype = 0
            core.print_out(
                "\nGradient: Custom function passed in without a defined dertype, assuming fd-energy based gradient.\n"
            )
        else:
            core.print_out("\nGradient: Custom function passed in with a dertype of %d\n" % user_dertype)
            dertype = user_dertype

        if dertype == 1:
            return name(gradient, kwargs.pop('label', 'custom function'), ptype='gradient', **kwargs)
        else:
            optstash = driver_util._set_convergence_criterion('energy', 'scf', 8, 10, 8, 10, 8)
            lowername = name

    elif gradient_type == 'nbody_gufunc':
        return driver_nbody.nbody_gufunc(gradient, name, ptype='gradient', **kwargs)

    elif gradient_type == 'cbs_wrapper':
        cbs_methods = driver_cbs._cbs_wrapper_methods(**kwargs)
        dertype = min([_find_derivative_type('gradient', method, user_dertype) for method in cbs_methods])
        if dertype == 1:
            # Bounce to CBS (directly) in pure-gradient mode if name is CBS and all parts have analytic grad. avail.
            return name(gradient, kwargs.pop('label', 'custom function'), ptype='gradient', **kwargs)
        else:
            optstash = driver_util._set_convergence_criterion('energy', cbs_methods[0], 8, 10, 8, 10, 8)
            lowername = name
            # Pass through to G by E

    elif gradient_type == 'cbs_gufunc':
        cbs_methods = driver_cbs._parse_cbs_gufunc_string(name.lower())[0]
        dertype = min([_find_derivative_type('gradient', method, user_dertype) for method in cbs_methods])
        lowername = name.lower()
        if dertype == 1:
            # Bounce to CBS in pure-gradient mode if "method/basis" name and all parts have analytic grad. avail.
            return driver_cbs._cbs_gufunc(gradient, name, ptype='gradient', **kwargs)
        else:
            # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
            optstash = driver_util._set_convergence_criterion('energy', cbs_methods[0], 8, 10, 8, 10, 8)

    else:
        # Allow specification of methods to arbitrary order
        lowername = name.lower()
        lowername, level = driver_util.parse_arbitrary_order(lowername)
        if level:
            kwargs['level'] = level

        # Prevent methods that do not have associated gradients
        if lowername in energy_only_methods:
            raise ValidationError("gradient('%s') does not have an associated gradient" % name)

        dertype = _find_derivative_type('gradient', lowername, user_dertype)

        # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
        optstash = driver_util._set_convergence_criterion('energy', lowername, 8, 10, 8, 10, 8)

    # Commit to procedures[] call hereafter
    return_wfn = kwargs.pop('return_wfn', False)
    core.clean_variables()

    # no analytic derivatives for scf_type cd
    if core.get_global_option('SCF_TYPE') == 'CD':
        if (dertype == 1):
            raise ValidationError("""No analytic derivatives for SCF_TYPE CD.""")

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    # Does dertype indicate an analytic procedure both exists and is wanted?
    if dertype == 1:
        core.print_out("""gradient() will perform analytic gradient computation.\n""")

        # Perform the gradient calculation
        wfn = procedures['gradient'][lowername](lowername, molecule=molecule, **kwargs)

    else:
        core.print_out("""gradient() will perform gradient computation by finite difference of analytic energies.\n""")

        opt_iter = kwargs.get('opt_iter', 1)
        if opt_iter is True:
            opt_iter = 1

        if opt_iter == 1:
            print('Performing finite difference calculations')

        # Obtain list of displacements
        findif_meta_dict = driver_findif.gradient_from_energy_geometries(molecule)
        ndisp = len(findif_meta_dict["displacements"]) + 1

        print(""" %d displacements needed ...""" % (ndisp), end='')

        wfn = _process_displacement(energy, lowername, molecule, findif_meta_dict["reference"], 1, ndisp,
                                    **kwargs)
        var_dict = core.variables()

        for n, displacement in enumerate(findif_meta_dict["displacements"].values(), start=2):
            _process_displacement(
                energy, lowername, molecule, displacement, n, ndisp, write_orbitals=False, **kwargs)

        # Reset variables
        for key, val in var_dict.items():
            core.set_variable(key, val)

        # Compute the gradient
        core.set_local_option('FINDIF', 'GRADIENT_WRITE', True)
        G = driver_findif.compute_gradient_from_energies(findif_meta_dict)
        grad_psi_matrix = core.Matrix.from_array(G)
        grad_psi_matrix.print_out()
        wfn.set_gradient(grad_psi_matrix)

        # Explicitly set the current energy..
        core.set_variable('CURRENT ENERGY', findif_meta_dict["reference"]["energy"])

    optstash.restore()

    if core.get_option('FINDIF', 'GRADIENT_WRITE'):
        filename = core.get_writer_file_prefix(wfn.molecule().name()) + ".grad"
        qcdb.gradparse.to_string(np.asarray(wfn.gradient()), filename, dtype='GRD', mol=molecule, energy=wfn.energy())

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

    :type name: string
    :param name: ``'ccsd'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type properties: array of strings
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

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()
    kwargs['molecule'] = molecule

    # Allow specification of methods to arbitrary order
    lowername = args[0].lower()
    lowername, level = driver_util.parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    if "/" in lowername:
        return driver_cbs._cbs_gufunc(properties, lowername, ptype='properties', **kwargs)

    return_wfn = kwargs.pop('return_wfn', False)
    props = kwargs.get('properties', ['dipole', 'quadrupole'])

    if len(args) > 1:
        props += args[1:]

    kwargs['properties'] = p4util.drop_duplicates(props)
    optstash = driver_util._set_convergence_criterion('properties', lowername, 6, 10, 6, 10, 8)
    wfn = procedures['properties'][lowername](lowername, **kwargs)

    optstash.restore()

    if return_wfn:
        return (core.variable('CURRENT ENERGY'), wfn)
    else:
        return core.variable('CURRENT ENERGY')


def optimize(name, **kwargs):
    r"""Function to perform a geometry optimization.

    :aliases: opt()

    :returns: *float* |w--w| Total electronic energy of optimized structure in Hartrees.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| energy and wavefunction when **return_wfn** specified.

    :raises: psi4.OptimizationConvergenceError if |optking__geom_maxiter| exceeded without reaching geometry convergence.

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`CURRENT ENERGY <CURRENTENERGY>`

    :type name: string
    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the database. May be any valid argument to
        :py:func:`~driver.energy`.

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

    :type hessian_with: string
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
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT) :ref:`[manual] <sec:scf>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | hf                      | HF self consistent field (SCF) :ref:`[manual] <sec:scf>`                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dcft                    | density cumulant functional theory :ref:`[manual] <sec:dcft>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp2                     | 2nd-order |MollerPlesset| perturbation theory (MP2) :ref:`[manual] <sec:dfmp2>` :ref:`[details] <tlmp2>`      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp3                     | 3rd-order |MollerPlesset| perturbation theory (MP3) :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tlmp3>`  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp2.5                   | average of MP2 and MP3 :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tlmp25>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp2                    | orbital-optimized second-order MP perturbation theory :ref:`[manual] <sec:occ_oo>`                            |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp3                    | orbital-optimized third-order MP perturbation theory :ref:`[manual] <sec:occ_oo>`                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp2.5                  | orbital-optimized MP2.5 :ref:`[manual] <sec:occ_oo>`                                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | lccd                    | Linear CCD :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tllccd>`                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | olccd                   | orbital optimized LCCD :ref:`[manual] <sec:occ_oo>`                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccd                     | coupled cluster doubles  (CCD) :ref:`[manual] <sec:occ_nonoo>`                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsd>`                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd(t)                 | CCSD with perturbative triples (CCSD(T)) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsdt>`                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | eom-ccsd                | equation of motion (EOM) CCSD :ref:`[manual] <sec:eomcc>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+

    .. _`table:grad_scf`:


    .. include:: ../autodoc_dft_opt.rst

    .. include:: ../cfour_table_grad.rst


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

    if hasattr(name, '__call__'):
        lowername = name
        custom_gradient = True
    else:
        lowername = name.lower()
        custom_gradient = False

    return_wfn = kwargs.pop('return_wfn', False)

    return_history = kwargs.pop('return_history', False)
    if return_history:
        # Add wfn once the deep copy issues are worked out
        step_energies = []
        step_gradients = []
        step_coordinates = []

    # For CBS and nbody wrappers, need to set retention on INTCO file
    if custom_gradient or ('/' in lowername) or kwargs.get('bsse_type', None) is not None:
        core.IOManager.shared_object().set_specific_retention(1, True)

    full_hess_every = core.get_option('OPTKING', 'FULL_HESS_EVERY')
    steps_since_last_hessian = 0

    if custom_gradient and core.has_option_changed('OPTKING', 'FULL_HESS_EVERY'):
        raise ValidationError("Optimize: Does not support custom Hessian's yet.")
    else:
        hessian_with_method = kwargs.get('hessian_with', lowername)

    optstash = p4util.OptionsState(
        ['OPTKING', 'INTRAFRAG_STEP_LIMIT'],
        ['FINDIF', 'HESSIAN_WRITE'],
        ['OPTKING', 'CART_HESS_READ'],
        ['SCF', 'GUESS_PERSIST'],  # handle on behalf of cbs()
        ['SCF', 'GUESS'])

    n = kwargs.get('opt_iter', 1)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())

    # If we are freezing cartesian, do not orient or COM
    if core.get_local_option("OPTKING", "FROZEN_CARTESIAN"):
        molecule.fix_orientation(True)
        molecule.fix_com(True)
    molecule.update_geometry()

    # Shifting the geometry so need to copy the active molecule
    moleculeclone = molecule.clone()

    initial_sym = moleculeclone.schoenflies_symbol()
    while n <= core.get_option('OPTKING', 'GEOM_MAXITER'):
        current_sym = moleculeclone.schoenflies_symbol()
        if initial_sym != current_sym:
            raise ValidationError("""Point group changed! (%s <-- %s) You should restart """
                                  """using the last geometry in the output, after """
                                  """carefully making sure all symmetry-dependent """
                                  """input, such as DOCC, is correct.""" % (current_sym, initial_sym))
        kwargs['opt_iter'] = n

        # Use orbitals from previous iteration as a guess
        #   set within loop so that can be influenced by fns to optimize (e.g., cbs)
        if (n > 1) and (not core.get_option('SCF', 'GUESS_PERSIST')):
            core.set_local_option('SCF', 'GUESS', 'READ')

        # Before computing gradient, save previous molecule and wavefunction if this is an IRC optimization
        if (n > 1) and (core.get_option('OPTKING', 'OPT_TYPE') == 'IRC'):
            old_thisenergy = core.variable('CURRENT ENERGY')

        # Compute the gradient - preserve opt data despite core.clean calls in gradient
        core.IOManager.shared_object().set_specific_retention(1, True)
        G, wfn = gradient(lowername, return_wfn=True, molecule=moleculeclone, **kwargs)
        thisenergy = core.variable('CURRENT ENERGY')

        # above, used to be getting energy as last of energy list from gradient()
        # thisenergy below should ultimately be testing on wfn.energy()

        # Record optimization steps
        # Add wavefunctions later
        if return_history:
            step_energies.append(thisenergy)
            step_coordinates.append(moleculeclone.geometry())
            step_gradients.append(G.clone())

        core.set_legacy_gradient(G)

        # opt_func = kwargs.get('opt_func', kwargs.get('func', energy))
        # if opt_func.__name__ == 'complete_basis_set':
        #     core.IOManager.shared_object().set_specific_retention(1, True)

        if full_hess_every > -1:
            core.set_global_option('HESSIAN_WRITE', True)

        # compute Hessian as requested; frequency wipes out gradient so stash it
        if ((full_hess_every > -1) and (n == 1)) or (steps_since_last_hessian + 1 == full_hess_every):
            G = core.get_legacy_gradient()  # TODO
            core.IOManager.shared_object().set_specific_retention(1, True)
            core.IOManager.shared_object().set_specific_path(1, './')
            frequencies(hessian_with_method, molecule=moleculeclone, **kwargs)
            steps_since_last_hessian = 0
            core.set_legacy_gradient(G)
            core.set_global_option('CART_HESS_READ', True)
        elif (full_hess_every == -1) and core.get_global_option('CART_HESS_READ') and (n == 1):
            pass
            # Do nothing; user said to read existing hessian once
        else:
            core.set_global_option('CART_HESS_READ', False)
            steps_since_last_hessian += 1

        # Take step. communicate to/from/within optking through legacy_molecule
        core.set_legacy_molecule(moleculeclone)
        optking_rval = core.optking()
        moleculeclone = core.get_legacy_molecule()
        moleculeclone.update_geometry()
        if optking_rval == core.PsiReturnType.EndLoop:
            # if this is the end of an IRC run, set wfn, energy, and molecule to that
            # of the last optimized IRC point
            if core.get_option('OPTKING', 'OPT_TYPE') == 'IRC':
                thisenergy = old_thisenergy
            print('Optimizer: Optimization complete!')
            core.print_out('\n    Final optimized geometry and variables:\n')
            moleculeclone.print_in_input_format()
            # Mark the optimization data as disposable now that the optimization is done.
            core.IOManager.shared_object().set_specific_retention(1, False)
            # Check if user wants to see the intcos; if so, don't delete them.
            if core.get_option('OPTKING', 'INTCOS_GENERATE_EXIT') == False:
                if core.get_option('OPTKING', 'KEEP_INTCOS') == False:
                    core.opt_clean()
            # Changing environment to optimized geometry as expected by user
            molecule.set_geometry(moleculeclone.geometry())
            for postcallback in hooks['optimize']['post']:
                postcallback(lowername, wfn=wfn, **kwargs)
            core.clean()

            # Cleanup binary file 1
            if custom_gradient or ('/' in lowername) or kwargs.get('bsse_type', None) is not None:
                core.IOManager.shared_object().set_specific_retention(1, False)

            optstash.restore()

            if return_history:
                history = {
                    'energy': step_energies,
                    'gradient': step_gradients,
                    'coordinates': step_coordinates,
                }

            if return_wfn and return_history:
                return (thisenergy, wfn, history)
            elif return_wfn and not return_history:
                return (thisenergy, wfn)
            elif return_history and not return_wfn:
                return (thisenergy, history)
            else:
                return thisenergy

        elif optking_rval == core.PsiReturnType.Failure:
            print('Optimizer: Optimization failed!')
            # Mark the optimization data as disposable now that the optimization is done.
            core.IOManager.shared_object().set_specific_retention(1, False)
            if (core.get_option('OPTKING', 'KEEP_INTCOS') == False):
                core.opt_clean()
            molecule.set_geometry(moleculeclone.geometry())
            core.clean()
            optstash.restore()
            raise OptimizationConvergenceError("""geometry optimization""", n - 1, wfn)
            return thisenergy

        core.print_out('\n    Structure for next step:\n')
        moleculeclone.print_in_input_format()

        n += 1

    if core.get_option('OPTKING', 'INTCOS_GENERATE_EXIT') == False:
        if core.get_option('OPTKING', 'KEEP_INTCOS') == False:
            core.opt_clean()

    optstash.restore()
    raise OptimizationConvergenceError("""geometry optimization""", n - 1, wfn)


def hessian(name, **kwargs):
    r"""Function complementary to :py:func:`~frequency`. Computes force
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
    kwargs = p4util.kwargs_lower(kwargs)

    # Figure out what kind of gradient this is
    if hasattr(name, '__call__'):
        if name.__name__ in ['cbs', 'complete_basis_set']:
            gradient_type = 'cbs_wrapper'
        else:
            # Bounce to name if name is non-CBS function
            gradient_type = 'custom_function'

    elif kwargs.get('bsse_type', None) is not None:
        gradient_type = 'nbody_gufunc'
    elif '/' in name:
        gradient_type = 'cbs_gufunc'
    else:
        gradient_type = 'conventional'

    # Call appropriate wrappers
    if gradient_type == 'nbody_gufunc':
        return driver_nbody.nbody_gufunc(hessian, name.lower(), ptype='hessian', **kwargs)
    # Check if this is a CBS extrapolation
    elif gradient_type == "cbs_gufunc":
        return driver_cbs._cbs_gufunc(hessian, name.lower(), **kwargs, ptype="hessian")
    elif gradient_type == "cbs_wrapper":
        return driver_cbs.cbs(hessian, "cbs", **kwargs, ptype="hessian")
    elif gradient_type != "conventional":
        raise ValidationError("Hessian: Does not yet support custom functions.")
    else:
        lowername = name.lower()

    return_wfn = kwargs.pop('return_wfn', False)
    core.clean_variables()
    dertype = 2

    # Prevent methods that do not have associated energies
    if lowername in energy_only_methods:
        raise ValidationError("hessian('%s') does not have an associated hessian" % name)

    optstash = p4util.OptionsState(
        ['FINDIF', 'HESSIAN_WRITE'],
        ['FINDIF', 'FD_PROJECT'],
    )

    # Allow specification of methods to arbitrary order
    lowername, level = driver_util.parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    dertype = _find_derivative_type('hessian', lowername, kwargs.pop('freq_dertype', kwargs.pop('dertype', None)))

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
    optstash_conv = driver_util._set_convergence_criterion('energy', lowername, 8, 10, 8, 10, 8)

    # Select certain irreps
    irrep = kwargs.get('irrep', -1)
    if irrep == -1:
        pass  # do all irreps
    else:
        irrep = driver_util.parse_cotton_irreps(irrep, molecule.schoenflies_symbol())
        irrep -= 1  # A1 irrep is externally 1, internally 0
        if dertype == 2:
            core.print_out(
                """hessian() switching to finite difference by gradients for partial Hessian calculation.\n""")
            dertype = 1

    # At stationary point?
    if 'ref_gradient' in kwargs:
        core.print_out("""hessian() using ref_gradient to assess stationary point.\n""")
        G0 = kwargs['ref_gradient']
    else:
        G0 = gradient(lowername, molecule=molecule, **kwargs)
    translations_projection_sound, rotations_projection_sound = _energy_is_invariant(G0)
    core.print_out(
        '\n  Based on options and gradient (rms={:.2E}), recommend {}projecting translations and {}projecting rotations.\n'
        .format(G0.rms(), '' if translations_projection_sound else 'not ',
                '' if rotations_projection_sound else 'not '))
    if not core.has_option_changed('FINDIF', 'FD_PROJECT'):
        core.set_local_option('FINDIF', 'FD_PROJECT', rotations_projection_sound)

    # Does an analytic procedure exist for the requested method?
    if dertype == 2:
        core.print_out("""hessian() will perform analytic frequency computation.\n""")

        # We have the desired method. Do it.
        wfn = procedures['hessian'][lowername](lowername, molecule=molecule, **kwargs)
        wfn.set_gradient(G0)
        optstash.restore()
        optstash_conv.restore()

        # TODO: check that current energy's being set to the right figure when this code is actually used
        core.set_variable('CURRENT ENERGY', wfn.energy())

    elif dertype == 1:
        core.print_out(
            """hessian() will perform frequency computation by finite difference of analytic gradients.\n""")

        # Obtain list of displacements
        findif_meta_dict = driver_findif.hessian_from_gradient_geometries(molecule, irrep)

        # Record undisplaced symmetry for projection of displaced point groups
        core.set_parent_symmetry(molecule.schoenflies_symbol())

        ndisp = len(findif_meta_dict["displacements"]) + 1

        print(""" %d displacements needed.""" % ndisp)

        wfn = _process_displacement(gradient, lowername, molecule, findif_meta_dict["reference"], 1, ndisp,
                                    **kwargs)
        var_dict = core.variables()

        for n, displacement in enumerate(findif_meta_dict["displacements"].values(), start=2):
            _process_displacement(
                gradient, lowername, molecule, displacement, n, ndisp, write_orbitals=False, **kwargs)

        # Reset variables
        for key, val in var_dict.items():
            core.set_variable(key, val)

        # Assemble Hessian from gradients
        #   Final disp is undisp, so wfn has mol, G, H general to freq calc
        H = driver_findif.compute_hessian_from_gradients(findif_meta_dict, irrep)
        wfn.set_hessian(core.Matrix.from_array(H))
        wfn.set_gradient(G0)

        # Explicitly set the current energy..
        core.set_variable('CURRENT ENERGY', findif_meta_dict["reference"]["energy"])
        wfn.set_variable('CURRENT ENERGY', findif_meta_dict["reference"]["energy"])

        core.set_parent_symmetry('')
        optstash.restore()
        optstash_conv.restore()

    else:
        core.print_out("""hessian() will perform frequency computation by finite difference of analytic energies.\n""")

        # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
        optstash.restore()
        optstash_conv.restore()
        optstash_conv = driver_util._set_convergence_criterion('energy', lowername, 10, 11, 10, 11, 10)

        # Obtain list of displacements
        findif_meta_dict = driver_findif.hessian_from_energy_geometries(molecule, irrep)

        # Record undisplaced symmetry for projection of diplaced point groups
        core.set_parent_symmetry(molecule.schoenflies_symbol())

        ndisp = len(findif_meta_dict["displacements"]) + 1

        print(' %d displacements needed.' % ndisp)

        wfn = _process_displacement(energy, lowername, molecule, findif_meta_dict["reference"], 1, ndisp,
                                    **kwargs)
        var_dict = core.variables()

        for n, displacement in enumerate(findif_meta_dict["displacements"].values(), start=2):
            _process_displacement(
                energy, lowername, molecule, displacement, n, ndisp, write_orbitals=False, **kwargs)

        # Reset variables
        for key, val in var_dict.items():
            core.set_variable(key, val)

        # Assemble Hessian from energies
        H = driver_findif.compute_hessian_from_energies(findif_meta_dict, irrep)
        wfn.set_hessian(core.Matrix.from_array(H))
        wfn.set_gradient(G0)

        # Explicitly set the current energy..
        core.set_variable('CURRENT ENERGY', findif_meta_dict["reference"]["energy"])
        wfn.set_variable('CURRENT ENERGY', findif_meta_dict["reference"]["energy"])

        core.set_parent_symmetry('')
        optstash.restore()
        optstash_conv.restore()

    _hessian_write(wfn)

    if return_wfn:
        return (wfn.hessian(), wfn)
    else:
        return wfn.hessian()


def frequency(name, **kwargs):
    r"""Function to compute harmonic vibrational frequencies.

    :aliases: frequencies(), freq()

    :returns: *float* |w--w| Total electronic energy in Hartrees.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| energy and wavefunction when **return_wfn** specified.

    :type name: string
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

    :type irrep: int or string
    :param irrep: |dl| ``-1`` |dr| || ``1`` || ``'b2'`` || ``'App'`` || etc.

        Indicates which symmetry block (:ref:`Cotton <table:irrepOrdering>` ordering) of vibrational
        frequencies to be computed. ``1``, ``'1'``, or ``'a1'`` represents
        :math:`a_1`, requesting only the totally symmetric modes.
        ``-1`` indicates a full frequency calculation.

    .. note:: Analytic hessians are only available for RHF. For all other methods, Frequencies will
        proceed through finite differences according to availability of gradients or energies.

    .. _`table:freq_gen`:

    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | name                    | calls method                                                                                                  |
    +=========================+===============================================================================================================+
    | scf                     | Hartree--Fock (HF) :ref:`[manual] <sec:scf>`                                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+

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
    translations_projection_sound, rotations_projection_sound = _energy_is_invariant(wfn.gradient())
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


def vibanal_wfn(wfn, hess=None, irrep=None, molecule=None, project_trans=True, project_rot=True):
    """Function to perform analysis of a hessian or hessian block, specifically...
    calling for and printing vibrational and thermochemical analysis, setting thermochemical variables,
    and writing the vibrec and normal mode files.

    Parameters
    ----------
    wfn : :py:class:`~psi4.core.Wavefunction`
        The wavefunction which had its Hessian computed.
    hess : ndarray of float, optional
        Hessian to analyze, if not the hessian in wfn.
        (3*nat, 3*nat) non-mass-weighted Hessian in atomic units, [Eh/a0/a0].
    irrep : int or string
        The irrep for which frequencies are calculated. Thermochemical analysis is skipped if this is given,
        as only one symmetry block of the hessian has been computed.
    molecule : :py:class:`~psi4.core.Molecule` or qcdb.Molecule, optional
        The molecule to pull information from, if not the molecule in wfn. Must at least have similar
        geometry to the molecule in wfn.
    project_trans : boolean
        Should translations be projected in the harmonic analysis?
    project_rot : boolean
        Should rotations be projected in the harmonic analysis?

    Returns
    -------
    vibinfo : dict
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

    vibinfo, vibtext = qcdb.vib.harmonic_analysis(
        nmwhess, geom, m, wfn.basisset(), irrep_labels, dipder=dipder, project_trans=project_trans, project_rot=project_rot)
    vibrec.update({k: qca.to_dict() for k, qca in vibinfo.items()})

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
        vibrec.update({k: qca.to_dict() for k, qca in therminfo.items()})

        core.set_variable("ZPVE", therminfo['ZPE_corr'].data)
        core.set_variable("THERMAL ENERGY CORRECTION", therminfo['E_corr'].data)
        core.set_variable("ENTHALPY CORRECTION", therminfo['H_corr'].data)
        core.set_variable("GIBBS FREE ENERGY CORRECTION", therminfo['G_corr'].data)

        core.set_variable("ZERO K ENTHALPY", therminfo['ZPE_tot'].data)
        core.set_variable("THERMAL ENERGY", therminfo['E_tot'].data)
        core.set_variable("ENTHALPY", therminfo['H_tot'].data)
        core.set_variable("GIBBS FREE ENERGY", therminfo['G_tot'].data)

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


def _hessian_write(wfn):
    if core.get_option('FINDIF', 'HESSIAN_WRITE'):
        filename = core.get_writer_file_prefix(wfn.molecule().name()) + ".hess"
        with open(filename, 'wb') as handle:
            qcdb.hessparse.to_string(np.asarray(wfn.hessian()), handle, dtype='psi4')


def gdma(wfn, datafile=""):
    """Function to use wavefunction information in *wfn* and, if specified,
    additional commands in *filename* to run GDMA analysis.

    .. include:: ../autodoc_abbr_options_c.rst

    .. versionadded:: 0.6

    :returns: None

    :type wfn: :py:class:`~psi4.core.Wavefunction`
    :param wfn: set of molecule, basis, orbitals from which to generate DMA analysis

    :type datafile: string
    :param datafile: optional control file (see GDMA manual) to peform more complicated DMA
                     analyses.  If this option is used, the File keyword must be set to read
                     a filename.fchk, where filename is provided by |globals__writer_file_label| .

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
                except:
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


def fchk(wfn, filename):
    """Function to write wavefunction information in *wfn* to *filename* in
    Gaussian FCHK format.

    .. versionadded:: 0.6

    :returns: None

    :type filename: string
    :param filename: destination file name for FCHK file

    :type wfn: :py:class:`~psi4.core.Wavefunction`
    :param wfn: set of molecule, basis, orbitals from which to generate fchk file

    Notes
    -----
    * A description of the FCHK format is http://wild.life.nctu.edu.tw/~jsyu/compchem/g09/g09ur/f_formchk.htm
    * The allowed headers for methods are general and limited, i.e., "Total SCF|MP2|CI|CC Density",
      so "CC" is always used for the post-HF case.

    :examples:

    >>> # [1] FCHK file for DFT calculation
    >>> E, wfn = energy('b3lyp', return_wfn=True)
    >>> fchk(wfn, 'mycalc.fchk')

    """
    fw = core.FCHKWriter(wfn)
    fw.write(filename)


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

    :type filename: string
    :param filename: destination file name for MOLDEN file (optional)

    :type density_a: :py:class:`~psi4.core.Matrix`
    :param density_a: density in the MO basis to build alpha NO's from (optional)

    :type density_b: :py:class:`~psi4.core.Matrix`
    :param density_b: density in the MO basis to build beta NO's from, assumes restricted if not supplied (optional)

    :type dovirtual: bool
    :param dovirtual: do write all the MOs to the MOLDEN file (true) or discard the unoccupied MOs, not valid for NO's (false) (optional)

    :examples:

    >>> # [1] Molden file for DFT calculation
    >>> E, wfn = energy('b3lyp', return_wfn=True)
    >>> molden(wfn, 'mycalc.molden')

    >>> # [2] Molden file for CI/MCSCF computation using NO roots
    >>> E, wfn = energy('ci', return_wfn=True)
    >>> molden(wfn, 'no_root1.molden', density_a=wfn.opdm(0, 0, "A", True))

    >>> # [3] The following does NOT work, please see below
    >>> E, wfn = energy('ccsd', return_wfn=True)
    >>> molden(wfn, 'ccsd_no.molden', density_a=wfn.Da())

    >>> # [4] This WILL work, note the transformation of Da (SO->MO)
    >>> E, wfn = properties('ccsd', properties=['dipole'], return_wfn=True)
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


# Aliases
opt = optimize
freq = frequency
frequencies = frequency
prop = properties
