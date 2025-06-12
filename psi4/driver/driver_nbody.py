#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
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
"""Plan, run, and assemble QC tasks to obtain many-body expansion and basis-set superposition error treatments.

=============
ManyBody Flow
=============
Bullet points are major actions
Lines of dashes denote function calls
e/d/dd=dg/g/h := energy, dipole, dipole derivative = dipole gradient, gradient, Hessian
`mc_(frag, bas)` := a modelchem index, mc; indices of real fragments, frag; set(bas - frag) are indices of ghost fragments. see "intermediates_energy" in big table below for example.
note that there's a lot of natural 1-indexing (1, 2, 3) rather than 0-indexing (0, 1, 2) in manybody. e.g., 2-body energy, Molecule.extract_subsets(1, (1, 2))
note that a "level" can be n-body level (how many real molecular fragments) or a modelchem level (`mc_`; e.g., CC on 1-bodies, MP2 on 2-bodies; "multilevel")

---------------------------
ManyBodyComputer.__init__()
---------------------------
* not an explicit function but pydantic handles some defaults and validation
* fields molecule, nfragments, bsse_type, return_total_data, and initial max_nbody set
* BaseComputer.__init__()

task_planner.py::task_planner()
-------------------------------
* computer gets modified from task_planner outside this file!
* modelchem (method and basis) treatment levels for each n-body level determined from user levels kwarg. fields nbodies_per_mc_level set and max_nbody reset
* for each modelchem treatment level, call build_tasks() below via one of four routes, depending on simple MB or layered MB(FD), MB(CBS), or MB(FD(CBS))

    ------------------------------
    ManyBodyComputer.build_tasks()
    ------------------------------
    * if supersystem requested as a modelchem level, request (frag, bas) indices for full nbody range of nocp treatment from build_nbody_compute_list()
    * otherwise, request (frag, bas) indices for specified nbody range covering specified bsse treatments from build_nbody_compute_list()

        build_nbody_compute_list()
        --------------------------
        * initializes dicts for each of nocp, cp, vmfc (2 for this one) with keys requested n-body levels and values empty sets
        * use combinatorics formulas to fill each key with (frag, bas) indices (what fragments are active and what fragments have basis functions)
          needed to compute the requested bsse treatments at the requested n-body levels.
        * merge by n-body level the sets of indices for each bsse treatment into an "all" dict. return this and all the per-bsse dicts.

    * merge (from different bsse_types) all the requested indices, prepend a modelchem treatment index to form `mc_(frag, bas)`
    * construct a molecule appropriately real/ghosted from active-fragment info in (frag, bas)
    * if embedding_charges active, prepare external_potentials array for atoms not in bas fragments
    * for any new `mc_(frag, bas)` index, append a new computer to self.task_list

--------------------------
ManyBodyComputer.compute()
--------------------------
* compute() for each job in self.task_list

----------------------------------
ManyBodyComputer.get_psi_results()
----------------------------------

    Computer.get_results()
    ----------------------
    * call get_results() for each job in task list
    * assemble all the computed energies, all the computed gradients, and all the computed hessians
    * call qcmb's core interface analyze() to sum assembled component results into MBE results
    * call qcmb's high-lvl interface get_results() to transform MBE results into QCSchema model
    * return MBE QCSchema model

* collect ManyBodyResult schema from self.get_results() (prior to v1.10, this was a ManyBody-flavored AtomicResult
* collect text summary table from schema. Print and logs formatted energy output separately for cp, nocp, vmfc
* build wfn from nbody mol and basis (always def2-svp)
* collect MBE properties from schema, translate them to qcvars, and push to P::e and wfn
* collect subsystem components properties, and write them to qcvars
* convert result to psi4.core.Matrix (non-energy) and set g/h on wfn
* return e/g/h and wfn

"""

__all__ = [
    "BsseEnum",
    "ManyBodyComputer",
    "nbody",
]

import copy
import logging
import pprint
# v2: from typing import TYPE_CHECKING, Any, ClassVar, Dict, Tuple, Union, Optional
from typing import TYPE_CHECKING, Any, Dict, Tuple, Union, Optional

from psi4 import core

from . import p4util
from .constants import constants, pp
from .driver_cbs import CompositeComputer
from .driver_findif import FiniteDifferenceComputer
from .p4util.exceptions import *
from .task_base import AtomicComputer, EnergyGradientHessianWfnReturn

from pydantic.v1 import validator

from qcelemental.models import Molecule
from qcmanybody.models.v1 import AtomicSpecification, ManyBodyInput, ManyBodyResult, BsseEnum
from qcmanybody import ManyBodyCore, ManyBodyComputer as ManyBodyComputerQCNG
from qcmanybody.utils import delabeler, labeler, translate_qcvariables, modelchem_labels

if TYPE_CHECKING:
    import qcportal

logger = logging.getLogger(__name__)

FragBasIndex = Tuple[Tuple[int], Tuple[int]]

SubTaskComputers = Union[AtomicComputer, CompositeComputer, FiniteDifferenceComputer]


# nbody function is here for the docstring
def nbody():
    """
    Computes the nbody interaction energy, gradient, or Hessian depending on input.
    This is a generalized universal function for computing interaction and total quantities.

    :returns: *return type of func* |w--w| The data.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| data and wavefunction with energy/gradient/hessian set appropriately when **return_wfn** specified.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type return_wfn: :ref:`boolean <op_py_boolean>`
    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|

        Indicate to additionally return the :py:class:`~psi4.core.Wavefunction`
        calculation result as the second element of a tuple.

    :type bsse_type: str or list
    :param bsse_type: ``'cp'`` || ``['nocp', 'vmfc']`` || |dl| ``None`` |dr| || etc.

        Type of BSSE correction to compute: CP for counterpoise correction, NoCP
        for plain supramolecular interaction energy, or VMFC for Valiron-Mayer
        Function Counterpoise correction. If a list is provided, the first string in
        the list determines which interaction or total energies/gradients/Hessians are
        returned by this function. By default, many-body treatments are inactive.

    :type max_nbody: int
    :param max_nbody: ``3`` || etc.

        Maximum n-body to compute, cannot exceed the number of fragments in the molecule.

    :type return_total_data: :ref:`boolean <op_py_boolean>`
    :param return_total_data: ``'on'`` || |dl| ``'off'`` |dr|

        If True returns the total data (energy/gradient/Hessian) of the system,
        otherwise returns interaction data. Default is ``'off'`` for energies,
        ``'on'`` for gradients and Hessians. Note that the calculation of total
        counterpoise corrected energies implies the calculation of the energies of
        monomers in the monomer basis, hence specifying ``return_total_data = True``
        may carry out more computations than ``return_total_data = False``.
        For gradients and Hessians, ``return_total_data = False`` is rarely useful.

    :type supersystem_ie_only: :ref:`boolean <op_py_boolean>`
    :param supersystem_ie_only: ``'on'`` || |dl| ``'off'`` |dr|

        Target the supersystem total/interaction energy (IE) data over the many-body expansion (MBE)
        analysis, thereby omitting intermediate-body calculations. When False (default), compute each n-body level
        in the MBE up through ``max_nbody``. When True (only allowed for ``max_nbody = nfragments`` ), only compute
        enough for the overall interaction/total energy: max_nbody-body and 1-body. When True, properties
        ``INTERACTION {driver} THROUGH {max_nbody}-BODY`` will always be available;
        ``TOTAL {driver} THROUGH {max_nbody}-BODY`` will be available depending on ``return_total_data`` ; and
        ``{max_nbody}-BODY CONTRIBUTION TO {driver}`` won't be available (except for dimers). This keyword produces
        no savings for a two-fragment molecule. But for the interaction energy of a three-fragment molecule, for
        example, 2-body subsystems can be skipped with ``supersystem_ie_only=True``. Do not use with ``vmfc`` in
        ``bsse_type`` as it cannot produce savings.

    :type levels: dict
    :param levels: ``{1: 'ccsd(t)', 2: 'mp2', 'supersystem': 'scf'}`` || ``{1: 2, 2: 'ccsd(t)', 3: 'mp2'}`` || etc

        Dictionary of different levels of theory for different levels of expansion
        Note that method_string is not used in this case. ``supersystem`` computes
        all higher order n-body effects up to the number of fragments.

    :type embedding_charges: dict
    :param embedding_charges: ``{1: [-0.834, 0.417, 0.417], ..}``

        Dictionary of atom-centered point charges. keys: 1-based index of fragment, values: list of charges for each fragment.
        Add atom-centered point charges for fragments whose basis sets are not included in the computation.

    Potential QCVariables set are:

        ptype_size = (1,)/(nat, 3)/(3 * nat, 3 * nat)
        e/g/h := energy or gradient or Hessian
        rtd := return_total_data

        .. |em| unicode:: U+02003 .. em space

        .. _`table:nbody_return`:

        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | item                                                          | size                 | present / zeroed                                                   | contents / interpretation                                                                                          |
        +===============================================================+======================+====================================================================+====================================================================================================================+
        | ret_ptype                                                     | ptype_size           | always                                                             | interaction data requested: IE or total (depending on return_total_data) e/g/h (depending on driver)               |
        |                                                               |                      |                                                                    |   with cp/nocp/vmfc treatment (depending on 1st of bsse_type)                                                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | ret_energy                                                    | 1                    | always                                                             | interaction energy: IE or total (depending on return_total_data) w/ cp/nocp/vmfc treat. (dep. on 1st of bsse_type) |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | ret_gradient                                                  | (nat, 3)             | when driver is g/h                                                 | interaction gradient: IE or total (depending on return_total_data) w/ cp/nocp/vmfc treat. (dep. on 1st of bsse_type|
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | ret_hessian                                                   | (nat * 3, nat * 3)   | when driver is h                                                   | interaction Hessian: IE or total (depending on return_total_data) w/ cp/nocp/vmfc treat. (dep. on 1st of bsse_type)|
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | nbody                                                         | >=1                  | always                                                             | energy n-body QCVariables to be set                                                                                |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY               |  |em| 1              | when cp in bsse_type                                               | MBE sum of subsystems of 1-body. summed are total energies with cp treatment                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY               |  |em| 1              | when cp in bsse_type & max_nbody>=2                                | MBE sum of subsystems of 2-body or fewer (cumulative); summed are total energies with cp treatment                 |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY            |  |em| 1              | when cp in bsse_type                                               | MBE sum of subsystems of {max_nbody}-body or fewer (cumulative); summed are total energies w/ cp treatment         |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED TOTAL ENERGY                              |  |em| 1              | when cp in bsse_type & rtd=T                                       | best available total energy with cp treatment: CP-CORRECTED TOTAL ENERGY THROUGH {max_nbody}-BODY                  |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY         |  |em| 1              | when cp in bsse_type & max_nbody>=2                                | 2-body total data less 1-body total data for cumulative IE; inputs are total energies with cp treatment            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY      |  |em| 1              | when cp in bsse_type                                               | {max_nbody}-body total data less 1-body total data for cumulative IE; inputs are total energies with cp treatment  |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED INTERACTION ENERGY                        |  |em| 1              | when cp in bsse_type                                               | best available interaction energy with cp treatment: CP-CORRECTED INTERACTION ENERGY THROUGH {max_nbody}-BODY      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY             |  |em| 1              | when cp in bsse_type & max_nbody>=2                                | 2-body total data less (2-1)-body total data for partial IE; inputs are total energies w/ cp treatment             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| CP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY          |  |em| 1              | when cp in bsse_type                                               | {max_nbody}-body total data less ({max_nbody}-1)-body data for partial IE; inputs are total energies w/ cp treat.  |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY             |  |em| 1              | when nocp in bsse_type                                             | MBE sum of subsystems of 1-body. summed are total energies without cp treatment                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY             |  |em| 1              | when nocp in bsse_type & max_nbody>=2                              | MBE sum of subsystems of 2-body or fewer (cumulative); summed are total energies without cp treatment              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY          |  |em| 1              | when nocp in bsse_type                                             | MBE sum of subsystems of {max_nbody}-body or fewer (cumulative); summed are total energies w/o cp treatment        |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED TOTAL ENERGY                            |  |em| 1              | when nocp in bsse_type                                             | best available total energy without cp treatment: NOCP-CORRECTED TOTAL ENERGY THROUGH {max_nbody}-BODY             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY       |  |em| 1              | when nocp in bsse_type & max_nbody>=2                              | 2-body total data less 1-body total data for cumulative IE; inputs are total energies w/o cp treatment             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY    |  |em| 1              | when nocp in bsse_type                                             | {max_nbody}-body total data less 1-body total data for cumulative IE; inputs are total energies w/o cp treatment   |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED INTERACTION ENERGY                      |  |em| 1              | when nocp in bsse_type                                             | best available interaction energy without cp treatment: NOCP-CORRECTED INTERACTION ENERGY THROUGH {max_nbody}-BODY |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY           |  |em| 1              | when nocp in bsse_type & max_nbody>=2                              | 2-body total data less (2-1)-body total data for partial IE; inputs are total energies w/o cp treatment            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| NOCP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY        |  |em| 1              | when nocp in bsse_type                                             | {max_nbody}-body total data less ({max_nbody}-1)-body data for partial IE; inputs are total energies w/o cp treat. |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY             |  |em| 1              | when vmfc in bsse_type                                             | MBE sum of subsystems of 1-body. summed are total energies with vmfc treatment                                     |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY             |  |em| 1              | when vmfc in bsse_type & max_nbody>=2                              | MBE sum of subsystems of 2-body or fewer (cumulative); summed are total energies with vmfc treatment               |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY          |  |em| 1              | when vmfc in bsse_type                                             | MBE sum of subsystems of {max_nbody}-body or fewer (cumulative); summed are total energies w/ vmfc treatment       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED TOTAL ENERGY                            |  |em| 1              | when vmfc in bsse_type                                             | best available total energy with vmfc treatment: VMFC-CORRECTED TOTAL ENERGY THROUGH {max_nbody}-BODY              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY       |  |em| 1              | when vmfc in bsse_type & max_nbody>=2                              | 2-body total data less 1-body total data for cumulative IE; inputs are total energies w/ vmfc treatment            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY    |  |em| 1              | when vmfc in bsse_type                                             | {max_nbody}-body total data less 1-body total data for cumulative IE; inputs are total energies w/ vmfc treatment  |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED INTERACTION ENERGY                      |  |em| 1              | when vmfc in bsse_type                                             | best available interaction energy with vmfc treatment: VMFC-CORRECTED INTERACTION ENERGY THROUGH {max_nbody}-BODY  |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY           |  |em| 1              | when vmfc in bsse_type & max_nbody>=2                              | 2-body total data less (2-1)-body total data for partial IE; inputs are total energies w/ vmfc treatment           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| VMFC-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY        |  |em| 1              | when vmfc in bsse_type                                             | {max_nbody}-body total data less ({max_nbody}-1)-body data for partial IE; inputs are total energies w/ vmfc treat.|
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | energy_body_dict                                              | max_nbody            | always                                                             | total energies at each n-body level                                                                                |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 1                                                       |  |em| 1              | always; zeroed if cp & rtd=F                                       | cumulative through 1-body total energies w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 2                                                       |  |em| 1              | max_nbody>=2                                                       | cumulative through 2-body total energies w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| {max_nbody}                                             |  |em| 1              | always                                                             | cumulative through {max_nbody}-body total energies w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | gradient_body_dict                                            | max_nbody            | when driver is g/h                                                 |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 1                                                       |  |em| (nat, 3)       | when driver is g/h                                                 | cumulative through 1-body total gradients with cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                   |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 2                                                       |  |em| (nat, 3)       | when driver is g/h & max_nbody>=2                                  | cumulative through 2-body total gradients with cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                   |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| {max_nbody}                                             |  |em| (nat, 3)       | when driver is g/h                                                 | cumulative through {max_nbody}-body total gradients w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | hessian_body_dict                                             | max_nbody            | when driver is h                                                   |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 1                                                       |  |em| (nat*3, nat*3) | when driver is h                                                   | cumulative through 1-body total Hessians w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| 2                                                       |  |em| (nat*3, nat*3) | when driver is h & max_nbody>=2                                    | cumulative through 2-body total Hessians w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |  |em| {max_nbody}                                             |  |em| (nat*3, nat*3) | when driver is h                                                   | cumulative through {max_nbody}-body total Hessians w/ cp/nocp/vmfc treatment (dep. on 1st of bsse_type)            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | cp_energy_body_dict                                           | max_nbody            | always; zeroed if cp not in bsse_type                              | total energies at each n-body level with cp treatment                                                              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1cp                                                    |  |em| 1              | always; zeroed if cp not in bsse_type or rtd=F                     | cumulative through 1-body total energies with cp treatment                                                         |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2cp                                                    |  |em| 1              | when max_nbody>=2; zeroed if cp not in bsse_type                   | cumulative through 2-body total energies with cp treatment                                                         |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}cp                                          |  |em| 1              | always; zeroed if cp not in bsse_type                              | cumulative through {max_nbody}-body total energies with cp treatment                                               |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | cp_gradient_body_dict                                         | max_nbody            | when driver is g/h; zeroed if cp not in bsse_type                  | total gradients at each n-body level with cp treatment                                                             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1cp                                                    |  |em| (nat, 3)       | when driver is g/h; zeroed if cp not in bsse_type or rtd=F         | cumulative through 1-body total gradients with cp treatment                                                        |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2cp                                                    |  |em| (nat, 3)       | when driver is g/h & max_nbody>=2; zeroed if cp not in bsse_type   | cumulative through 2-body total gradients with cp treatment                                                        |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}cp                                          |  |em| (nat, 3)       | when driver is g/h; zeroed if cp not in bsse_type                  | cumulative through {max_nbody}-body total gradients with cp treatment                                              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | cp_hessian_body_dict                                          | max_nbody            | when driver is h; zeroed if cp not in bsse_type                    | total Hessians at each n-body level with cp treatment                                                              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1cp                                                    |  |em| (nat*3, nat*3) | when driver is h; zeroed if cp not in bsse_type or rtd=F           | cumulative through 1-body total Hessians with cp treatment                                                         |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2cp                                                    |  |em| (nat*3, nat*3) | when driver is h & max_nbody>=2; zeroed if cp not in bsse_type     | cumulative through 2-body total Hessians with cp treatment                                                         |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}cp                                          |  |em| (nat*3, nat*3) | when driver is h; zeroed if cp not in bsse_type                    | cumulative through {max_nbody}-body total Hessians with cp treatment                                               |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | nocp_energy_body_dict                                         | max_nbody            | always; zeroed if nocp not in bsse_type                            | total energies at each n-body level with nocp treatment                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1nocp                                                  |  |em| 1              | always; zeroed if nocp not in bsse_type                            | cumulative through 1-body total energies with nocp treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2nocp                                                  |  |em| 1              | when max_nbody>=2; zeroed if nocp not in bsse_type                 | cumulative through 2-body total energies with nocp treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}nocp                                        |  |em| 1              | always; zeroed if nocp not in bsse_type                            | cumulative through {max_nbody}-body total energies with nocp treatment                                             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | nocp_gradient_body_dict                                       | max_nbody            | when driver is g/h; zeroed if nocp not in bsse_type                | total gradients at each n-body level with nocp treatment                                                           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1nocp                                                  |  |em| (nat, 3)       | when driver is g/h; zeroed if nocp not in bsse_type                | cumulative through 1-body total gradients with nocp treatment                                                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2nocp                                                  |  |em| (nat, 3)       | when driver is g/h & max_nbody>=2; zeroed if nocp not in bsse_type | cumulative through 2-body total gradients with nocp treatment                                                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}nocp                                        |  |em| (nat, 3)       | when driver is g/h; zeroed if nocp not in bsse_type                | cumulative through {max_nbody}-body total gradients with nocp treatment                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | nocp_hessian_body_dict                                        | max_nbody            | when driver is h; zeroed if nocp not in bsse_type                  | total Hessians at each n-body level with nocp treatment                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1nocp                                                  |  |em| (nat*3, nat*3) | when driver is h; zeroed if nocp not in bsse_type                  | cumulative through 1-body total Hessians with nocp treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2nocp                                                  |  |em| (nat*3, nat*3) | when driver is h & max_nbody>=2; zeroed if nocp not in bsse_type   | cumulative through 2-body total Hessians with nocp treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}nocp                                        |  |em| (nat*3, nat*3) | when driver is h; zeroed if nocp not in bsse_type                  | cumulative through {max_nbody}-body total Hessians with nocp treatment                                             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | vmfc_energy_body_dict                                         | max_nbody            | always; zeroed if vmfc not in bsse_type                            | total energies at each n-body level with vmfc treatment                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1vmfc                                                  |  |em| 1              | always; zeroed if vmfc not in bsse_type                            | cumulative through 1-body total energies with vmfc treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2vmfc                                                  |  |em| 1              | when max_nbody>=2; zeroed if vmfc not in bsse_type                 | cumulative through 2-body total energies with vmfc treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}vmfc                                        |  |em| 1              | always; zeroed if vmfc not in bsse_type                            | cumulative through {max_nbody}-body total energies with vmfc treatment                                             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | vmfc_gradient_body_dict                                       | max_nbody            | when driver is g/h; zeroed if vmfc not in bsse_type                | total gradients at each n-body level with vmfc treatment                                                           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1vmfc                                                  |  |em| (nat, 3)       | when driver is g/h; zeroed if vmfc not in bsse_type                | cumulative through 1-body total gradients with vmfc treatment                                                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2vmfc                                                  |  |em| (nat, 3)       | when driver is g/h & max_nbody>=2; zeroed if vmfc not in bsse_type | cumulative through 2-body total gradients with vmfc treatment                                                      |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}vmfc                                        |  |em| (nat, 3)       | when driver is g/h; zeroed if vmfc not in bsse_type                | cumulative through {max_nbody}-body total gradients with vmfc treatment                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | vmfc_hessian_body_dict                                        | max_nbody            | when driver is h; zeroed if vmfc not in bsse_type                  | total Hessians at each n-body level with vmfc treatment                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1vmfc                                                  |  |em| (nat*3, nat*3) | when driver is h; zeroed if vmfc not in bsse_type                  | cumulative through 1-body total Hessians with vmfc treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2vmfc                                                  |  |em| (nat*3, nat*3) | when driver is h & max_nbody>=2; zeroed if vmfc not in bsse_type   | cumulative through 2-body total Hessians with vmfc treatment                                                       |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| {max_nbody}vmfc                                        |  |em| (nat*3, nat*3) | when driver is h; zeroed if vmfc not in bsse_type                  | cumulative through {max_nbody}-body total Hessians with vmfc treatment                                             |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | intermediates                                                 | ntasks               | always                                                             | all individual energies with nice labels                                                                           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| N-BODY (1, 2)@(1, 2) TOTAL ENERGY                      |  |em| 1              | always                                                             | total energy for 1st modelchem, 1st & 2nd fragments in basis of 1st & 2nd fragments                                |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| N-BODY (3)@(2, 3) TOTAL ENERGY                         |  |em| 1              | always                                                             | total energy for 2nd modelchem, 3rd fragment in basis of 2nd and 3rd fragments                                     |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| ...                                                    |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | intermediates_energy                                          | ntasks               | always                                                             | all individual energies                                                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1_((1, 2), (1, 2))                                     |  |em| 1              | always                                                             | total energy for 1st modelchem, 1st & 2nd fragments in basis of 1st & 2nd fragments                                |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2_((3,), (2, 3))                                       |  |em| 1              | always                                                             | total energy for 2nd modelchem, 3rd fragment in basis of 2nd and 3rd fragments                                     |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| ...                                                    |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | intermediates_gradient                                        | ntasks               | when driver is g/h                                                 | all individual gradients                                                                                           |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1_((1, 2), (1, 2))                                     |  |em| (nat, 3)       | when driver is g/h                                                 | total gradient for 1st modelchem, 1st & 2nd fragments in basis of 1st & 2nd fragments                              |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2_((3,), (2, 3))                                       |  |em| (nat, 3)       | when driver is g/h                                                 | total gradient for 2nd modelchem, 3rd fragment in basis of 2nd and 3rd fragments                                   |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| ...                                                    |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |                                                               |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        | intermediates_hessian                                         | ntasks               | when driver is h                                                   | all individual Hessians                                                                                            |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 1_((1, 2), (1, 2))                                     |  |em| (nat*3, nat*3) | when driver is h                                                   | total Hessian for 1st modelchem, 1st & 2nd fragments in basis of 1st & 2nd fragments                               |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| 2_((3,), (2, 3))                                       |  |em| (nat*3, nat*3) | when driver is h                                                   | total Hessian for 2nd modelchem, 3rd fragment in basis of 2nd and 3rd fragments                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
        |   |em| ...                                                    |                      |                                                                    |                                                                                                                    |
        +---------------------------------------------------------------+----------------------+--------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+

    """
    pass

# TODO questions to check:
# * can work with supersystem and embedding_charges?
# * can levels work with same method, different basis?


class ManyBodyComputer(ManyBodyComputerQCNG):

    # v2, probably: @field_validator("molecule", mode="before")
    @validator("molecule", pre=True)
    @classmethod
    def set_molecule(cls, v: Any) -> Molecule:
        if isinstance(v, core.Molecule):
            v = v.to_schema(dtype=2)

        return v

    @classmethod
    def from_psi4_task_planner(cls, *, molecule, method, basis, driver, keywords, **mbkwargs):

        atomic_spec = AtomicSpecification(
            model={"method": method, "basis": basis},
            program="psi4",
            driver=driver,
            keywords=keywords,
        )

        if mbkwargs.get("levels"):
            specification = {v: atomic_spec.copy(deep=True) for v in mbkwargs["levels"].values()}
        else:
            specification = atomic_spec  # will turn into {"(auto)": atomic_spec}

        input_model = ManyBodyInput(
            specification={
                "specification": specification,
                "keywords": mbkwargs,
                "driver": driver,
            },
            molecule=molecule.to_schema(dtype=2),
        )

        computer_model = cls(
            molecule=input_model.molecule,
            driver=input_model.specification.driver,
            # v2: **input_model.specification.keywords.model_dump(exclude_unset=True),
            **input_model.specification.keywords.dict(exclude_unset=True),
            input_data=input_model,  # storage, to reconstitute ManyBodyResult
        )
        nb_per_mc = computer_model.nbodies_per_mc_level

        # print("\n<<<  (ZZZ 1) Psi4 harness ManyBodyComputer.from_psi4_task_planner  >>>")
        # v2: pprint.pprint(computer_model.model_dump(), width=200)
        # pprint.pprint(computer_model.dict(), width=200)
        # print(f"nbodies_per_mc_level={nb_per_mc}")

        comp_levels = {}
        for mc_level_idx, mtd in enumerate(computer_model.levels.values()):
            for lvl1 in nb_per_mc[mc_level_idx]:
                key = "supersystem" if lvl1 == "supersystem" else int(lvl1)
                comp_levels[key] = mtd

        specifications = {}
        for mtd, spec in computer_model.input_data.specification.specification.items():
            spec = spec.dict()
            specifications[mtd] = {}
            specifications[mtd]["program"] = spec.pop("program")
            specifications[mtd]["specification"] = spec
            specifications[mtd]["specification"].pop("schema_name", None)
            specifications[mtd]["specification"].pop("protocols", None)

        computer_model.qcmb_core = ManyBodyCore(
            computer_model.molecule,
            computer_model.bsse_type,
            comp_levels,
            return_total_data=computer_model.return_total_data,
            supersystem_ie_only=computer_model.supersystem_ie_only,
            embedding_charges=computer_model.embedding_charges,
        )

        info, dcount = computer_model.qcmb_core.format_calc_plan()
        logger.info(info)
        core.print_out("\n" + p4util.banner(f" ManyBody Setup: N-Body Levels {computer_model.max_nbody}", strNotOutfile=True) + "\n")
        core.print_out(info)

        #    if not result.success:
        #        print(result.error.error_message)
        #        raise RuntimeError("Calculation did not succeed! Error:\n" + result.error.error_message)

        return computer_model

    def build_tasks(
        self,
        mb_computer: AtomicComputer,  # SubTaskComputers,
        mc_level_idx: int,
        **kwargs: Dict[str, Any],
    ) -> int:
        """Adds to the task_list as many new unique tasks as necessary to treat a single model chemistry level at one
        or several n-body levels. New tasks are of type *mb_computer* with model chemistry level specified in *kwargs*
        and n-body levels accessed through *mc_level_idx*.

        Parameters
        ----------
        mb_computer
            Class of task computers to instantiate and add to self.task_list. Usually :class:`~psi4.driver.AtomicComputer` but may be other when wrappers are layered.
        mc_level_idx
            Position in field self.nbodies_per_mc_level used to obtain ``nbodies``, the list of n-body
            levels (e.g., `[1]` or `[1, 2]` or `["supersystem"]`) to which the modelchem specified in **kwargs** applies.
            That is, `nbodies = self.nbodies_per_mc_level[mc_level_idx]`.
            Note the natural 1-indexing of ``nbodies`` _contents_, so `[1]` covers one-body contributions.
            The corresponding user label is the 1-indexed counterpart, `mc_level_lbl = mc_level_idx + 1`
            Formerly nlevel as in `nbody = self.nbody_list[nbody_level=nlevel]`. (edit)
        kwargs
            Other arguments for initializing **mb_computer**. In particular, specifies model chemistry. (edit)

        Returns
        -------
        count : int
            Number of new tasks planned by this call including any supersystem.

        """
        # TODO method not coming from levels right

        # Get the n-body orders for this level. e.g., [1] or [2, 3] or ["supersystem"]
        nbodies = self.nbodies_per_mc_level[mc_level_idx]
        # print(f"{self.nbodies_per_mc_level=} {nbodies=} {mc_level_idx=}")
        # print(f"{self.levels=}")
        # print(f"{self.qcmb_core.nbodies_per_mc_level=}")

        for spec_key, nb_lst in self.qcmb_core.nbodies_per_mc_level.items():
            if nb_lst == nbodies:
                filter_key = spec_key
        #filter_key = rev_nbpermclvl[self.nbodies_per_mc_level[mc_level_idx]]
        # print(f"{filter_key=}")

        for kwg in ['dft_functional']:
            if kwg in kwargs:
                kwargs['keywords']['function_kwargs'][kwg] = kwargs.pop(kwg)

        count = 0
        template = copy.deepcopy(kwargs)
        logger.info(f"TEMPLATE {template=}")

        for chem, label, imol in self.qcmb_core.iterate_molecules():
            mckey, frag, bas = delabeler(label)
            # print(f"{chem=} {label=}, {imol=}      {mckey=} {frag=} {bas=}    {imol.extras=}")
            if mckey != filter_key:
                continue
            # TODO probably haven't filtered by body
            if label in self.task_list:
                continue

            data = copy.deepcopy(template)
            data["molecule"] = core.Molecule.from_schema(imol.dict()) #nonphysical=True
            if self.embedding_charges:
                charges = imol.extras.pop("embedding_charges")
                data['keywords']['function_kwargs'].update({'external_potentials': charges})

            self.task_list[label] = mb_computer(**data)
            count += 1

        return count

    def compute(self, client: Optional["qcportal.client.PortalClient"] = None):
        """Run quantum chemistry.

        Parameters
        ----------
        client
            QCFractal client if using QCArchive for distributed compute.

        """
        info = "\n" + p4util.banner(f" ManyBody Computations ", strNotOutfile=True) + "\n"
        logger.info(info)

        with p4util.hold_options_state():
            for t in self.task_list.values():
                t.compute(client=client)

    def get_results(
        self,
        client: Optional["qcportal.client.PortalClient"] = None) -> ManyBodyResult:
        """Process the QC results from all n-body component molecular systems
        and model chemistry levels into final quantities and form QCSchema model.

        Parameters
        ----------
        client
            QCFractal client if using QCArchive for distributed compute.

        Returns
        -------
        mbres
            All MBE results collected into a QCSchema model.

        """
        component_results = {k: v.get_results(client=client) for k, v in self.task_list.items()}

#            if not result.success:
#                print(result.error.error_message)
#                raise RuntimeError("Calculation did not succeed! Error:\n" + result.error.error_message)

        component_properties = {}
        props = {"energy", "gradient", "hessian"}

        for label, result in component_results.items():
            component_properties[label] = {}

            for egh in props:
                v = getattr(result.properties, f"return_{egh}", None)
                if v is not None:
                    component_properties[label][egh] = v

        findif_number = {label: res.extras["qcvars"].get("FINDIF NUMBER") for label, res in component_results.items()}

        # print("\n<<<  (RRR 2) Psi4 ManyBodyComputer.get_results component_properties  >>>")
        # pprint.pprint(component_properties, width=200)

        analyze_back = self.qcmb_core.analyze(component_properties)
        analyze_back["nbody_number"] = len(component_properties)

        nbody_model = super().get_results(external_results=analyze_back, component_results=component_results, client=client)

        return nbody_model

    def get_psi_results(
        self,
        client: Optional["qcportal.client.PortalClient"] = None,
        *,
        return_wfn: bool = False) -> EnergyGradientHessianWfnReturn:
        """Called by driver to assemble results into ManyBody-flavored QCSchema,
        then reshape and return them in the customary Psi4 driver interface: ``(e/g/h, wfn)``.

        Parameters
        ----------
        return_wfn
            Whether to additionally return the dummy :py:class:`~psi4.core.Wavefunction`
            calculation result as the second element of a tuple. Contents are:

            - supersystem molecule
            - dummy basis, def2-svp
            - e/g/h member data
            - QCVariables

        Returns
        -------
        ret
            Energy, gradient, or Hessian according to self.driver.
        wfn
            Wavefunction described above when *return_wfn* specified.

        """
        nbody_model = self.get_results(client=client)

        # print("\n<<<  (RRR 4) Psi4 ManyBodyComputer.get_psi_results nbody_model  >>>")
        # pprint.pprint(nbody_model.dict(), width=200)

        info = "\n" + p4util.banner(f" ManyBody Results ", strNotOutfile=True) + "\n"
        core.print_out(info)
        core.print_out(nbody_model.stdout)
        logger.info('\nQCManyBody:\n' + nbody_model.stdout)

        # v2: nbody_model.model_dump()
        logger.debug('\nManyBodyResult QCSchema:\n' + pp.pformat(nbody_model.dict()))

        qmol = core.Molecule.from_schema(self.molecule.dict())  # nonphy
        wfn = core.Wavefunction.build(qmol, "def2-svp", quiet=True)

        qcvars = translate_qcvariables(nbody_model.properties.dict())
        for qcv, val in qcvars.items():
            for obj in [core, wfn]:
                obj.set_variable(qcv, val)

        mc_labels = modelchem_labels(self.qcmb_core.nbodies_per_mc_level)
        full_to_ordinal_mc_lbl = {v[0]: v[1] for v in mc_labels.values()}

        for qcmb_label, dval in nbody_model.component_properties.items():
            mckey, frag, bas = delabeler(qcmb_label)
            if len(full_to_ordinal_mc_lbl) == 1:
                mc_lbl = None
            else:
                mc_lbl = full_to_ordinal_mc_lbl[mckey][1:]  # avoid §§
            viz_label = labeler(mc_lbl, frag, bas, opaque=False)

            available_properties = dval.dict().keys()
            for obj in [core, wfn]:
                if "return_energy" in available_properties:
                    obj.set_variable(f"N-BODY {viz_label} TOTAL ENERGY", dval.return_energy)
                if "return_gradient" in available_properties:
                    obj.set_variable(f"N-BODY {viz_label} TOTAL GRADIENT", dval.return_gradient)
                if "return_hessian" in available_properties:
                    obj.set_variable(f"N-BODY {viz_label} TOTAL HESSIAN", dval.return_hessian)

        if self.driver == "energy":
            ret = nbody_model.return_result
        elif self.driver == "gradient":
            ret = core.Matrix.from_array(nbody_model.return_result)
            wfn.set_gradient(ret)
        elif self.driver == "hessian":
            ret = core.Matrix.from_array(nbody_model.return_result)
            grad = core.Matrix.from_array(nbody_model.properties.return_gradient)
            wfn.set_hessian(ret)
            wfn.set_gradient(grad)

        if return_wfn:
            return (ret, wfn)
        else:
            return ret
