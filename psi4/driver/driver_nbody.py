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

        Computer._prepare_results()
        ---------------------------
        * if multiple modelchems (multilevel):

            multilevel.prepare_results()
            ----------------------------
            * from the pool of calcs, partition them by modelchem treatment and call _prepare_results on each subpool
            * sums modelchem levels and returns small dict back to get_results()

        * call get_results() for each job in task list
        * assemble all the computed energies, all the computed gradients, and all the computed hessians
        * for each available derivative, call:

            assemble_nbody_components()
            ---------------------------
            * re-call build_nbody_compute_list to get the cp/nocp/vmfc lists again

                build_nbody_compute_list()
                --------------------------

            * slice up the supersystem mol into fragment atom ranges to help g/h arrays build piecemeal
            * prepare empty {bsse_type}_by_level and {bsse_type}_body_dict structs. the former have different contents for vmfc
            * for cp and nocp, resort the build_nbody_compute_list returns into per-body lists suitable for summing
            * note that nb loops often run over more than active nbodies_per_mc_level item due to 1-body for subtraction and multilevel complications
            * for each possibly active n-body level and each active bsse_type, call _sum_cluster_ptype_data to build by_level structs

                _sum_cluster_ptype_data()
                -------------------------
                * sum up ene, grad, or Hess in per-fragment pieces based on list of (frag, bas) subjobs active for that bsse treatment

            * compute special case of monomers in monomer basis
            * for each of cp/nocp/vmfc, apply appropriate formula to build each n-body level of cumulative total energy into body_dict
            * for driver=energy, set several qcvars and call:

                _print_nbody_energy()
                ---------------------
                * prints and logs formatted energy output. called separately for cp, nocp, vmfc

            * collect qcvars and summed levels into a return dictionary with some extra aliases for target bsse_type and target driver

        * merge all the assemble_nbody_components return dictionaries
        * in struct["intermediates"], store dict of `"N-BODY (?)@(?) TOTAL ENERGY" = return_energy` for all in task_list or results kwarg
        * in struct["intermediates_{ptype}"], store dict of `task_list key = return_{ptype}` for all in task_list or results kwarg. ptype=e/g/h
          always for ptype=energy, as available for higher derivatives when driver=g/h

    * form nbody qcvars and properties, inc'l number, current e/g/h as available
    * pull results (incl dicts!) into qcvars
    * form model, including copy of class with mols converted to qcsk at atomicresult.extras["component_results"]

* collect ManyBody-flavored AtomicResult from self.get_results()
* build wfn from nbody mol and basis (always def2-svp)
* push qcvars to P::e and wfn. push various internal dicts to qcvars, too
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
from ast import literal_eval
# v2: from typing import TYPE_CHECKING, Any, ClassVar, Dict, List, Tuple, Union, Optional
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Tuple, Union, Optional

# printing and logging formatting niceties
import pprint
from functools import partial
import numpy as np
pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
nppp = partial(np.array_str, max_line_width=120, precision=8, suppress_small=True)
nppp10 = partial(np.array_str, max_line_width=120, precision=10, suppress_small=True)

from psi4 import core

from . import p4util
from .constants import constants, pp
from .driver_cbs import CompositeComputer
from .driver_findif import FiniteDifferenceComputer
from .p4util.exceptions import *
from .task_base import AtomicComputer, BaseComputer, EnergyGradientHessianWfnReturn

try:
    from pydantic.v1 import validator
except ImportError:
    from pydantic import validator

from qcelemental.models import FailedOperation, Molecule, DriverEnum, ProtoModel, AtomicResult, AtomicInput
import qcengine as qcng

import qcmanybody as qcmb
from qcmanybody.models import AtomicSpecification, BsseEnum, ManyBodyKeywords, ManyBodyInput, ManyBodyResult, ManyBodyResultProperties
from qcmanybody import ManyBodyCalculator
from qcmanybody.computer import ManyBodyComputer as ManyBodyComputerQCNG
from qcmanybody.utils import delabeler, provenance_stamp as qcmb_provenance_stamp

if TYPE_CHECKING:
    import qcportal

logger = logging.getLogger(__name__)

FragBasIndex = Tuple[Tuple[int], Tuple[int]]

SubTaskComputers = Union[AtomicComputer, CompositeComputer, FiniteDifferenceComputer]

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

    :type levels: dict
    :param levels: ``{1: 'ccsd(t)', 2: 'mp2', 'supersystem': 'scf'}`` || ``{1: 2, 2: 'ccsd(t)', 3: 'mp2'}`` || etc

        Dictionary of different levels of theory for different levels of expansion
        Note that method_string is not used in this case. ``supersystem`` computes
        all higher order n-body effects up to the number of fragments.

    :type embedding_charges: dict
    :param embedding_charges: ``{1: [-0.834, 0.417, 0.417], ..}``

        Dictionary of atom-centered point charges. keys: 1-based index of fragment, values: list of charges for each fragment.
        Add atom-centered point charges for fragments whose basis sets are not included in the computation.

    """
    pass


hide_stuff = '''

    def plan(self):
        # uncalled function
        return [t.plan() for t in self.task_list.values()]

    def prepare_results(
        self,
        results: Optional[Dict[str, SubTaskComputers]] = None,
        client: Optional["qcportal.client.PortalClient"] = None,
    ) -> Dict[str, Any]:
        """Process the results from all n-body component molecular systems and model chemistry levels into final quantities.

        Parameters
        ----------
        results
            A set of tasks to process instead of self.task_list. Used in multilevel processing to pass a subset of
            self.task_list filtered to only one modelchem level.
        client
            QCFractal client if using QCArchive for distributed compute.

        Returns
        -------
        nbody_results
            When the ManyBodyComputer specifies a single model chemistry level (see self.nbodies_per_mc_level), the
            return is a dictionary, nbody_results, described in the table below. Many of the items are actually filled
            by successive calls to assemble_nbody_components(). When multiple model chemistry levels are specified, this
            function diverts its return to driver_nbody_multilevel.prepare_results() wherein each mc level calls this
            function again and collects separate nbody_results dictionaries and processes them into a final return that
            is a small subset of the table below.


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
        if results is None:
            results = {}

        # formerly nlevels
        mc_level_labels = {i.split("_")[0] for i in self.task_list}
        if len(mc_level_labels) > 1 and not results:
            return driver_nbody_multilevel.prepare_results(self, client)

        results_list = {k: v.get_results(client=client) for k, v in (results.items() or self.task_list.items())}
        trove = {  # AtomicResult.properties return None if missing
            "energy": {k: v.properties.return_energy for k, v in results_list.items()},
            "gradient": {k: v.properties.return_gradient for k, v in results_list.items()},
            "hessian": {k: v.properties.return_hessian for k, v in results_list.items()},
        }

        # TODO: make assemble_nbody_components and driver_nbody_multilevel.prepare_results into class functions.
        #   note that the former uses metadata as read-only (except for one solveable case) while the latter overwrites self (!).
        metadata = {
            "quiet": self.quiet,
            "nbodies_per_mc_level": self.nbodies_per_mc_level,
            "bsse_type": self.bsse_type,
            "nfragments": self.nfragments,
            "return_total_data": self.return_total_data,
            "molecule": self.molecule,
            "embedding_charges": bool(self.embedding_charges),
            "max_nbody": self.max_nbody,
        }
        if self.driver.name == "energy":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())

        elif self.driver.name == "gradient":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())
            nbody_results.update(assemble_nbody_components("gradient", trove["gradient"], metadata.copy()))

        elif self.driver.name == "hessian":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())
            nbody_results.update(assemble_nbody_components("gradient", trove["gradient"], metadata.copy()))
            nbody_results.update(assemble_nbody_components("hessian", trove["hessian"], metadata.copy()))

        def delabeler(item: str, return_obj: bool = False) -> Union[Tuple[str, str, str], Tuple[int, Tuple[int], Tuple[int]]]:
            """Transform labels like string "1_((2,), (1, 2))" into string tuple ("1", "2", "1, 2") or object tuple (1, (2,), (1, 2))."""

            mc, _, fragbas = item.partition("_")
            frag, bas = literal_eval(fragbas)

            if return_obj:
                return int(mc), frag, bas
            else:
                return mc, ", ".join(map(str, frag)), ", ".join(map(str, bas))

        # save some mc_(frag, bas) component results
        # * formerly, intermediates_energy was intermediates2
        # * formerly, intermediates_gradient was intermediates_ptype
        # * formerly, intermediates_hessian was intermediates_ptype

        nbody_results["intermediates"] = {}
        for idx, task in results_list.items():
            mc, frag, bas = delabeler(idx)
            nbody_results["intermediates"][f"N-BODY ({frag})@({bas}) TOTAL ENERGY"] = task.properties.return_energy

        nbody_results["intermediates_energy"] = trove["energy"]

        if not all(x is None for x in trove["gradient"].values()):
            nbody_results["intermediates_gradient"] = trove["gradient"]

        if not all(x is None for x in trove["hessian"].values()):
            nbody_results["intermediates_hessian"] = trove["hessian"]

        debug = False
        if debug:
            for k, v in nbody_results.items():
                if isinstance(v, np.ndarray):
                    print(f"CLS-prepared results >>> {k} {v.size}")
                elif isinstance(v, dict):
                    print(f"CLS-prepared results >>> {k} {len(v)}")
                    for k2, v2 in v.items():
                        if isinstance(v2, np.ndarray):
                            print(f"CLS-prepared results      >>> {k2} {v2.size}")
                        else:
                            print(f"CLS-prepared results      >>> {k2} {v2}")
                else:
                    print(f"CLS-prepared results >>> {k} {v}")

        return nbody_results

    def get_results(self, client: Optional["qcportal.client.PortalClient"] = None) -> AtomicResult:
        """Return results as ManyBody-flavored QCSchema."""

        info = "\n" + p4util.banner(f" ManyBody Results ", strNotOutfile=True) + "\n"
        core.print_out(info)
        logger.info(info)

        results = self.prepare_results(client=client)
        ret_energy = results.pop("ret_energy")
        ret_ptype = results.pop("ret_ptype")
        ret_gradient = results.pop("ret_gradient", None)

        # load QCVariables
        qcvars = {
            'NUCLEAR REPULSION ENERGY': self.molecule.nuclear_repulsion_energy(),
            'NBODY NUMBER': len(self.task_list),
        }

        properties = {
            "calcinfo_natom": self.molecule.natom(),
            "nuclear_repulsion_energy": self.molecule.nuclear_repulsion_energy(),
            "return_energy": ret_energy,
        }

        for k, val in results.items():
            qcvars[k] = val

        qcvars['CURRENT ENERGY'] = ret_energy
        if self.driver == 'gradient':
            qcvars['CURRENT GRADIENT'] = ret_ptype
            properties["return_gradient"] = ret_ptype
        elif self.driver == 'hessian':
            qcvars['CURRENT GRADIENT'] = ret_gradient
            qcvars['CURRENT HESSIAN'] = ret_ptype
            properties["return_gradient"] = ret_gradient
            properties["return_hessian"] = ret_ptype

        component_results = self.dict()['task_list']
        for k, val in component_results.items():
            val['molecule'] = val['molecule'].to_schema(dtype=2)

        nbody_model = AtomicResult(
            **{
                'driver': self.driver,
                'model': {
                    'method': self.method,
                    'basis': self.basis,
                },
                'molecule': self.molecule.to_schema(dtype=2),
                'properties': properties,
                'provenance': p4util.provenance_stamp(__name__),
                'extras': {
                    'qcvars': qcvars,
                    'component_results': component_results,
                },
                'return_result': ret_ptype,
                'success': True,
            })

        logger.debug('\nNBODY QCSchema:\n' + pp.pformat(nbody_model.dict()))

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
        ret = nbody_model.return_result

        wfn = core.Wavefunction.build(self.molecule, "def2-svp", quiet=True)

        # TODO all besides nbody may be better candidates for extras than qcvars. energy/gradient/hessian_body_dict in particular are too simple for qcvars (e.g., "2")
        dicts = [
            #"energies",  # retired
            #"ptype",     # retired
            "intermediates",
            "intermediates_energy",  #"intermediates2",
            "intermediates_gradient",  #"intermediates_ptype",
            "intermediates_hessian",  #"intermediates_ptype",
            "energy_body_dict",
            "gradient_body_dict",  # ptype_body_dict
            "hessian_body_dict",  # ptype_body_dict
            "nbody",
            "cp_energy_body_dict",
            "nocp_energy_body_dict",
            "vmfc_energy_body_dict",
            "cp_gradient_body_dict",
            "nocp_gradient_body_dict",
            "vmfc_gradient_body_dict",
            "cp_hessian_body_dict",
            "nocp_hessian_body_dict",
            "vmfc_hessian_body_dict",
        ]

        for qcv, val in nbody_model.extras['qcvars'].items():
            if isinstance(val, dict):
                if qcv in dicts:
                    for qcv2, val2 in val.items():
                        for obj in [core, wfn]:
                            try:
                                obj.set_variable(str(qcv2), val2)
                            except ValidationError:
                                obj.set_variable(f"{self.driver.name} {qcv2}", val2)
            else:
                for obj in [core, wfn]:
                    obj.set_variable(qcv, val)

        if self.driver == 'gradient':
            ret = core.Matrix.from_array(ret)
            wfn.set_gradient(ret)
        elif self.driver == 'hessian':
            ret = core.Matrix.from_array(ret)
            grad = core.Matrix.from_array(nbody_model.properties.return_gradient)
            wfn.set_hessian(ret)
            wfn.set_gradient(grad)

        if return_wfn:
            return (ret, wfn)
        else:
            return ret


# TODO questions to check:
# * can work with supersystem and embedding_charges?
# * can levels work with same method, different basis?

#####################################################################################################



'''


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

        #if "levels" in mbkwargs:
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
        # print(f"{comp_levels=}")

        specifications = {}
        for mtd, spec in computer_model.input_data.specification.specification.items():
            spec = spec.dict()
            specifications[mtd] = {}
            specifications[mtd]["program"] = spec.pop("program")
            specifications[mtd]["specification"] = spec
            specifications[mtd]["specification"].pop("schema_name", None)
            specifications[mtd]["specification"].pop("protocols", None)
        # print("POST specifications")
        # pprint.pprint(specifications)

        calculator_cls = ManyBodyCalculator(
            computer_model.molecule,
            computer_model.bsse_type,
            comp_levels,
            computer_model.return_total_data,
            computer_model.supersystem_ie_only,
            computer_model.embedding_charges,
        )
        computer_model.qcmb_calculator = calculator_cls

        # print("\n<<<  (ZZZ 2) QCManyBody module ManyBodyCalculator  >>>")
        # print(dir(calculator_cls))

        info, dcount = calculator_cls.format_calc_plan()
        # print(info)
        # pprint.pprint(dcount)
        logger.info(info)
        core.print_out("\n" + p4util.banner(f" ManyBody Setup: N-Body Levels {computer_model.max_nbody}", strNotOutfile=True) + "\n")
        core.print_out(info)

        component_results = {}

        for chem, label, imol in calculator_cls.iterate_molecules():
            # print(f"{chem=} {label=}, {imol=} {imol.extras=}")
            pass
        #    inp = AtomicInput(molecule=imol, **specifications[chem]["specification"])

        #    _, real, bas = delabeler(label)
        #    result = qcng.compute(inp, specifications[chem]["program"])

        #    if not result.success:
        #        print(result.error.error_message)
        #        raise RuntimeError("Calculation did not succeed! Error:\n" + result.error.error_message)

        #    # pull out stuff
        #    props = {"energy", "gradient", "hessian"}

        #    component_results[label] = {}

        #    for p in props:
        #        if hasattr(result.properties, f"return_{p}"):
        #            v = getattr(result.properties, f"return_{p}")
        #            # print(f"  {label} {p}: {v}")
        #            if v is not None:
        #                component_results[label][p] = v

        #print("\n<<<  (ZZZ 2) Psi4 harness ManyBodyComputer.from_psi4_task_planner component_results  >>>")
        #pprint.pprint(component_results, width=200)

        #print("start to analyze")
        #analyze_back = calculator_cls.analyze(component_results)
        #analyze_back["nbody_number"] = len(component_results)
        #print("\n<<<  (ZZZ 3) Psi4 harness ManyBodyComputer.from_psi4_task_planner analyze_back  >>>")
        #pprint.pprint(analyze_back, width=200)

        #return computer_model.get_results(external_results=analyze_back)
        return computer_model

    def build_tasks(
        self,
        mb_computer: AtomicComputer, #MBETaskComputers,
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
            Number of new tasks planned by this call.
            Formerly, didn't include supersystem in count.

        """
# qcng:        from psi4.driver.driver_nbody import build_nbody_compute_list

        # TODO method not coming from levels right

        # Get the n-body orders for this level. e.g., [1] or [2, 3] or ["supersystem"]
        nbodies = self.nbodies_per_mc_level[mc_level_idx]
        # print(f"{self.nbodies_per_mc_level=} {nbodies=} {mc_level_idx=}")
        # print(f"{self.levels=}")
        # print(f"{self.qcmb_calculator.nbodies_per_mc_level=}")

#self.nbodies_per_mc_level=[[1], [2]] nbodies=[2] mc_level_idx=1
#self.levels={1: 'mp2/sto-3g', 2: 'scf/sto-3g'}
#self.qcmb_calculator.nbodies_per_mc_level={'scf/sto-3g': [2], 'mp2/sto-3g': [1]}
        for spec_key, nb_lst in self.qcmb_calculator.nbodies_per_mc_level.items():
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

        # FROM CPTR # Get compute list
        # FROM CPTR if nbodies == ["supersystem"]:
        # FROM CPTR     # Add supersystem computation if requested -- always nocp
        # FROM CPTR     data = template
        # FROM CPTR     data["molecule"] = self.molecule
        # FROM CPTR     key = f"supersystem_{self.nfragments}"
        # FROM CPTR     self.task_list[key] = mb_computer(**data)
        # FROM CPTR     count += 1

        # FROM CPTR     compute_dict = build_nbody_compute_list(
        # FROM CPTR         ["nocp"], list(range(1, self.max_nbody + 1)),
        # FROM CPTR         self.nfragments, self.return_total_data, self.supersystem_ie_only)
        # FROM CPTR else:
        # FROM CPTR     compute_dict = build_nbody_compute_list(
        # FROM CPTR         self.bsse_type, nbodies,
        # FROM CPTR         self.nfragments, self.return_total_data, self.supersystem_ie_only 

        def lab_labeler_alt(item) -> str:
            mckey, frag, bas = qcmb.delabeler(label)
            # print("parts", mckey, frag, bas)
            # note 0-index to 1-index shift for label
            return f"{mc_level_idx + 1}_{item}"

        #print("HHHH", compute_dict)
        for chem, label, imol in self.qcmb_calculator.iterate_molecules():
            mckey, frag, bas = qcmb.delabeler(label)
            # print(f"{chem=} {label=}, {imol=}      {mckey=} {frag=} {bas=}    {imol.extras=}")
            lbl2 = lab_labeler_alt(label)
            # print(f"{lbl2=}")
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

            # print(f"ADDING TO TASKS {label=} {lbl2=} {data=} {template=}")
            self.task_list[label] = mb_computer(**data)
            count += 1

    def compute(self, client: Optional["qcportal.client.PortalClient"] = None):
        """Run quantum chemistry."""

        info = "\n" + p4util.banner(f" ManyBody Computations ", strNotOutfile=True) + "\n"
        #core.print_out(info)
        logger.info(info)

        with p4util.hold_options_state():
            for t in self.task_list.values():
                t.compute(client=client)

    def get_results(
        self,
        client: Optional["qcportal.client.PortalClient"] = None) -> ManyBodyResult:
        """
        return back the model.
        """
        results_list = {k: v.get_results(client=client) for k, v in self.task_list.items()}
        # print("RRRRR RESULTS_LIST")
        # pprint.pprint(results_list, width=200)

#        trove = {  # AtomicResult.properties return None if missing
#            "energy": {k: v.properties.return_energy for k, v in results_list.items()},
#            "gradient": {k: v.properties.return_gradient for k, v in results_list.items()},
#            "hessian": {k: v.properties.return_hessian for k, v in results_list.items()},
#        }
        #findif_number = [for k, v in results_list.items()

#        for chem, label, imol in calculator_cls.iterate_molecules():
#            inp = AtomicInput(molecule=imol, **specifications[chem]["specification"])
#
#            _, real, bas = delabeler(label)
#            result = qcng.compute(inp, specifications[chem]["program"])
#
#            if not result.success:
#                print(result.error.error_message)
#                raise RuntimeError("Calculation did not succeed! Error:\n" + result.error.error_message)
#
        component_properties = {}
        # pull out stuff
        props = {"energy", "gradient", "hessian"}

        for label, result in results_list.items():
            component_properties[label] = {}

            for egh in props:
                if hasattr(result.properties, f"return_{egh}"):
                    v = getattr(result.properties, f"return_{egh}")
                    # print(f"  {label} {egh}: {v}")
                    if v is not None:
                        component_properties[label][egh] = v

        #findif_number = set(res.extras["qcvars"].get("FINDIF NUMBER") for label, res in results_list.items())
        findif_number = {label: res.extras["qcvars"].get("FINDIF NUMBER") for label, res in results_list.items()}
        # print(f"{findif_number=}")

        # print("\n<<<  (RRR 2) Psi4 ManyBodyComputer.get_psi_results component_properties  >>>")
        # pprint.pprint(component_properties, width=200)

        # print("start to analyze")
        analyze_back = self.qcmb_calculator.analyze(component_properties)
        analyze_back["nbody_number"] = len(component_properties)
        # print("\n<<<  (RRR 3) Psi4 ManyBodyComputer.get_psi_results analyze_back  >>>")
        # pprint.pprint(analyze_back, width=200)

        nbody_model = super().get_results(external_results=analyze_back, component_results=results_list, client=client)

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
#SEP        results_list = {k: v.get_results(client=client) for k, v in self.task_list.items()}
#SEP        print("RRRRR RESULTS_LIST")
#SEP        pprint.pprint(results_list, width=200)
#SEP
#SEP#        trove = {  # AtomicResult.properties return None if missing
#SEP#            "energy": {k: v.properties.return_energy for k, v in results_list.items()},
#SEP#            "gradient": {k: v.properties.return_gradient for k, v in results_list.items()},
#SEP#            "hessian": {k: v.properties.return_hessian for k, v in results_list.items()},
#SEP#        }
#SEP        #findif_number = [for k, v in results_list.items()
#SEP
#SEP#        for chem, label, imol in calculator_cls.iterate_molecules():
#SEP#            inp = AtomicInput(molecule=imol, **specifications[chem]["specification"])
#SEP#
#SEP#            _, real, bas = delabeler(label)
#SEP#            result = qcng.compute(inp, specifications[chem]["program"])
#SEP#
#SEP#            if not result.success:
#SEP#                print(result.error.error_message)
#SEP#                raise RuntimeError("Calculation did not succeed! Error:\n" + result.error.error_message)
#SEP#
#SEP        component_properties = {}
#SEP        # pull out stuff
#SEP        props = {"energy", "gradient", "hessian"}
#SEP
#SEP        for label, result in results_list.items():
#SEP            component_properties[label] = {}
#SEP
#SEP            for egh in props:
#SEP                if hasattr(result.properties, f"return_{egh}"):
#SEP                    v = getattr(result.properties, f"return_{egh}")
#SEP                    # print(f"  {label} {egh}: {v}")
#SEP                    if v is not None:
#SEP                        component_properties[label][egh] = v
#SEP
#SEP        #findif_number = set(res.extras["qcvars"].get("FINDIF NUMBER") for label, res in results_list.items())
#SEP        findif_number = {label: res.extras["qcvars"].get("FINDIF NUMBER") for label, res in results_list.items()}
#SEP        print(f"{findif_number=}")
#SEP
#SEP        print("\n<<<  (RRR 2) Psi4 ManyBodyComputer.get_psi_results component_properties  >>>")
#SEP        pprint.pprint(component_properties, width=200)
#SEP
#SEP        print("start to analyze")
#SEP        analyze_back = self.qcmb_calculator.analyze(component_properties)
#SEP        analyze_back["nbody_number"] = len(component_properties)
#SEP        print("\n<<<  (RRR 3) Psi4 ManyBodyComputer.get_psi_results analyze_back  >>>")
#SEP        pprint.pprint(analyze_back, width=200)
#SEP
#SEP        nbody_model = self.get_results(external_results=analyze_back, component_results=results_list, client=client)

        nbody_model = self.get_results(client=client)

        # print("\n<<<  (RRR 4) Psi4 ManyBodyComputer.get_psi_results nbody_model  >>>")
        # pprint.pprint(nbody_model.dict(), width=200)
#SEP        print(f"{findif_number=}")

        core.print_out(nbody_model.stdout)
        logger.info('\nQCManyBody:\n' + nbody_model.stdout)
        # v2: logger.debug('\nManyBodyResult QCSchema:\n' + pp.pformat(nbody_model.model_dump()))
        logger.debug('\nManyBodyResult QCSchema:\n' + pp.pformat(nbody_model.dict()))

        qmol = core.Molecule.from_schema(self.molecule.dict())  # nonphy
        wfn = core.Wavefunction.build(qmol, "def2-svp", quiet=True)

#        # TODO all besides nbody may be better candidates for extras than qcvars. energy/gradient/hessian_body_dict in particular are too simple for qcvars (e.g., "2")
        dicts = [
#            #"energies",  # retired
#            #"ptype",     # retired
#            "intermediates",
#            "intermediates_energy",  #"intermediates2",
#            "intermediates_gradient",  #"intermediates_ptype",
#            "intermediates_hessian",  #"intermediates_ptype",
#            "energy_body_dict",
#            "gradient_body_dict",  # ptype_body_dict
#            "hessian_body_dict",  # ptype_body_dict
            "nbody",
#            "cp_energy_body_dict",
#            "nocp_energy_body_dict",
#            "vmfc_energy_body_dict",
#            "cp_gradient_body_dict",
#            "nocp_gradient_body_dict",
#            "vmfc_gradient_body_dict",
#            "cp_hessian_body_dict",
#            "nocp_hessian_body_dict",
#            "vmfc_hessian_body_dict",
        ]

        for qcv, val in nbody_model.extras['qcvars'].items():
            if isinstance(val, dict):
                if qcv in dicts:
                    for qcv2, val2 in val.items():
                        for obj in [core, wfn]:
#                            try:
                                obj.set_variable(str(qcv2), val2)
#                            except ValidationError:
#                                obj.set_variable(f"{self.driver.name} {qcv2}", val2)
            else:
                for obj in [core, wfn]:
                    obj.set_variable(qcv, val)

        # print(f"{self.qcmb_calculator.levels=}")
        mc_ordered_list = list(set(self.qcmb_calculator.levels.values()))

        for mcfragbas, dval in nbody_model.component_properties.items():
            # print(f"{mcfragbas=} {delabeler(mcfragbas)=}")
            mckey, frag, bas = delabeler(mcfragbas)
            # print("DVAL", mckey, frag, bas, dval.dict())
            psi_mcfragbas = new_lab_labeler(mc_ordered_list.index(mckey), frag, bas)
            psi_mcfragbas2 = new_lab_labeler2(mc_ordered_list.index(mckey), frag, bas)  # TODO CHOOSE


            available_properties = dval.dict().keys()
            for obj in [core, wfn]:
                if "return_energy" in available_properties:
                    obj.set_variable(f"N-BODY {psi_mcfragbas2} TOTAL ENERGY", dval.return_energy)
                if "return_gradient" in available_properties:
                    obj.set_variable(f"N-BODY {psi_mcfragbas2} TOTAL GRADIENT", dval.return_gradient)
                if "return_hessian" in available_properties:
                    obj.set_variable(f"N-BODY {psi_mcfragbas2} TOTAL HESSIAN", dval.return_hessian)

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

        # print('WFN')
        # pprint.pprint(wfn.variables(), width=200)
        # print('CORE')
        # pprint.pprint(core.variables(), width=200)

        if return_wfn:
            return (ret, wfn)
        else:
            return ret


def new_lab_labeler(mc_level_idx: str, frag: Tuple, bas:Tuple) -> str:
        return str(mc_level_idx + 1) + "_" + str((tuple(frag), tuple(bas)))

def new_lab_labeler2(mc_level_idx: str, frag: Tuple, bas:Tuple) -> str:
    # returns "1_(1)@(1, 2)"
        #        return mc, ", ".join(map(str, frag)), ", ".join(map(str, bas))
        return f"{mc_level_idx + 1}_({', '.join(map(str, frag))})@({', '.join(map(str, bas))})"


qcvars_to_manybodyproperties = {}
# v2: for skprop in ManyBodyResultProperties.model_fields.keys():
for skprop in ManyBodyResultProperties.__fields__.keys():
    qcvar = skprop.replace("_body", "-body").replace("_corr", "-corr").replace("_", " ").upper()
    qcvars_to_manybodyproperties[qcvar] = skprop
qcvars_to_manybodyproperties["CURRENT ENERGY"] = "return_energy"
qcvars_to_manybodyproperties["CURRENT GRADIENT"] = "return_gradient"
qcvars_to_manybodyproperties["CURRENT HESSIAN"] = "return_hessian"


# TODO where does this live?
def build_manybodyproperties(qcvars: Mapping) -> ManyBodyResultProperties:
    """For results extracted from QC output in QCDB terminology, translate to QCSchema terminology.

    Parameters
    ----------
    qcvars : PreservingDict
        Dictionary of calculation information in QCDB QCVariable terminology.

    Returns
    -------
    atprop : ManyBodyResultProperties
        Object of calculation information in QCSchema ManyBodyResultProperties terminology.

    """
    atprop = {}
    for pv, dpv in qcvars.items():
        if pv in qcvars_to_manybodyproperties:
            atprop[qcvars_to_manybodyproperties[pv]] = dpv

    return ManyBodyResultProperties(**atprop)
