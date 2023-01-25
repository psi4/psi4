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
"""Plan, run, and assemble QC tasks to obtain derivatives by finite difference of lesser derivatives.

===========
FINDIF Flow
===========
Bullet points are major actions
Lines of dashes denote function calls
e/d/dd=dg/g/h := energy, dipole, dipole derivative = dipole gradient, gradient, Hessian

-----------------------------------
FiniteDifferenceComputer.__init__()
-----------------------------------
* collect findif_stencil_size, findif_step_size from initializer kwargs
* BaseComputer.__init__()
* negotiate safety and user wishes on translation and rotation projection

    gradient_from_energies_geometries()
    -----------------------------------
    hessian_from_gradients_geometries()
    -----------------------------------
    hessian_from_energies_geometries()
    ----------------------------------

        _geom_generator()
        -----------------

            _initialize_findif()
            --------------------
            * initialize CdSalcs, partition them per irrep, apply user irreps

        * start the governing dict findifrec with parameters, size, mol
        * for each irrep, for each relevant salc ...

            _displace_cart()
            ----------------
            * form new geometry by linear combination

        * ... and collect geometry into a field of findifrec["displacements"].<label>
        * for (2, 0) also collect off-diagonal displacements
        * also collect undisplaced geometry into field reference
        * return findifrec

* form AtomicComputers for each displacement, particularly changing mol and driver, and possibly relaxing disp symm
* form dict task_list with keys findifrec labels and vals AtomicComputers

--------------------------------------
FiniteDifferenceComputer.build_tasks()
--------------------------------------
* pass

----------------------------------
FiniteDifferenceComputer.compute()
----------------------------------
* compute() for each job in task list

------------------------------------------
FiniteDifferenceComputer.get_psi_results()
------------------------------------------

    Computer.get_results()
    ----------------------

        Computer._prepare_results()
        ---------------------------
        * get_results() for each job in task list
        * arrange atomicresult data into e/d/g/h fields as available on each of reference and displacements entries

            assemble_hessian_from_energies()
            --------------------------------
            assemble_hessian_from_gradients()
            ---------------------------------

                _process_hessian_symmetry_block()
                ---------------------------------
                * resymmetrize each H block

                _process_hessian()
                ------------------
                * transform H to Cartesians and unmasswt

            assemble_gradient_from_energies()
            ---------------------------------
            assemble_dipder_from_dipoles()
            ------------------------------

            * form DD, G, H from lower derivative points

        * place as many of DD, G, H as available onto reference entry

    * pull qcvars off reference job
    * from reference job, set add'l mol, DD, G, H as available
    * form model, including detailed dict at atomicresult.extras["findif_record"]

* convert result to psi4.core.Matrix

    _findif_schema_to_wfn()
    -----------------------
    * build wfn from findif mol and basis (if singular) and module (if singular)
    * push qcvars to P::e and wfn

    gradient_write()
    ----------------
    * write .grad file if requested

    hessian_write()
    ---------------
    * write .hess file if requested

* return dd/g/h and wfn

"""

import copy
import logging
from functools import partial
from typing import Any, Callable, Dict, Iterator, List, Optional, Tuple, Union, TYPE_CHECKING

import numpy as np
from pydantic import Field, validator
from qcelemental.models import DriverEnum, AtomicResult
from qcelemental import constants

from psi4 import core
from psi4.driver import p4util, pp, qcdb, nppp10
from psi4.driver.p4util.exceptions import ValidationError
from psi4.driver.task_base import AtomicComputer, BaseComputer, EnergyGradientHessianWfnReturn

if TYPE_CHECKING:
    import qcportal

logger = logging.getLogger(__name__)

# CONVENTIONS:
# n_ at the start of a variable name is short for "number of."
# _pi at the end of a variable name is short for "per irrep."
# h is the index of an irrep.


def _displace_cart(mass: np.ndarray, geom: np.ndarray, salc_list: core.CdSalcList, i_m: Iterator[Tuple], step_size: float) -> Tuple[np.ndarray, str]:
    """Displace a geometry along the specified displacement SALCs.

    Parameters
    ----------
    mass
        (nat, ) masses [u] of atoms of the molecule (const).
    geom
        (nat, 3) reference geometry [a0] of the molecule (const).
    salc_list
        A list of Cartesian displacement SALCs
    i_m
        An iterator containing tuples. Each tuple has the index of a salc in
        salc_list and the number of steps (positive or negative) to displace
        the salc at that index.
    step_size
        The size of a single "step," i.e., the stencil size.

    Returns
    -------
    disp_geom
        (nat, 3) Displaced geometry.
    label
        Displacement label for the metadata dictionary.

    """
    label = []
    disp_geom = np.copy(geom)
    # This for loop and tuple unpacking is why the function can handle
    # an arbitrary number of SALCs.
    for salc_index, disp_steps in i_m:
        # * Python error if iterate through `salc_list`
        for i in range(len(salc_list[salc_index])):
            component = salc_list[salc_index][i]
            disp_geom[component.atom, component.xyz] += disp_steps * step_size * component.coef / np.sqrt(mass[component.atom])
        label.append(f"{salc_index}: {disp_steps}")

    # salc_index is in descending order. We want the label in ascending order, so...
    # ...add the new label part from the left of the string, not the right.
    label = ', '.join(reversed(label))
    return disp_geom, label


def _initialize_findif(mol: Union["qcdb.Molecule", core.Molecule],
                       freq_irrep_only: int,
                       mode: str,
                       stencil_size: int,
                       step_size: float,
                       initialize_string: Callable,
                       t_project: bool,
                       r_project: bool,
                       initialize: bool,
                       verbose: int = 0) -> Dict:
    """Perform initialization tasks needed by all primary functions.

    Parameters
    ----------
    mol
        The molecule to displace
    freq_irrep_only
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.
    mode
        {"1_0", "2_0", "2_1"}
         The first number specifies the derivative level determined from
         displacements, and the second number is the level determined at.
    stencil_size
        {3, 5}
        Number of points to evaluate for each displacement basis vector inclusive of central reference geometry.
    step_size
        [a0]
    initialize_string
         A function that returns the string to print to show the caller was entered.
         The string is both caller-specific and dependent on values determined
         in this function.
    initialize
        For printing, whether call is from generator or assembly stages.
    verbose
         Set to 0 to silence extra print information, regardless of the print level.
         Used so the information is printed only during geometry generation, and not
         during the derivative computation as well.

    Returns
    -------
    data
        Miscellaneous information required by callers.
    """

    info = """
         ----------------------------------------------------------
                                   FINDIF
                     R. A. King and Jonathon Misiewicz
         ----------------------------------------------------------

"""
    if initialize:
        core.print_out(info)
        logger.info(info)

    print_lvl = core.get_option("FINDIF", "PRINT")

    data = {"print_lvl": print_lvl, "stencil_size": stencil_size, "step_size": step_size}

    if print_lvl:
        info = initialize_string(data)
        core.print_out(info)
        logger.info(info)

    # Get settings for CdSalcList, then get the CdSalcList.
    method_allowed_irreps = 0x1 if mode == "1_0" else 0xFF
    # core.get_option returns an int, but CdSalcList expect a bool, so re-cast
    salc_list = core.CdSalcList(mol, method_allowed_irreps, t_project, r_project)

    n_atom = mol.natom()
    n_irrep = salc_list.nirrep()
    n_salc = salc_list.ncd()

    if print_lvl and verbose:
        info = f"    Number of atoms is {n_atom}.\n"
        if method_allowed_irreps != 0x1:
            info += f"    Number of irreps is {n_irrep}.\n"
        info += "    Number of {!s}SALCs is {:d}.\n".format("" if method_allowed_irreps != 0x1 else "symmetric ",
                                                            n_salc)
        info += f"    Translations projected? {t_project:d}. Rotations projected? {r_project:d}.\n"
        core.print_out(info)
        logger.info(info)

    # TODO: Replace with a generator from a stencil to a set of points.
    # Diagonal displacements differ between the totally symmetric irrep, compared to all others.
    # Off-diagonal displacements are the same for both.
    pts_dict = {
        3: {
            "sym_irr": ((-1, ), (1, )),
            "asym_irr": ((-1, ), ),
            "off": ((1, 1), (-1, -1))
        },
        5: {
            "sym_irr": ((-2, ), (-1, ), (1, ), (2, )),
            "asym_irr": ((-2, ), (-1, )),
            "off": ((-1, -2), (-2, -1), (-1, -1), (1, -1), (-1, 1), (1, 1), (2, 1), (1, 2))
        }
    }

    try:
        disps = pts_dict[stencil_size]
    except KeyError:
        raise ValidationError(f"FINDIF: Number of points ({stencil_size}) not among {pts_dict.keys()}!")

    # Convention: x_pi means x_per_irrep. The ith element is x for irrep i, with Cotton ordering.
    salc_indices_pi = [[] for h in range(n_irrep)]

    # Validate that we have an irrep matching the user-specified irrep, if any.
    try:
        salc_indices_pi[freq_irrep_only]
    except (TypeError, IndexError):
        if freq_irrep_only != -1:
            raise ValidationError(
                f"FINDIF: 0-indexed Irrep value ({freq_irrep_only}) not in valid range: <{len(salc_indices_pi)}.")

    # Populate salc_indices_pi for all irreps.
    # * Python error if iterate through `salc_list`
    for i in range(len(salc_list)):
        salc_indices_pi[salc_list[i].irrep_index()].append(i)

    # If the method allows more than one irrep, print how the irreps partition the SALCS.
    if print_lvl and method_allowed_irreps != 0x1 and verbose:
        info = "    Index of SALCs per irrep:\n"
        for h in range(n_irrep):
            if print_lvl > 1 or freq_irrep_only in {h, -1}:
                tmp = (" {:d} " * len(salc_indices_pi[h])).format(*salc_indices_pi[h])
                info += "      {:d} : ".format(h + 1) + tmp + "\n"
        info += "    Number of SALCs per irrep:\n"
        for h in range(n_irrep):
            if print_lvl > 1 or freq_irrep_only in {h, -1}:
                info += "      Irrep {:d}: {:d}\n".format(h + 1, len(salc_indices_pi[h]))
        core.print_out(info)
        logger.info(info)

    # Now that we've printed the SALCs, clear any that are not of user-specified symmetry.
    if freq_irrep_only != -1:
        for h in range(n_irrep):
            if h != freq_irrep_only:
                salc_indices_pi[h].clear()

    n_disp_pi = []

    for irrep, indices in enumerate(salc_indices_pi):
        n_disp = len(indices) * len(disps["asym_irr" if irrep != 0 else "sym_irr"])
        if mode == "2_0":
            # Either len(indices) or len(indices)-1 is even, so dividing by two is safe.
            n_disp += len(indices) * (len(indices) - 1) // 2 * len(disps["off"])
        n_disp_pi.append(n_disp)

    # Let's print out the number of geometries, the displacement multiplicity, and the CdSALCs!
    if print_lvl and verbose:
        info = f"    Number of geometries (including reference) is {sum(n_disp_pi) + 1}.\n"
        if method_allowed_irreps != 0x1:
            info += "    Number of displacements per irrep:\n"
            for i, ndisp in enumerate(n_disp_pi, start=1):
                info += f"      Irrep {i}: {ndisp}\n"
        core.print_out(info)
        logger.info(info)

    if print_lvl > 1 and verbose:
        for i in range(len(salc_list)):
            salc_list[i].print_out()

    data.update({
        "n_disp_pi": n_disp_pi,
        "n_irrep": n_irrep,
        "n_salc": n_salc,
        "n_atom": n_atom,
        "salc_list": salc_list,
        "salc_indices_pi": salc_indices_pi,
        "disps": disps,
        "project_translations": t_project,
        "project_rotations": r_project
    })

    return data


def _geom_generator(mol: Union["qcdb.Molecule", core.Molecule], freq_irrep_only: int, mode: str, *, t_project: bool = True, r_project: bool = True, stencil_size: int = 3, step_size: float = 0.005) -> Dict:
    """
    Generate geometries for the specified molecule and derivative levels.
    You probably want to instead use one of the convenience functions:
    gradient_from_energies_geometries, hessian_from_energies_geometries,
    hessian_from_gradients_geometries.

    Parameters
    ----------
    mol
        The molecule on which to perform a finite difference calculation.
    freq_irrep_only
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps. Irrelevant for "1_0".
    mode
        {"1_0", "2_0", "2_1"}
        The first number specifies the targeted derivative level. The
        second number is the compute derivative level. E.g., "2_0"
        is hessian from energies.
    stencil_size
        {3, 5}
        Number of points to evaluate for each displacement basis vector inclusive of central reference geometry.
    step_size
        Displacement size [a0].

    Returns
    -------
    findifrec
        Dictionary of finite difference data, specified below.
        The dictionary makes findifrec _extensible_. If you need a new field
        in the record, just add it.
        All fields should be present at all times, with two exceptions:
            1. Fields for computed quantities will not be available until
               after they are computed.
            2. Displacement specific overrides for globals will not be
               available unless the user specified the overrides.
               (Such overrides are not implemented at time of writing. An example
               is giving a displacement its own step dict.)

    step : dict
        A descriptor for the finite difference step.
        In future, this can be overriden by step fields for individual displacements.

        units : {'Bohr'}
            The units for the displacement. The code currently assumes "bohr," per MolSSI standards.
        size : float
            The step size for the displacement.

    stencil_size : {3, 5}
        Number of points to evaluate at for each displacement basis vector. Count
        includes the central reference point.

    displacement_space : {'CdSalc'}
        A string specifying the vector space in which displacements are performed.
        Currently, only CdSalc is supported.

    project_translations : bool
        Whether translations are to be projected out of the displacements.

    project_rotations : bool
        Whether rotations are to be projected out of the displacements.

    molecule : dict
        The reference molecule, in MolSSI schema. See
        https://molssi-qc-schema.readthedocs.io/en/latest/auto_topology.html

    displacements : dict
        A dictionary mapping labels specifying the displacement to data about
        the geometry. Labels are of the form "A: a, B: b" where A and B index the
        basis vector in displacement space and A < B, and a and b index the step
        magnitude. For instance, "0: 1, 1: -1" specifies displacing +1 in
        displacement vector 0 and -1 in displacement vector 1. "1: -1, 0: 1" is
        forbidden for breaking ordering. Generalizes to arbitrary numbers of
        simultaneous displacements in the obvious way.

        The possible geometry data is as follows:

        geometry: list of floats
            (3 * nat) The molecular geometry as a flat list in bohr. All coordinates
            are given for one atom before proceeding to the next atom.

        energy: int
            The last computed electronic energy at the geometry.

        gradient: list of floats
            (3 * nat) The last computed gradient of energy with respect to changes in
            geometry at the geometry, as a flat list. All coordinates are given for
            displacing one atom before proceeding to the next atom.

    reference : dict
         A geometry data dict, as described above, for the reference geometry.
    """

    msg_dict = {
        "1_0":
        "energies to determine gradients",
        "2_1":
        "gradients to determine vibrational frequencies and \n"
        "  normal modes. Resulting frequencies are only valid at stationary points",
        "2_0":
        "gradients to determine vibrational frequencies and \n"
        "  normal modes. Resulting frequencies are only valid at stationary points"
    }

    try:
        print_msg = msg_dict[mode]
    except KeyError:
        raise ValidationError("FINDIF: Mode {} not recognized.".format(mode))

    def init_string(data):
        return f"""  Using finite-differences of {print_msg}.
    Generating geometries for use with {data["stencil_size"]}-point formula.
    Displacement size will be {data["step_size"]:6.2e}.\n"""

    # Genuine support for qcdb molecules would be nice. But that requires qcdb CdSalc tech.
    # Until then, silently swap the qcdb molecule out for a psi4.core.Molecule.
    if isinstance(mol, qcdb.Molecule):
        mol = core.Molecule.from_dict(mol.to_dict())

    data = _initialize_findif(mol, freq_irrep_only, mode, stencil_size, step_size, init_string, t_project, r_project,
                              True, 1)

    # We can finally start generating displacements.
    ref_geom = np.array(mol.geometry())

    # Now we generate the metadata...
    findifrec = {
        "step": {
            "units": "bohr",
            "size": data["step_size"]
        },
        "stencil_size": data["stencil_size"],
        "displacement_space": "CdSALC",
        "project_translations": data["project_translations"],
        "project_rotations": data["project_rotations"],
        "molecule": mol.to_schema(dtype=2, units='Bohr'),
        "displacements": {},
        "reference": {}
    }

    def append_geoms(indices, steps):
        """Given a list of indices and a list of steps to displace each, append the corresponding geometry to the list."""

        # Next, to make this salc/magnitude composite.
        disp_geom, label = _displace_cart(findifrec['molecule']['masses'], ref_geom, data["salc_list"],
                                          zip(indices, steps), data["step_size"])
        if data["print_lvl"] > 2:
            info = "\nDisplacement '{}'\n{}\n".format(label, nppp10(disp_geom))
            core.print_out(info)
            logger.info(info)
        findifrec["displacements"][label] = {"geometry": disp_geom}

    for h in range(data["n_irrep"]):
        active_indices = data["salc_indices_pi"][h]

        for index in active_indices:
            # Displace along the diagonal.
            # Remember that the totally symmetric irrep has special displacements.
            for val in data["disps"]["sym_irr" if h == 0 else "asym_irr"]:
                append_geoms((index, ), val)

        # Hessian from energies? We have off-diagonal displacements to worry about.
        if mode == "2_0":
            # i indexes SALC indices of the current irrep.
            for i, index in enumerate(active_indices):
                for index2 in active_indices[:i]:
                    for val in data["disps"]["off"]:
                        append_geoms((index, index2), val)

    if data["print_lvl"] > 2:
        logger.info("\nReference\n{}\n".format(nppp10(ref_geom)))
    findifrec["reference"]["geometry"] = ref_geom

    if data["print_lvl"] > 1:
        logger.info("\n-------------------------------------------------------------")

    return findifrec


_der_from_lesser_docstring = """
    Parameters
    ----------
    mol
        The molecule on which to perform a finite difference calculation.
    freq_irrep_only
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps. Irrelevant for "1_0".
    stencil_size
        {3, 5}
        Number of points to evaluate for each displacement basis vector inclusive of central reference geometry.
    step_size
        Displacement size [a0].

    Returns
    -------
    findifrec : dict
        Dictionary of finite difference data, specified in _geom_generator docstring.

"""


gradient_from_energies_geometries = partial(_geom_generator, freq_irrep_only=-1, mode="1_0")
hessian_from_gradients_geometries = partial(_geom_generator, mode="2_1")
hessian_from_energies_geometries = partial(_geom_generator, mode="2_0")

gradient_from_energies_geometries.__doc__ = "Generate geometries for a gradient by finite difference of energies." + _der_from_lesser_docstring
hessian_from_gradients_geometries.__doc__ = "Generate geometries for a Hessian by finite difference of gradients." + _der_from_lesser_docstring
hessian_from_energies_geometries.__doc__ = "Generate geometries for a Hessian by finite difference of energies." + _der_from_lesser_docstring


def assemble_gradient_from_energies(findifrec: Dict) -> np.ndarray:
    """Compute the gradient by finite difference of energies.

    Parameters
    ----------
    findifrec
        Dictionary of finite difference data, specified in _geom_generator docstring.

    Returns
    -------
    gradient
        (nat, 3) Cartesian gradient [Eh/a0].
    """

    # This *must* be a Psi molecule at present - CdSalcList generation panics otherwise
    mol = core.Molecule.from_schema(findifrec["molecule"], nonphysical=True, verbose=0)

    def init_string(data):
        return f"""  Computing gradient from energies.
    Using {findifrec["stencil_size"]}-point formula.
    Energy without displacement: {findifrec["reference"]["energy"]:15.10f}
    Check energies below for precision!
    Forces are for mass-weighted, symmetry-adapted cartesians [a0].\n"""

    data = _initialize_findif(mol, -1, "1_0", findifrec['stencil_size'], findifrec['step']['size'], init_string,
                              findifrec['project_translations'], findifrec['project_rotations'], False)
    salc_indices = data["salc_indices_pi"][0]

    # Extract the energies, and turn then into an ndarray for easy manipulating
    # E(i, j) := Energy on displacing the ith SALC we care about in the jth step
    # Steps are ordered, for example, -2, -1, 1, 2
    max_disp = (findifrec["stencil_size"] - 1) // 2  # The numerator had better be divisible by two.
    e_per_salc = 2 * max_disp
    E = np.zeros((len(salc_indices), e_per_salc))

    for i, salc_index in enumerate(salc_indices):
        for j in range(1, max_disp + 1):
            E[i, max_disp - j] = findifrec["displacements"][f"{salc_index}: {-j}"]["energy"]
            E[i, max_disp + j - 1] = findifrec["displacements"][f"{salc_index}: {j}"]["energy"]

    # Perform the finite difference.
    if findifrec["stencil_size"] == 3:
        g_q = (E[:, 1] - E[:, 0]) / (2.0 * findifrec["step"]["size"])
    elif findifrec["stencil_size"] == 5:
        g_q = (E[:, 0] - 8.0 * E[:, 1] + 8.0 * E[:, 2] - E[:, 3]) / (12.0 * findifrec["step"]["size"])
    else:  # This error SHOULD have already been caught, but just in case...
        raise ValidationError("FINDIF: {} is an invalid number of points.".format(findifrec["stencil_size"]))
    g_q = np.asarray(g_q)

    if data["print_lvl"]:
        energy_string = ""
        for i in range(1, max_disp + 1):
            energy_string = f"Energy(-{i})        " + energy_string + f"Energy(+{i})        "
        info = "\n     Coord      " + energy_string + "    Force"
        for salc in range(data["n_salc"]):
            print_str = "\n    {:5d}" + " {:17.10f}" * (e_per_salc) + " {force:17.10f}"
            energies = E[salc]
            info += print_str.format(salc, force=g_q[salc], *energies)
        core.print_out(info)
        logger.info(info)

    # Transform the gradient from mass-weighted SALCs to non-mass-weighted Cartesians
    B = data["salc_list"].matrix()
    g_cart = np.dot(g_q, B)
    g_cart = g_cart.reshape(data["n_atom"], 3)
    massweighter = np.array([mol.mass(a) for a in range(data["n_atom"])])**(0.5)
    g_cart = (g_cart.T * massweighter).T

    if data["print_lvl"]:
        info = "\n         -------------------------------------------------------------\n"
        core.print_out(info)
        logger.info(info)

    return g_cart


def _process_hessian_symmetry_block(H_block: np.ndarray, B_block: np.ndarray, massweighter: np.ndarray, irrep: str, print_lvl: int) -> np.ndarray:
    """Perform post-construction processing for a symmetry block of the Hessian.
       Statements need to be printed, and the Hessian must be made orthogonal.

    Parameters
    ----------
    H_block
        A block of the Hessian for an irrep, in mass-weighted salcs.
        (nsalc, nsalc)
    B_block
        A block of the B matrix for an irrep, which transforms CdSalcs to Cartesians.
        (nsalc, 3 * nat)
    massweighter
        The mass associated with each atomic coordinate.
        (3 * nat, ) Due to x, y, z, values appear in groups of three.
    irrep
        A string identifying the irrep H_block and B_block are of.
    print_lvl
        The level of printing information requested by the user.

    Returns
    -------
    H_block
        H_block, but made into an orthogonal array.
    """

    # Symmetrize our Hessian block.
    # The symmetric structure is lost due to errors in the computation
    H_block = (H_block + H_block.T) / 2.0

    if print_lvl >= 3:
        core.print_out(f"Force Constants for irrep {irrep} in mass-weighted, symmetry-adapted Cartesian coordinates.")
        core.print_out("\n{}\n".format(nppp10(H_block)))

    evals, evects = np.linalg.eigh(H_block)
    # Get our eigenvalues and eigenvectors in descending order.
    idx = evals.argsort()[::-1]
    evals = evals[idx]
    evects = evects[:, idx]

    normal_irr = np.dot((B_block * massweighter).T, evects)

    if print_lvl >= 2:
        core.print_out("\n    Normal coordinates (non-mass-weighted) for irrep {}:\n".format(irrep))
        core.print_out("\n{}\n".format(nppp10(normal_irr)))

    return H_block


def _process_hessian(H_blocks: List[np.ndarray], B_blocks: List[np.ndarray], massweighter: np.ndarray, print_lvl: int) -> np.ndarray:
    """Perform post-construction processing for the Hessian.
       Statements need to be printed, and the Hessian must be transformed.

    Parameters
    ----------
    H_blocks
        A list of blocks of the Hessian per irrep, in mass-weighted salcs.
        Each is (nsalc_in_irrep, nsalc_in_irrep)
    B_blocks
        A block of the B matrix per irrep, which transforms CdSalcs to Cartesians.
        Each is (nsalc_in_irrep, 3 * nat)
    massweighter
        The mass associated with each atomic coordinate.
        (3 * nat, ) Due to x, y, z, values appear in groups of three.
    print_lvl
        The level of printing information requested by the user.

    Returns
    -------
    Hx
        The Hessian in non-mass weighted cartesians.
    """

    # Handle empty case (atom)
    if not H_blocks and not B_blocks:
        nat3 = massweighter.size
        return np.zeros((nat3, nat3), dtype=np.float64)

    # We have the Hessian in each irrep! The final task is to perform coordinate transforms.
    H = p4util.block_diagonal_array(*H_blocks)
    B = np.vstack(B_blocks)

    if print_lvl >= 3:
        core.print_out("\n    Force constant matrix for all computed irreps in mass-weighted SALCS.\n")
        core.print_out("\n{}\n".format(nppp10(H)))

    # Transform the massweighted Hessian from the CdSalc basis to Cartesians.
    # The Hessian is the matrix not of a linear transformation, but of a (symmetric) bilinear form
    # As such, the change of basis is formula A' = Xt A X, no inverses!
    # More conceptually, it's A'_kl = A_ij X_ik X_jl; Each index transforms linearly.
    Hx = np.dot(np.dot(B.T, H), B)
    if print_lvl >= 3:
        core.print_out("\n    Force constants in mass-weighted Cartesian coordinates.\n")
        core.print_out("\n{}\n".format(nppp10(Hx)))

    # Un-massweight the Hessian.
    Hx = np.transpose(Hx / massweighter) / massweighter

    if print_lvl >= 3:
        core.print_out("\n    Force constants in Cartesian coordinates.\n")
        core.print_out("\n{}\n".format(nppp10(Hx)))

    if print_lvl:
        core.print_out("\n-------------------------------------------------------------\n")

    return Hx


def assemble_dipder_from_dipoles(findifrec: Dict, freq_irrep_only: int) -> np.ndarray:
    """Compute the dipole derivatives by finite difference of dipoles.

    Parameters
    ----------
    findifrec
        Dictionary of finite difference data, specified in _geom_generator docstring.
    freq_irrep_only
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.

    Returns
    -------
    dipder
        (3 * nat, 3) Cartesian Dipole Derivatives [Eh/a0^2]

    """

    # This *must* be a Psi molecule at present - CdSalcList generation panics otherwise
    mol = core.Molecule.from_schema(findifrec["molecule"], nonphysical=True, verbose=0)

    pg = mol.point_group()
    ct = pg.char_table()
    order = pg.order()

    displacements = findifrec["displacements"]

    def init_string(data):
        return ("")

    data = _initialize_findif(mol, freq_irrep_only, "2_1", findifrec['stencil_size'], findifrec['step']['size'],
                              init_string, findifrec['project_translations'], findifrec['project_rotations'], False)
    salc_indices = data["salc_indices_pi"][0]
    max_disp = (findifrec["stencil_size"] - 1) // 2  # The numerator had better be divisible by two.
    d_per_salc = 2 * max_disp

    # Populating with positive and negative displacements for the identity point group
    dipole = np.zeros(shape=(data['n_salc'], d_per_salc, 3))
    for salc_index in salc_indices:
        for j in range(1, max_disp + 1):
            dipole[salc_index, max_disp - j] = displacements[f"{salc_index}: {-j}"]["dipole"]
            dipole[salc_index, max_disp + j - 1] = displacements[f"{salc_index}: {j}"]["dipole"]

    for h in range(1, data["n_irrep"]):
        # Find the group operation that converts + to - displacements.
        gamma = ct.gamma(h)
        for group_op in range(order):
            if gamma.character(group_op) == -1:
                break
        else:
            raise ValidationError("A symmetric dipole passed for a non-symmetric one.")

        sym_op = np.array(ct.symm_operation(group_op).matrix())
        salc_indices = data["salc_indices_pi"][h]

        # Creating positive displacements and populating for the other point groups
        for salc_index in salc_indices:
            for j in range(1, max_disp + 1):
                pos_disp_dipole = np.dot(sym_op, displacements[f"{salc_index}: {-j}"]["dipole"].T)
                dipole[salc_index, max_disp - j] = displacements[f"{salc_index}: {-j}"]["dipole"]
                dipole[salc_index, max_disp + j - 1] = pos_disp_dipole

    # Computing the dipole derivative by finite differnce
    if findifrec["stencil_size"] == 3:
        dipder_q = (dipole[:, 1] - dipole[:, 0]) / (2.0 * findifrec["step"]["size"])
    elif findifrec["stencil_size"] == 5:
        dipder_q = (dipole[:, 0] - 8.0 * dipole[:, 1] + 8.0 * dipole[:, 2] -
                    dipole[:, 3]) / (12.0 * findifrec["step"]["size"])

    # Transform the dipole derivates from mass-weighted SALCs to non-mass-weighted Cartesians
    B = np.asarray(data["salc_list"].matrix())
    dipder_cart = np.dot(dipder_q.T, B)
    dipder_cart = dipder_cart.T.reshape(data["n_atom"], 9)

    massweighter = np.array([mol.mass(a) for a in range(data["n_atom"])])**(0.5)
    dipder_cart = (dipder_cart.T * massweighter).T

    dipder_cart = dipder_cart.reshape(3 * data["n_atom"], 3)
    return dipder_cart


def assemble_hessian_from_gradients(findifrec: Dict, freq_irrep_only: int) -> np.ndarray:
    """Compute the Hessian by finite difference of gradients.

    Parameters
    ----------
    findifrec
        Dictionary of finite difference data, specified in _geom_generator docstring.
    freq_irrep_only
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.

    Returns
    -------
    hessian
        (3 * nat, 3 * nat) Cartesian Hessian [Eh/a0^2]
    """

    # This *must* be a Psi molecule at present - CdSalcList generation panics otherwise
    mol = core.Molecule.from_schema(findifrec["molecule"], nonphysical=True, verbose=0)

    displacements = findifrec["displacements"]

    def init_string(data):
        return ("  Computing second-derivative from gradients using projected, \n"
                "  symmetry-adapted, cartesian coordinates.\n\n"
                "  {:d} gradients passed in, including the reference geometry.\n".format(len(displacements) + 1))

    data = _initialize_findif(mol, freq_irrep_only, "2_1", findifrec['stencil_size'], findifrec['step']['size'],
                              init_string, findifrec['project_translations'], findifrec['project_rotations'], False)

    # For non-totally symmetric CdSALCs, a symmetry operation can convert + and - displacements.
    # Good News: By taking advantage of that, we (potentially) ran less computations.
    # Bad News: We need to find the - displacements from the + computations now.
    # The next ~80 lines of code are dedicated to that task.
    if data["print_lvl"]:
        core.print_out("  Generating complete list of displacements from unique ones.\n\n")

    pg = mol.point_group()
    ct = pg.char_table()
    order = pg.order()

    # Determine what atoms map to what other atoms under the point group operations.
    # The py-side compute_atom_map will work whether mol is a Py-side or C-side object.
    atom_map = qcdb.compute_atom_map(mol)
    if data["print_lvl"] >= 3:
        core.print_out("    The atom map:\n")
        for atom, sym_image_list in enumerate(atom_map):
            core.print_out(f"     {atom + 1:d} : ")
            for image_atom in sym_image_list:
                core.print_out(f"{image_atom + 1:4d}")
            core.print_out("\n")
        core.print_out("\n")

    # A list of lists of gradients, per irrep
    gradients_pi = [[]]
    # Extract and print the symmetric gradients. These need no additional processing.
    max_disp = (findifrec["stencil_size"] - 1) // 2  # The numerator had better be divisible by two.
    for i in data["salc_indices_pi"][0]:
        for n in range(-max_disp, 0):
            grad_raw = displacements[f"{i}: {n}"]["gradient"]
            gradients_pi[0].append(np.reshape(grad_raw, (-1, 3)))
        for n in range(1, max_disp + 1):
            grad_raw = displacements[f"{i}: {n}"]["gradient"]
            gradients_pi[0].append(np.reshape(grad_raw, (-1, 3)))

    if data["print_lvl"] >= 3:
        core.print_out("    Symmetric gradients\n")
        for gradient in gradients_pi[0]:
            core.print_out("\n{}\n".format(nppp10(gradient)))

    # Asymmetric gradient. There's always SOME operation that transforms a positive
    # into a negative displacement.By doing extra things here, we can find the
    # gradients at the positive displacements.
    for h in range(1, data["n_irrep"]):

        # If there are no CdSALCs in this irrep, let's skip it.
        if not data["n_disp_pi"][h]:
            gradients_pi.append([])
            continue

        gamma = ct.gamma(h)
        if data["print_lvl"] >= 3:
            core.print_out(f"Characters for irrep {h}\n")
            for group_op in range(order):
                core.print_out(" {:5.1f}".format(gamma.character(group_op)))
            core.print_out("\n")

        # Find the group operation that converts + to - displacements.
        for group_op in range(order):
            if gamma.character(group_op) == -1:
                break
        else:
            raise ValidationError("A symmetric gradient passed for a non-symmetric one.")
        if data["print_lvl"]:
            core.print_out("    Operation {} takes plus displacements of irrep {} to minus ones.\n".format(
                group_op + 1, gamma.symbol()))

        sym_op = np.array(ct.symm_operation(group_op).matrix())
        gradients = []

        def recursive_gradients(i, n):
            """Populate gradients, with step -n, -n+1, ... -1, 1, ... n. Positive displacements are computed."""

            grad_raw = displacements[f"{i}: {-n}"]["gradient"]
            gradients.append(np.reshape(grad_raw, (-1, 3)))
            new_grad = np.zeros((data["n_atom"], 3))
            for atom, image in enumerate(atom_map):
                atom2 = image[group_op]
                new_grad[atom2] = np.einsum("xy,y->x", sym_op, gradients[-1][atom])
            if n > 1:
                recursive_gradients(i, n - 1)
            gradients.append(new_grad)

        for i in data["salc_indices_pi"][h]:
            recursive_gradients(i, max_disp)
        gradients_pi.append(gradients)

    # Massweight all gradients.
    # Remember, the atom currently corresponds to our 0 axis, hence these transpose tricks.
    massweighter = np.asarray([mol.mass(a) for a in range(data["n_atom"])])**(-0.5)
    gradients_pi = [[(grad.T * massweighter).T for grad in gradients] for gradients in gradients_pi]

    if data["print_lvl"] >= 3:
        core.print_out("    All mass-weighted gradients\n")
        for gradients in gradients_pi:
            for grad in gradients:
                core.print_out("\n{}\n".format(nppp10(grad)))

    # We have all our gradients generated now!
    # Next, time to get our Hessian.

    H_pi = []
    B_pi = []
    irrep_lbls = mol.irrep_labels()
    massweighter = np.repeat(massweighter, 3)

    for h in range(data["n_irrep"]):
        n_disp = data["n_disp_pi"][h]
        Nindices = len(data["salc_indices_pi"][h])
        gradients = gradients_pi[h]

        if not Nindices:
            continue

        # Flatten each gradient, and turn it into a COLUMN of the matrix.
        gradient_matrix = np.array([grad.flatten() for grad in gradients]).T
        # Transform disps from Cartesian to CdSalc coordinates.
        # For future convenience, we transpose.
        # Rows are gradients and columns are coordinates with respect to a particular CdSALC.
        B_pi.append(data["salc_list"].matrix_irrep(h))
        grads_adapted = np.dot(B_pi[-1], gradient_matrix).T

        if data["print_lvl"] >= 3:
            core.print_out("Gradients in B-matrix coordinates\n")
            for disp in range(n_disp):
                core.print_out(f" disp {disp}: ")
                for salc in grads_adapted[disp]:
                    core.print_out(f"{salc:15.10f}")
                core.print_out("\n")

        H_pi.append(np.empty([Nindices, Nindices]))

        if findifrec["stencil_size"] == 3:
            H_pi[-1] = (grads_adapted[1::2] - grads_adapted[::2]) / (2.0 * findifrec["step"]["size"])
        elif findifrec["stencil_size"] == 5:
            H_pi[-1] = (grads_adapted[::4] - 8 * grads_adapted[1::4] + 8 * grads_adapted[2::4] -
                        grads_adapted[3::4]) / (12.0 * findifrec["step"]["size"])

        H_pi[-1] = _process_hessian_symmetry_block(H_pi[-1], B_pi[-1], massweighter, irrep_lbls[h], data["print_lvl"])

    # All blocks of the Hessian are now constructed!
    return _process_hessian(H_pi, B_pi, massweighter, data["print_lvl"])


def assemble_hessian_from_energies(findifrec: Dict, freq_irrep_only: int) -> np.ndarray:
    """Compute the Hessian by finite difference of energies.

    Parameters
    ----------
    findifrec
        Dictionary of finite difference data, specified in _geom_generator docstring.
    freq_irrep_only
        The 0-indexed Cotton ordered irrep to get frequencies for. Choose -1 for all irreps.

    Returns
    -------
    hessian
        (3 * nat, 3 * nat) Cartesian Hessian [Eh/a0^2].
    """

    # This *must* be a Psi molecule at present - CdSalcList generation panics otherwise
    mol = core.Molecule.from_schema(findifrec["molecule"], nonphysical=True, verbose=0)

    displacements = findifrec["displacements"]
    ref_energy = findifrec["reference"]["energy"]

    def init_string(data):
        max_label_len = str(max([len(label) for label in displacements], default=3))
        out_str = ""
        for label, disp_data in displacements.items():
            out_str += ("    {:" + max_label_len + "s} : {:20.10f}\n").format(label, disp_data["energy"])
        return ("  Computing second-derivative from energies using projected, \n"
                "  symmetry-adapted, cartesian coordinates.\n\n"
                "  {:d} energies passed in, including the reference geometry.\n"
                "    Using {:d}-point formula.\n"
                "    Energy without displacement: {:15.10f}\n"
                "    Check energies below for precision!\n{}".format(
                    len(displacements) + 1, findifrec["stencil_size"], ref_energy, out_str))

    data = _initialize_findif(mol, freq_irrep_only, "2_0", findifrec['stencil_size'], findifrec['step']['size'],
                              init_string, findifrec['project_translations'], findifrec['project_rotations'], False)

    massweighter = np.repeat([mol.mass(a) for a in range(data["n_atom"])], 3)**(-0.5)
    B_pi = []
    H_pi = []
    irrep_lbls = mol.irrep_labels()
    max_disp = (findifrec["stencil_size"] - 1) // 2
    e_per_diag = 2 * max_disp

    # Unlike in the gradient case, we have no symmetry transformations to worry about.
    # We get to the task directly: assembling the force constants in each irrep block.
    for h in range(data["n_irrep"]):
        salc_indices = data["salc_indices_pi"][h]
        if not salc_indices: continue

        n_salcs = len(salc_indices)
        E = np.zeros((len(salc_indices), e_per_diag))

        # Step One: Diagonals
        # For asymmetric irreps, the energy at a + disp is the same as at a - disp
        # Just reuse the - disp energy for the + disp energy

        for i, salc_index in enumerate(salc_indices):
            for j in range(1, max_disp + 1):
                E[i, max_disp - j] = displacements[f"{salc_index}: {-j}"]["energy"]
                k = -j if h else j  # Because of the +- displacement trick
                E[i, max_disp + j - 1] = displacements[f"{salc_index}: {k}"]["energy"]
        # Now determine all diagonal force constants for this irrep.
        if findifrec["stencil_size"] == 3:
            diag_fcs = E[:, 0] + E[:, 1]
            diag_fcs -= 2 * ref_energy
            diag_fcs /= (findifrec["step"]["size"]**2)
        elif findifrec["stencil_size"] == 5:
            diag_fcs = -E[:, 0] + 16 * E[:, 1] + 16 * E[:, 2] - E[:, 3]
            diag_fcs -= 30 * ref_energy
            diag_fcs /= (12 * findifrec["step"]["size"]**2)
        H_irr = np.diag(diag_fcs)

        # TODO: It's a bit ugly to use the salc indices to grab the off-diagonals but the indices
        # within the irrep to grab the diagonals. Is there a better way to do this?

        # Step Two: Off-diagonals
        # We need off-diagonal energies, diagonal energies, AND the reference energy
        # Grabbing off-diagonal energies is a pain, so once we know our SALCs...
        # ...define offdiag_en to do that for us.
        for i, salc in enumerate(salc_indices):
            for j, salc2 in enumerate(salc_indices[:i]):
                offdiag_en = lambda index: displacements["{l}: {}, {k}: {}".format(
                    k=salc, l=salc2, *data["disps"]["off"][index])]["energy"]
                if findifrec["stencil_size"] == 3:
                    fc = (+offdiag_en(0) + offdiag_en(1) + 2 * ref_energy - E[i][0] - E[i][1] - E[j][0] -
                          E[j][1]) / (2 * findifrec["step"]["size"]**2)
                elif findifrec["stencil_size"] == 5:
                    fc = (-offdiag_en(0) - offdiag_en(1) + 9 * offdiag_en(2) - offdiag_en(3) - offdiag_en(4) +
                          9 * offdiag_en(5) - offdiag_en(6) - offdiag_en(7) + E[i][0] - 7 * E[i][1] - 7 * E[i][2] +
                          E[i][3] + E[j][0] - 7 * E[j][1] - 7 * E[j][2] + E[j][3] +
                          12 * ref_energy) / (12 * findifrec["step"]["size"]**2)
                H_irr[i, j] = fc
                H_irr[j, i] = fc

        B_pi.append(data["salc_list"].matrix_irrep(h))
        H_pi.append(_process_hessian_symmetry_block(H_irr, B_pi[-1], massweighter, irrep_lbls[h], data["print_lvl"]))

    # All blocks of the Hessian are now constructed!
    return _process_hessian(H_pi, B_pi, massweighter, data["print_lvl"])


class FiniteDifferenceComputer(BaseComputer):

    molecule: Any
    driver: DriverEnum
    metameta: Dict[str, Any] = {}
    task_list: Dict[str, BaseComputer] = {}
    findifrec: Dict[str, Any] = {}
    computer: BaseComputer = AtomicComputer
    method: str

    @validator('driver')
    def set_driver(cls, driver):
        egh = ['energy', 'gradient', 'hessian']
        if driver not in egh:
            raise ValidationError(f"""Wrapper is unhappy to be calling function ({driver}) not among {egh}.""")

        return driver

    @validator('molecule')
    def set_molecule(cls, mol):
        mol.update_geometry()
        mol.fix_com(True)
        mol.fix_orientation(True)
        return mol

    def __init__(self, **data):
        """Initialize FiniteDifference class.

        data keywords include
        * general AtomicInput keys like molecule, driver, method, basis, and keywords.
        * specialized findif keys like findif_mode, findif_irrep, and those converted from keywords to kwargs:
          findif_stencil_size, findif_step_size, and findif_verbose.
        * TODO hangers-on keys present at class initiation get automatically attached to class since `extra = "allow"` but should be pruned

        """
        findif_stencil_size = data.pop('findif_stencil_size')
        findif_step_size = data.pop('findif_step_size')

        BaseComputer.__init__(self, **data)

        translations_projection_sound = (not "external_potentials" in data['keywords']['function_kwargs']
                                         and not core.get_option('SCF', 'PERTURB_H')
                                         and not hasattr(self.molecule, 'EFP'))
        if 'ref_gradient' in data:
            logger.info("""hessian() using ref_gradient to assess stationary point.""")
            stationary_criterion = 1.e-2  # pulled out of a hat
            stationary_point = _rms(data['ref_gradient']) < stationary_criterion
        else:
            stationary_point = False  # unknown, so F to be safe
        rotations_projection_sound_grad = translations_projection_sound
        rotations_projection_sound_hess = translations_projection_sound and stationary_point
        if core.has_option_changed('FINDIF', 'FD_PROJECT'):
            r_project_grad = core.get_option('FINDIF', 'FD_PROJECT')
            r_project_hess = core.get_option('FINDIF', 'FD_PROJECT')
        else:
            r_project_grad = rotations_projection_sound_grad
            r_project_hess = rotations_projection_sound_hess

        for kwg in ['dft_functional']:
            if kwg in data:
                data['keywords']['function_kwargs'][kwg] = data.pop(kwg)
        # I have the feeling the keywords.function_kwargs should be all left over in data
        #   after the findif control ones are removed, not this by-name procedure
        data['keywords']['PARENT_SYMMETRY'] = self.molecule.point_group().full_name()

        self.method = data['method']

        self.metameta['mode'] = str(data['findif_mode'][0]) + '_' + str(data['findif_mode'][1])
        self.metameta['irrep'] = data.pop('findif_irrep', -1)

        if self.metameta['mode'] == '1_0':
            self.metameta['proxy_driver'] = 'energy'
            self.findifrec = gradient_from_energies_geometries(self.molecule,
                                                               stencil_size=findif_stencil_size,
                                                               step_size=findif_step_size,
                                                               t_project=translations_projection_sound,
                                                               r_project=r_project_grad)

        elif self.metameta['mode'] == '2_1':
            self.metameta['proxy_driver'] = 'gradient'
            self.findifrec = hessian_from_gradients_geometries(self.molecule,
                                                               freq_irrep_only=self.metameta['irrep'],
                                                               stencil_size=findif_stencil_size,
                                                               step_size=findif_step_size,
                                                               t_project=translations_projection_sound,
                                                               r_project=r_project_hess)

        elif self.metameta['mode'] == '2_0':
            self.metameta['proxy_driver'] = 'energy'
            self.findifrec = hessian_from_energies_geometries(self.molecule,
                                                              freq_irrep_only=self.metameta['irrep'],
                                                              stencil_size=findif_stencil_size,
                                                              step_size=findif_step_size,
                                                              t_project=translations_projection_sound,
                                                              r_project=r_project_hess)

        ndisp = len(self.findifrec["displacements"]) + 1
        info = f""" {ndisp} displacements needed ...\n"""
        core.print_out(info)
        logger.debug(info)

        # var_dict = core.variables()
        packet = {
            "molecule": self.molecule,
            "driver": self.metameta['proxy_driver'],
            "method": self.method,
            "basis": data["basis"],
            "keywords": data["keywords"] or {},
        }
        if 'cbs_metadata' in data:
            packet['cbs_metadata'] = data['cbs_metadata']
        passalong = {k: v for k, v in data.items() if k not in packet}
        passalong.pop('ptype', None)

        self.task_list["reference"] = self.computer(**packet, **passalong)

        parent_group = self.molecule.point_group()
        for label, displacement in self.findifrec["displacements"].items():
            clone = self.molecule.clone()
            clone.reinterpret_coordentry(False)
            #clone.fix_orientation(True)

            # Load in displacement into the active molecule
            clone.set_geometry(core.Matrix.from_array(displacement["geometry"]))

            # If the user insists on symmetry, weaken it if some is lost when displacing.
            # or 'fix_symmetry' in self.findifrec.molecule
            logger.debug(f'SYMM {clone.schoenflies_symbol()}')
            if self.molecule.symmetry_from_input():
                disp_group = clone.find_highest_point_group()
                new_bits = parent_group.bits() & disp_group.bits()
                new_symm_string = qcdb.PointGroup.bits_to_full_name(new_bits)
                clone.reset_point_group(new_symm_string)

            packet = {
                "molecule": clone,
                "driver": self.metameta['proxy_driver'],
                "method": self.method,
                "basis": data["basis"],
                "keywords": data["keywords"] or {},
            }
            # Displacements can run in lower symmetry. Don't overwrite orbitals from reference geom
            packet['keywords']['function_kwargs'].update({"write_orbitals": False})
            if 'cbs_metadata' in data:
                packet['cbs_metadata'] = data['cbs_metadata']

            self.task_list[label] = self.computer(**packet, **passalong)


#        for n, displacement in enumerate(findif_meta_dict["displacements"].values(), start=2):
#            _process_displacement(energy, lowername, molecule, displacement, n, ndisp, write_orbitals=False, **kwargs)

    def build_tasks(self, obj, **kwargs):
        # permanently a dummy function
        pass

    def plan(self):
        # uncalled function
        return [t.plan() for t in self.task_list.values()]

    def compute(self, client: Optional["qcportal.FractalClient"] = None):
        """Run each job in task list."""
        instructions = "\n" + p4util.banner(f" FiniteDifference Computations", strNotOutfile=True) + "\n"
        logger.debug(instructions)
        core.print_out(instructions)

        with p4util.hold_options_state():
            for t in self.task_list.values():
                t.compute(client=client)

    def _prepare_results(self, client: Optional["qcportal.FractalClient"] = None):
        results_list = {k: v.get_results(client=client) for k, v in self.task_list.items()}

        # load AtomicComputer results into findifrec[reference]
        reference = self.findifrec["reference"]
        task = results_list["reference"]
        response = task.return_result
        reference["module"] = getattr(task.provenance, "module", None)

        if task.driver == 'energy':
            reference['energy'] = response

        elif task.driver == 'gradient':
            reference['gradient'] = response
            reference['energy'] = task.extras['qcvars']['CURRENT ENERGY']

        elif task.driver == 'hessian':
            reference['hessian'] = response
            reference['energy'] = task.extras['qcvars']['CURRENT ENERGY']
            if 'CURRENT GRADIENT' in task.extras['qcvars']:
                reference['gradient'] = task.extras['qcvars']['CURRENT GRADIENT']

        dipole_available = False
        if 'CURRENT DIPOLE' in task.extras['qcvars']:
            reference['dipole'] = task.extras['qcvars']['CURRENT DIPOLE']
            dipole_available = True

        # load AtomicComputer results into findifrec[displacements]
        for label, displacement in self.findifrec["displacements"].items():
            task = results_list[label]
            response = task.return_result

            if task.driver == 'energy':
                displacement['energy'] = response

            elif task.driver == 'gradient':
                displacement['gradient'] = response
                displacement['energy'] = task.extras['qcvars']['CURRENT ENERGY']

            elif task.driver == 'hessian':
                displacement['hessian'] = response
                displacement['energy'] = task.extras['qcvars']['CURRENT ENERGY']
                if 'CURRENT GRADIENT' in task.extras['qcvars']:
                    displacement['gradient'] = task.extras['qcvars']['CURRENT GRADIENT']

            if 'CURRENT DIPOLE' in task.extras['qcvars']:
                displacement['dipole'] = task.extras['qcvars']['CURRENT DIPOLE']

        # apply finite difference formulas and load derivatives into findifrec[reference]
        if self.metameta['mode'] == '1_0':
            G0 = assemble_gradient_from_energies(self.findifrec)
            self.findifrec["reference"][self.driver.name] = G0

        elif self.metameta['mode'] == '2_1':
            if dipole_available:
                DD0 = assemble_dipder_from_dipoles(self.findifrec, self.metameta['irrep'])
                self.findifrec["reference"]["dipole derivative"] = DD0

            H0 = assemble_hessian_from_gradients(self.findifrec, self.metameta['irrep'])
            self.findifrec["reference"][self.driver.name] = H0

        elif self.metameta['mode'] == '2_0':
            try:
                G0 = assemble_gradient_from_energies(self.findifrec)
            except KeyError:
                core.print_out("Unable to construct reference gradient from Hessian displacements.")
                # TODO: this happens properly when the requested symmetry block
                #  of displacements don't have the totally symmetric displacements
                #  needed for gradient. For both this case
                #  and distributed computing are-we-there-yet? queries,
                #  should have a probe as to whether all the
                #  findif[displacement] labels are present and whether
                #  all the findif[displacement][energy-or-gradient] values
                #  are ready. Not sure what type of error/query is best,
                #  so deferring for now. Also, possibly need to check if
                #  step size matches before using values from one findifrec
                #  to construct another quantity.
            else:
                self.findifrec["reference"]["gradient"] = G0

            if dipole_available:
                DD0 = assemble_dipder_from_dipoles(self.findifrec, self.metameta['irrep'])
                self.findifrec["reference"]["dipole derivative"] = DD0

            H0 = assemble_hessian_from_energies(self.findifrec, self.metameta['irrep'])
            self.findifrec["reference"][self.driver.name] = H0

    def get_results(self, client: Optional["qcportal.FractalClient"] = None) -> AtomicResult:
        """Return results as FiniteDifference-flavored QCSchema."""

        instructions = "\n" + p4util.banner(f" FiniteDifference Results", strNotOutfile=True) + "\n"
        core.print_out(instructions)

        self._prepare_results(client=client)  # assembled_results

        # load QCVariables & properties
        qcvars = self.task_list['reference'].get_results().extras['qcvars']
        E0 = self.findifrec['reference']['energy']

        properties = {
            "calcinfo_natom": self.molecule.natom(),
            "nuclear_repulsion_energy": self.molecule.nuclear_repulsion_energy(),
            "return_energy": E0,
        }

        qcvars['FINDIF NUMBER'] = len(self.task_list)
        qcvars['NUCLEAR REPULSION ENERGY'] = self.molecule.nuclear_repulsion_energy()
        qcvars['CURRENT ENERGY'] = E0

        DD0 = self.findifrec['reference'].get('dipole derivative')
        if DD0 is not None:
            qcvars['CURRENT DIPOLE GRADIENT'] = DD0
            qcvars[f"{self.method.upper()} DIPOLE GRADIENT"] = DD0

        G0 = self.findifrec['reference'].get('gradient')
        if G0 is not None:
            qcvars['CURRENT GRADIENT'] = G0
            qcvars[f"{self.method.upper()} TOTAL GRADIENT"] = G0
            properties["return_gradient"] = G0

        H0 = self.findifrec['reference'].get('hessian')
        if H0 is not None:
            qcvars['CURRENT HESSIAN'] = H0
            qcvars[f"{self.method.upper()} TOTAL HESSIAN"] = H0
            properties["return_hessian"] = H0

#        if isinstance(lowername, str) and lowername in procedures['energy']:
#            # this correctly filters out cbs fn and "hf/cc-pvtz"
#            # it probably incorrectly filters out mp5, but reconsider in DDD

        findif_model = AtomicResult(
            **{
                'driver': self.driver,
                'model': {
                    "basis": self.basis,
                    'method': self.method,
                },
                'molecule': self.molecule.to_schema(dtype=2),
                'properties': properties,
                'provenance': p4util.provenance_stamp(__name__, module=self.findifrec["reference"]["module"]),
                'extras': {
                    'qcvars': qcvars,
                    'findif_record': copy.deepcopy(self.findifrec),
                },
                'return_result': self.findifrec['reference'][self.driver.name],
                'success': True,
            })

        logger.debug('\nFINDIF QCSchema:\n' + pp.pformat(findif_model))

        return findif_model

    def get_psi_results(self, return_wfn: bool = False) -> EnergyGradientHessianWfnReturn:
        """Called by driver to assemble results into FiniteDifference-flavored QCSchema,
        then reshape and return them in the customary Psi4 driver interface: ``(e/g/h, wfn)``.

        Parameters
        ----------
        return_wfn
            Whether to additionally return the dummy :py:class:`~psi4.core.Wavefunction`
            calculation result as the second element of a tuple. Contents are:

            - undisplaced molecule
            - compute basis if simple, else dummy basis def2-svp
            - e/g/h member data
            - QCVariables
            - module

        Returns
        -------
        ret
            Gradient or Hessian according to self.driver.
        wfn
            Wavefunction described above when *return_wfn* specified.

        """
        findif_model = self.get_results()

        ret_ptype = core.Matrix.from_array(findif_model.return_result)
        wfn = _findif_schema_to_wfn(findif_model)

        gradient_write(wfn)
        hessian_write(wfn)

        if return_wfn:
            return (ret_ptype, wfn)
        else:
            return ret_ptype


def _findif_schema_to_wfn(findif_model: AtomicResult) -> core.Wavefunction:
    """Helper function to produce Wavefunction and Psi4 files from a FiniteDifference-flavored AtomicResult."""

    # new skeleton wavefunction w/mol, highest-SCF basis (just to choose one), & not energy
    mol = core.Molecule.from_schema(findif_model.molecule.dict(), nonphysical=True)
    sbasis = "def2-svp" if (findif_model.model.basis == "(auto)") else findif_model.model.basis
    basis = core.BasisSet.build(mol, "ORBITAL", sbasis, quiet=True)
    wfn = core.Wavefunction(mol, basis)
    if hasattr(findif_model.provenance, "module"):
        wfn.set_module(findif_model.provenance.module)

    # setting CURRENT E/G/H on wfn below catches Wfn.energy_, gradient_, hessian_
    # setting CURRENT E/G/H on core below is authoritative P::e record
    for qcv, val in findif_model.extras["qcvars"].items():
        for obj in [core, wfn]:
            obj.set_variable(qcv, val)

    return wfn


def hessian_write(wfn: core.Wavefunction):
    if core.get_option('FINDIF', 'HESSIAN_WRITE'):
        filename = core.get_writer_file_prefix(wfn.molecule().name()) + ".hess"
        with open(filename, 'wb') as handle:
            qcdb.hessparse.to_string(np.asarray(wfn.hessian()), handle, dtype='psi4')


def gradient_write(wfn: core.Wavefunction):
    if core.get_option('FINDIF', 'GRADIENT_WRITE'):
        filename = core.get_writer_file_prefix(wfn.molecule().name()) + ".grad"
        qcdb.gradparse.to_string(np.asarray(wfn.gradient()),
                                 filename,
                                 dtype='GRD',
                                 mol=wfn.molecule(),
                                 energy=wfn.energy())


def _rms(arr: Union[core.Matrix, np.ndarray]) -> float:
    """Compute root-mean-square of array, be it Psi4 or NumPy array."""
    if isinstance(arr, np.ndarray):
        return np.sqrt(np.mean(np.square(arr)))
    else:
        return arr.rms()
