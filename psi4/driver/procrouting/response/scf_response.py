#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2021 The Psi4 Developers.
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

from typing import Union, List
try:
    from dataclasses import dataclass
except ImportError:
    from pydantic.dataclasses import dataclass

import numpy as np

from psi4 import core
from psi4.driver import constants
from psi4.driver.p4util import solvers
from psi4.driver.p4util.exceptions import *
from psi4.driver.procrouting.response.scf_products import (TDRSCFEngine, TDUSCFEngine)

dipole = {
    'name': 'Dipole polarizabilities',
    'printout_labels': ['X', 'Y', 'Z'],
    'mints_function': core.MintsHelper.ao_dipole,
    'vector names': ['AO Mux', 'AO Muy', 'AO Muz']
}

quadrupole = {
    'name': 'Quadrupole polarizabilities',
    'printout_labels': ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'],
    'mints_function': core.MintsHelper.ao_quadrupole,
}
quadrupole['vector names'] = ["AO Quadrupole " + x for x in quadrupole["printout_labels"]]

traceless_quadrupole = {
    'name': 'Traceless quadrupole polarizabilities',
    'printout_labels': ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'],
    'mints_function': core.MintsHelper.ao_traceless_quadrupole,
}
traceless_quadrupole['vector names'] = [
    "AO Traceless Quadrupole " + x for x in traceless_quadrupole["printout_labels"]
]

property_dicts = {
    'DIPOLE_POLARIZABILITIES': dipole,
    'QUADRUPOLE_POLARIZABILITIES': quadrupole,
    'TRACELESS_QUADRUPOLE_POLARIZABILITIES': traceless_quadrupole
}


def cpscf_linear_response(wfn, *args, **kwargs):
    """
    Compute the static properties from a reference wavefunction. The currently implemented properties are
      - dipole polarizability
      - quadrupole polarizability

    Parameters
    ----------
    wfn : psi4 wavefunction
        The reference wavefunction.
    args : list
        The list of arguments. For each argument, such as ``dipole polarizability``, will return the corresponding
        response. The user may also choose to pass a list or tuple of custom vectors.
    kwargs : dict
        Options that control how the response is computed. The following options are supported (with default values):
          - ``conv_tol``: 1e-5
          - ``max_iter``: 10
          - ``print_lvl``: 2

    Returns
    -------
    responses : list
        The list of responses.
    """
    mints = core.MintsHelper(wfn.basisset())

    # list of dictionaries to control response calculations, count how many user-supplied vectors we have
    complete_dict = []
    n_user = 0

    for arg in args:

        # for each string keyword, append the appropriate dictionary (vide supra) to our list
        if isinstance(arg, str):
            ret = property_dicts.get(arg)
            if ret:
                complete_dict.append(ret)
            else:
                raise ValidationError('Do not understand {}. Abort.'.format(arg))

        # the user passed a list of vectors. absorb them into a dictionary
        elif isinstance(arg, tuple) or isinstance(arg, list):
            complete_dict.append({
                'name': 'User Vectors',
                'length': len(arg),
                'vectors': arg,
                'vector names': ['User Vector {}_{}'.format(n_user, i) for i in range(len(arg))]
            })
            n_user += len(arg)

        # single vector passed. stored in a dictionary as a list of length 1 (can be handled as the case above that way)
        # note: the length is set to '0' to designate that it was not really passed as a list
        else:
            complete_dict.append({
                'name': 'User Vector',
                'length': 0,
                'vectors': [arg],
                'vector names': ['User Vector {}'.format(n_user)]
            })
            n_user += 1

    # vectors will be passed to the cphf solver, vector_names stores the corresponding names
    vectors = []
    vector_names = []

    # construct the list of vectors. for the keywords, fetch the appropriate tensors from MintsHelper
    for prop in complete_dict:
        if 'User' in prop['name']:
            for name, vec in zip(prop['vector names'], prop['vectors']):
                vectors.append(vec)
                vector_names.append(name)

        else:
            tmp_vectors = prop['mints_function'](mints)
            for tmp in tmp_vectors:
                tmp.scale(-2.0)  # RHF only
                vectors.append(tmp)
                vector_names.append(tmp.name)

    # do we have any vectors to work with?
    if len(vectors) == 0:
        raise ValidationError('I have no vectors to work with. Aborting.')

    # print information on module, vectors that will be used
    _print_header(complete_dict, n_user)

    # fetch wavefunction information
    nmo = wfn.nmo()
    ndocc = wfn.nalpha()
    nvirt = nmo - ndocc

    c_occ = wfn.Ca_subset("AO", "OCC")
    c_vir = wfn.Ca_subset("AO", "VIR")
    nbf = c_occ.shape[0]

    # the vectors need to be in the MO basis. if they have the shape nbf x nbf, transform.
    for i in range(len(vectors)):
        shape = vectors[i].shape

        if shape == (nbf, nbf):
            vectors[i] = core.triplet(c_occ, vectors[i], c_vir, True, False, False)

        # verify that this vector already has the correct shape
        elif shape != (ndocc, nvirt):
            raise ValidationError('ERROR: "{}" has an unrecognized shape ({}, {}). Must be either ({}, {}) or ({}, {})'.format(
                vector_names[i], shape[0], shape[1], nbf, nbf, ndocc, nvirt))

    # compute response vectors for each input vector
    params = [kwargs.pop("conv_tol", 1.e-5), kwargs.pop("max_iter", 10), kwargs.pop("print_lvl", 2)]

    responses = wfn.cphf_solve(vectors, *params)

    # zip vectors, responses for easy access
    vectors = {k: v for k, v in zip(vector_names, vectors)}
    responses = {k: v for k, v in zip(vector_names, responses)}

    # compute response values, format output
    output = []
    for prop in complete_dict:

        # try to replicate the data structure of the input
        if 'User' in prop['name']:
            if prop['length'] == 0:
                output.append(responses[prop['vector names'][0]])
            else:
                buf = []
                for name in prop['vector names']:
                    buf.append(responses[name])
                output.append(buf)

        else:
            names = prop['vector names']
            dim = len(names)

            buf = np.zeros((dim, dim))

            for i, i_name in enumerate(names):
                for j, j_name in enumerate(names):
                    buf[i, j] = -1.0 * vectors[i_name].vector_dot(responses[j_name])

            output.append(buf)

    _print_output(complete_dict, output)

    return output


def _print_header(complete_dict, n_user):
    core.print_out('\n\n         ---------------------------------------------------------\n'
                   '         {:^57}\n'.format('CPSCF Linear Response Solver') +
                   '         {:^57}\n'.format('by Marvin Lechner and Daniel G. A. Smith') +
                   '         ---------------------------------------------------------\n')

    core.print_out('\n   ==> Requested Responses <==\n\n')

    for prop in complete_dict:
        if 'User' not in prop['name']:
            core.print_out('    {}\n'.format(prop['name']))

    if n_user != 0:
        core.print_out('    {} user-supplied vector(s)\n'.format(n_user))


def _print_matrix(descriptors, content, title):
    length = len(descriptors)

    matrix_header = '         ' + ' {:^10}' * length + '\n'
    core.print_out(matrix_header.format(*descriptors))
    core.print_out('    -----' + ' ----------' * length + '\n')

    for i, desc in enumerate(descriptors):
        core.print_out('    {:^5}'.format(desc))
        for j in range(length):
            core.print_out(' {:>10.5f}'.format(content[i, j]))

            # Set the name
            var_name = title + " " + descriptors[i] + descriptors[j]
            core.set_variable(var_name, content[i, j])
        core.print_out('\n')


def _print_output(complete_dict, output):
    core.print_out('\n   ==> Response Properties <==\n')

    for i, prop in enumerate(complete_dict):
        if not 'User' in prop['name']:
            core.print_out('\n    => {} <=\n\n'.format(prop['name']))
            directions = prop['printout_labels']
            var_name = prop['name'].upper().replace("IES", "Y")
            _print_matrix(directions, output[i], var_name)


def _print_tdscf_header(*, r_convergence: float, guess_type: str, restricted: bool, ptype: str):
    core.print_out("\n\n         ---------------------------------------------------------\n"
                   f"         {'TDSCF excitation energies':^57}\n" +
                   f"         {'by Andrew M. James and Daniel G. A. Smith':^57}\n" +
                   "         ---------------------------------------------------------\n")

    core.print_out("\n  ==> Options <==\n\n")
    core.print_out(f"     {'Residual threshold':<20s}: {r_convergence:.4e}\n")
    core.print_out(f"     {'Initial guess':20s}: {guess_type.lower()}\n")
    reference = 'RHF' if restricted else 'UHF'
    core.print_out(f"     {'Reference':20s}: {reference}\n")
    solver_type = 'Hamiltonian' if ptype == "RPA" else "Davidson"
    core.print_out(f"     {'Solver type':20s}: {ptype} ({solver_type})\n")

    core.print_out("\n")


@dataclass
class _TDSCFResults:
    E_ex_au: float
    irrep_GS: str
    irrep_ES: str
    irrep_trans: str
    edtm_length: np.ndarray
    f_length: float
    edtm_velocity: np.ndarray
    f_velocity: float
    mdtm: np.ndarray
    R_length: float
    R_velocity: float
    spin_mult: str
    R_eigvec: Union[core.Matrix, List[core.Matrix]]
    L_eigvec: Union[core.Matrix, List[core.Matrix]]


def _solve_loop(wfn,
                ptype,
                solve_function,
                states_per_irrep: List[int],
                maxiter: int,
                restricted: bool = True,
                spin_mult: str = "singlet") -> List[_TDSCFResults]:
    """

    References
    ----------
    For the expression of the transition moments in length and velocity gauges:

    - T. B. Pedersen, A. E. Hansen, "Ab Initio Calculation and Display of the
    Rotary Strength Tensor in the Random Phase Approximation. Method and Model
    Studies." Chem. Phys. Lett., 246, 1 (1995)
    - P. J. Lestrange, F. Egidi, X. Li, "The Consequences of Improperly
    Describing Oscillator Strengths beyond the Electric Dipole Approximation."
    J. Chem. Phys., 143, 234103 (2015)
    """

    core.print_out("\n  ==> Requested Excitations <==\n\n")
    for nstate, state_sym in zip(states_per_irrep, wfn.molecule().irrep_labels()):
        core.print_out(f"      {nstate} {spin_mult} states with {state_sym} symmetry\n")

    # construct the engine
    if restricted:
        if spin_mult == "triplet":
            engine = TDRSCFEngine(wfn, ptype=ptype.lower(), triplet=True)
        else:
            engine = TDRSCFEngine(wfn, ptype=ptype.lower(), triplet=False)
    else:
        engine = TDUSCFEngine(wfn, ptype=ptype.lower())

    # collect results and compute some spectroscopic observables
    mints = core.MintsHelper(wfn.basisset())
    results = []
    irrep_GS = wfn.molecule().irrep_labels()[engine.G_gs]
    for state_sym, nstates in enumerate(states_per_irrep):
        if nstates == 0:
            continue
        irrep_ES = wfn.molecule().irrep_labels()[state_sym]
        core.print_out(f"\n\n  ==> Seeking the lowest {nstates} {spin_mult} states with {irrep_ES} symmetry")
        engine.reset_for_state_symm(state_sym)
        guess_ = engine.generate_guess(nstates * 4)

        # ret = {"eigvals": ee, "eigvecs": (rvecs, rvecs), "stats": stats} (TDA)
        # ret = {"eigvals": ee, "eigvecs": (rvecs, lvecs), "stats": stats} (RPA)
        ret = solve_function(engine, nstates, guess_, maxiter)

        # check whether all roots converged
        if not ret["stats"][-1]["done"]:
            # raise error
            raise TDSCFConvergenceError(maxiter, wfn, f"singlet excitations in irrep {irrep_ES}", ret["stats"][-1])

        # flatten dictionary: helps with sorting by energy
        # also append state symmetry to return value
        for e, (R, L) in zip(ret["eigvals"], ret["eigvecs"]):
            irrep_trans = wfn.molecule().irrep_labels()[engine.G_gs ^ state_sym]

            # length-gauge electric dipole transition moment
            edtm_length = engine.residue(R, mints.so_dipole())
            # length-gauge oscillator strength
            f_length = ((2 * e) / 3) * np.sum(edtm_length**2)
            # velocity-gauge electric dipole transition moment
            edtm_velocity = engine.residue(L, mints.so_nabla())
            ## velocity-gauge oscillator strength
            f_velocity = (2 / (3 * e)) * np.sum(edtm_velocity**2)
            # length gauge magnetic dipole transition moment
            # 1/2 is the Bohr magneton in atomic units
            mdtm = 0.5 * engine.residue(L, mints.so_angular_momentum())
            # NOTE The signs for rotatory strengths are opposite WRT the cited paper.
            # This is becasue Psi4 defines length-gauge dipole integral to include the electron charge (-1.0)
            # length gauge rotatory strength
            R_length = np.einsum("i,i", edtm_length, mdtm)
            # velocity gauge rotatory strength
            R_velocity = -np.einsum("i,i", edtm_velocity, mdtm) / e

            results.append(
                _TDSCFResults(e, irrep_GS, irrep_ES, irrep_trans, edtm_length, f_length, edtm_velocity, f_velocity,
                              mdtm, R_length, R_velocity, spin_mult, R, L))

    return results


def _states_per_irrep(states, nirrep):
    """Distributes states into nirrep"""
    spi = [states // nirrep] * nirrep
    for i in range(states % nirrep):
        spi[i] += 1

    return spi


def _validate_tdscf(*, wfn, states, triplets, guess) -> None:

    # validate states
    if not isinstance(states, (int, list)):
        raise ValidationError("TDSCF: Number of states must be either an integer or a list of integers")

    # list of states per irrep given, validate it
    if isinstance(states, list):
        if len(states) != wfn.nirrep():
            raise ValidationError(f"TDSCF: States requested ({states}) do not match number of irreps ({wfn.nirrep()})")

    # do triplets?
    if triplets not in ["NONE", "ALSO", "ONLY"]:
        raise ValidationError(
            f"TDSCF: Triplet option ({triplets}) unrecognized. Must be one of 'NONE', 'ALSO' or 'ONLY'")

    restricted = wfn.same_a_b_orbs()
    do_triplets = False if triplets == "NONE" else True
    if (not restricted) and do_triplets:
        raise ValidationError("TDSCF: Cannot compute triplets with an unrestricted reference")

    # determine how many states per irrep to seek and apportion them between singlets/triplets and irreps.

    # validate calculation
    if restricted and wfn.functional().needs_xc() and do_triplets:
        raise ValidationError("TDSCF: Restricted Vx kernel only spin-adapted for singlets")

    not_lda = wfn.functional().is_gga() or wfn.functional().is_meta()
    if (not restricted) and not_lda:
        raise ValidationError("TDSCF: Unrestricted Kohn-Sham Vx kernel currently limited to SVWN functional")

    if guess != "DENOMINATORS":
        raise ValidationError(f"TDSCF: Guess type {guess} is not valid")


def tdscf_excitations(wfn,
                      *,
                      states: Union[int, List[int]],
                      triplets: str = "NONE",
                      tda: bool = False,
                      r_convergence: float = 1.0e-4,
                      maxiter: int = 60,
                      guess: str = "DENOMINATORS",
                      verbose: int = 1):
    """Compute excitations from a SCF(HF/KS) wavefunction

    Parameters
    -----------
    wfn : :py:class:`psi4.core.Wavefunction`
       The reference wavefunction
    states : Union[int, List[int]]
       How many roots (excited states) should the solver seek to converge?
       This function accepts either an integer or a list of integers:
         - The list has :math:`n_{\mathrm{irrep}}` elements and is only
           acceptable if the system has symmetry. It tells the solver how many
           states per irrep to calculate.
         - If an integer is given _and_ the system has symmetry, the states
           will be distributed among irreps.
           For example, ``states = 10`` for a D2h system will compute 10 states
           distributed as ``[2, 2, 1, 1, 1, 1, 1, 1]`` among irreps.
    triplets : {"NONE", "ONLY", "ALSO"}
       Should the solver seek to converge states of triplet symmetry?
       Default is `none`: do not seek to converge triplets.
       Valid options are:
         - `NONE`. Do not seek to converge triplets.
         - `ONLY`. Only seek to converge triplets.
         - `ALSO`. Seek to converge both triplets and singlets. This choice is
           only valid for restricted reference wavefunction.
           The number of states given will be apportioned roughly 50-50 between
           singlet and triplet states, preferring the former. For example:
           given ``state = 5, triplets = "ALSO"``, the solver will seek to
           converge 3 states of singlet spin symmetry and 2 of triplet spin
           symmetry. When asking for ``states = [3, 3, 3, 3], triplets =
           "ALSO"`` states (C2v symmetry), ``[2, 2, 2, 2]`` will be of singlet
           spin symmetry and ``[1, 1, 1, 1]``` will be of triplet spin
           symmetry.
    tda :  bool, optional.
       Should the solver use the Tamm-Dancoff approximation (TDA) or the
       random-phase approximation (RPA)?
       Default is ``False``: use RPA.
       Note that TDA is equivalent to CIS for HF references.
    r_convergence : float, optional.
       The convergence threshold for the norm of the residual vector.
       Default: 1.0e-4
       Using a tighter convergence threshold here requires tighter SCF ground
       state convergence threshold. As a rule of thumb, with the SCF ground
       state density converged to :math:`10^{-N}` (``D_CONVERGENGE = 1.0e-N``),
       you can afford converging a corresponding TDSCF calculation to
       :math:`10^{-(N-2)}`.
       The default value is consistent with the default value for
       ``D_CONVERGENCE``.
    maxiter : int, optional
       Maximum number of iterations.
       Default: 60
    guess : str, optional.
       How should the starting trial vectors be generated?
       Default: `DENOMINATORS`, i.e. use orbital energy differences to generate
       guess vectors.
    verbose : int, optional.
       How verbose should the solver be?
       Default: 1

    Notes
    -----
    The algorithm employed to solve the non-Hermitian eigenvalue problem (``tda = False``)
    will fail when the SCF wavefunction has a triplet instability.

    This function can be used for:
      - restricted singlets: RPA or TDA, any functional
      - restricted triplets: RPA or TDA, Hartree-Fock only
      - unresctricted: RPA or TDA, Hartre-Fock and LDA only

    Tighter convergence thresholds will require a larger iterative subspace.
    The maximum size of the iterative subspace is calculated based on `r_convergence`:

       max_vecs_per_root = -np.log10(r_convergence) * 50

    for the default converegence threshold this gives 200 trial vectors per root and a maximum subspace size
    of:

       max_ss_size = max_vecs_per_root * n

    where `n` are the number of roots to seek in the given irrep.
    For each irrep, the algorithm will store up to `max_ss_size` trial vectors
    before collapsing (restarting) the iterations from the `n` best
    approximations.
    """

    # validate input parameters
    triplets = triplets.upper()
    guess = guess.upper()
    _validate_tdscf(wfn=wfn, states=states, triplets=triplets, guess=guess)

    restricted = wfn.same_a_b_orbs()

    # determine how many states per irrep to seek and apportion them between singlets/triplets and irreps.
    singlets_per_irrep = []
    triplets_per_irrep = []
    if isinstance(states, list):
        if triplets == "ONLY":
            triplets_per_irrep = states
        elif triplets == "ALSO":
            singlets_per_irrep = [(s // 2) + (s % 2) for s in states]
            triplets_per_irrep = [(s // 2) for s in states]
        else:
            singlets_per_irrep = states
    else:
        # total number of states given
        # first distribute them among singlets and triplets, preferring the
        # former then distribute them among irreps
        if triplets == "ONLY":
            triplets_per_irrep = _states_per_irrep(states, wfn.nirrep())
        elif triplets == "ALSO":
            spi = (states // 2) + (states % 2)
            singlets_per_irrep = _states_per_irrep(spi, wfn.nirrep())
            tpi = states - spi
            triplets_per_irrep = _states_per_irrep(tpi, wfn.nirrep())
        else:
            singlets_per_irrep = _states_per_irrep(states, wfn.nirrep())

    # tie maximum number of vectors per root to requested residual tolerance
    # This gives 200 vectors per root with default tolerance
    max_vecs_per_root = int(-np.log10(r_convergence) * 50)

    def rpa_solver(e, n, g, m):
        return solvers.hamiltonian_solver(engine=e,
                                          nroot=n,
                                          guess=g,
                                          r_convergence=r_convergence,
                                          max_ss_size=max_vecs_per_root * n,
                                          verbose=verbose)

    def tda_solver(e, n, g, m):
        return solvers.davidson_solver(engine=e,
                                       nroot=n,
                                       guess=g,
                                       r_convergence=r_convergence,
                                       max_ss_size=max_vecs_per_root * n,
                                       verbose=verbose)

    # determine which solver function to use: Davidson for TDA or Hamiltonian for RPA?
    if tda:
        ptype = "TDA"
        solve_function = tda_solver
    else:
        ptype = "RPA"
        solve_function = rpa_solver

    _print_tdscf_header(r_convergence=r_convergence, guess_type=guess, restricted=restricted, ptype=ptype)

    # collect solver results into a list
    _results = []

    # singlets solve loop
    if triplets == "NONE" or triplets == "ALSO":
        res_1 = _solve_loop(wfn, ptype, solve_function, singlets_per_irrep, maxiter, restricted, "singlet")
        _results.extend(res_1)

    # triplets solve loop
    if triplets == "ALSO" or triplets == "ONLY":
        res_3 = _solve_loop(wfn, ptype, solve_function, triplets_per_irrep, maxiter, restricted, "triplet")
        _results.extend(res_3)

    # sort by energy
    _results = sorted(_results, key=lambda x: x.E_ex_au)

    core.print_out("\n{}\n".format("*"*90) +
                   "{}{:^70}{}\n".format("*"*10, "WARNING", "*"*10) +
                   "{}{:^70}{}\n".format("*"*10, "Length-gauge rotatory strengths are **NOT** gauge-origin invariant", "*"*10) +
                   "{}\n\n".format("*"*90)) #yapf: disable

    # print results
    core.print_out("        " + (" " * 20) + " " + "Excitation Energy".center(31) + f" {'Total Energy':^15}" +
                   "Oscillator Strength".center(31) + "Rotatory Strength".center(31) + "\n")
    core.print_out(
        f"    {'#':^4} {'Sym: GS->ES (Trans)':^20} {'au':^15} {'eV':^15} {'au':^15} {'au (length)':^15} {'au (velocity)':^15} {'au (length)':^15} {'au (velocity)':^15}\n"
    )
    core.print_out(
        f"    {'-':->4} {'-':->20} {'-':->15} {'-':->15} {'-':->15} {'-':->15} {'-':->15} {'-':->15} {'-':->15}\n")

    # collect results
    solver_results = []
    for i, x in enumerate(_results):
        sym_descr = f"{x.irrep_GS}->{x.irrep_ES} ({1 if x.spin_mult== 'singlet' else 3} {x.irrep_trans})"

        E_ex_ev = constants.conversion_factor('hartree', 'eV') * x.E_ex_au

        E_tot_au = wfn.energy() + x.E_ex_au

        # prepare return dictionary for this root
        solver_results.append({
            "EXCITATION ENERGY": x.E_ex_au,
            "ELECTRIC DIPOLE TRANSITION MOMENT (LEN)": x.edtm_length,
            "OSCILLATOR STRENGTH (LEN)": x.f_length,
            "ELECTRIC DIPOLE TRANSITION MOMENT (VEL)": x.edtm_velocity,
            "OSCILLATOR STRENGTH (VEL)": x.f_velocity,
            "MAGNETIC DIPOLE TRANSITION MOMENT": x.mdtm,
            "ROTATORY STRENGTH (LEN)": x.R_length,
            "ROTATORY STRENGTH (VEL)": x.R_velocity,
            "SYMMETRY": x.irrep_trans,
            "SPIN": x.spin_mult,
            "RIGHT EIGENVECTOR ALPHA": x.R_eigvec if restricted else x.R_eigvec[0],
            "LEFT EIGENVECTOR ALPHA": x.L_eigvec if restricted else x.L_eigvec[0],
            "RIGHT EIGENVECTOR BETA": x.R_eigvec if restricted else x.R_eigvec[1],
            "LEFT EIGENVECTOR BETA": x.L_eigvec if restricted else x.L_eigvec[1],
        })

        # stash in psivars/wfnvars
        ssuper_name = wfn.functional().name()
        # wfn.set_variable("TD-fctl ROOT n TOTAL ENERGY - h SYMMETRY")  # P::e SCF
        # wfn.set_variable("TD-fctl ROOT 0 -> ROOT m EXCITATION ENERGY - h SYMMETRY")  # P::e SCF
        # wfn.set_variable("TD-fctl ROOT 0 -> ROOT m OSCILLATOR STRENGTH (LEN) - h SYMMETRY")  # P::e SCF
        # wfn.set_variable("TD-fctl ROOT 0 -> ROOT m OSCILLATOR STRENGTH (VEL) - h SYMMETRY")  # P::e SCF
        # wfn.set_variable("TD-fctl ROOT 0 -> ROOT m ROTATORY STRENGTH (LEN) - h SYMMETRY")  # P::e SCF
        # wfn.set_variable("TD-fctl ROOT 0 -> ROOT m ROTATORY STRENGTH (VEL) - h SYMMETRY")  # P::e SCF
        wfn.set_variable(f"TD-{ssuper_name} ROOT {i+1} TOTAL ENERGY - {x.irrep_ES} SYMMETRY", E_tot_au)
        wfn.set_variable(f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} EXCITATION ENERGY - {x.irrep_ES} SYMMETRY", x.E_ex_au)
        wfn.set_variable(f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} OSCILLATOR STRENGTH (LEN) - {x.irrep_ES} SYMMETRY",
                         x.f_length)
        wfn.set_variable(f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} OSCILLATOR STRENGTH (VEL) - {x.irrep_ES} SYMMETRY",
                         x.f_velocity)
        wfn.set_variable(f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} ROTATORY STRENGTH (LEN) - {x.irrep_ES} SYMMETRY",
                         x.R_length)
        wfn.set_variable(f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} ROTATORY STRENGTH (VEL) - {x.irrep_ES} SYMMETRY",
                         x.R_velocity)
        wfn.set_array_variable(
            f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} ELECTRIC TRANSITION DIPOLE MOMENT (LEN) - {x.irrep_ES} SYMMETRY",
            core.Matrix.from_array(x.edtm_length.reshape((1, 3))))
        wfn.set_array_variable(
            f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} ELECTRIC TRANSITION DIPOLE MOMENT (VEL) - {x.irrep_ES} SYMMETRY",
            core.Matrix.from_array(x.edtm_velocity.reshape((1, 3))))
        wfn.set_array_variable(
            f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} MAGNETIC TRANSITION DIPOLE MOMENT - {x.irrep_ES} SYMMETRY",
            core.Matrix.from_array(x.mdtm.reshape((1, 3))))
        wfn.set_array_variable(
            f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} RIGHT EIGENVECTOR ALPHA - {x.irrep_ES} SYMMETRY",
            x.R_eigvec if restricted else x.R_eigvec[0])
        wfn.set_array_variable(f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} LEFT EIGENVECTOR ALPHA - {x.irrep_ES} SYMMETRY",
                               x.L_eigvec if restricted else x.L_eigvec[0])
        wfn.set_array_variable(f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} RIGHT EIGENVECTOR BETA - {x.irrep_ES} SYMMETRY",
                               x.R_eigvec if restricted else x.R_eigvec[1])
        wfn.set_array_variable(f"TD-{ssuper_name} ROOT 0 -> ROOT {i+1} LEFT EIGENVECTOR ALPHA - {x.irrep_ES} SYMMETRY",
                               x.L_eigvec if restricted else x.L_eigvec[1])

        core.print_out(
            f"    {i+1:^4} {sym_descr:^20} {x.E_ex_au:< 15.5f} {E_ex_ev:< 15.5f} {E_tot_au:< 15.5f} {x.f_length:< 15.4f} {x.f_velocity:< 15.4f} {x.R_length:< 15.4f} {x.R_velocity:< 15.4f}\n"
        )

    core.print_out("\n")

    return solver_results
