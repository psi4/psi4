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
    nbf = wfn.nmo()
    ndocc = wfn.nalpha()
    nvirt = nbf - ndocc

    c_occ = wfn.Ca_subset("AO", "OCC")
    c_vir = wfn.Ca_subset("AO", "VIR")

    # the vectors need to be in the MO basis. if they have the shape nbf x nbf, transform.
    for i in range(len(vectors)):
        shape = vectors[i].shape

        if shape == (nbf, nbf):
            vectors[i] = core.triplet(c_occ, vectors[i], c_vir, True, False, False)

        # verify that this vector already has the correct shape
        elif shape != (ndocc, nvirt):
            raise ValidationError('ERROR: "{}" has an unrecognized shape. Must be either ({}, {}) or ({}, {})'.format(
                vector_names[i], nbf, nbf, ndocc, nvirt))

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


def _print_tdscf_header(**options):
    core.print_out('\n\n         ---------------------------------------------------------\n'
                   '         {:^57}\n'.format('TDSCF excitation energies') +
                   '         {:^57}\n'.format('by Andrew M. James and Daniel G. A. Smith') +
                   '         ---------------------------------------------------------\n')

    core.print_out("{}\n".format("*"*80) +
                   "{}{:^60}{}\n".format("*"*10, "WARNING", "*"*10) +
                   "{}{:^60}{}\n".format("*"*10, "TDSCF is experimental results may be inaccurate", "*"*10) +
                   "{}\n".format("*"*80)) #yapf: disable

    core.print_out("\n  ==> Requested Excitations <==\n\n")
    state_info = options.pop('states')
    for nstate, state_sym in state_info:
        core.print_out("      {} states with {} symmetry\n".format(nstate, state_sym))

    core.print_out("\n  ==> Options <==\n\n")
    for k, v in options.items():
        core.print_out("     {:<10s}:              {}\n".format(k, v))

    core.print_out("\n")


def tdscf_excitations(wfn, **kwargs):
    """Compute excitations from a scf(HF/KS) wavefunction:

    Parameters
    -----------
    wfn : :py:class:`psi4.core.Wavefunction`
       The reference wavefunction
    states_per_irrep : list (int), {optional}
       The solver will find this many lowest excitations for each irreducible representation of the computational point group.
       The default is to find the lowest excitation of each symmetry. If this option is provided it must have the same number of elements
       as the number of irreducible representations as the computational point group.
    triplets : str {optional, ``none``, ``only``, ``also``}
       The default ``none`` will solve for no triplet states, ``only`` will solve for triplet states only, and ``also``
       will solve for the requested number of states of both triplet and singlet. This option is only valid for restricted references,
       and is ignored otherwise. The triplet and singlet solutions are found separately so using the ``also`` option will roughly double
       the computational cost of the calculation.
    tda :  bool {optional ``False``}
       If true the Tamm-Dancoff approximation (TDA) will be employed. For HF references this is equivalent to CIS.
    e_tol : float, {optional, 1.0e-6}
       The convergence threshold for the excitation energy
    r_tol : float, {optional, 1.0e-8}
       The convergence threshold for the norm of the residual vector
    max_ss_vectors: int {optional}
       The maximum number of ss vectors that will be stored before a collapse is done.
    guess : str
       If string the guess that will be used. Allowed choices:
       - ``denominators``: {default} uses orbital energy differences to generate guess vectors.


    ..note:: The algorithm employed to solve the non-Hermitian eigenvalue problem
             (when ``tda`` is False) will fail when the SCF wavefunction has a triplet instability.
    """
    # gather arguments
    e_tol = kwargs.pop('e_tol', 1.0e-6)
    r_tol = kwargs.pop('r_tol', 1.0e-8)
    max_ss_vec = kwargs.pop('max_ss_vectors', 50)
    verbose = kwargs.pop('print_lvl', 0)

    # how many states
    passed_spi = kwargs.pop('states_per_irrep', [0 for _ in range(wfn.nirrep())])

    #TODO:states_per_irrep = _validate_states_args(wfn, passed_spi, nstates)
    states_per_irrep = passed_spi

    #TODO: guess types, user guess
    guess_type = kwargs.pop("guess", "denominators")
    if guess_type != "denominators":
        raise ValidationError("Guess type {} is not valid".format(guess_type))

    # which problem
    ptype = 'rpa'
    solve_function = solvers.hamiltonian_solver
    if kwargs.pop('tda', False):
        ptype = 'tda'
        solve_function = solvers.davidson_solver

    restricted = wfn.same_a_b_orbs()
    if restricted:
        triplet = kwargs.pop('triplet', False)
    else:
        triplet = None

    _print_tdscf_header(
        etol=e_tol,
        rtol=r_tol,
        states=[(count, label) for count, label in zip(states_per_irrep,
                                                       wfn.molecule().irrep_labels())],
        guess_type=guess_type,
        restricted=restricted,
        triplet=triplet,
        ptype=ptype)

    # construct the engine
    if restricted:
        engine = TDRSCFEngine(wfn, triplet=triplet, ptype=ptype)
    else:
        engine = TDUSCFEngine(wfn, ptype=ptype)

    # just energies for now
    solver_results = []
    for state_sym, nstates in enumerate(states_per_irrep):
        if nstates == 0:
            continue
        engine.reset_for_state_symm(state_sym)
        guess_ = engine.generate_guess(nstates * 2)

        vecs_per_root = max_ss_vec // nstates

        # ret = (ee, rvecs, stats) (TDA)
        # ret = (ee, rvecs, lvecs, stats) (full TDSCF)
        ret = solve_function(
            engine=engine,
            e_tol=e_tol,
            r_tol=r_tol,
            max_vecs_per_root=vecs_per_root,
            nroot=nstates,
            guess=guess_,
            verbose=verbose)

        # store excitation energies tagged with final state symmetry (for printing)
        # TODO: handle R eigvecs (TDA) R/L eigvecs(full TDSCF): solver maybe should return dicts
        for ee in ret[0]:
            solver_results.append((ee, state_sym))

    # sort by energy symmetry is just meta data
    solver_results.sort(key=lambda x: x[0])

    # print excitation energies
    core.print_out("\n\nFinal Energetic Summary:\n")
    core.print_out("        " + (" " * 20) + " " + "Excitation Energy".center(31) + " {:^15}\n".format("Total Energy"))
    core.print_out("    {:^4} {:^20} {:^15} {:^15} {:^15}\n".format("#", "Sym: GS->ES (Trans)", "[au]", "[eV]",
                                                                    "(au)"))
    core.print_out("    {:->4} {:->20} {:->15} {:->15} {:->15}\n".format("-", "-", "-", "-", "-"))

    irrep_GS = wfn.molecule().irrep_labels()[engine.G_gs]
    for i, (E_ex_au, final_sym) in enumerate(solver_results):
        irrep_ES = wfn.molecule().irrep_labels()[final_sym]
        irrep_trans = wfn.molecule().irrep_labels()[engine.G_gs ^ final_sym]
        sym_descr = "{}->{} ({})".format(irrep_GS, irrep_ES, irrep_trans)

        #TODO: psivars/wfnvars

        E_ex_ev = constants.conversion_factor('hartree', 'eV') * E_ex_au

        E_tot_au = wfn.energy() + E_ex_au
        core.print_out("    {:^4} {:^20} {:< 15.5f} {:< 15.5f} {:< 15.5f}\n".format(
            i + 1, sym_descr, E_ex_au, E_ex_ev, E_tot_au))

    core.print_out("\n")

    #TODO: output table

    #TODO: oscillator strengths

    #TODO: check/handle convergence failures

    return solver_results
