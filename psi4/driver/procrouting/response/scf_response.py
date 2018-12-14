#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2018 The Psi4 Developers.
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
from psi4.driver.p4util import solvers
from psi4.driver.p4util.exceptions import *
from psi4.driver.procrouting.response.scf_products import (TDRSCFEngine,
                                                           TDUSCFEngine)

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


def _single_random_fill(mat):
    for h in range(mat.nirrep()):
        shape = mat.nph[h].shape
        mat.nph[h][:, :] = np.random.randn(*shape)
    return mat


def _denom_guess_uhf(engine, nstates):
    guess_vecs = []
    for ex_spin in range(2):
        for irr_occ in range(engine.wfn.nirrep()):
            irr_vir = irr_occ ^ engine.Gx
            begin_occ = max(engine.occpi[ex_spin][irr_occ] - nstates * 2, 0)
            end_vir = min(engine.virpi[ex_spin][irr_vir], nstates * 2)
            for i in range(begin_occ, engine.occpi[ex_spin][irr_occ]):
                for a in range(0, end_vir):
                    new_vec = engine.new_vector("Guess vector")
                    new_vec[ex_spin].set(irr_occ, i, a, 1.0)
                    guess_vecs.append(new_vec)
    return guess_vecs


def _denom_guess_rhf(engine, nstates):
    guess_vecs = []
    for irr_occ in range(engine.wfn.nirrep()):
        irr_vir = irr_occ ^ engine.Gx
        begin_occ = max(engine.occpi[irr_occ] - nstates * 2, 0)
        end_vir = min(engine.virpi[irr_vir], nstates * 2)
        for i in range(begin_occ, engine.occpi[irr_occ]):
            for a in range(0, end_vir):
                new_vec = engine.new_vector("Guess vector")
                new_vec.set(irr_occ, i, a, 1.0)
                guess_vecs.append(new_vec)
    return guess_vecs


def tdscf_excitations(wfn, **kwargs):
    """Compute excitations from a scf(HF/KS) wavefunction:

    Parameters
    -----------
    wfn : :py:class:`psi4.core.Wavefunction`
       The reference wavefunction
    nstates : int {optional}
       The number of states to find, if the system has symmetry this is the same as passing ``states_per_irrep`` with ``nstates``
       as every element. (Convenient for C1 systems)
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
       - ``random``: generates random vectors


    ..note:: The algorithm employed to solve the non-Hermitian eigenvalue problem
             (when ``tda`` is False) will fail when the SCF wavefunction has a triplet instability.
    """
    # gather arguments
    etol = kwargs.pop('e_tol', 1.0e-6)
    rtol = kwargs.pop('r_tol', 1.0e-8)
    max_ss_vec = kwargs.pop('max_ss_vectors', 50)

    # how many
    passed_spi = kwargs.pop('states_per_irrep', [])
    nstates = kwargs.pop('nstates', 1)
    #states_per_irrep = _validate_states_args(wfn, passed_spi, nstates)
    states_per_irrep = passed_spi

    # guess
    passed_guess = kwargs.pop('guess', 'random')
    guess_type = None
    if isinstance(passed_guess, (list, tuple)):
        guess_vectors = _validate_user_guess(wfn, states_per_irrep, guess)
        guess_type = 'user'
    else:
        guess_vectors = []
        guess_type = 'random'

    # which problem
    ptype = 'rpa'
    solve_function = solvers.hamiltonian_solver
    if kwargs.pop('tda', False):
        ptype = 'tda'
        solve_function = solvers.davidson_solver

    # construct the engine
    if wfn.same_a_b_orbs():
        guess_gen = _denom_guess_rhf
        engine = TDRSCFEngine(wfn, triplet=kwargs.pop('triplet', False), ptype=ptype)
    else:
        guess_gen = _denom_guess_uhf
        engine = TDUSCFEngine(wfn, ptype=ptype)

    solver_results = []
    for state_sym, nstates in enumerate(states_per_irrep):
        if nstates == 0:
            solver_results.append(None)
            continue
        engine.reset_symmetry(state_sym)
        guess_ = []
        if guess_type == 'random':
            guess_ = guess_gen(engine, nstates)
        elif guess_type == 'user':
            guess_ = guess_vectors[state_sym]

        vecs_per_root = max_ss_vec // nstates
        solver_results.append(
            solve_function(
                engine=engine,
                e_tol=etol,
                r_tol=rtol,
                max_vecs_per_root=vecs_per_root,
                nroot=nstates,
                guess=guess_,
                verbose=2))

    return solver_results
