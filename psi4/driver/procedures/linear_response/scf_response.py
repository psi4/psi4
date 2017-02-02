import numpy as np
from psi4 import core

dipole = {'name': 'Dipole polarizabilities', 'printout labels': ['X', 'Y', 'Z'],
          'mints function': core.MintsHelper.so_dipole, 'vector names': ["SO Dipole x", "SO Dipole y", "SO Dipole z"]}

quadrupole = {'name': 'Quadrupole polarizabilities', 'printout labels': ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'],
              'mints function': core.MintsHelper.so_quadrupole,
              'vector names': ['SO Quadrupole x2', 'SO Quadrupole xy', 'SO Quadrupole xz', 'SO Quadrupole y2',
                               'SO Quadrupole yz', 'SO Quadrupole z2']}

traceless_quadrupole = {'name': 'Traceless quadrupole polarizabilities',
                        'printout labels': ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'],
                        'mints function': core.MintsHelper.so_traceless_quadrupole,
                        'vector names': ['SO Traceless Quadrupole x2', 'SO Traceless Quadrupole xy',
                                         'SO Traceless Quadrupole xz', 'SO Traceless Quadrupole y2',
                                         'SO Traceless Quadrupole yz', 'SO Traceless Quadrupole z2']}

property_dicts = {'dipole_polarizabilities': dipole, 'quadrupole_polarizabilities': quadrupole,
                  'traceless_quadrupole_polarizabilities': traceless_quadrupole}


def static_solve(wfn, *args, **kwargs):
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
        argtype = type(arg)

        # for each string keyword, append the appropriate dictionary (vide supra) to our list
        if argtype is str:
            ret = property_dicts.get(arg)
            if ret:
                complete_dict.append(ret)
            else:
                print('Do not understand {}. Skipping.'.format(arg))

        # the user passed a list of vectors. absorb them into a dictionary
        elif argtype is tuple or argtype is list:
            complete_dict.append(
                {'name': 'User Vectors', 'length': len(arg), 'vectors': arg,
                 'vector names': ['User Vector {}_{}'.format(n_user, i) for i in range(len(arg))]})
            n_user += len(arg)

        # single vector passed. stored in a dictionary as a list of length 1 (can be handled as the case above that way)
        # note: the length is set to '0' to designate that it was not really passed as a list
        else:
            complete_dict.append({'name': 'User Vector', 'length': 0, 'vectors': [arg],
                                  'vector names': ['User Vector {}'.format(n_user)]})
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
            tmp_vectors = prop['mints function'](mints)
            for tmp in tmp_vectors:
                tmp.scale(-2.0)
                vectors.append(tmp)
                vector_names.append(tmp.name)

    # do we have any vectors to work with?
    if len(vectors) == 0:
        print('I have no vectors to work with. Aborting.')
        return None

    # print information on module, vectors that will be used
    _print_header(complete_dict, n_user)

    # fetch wavefunction information
    nbf = wfn.nmo()
    ndocc = wfn.nalpha()
    nvirt = nbf - ndocc

    c_occ = wfn.Ca_subset("SO", "OCC")
    c_vir = wfn.Ca_subset("SO", "VIR")

    # the vectors need to be in the MO basis. if they have the shape nbf x nbf, transform.
    for i in range(len(vectors)):
        shape = vectors[i].shape

        if shape == (nbf, nbf):
            vectors[i] = core.Matrix.triplet(c_occ, vectors[i], c_vir, True, False, False)

        # verify that this vector already has the correct shape
        elif shape != (ndocc, nvirt):
            core.print_out('ERROR: "{}" has an unrecognized shape. Must be either ({}, {}) or ({}, {})\n'.format(
                vector_names[i], nbf, nbf, ndocc, nvirt
            ))

    # compute response vectors for each input vector
    params = _get_cphf_solver_params(kwargs)        # TODO get the options from Psi4?
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
                   '         {:^57}\n'.format('CPHF Solver') +
                   '         {:^57}\n'.format('by Marvin Lechner and Daniel Smith') +
                   '         ---------------------------------------------------------\n')

    core.print_out('\n   ==> Requested Responses <==\n\n')

    for prop in complete_dict:
        if 'User' not in prop['name']:
            core.print_out('    {}\n'.format(prop['name']))

    if n_user != 0:
        core.print_out('    {} user-supplied vector(s)\n'.format(n_user))


def _print_matrix(descriptors, content):
    length = len(descriptors)

    matrix_header = '         ' + ' {:^10}' * length + '\n'
    core.print_out(matrix_header.format(*descriptors))
    core.print_out('    -----' + ' ----------' * length + '\n')

    for i, desc in enumerate(descriptors):
        core.print_out('    {:^5}'.format(desc))
        for j in range(length):
            core.print_out(' {:>10.5f}'.format(content[i, j]))
        core.print_out('\n')


def _print_output(complete_dict, output):
    core.print_out('\n   ==> Response Properties <==\n')

    for i, prop in enumerate(complete_dict):
        if not 'User' in prop['name']:
            core.print_out('\n    => {} <=\n\n'.format(prop['name']))
            directions = prop['printout labels']
            _print_matrix(directions, output[i])


def _get_cphf_solver_params(kwargs):
    """
    Extract the parameters for the cphf solver from the keyword arguments.

    Parameters
    ----------
    kwargs : dict
        The keyword parameters passed by the user.

    Returns
    -------
    params : list
        The cphf solver parameters ``conv_tol``, ``max_iter``, ``print_lvl``.
    """
    params = []
    needed_params = ['conv_tol', 'max_iter', 'print_lvl']
    default_values = [1e-5, 10, 2]

    for i, param in enumerate(needed_params):
        ret = kwargs.get(param)

        if ret is None:
            params.append(default_values[i])
        else:
            params.append(ret)

    return params
