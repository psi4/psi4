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

import math
import itertools

import numpy as np

from psi4 import core
from psi4.driver import p4util
from psi4.driver import constants
from psi4.driver.p4util.exceptions import *
from psi4.driver import driver_nbody_helper

### Math helper functions


def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n - r)


### Begin CBS gufunc data


def _sum_cluster_ptype_data(ptype,
                            ptype_dict,
                            compute_list,
                            fragment_slice_dict,
                            fragment_size_dict,
                            ret,
                            vmfc=False,
                            n=0):
    """
    Sums gradient and hessian data from compute_list.

    compute_list comes in as a tuple(frag, basis)
    """

    if len(compute_list) == 0:
        return

    sign = 1

    # Do ptype
    if ptype == 'gradient':
        for fragn, basisn in compute_list:
            start = 0
            grad = np.asarray(ptype_dict[(fragn, basisn)])

            if vmfc:
                sign = ((-1)**(n - len(fragn)))

            for bas in basisn:
                end = start + fragment_size_dict[bas]
                ret[fragment_slice_dict[bas]] += sign * grad[start:end]
                start += fragment_size_dict[bas]

    elif ptype == 'hessian':
        for fragn, basisn in compute_list:
            hess = np.asarray(ptype_dict[(fragn, basisn)])

            if vmfc:
                sign = ((-1)**(n - len(fragn)))

            # Build up start and end slices
            abs_start, rel_start = 0, 0
            abs_slices, rel_slices = [], []
            for bas in basisn:
                rel_end = rel_start + 3 * fragment_size_dict[bas]
                rel_slices.append(slice(rel_start, rel_end))
                rel_start += 3 * fragment_size_dict[bas]

                tmp_slice = fragment_slice_dict[bas]
                abs_slices.append(slice(tmp_slice.start * 3, tmp_slice.stop * 3))

            for abs_sl1, rel_sl1 in zip(abs_slices, rel_slices):
                for abs_sl2, rel_sl2 in zip(abs_slices, rel_slices):
                    ret[abs_sl1, abs_sl2] += hess[rel_sl1, rel_sl2]

    else:
        raise KeyError("ptype can only be gradient or hessian How did you end up here?")


def _print_nbody_energy(energy_body_dict, header, embedding=False):
    core.print_out("""\n   ==> N-Body: %s energies <==\n\n""" % header)
    core.print_out("""   n-Body     Total Energy [Eh]       I.E. [kcal/mol]      Delta [kcal/mol]\n""")
    previous_e = energy_body_dict[1]
    nbody_range = list(energy_body_dict)
    nbody_range.sort()
    for n in nbody_range:
        delta_e = (energy_body_dict[n] - previous_e)
        delta_e_kcal = delta_e * constants.hartree2kcalmol
        int_e_kcal = (
            energy_body_dict[n] - energy_body_dict[1]) * constants.hartree2kcalmol if not embedding else np.nan
        core.print_out(
            """     %4s  %20.12f  %20.12f  %20.12f\n""" % (n, energy_body_dict[n], int_e_kcal, delta_e_kcal))
        previous_e = energy_body_dict[n]
    core.print_out("\n")


def nbody_gufunc(func, method_string, **kwargs):
    """
    Computes the nbody interaction energy, gradient, or Hessian depending on input.
    This is a generalized univeral function for computing interaction and total quantities.

    :returns: *return type of func* |w--w| The data.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| data and wavefunction with energy/gradient/hessian set appropriately when **return_wfn** specified.

    :type func: function
    :param func: ``energy`` || etc.

        Python function that accepts method_string and a molecule. Returns a
        energy, gradient, or Hessian as requested.

    :type method_string: string
    :param method_string: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

        First argument, lowercase and usually unlabeled. Indicates the computational
        method to be passed to func.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type return_wfn: :ref:`boolean <op_py_boolean>`
    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|

        Indicate to additionally return the :py:class:`~psi4.core.Wavefunction`
        calculation result as the second element of a tuple.

    :type bsse_type: string or list
    :param bsse_type: ``'cp'`` || ``['nocp', 'vmfc']`` || |dl| ``None`` |dr| || etc.

        Type of BSSE correction to compute: CP, NoCP, or VMFC. The first in this
        list is returned by this function. By default, this function is not called.

    :type max_nbody: int
    :param max_nbody: ``3`` || etc.

        Maximum n-body to compute, cannot exceed the number of fragments in the moleucle.

    :type ptype: string
    :param ptype: ``'energy'`` || ``'gradient'`` || ``'hessian'``

        Type of the procedure passed in.

    :type return_total_data: :ref:`boolean <op_py_boolean>`
    :param return_total_data: ``'on'`` || |dl| ``'off'`` |dr|

        If True returns the total data (energy/gradient/etc) of the system,
        otherwise returns interaction data.

    :type levels: dict
    :param levels: ``{1: 'ccsd(t)', 2: 'mp2', 'supersystem': 'scf'}`` || ``{1: 2, 2: 'ccsd(t)', 3: 'mp2'}`` || etc

        Dictionary of different levels of theory for different levels of expansion
        Note that method_string is not used in this case.
        supersystem computes all higher order n-body effects up to nfragments.

    :type embedding_charges: dict
    :param embedding_charges: ``{1: [-0.834, 0.417, 0.417], ..}``

        Dictionary of atom-centered point charges. keys: 1-based index of fragment, values: list of charges for each fragment.

    :type charge_method: string
    :param charge_method: ``scf/6-31g`` || ``b3lyp/6-31g*`` || etc

        Method to compute point charges for monomers. Overridden by embedding_charges if both are provided.

    :type charge_type: string
    :param charge_type: ``MULLIKEN_CHARGES`` || ``LOWDIN_CHARGES`` 

        Default is ``MULLIKEN_CHARGES``
    """

    # Initialize dictionaries for easy data passing
    metadata, component_results, nbody_results = {}, {}, {}

    # Parse some kwargs
    kwargs = p4util.kwargs_lower(kwargs)
    if kwargs.get('levels', False):
        return driver_nbody_helper.multi_level(func, **kwargs)
    metadata['ptype'] = kwargs.pop('ptype', None)
    metadata['return_wfn'] = kwargs.pop('return_wfn', False)
    metadata['return_total_data'] = kwargs.pop('return_total_data', False)
    metadata['molecule'] = kwargs.pop('molecule', core.get_active_molecule())
    metadata['molecule'].update_geometry()
    metadata['molecule'].fix_com(True)
    metadata['molecule'].fix_orientation(True)
    metadata['embedding_charges'] = kwargs.get('embedding_charges', False)
    metadata['kwargs'] = kwargs
    core.clean_variables()

    if metadata['ptype'] not in ['energy', 'gradient', 'hessian']:
        raise ValidationError("""N-Body driver: The ptype '%s' is not regonized.""" % metadata['ptype'])

    # Parse bsse_type, raise exception if not provided or unrecognized
    metadata['bsse_type_list'] = kwargs.pop('bsse_type')
    if metadata['bsse_type_list'] is None:
        raise ValidationError("N-Body GUFunc: Must pass a bsse_type")
    if not isinstance(metadata['bsse_type_list'], list):
        metadata['bsse_type_list'] = [metadata['bsse_type_list']]

    for num, btype in enumerate(metadata['bsse_type_list']):
        metadata['bsse_type_list'][num] = btype.lower()
        if btype.lower() not in ['cp', 'nocp', 'vmfc']:
            raise ValidationError("N-Body GUFunc: bsse_type '%s' is not recognized" % btype.lower())

    metadata['max_nbody'] = kwargs.get('max_nbody', -1)
    metadata['max_frag'] = metadata['molecule'].nfragments()
    if metadata['max_nbody'] == -1:
        metadata['max_nbody'] = metadata['molecule'].nfragments()
    else:
        metadata['max_nbody'] = min(metadata['max_nbody'], metadata['max_frag'])

    # Flip this off for now, needs more testing
    # If we are doing CP lets save them integrals
    #if 'cp' in bsse_type_list and (len(bsse_type_list) == 1):
    #    # Set to save RI integrals for repeated full-basis computations
    #    ri_ints_io = core.get_global_option('DF_INTS_IO')

    #    # inquire if above at all applies to dfmp2 or just scf
    #    core.set_global_option('DF_INTS_IO', 'SAVE')
    #    psioh = core.IOManager.shared_object()
    #    psioh.set_specific_retention(97, True)

    bsse_str = metadata['bsse_type_list'][0]
    if len(metadata['bsse_type_list']) > 1:
        bsse_str = str(metadata['bsse_type_list'])
    core.print_out("\n\n")
    core.print_out("   ===> N-Body Interaction Abacus <===\n")
    core.print_out("        BSSE Treatment:                     %s\n" % bsse_str)

    # Get compute list
    metadata = build_nbody_compute_list(metadata)

    # Compute N-Body components
    component_results = compute_nbody_components(func, method_string, metadata)

    # Assemble N-Body quantities
    nbody_results = assemble_nbody_components(metadata, component_results)

    # Build wfn and bind variables
    wfn = core.Wavefunction.build(metadata['molecule'], 'def2-svp')
    dicts = [
        'energies', 'ptype', 'intermediates', 'energy_body_dict', 'gradient_body_dict', 'hessian_body_dict', 'nbody',
        'cp_energy_body_dict', 'nocp_energy_body_dict', 'vmfc_energy_body_dict'
    ]
    if metadata['ptype'] == 'gradient':
        wfn.set_gradient(nbody_results['ret_ptype'])
        nbody_results['gradient_body_dict'] = nbody_results['ptype_body_dict']
    elif metadata['ptype'] == 'hessian':
        nbody_results['hessian_body_dict'] = nbody_results['ptype_body_dict']
        wfn.set_hessian(nbody_results['ret_ptype'])
        component_results_gradient = component_results.copy()
        component_results_gradient['ptype'] = component_results_gradient['gradients']
        metadata['ptype'] = 'gradient'
        nbody_results_gradient = assemble_nbody_components(metadata, component_results_gradient)
        wfn.set_gradient(nbody_results_gradient['ret_ptype'])
        nbody_results['gradient_body_dict'] = nbody_results_gradient['ptype_body_dict']

    for r in [component_results, nbody_results]:
        for d in r:
            if d in dicts:
                for var, value in r[d].items():
                    try:
                        wfn.set_scalar_variable(str(var), value)
                        core.set_scalar_variable(str(var), value)
                    except:
                        wfn.set_array_variable(d.split('_')[0].upper() + ' ' + str(var), core.Matrix.from_array(value))

    core.set_variable("CURRENT ENERGY", nbody_results['ret_energy'])
    wfn.set_variable("CURRENT ENERGY", nbody_results['ret_energy'])

    if metadata['return_wfn']:
        return (nbody_results['ret_ptype'], wfn)
    else:
        return nbody_results['ret_ptype']


def build_nbody_compute_list(metadata):
    """Generates the list of N-Body computations to be performed for a given BSSE type.

    Parameters
    ----------
    metadata : dict of str
        Dictionary containing N-body metadata.

        Required ``'key': value`` pairs:
        ``'bsse_type_list'``: list of str
            List of requested BSSE treatments.  Possible values include lowercase ``'cp'``, ``'nocp'``,
            and ``'vmfc'``.
        ``'max_nbody'``: int
            Maximum number of bodies to include in the N-Body treatment.
            Possible: `max_nbody` <= `max_frag`
            Default: `max_nbody` = `max_frag`
        ``'max_frag'``: int
            Number of distinct fragments comprising full molecular supersystem.

    Returns
    -------
    metadata : dict of str
        Dictionary containing N-body metadata.
        
        New ``'key': value`` pair:
        ``'compute_dict'`` : dict of str: dict
            Dictionary containing subdicts enumerating compute lists for each possible BSSE treatment.

            Contents:
            ``'all'``: dict of int: set
                Set containing full list of computations required
            ``'cp'``: dict of int: set
                Set containing list of computations required for CP procedure
            ``'nocp'``: dict of int: set
                Set containing list of computations required for non-CP procedure
            ``'vmfc_compute'``: dict of int: set
                Set containing list of computations required for VMFC procedure
            ``'vmfc_levels'``: dict of int: set
                Set containing list of levels required for VMFC procedure
    """
    # What levels do we need?
    nbody_range = range(1, metadata['max_nbody'] + 1)
    fragment_range = range(1, metadata['max_frag'] + 1)

    cp_compute_list = {x: set() for x in nbody_range}
    nocp_compute_list = {x: set() for x in nbody_range}
    vmfc_compute_list = {x: set() for x in nbody_range}
    vmfc_level_list = {x: set() for x in nbody_range}  # Need to sum something slightly different

    # Verify proper passing of bsse_type_list
    bsse_type_remainder = set(metadata['bsse_type_list']) - {'cp', 'nocp', 'vmfc'}
    if bsse_type_remainder:
        raise ValidationError("""Unrecognized BSSE type(s): %s
Possible values are 'cp', 'nocp', and 'vmfc'.""" % ', '.join(str(i) for i in bsse_type_remainder))

    # Build up compute sets
    if 'cp' in metadata['bsse_type_list']:
        # Everything is in dimer basis
        basis_tuple = tuple(fragment_range)
        for nbody in nbody_range:
            for x in itertools.combinations(fragment_range, nbody):
                if metadata['max_nbody'] == 1: break
                cp_compute_list[nbody].add((x, basis_tuple))
        # Add monomers in monomer basis
        for x in fragment_range:
            cp_compute_list[1].add(((x, ), (x, )))

    if 'nocp' in metadata['bsse_type_list']:
        # Everything in monomer basis
        for nbody in nbody_range:
            for x in itertools.combinations(fragment_range, nbody):
                nocp_compute_list[nbody].add((x, x))

    if 'vmfc' in metadata['bsse_type_list']:
        # Like a CP for all combinations of pairs or greater
        for nbody in nbody_range:
            for cp_combos in itertools.combinations(fragment_range, nbody):
                basis_tuple = tuple(cp_combos)
                for interior_nbody in nbody_range:
                    for x in itertools.combinations(cp_combos, interior_nbody):
                        combo_tuple = (x, basis_tuple)
                        vmfc_compute_list[interior_nbody].add(combo_tuple)
                        vmfc_level_list[len(basis_tuple)].add(combo_tuple)

    # Build a comprehensive compute_range
    compute_list = {x: set() for x in nbody_range}
    for n in nbody_range:
        compute_list[n] |= cp_compute_list[n]
        compute_list[n] |= nocp_compute_list[n]
        compute_list[n] |= vmfc_compute_list[n]
        core.print_out("        Number of %d-body computations:     %d\n" % (n, len(compute_list[n])))

    metadata['compute_dict'] = {
        'all': compute_list,
        'cp': cp_compute_list,
        'nocp': nocp_compute_list,
        'vmfc_compute': vmfc_compute_list,
        'vmfc_levels': vmfc_level_list
    }

    return metadata


def compute_nbody_components(func, method_string, metadata):
    """Computes requested N-body components.

    Performs requested computations for psi4::Molecule object `molecule` according to
    `compute_list` with function `func` at `method_string` level of theory.

    Parameters
    ----------
    func : {'energy', 'gradient', 'hessian'}
        Function object to be called within N-Body procedure.
    method_string : str
        Indicates level of theory to be passed to function `func`.
    metadata : dict of str
        Dictionary of N-body metadata.

        Required ``'key': value`` pairs:
        ``'compute_list'``: dict of int: set
            List of computations to perform.  Keys indicate body-levels, e.g,. `compute_list[2]` is the
            list of all 2-body computations required.
        ``'kwargs'``: dict
            Arbitrary keyword arguments to be passed to function `func`.

    Returns
    -------
    dict of str: dict
        Dictionary containing computed N-body components.

        Contents:
        ``'energies'``: dict of set: float64
               Dictionary containing all energy components required for given N-body procedure.
        ``'ptype'``: dict of set: float64 or dict of set: psi4.Matrix
               Dictionary of returned quantities from calls of function `func` during N-body computations
        ``'intermediates'``: dict of str: float64
               Dictionary of psivars for intermediate N-body computations to be set at the end of the
               N-body procedure.
    """
    # Get required metadata
    kwargs = metadata['kwargs']
    molecule = metadata['molecule']
    #molecule = core.get_active_molecule()
    compute_list = metadata['compute_dict']['all']

    # Now compute the energies
    energies_dict = {}
    gradients_dict = {}
    ptype_dict = {}
    intermediates_dict = {}
    if kwargs.get('charge_method', False) and not metadata['embedding_charges']:
        metadata['embedding_charges'] = driver_nbody_helper.compute_charges(kwargs['charge_method'],
                                            kwargs.get('charge_type', 'MULLIKEN_CHARGES').upper(), molecule)
    for count, n in enumerate(compute_list.keys()):
        core.print_out("\n   ==> N-Body: Now computing %d-body complexes <==\n\n" % n)
        total = len(compute_list[n])
        for num, pair in enumerate(compute_list[n]):
            core.print_out(
                "\n       N-Body: Computing complex (%d/%d) with fragments %s in the basis of fragments %s.\n\n" %
                (num + 1, total, str(pair[0]), str(pair[1])))
            ghost = list(set(pair[1]) - set(pair[0]))

            current_mol = molecule.extract_subsets(list(pair[0]), ghost)
            current_mol.set_name("%s_%i_%i" % (current_mol.name(), count, num))
            if metadata['embedding_charges']: driver_nbody_helper.electrostatic_embedding(metadata, pair=pair)
            # Save energies info
            ptype_dict[pair], wfn = func(method_string, molecule=current_mol, return_wfn=True, **kwargs)
            core.set_global_option_python('EXTERN', None)
            energies_dict[pair] = core.variable("CURRENT ENERGY")
            gradients_dict[pair] = wfn.gradient()
            var_key = "N-BODY (%s)@(%s) TOTAL ENERGY" % (', '.join([str(i) for i in pair[0]]), ', '.join(
                [str(i) for i in pair[1]]))
            intermediates_dict[var_key] = core.variable("CURRENT ENERGY")
            core.print_out("\n       N-Body: Complex Energy (fragments = %s, basis = %s: %20.14f)\n" % (str(
                pair[0]), str(pair[1]), energies_dict[pair]))
            # Flip this off for now, needs more testing
            #if 'cp' in bsse_type_list and (len(bsse_type_list) == 1):
            #    core.set_global_option('DF_INTS_IO', 'LOAD')

            core.clean()

    return {
        'energies': energies_dict,
        'gradients': gradients_dict,
        'ptype': ptype_dict,
        'intermediates': intermediates_dict
    }


def assemble_nbody_components(metadata, component_results):
    """Assembles N-body components into interaction quantities according to requested BSSE procedure(s).

    Parameters
    -----------
    metadata : dict of str
        Dictionary of N-body metadata.

        Required ``'key': value`` pairs:
            ``'ptype'``: {'energy', 'gradient', 'hessian'}
                   Procedure which has generated the N-body components to be combined.
            ``'bsse_type_list'``: list of str
                   List of requested BSSE treatments.  Possible values include lowercase ``'cp'``, ``'nocp'``,
                   and ``'vmfc'``.
            ``'max_nbody'``: int
                   Maximum number of bodies to include in the N-Body treatment.
                   Possible: `max_nbody` <= `max_frag`
                   Default: `max_nbody` = `max_frag`
            ``'max_frag'``: int
                   Number of distinct fragments comprising full molecular supersystem.
            ``'energies_dict'``: dict of set: float64
                   Dictionary containing all energy components required for given N-body procedure.
            ``'ptype_dict'``: dict of set: float64 or dict of set: psi4.Matrix
                   Dictionary of returned quantities from calls of function `func` during N-body computations
            ``'compute_dict'``: dict of str: dict
                   Dictionary containing {int: set} subdicts enumerating compute lists for each possible
                   BSSE treatment.
            ``'kwargs'``: dict
                   Arbitrary keyword arguments.
    component_results : dict of str: dict
        Dictionary containing computed N-body components.

        Required ``'key': value`` pairs:
        ``'energies'``: dict of set: float64
               Dictionary containing all energy components required for given N-body procedure.
        ``'ptype'``: dict of set: float64 or dict of set: psi4.Matrix
               Dictionary of returned quantities from calls of function `func` during N-body computations
        ``'intermediates'``: dict of str: float64
               Dictionary of psivars for intermediate N-body computations to be set at the end of the
               N-body procedure.
    Returns
    -------
    results : dict of str
        Dictionary of all N-body results.

        Contents:
        ``'ret_energy'``: float64
            Interaction data requested.  If multiple BSSE types requested in `bsse_type_list`, the
            interaction data associated with the *first* BSSE type in the list is returned.
        ``'nbody_dict'``: dict of str: float64
            Dictionary of relevant N-body psivars to be set
        ``'energy_body_dict'``: dict of int: float64
            Dictionary of total energies at each N-body level, i.e., ``results['energy_body_dict'][2]`` 
            is the sum of all 2-body total energies for the supersystem.
        ``'ptype_body_dict'``: dict or dict of int: array_like
            Empty dictionary if `ptype is ``'energy'``, or dictionary of total ptype arrays at each
            N-body level; i.e., ``results['ptype_body_dict'][2]`` for `ptype` ``'gradient'``is the 
            total 2-body gradient.
    """
    # Unpack metadata
    kwargs = metadata['kwargs']

    nbody_range = range(1, metadata['max_nbody'] + 1)

    # Unpack compute list metadata
    compute_list = metadata['compute_dict']['all']
    cp_compute_list = metadata['compute_dict']['cp']
    nocp_compute_list = metadata['compute_dict']['nocp']
    vmfc_compute_list = metadata['compute_dict']['vmfc_compute']
    vmfc_level_list = metadata['compute_dict']['vmfc_levels']

    # Build size and slices dictionaries
    fragment_size_dict = {
        frag: metadata['molecule'].extract_subsets(frag).natom()
        for frag in range(1, metadata['max_frag'] + 1)
    }
    start = 0
    fragment_slice_dict = {}
    for k, v in fragment_size_dict.items():
        fragment_slice_dict[k] = slice(start, start + v)
        start += v

    molecule_total_atoms = sum(fragment_size_dict.values())

    # Final dictionaries
    cp_energy_by_level = {n: 0.0 for n in nbody_range}
    nocp_energy_by_level = {n: 0.0 for n in nbody_range}

    cp_energy_body_dict = {n: 0.0 for n in nbody_range}
    nocp_energy_body_dict = {n: 0.0 for n in nbody_range}
    vmfc_energy_body_dict = {n: 0.0 for n in nbody_range}

    # Build out ptype dictionaries if needed
    if metadata['ptype'] != 'energy':
        if metadata['ptype'] == 'gradient':
            arr_shape = (molecule_total_atoms, 3)
        elif metadata['ptype'] == 'hessian':
            arr_shape = (molecule_total_atoms * 3, molecule_total_atoms * 3)
        else:
            raise KeyError("N-Body: ptype '%s' not recognized" % ptype)

        cp_ptype_by_level = {n: np.zeros(arr_shape) for n in nbody_range}
        nocp_ptype_by_level = {n: np.zeros(arr_shape) for n in nbody_range}
        vmfc_ptype_by_level = {n: np.zeros(arr_shape) for n in nbody_range}

        cp_ptype_body_dict = {n: np.zeros(arr_shape) for n in nbody_range}
        nocp_ptype_body_dict = {n: np.zeros(arr_shape) for n in nbody_range}
        vmfc_ptype_body_dict = {n: np.zeros(arr_shape) for n in nbody_range}
    else:
        cp_ptype_by_level, cp_ptype_body_dict = {}, {}
        nocp_ptype_by_level, nocp_ptype_body_dict = {}, {}
        vmfc_ptype_body_dict = {}

    # Sum up all of the levels
    nbody_dict = {}
    for n in nbody_range:

        # Energy
        # Extract energies for monomers in monomer basis for CP total data
        if n == 1:
            cp_monomers_in_monomer_basis = [v for v in cp_compute_list[1] if len(v[1]) == 1]
            cp_monomer_energies = 0.0
            cp_monomer_energy_list = []
            for i in cp_monomers_in_monomer_basis:
                cp_monomer_energy_list.append(component_results['energies'][i])
                cp_monomer_energies += component_results['energies'][i]
                cp_compute_list[1].remove(i)

        cp_energy_by_level[n] = sum(component_results['energies'][v] for v in cp_compute_list[n])
        nocp_energy_by_level[n] = sum(component_results['energies'][v] for v in nocp_compute_list[n])

        # Special vmfc case
        if n > 1:
            vmfc_energy_body_dict[n] = vmfc_energy_body_dict[n - 1]
        for tup in vmfc_level_list[n]:
            vmfc_energy_body_dict[n] += ((-1)**(n - len(tup[0]))) * component_results['energies'][tup]

        # Do ptype
        if metadata['ptype'] != 'energy':
            _sum_cluster_ptype_data(metadata['ptype'], component_results['ptype'], cp_compute_list[n],
                                    fragment_slice_dict, fragment_size_dict, cp_ptype_by_level[n])
            _sum_cluster_ptype_data(metadata['ptype'], component_results['ptype'], nocp_compute_list[n],
                                    fragment_slice_dict, fragment_size_dict, nocp_ptype_by_level[n])
            _sum_cluster_ptype_data(
                metadata['ptype'],
                component_results['ptype'],
                vmfc_level_list[n],
                fragment_slice_dict,
                fragment_size_dict,
                vmfc_ptype_by_level[n],
                vmfc=True,
                n=n)

        # Add extracted monomers back.
        for i, j in enumerate(cp_monomers_in_monomer_basis):
            cp_compute_list[1].add(j)

    if metadata['ptype'] != 'energy':
        # Extract ptype data for monomers in monomer basis for CP total data
        cp_monomer_ptype = np.zeros(arr_shape)
        _sum_cluster_ptype_data(metadata['ptype'], component_results['ptype'], cp_monomers_in_monomer_basis,
                                fragment_slice_dict, fragment_size_dict, cp_monomer_ptype)

    # Compute cp energy and ptype
    if 'cp' in metadata['bsse_type_list']:
        for n in nbody_range:
            if n == metadata['max_frag']:
                cp_energy_body_dict[n] = cp_energy_by_level[n] - bsse
                if metadata['ptype'] != 'energy':
                    cp_ptype_body_dict[n][:] = cp_ptype_by_level[n] - bsse_ptype
                continue

            for k in range(1, n + 1):
                take_nk = nCr(metadata['max_frag'] - k - 1, n - k)
                sign = ((-1)**(n - k))
                value = cp_energy_by_level[k]
                cp_energy_body_dict[n] += take_nk * sign * value

                if metadata['ptype'] != 'energy':
                    value = cp_ptype_by_level[k]
                    cp_ptype_body_dict[n] += take_nk * sign * value

            if n == 1:
                bsse = cp_energy_body_dict[n] - cp_monomer_energies
                cp_energy_body_dict[n] = cp_monomer_energies
                if metadata['ptype'] != 'energy':
                    bsse_ptype = cp_ptype_body_dict[n] - cp_monomer_ptype
                    cp_ptype_body_dict[n] = cp_monomer_ptype.copy()

            else:
                cp_energy_body_dict[n] -= bsse
                if metadata['ptype'] != 'energy':
                    cp_ptype_body_dict[n] -= bsse_ptype

        cp_interaction_energy = cp_energy_body_dict[metadata['max_nbody']] - cp_energy_body_dict[1]
        nbody_dict['Counterpoise Corrected Interaction Energy'] = cp_interaction_energy

        for n in nbody_range[1:]:
            var_key = 'CP-CORRECTED %d-BODY INTERACTION ENERGY' % n
            nbody_dict[var_key] = cp_energy_body_dict[n] - cp_energy_body_dict[1]

        _print_nbody_energy(cp_energy_body_dict, "Counterpoise Corrected (CP)", metadata['embedding_charges'])
        cp_interaction_energy = cp_energy_body_dict[metadata['max_nbody']] - cp_energy_body_dict[1]
        nbody_dict['Counterpoise Corrected Total Energy'] = cp_energy_body_dict[metadata['max_nbody']]
        nbody_dict['Counterpoise Corrected Interaction Energy'] = cp_interaction_energy

    # Compute nocp energy and ptype
    if 'nocp' in metadata['bsse_type_list']:
        for n in nbody_range:
            if n == metadata['max_frag']:
                nocp_energy_body_dict[n] = nocp_energy_by_level[n]
                if metadata['ptype'] != 'energy':
                    nocp_ptype_body_dict[n][:] = nocp_ptype_by_level[n]
                continue

            for k in range(1, n + 1):
                take_nk = nCr(metadata['max_frag'] - k - 1, n - k)
                sign = ((-1)**(n - k))
                value = nocp_energy_by_level[k]
                nocp_energy_body_dict[n] += take_nk * sign * value

                if metadata['ptype'] != 'energy':
                    value = nocp_ptype_by_level[k]
                    nocp_ptype_body_dict[n] += take_nk * sign * value

        _print_nbody_energy(nocp_energy_body_dict, "Non-Counterpoise Corrected (NoCP)", metadata['embedding_charges'])
        nocp_interaction_energy = nocp_energy_body_dict[metadata['max_nbody']] - nocp_energy_body_dict[1]
        nbody_dict['Non-Counterpoise Corrected Total Energy'] = nocp_energy_body_dict[metadata['max_nbody']]
        nbody_dict['Non-Counterpoise Corrected Interaction Energy'] = nocp_interaction_energy

        for n in nbody_range[1:]:
            var_key = 'NOCP-CORRECTED %d-BODY INTERACTION ENERGY' % n
            nbody_dict[var_key] = nocp_energy_body_dict[n] - nocp_energy_body_dict[1]

    # Compute vmfc ptype
    if 'vmfc' in metadata['bsse_type_list']:
        if metadata['ptype'] != 'energy':
            for n in nbody_range:
                if n > 1:
                    vmfc_ptype_body_dict[n] = vmfc_ptype_by_level[n-1]
                vmfc_ptype_body_dict[n] += vmfc_ptype_by_level[n]

        _print_nbody_energy(vmfc_energy_body_dict, "Valiron-Mayer Function Couterpoise (VMFC)",
                            metadata['embedding_charges'])
        vmfc_interaction_energy = vmfc_energy_body_dict[metadata['max_nbody']] - vmfc_energy_body_dict[1]
        nbody_dict['Valiron-Mayer Function Couterpoise Total Energy'] = vmfc_energy_body_dict[metadata['max_nbody']]
        nbody_dict['Valiron-Mayer Function Couterpoise Interaction Energy'] = vmfc_interaction_energy

        for n in nbody_range[1:]:
            var_key = 'VMFC-CORRECTED %d-BODY INTERACTION ENERGY' % n
            nbody_dict[var_key] = vmfc_energy_body_dict[n] - vmfc_energy_body_dict[1]

    # Returns
    results = {}
    results['nbody'] = nbody_dict
    for b in ['cp', 'nocp', 'vmfc']:
        results['%s_energy_body_dict' % b] = eval('%s_energy_body_dict' % b)
        results['%s_energy_body_dict' % b] = {str(i) + b: j for i, j in results['%s_energy_body_dict' % b].items()}

    # Figure out and build return types
    return_method = metadata['bsse_type_list'][0]

    if return_method == 'cp':
        results['ptype_body_dict'] = cp_ptype_body_dict
        results['energy_body_dict'] = cp_energy_body_dict
    elif return_method == 'nocp':
        results['ptype_body_dict'] = nocp_ptype_body_dict
        results['energy_body_dict'] = nocp_energy_body_dict
    elif return_method == 'vmfc':
        results['ptype_body_dict'] = vmfc_ptype_body_dict
        results['energy_body_dict'] = vmfc_energy_body_dict
    else:
        raise ValidationError(
            "N-Body Wrapper: Invalid return type. Should never be here, please post this error on github.")

    if metadata['return_total_data']:
        results['ret_energy'] = results['energy_body_dict'][metadata['max_nbody']]
    else:
        results['ret_energy'] = results['energy_body_dict'][metadata['max_nbody']]
        results['ret_energy'] -= results['energy_body_dict'][1]

    if metadata['ptype'] != 'energy':
        if metadata['return_total_data']:
            np_final_ptype = results['ptype_body_dict'][metadata['max_nbody']].copy()
        else:
            np_final_ptype = results['ptype_body_dict'][metadata['max_nbody']].copy()
            np_final_ptype -= results['ptype_body_dict'][1]

        results['ret_ptype'] = core.Matrix.from_array(np_final_ptype)
    else:
        results['ret_ptype'] = results['ret_energy']

    return results
