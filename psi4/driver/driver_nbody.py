#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

from __future__ import print_function
from __future__ import absolute_import

import math
import numpy as np
import itertools as it

from psi4 import core

# Import driver helpers
from psi4.driver import p4util
from psi4.driver import constants

from psi4.driver.p4util.exceptions import *

### Math helper functions

def nCr(n, r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

### Begin CBS gufunc data

def _sum_cluster_ptype_data(ptype, ptype_dict, compute_list, fragment_slice_dict, fragment_size_dict, ret, vmfc=False):
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
                sign = ((-1) ** (n - len(fragn)))

            for bas in basisn:
                end = start + fragment_size_dict[bas]
                ret[fragment_slice_dict[bas]] += sign * grad[start:end]
                start += fragment_size_dict[bas]

    elif ptype == 'hessian':
        for fragn, basisn in compute_list:
            hess = np.asarray(ptype_dict[(fragn, basisn)])
            if vmfc:
                raise Exception("VMFC for hessian NYI")

            # Build up start and end slices
            abs_start, rel_start = 0, 0
            abs_slices, rel_slices = [], []
            for bas in basisn:
                rel_end = rel_start + 3 * fragment_size_dict[bas]
                rel_slices.append(slice(rel_start, rel_end))
                rel_start += 3 * fragment_size_dict[bas]

                tmp_slice = fragment_slice_dict[bas]
                abs_slices.append(slice(tmp_slice.start * 3, tmp_slice.end * 3))

            for abs_sl1, rel_sl1 in zip(abs_slices, rel_slices):
                for abs_sl2, rel_sl2 in zip(abs_slices, rel_slices):
                    ret[abs_sl1, abs_sl2] += hess[rel_sl1, rel_sl2]

    else:
        raise KeyError("ptype can only be gradient or hessian How did you end up here?")

def _print_nbody_energy(energy_body_dict, header):
        core.print_out("""\n   ==> N-Body: %s  energies <==\n\n""" % header)
        core.print_out("""   n-Body     Total Energy [Eh]       I.E. [kcal/mol]      Delta [kcal/mol]\n""")
        previous_e = energy_body_dict[1]
        nbody_range = list(energy_body_dict)
        nbody_range.sort()
        for n in nbody_range:
            delta_e = (energy_body_dict[n] - previous_e)
            delta_e_kcal = delta_e * constants.hartree2kcalmol
            int_e_kcal = (energy_body_dict[n] - energy_body_dict[1]) * constants.hartree2kcalmol
            core.print_out("""     %4s  %20.12f  %20.12f  %20.12f\n""" %
                                        (n, energy_body_dict[n], int_e_kcal, delta_e_kcal))
            previous_e = energy_body_dict[n]
        core.print_out("\n")

def nbody_gufunc(func, method_string, **kwargs):
    """
    Computes the nbody interaction energy, gradient, or Hessian depending on input.
    This is a generalized univeral function for computing interaction quantities.

    :returns: *return type of func* |w--w| The interaction data.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| interaction data and wavefunction with energy/gradient/hessian set appropriately when **return_wfn** specified.

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
    """

    ### ==> Parse some kwargs <==
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    ptype = kwargs.pop('ptype', None)
    return_total_data = kwargs.pop('return_total_data', False)
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()
    core.clean_variables()

    if ptype not in ['energy', 'gradient', 'hessian']:
        raise ValidationError("""N-Body driver: The ptype '%s' is not regonized.""" % ptype)

    # Figure out BSSE types
    do_cp = False
    do_nocp = False
    do_vmfc = False
    return_method = False

    # Must be passed bsse_type
    bsse_type_list = kwargs.pop('bsse_type')
    if bsse_type_list is None:
        raise ValidationError("N-Body GUFunc: Must pass a bsse_type")
    if not isinstance(bsse_type_list, list):
        bsse_type_list = [bsse_type_list]

    for num, btype in enumerate(bsse_type_list):
        if btype.lower() == 'cp':
            do_cp = True
            if (num == 0): return_method = 'cp'
        elif btype.lower() == 'nocp':
            do_nocp = True
            if (num == 0): return_method = 'nocp'
        elif btype.lower() == 'vmfc':
            do_vmfc = True
            if (num == 0): return_method = 'vmfc'
        else:
            raise ValidationError("N-Body GUFunc: bsse_type '%s' is not recognized" % btype.lower())

    max_nbody = kwargs.get('max_nbody', -1)
    max_frag = molecule.nfragments()
    if max_nbody == -1:
        max_nbody = molecule.nfragments()
    else:
        max_nbody = min(max_nbody, max_frag)

    # What levels do we need?
    nbody_range = range(1, max_nbody + 1)
    fragment_range = range(1, max_frag + 1)

    # Flip this off for now, needs more testing
    # If we are doing CP lets save them integrals
    #if 'cp' in bsse_type_list and (len(bsse_type_list) == 1):
    #    # Set to save RI integrals for repeated full-basis computations
    #    ri_ints_io = core.get_global_option('DF_INTS_IO')

    #    # inquire if above at all applies to dfmp2 or just scf
    #    core.set_global_option('DF_INTS_IO', 'SAVE')
    #    psioh = core.IOManager.shared_object()
    #    psioh.set_specific_retention(97, True)


    bsse_str = bsse_type_list[0]
    if len(bsse_type_list) >1:
        bsse_str =  str(bsse_type_list)
    core.print_out("\n\n")
    core.print_out("   ===> N-Body Interaction Abacus <===\n")
    core.print_out("        BSSE Treatment:                     %s\n" % bsse_str)


    cp_compute_list = {x:set() for x in nbody_range}
    nocp_compute_list = {x:set() for x in nbody_range}
    vmfc_compute_list = {x:set() for x in nbody_range}
    vmfc_level_list = {x:set() for x in nbody_range} # Need to sum something slightly different

    # Build up compute sets
    if do_cp:
        # Everything is in dimer basis
        basis_tuple = tuple(fragment_range)
        for nbody in nbody_range:
            for x in it.combinations(fragment_range, nbody):
                cp_compute_list[nbody].add( (x, basis_tuple) )

    if do_nocp:
        # Everything in monomer basis
        for nbody in nbody_range:
            for x in it.combinations(fragment_range, nbody):
                nocp_compute_list[nbody].add( (x, x) )

    if do_vmfc:
        # Like a CP for all combinations of pairs or greater
        for nbody in nbody_range:
            for cp_combos in it.combinations(fragment_range, nbody):
                basis_tuple = tuple(cp_combos)
                for interior_nbody in nbody_range:
                    for x in it.combinations(cp_combos, interior_nbody):
                        combo_tuple = (x, basis_tuple)
                        vmfc_compute_list[interior_nbody].add( combo_tuple )
                        vmfc_level_list[len(basis_tuple)].add( combo_tuple )

    # Build a comprehensive compute_range
    compute_list = {x:set() for x in nbody_range}
    for n in nbody_range:
        compute_list[n] |= cp_compute_list[n]
        compute_list[n] |= nocp_compute_list[n]
        compute_list[n] |= vmfc_compute_list[n]
        core.print_out("        Number of %d-body computations:     %d\n" % (n, len(compute_list[n])))


    # Build size and slices dictionaries
    fragment_size_dict = {frag: molecule.extract_subsets(frag).natom() for
                                           frag in range(1, max_frag+1)}

    start = 0
    fragment_slice_dict = {}
    for k, v in fragment_size_dict.items():
        fragment_slice_dict[k] = slice(start, start + v)
        start += v

    molecule_total_atoms = sum(fragment_size_dict.values())

    # Now compute the energies
    energies_dict = {}
    ptype_dict = {}
    for n in compute_list.keys():
        core.print_out("\n   ==> N-Body: Now computing %d-body complexes <==\n\n" % n)
        total = len(compute_list[n])
        for num, pair in enumerate(compute_list[n]):
            core.print_out("\n       N-Body: Computing complex (%d/%d) with fragments %s in the basis of fragments %s.\n\n" %
                                                                    (num + 1, total, str(pair[0]), str(pair[1])))
            ghost = list(set(pair[1]) - set(pair[0]))

            current_mol = molecule.extract_subsets(list(pair[0]), ghost)
            ptype_dict[pair] = func(method_string, molecule=current_mol, **kwargs)
            energies_dict[pair] = core.get_variable("CURRENT ENERGY")
            core.print_out("\n       N-Body: Complex Energy (fragments = %s, basis = %s: %20.14f)\n" %
                                                                (str(pair[0]), str(pair[1]), energies_dict[pair]))

            # Flip this off for now, needs more testing
            #if 'cp' in bsse_type_list and (len(bsse_type_list) == 1):
            #    core.set_global_option('DF_INTS_IO', 'LOAD')

            core.clean()

    # Final dictionaries
    cp_energy_by_level   = {n: 0.0 for n in nbody_range}
    nocp_energy_by_level = {n: 0.0 for n in nbody_range}

    cp_energy_body_dict =   {n: 0.0 for n in nbody_range}
    nocp_energy_body_dict = {n: 0.0 for n in nbody_range}
    vmfc_energy_body_dict = {n: 0.0 for n in nbody_range}

    # Build out ptype dictionaries if needed
    if ptype != 'energy':
        if ptype == 'gradient':
            arr_shape = (molecule_total_atoms, 3)
        elif ptype == 'hessian':
            arr_shape = (molecule_total_atoms * 3, molecule_total_atoms * 3)
        else:
            raise KeyError("N-Body: ptype '%s' not recognized" % ptype)

        cp_ptype_by_level   =  {n: np.zeros(arr_shape) for n in nbody_range}
        nocp_ptype_by_level =  {n: np.zeros(arr_shape) for n in nbody_range}
        vmfc_ptype_by_level = {n: np.zeros(arr_shape) for n in nbody_range}

        cp_ptype_body_dict   = {n: np.zeros(arr_shape) for n in nbody_range}
        nocp_ptype_body_dict = {n: np.zeros(arr_shape) for n in nbody_range}
        vmfc_ptype_body_dict = {n: np.zeros(arr_shape) for n in nbody_range}
    else:
        cp_ptype_by_level, cp_ptype_body_dict = None, None
        nocp_ptype_by_level, nocp_ptype_body_dict = None, None
        vmfc_ptype_body_dict = None


    # Sum up all of the levels
    for n in nbody_range:

        # Energy
        cp_energy_by_level[n]   = sum(energies_dict[v] for v in cp_compute_list[n])
        nocp_energy_by_level[n] = sum(energies_dict[v] for v in nocp_compute_list[n])

        # Special vmfc case
        if n > 1:
            vmfc_energy_body_dict[n] = vmfc_energy_body_dict[n - 1]
        for tup in vmfc_level_list[n]:
            vmfc_energy_body_dict[n] += ((-1) ** (n - len(tup[0]))) * energies_dict[tup]


        # Do ptype
        if ptype != 'energy':
            _sum_cluster_ptype_data(ptype, ptype_dict, cp_compute_list[n],
                                      fragment_slice_dict, fragment_size_dict,
                                      cp_ptype_by_level[n])
            _sum_cluster_ptype_data(ptype, ptype_dict, nocp_compute_list[n],
                                      fragment_slice_dict, fragment_size_dict,
                                      nocp_ptype_by_level[n])
            _sum_cluster_ptype_data(ptype, ptype_dict, vmfc_level_list[n],
                                      fragment_slice_dict, fragment_size_dict,
                                      vmfc_ptype_by_level[n], vmfc=True)

    # Compute cp energy and ptype
    if do_cp:
        for n in nbody_range:
            if n == max_frag:
                cp_energy_body_dict[n] = cp_energy_by_level[n]
                if ptype != 'energy':
                    cp_ptype_body_dict[n][:] = cp_ptype_by_level[n]
                continue

            for k in range(1, n + 1):
                take_nk =  nCr(max_frag - k - 1, n - k)
                sign = ((-1) ** (n - k))
                value = cp_energy_by_level[k]
                cp_energy_body_dict[n] += take_nk * sign * value

                if ptype != 'energy':
                    value = cp_ptype_by_level[k]
                    cp_ptype_body_dict[n] += take_nk * sign * value

        _print_nbody_energy(cp_energy_body_dict, "Counterpoise Corrected (CP)")
        cp_interaction_energy = cp_energy_body_dict[max_nbody] - cp_energy_body_dict[1]
        core.set_variable('Counterpoise Corrected Total Energy', cp_energy_body_dict[max_nbody])
        core.set_variable('Counterpoise Corrected Interaction Energy', cp_interaction_energy)

        for n in nbody_range[1:]:
            var_key = 'CP-CORRECTED %d-BODY INTERACTION ENERGY' % n
            core.set_variable(var_key, cp_energy_body_dict[n] - cp_energy_body_dict[1])

    # Compute nocp energy and ptype
    if do_nocp:
        for n in nbody_range:
            if n == max_frag:
                nocp_energy_body_dict[n] = nocp_energy_by_level[n]
                if ptype != 'energy':
                    nocp_ptype_body_dict[n][:] = nocp_ptype_by_level[n]
                continue

            for k in range(1, n + 1):
                take_nk =  nCr(max_frag - k - 1, n - k)
                sign = ((-1) ** (n - k))
                value = nocp_energy_by_level[k]
                nocp_energy_body_dict[n] += take_nk * sign * value

                if ptype != 'energy':
                    value = nocp_ptype_by_level[k]
                    nocp_ptype_body_dict[n] += take_nk * sign * value

        _print_nbody_energy(nocp_energy_body_dict, "Non-Counterpoise Corrected (NoCP)")
        nocp_interaction_energy = nocp_energy_body_dict[max_nbody] - nocp_energy_body_dict[1]
        core.set_variable('Non-Counterpoise Corrected Total Energy', nocp_energy_body_dict[max_nbody])
        core.set_variable('Non-Counterpoise Corrected Interaction Energy', nocp_interaction_energy)

        for n in nbody_range[1:]:
            var_key = 'NOCP-CORRECTED %d-BODY INTERACTION ENERGY' % n
            core.set_variable(var_key, nocp_energy_body_dict[n] - nocp_energy_body_dict[1])


    # Compute vmfc energy and ptype
    if do_vmfc:
        _print_nbody_energy(vmfc_energy_body_dict, "Valiron-Mayer Function Couterpoise (VMFC)")
        vmfc_interaction_energy = vmfc_energy_body_dict[max_nbody] - vmfc_energy_body_dict[1]
        core.set_variable('Valiron-Mayer Function Couterpoise Total Energy', vmfc_energy_body_dict[max_nbody])
        core.set_variable('Valiron-Mayer Function Couterpoise Interaction Energy', vmfc_interaction_energy)

        for n in nbody_range[1:]:
            var_key = 'VMFC-CORRECTED %d-BODY INTERACTION ENERGY' % n
            core.set_variable(var_key, vmfc_energy_body_dict[n] - vmfc_energy_body_dict[1])

    if return_method == 'cp':
        ptype_body_dict = cp_ptype_body_dict
        energy_body_dict = cp_energy_body_dict
    elif return_method == 'nocp':
        ptype_body_dict = nocp_ptype_body_dict
        energy_body_dict = nocp_energy_body_dict
    elif return_method == 'vmfc':
        ptype_body_dict = vmfc_ptype_body_dict
        energy_body_dict = vmfc_energy_body_dict
    else:
        raise ValidationError("N-Body Wrapper: Invalid return type. Should never be here, please post this error on github.")


    # Figure out and build return types
    if return_total_data:
        ret_energy = energy_body_dict[max_nbody]
    else:
        ret_energy = energy_body_dict[max_nbody]
        ret_energy -= energy_body_dict[1]


    if ptype != 'energy':
        if return_total_data:
            np_final_ptype = ptype_body_dict[max_nbody].copy()
        else:
            np_final_ptype = ptype_body_dict[max_nbody].copy()
            np_final_ptype -= ptype_body_dict[1]

            ret_ptype = core.Matrix.from_array(np_final_ptype)
    else:
        ret_ptype = ret_energy

    # Build and set a wavefunction
    wfn = core.Wavefunction.build(molecule, 'sto-3g')
    wfn.nbody_energy = energies_dict
    wfn.nbody_ptype = ptype_dict
    wfn.nbody_body_energy = energy_body_dict
    wfn.nbody_body_ptype = ptype_body_dict

    if ptype == 'gradient':
        wfn.set_gradient(ret_ptype)
    elif ptype == 'hessian':
        wfn.set_hessian(ret_ptype)

    core.set_variable("CURRENT ENERGY", ret_energy)

    if return_wfn:
        return (ret_ptype, wfn)
    else:
        return ret_ptype


