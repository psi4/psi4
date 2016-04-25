#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
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
#@END LICENSE
#

from __future__ import print_function
from __future__ import absolute_import

import psi4

import numpy as np
import itertools as it

# Import driver helpers
import driver_util
import p4util

from qmmm import *
from procedures import *
from p4util.exceptions import *

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
        print("""\n   ==> N-Body: %s  energies <==""" % header)
        print("""   n-Body        Total Energy [Eh]          I.E. [kcal/mol]         Delta [kcal/mol]""")
        previous_e = energy_body_dict[1]
        nbody_range = energy_body_dict.keys()
        nbody_range.sort()
        for n in nbody_range:
            delta_e = (energy_body_dict[n] - previous_e)
            delta_e_kcal = delta_e * p4const.psi_hartree2kcalmol
            int_e_kcal = (energy_body_dict[n] - energy_body_dict[1]) * p4const.psi_hartree2kcalmol
            print("""     %4s        % 16.14f        % 16.14f        % 16.14f""" %
                                        (n, energy_body_dict[n], int_e_kcal, delta_e_kcal))
            previous_e = energy_body_dict[n]
        print("\n")

def nCr(n, r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

def _nbody_gufunc(func, method_string, **kwargs):
    """
    Computes the nbody interaction energy, gradient, or Hessian depending on input.

    Parameters
    ----------
    func : python function
        Python function that accepts method_string and a molecule and returns a energy, gradient, or Hessian.
    method_string : str
        Lowername to be passed to function
    molecule : psi4.Molecule (default: Global Molecule)
        Molecule to use in all computations
    return_wfn : bool (default: False)
        Return a wavefunction or not
    bsse_type : str or list (default: None, this function is not called)
        Type of BSSE correction to compute: CP, NoCP, or VMFC. The first in this list is returned by this function.
    max_nbody : int
        Maximum n-body to compute, cannot exceede the number of fragments in the moleucle

    Returns
    -------
    data : return type of func
        The interaction data
    wfn : psi4.Wavefunction (optional)
        A wavefunction with... something in it?

    Notes
    -----
    This is a generalized univeral function for compute interaction quantities.

    Examples
    --------
    """

    ### ==> Parse some kwargs <==
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    ptype = kwars.pop('ptype', None)
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()
    psi4.clean_variables()

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

    # If we are doing CP lets save them integrals
    if 'cp' in bsse_type_list and (len(bsse_type_list) == 1):
        # Set to save RI integrals for repeated full-basis computations
        ri_ints_io = psi4.get_global_option('DF_INTS_IO')

        # inquire if above at all applies to dfmp2 or just scf
        psi4.set_global_option('DF_INTS_IO', 'SAVE')
        psioh = psi4.IOManager.shared_object()
        psioh.set_specific_retention(97, True)


    print("\n\n")
    print("   ===> N-Body Interaction Abacus <===\n")
    print("        BSSE Treatment:          %s")


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
        print("        Number of %d-body computations:          %d" % (n, len(compute_list[n])))


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
        print("\n   ==> N-Body: Now computing %d-body complexes <==\n" % n)
        total = len(compute_list[n])
        for num, pair in enumerate(compute_list[n]):
            print("\n       N-Body: Computing complex (%d/%d) with fragments %s in the basis of fragments %s.\n" %
                                                                    (num + 1, total, str(pair[0]), str(pair[1])))
            ghost = list(set(pair[1]) - set(pair[0]))
            current_mol = molecule.extract_subsets(list(pair[0]), ghost)
            ptype_dict[pair] = func(method_string, molecule=current_mol, **kwargs)
            energies_dict[pair] = psi4.get_variable("CURRENT ENERGY")
            psi4.clean()

    # Figure out if we are dealing with a matrix or an array
    ptype = None
    try:
        shape = ptype_dict[ptype_dict.keys()[0]].__array_interface__['shape']
        if len(shape) > 2:
            raise ValueError("N-Body: Return type of passed function has too many dimensions!")
        elif shape[1] == 3:
            ptype = 'gradient'
        elif shape[0] == shape[1]:
            ptype = 'hessian'
        else:
            raise ValueError("N-Body: Return type of passed function has incorrect dimensions.")
    except:
        if not isinstance(ptype_dict[ptype_dict.keys()[0]], float):
            raise ValueError("N-Body: Return type of passed function is not understood: '%s'." %
                             type(ptype_dict[ptype_dict.keys()[0]]))
        ptype = 'energy'


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

        cp_ptype_body_dict   = {n: np.zeros(arr_shape) for n in nbody_range}
        nocp_ptype_body_dict = {n: np.zeros(arr_shape) for n in nbody_range}
        vmfc_ptype_body_dict = {n: np.zeros(arr_shape) for n in nbody_range}
    else:
        cp_ptype_by_level, cp_ptype_body_dict = None, None
        nocp_ptype_by_level, nocp_ptype_body_dict = None, None
        vmfc_ptype_by_level= None


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
        cp_final_energy = cp_energy_body_dict[n] - cp_energy_body_dict[1]
        psi4.set_variable('Counterpoise Corrected Interaction Energy', cp_final_energy)

        for n in nbody_range[1:]:
            var_key = 'CP-CORRECTED %d-BODY INTERACTION ENERGY' % n
            psi4.set_variable(var_key, cp_energy_body_dict[n] - cp_energy_body_dict[1])

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
        nocp_final_energy = nocp_energy_body_dict[n] - nocp_energy_body_dict[1]
        psi4.set_variable('Non-Counterpoise Corrected Interaction Energy', nocp_final_energy)

        for n in nbody_range[1:]:
            var_key = 'NOCP-CORRECTED %d-BODY INTERACTION ENERGY' % n
            psi4.set_variable(var_key, nocp_energy_body_dict[n] - nocp_energy_body_dict[1])


    # Compute vmfc energy and ptype
    if do_vmfc:
        _print_nbody_energy(vmfc_energy_body_dict, "Valiron-Mayer Function Couterpoise (VMFC)")
        vmfc_final_energy = vmfc_energy_body_dict[n] - vmfc_energy_body_dict[1]
        psi4.set_variable('Valiron-Mayer Function Couterpoise Interaction Energy', vmfc_final_energy)

        for n in nbody_range[1:]:
            var_key = 'VMFC-CORRECTED %d-BODY INTERACTION ENERGY' % n
            psi4.set_variable(var_key, vmfc_energy_body_dict[n] - vmfc_energy_body_dict[1])


    if return_method == 'cp':
        ret_energy = cp_final_energy
        ptype_body_dict = cp_ptype_body_dict
        total_energy = cp_energy_body_dict[-1]
    elif return_method == 'nocp':
        ret_energy = nocp_final_energy
        ptype_body_dict = nocp_ptype_body_dict
        total_energy = nocp_energy_body_dict[-1]
    elif return_method == 'vmfc':
        ret_energy = vmfc_final_energy
        ptype_body_dict = vmfc_ptype_body_dict
        total_energy = vmfc_energy_body_dict[-1]
    else:
        raise ValidationError("N-Body Wrapper: Invalid return type. Should never be here, please post this error on github.")


    if ptype != 'energy':
        np_final_ptype = ptype_body_dict[max_nbody].copy()
        np_final_ptype -= ptype_body_dict[1]

        ret_ptype = psi4.Matrix(*np_cp_final_ptype.shape)
        ret_ptype_view = np.asarray(final_ptype)
        ret_ptype_view[:] = np_final_ptype

    # Build a wavefunction
    wfn = psi4.new_wavefunction(molecule, 'sto-3g')
    wfn.energy_dict = energies_dict
    wfn.ptype_dict = ptype_dict
    if ptype == 'gradient':
        wfn.set_gradient(ret_ptype)
    elif ptype == 'hessian':
        wfn.set_hessian(ret_ptype)

    psi4.set_variable("CURRENT ENERGY", ret_energy)

    if return_wfn:
        return (ptype_value, wfn)
    else:
        return ptype_value


