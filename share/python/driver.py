
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
"""Module with a *procedures* dictionary specifying available quantum
chemical methods and functions driving the main quantum chemical
functionality, namely single-point energies, geometry optimizations,
properties, and vibrational frequency calculations.

"""
from __future__ import absolute_import
import sys
import re
from proc import *
from interface_cfour import *
import driver_util
import numpy as np
import itertools as it
import math
# never import wrappers or aliases into this file



# Procedure lookup tables
procedures = {
        'energy': {
            'hf'            : run_scf,
            'scf'           : run_scf,
            'mcscf'         : run_mcscf,
            'dcft'          : run_dcft,
            'mp3'           : select_mp3,
            'mp2.5'         : select_mp2p5,
            'mp2'           : select_mp2,
            'omp2'          : select_omp2,
            'scs-omp2'      : run_occ,
            'scs(n)-omp2'   : run_occ,
            'scs-omp2-vdw'  : run_occ,
            'sos-omp2'      : run_occ,
            'sos-pi-omp2'   : run_occ,
            'omp3'          : select_omp3,
            'scs-omp3'      : run_occ,
            'scs(n)-omp3'   : run_occ,
            'scs-omp3-vdw'  : run_occ,
            'sos-omp3'      : run_occ,
            'sos-pi-omp3'   : run_occ,
            'olccd'         : select_olccd,
            'omp2.5'        : select_omp2p5,
            'dfocc'         : run_dfocc,
            'qchf'          : run_qchf,
            'ccd'           : run_dfocc,
            'sapt0'         : run_sapt,
            'sapt2'         : run_sapt,
            'sapt2+'        : run_sapt,
            'sapt2+(3)'     : run_sapt,
            'sapt2+3'       : run_sapt,
            'sapt2+(ccd)'   : run_sapt,
            'sapt2+(3)(ccd)': run_sapt,
            'sapt2+3(ccd)'  : run_sapt,
            'sapt2+dmp2'    : run_sapt,
            'sapt2+(3)dmp2' : run_sapt,
            'sapt2+3dmp2'   : run_sapt,
            'sapt2+(ccd)dmp2' : run_sapt,
            'sapt2+(3)(ccd)dmp2' : run_sapt,
            'sapt2+3(ccd)dmp2' : run_sapt,
            'sapt0-ct'      : run_sapt_ct,
            'sapt2-ct'      : run_sapt_ct,
            'sapt2+-ct'     : run_sapt_ct,
            'sapt2+(3)-ct'  : run_sapt_ct,
            'sapt2+3-ct'    : run_sapt_ct,
            'sapt2+(ccd)-ct'     : run_sapt_ct,
            'sapt2+(3)(ccd)-ct'  : run_sapt_ct,
            'sapt2+3(ccd)-ct'    : run_sapt_ct,
            'fisapt0'       : run_fisapt,
            'ccenergy'      : run_ccenergy,  # full control over ccenergy
            'ccsd'          : select_ccsd,
            'ccsd(t)'       : select_ccsd_t_,
            'ccsd(at)'      : select_ccsd_at_,
            'cc2'           : run_ccenergy,
            'cc3'           : run_ccenergy,
            'mrcc'          : run_mrcc,  # interface to Kallay's MRCC program
            'bccd'          : run_bccd,
            'bccd(t)'       : run_bccd,
            'eom-ccsd'      : run_eom_cc,
            'eom-cc2'       : run_eom_cc,
            'eom-cc3'       : run_eom_cc,
            'detci'         : run_detci,  # full control over detci
            'mp'            : run_detci,  # arbitrary order mp(n)
            'zapt'          : run_detci,  # arbitrary order zapt(n)
            'cisd'          : select_cisd,
            'cisdt'         : run_detci,
            'cisdtq'        : run_detci,
            'ci'            : run_detci,  # arbitrary order ci(n)
            'fci'           : run_detci,
            'casscf'        : run_detcas,
            'rasscf'        : run_detcas,
            'adc'           : run_adc,
#            'cphf'          : run_libfock,
#            'cis'           : run_libfock,
#            'tdhf'          : run_libfock,
#            'cpks'          : run_libfock,
#            'tda'           : run_libfock,
#            'tddft'         : run_libfock,
            'psimrcc'       : run_psimrcc,
            'psimrcc_scf'   : run_psimrcc_scf,
            'qcisd'         : run_fnocc,
            'qcisd(t)'      : run_fnocc,
            'mp4'           : select_mp4,
            'mp4(sdq)'      : run_fnocc,
            'fno-ccsd'      : select_fnoccsd,
            'fno-ccsd(t)'   : select_fnoccsd_t_,
            'fno-qcisd'     : run_fnocc,
            'fno-qcisd(t)'  : run_fnocc,
            'fno-mp3'       : run_fnocc,
            'fno-mp4(sdq)'  : run_fnocc,
            'fno-mp4'       : run_fnocc,
            'fno-lccd'      : run_cepa,
            'fno-lccsd'     : run_cepa,
            'fno-cepa(0)'   : run_cepa,
            'fno-cepa(1)'   : run_cepa,
            'fno-cepa(3)'   : run_cepa,
            'fno-acpf'      : run_cepa,
            'fno-aqcc'      : run_cepa,
            'fno-cisd'      : run_cepa,
            'lccd'          : select_lccd,
            'lccsd'         : run_cepa,
            'cepa(0)'       : run_cepa,
            'cepa(1)'       : run_cepa,
            'cepa(3)'       : run_cepa,
            'acpf'          : run_cepa,
            'aqcc'          : run_cepa,
            'efp'           : run_efp,
            'dmrgscf'       : run_dmrgscf,
            'dmrgci'        : run_dmrgci,
            # Upon adding a method to this list, add it to the docstring in energy() below
            # Aliases are discouraged. If you must add an alias to this list (e.g.,
            #    lccsd/cepa(0)), please search the whole driver to find uses of
            #    name in return values and psi variables and extend the logic to
            #    encompass the new alias.
        },
        'gradient' : {
            'hf'            : run_scf_gradient,
            'scf'           : run_scf_gradient,
            'ccsd'          : select_ccsd_gradient,
            'ccsd(t)'       : select_ccsd_t__gradient,
            'mp2'           : select_mp2_gradient,
            'eom-ccsd'      : run_eom_cc_gradient,
            'dcft'          : run_dcft_gradient,
            'omp2'          : select_omp2_gradient,
            'omp3'          : select_omp3_gradient,
            'mp3'           : select_mp3_gradient,
            'mp2.5'         : select_mp2p5_gradient,
            'omp2.5'        : select_omp2p5_gradient,
            'lccd'          : select_lccd_gradient,
            'olccd'         : select_olccd_gradient,
            'ccd'           : run_dfocc_gradient,
            # Upon adding a method to this list, add it to the docstring in optimize() below
        },
        'hessian' : {
            # Upon adding a method to this list, add it to the docstring in frequency() below
        },
        'property' : {
            'hf'       : run_scf_property,
            'scf'      : run_scf_property,
            'mp2'      : select_mp2_property,
            'cc2'      : run_cc_property,
            'ccsd'     : run_cc_property,
            'eom-cc2'  : run_cc_property,
            'eom-ccsd' : run_cc_property,
            'detci'    : run_detci_property,  # full control over detci
            'cisd'     : run_detci_property,
            'cisdt'    : run_detci_property,
            'cisdtq'   : run_detci_property,
            'ci'       : run_detci_property,  # arbitrary order ci(n)
            'fci'      : run_detci_property,
            'rasscf'   : run_detci_property,
            'casscf'   : run_detci_property,
            # Upon adding a method to this list, add it to the docstring in property() below
        }}

# Will only allow energy to be run for the following methods
energy_only_methods = [x for x in procedures['energy'].keys() if 'sapt' in x]
energy_only_methods += ['adc', 'efp', 'cphf', 'tdhf', 'cis']

# dictionary to register pre- and post-compute hooks for driver routines
hooks = dict((k1, dict((k2, []) for k2 in ['pre', 'post'])) for k1 in ['energy', 'optimize', 'frequency'])

# Integrate DFT with driver routines
for ssuper in superfunctional_list():
    procedures['energy'][ssuper.name().lower()] = run_dft
    if not ssuper.is_c_hybrid():
        procedures['property'][ssuper.name().lower()] = run_dft_property

for ssuper in superfunctional_list():
    if ((not ssuper.is_c_hybrid()) and (not ssuper.is_c_lrc()) and (not ssuper.is_x_lrc())):
        procedures['gradient'][ssuper.name().lower()] = run_dft_gradient

# Integrate CFOUR with driver routines
for ssuper in cfour_list():
    procedures['energy'][ssuper.lower()] = run_cfour

for ssuper in cfour_gradient_list():
    procedures['gradient'][ssuper.lower()] = run_cfour

### Helper functions

def _method_exists(ptype, method_name):
    r"""
    Quick check to see if this method exists, if it does not exist we raise a convenient flag.
    """
    if method_name not in procedures[ptype].keys():
        alternatives = ""
        alt_method_name = p4util.text.find_approximate_string_matches(method_name,
                                                                procedures[ptype].keys(), 2)
        if len(alt_method_name) > 0:
            alternatives = " Did you mean? %s" % (" ".join(alt_method_name))
        Cptype = ptype[0].upper() + ptype[1:]
        raise ValidationError('%s method "%s" is not available.%s' % (Cptype, method_name, alternatives))

def _set_convergence_criterion(ptype, method_name, scf_Ec, pscf_Ec, scf_Dc, pscf_Dc, gen_Ec):
    r"""

    This function will set local SCF and global energy convergence criterion
    to the defaults listed at:
    http://www.psicode.org/psi4manual/master/scf.html#convergence-and-
    algorithm-defaults. SCF will be convergence more tightly if a post-SCF
    method is select (pscf_Ec, and pscf_Dc) else the looser (scf_Ec, and
    scf_Dc convergence criterion will be used).

    ptype -         Procedure type (energy, gradient, etc)
    method_name -   Name of the method
    scf_Ec -        SCF E convergence criterion
    pscf_Ec -       Post-SCF E convergence criterion
    scf_Dc -        SCF D convergence criterion
    pscf_Dc -       Post-SCF D convergence criterion
    gen_Ec -        General E convergence for all methods
    """

    # Kind of want to move this out of here
    _method_exists(ptype, method_name)

    # Set method-dependent scf convergence criteria, check against energy routines
    if not psi4.has_option_changed('SCF', 'E_CONVERGENCE'):
        if procedures['energy'][method_name] in [run_scf, run_dft]:
            psi4.set_local_option('SCF', 'E_CONVERGENCE', scf_Ec)
        else:
            psi4.set_local_option('SCF', 'E_CONVERGENCE', pscf_Ec)

    if not psi4.has_option_changed('SCF', 'D_CONVERGENCE'):
        if procedures['energy'][method_name] in [run_scf, run_dft]:
            psi4.set_local_option('SCF', 'D_CONVERGENCE', scf_Dc)
        else:
            psi4.set_local_option('SCF', 'D_CONVERGENCE', pscf_Dc)

    # Set post-scf convergence criteria (global will cover all correlated modules)
    if not psi4.has_global_option_changed('E_CONVERGENCE'):
        if procedures['energy'][method_name] not in [run_scf, run_dft]:
            psi4.set_global_option('E_CONVERGENCE', gen_Ec)

    # We passed!
    return True


def _find_derivative_type(ptype, method_name, user_dertype):
    r"""
    Figures out the derivate type (0, 1, 2) for a given method_name. Will
    first use user default and then the highest available derivative type for
    a given method.
    """

    if ptype not in ['gradient', 'hessian']:
        raise KeyError("_find_derivative_type: ptype must either be gradient or hessian.")

    dertype = "(auto)"

    # If user type is None, try to find the highest derivative
    if user_dertype is None:
        if (ptype == 'hessian') and (method_name in procedures['hessian']):
            dertype = 2
            # Will need special logic if we ever have managed Hessians
        elif method_name in procedures['gradient']:
            dertype = 1
            if procedures['gradient'][method_name].__name__.startswith('select_'):
                try:
                    procedures['gradient'][method_name](method_name, probe=True)
                except ManagedMethodError:
                    dertype = 0
        elif method_name in procedures['energy']:
            dertype = 0
    else:
        # Quick sanity check. Only *should* be able to be None or int, but hey, kids today...
        if not isinstance(user_dertype, int):
            raise TypeError("_find_derivative_type: user_dertype should only be None or int!")
        dertype = user_dertype

    # Summary validation
    if (dertype == 2) and (method_name in procedures['hessian']):
        method = hessian
    elif (dertype == 1) and (method_name in procedures['gradient']):
        method = gradient
    elif (dertype == 0) and (method_name in procedures['energy']):
        method = energy
    else:
        alternatives = ''
        alt_method_name = p4util.text.find_approximate_string_matches(method_name, procedures['energy'].keys(), 2)
        if len(alt_method_name) > 0:
            alternatives = """ Did you mean? %s""" % (' '.join(alt_method_name))

        raise ValidationError("""Derivative method 'name' %s and derivative level 'dertype' %s are not available.%s"""
            % (method_name, str(dertype), alternatives))

    return (dertype, method)

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
    psi4.clean_variables()

    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()
#    moleucle.

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


# End CBS gufunc data

def _cbs_gufunc(ptype, total_method_name, **kwargs):

    # Catch kwarg issues
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    psi4.clean_variables()
    user_dertype = kwargs.pop('dertype', None)
    cbs_verbose = kwargs.pop('cbs_verbose', True)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # Same some global variables so we can reset them later
    optstash = p4util.OptionsState(
        ['BASIS'],
    )

    # Sanitize total_method_name
    total_method_name = total_method_name.lower()
    replacement_list = [(' ', ''), ('D-', 'D:'), ('{', '['), ('}', ']')]
    for char, rep in replacement_list:
        total_method_name = total_method_name.replace(char, rep)

    # Generate list of list [[method_name, [bas1, bas2, ...], ...]
    method_list = []

    # Split into components
    total_method_name_list = re.split(r'/\+(?!([^\(]*\)|[^\[]*\]))$', total_method_name)
    for method_str in total_method_name_list:
        if (method_str.count("[") > 1) or (method_str.count("]") > 1):
            raise ValidationError("""CBS gufunc: Too many brakcets given! %s """ % method_str)

        if method_str.count('/') != 1:
            raise ValidationError("""CBS gufunc: All methods must specify a basis with '/'. %s""" % method_str)

        # Find method and basis
        method, basis_str = method_str.split('/')
        isDelta = False
        if 'D:' in method:
            if 'D:' in method[2:]:
                raise ValidationError("""CBS gufunc: Invalid Delta position. %s""" % method)
            isDelta = True
            method = method[2:]

        # Expand basis set (aug-cc-pv[d,t]z -> [aug-cc-pvdz, aug-cc-pvtz])
        molstr = molecule.create_psi4_string_from_molecule()
        basissets, dunning_num = driver_util.expand_bracketed_basis(basis_str, molecule=molstr)

        # check method, need to be careful with gradients
        _method_exists(ptype, method)

        if (method in energy_only_methods) and (len(bassisets) > 2):
            raise ValidationError("CBS gufunc: Method '%s' cannot be extrapolated" % method)

        # Build up dictionary of information
        method_dict = {}
        method_dict['basissets'] = basissets
        method_dict['basis_zetas'] = dunning_num
        method_dict['method_name'] = method
        method_dict['isSCF'] = (procedures['energy'][method] in [run_scf, run_dft])
        method_dict['isDelta'] = isDelta
        method_dict['dertype'] = None
        if ptype != 'energy':
            method_dict['dertype'] =  _find_derivative_type(ptype, method, user_dertype)[0]

        # Figure out extrapolation method
        method_dict['extrapolation_type'] = None
        if method_dict['isSCF']:
            if len(basissets) == 1:
                method_dict['extrapolation_type'] = None
            elif len(basissets) == 2:
                method_dict['extrapolation_type'] = driver_util.scf_xtpl_helgaker_2
            elif len(basissets) == 3:
                method_dict['extrapolation_type'] = driver_util.scf_xtpl_helgaker_3
            else:
                raise ValidationError("CBS Gufunc: SCF-type method '%s' can only be supplied 1, 2, or"
                                      " 3 basis sets (given %d)." % (method, len(basissets)))
        else:
            if len(basissets) == 1:
                method_dict['extrapolation_type'] = None
            elif len(basissets) == 2:
                method_dict['extrapolation_type'] = driver_util.corl_xtpl_helgaker_2
            else:
                raise ValidationError("CBS Gufunc: post-SCF-type method '%s' can only be supplied 1 or 2"
                                      " basis sets (given %d)." % (method, len(basissets)))

        method_list.append(method_dict)

    # SCF/cc-pv[DTQ]Z + D:MP2/cc-pv[DT]Z + D:CCSD(T)
    # MP2/cc-pv[DT]Z + D:CCSD(T)

    # SCF/cc-pv[DT]Z
    # scf_basis = [...]
    # corl_wfn
    # corl_basis
    # delta_wfn
    # delta_basis 
    # delta2_wfn
    # delta2_basis 

    # Validate the method_list

    # Pick the correct functin to call
    if ptype == 'energy':
        pfunc = energy
    elif ptype == 'gradient':
        pfunc = gradient
    elif ptype == 'hessian':
        pfunc = hessian
    else:
        raise ValidtionError("""CBS gufunc: ptype '%s' is not recognized""" % ptype)

    # We also need energy for gradient and Hessian methods
    energy_list = []
    ptype_list = []


    ### DGAS compute loop
    # Loop over methods
    for method_dict in method_list:
        scf_ptype_list = []
        corr_ptype_list = []

        scf_energy_list = []
        corr_energy_list = []

        # Loop over basis sets
        method_name = method_dict['method_name']
        for basis in method_dict['basissets']:
            #print(ptype, method_name, basis)
            psi4.set_global_option('BASIS', basis)

            # Just need a single energy list
            if method_dict['isSCF'] or (method_dict['extrapolation_type'] is None):
                pvalue, pwfn = pfunc(method, molecule=molecule, return_wfn=True)
                scf_ptype_list.append(pvalue)
                scf_energy_list.append(psi4.get_variable('CURRENT ENERGY'))

            # Need to seperate into scf and corr values
            else:
                scf_pvalue, scf_wfn = pfunc('SCF', molecule=molecule, return_wfn=True)
                scf_energy = psi4.get_variable('CURRENT ENERGY')

                scf_ptype_list.append(scf_pvalue)
                scf_energy_list.append(scf_energy)

                pvalue, pwfn = pfunc(method, ref_wfn=scf_wfn, return_wfn=True)
                corr_energy = psi4.get_variable('CURRENT ENERGY')

                # Correlation only part
                if ptype == 'energy':
                    pvalue = pvalue - scf_pvalue
                else:
                    pvalue.subtract(scf_pvalue)

                corr_ptype_list.append(pvalue)
                corr_energy_list.append(corr_energy - scf_energy)

            psi4.clean()

        # Now build final values
        method_total_energy = 0
        method_ptype_data = None

        if method_dict['extrapolation_type'] is None:
            method_total_energy = scf_energy_list[0]
            method_ptype_data = scf_ptype_list[0]

        else:
            xt_type = method_dict['extrapolation_type']
            mt_name = method_dict['method_name']

            # Pick which list to extrapolate
            if method_dict['isSCF']:
                e_list = scf_energy_list
                p_list = scf_ptype_list
            else:
                e_list = corr_energy_list
                p_list = corr_ptype_list

            # Extrapolate the list for energy
            extrap_data = []
            for z, data in zip(method_dict['basis_zetas'], e_list):
                extrap_data.append(z)
                extrap_data.append(data)

            if method_dict['isSCF']:
                method_total_energy = xt_type(mt_name,
                                              *extrap_data, verbose=cbs_verbose)
            else:
                method_total_energy = xt_type(mt_name, scf_energy_list[-1],
                                              *extrap_data, verbose=cbs_verbose)

            # Extrapolate the list of ptype data
            if ptype == 'energy':
                # We already have energy, skip extrapolation
                method_ptype_data = method_total_energy
            else:
                extrap_data = []
                for z, data in zip(method_dict['basis_zetas'], p_list):
                    extrap_data.append(z)
                    extrap_data.append(data)

                if method_dict['isSCF']:
                    method_ptype_data = xt_type(mt_name,
                                                *extrap_data, verbose=cbs_verbose)
                else:
                    method_ptype_data = xt_type(mt_name, scf_ptype_list[-1],
                                                *extrap_data, verbose=cbs_verbose)


        # Append to total list
        energy_list.append(method_total_energy)
        ptype_list.append(method_ptype_data)


    total_energy = sum(energy_list)
    psi4.set_variable('CURRENT ENERGY', total_energy)
    if ptype == 'energy':
        ptype_value = sum(ptype_list)
    else:
        ptype_value = ptype_list[0].clone()
        for val in ptype_list[1:]:
            ptype_value.add(val)


    wfn = None

    # Reset modified global variables
    optstash.restore()

    if return_wfn:
        return (ptype_value, wfn)
    else:
        return ptype_value

def energy(name, **kwargs):
    r"""Function to compute the single-point electronic energy.

    :returns: *float* |w--w| Total electronic energy in Hartrees. SAPT & EFP return interaction energy.

    :returns: (*float*, :ref:`Wavefunction<sec:psimod_Wavefunction>`) |w--w| energy and wavefunction when **return_wfn** specified.

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`CURRENT ENERGY <CURRENTENERGY>`
       * :psivar:`CURRENT REFERENCE ENERGY <CURRENTREFERENCEENERGY>`
       * :psivar:`CURRENT CORRELATION ENERGY <CURRENTCORRELATIONENERGY>`

    :type name: string
    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type return_wfn: :ref:`boolean <op_py_boolean>`
    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|

        Indicate to additionally return the :ref:`Wavefunction<sec:psimod_Wavefunction>`
        calculation result as the second element (after *float* energy) of a tuple.

    :type restart_file: string
    :param restart_file: ``['file.1, file.32]`` || ``./file`` || etc.

        Binary data files to be renamed for calculation restart.


    .. _`table:energy_gen`:

    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | name                    | calls method                                                                                                  |
    +=========================+===============================================================================================================+
    | efp                     | effective fragment potential (EFP) :ref:`[manual] <sec:libefp>`                                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT) :ref:`[manual] <sec:scf>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | hf                      | HF self consistent field (SCF) :ref:`[manual] <sec:scf>`                                                      |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dcft                    | density cumulant functional theory :ref:`[manual] <sec:dcft>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mcscf                   | multiconfigurational self consistent field (SCF)                                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp2                     | 2nd-order Moller-Plesset perturbation theory (MP2) :ref:`[manual] <sec:dfmp2>` :ref:`[details] <tlmp2>`       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp3                     | 3rd-order Moller-Plesset perturbation theory (MP3) :ref:`[details] <tlmp3>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-mp3                 | MP3 with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp2.5                   | average of MP2 and MP3 :ref:`[details] <tlmp25>`                                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp4(sdq)                | 4th-order MP perturbation theory (MP4) less triples :ref:`[manual] <sec:fnompn>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-mp4(sdq)            | MP4 (less triples) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp4                     | full MP4 :ref:`[manual] <sec:fnompn>` :ref:`[details] <tlmp4>`                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-mp4                 | full MP4 with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp\ *n*                 | *n*\ th-order Moller--Plesset (MP) perturbation theory :ref:`[manual] <sec:arbpt>`                            |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | zapt\ *n*               | *n*\ th-order z-averaged perturbation theory (ZAPT) :ref:`[manual] <sec:arbpt>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp2                    | orbital-optimized second-order MP perturbation theory :ref:`[manual] <sec:occ_oo>`                            |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs-omp2                | spin-component scaled OMP2 :ref:`[manual] <sec:occ_oo>`                                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs(n)-omp2             | a special version of SCS-OMP2 for nucleobase interactions :ref:`[manual] <sec:occ_oo>`                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs-omp2-vdw            | a special version of SCS-OMP2 (from ethene dimers) :ref:`[manual] <sec:occ_oo>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sos-omp2                | spin-opposite scaled OMP2 :ref:`[manual] <sec:occ_oo>`                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sos-pi-omp2             | A special version of SOS-OMP2 for pi systems :ref:`[manual] <sec:occ_oo>`                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp3                    | orbital-optimized third-order MP perturbation theory :ref:`[manual] <sec:occ_oo>`                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs-omp3                | spin-component scaled OMP3 :ref:`[manual] <sec:occ_oo>`                                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs(n)-omp3             | a special version of SCS-OMP3 for nucleobase interactions :ref:`[manual] <sec:occ_oo>`                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | scs-omp3-vdw            | a special version of SCS-OMP3 (from ethene dimers) :ref:`[manual] <sec:occ_oo>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sos-omp3                | spin-opposite scaled OMP3 :ref:`[manual] <sec:occ_oo>`                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sos-pi-omp3             | A special version of SOS-OMP3 for pi systems :ref:`[manual] <sec:occ_oo>`                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp2.5                  | orbital-optimized MP2.5 :ref:`[manual] <sec:occ_oo>`                                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ocepa(0)                | orbital-optimized coupled electron pair approximation :ref:`[manual] <sec:occ_oo>`                            |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cepa(0)                 | coupled electron pair approximation variant 0 :ref:`[manual] <sec:fnocepa>` :ref:`[details] <tllccsd>`        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-cepa(0)             | CEPA(0) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cepa(1)                 | coupled electron pair approximation variant 1 :ref:`[manual] <sec:fnocepa>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-cepa(1)             | CEPA(1) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cepa(3)                 | coupled electron pair approximation variant 3 :ref:`[manual] <sec:fnocepa>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-cepa(3)             | CEPA(3) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | acpf                    | averaged coupled-pair functional :ref:`[manual] <sec:fnocepa>`                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-acpf                | ACPF with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | aqcc                    | averaged quadratic coupled cluster :ref:`[manual] <sec:fnocepa>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-aqcc                | AQCC with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | qcisd                   | quadratic CI singles doubles (QCISD) :ref:`[manual] <sec:fnocc>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-qcisd               | QCISD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cc2                     | approximate coupled cluster singles and doubles (CC2) :ref:`[manual] <sec:cc>`                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccd                     | coupled cluster doubles  (CCD) :ref:`[manual] <sec:cc>`                                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD) :ref:`[manual] <sec:cc>` :ref:`[details] <tllccsd>`                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | bccd                    | Brueckner coupled cluster doubles (BCCD) :ref:`[manual] <sec:cc>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | lccd                    | Linear CCD :ref:`[manual] <sec:cc>` :ref:`[details] <tllccd>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | olccd                   | orbital optimized LCCD :ref:`[manual] <sec:occ_oo>`                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-lccd                | LCCD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-ccsd                | CCSD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | qcisd(t)                | QCISD with perturbative triples :ref:`[manual] <sec:fnocc>`                                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-qcisd(t)            | QCISD(T) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd(t)                 | CCSD with perturbative triples (CCSD(T)) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsdt>`                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd(at)                | CCSD with asymmetric perturbative triples (CCSD(AT)) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsdat>`     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | bccd(t)                 | BCCD with perturbative triples :ref:`[manual] <sec:cc>`                                                       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-ccsd(t)             | CCSD(T) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cc3                     | approximate CC singles, doubles, and triples (CC3) :ref:`[manual] <sec:cc>`                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccenergy                | **expert** full control over ccenergy module                                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dfocc                   | orbital optimized CC with density fitting :ref:`[manual] <sec:occ_oo>`                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cisd                    | configuration interaction (CI) singles and doubles (CISD) :ref:`[manual] <sec:ci>` :ref:`[details] <tlcisd>`  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fno-cisd                | CISD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cisdt                   | CI singles, doubles, and triples (CISDT) :ref:`[manual] <sec:ci>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cisdtq                  | CI singles, doubles, triples, and quadruples (CISDTQ) :ref:`[manual] <sec:ci>`                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ci\ *n*                 | *n*\ th-order CI :ref:`[manual] <sec:ci>`                                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | fci                     | full configuration interaction (FCI) :ref:`[manual] <sec:ci>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | detci                   | **expert** full control over detci module                                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | casscf                  | complete active space self consistent field (CASSCF)  :ref:`[manual] <sec:ci>`                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | rasscf                  | restricted active space self consistent field (RASSCF)  :ref:`[manual] <sec:ci>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt0                   | 0th-order symmetry adapted perturbation theory (SAPT) :ref:`[manual] <sec:sapt>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2                   | 2nd-order SAPT, traditional definition :ref:`[manual] <sec:sapt>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+                  | SAPT including all 2nd-order terms :ref:`[manual] <sec:sapt>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)               | SAPT including perturbative triples :ref:`[manual] <sec:sapt>`                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3                 | SAPT including all 3rd-order terms :ref:`[manual] <sec:sapt>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(ccd)             | SAPT2+ with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                    |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)          | SAPT2+(3) with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3(ccd)            | SAPT2+3 with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+dmp2              | SAPT including all 2nd-order terms and MP2 correction :ref:`[manual] <sec:sapt>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)dmp2           | SAPT including perturbative triples and MP2 correction :ref:`[manual] <sec:sapt>`                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3dmp2             | SAPT including all 3rd-order terms and MP2 correction :ref:`[manual] <sec:sapt>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(ccd)dmp2         | SAPT2+ with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)dmp2      | SAPT2+(3) with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3(ccd)dmp2        | SAPT2+3 with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt0-ct                | 0th-order SAPT plus charge transfer (CT) calculation :ref:`[manual] <sec:saptct>`                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2-ct                | SAPT2 plus CT :ref:`[manual] <sec:saptct>`                                                                    |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+-ct               | SAPT2+ plus CT :ref:`[manual] <sec:saptct>`                                                                   |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)-ct            | SAPT2+(3) plus CT :ref:`[manual] <sec:saptct>`                                                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3-ct              | SAPT2+3 plus CT :ref:`[manual] <sec:saptct>`                                                                  |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(ccd)-ct          | SAPT2+(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+(3)(ccd)-ct       | SAPT2+(3)(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | sapt2+3(ccd)-ct         | SAPT2+3(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | adc                     | 2nd-order algebraic diagrammatic construction (ADC) :ref:`[manual] <sec:adc>`                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | eom-cc2                 | EOM-CC2 :ref:`[manual] <sec:eomcc>`                                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | eom-ccsd                | equation of motion (EOM) CCSD :ref:`[manual] <sec:eomcc>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | eom-cc3                 | EOM-CC3 :ref:`[manual] <sec:eomcc>`                                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dmrgscf                 | density matrix renormalization group SCF :ref:`[manual] <sec:dmrg>`                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dmrgci                  | density matrix renormalization group CI :ref:`[manual] <sec:dmrg>`                                            |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+

    .. include:: autodoc_dft_energy.rst

    .. include:: mrcc_table_energy.rst

    .. include:: cfour_table_energy.rst

    :examples:

    >>> # [1] Coupled-cluster singles and doubles calculation with psi code
    >>> energy('ccsd')

    >>> # [2] Charge-transfer SAPT calculation with scf projection from small into
    >>> #     requested basis, with specified projection fitting basis
    >>> set basis_guess true
    >>> set df_basis_guess jun-cc-pVDZ-JKFIT
    >>> energy('sapt0-ct')

    >>> # [3] Arbitrary-order MPn calculation
    >>> energy('mp7')

    >>> # [4] Converge scf as singlet, then run detci as triplet upon singlet reference
    >>> # Note that the integral transformation is not done automatically when detci is run in a separate step.
    >>> molecule H2 {\n0 1\nH\nH 1 0.74\n}
    >>> set global basis cc-pVDZ
    >>> set global reference rohf
    >>> scf_e, scf_wfn = energy('scf', return_wfn = True)
    >>> H2.set_multiplicity(3)
    >>> psi4.MintsHelper(scf_wfn.basisset()).integrals()
    >>> energy('detci', ref_wfn=scf_wfn)

    >>> # [5] Run two CI calculations, keeping the integrals generated in the first one.
    >>> molecule ne {\nNe\n}
    >>> set globals  basis cc-pVDZ
    >>> cisd_e, cisd_wfn = energy('cisd', return_wfn = True)
    >>> energy('fci', ref_wfn=cisd_wfn)

    """

    if hasattr(name, '__call__'):
        return name(energy, **kwargs)

    lowername = name.lower()

    # Do a cp thing
    if kwargs.get('bsse_type', None) is not None:
        return _nbody_gufunc(energy, lowername, **kwargs)

    # Check if this is a CBS extrapolation
    if "/" in lowername:
        return _cbs_gufunc('energy', lowername, **kwargs)

    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    psi4.clean_variables()

    optstash = p4util.OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['E_CONVERGENCE'])

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # Allow specification of methods to arbitrary order
    lowername, level = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    for precallback in hooks['energy']['pre']:
        precallback(lowername, **kwargs)

    _set_convergence_criterion('energy', lowername, 6, 8, 6, 8, 6)

    # Before invoking the procedure, we rename any file that should be read.
    # This is a workaround to do restarts with the current PSI4 capabilities
    # before actual, clean restarts are put in there
    # Restartfile is always converted to a single-element list if
    # it contains a single string
    if 'restart_file' in kwargs:
        restartfile = kwargs['restart_file']  # Option still available for procedure-specific action
        if restartfile != list(restartfile):
            restartfile = [restartfile]
        # Rename the files to be read to be consistent with psi4's file system
        for item in restartfile:
            name_split = re.split(r'\.', item)
            filenum = name_split[len(name_split) - 1]
            try:
                filenum = int(filenum)
            except ValueError:
                filenum = 32  # Default file number is the checkpoint one
            psioh = psi4.IOManager.shared_object()
            psio = psi4.IO.shared_object()
            filepath = psioh.get_file_path(filenum)
            namespace = psio.get_default_namespace()
            pid = str(os.getpid())
            prefix = 'psi'
            targetfile = filepath + prefix + '.' + pid + '.' + namespace + '.' + str(filenum)
            shutil.copy(item, targetfile)

    wfn = procedures['energy'][lowername](lowername, molecule=molecule, **kwargs)

    for postcallback in hooks['energy']['post']:
        postcallback(lowername, **kwargs)

    optstash.restore()
    if return_wfn:  # TODO current energy safer than wfn.energy() for now, but should be revisited

        # TODO place this with the associated call, very awkward to call this in other areas at the moment
        if lowername in ['EFP', 'MRCC', 'DMRG', 'PSIMRCC']:
            psi4.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
            psi4.print_out("The returned wavefunction is the incoming reference wavefunction.\n\n")
        elif 'sapt' in lowername:
            psi4.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
            psi4.print_out("The returned wavefunction is the dimer SCF wavefunction.\n\n")

        return (psi4.get_variable('CURRENT ENERGY'), wfn)
    else:
        return psi4.get_variable('CURRENT ENERGY')


def gradient(name, **kwargs):
    r"""Function complementary to :py:func:~driver.optimize(). Carries out one gradient pass,
    deciding analytic or finite difference.

    :returns: :ref:`Matrix<sec:psimod_Matrix>` |w--w| Total electronic gradient in Hartrees/Bohr.

    :returns: (:ref:`Matrix<sec:psimod_Matrix>`, :ref:`Wavefunction<sec:psimod_Wavefunction>`) |w--w| gradient and wavefunction when **return_wfn** specified.

    :examples:

    >>> # [1] Single-point dft gradient getting the gradient
    >>> #     in file, psi4.Matrix, and np.array forms
    >>> set gradient_write on
    >>> G, wfn = gradient('b3lyp-d', return_wfn=True)
    >>> wfn.gradient().print_out()
    >>> np.array(G)

    """
    lowername = name.lower()

    # Check if this is a CBS extrapolation
    if "/" in lowername:
        return _cbs_gufunc('gradient', lowername, **kwargs)

    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    psi4.clean_variables()
    dertype = 1

    # Prevent methods that do not have associated energies
    if lowername in energy_only_methods:
    	raise ValidationError("gradient('%s') does not have an associated gradient" % name)

    optstash = p4util.OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['E_CONVERGENCE'])

    # Allow specification of methods to arbitrary order
    lowername, level = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    dertype, func = _find_derivative_type('gradient', lowername, kwargs.pop('dertype', None))

    # no analytic derivatives for scf_type cd
    if psi4.get_option('SCF', 'SCF_TYPE') == 'CD':
        if (dertype == 1):
            raise ValidationError("""No analytic derivatives for SCF_TYPE CD.""")

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # S/R: Mode of operation- whether finite difference opt run in one job or files farmed out
    opt_mode = kwargs.get('mode', 'continuous').lower()
    if opt_mode == 'continuous':
        pass
    elif opt_mode == 'sow':
        if dertype == 1:
            raise ValidationError("""Optimize execution mode 'sow' not valid for analytic gradient calculation.""")
    elif opt_mode == 'reap':
        opt_linkage = kwargs.get('linkage', None)
        if opt_linkage is None:
            raise ValidationError("""Optimize execution mode 'reap' requires a linkage option.""")
    else:
        raise ValidationError("""Optimize execution mode '%s' not valid.""" % (opt_mode))

    # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
    _set_convergence_criterion('energy', lowername, 8, 10, 8, 10, 8)

    # Does dertype indicate an analytic procedure both exists and is wanted?
    if dertype == 1:
        psi4.print_out("""gradient() will perform analytic gradient computation.\n""")

        # Perform the gradient calculation
        wfn = procedures['gradient'][lowername](lowername, molecule=molecule, **kwargs)

        optstash.restore()
        if return_wfn:
            return (wfn.gradient(), wfn)
        else:
            return wfn.gradient()

    else:
        psi4.print_out("""gradient() will perform gradient computation by finite difference of analytic energies.\n""")

        opt_iter = kwargs.get('opt_iter', 1)

        if opt_iter == 1:
            print('Performing finite difference calculations')

        # Shifting the geometry so need to copy the active molecule
        moleculeclone = molecule.clone()

        # Obtain list of displacements
        displacements = psi4.fd_geoms_1_0(moleculeclone)
        ndisp = len(displacements)

        # This version is pretty dependent on the reference geometry being last (as it is now)
        print(""" %d displacements needed ...""" % (ndisp), end='')
        energies = []

        # S/R: Write instructions for sow/reap procedure to output file and reap input file
        if opt_mode == 'sow':
            instructionsO = """\n    The optimization sow/reap procedure has been selected through mode='sow'. In addition\n"""
            instructionsO += """    to this output file (which contains no quantum chemical calculations), this job\n"""
            instructionsO += """    has produced a number of input files (OPT-%s-*.in) for individual components\n""" % (str(opt_iter))
            instructionsO += """    and a single input file (OPT-master.in) with an optimize(mode='reap') command.\n"""
            instructionsO += """    These files may look very peculiar since they contain processed and pickled python\n"""
            instructionsO += """    rather than normal input. Follow the instructions in OPT-master.in to continue.\n\n"""
            instructionsO += """    Alternatively, a single-job execution of the gradient may be accessed through\n"""
            instructionsO += """    the optimization wrapper option mode='continuous'.\n\n"""
            psi4.print_out(instructionsO)

            instructionsM = """\n#    Follow the instructions below to carry out this optimization cycle.\n#\n"""
            instructionsM += """#    (1)  Run all of the OPT-%s-*.in input files on any variety of computer architecture.\n""" % (str(opt_iter))
            instructionsM += """#       The output file names must be as given below.\n#\n"""
            for rgt in range(ndisp):
                pre = 'OPT-' + str(opt_iter) + '-' + str(rgt + 1)
                instructionsM += """#             psi4 -i %-27s -o %-27s\n""" % (pre + '.in', pre + '.out')
            instructionsM += """#\n#    (2)  Gather all the resulting output files in a directory. Place input file\n"""
            instructionsM += """#         OPT-master.in into that directory and run it. The job will be minimal in\n"""
            instructionsM += """#         length and give summary results for the gradient step in its output file.\n#\n"""
            if opt_iter == 1:
                instructionsM += """#             psi4 -i %-27s -o %-27s\n#\n""" % ('OPT-master.in', 'OPT-master.out')
            else:
                instructionsM += """#             psi4 -a -i %-27s -o %-27s\n#\n""" % ('OPT-master.in', 'OPT-master.out')
            instructionsM += """#    After each optimization iteration, the OPT-master.in file is overwritten so return here\n"""
            instructionsM += """#    for new instructions. With the use of the psi4 -a flag, OPT-master.out is not\n"""
            instructionsM += """#    overwritten and so maintains a history of the job. To use the (binary) optimizer\n"""
            instructionsM += """#    data file to accelerate convergence, the OPT-master jobs must run on the same computer.\n\n"""

            with open('OPT-master.in', 'wb') as fmaster:
                fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n'.encode('utf-8'))
                fmaster.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
                fmaster.write(p4util.format_options_for_input().encode('utf-8'))
                p4util.format_kwargs_for_input(fmaster, lmode=2, return_wfn=True, **kwargs)
                fmaster.write(("""retE, retwfn = %s('%s', **kwargs)\n\n""" % (optimize.__name__, lowername)).encode('utf-8'))
                fmaster.write(instructionsM.encode('utf-8'))

        for n, displacement in enumerate(displacements):
            rfile = 'OPT-%s-%s' % (opt_iter, n + 1)

            # Build string of title banner
            banners = ''
            banners += """psi4.print_out('\\n')\n"""
            banners += """p4util.banner(' Gradient %d Computation: Displacement %d ')\n""" % (opt_iter, n + 1)
            banners += """psi4.print_out('\\n')\n\n"""

            if opt_mode == 'continuous':

                # print progress to file and screen
                psi4.print_out('\n')
                p4util.banner('Loading displacement %d of %d' % (n + 1, ndisp))
                print(""" %d""" % (n + 1), end=('\n' if (n + 1 == ndisp) else ''))
                sys.stdout.flush()

                # Load in displacement into the active molecule
                moleculeclone.set_geometry(displacement)

                # Perform the energy calculation
                E, wfn = func(lowername, return_wfn=True, molecule=moleculeclone, **kwargs)
                energies.append(psi4.get_variable('CURRENT ENERGY'))

            # S/R: Write each displaced geometry to an input file
            elif opt_mode == 'sow':
                moleculeclone.set_geometry(displacement)

                # S/R: Prepare molecule, options, and kwargs
                with open('%s.in' % (rfile), 'wb') as freagent:
                    freagent.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n'.encode('utf-8'))
                    freagent.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
                    freagent.write(p4util.format_options_for_input().encode('utf-8'))
                    p4util.format_kwargs_for_input(freagent, **kwargs)

                    # S/R: Prepare function call and energy save
                    freagent.write(("""electronic_energy = %s('%s', **kwargs)\n\n""" % (func.__name__, lowername)).encode('utf-8'))
                    freagent.write(("""psi4.print_out('\\nGRADIENT RESULT: computation %d for item %d """ % (os.getpid(), n + 1)).encode('utf-8'))
                    freagent.write("""yields electronic energy %20.12f\\n' % (electronic_energy))\n\n""".encode('utf-8'))

            # S/R: Read energy from each displaced geometry output file and save in energies array
            elif opt_mode == 'reap':
                exec(banners)
                psi4.set_variable('NUCLEAR REPULSION ENERGY', moleculeclone.nuclear_repulsion_energy())
                energies.append(p4util.extract_sowreap_from_output(rfile, 'GRADIENT', n, opt_linkage, True))

        # S/R: Quit sow after writing files. Initialize skeleton wfn to receive grad for reap
        if opt_mode == 'sow':
            optstash.restore()
            if return_wfn:
                return (None, None)  # any point to building a dummy wfn here?
            else:
                return None
        elif opt_mode == 'reap':
            psi4.set_variable('CURRENT ENERGY', energies[-1])
            wfn = psi4.new_wavefunction(molecule, psi4.get_global_option('BASIS'))

        # Compute the gradient; last item in 'energies' is undisplaced
        psi4.set_local_option('FINDIF', 'GRADIENT_WRITE', True)
        G = psi4.fd_1_0(molecule, energies)
        G.print_out()
        wfn.set_gradient(G)

        optstash.restore()

        if return_wfn:
            return (wfn.gradient(), wfn)
        else:
            return wfn.gradient()


def property(name, **kwargs):
    r"""Function to compute various properties.

    :aliases: prop()

    :returns: none.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - This function at present has a limited functionality.
         Consult the keywords sections of other modules for further property capabilities.

    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | Name               | Calls Method                                  | Reference      | Supported Properties                                          |
    +====================+===============================================+================+===============================================================+
    | scf                | Self-consistent field method(s)               | RHF/ROHF/UHF   | Listed :ref:`here <sec:oeprop>`                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | hf                 | HF Self-consistent field method(s)            | RHF/ROHF/UHF   | Listed :ref:`here <sec:oeprop>`                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | mp2                | MP2 ('mp2_type df' only)                      | RHF            | Listed :ref:`here <sec:oeprop>`                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | cc2                | 2nd-order approximate CCSD                    | RHF            | dipole, quadrupole, polarizability, rotation, roa             |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | ccsd               | Coupled cluster singles and doubles (CCSD)    | RHF            | dipole, quadrupole, polarizability, rotation, roa             |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | eom-cc2            | 2nd-order approximate EOM-CCSD                | RHF            | oscillator_strength, rotational_strength                      |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | eom-ccsd           | Equation-of-motion CCSD (EOM-CCSD)            | RHF            | oscillator_strength, rotational_strength                      |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | cisd, cisdt,       | Configuration interaction                     | RHF/ROHF       | Listed :ref:`here <sec:oeprop>`, transition_dipole,           |
    | cisdt, cisdtq,     |                                               |                | transition_quadrupole                                         |
    | ci5, ..., fci      |                                               |                |                                                               |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
    | casscf, rasscf     | Multi-configurational SCF                     | RHF/ROHF       | Listed :ref:`here <sec:oeprop>`, transition_dipole,           |
    |                    |                                               |                | transition_quadrupole                                         |
    +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+

    :type name: string
    :param name: ``'ccsd'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type properties: array of strings
    :param properties: |dl| ``[]`` |dr| || ``['rotation', 'polarizability', 'oscillator_strength', 'roa']`` || etc.

        Indicates which properties should be computed. Defaults to dipole and quadrupole.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :examples:

    >>> # [1] Optical rotation calculation
    >>> property('cc2', properties=['rotation'])

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)

    optstash = p4util.OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['E_CONVERGENCE'])

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # Allow specification of methods to arbitrary order
    lowername, level = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    properties = kwargs.get('properties', ['dipole', 'quadrupole'])
    kwargs['properties'] = p4util.drop_duplicates(properties)
    _set_convergence_criterion('property', lowername, 6, 10, 6, 10, 8)
    wfn = procedures['property'][lowername](lowername, **kwargs)

    optstash.restore()

    if return_wfn:
        return (psi4.get_variable('CURRENT ENERGY'), wfn)
    else:
        return psi4.get_variable('CURRENT ENERGY')


def optimize(name, **kwargs):
    r"""Function to perform a geometry optimization.

    :aliases: opt()

    :returns: *float* |w--w| Total electronic energy of optimized structure in Hartrees.

    :returns: (*float*, :ref:`Wavefunction<sec:psimod_Wavefunction>`) |w--w| energy and wavefunction when **return_wfn** specified.

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`CURRENT ENERGY <CURRENTENERGY>`

    .. note:: Analytic gradients area available for all methods in the table
        below. Optimizations with other methods in the energy table proceed
        by finite differences.

    .. _`table:grad_gen`:

    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | name                    | calls method                                                                                                  |
    +=========================+===============================================================================================================+
    | scf                     | Hartree--Fock (HF) or density functional theory (DFT) :ref:`[manual] <sec:scf>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | hf                      | Hartree--Fock (HF)  :ref:`[manual] <sec:scf>`                                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | dcft                    | density cumulant functional theory :ref:`[manual] <sec:dcft>`                                                 |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp2                     | 2nd-order Moller-Plesset perturbation theory (MP2) :ref:`[manual] <sec:dfmp2>` :ref:`[details] <tlmp2>`       |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp2.5                   | MP2.5 :ref:`[manual] <sec:convocc>`                                                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | mp3                     | third-order MP perturbation theory :ref:`[manual] <sec:convocc>`                                              |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp2                    | orbital-optimized second-order MP perturbation theory :ref:`[manual] <sec:occ>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp2.5                  | orbital-optimized MP2.5 :ref:`[manual] <sec:occ>`                                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | omp3                    | orbital-optimized third-order MP perturbation theory :ref:`[manual] <sec:occ>`                                |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ocepa                   | orbital-optimized coupled electron pair approximation :ref:`[manual] <sec:occ>`                               |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | cepa0                   | coupled electron pair approximation(0) :ref:`[manual] <sec:convocc>`                                          |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd                    | coupled cluster singles and doubles (CCSD) :ref:`[manual] <sec:cc>`                                           |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | ccsd(t)                 | CCSD with perturbative triples (CCSD(T)) :ref:`[manual] <sec:cc>`                                             |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | eom-ccsd                | equation of motion (EOM) CCSD :ref:`[manual] <sec:eomcc>`                                                     |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+
    | efp                     | efp-only optimizations                                                                                        |
    +-------------------------+---------------------------------------------------------------------------------------------------------------+

    .. _`table:grad_scf`:


    .. include:: autodoc_dft_opt.rst

    .. include:: cfour_table_grad.rst

    .. warning:: Optimizations where the molecule is specified in Z-matrix format
       with dummy atoms will result in the geometry being converted to a Cartesian representation.

    :type name: string
    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the database. May be any valid argument to
        :py:func:`~driver.energy`.

    :type func: :ref:`function <op_py_function>`
    :param func: |dl| ``gradient`` |dr| || ``energy`` || ``cbs``

        Indicates the type of calculation to be performed on the molecule.
        The default dertype accesses ``'gradient'`` or ``'energy'``, while
        ``'cbs'`` performs a multistage finite difference calculation.
        If a nested series of python functions is intended (see :ref:`sec:intercalls`),
        use keyword ``opt_func`` instead of ``func``.

    :type mode: string
    :param mode: |dl| ``'continuous'`` |dr| || ``'sow'`` || ``'reap'``

        For a finite difference of energies optimization, indicates whether
        the calculations required to complete the
        optimization are to be run in one file (``'continuous'``) or are to be
        farmed out in an embarrassingly parallel fashion
        (``'sow'``/``'reap'``). For the latter, run an initial job with
        ``'sow'`` and follow instructions in its output file. For maximum
        flexibility, ``return_wfn`` is always on in ``'reap'`` mode.

    :type dertype: :ref:`dertype <op_py_dertype>`
    :param dertype: ``'gradient'`` || ``'energy'``

        Indicates whether analytic (if available) or finite difference
        optimization is to be performed.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type hessian_with: string
    :param hessian_with: ``'scf'`` || ``'mp2'`` || etc.

        Indicates the computational method with which to perform a hessian
        analysis to guide the geometry optimization.

    :examples:

    >>> # [1] Analytic hf optimization
    >>> optimize('hf')

    >>> # [2] Finite difference mp5 optimization with gradient
    >>> #     printed to output file
    >>> e, wfn = opt('mp5', return_wfn='yes')
    >>> wfn.gradient().print_out()

    >>> # [3] Forced finite difference hf optimization run in
    >>> #     embarrassingly parallel fashion
    >>> optimize('hf', dertype='energy', mode='sow')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)

    full_hess_every = psi4.get_option('OPTKING', 'FULL_HESS_EVERY')
    steps_since_last_hessian = 0
    hessian_with_method = kwargs.get('hessian_with', lowername)

    # are we in sow/reap mode?
    opt_mode = kwargs.get('mode', 'continuous').lower()
    if opt_mode not in ['continuous', 'sow', 'reap']:
        raise ValidationError("""Optimize execution mode '%s' not valid.""" % (opt_mode))

    optstash = p4util.OptionsState(
        ['OPTKING', 'INTRAFRAG_STEP_LIMIT'],
        ['FINDIF', 'HESSIAN_WRITE'],
        ['OPTKING', 'CART_HESS_READ'],
        ['SCF', 'GUESS_PERSIST'],  # handle on behalf of cbs()
        ['SCF', 'GUESS'])

    n = kwargs.get('opt_iter', 1)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())

    # If we are feezing cartesian, do not orient or COM
    if psi4.get_local_option("OPTKING", "FROZEN_CARTESIAN"):
        molecule.fix_orientation(True)
        molecule.fix_com(True)   
    molecule.update_geometry()

    # Shifting the geometry so need to copy the active molecule
    moleculeclone = molecule.clone()

    initial_sym = moleculeclone.schoenflies_symbol()
    while n <= psi4.get_option('OPTKING', 'GEOM_MAXITER'):
        current_sym = moleculeclone.schoenflies_symbol()
        if initial_sym != current_sym:
            raise ValidationError("""Point group changed! You should restart """
                                  """using the last geometry in the output, after """
                                  """carefully making sure all symmetry-dependent """
                                  """input, such as DOCC, is correct.""")
        kwargs['opt_iter'] = n

        # Use orbitals from previous iteration as a guess
        #   set within loop so that can be influenced by fns to optimize (e.g., cbs)
        if (n > 1) and (opt_mode == 'continuous') and (not psi4.get_option('SCF', 'GUESS_PERSIST')):
            psi4.set_local_option('SCF', 'GUESS', 'READ')

        # Before computing gradient, save previous molecule and wavefunction if this is an IRC optimization
        if (n > 1) and (psi4.get_option('OPTKING', 'OPT_TYPE') == 'IRC'):
            old_thisenergy = psi4.get_variable('CURRENT ENERGY')

        # Compute the gradient
        G, wfn = gradient(lowername, return_wfn=True, molecule=moleculeclone, **kwargs)
        thisenergy = psi4.get_variable('CURRENT ENERGY')

        # above, used to be getting energy as last of energy list from gradient()
        # thisenergy below should ultimately be testing on wfn.energy()

        # S/R: Quit after getting new displacements or if forming gradient fails
        if opt_mode  == 'sow':
            return (0.0, None)
        elif opt_mode == 'reap' and thisenergy == 0.0:
            return (0.0, None)

        psi4.set_gradient(G)

        # S/R: Move opt data file from last pass into namespace for this pass
        if opt_mode == 'reap' and n != 0:
            psi4.IOManager.shared_object().set_specific_retention(1, True)
            psi4.IOManager.shared_object().set_specific_path(1, './')
            if 'opt_datafile' in kwargs:
                restartfile = kwargs.pop('opt_datafile')
                #if psi4.me() == 0:  TODO ask Ryan
                shutil.copy(restartfile, p4util.get_psifile(1))

        opt_func = kwargs.get('opt_func', kwargs.get('func', energy))
        if opt_func.__name__ == 'complete_basis_set':
            psi4.IOManager.shared_object().set_specific_retention(1, True)

        if full_hess_every > -1:
            psi4.set_global_option('HESSIAN_WRITE', True)

        # compute Hessian as requested; frequency wipes out gradient so stash it
        if ((full_hess_every > -1) and (n == 1)) or (steps_since_last_hessian + 1 == full_hess_every):
            G = psi4.get_gradient()  # TODO
            psi4.IOManager.shared_object().set_specific_retention(1, True)
            psi4.IOManager.shared_object().set_specific_path(1, './')
            frequencies(hessian_with_method, **kwargs)
            steps_since_last_hessian = 0
            psi4.set_gradient(G)
            psi4.set_global_option('CART_HESS_READ', True)
        elif (full_hess_every == -1) and psi4.get_global_option('CART_HESS_READ') and (n == 1):
            pass
            # Do nothing; user said to read existing hessian once
        else:
            psi4.set_global_option('CART_HESS_READ', False)
            steps_since_last_hessian += 1

        # Take step. communicate to/from/within optking through legacy_molecule
        psi4.set_legacy_molecule(moleculeclone)
        optking_rval = psi4.optking()
        moleculeclone = psi4.get_legacy_molecule()
        moleculeclone.update_geometry()
        if optking_rval == psi4.PsiReturnType.EndLoop:
            # if this is the end of an IRC run, set wfn, energy, and molecule to that
            # of the last optimized IRC point
            if psi4.get_option('OPTKING', 'OPT_TYPE') == 'IRC':
                thisenergy = old_thisenergy
            print('Optimizer: Optimization complete!')
            psi4.print_out('\n    Final optimized geometry and variables:\n')
            moleculeclone.print_in_input_format()
            # Check if user wants to see the intcos; if so, don't delete them.
            if psi4.get_option('OPTKING', 'INTCOS_GENERATE_EXIT') == False:
                if psi4.get_option('OPTKING', 'KEEP_INTCOS') == False:
                    psi4.opt_clean()
            # Changing environment to optimized geometry as expected by user
            molecule.set_geometry(moleculeclone.geometry())
            for postcallback in hooks['optimize']['post']:
                postcallback(lowername, **kwargs)
            psi4.clean()

            # S/R: Clean up opt input file
            if opt_mode == 'reap':
                with open('OPT-master.in', 'wb') as fmaster:
                    fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n'.encode('utf-8'))
                    fmaster.write('# Optimization complete!\n\n'.encode('utf-8'))

            if opt_func.__name__ == 'complete_basis_set':
                psi4.IOManager.shared_object().set_specific_retention(1, False)

            optstash.restore()

            if return_wfn:
                return (thisenergy, wfn)
            else:
                return thisenergy

        elif optking_rval == psi4.PsiReturnType.Failure:
            print('Optimizer: Optimization failed!')
            if (psi4.get_option('OPTKING', 'KEEP_INTCOS') == False):
                psi4.opt_clean()
            molecule.set_geometry(moleculeclone.geometry())
            psi4.clean()
            optstash.restore()
            return thisenergy

        psi4.print_out('\n    Structure for next step:\n')
        moleculeclone.print_in_input_format()

        # S/R: Preserve opt data file for next pass and switch modes to get new displacements
        if opt_mode == 'reap':
            kwargs['opt_datafile'] = p4util.get_psifile(1)
            kwargs['mode'] = 'sow'

        n += 1

    psi4.print_out('\tOptimizer: Did not converge!')
    if psi4.get_option('OPTKING', 'INTCOS_GENERATE_EXIT') == False:
        if psi4.get_option('OPTKING', 'KEEP_INTCOS') == False:
            psi4.opt_clean()

    optstash.restore()


def parse_arbitrary_order(name):
    r"""Function to parse name string into a method family like CI or MRCC and specific
    level information like 4 for CISDTQ or MRCCSDTQ.

    """
    # matches 'mrccsdt(q)'
    if name.startswith('mrcc'):

        # avoid undoing fn's good work when called twice
        if name == 'mrcc':
            return name, None

        # grabs 'sdt(q)'
        ccfullname = name[4:]

        # A negative order indicates perturbative method
        methods = {
            'sd'          : { 'method': 1, 'order':  2, 'fullname': 'CCSD'         },
            'sdt'         : { 'method': 1, 'order':  3, 'fullname': 'CCSDT'        },
            'sdtq'        : { 'method': 1, 'order':  4, 'fullname': 'CCSDTQ'       },
            'sdtqp'       : { 'method': 1, 'order':  5, 'fullname': 'CCSDTQP'      },
            'sdtqph'      : { 'method': 1, 'order':  6, 'fullname': 'CCSDTQPH'     },
            'sd(t)'       : { 'method': 3, 'order': -3, 'fullname': 'CCSD(T)'      },
            'sdt(q)'      : { 'method': 3, 'order': -4, 'fullname': 'CCSDT(Q)'     },
            'sdtq(p)'     : { 'method': 3, 'order': -5, 'fullname': 'CCSDTQ(P)'    },
            'sdtqp(h)'    : { 'method': 3, 'order': -6, 'fullname': 'CCSDTQP(H)'   },
            'sd(t)_l'     : { 'method': 4, 'order': -3, 'fullname': 'CCSD(T)_L'    },
            'sdt(q)_l'    : { 'method': 4, 'order': -4, 'fullname': 'CCSDT(Q)_L'   },
            'sdtq(p)_l'   : { 'method': 4, 'order': -5, 'fullname': 'CCSDTQ(P)_L'  },
            'sdtqp(h)_l'  : { 'method': 4, 'order': -6, 'fullname': 'CCSDTQP(H)_L' },
            'sdt-1a'      : { 'method': 5, 'order':  3, 'fullname': 'CCSDT-1a'     },
            'sdtq-1a'     : { 'method': 5, 'order':  4, 'fullname': 'CCSDTQ-1a'    },
            'sdtqp-1a'    : { 'method': 5, 'order':  5, 'fullname': 'CCSDTQP-1a'   },
            'sdtqph-1a'   : { 'method': 5, 'order':  6, 'fullname': 'CCSDTQPH-1a'  },
            'sdt-1b'      : { 'method': 6, 'order':  3, 'fullname': 'CCSDT-1b'     },
            'sdtq-1b'     : { 'method': 6, 'order':  4, 'fullname': 'CCSDTQ-1b'    },
            'sdtqp-1b'    : { 'method': 6, 'order':  5, 'fullname': 'CCSDTQP-1b'   },
            'sdtqph-1b'   : { 'method': 6, 'order':  6, 'fullname': 'CCSDTQPH-1b'  },
            '2'           : { 'method': 7, 'order':  2, 'fullname': 'CC2'          },
            '3'           : { 'method': 7, 'order':  3, 'fullname': 'CC3'          },
            '4'           : { 'method': 7, 'order':  4, 'fullname': 'CC4'          },
            '5'           : { 'method': 7, 'order':  5, 'fullname': 'CC5'          },
            '6'           : { 'method': 7, 'order':  6, 'fullname': 'CC6'          },
            'sdt-3'       : { 'method': 8, 'order':  3, 'fullname': 'CCSDT-3'      },
            'sdtq-3'      : { 'method': 8, 'order':  4, 'fullname': 'CCSDTQ-3'     },
            'sdtqp-3'     : { 'method': 8, 'order':  5, 'fullname': 'CCSDTQP-3'    },
            'sdtqph-3'    : { 'method': 8, 'order':  6, 'fullname': 'CCSDTQPH-3'   }
        }

        # looks for 'sdt(q)' in dictionary
        if ccfullname in methods:
            return 'mrcc', methods[ccfullname]
        else:
            raise ValidationError('MRCC method \'%s\' invalid.' % (name))

    elif re.match(r'^[a-z]+\d+$', name):
        decompose = re.compile(r'^([a-z]+)(\d+)$').match(name)
        namestump = decompose.group(1)
        namelevel = int(decompose.group(2))

        if namestump in ['mp', 'zapt', 'ci']:
            # Let mp2, mp3, mp4 pass through to select functions
            if namestump == 'mp' and namelevel in [2, 3, 4]:
                return name, None
            # Otherwise return method and order
            else:
                return namestump, namelevel
        else:
            return name, None
    else:
        return name, None


def hessian(name, **kwargs):
    r"""Function complementary to :py:func:`~frequency`. Computes force
    constants, deciding analytic, finite difference of gradients, or
    finite difference of energies.

    :returns: :ref:`Matrix<sec:psimod_Matrix>` |w--w| Total non-mass-weighted electronic Hessian in Hartrees/Bohr/Bohr.

    :returns: (:ref:`Matrix<sec:psimod_Matrix>`, :ref:`Wavefunction<sec:psimod_Wavefunction>`) |w--w| Hessian and wavefunction when **return_wfn** specified.

    :examples:

    >>> # [1] Frequency calculation without thermochemical analysis
    >>> hessian('mp3')

    >>> # [2] Frequency calc w/o thermo analysis getting the Hessian
    >>> #     in file, psi4.Matrix, and np.array forms
    >>> set hessian_write on
    >>> H, wfn = hessian('ccsd', return_wfn=True)
    >>> wfn.hessian().print_out()
    >>> np.array(H)

    """
    lowername = name.lower()

    # Check if this is a CBS extrapolation
    if "/" in lowername:
        return _cbs_gufunc('hessian', lowername, **kwargs)

    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    psi4.clean_variables()
    dertype = 2

    # Prevent methods that do not have associated energies
    if lowername in energy_only_methods:
	    raise ValidationError("hessian('%s') does not have an associated hessian" % name)

    optstash = p4util.OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['FINDIF', 'HESSIAN_WRITE'],
        ['E_CONVERGENCE'])

    # Allow specification of methods to arbitrary order
    lowername, level = parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    dertype, func = _find_derivative_type('hessian', lowername, kwargs.pop('dertype', None))

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # S/R: Mode of operation- whether finite difference freq run in one job or files farmed out
    freq_mode = kwargs.get('mode', 'continuous').lower()
    if freq_mode == 'continuous':
        pass
    elif freq_mode == 'sow':
        if dertype == 2:
            raise ValidationError("""Frequency execution mode 'sow' not valid for analytic Hessian calculation.""")
    elif freq_mode == 'reap':
        freq_linkage = kwargs.get('linkage', None)
        if freq_linkage is None:
            raise ValidationError("""Frequency execution mode 'reap' requires a linkage option.""")
    else:
        raise ValidationError("""Frequency execution mode '%s' not valid.""" % (freq_mode))

    # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
    _set_convergence_criterion('energy', lowername, 8, 10, 8, 10, 8)

    # Select certain irreps
    irrep = kwargs.get('irrep', -1)
    if irrep == -1:
        pass  # do all irreps
    else:
        irrep = parse_cotton_irreps(irrep, molecule.schoenflies_symbol())
        irrep -= 1  # A1 irrep is externally 1, internally 0

    # Does an analytic procedure exist for the requested method?
    if dertype == 2:
        psi4.print_out("""hessian() will perform analytic frequency computation.\n""")

        # We have the desired method. Do it.
        wfn = procedures['hessian'][lowername](lowername, molecule=molecule, **kwargs)
        optstash.restore()

        # TODO: check that current energy's being set to the right figure when this code is actually used
        psi4.set_variable('CURRENT ENERGY', wfn.energy())

        if return_wfn:
            return (wfn.hessian(), wfn)
        else:
            return wfn.hessian()

    elif dertype == 1:
        psi4.print_out("""hessian() will perform frequency computation by finite difference of analytic gradients.\n""")

        func = procedures['gradient'][lowername]

        # Shifting the geometry so need to copy the active molecule
        moleculeclone = molecule.clone()

        # Obtain list of displacements
        displacements = psi4.fd_geoms_freq_1(moleculeclone, irrep)
        moleculeclone.reinterpret_coordentry(False)
        moleculeclone.fix_orientation(True)

        # Record undisplaced symmetry for projection of diplaced point groups
        psi4.set_parent_symmetry(molecule.schoenflies_symbol())

        ndisp = len(displacements)
        print(""" %d displacements needed.""" % ndisp)
        gradients = []
        energies = []

        # S/R: Write instructions for sow/reap procedure to output file and reap input file
        if freq_mode == 'sow':
            instructionsO = """\n#    The frequency sow/reap procedure has been selected through mode='sow'. In addition\n"""
            instructionsO += """#    to this output file (which contains no quantum chemical calculations), this job\n"""
            instructionsO += """#    has produced a number of input files (FREQ-*.in) for individual components\n"""
            instructionsO += """#    and a single input file (FREQ-master.in) with a frequency(mode='reap') command.\n"""
            instructionsO += """#    These files may look very peculiar since they contain processed and pickled python\n"""
            instructionsO += """#    rather than normal input. Follow the instructions below (repeated in FREQ-master.in)\n"""
            instructionsO += """#    to continue.\n#\n"""
            instructionsO += """#    Alternatively, a single-job execution of the hessian may be accessed through\n"""
            instructionsO += """#    the frequency wrapper option mode='continuous'.\n#\n"""
            psi4.print_out(instructionsO)

            instructionsM = """\n#    Follow the instructions below to carry out this frequency computation.\n#\n"""
            instructionsM += """#    (1)  Run all of the FREQ-*.in input files on any variety of computer architecture.\n"""
            instructionsM += """#       The output file names must be as given below (these are the defaults when executed\n"""
            instructionsM += """#       as `psi4 FREQ-1.in`, etc.).\n#\n"""
            for rgt in range(ndisp):
                pre = 'FREQ-' + str(rgt + 1)
                instructionsM += """#             psi4 -i %-27s -o %-27s\n""" % (pre + '.in', pre + '.out')
            instructionsM += """#\n#    (2)  Gather all the resulting output files in a directory. Place input file\n"""
            instructionsM += """#         FREQ-master.in into that directory and run it. The job will be minimal in\n"""
            instructionsM += """#         length and give summary results for the frequency computation in its output file.\n#\n"""
            instructionsM += """#             psi4 -i %-27s -o %-27s\n#\n\n""" % ('FREQ-master.in', 'FREQ-master.out')

            with open('FREQ-master.in', 'wb') as fmaster:
                fmaster.write('# This is a psi4 input file auto-generated from the hessian() wrapper.\n\n'.encode('utf-8'))
                fmaster.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
                fmaster.write(p4util.format_options_for_input(moleculeclone, **kwargs))
                p4util.format_kwargs_for_input(fmaster, lmode=2, return_wfn=True, **kwargs)
                fmaster.write(("""retE, retwfn = %s('%s', **kwargs)\n\n""" % (frequency.__name__, lowername)).encode('utf-8'))
                fmaster.write(instructionsM.encode('utf-8'))
            psi4.print_out(instructionsM)

        for n, displacement in enumerate(displacements):
            rfile = 'FREQ-%s' % (n + 1)

            # Build string of title banner
            banners = ''
            banners += """psi4.print_out('\\n')\n"""
            banners += """p4util.banner(' Hessian Computation: Gradient Displacement %d ')\n""" % (n + 1)
            banners += """psi4.print_out('\\n')\n\n"""

            if freq_mode == 'continuous':

                # print progress to file and screen
                psi4.print_out('\n')
                p4util.banner('Loading displacement %d of %d' % (n + 1, ndisp))
                print(""" %d""" % (n + 1), end=('\n' if (n + 1 == ndisp) else ''))
                sys.stdout.flush()

                # Load in displacement into the active molecule (xyz coordinates only)
                moleculeclone.set_geometry(displacement)

                # Perform the gradient calculation
                wfn = func(lowername, molecule=moleculeclone, **kwargs)
                gradients.append(wfn.gradient())
                energies.append(psi4.get_variable('CURRENT ENERGY'))

                # clean may be necessary when changing irreps of displacements
                psi4.clean()

            # S/R: Write each displaced geometry to an input file
            elif freq_mode == 'sow':
                moleculeclone.set_geometry(displacement)

                # S/R: Prepare molecule, options, kwargs, function call and energy save
                #      forcexyz in molecule writer S/R enforcement of !reinterpret_coordentry above
                with open('%s.in' % (rfile), 'wb') as freagent:
                    freagent.write('# This is a psi4 input file auto-generated from the hessian() wrapper.\n\n')
                    freagent.write(p4util.format_molecule_for_input(moleculeclone, forcexyz=True).encode('utf-8'))
                    freagent.write(p4util.format_options_for_input(moleculeclone, **kwargs).encode('utf-8'))
                    p4util.format_kwargs_for_input(freagent, **kwargs)
                    freagent.write("""wfn = %s('%s', **kwargs)\n\n""" % (func.__name__, lowername))
                    freagent.write("""psi4.print_out('\\nHESSIAN RESULT: computation %d for item %d """ % (os.getpid(), n + 1))
                    freagent.write("""yields electronic gradient %r\\n' % (p4util.mat2arr(wfn.gradient())))\n\n""")
                    freagent.write("""psi4.print_out('\\nHESSIAN RESULT: computation %d for item %d """ % (os.getpid(), n + 1))
                    freagent.write("""yields electronic energy %20.12f\\n' % (get_variable('CURRENT ENERGY')))\n\n""")

            # S/R: Read energy from each displaced geometry output file and save in energies array
            elif freq_mode == 'reap':
                exec(banners)
                psi4.set_variable('NUCLEAR REPULSION ENERGY', moleculeclone.nuclear_repulsion_energy())
                pygrad = p4util.extract_sowreap_from_output(rfile, 'HESSIAN', n, freq_linkage, True, label='electronic gradient')
                p4mat = psi4.Matrix(moleculeclone.natom(), 3)
                p4mat.set(pygrad)
                p4mat.print_out()
                gradients.append(p4mat)
                energies.append(p4util.extract_sowreap_from_output(rfile, 'HESSIAN', n, freq_linkage, True))

        # S/R: Quit sow after writing files. Initialize skeleton wfn to receive grad for reap
        if freq_mode == 'sow':
            optstash.restore()
            if return_wfn:
                return (None, None)
            else:
                return None
        elif freq_mode == 'reap':
            wfn = psi4.new_wavefunction(molecule, psi4.get_global_option('BASIS'))

        # Assemble Hessian from gradients
        #   Final disp is undisp, so wfn has mol, G, H general to freq calc
        H = psi4.fd_freq_1(molecule, gradients, irrep)  # TODO or moleculeclone?
        wfn.set_hessian(H)
        wfn.set_frequencies(psi4.get_frequencies())

        # The last item in the list is the reference energy, return it
        psi4.set_variable('CURRENT ENERGY', energies[-1])

        psi4.set_parent_symmetry('')
        optstash.restore()

        if return_wfn:
            return (wfn.hessian(), wfn)
        else:
            return wfn.hessian()

    else:
        psi4.print_out("""hessian() will perform frequency computation by finite difference of analytic energies.\n""")

        # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
        optstash.restore()
        if not psi4.has_option_changed('SCF', 'E_CONVERGENCE'):
            if procedures['energy'][lowername] in [run_scf, run_dft]:
                psi4.set_local_option('SCF', 'E_CONVERGENCE', 10)
            else:
                psi4.set_local_option('SCF', 'E_CONVERGENCE', 11)
        if not psi4.has_option_changed('SCF', 'D_CONVERGENCE'):
            if procedures['energy'][lowername] in [run_scf, run_dft]:
                psi4.set_local_option('SCF', 'D_CONVERGENCE', 10)
            else:
                psi4.set_local_option('SCF', 'D_CONVERGENCE', 11)

        # Set post-scf convergence criteria (global will cover all correlated modules)
        if not psi4.has_global_option_changed('E_CONVERGENCE'):
            if procedures['energy'][lowername] not in [run_scf, run_dft]:
                psi4.set_global_option('E_CONVERGENCE', 10)

        # Shifting the geometry so need to copy the active molecule
        moleculeclone = molecule.clone()

        # Obtain list of displacements
        displacements = psi4.fd_geoms_freq_0(moleculeclone, irrep)
        moleculeclone.fix_orientation(True)
        moleculeclone.reinterpret_coordentry(False)

        # Record undisplaced symmetry for projection of diplaced point groups
        psi4.set_parent_symmetry(molecule.schoenflies_symbol())

        ndisp = len(displacements)

        # This version is pretty dependent on the reference geometry being last (as it is now)
        print(' %d displacements needed.' % ndisp)
        energies = []

        # S/R: Write instructions for sow/reap procedure to output file and reap input file
        if freq_mode == 'sow':
            instructionsO = """\n#    The frequency sow/reap procedure has been selected through mode='sow'. In addition\n"""
            instructionsO += """#    to this output file (which contains no quantum chemical calculations), this job\n"""
            instructionsO += """#    has produced a number of input files (FREQ-*.in) for individual components\n"""
            instructionsO += """#    and a single input file (FREQ-master.in) with a frequency(mode='reap') command.\n"""
            instructionsO += """#    These files may look very peculiar since they contain processed and pickled python\n"""
            instructionsO += """#    rather than normal input. Follow the instructions below (repeated in FREQ-master.in)\n"""
            instructionsO += """#    to continue.\n#\n"""
            instructionsO += """#    Alternatively, a single-job execution of the hessian may be accessed through\n"""
            instructionsO += """#    the frequency wrapper option mode='continuous'.\n#\n"""
            psi4.print_out(instructionsO)

            instructionsM = """\n#    Follow the instructions below to carry out this frequency computation.\n#\n"""
            instructionsM += """#    (1)  Run all of the FREQ-*.in input files on any variety of computer architecture.\n"""
            instructionsM += """#       The output file names must be as given below (these are the defaults when executed\n"""
            instructionsM += """#       as `psi4 FREQ-1.in`, etc.).\n#\n"""
            for rgt in range(ndisp):
                pre = 'FREQ-' + str(rgt + 1)
                instructionsM += """#             psi4 -i %-27s -o %-27s\n""" % (pre + '.in', pre + '.out')
            instructionsM += """#\n#    (2)  Gather all the resulting output files in a directory. Place input file\n"""
            instructionsM += """#         FREQ-master.in into that directory and run it. The job will be minimal in\n"""
            instructionsM += """#         length and give summary results for the frequency computation in its output file.\n#\n"""
            instructionsM += """#             psi4 -i %-27s -o %-27s\n#\n\n""" % ('FREQ-master.in', 'FREQ-master.out')

            with open('FREQ-master.in', 'wb') as fmaster:
                fmaster.write('# This is a psi4 input file auto-generated from the hessian() wrapper.\n\n'.encode('utf-8'))
                fmaster.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
                fmaster.write(p4util.format_options_for_input(moleculeclone, **kwargs))
                p4util.format_kwargs_for_input(fmaster, lmode=2, return_wfn=True, **kwargs)
                fmaster.write(("""retE, retwfn = %s('%s', **kwargs)\n\n""" % (frequency.__name__, lowername)).encode('utf-8'))
                fmaster.write(instructionsM.encode('utf-8'))
            psi4.print_out(instructionsM)

        for n, displacement in enumerate(displacements):
            rfile = 'FREQ-%s' % (n + 1)

            # Build string of title banner
            banners = ''
            banners += """psi4.print_out('\\n')\n"""
            banners += """p4util.banner(' Hessian Computation: Energy Displacement %d ')\n""" % (n + 1)
            banners += """psi4.print_out('\\n')\n\n"""

            if freq_mode == 'continuous':

                # print progress to file and screen
                psi4.print_out('\n')
                p4util.banner('Loading displacement %d of %d' % (n + 1, ndisp))
                print(""" %d""" % (n + 1), end=('\n' if (n + 1 == ndisp) else ''))
                sys.stdout.flush()

                # Load in displacement into the active molecule
                moleculeclone.set_geometry(displacement)

                # Perform the energy calculation
                E, wfn = func(lowername, return_wfn=True, molecule=moleculeclone, **kwargs)
                energies.append(psi4.get_variable('CURRENT ENERGY'))

                # clean may be necessary when changing irreps of displacements
                psi4.clean()

            # S/R: Write each displaced geometry to an input file
            elif freq_mode == 'sow':
                moleculeclone.set_geometry(displacement)

                # S/R: Prepare molecule, options, kwargs, function call and energy save
                with open('%s.in' % (rfile), 'wb') as freagent:
                    freagent.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n')
                    freagent.write(p4util.format_molecule_for_input(moleculeclone, forcexyz=True).encode('utf-8'))
                    freagent.write(p4util.format_options_for_input(moleculeclone, **kwargs).encode('utf-8'))
                    p4util.format_kwargs_for_input(freagent, **kwargs)
                    freagent.write("""electronic_energy = %s('%s', **kwargs)\n\n""" % (func.__name__, lowername))
                    freagent.write("""psi4.print_out('\\nHESSIAN RESULT: computation %d for item %d """ % (os.getpid(), n + 1))
                    freagent.write("""yields electronic energy %20.12f\\n' % (electronic_energy))\n\n""")

            # S/R: Read energy from each displaced geometry output file and save in energies array
            elif freq_mode == 'reap':
                exec(banners)
                psi4.set_variable('NUCLEAR REPULSION ENERGY', moleculeclone.nuclear_repulsion_energy())
                energies.append(p4util.extract_sowreap_from_output(rfile, 'HESSIAN', n, freq_linkage, True))

        # S/R: Quit sow after writing files. Initialize skeleton wfn to receive grad for reap
        if freq_mode == 'sow':
            optstash.restore()
            if return_wfn:
                return (None, None)
            else:
                return None
        elif freq_mode == 'reap':
        #    psi4.set_variable('CURRENT ENERGY', energies[-1])
            wfn = psi4.new_wavefunction(molecule, psi4.get_global_option('BASIS'))

        # Assemble Hessian from energies
        H = psi4.fd_freq_0(molecule, energies, irrep)
        wfn.set_hessian(H)
        wfn.set_frequencies(psi4.get_frequencies())

        # The last item in the list is the reference energy, return it
        psi4.set_variable('CURRENT ENERGY', energies[-1])

        psi4.set_parent_symmetry('')
        optstash.restore()

        if return_wfn:
            return (wfn.hessian(), wfn)
        else:
            return wfn.hessian()


def frequency(name, **kwargs):
    r"""Function to compute harmonic vibrational frequencies.

    :aliases: frequencies(), freq()

    :returns: *float* |w--w| Total electronic energy in Hartrees.

    :returns: (*float*, :ref:`Wavefunction<sec:psimod_Wavefunction>`) |w--w| energy and wavefunction when **return_wfn** specified.

    .. note:: Analytic hessians are not available. Frequencies will proceed through
        finite differences according to availability of gradients or energies.

    .. _`table:freq_gen`:

    :type name: string
    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the system.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :type return_wfn: :ref:`boolean <op_py_boolean>`
    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|

        Indicate to additionally return the :ref:`Wavefunction<sec:psimod_Wavefunction>`
        calculation result as the second element (after *float* energy) of a tuple.
        Arrays of frequencies and the Hessian can be accessed through the wavefunction.

    :type dertype: :ref:`dertype <op_py_dertype>`
    :param dertype: |dl| ``'hessian'`` |dr| || ``'gradient'`` || ``'energy'``

        Indicates whether analytic (if available- they're not), finite
        difference of gradients (if available) or finite difference of
        energies is to be performed.

    :type mode: string
    :param mode: |dl| ``'continuous'`` |dr| || ``'sow'`` || ``'reap'``

        For a finite difference of energies or gradients frequency, indicates
        whether the calculations required to complete the frequency are to be run
        in one file (``'continuous'``) or are to be farmed out in an
        embarrassingly parallel fashion (``'sow'``/``'reap'``)/ For the latter,
        run an initial job with ``'sow'`` and follow instructions in its output file.
        For maximum flexibility, ``return_wfn`` is always on in ``'reap'`` mode.

    :type irrep: int or string
    :param irrep: |dl| ``-1`` |dr| || ``1`` || ``'b2'`` || ``'App'`` || etc.

        Indicates which symmetry block (:ref:`Cotton <table:irrepOrdering>` ordering) of vibrational
        frequencies to be computed. ``1``, ``'1'``, or ``'a1'`` represents
        :math:`a_1`, requesting only the totally symmetric modes.
        ``-1`` indicates a full frequency calculation.

    :examples:

    >>> # [1] Frequency calculation for all modes through highest available derivatives
    >>> frequency('ccsd')

    >>> # [2] Frequency calculation for b2 modes through finite difference of gradients
    >>> #     printing lowest mode frequency to screen and Hessian to output
    >>> E, wfn = frequencies('scf', dertype=1, irrep=4, return_wfn=True)
    >>> print wfn.frequencies().get(0, 0)
    >>> wfn.hessian().print_out()

    >>> # [3] Frequency calculation at default conditions and Hessian reuse at STP
    >>> E, wfn = freq('mp2', return_wfn=True)
    >>> set t 273.15
    >>> set p 100000
    >>> thermo(wfn, wfn.frequencies())

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)

    # are we in sow/reap mode?
    freq_mode = kwargs.get('mode', 'continuous').lower()
    if freq_mode not in ['continuous', 'sow', 'reap']:
        raise ValidationError("""Frequency execution mode '%s' not valid.""" % (freq_mode))

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', psi4.get_active_molecule())
    molecule.update_geometry()

    # Compute the hessian
    H, wfn = hessian(lowername, return_wfn=True, molecule=molecule, **kwargs)

    # S/R: Quit after getting new displacements
    if freq_mode == 'sow':
        return 0.0

    wfn.frequencies().print_out()
    psi4.thermo(wfn, wfn.frequencies())

    for postcallback in hooks['frequency']['post']:
        postcallback(lowername, **kwargs)

    if return_wfn:
        return (psi4.get_variable('CURRENT ENERGY'), wfn)
    else:
        return psi4.get_variable('CURRENT ENERGY')


def gdma(wfn, datafile=""):
    """Function to write wavefunction information in *wfn* to *filename* in
    molden format.

    .. versionadded:: 0.6

    :returns: None

    :type datafile: string
    :param datafile: optional control file (see GDMA manual) to peform more complicated DMA
                     analyses.  If this option is used, the File keyword must be set to read
                     a filename.fchk, where filename is provided by |globals__writer_file_label| .

    :type wfn: :ref:`Wavefunction<sec:psimod_Wavefunction>`
    :param wfn: set of molecule, basis, orbitals from which to generate DMA analysis

    :examples:

    >>> # [1] DMA analysis from MP2 wavefunction.  N.B. gradient must be requested to generate MP2 density.
    >>> grad, wfn = gradient('mp2', return_wfn=True)
    >>> gdma(wfn)

    """

    # Start by writing a G* checkpoint file, for the GDMA code to read in
    fw = psi4.FCHKWriter(wfn)
    molname = wfn.molecule().name()
    prefix = psi4.get_writer_file_prefix(molname)
    fchkfile = prefix + '.fchk'
    fw.write(fchkfile)

    if datafile:
        commands = datafile
    else:
        densname = wfn.name()
        if densname == "DFT":
            densname = "SCF"
        commands = 'psi4_dma_datafile.dma'
        radii = psi4.get_option('GDMA', 'GDMA_RADIUS')
        origin = psi4.get_option('GDMA', 'GDMA_ORIGIN')
        with open(commands, 'w') as f:
            f.write("File %s Density %s\n" % (fchkfile, densname))
            f.write("Angstrom\n")
            f.write("%s\n" % psi4.get_option('GDMA', 'GDMA_MULTIPOLE_UNITS'))
            f.write("Multipoles\n")
            if origin:
                try:
                    f.write("Origin %f %f %f\n" % (float(origin[0]), float(origin[1]), float(origin[2])))
                except:
                    raise ValidationError("The GDMA origin array should contain three entries: x, y, and z.")
            f.write("Switch %f\n" % psi4.get_option('GDMA', 'GDMA_SWITCH'))
            if radii:
                f.write("Radius %s\n" % " ".join([str(r) for r in radii]))
            f.write("Limit %d\n" % psi4.get_option('GDMA', 'GDMA_LIMIT') )
            f.write("Start\n")
            f.write("Finish\n")
    psi4.run_gdma(wfn, commands)

    os.remove(fchkfile)
    # If we generated the DMA control file, we should clean up here
    if not datafile:
        os.remove(commands)


def fchk(wfn, filename):
    """Function to write wavefunction information in *wfn* to *filename* in
    Gaussian FCHK format.

    .. versionadded:: 0.6

    :returns: None

    :type filename: string
    :param filename: destination file name for FCHK file

    :type wfn: :ref:`Wavefunction<sec:psimod_Wavefunction>`
    :param wfn: set of molecule, basis, orbitals from which to generate fchk file

    :examples:

    >>> # [1] FCHK file for DFT calculation
    >>> E, wfn = energy('b3lyp', return_wfn=True)
    >>> fchk(wfn, 'mycalc.fchk')

    """
    fw = psi4.FCHKWriter(wfn)
    fw.write(filename)


def molden(wfn, filename):
    """Function to write wavefunction information in *wfn* to *filename* in
    molden format.

    .. versionadded:: 0.5
       *wfn* parameter passed explicitly

    :returns: None

    :type filename: string
    :param filename: destination file name for MOLDEN file

    :type wfn: :ref:`Wavefunction<sec:psimod_Wavefunction>`
    :param wfn: set of molecule, basis, orbitals from which to generate cube files

    :examples:

    >>> # [1] Molden file for DFT calculation
    >>> E, wfn = energy('b3lyp', return_wfn=True)
    >>> molden(wfn, 'mycalc.molden')

    """
    try:
        occa = wfn.occupation_a()
        occb = wfn.occupation_a()
    except AttributeError:
        psi4.print_out("\n!Molden warning: This wavefunction does not have occupation numbers.\n"
                       "Writing zero's for occupation numbers\n\n")
        occa = psi4.Vector(wfn.nmopi())
        occb = psi4.Vector(wfn.nmopi())

    # At this point occupation number will be difficult to build, lets set them to zero
    mw = psi4.MoldenWriter(wfn)
    mw.write(filename, wfn.Ca(), wfn.Cb(), wfn.epsilon_a(), wfn.epsilon_b(), occa, occb)

def parse_cotton_irreps(irrep, point_group):
    r"""Function to return validated Cotton ordering index for molecular
    *point_group* from string or integer irreducible representation *irrep*.

    """
    cotton = {
        'c1': {
            'a': 1,
            '1': 1
        },
        'ci': {
            'ag': 1,
            'au': 2,
            '1': 1,
            '2': 2
        },
        'c2': {
            'a': 1,
            'b': 2,
            '1': 1,
            '2': 2
        },
        'cs': {
            'ap': 1,
            'app': 2,
            '1': 1,
            '2': 2
        },
        'd2': {
            'a': 1,
            'b1': 2,
            'b2': 3,
            'b3': 4,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4
        },
        'c2v': {
            'a1': 1,
            'a2': 2,
            'b1': 3,
            'b2': 4,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4
        },
        'c2h': {
            'ag': 1,
            'bg': 2,
            'au': 3,
            'bu': 4,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4,
        },
        'd2h': {
            'ag': 1,
            'b1g': 2,
            'b2g': 3,
            'b3g': 4,
            'au': 5,
            'b1u': 6,
            'b2u': 7,
            'b3u': 8,
            '1': 1,
            '2': 2,
            '3': 3,
            '4': 4,
            '5': 5,
            '6': 6,
            '7': 7,
            '8': 8
        }
    }

    try:
        return cotton[point_group.lower()][str(irrep).lower()]
    except KeyError:
        raise ValidationError("""Irrep '%s' not valid for point group '%s'.""" % (str(irrep), point_group))


# Aliases
opt = optimize
freq = frequency
frequencies = frequency
prop = property
