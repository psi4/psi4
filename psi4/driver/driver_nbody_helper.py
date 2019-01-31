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
from psi4.driver.p4util.exceptions import *


def multi_level(func, **kwargs):
    """
    Use different levels of theory for different expansion levels
    See kwargs description in driver_nbody.nbody_gufunc

    :returns: *return type of func* |w--w| The data.

    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| data and wavefunction with energy/gradient/hessian set appropriately when **return_wfn** specified.

    """
    from psi4.driver.driver_nbody import nbody_gufunc
    from psi4.driver.driver_nbody import _print_nbody_energy

    ptype = kwargs['ptype']
    return_wfn = kwargs.get('return_wfn', False)
    kwargs['return_wfn'] = True
    levels = kwargs.pop('levels')
    for i in levels:
        if isinstance(i, str): levels[i.lower()] = levels.pop(i)
    supersystem = levels.pop('supersystem', False)
    molecule = kwargs.get('molecule', core.get_active_molecule())
    kwargs['bsse_type'] = [kwargs['bsse_type']] if isinstance(kwargs['bsse_type'], str) else kwargs['bsse_type']
    natoms = molecule.natom()

    # Initialize with zeros
    energy_result, gradient_result, hessian_result = 0, None, None
    energy_body_contribution = {b: {} for b in kwargs['bsse_type']}
    energy_body_dict = {b: {} for b in kwargs['bsse_type']}
    wfns = {}
    if ptype in ['gradient', 'hessian']:
        gradient_result = np.zeros((natoms, 3))
    if ptype == 'hessian':
        hessian_result = np.zeros((natoms * 3, natoms * 3))

    if kwargs.get('charge_method', False) and not kwargs.get('embedding_charges', False):
        kwargs['embedding_charges'] = compute_charges(kwargs['charge_method'],
                                      kwargs.get('charge_type', 'MULLIKEN_CHARGES').upper(), molecule)

    for n in sorted(levels)[::-1]:
        molecule.set_name('%i' %n)
        kwargs_copy = kwargs.copy()
        kwargs_copy['max_nbody'] = n
        energy_bsse_dict = {b: 0 for b in kwargs['bsse_type']}
        if isinstance(levels[n], str):
            # If a new level of theory is provided, compute contribution
            ret, wfn = nbody_gufunc(func, levels[n], **kwargs_copy)
            wfns[n] = wfn
        else:
            # For the n-body contribution, use available data from the higher order levels[n]-body
            wfn = wfns[levels[n]]

        for m in range(n - 1, n + 1):
            if m == 0: continue
            # Subtract the (n-1)-body contribution from the n-body contribution to get the n-body effect
            sign = (-1)**(1 - m // n)
            for b in kwargs['bsse_type']:
                energy_bsse_dict[b] += sign * wfn.variable('%i%s' % (m, b.lower()))
            if ptype in ['gradient', 'hessian']:
                gradient_result += sign * np.array(wfn.variable('GRADIENT ' + str(m)))
                # Keep 1-body contribution to compute interaction data
                if n == 1:
                    gradient1 = np.array(wfn.variable('GRADIENT ' + str(m)))
            if ptype == 'hessian':
                hessian_result += sign * np.array(wfn.variable('HESSIAN ' + str(m)))
                if n == 1:
                    hessian1 = np.array(wfn.variable('HESSIAN ' + str(m)))
        energy_result += energy_bsse_dict[kwargs['bsse_type'][0]]
        for b in kwargs['bsse_type']:
            energy_body_contribution[b][n] = energy_bsse_dict[b]

    if supersystem:
        # Super system recovers higher order effects at a lower level
        molecule.set_name('supersystem')
        kwargs_copy = kwargs.copy()
        kwargs_copy.pop('bsse_type')
        kwargs_copy.pop('ptype')
        ret, wfn_super = func(supersystem, **kwargs_copy)
        core.clean()
        kwargs_copy = kwargs.copy()
        kwargs_copy['bsse_type'] = 'nocp'
        kwargs_copy['max_nbody'] = max(levels)
        # Subtract lower order effects to avoid double counting
        ret, wfn = nbody_gufunc(func, supersystem, **kwargs_copy)
        energy_result += wfn_super.energy() - wfn.variable(str(max(levels)))
        for b in kwargs['bsse_type']:
            energy_body_contribution[b][molecule.nfragments()] = wfn_super.energy() - wfn.variable(
                str(max(levels)))

        if ptype in ['gradient', 'hessian']:
            gradient_result += np.array(wfn_super.gradient()) - np.array(wfn.variable('GRADIENT ' + str(max(levels))))
        if ptype == 'hessian':
            hessian_result += np.array(wfn_super.hessian()) - np.array(wfn.variable('HESSIAN ' + str(max(levels))))
        levels['supersystem'] = supersystem

    for b in kwargs['bsse_type']:
        for n in energy_body_contribution[b]:
            energy_body_dict[b][n] = sum(
                [energy_body_contribution[b][i] for i in range(1, n + 1) if i in energy_body_contribution[b]])

    is_embedded = kwargs.get('embedding_charges', False) or kwargs.get('charge_method', False)
    for b in kwargs['bsse_type']:
        _print_nbody_energy(energy_body_dict[b], '%s-corrected multilevel many-body expansion' % b.upper(),
                            is_embedded)

    if not kwargs['return_total_data']:
        # Remove monomer cotribution for interaction data
        energy_result -= energy_body_dict[kwargs['bsse_type'][0]][1]
        if ptype in ['gradient', 'hessian']:
            gradient_result -= gradient1
        if ptype == 'hessian':
            hessian_result -= hessian1
    wfn_out = core.Wavefunction.build(molecule, 'def2-svp')
    core.set_variable("CURRENT ENERGY", energy_result)
    wfn_out.set_variable("CURRENT ENERGY", energy_result)
    gradient_result = core.Matrix.from_array(gradient_result) if gradient_result is not None else None
    wfn_out.set_gradient(gradient_result)
    hessian_result = core.Matrix.from_array(hessian_result) if hessian_result is not None else None
    wfn_out.set_hessian(hessian_result)
    ptype_result = eval(ptype + '_result')
    for b in kwargs['bsse_type']:
        for i in energy_body_dict[b]:
            wfn_out.set_variable(str(i) + b, energy_body_dict[b][i])

    if kwargs['return_wfn']:
        return (ptype_result, wfn_out)
    else:
        return ptype_result


def compute_charges(charge_method, charge_type, molecule):
    """
    Compute charges for nbody fragments
    """
    from psi4.driver.driver import energy
    from psi4.driver.p4util.util import oeprop

    charges = {}
    molecule = molecule.clone()
    for i in range(1, molecule.nfragments() + 1):
        molecule.set_name('charges%i' %i)
        e, wfn = energy(charge_method, molecule=molecule.extract_subsets([i]), return_wfn=True)
        oeprop(wfn, charge_type)
        charges[i] = wfn.atomic_point_charges().np

    return charges


def electrostatic_embedding(metadata, pair):
    """
    Add atom-centered point charges for fragments whose basis sets are not included in the computation.
    """
    from psi4.driver import qmmm
    from psi4.driver import constants

    if not metadata['return_total_data']:
        raise Exception('Cannot return interaction data when using embedding scheme.')
    # Add embedding point charges
    Chrgfield = qmmm.QMMM()
    for p in metadata['embedding_charges']:
        if p in pair[1]: continue
        mol = metadata['molecule'].extract_subsets([p])
        for i in range(mol.natom()):
            geom = np.array([mol.x(i), mol.y(i), mol.z(i)])
            if mol.units() == 'Angstrom':
                geom *= constants.bohr2angstroms
            Chrgfield.extern.addCharge(metadata['embedding_charges'][p][i], geom[0], geom[1], geom[2])
    core.set_global_option_python('EXTERN', Chrgfield.extern)
