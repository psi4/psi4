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

def prepare_results(self):
    """
    Use different levels of theory for different expansion levels
    See kwargs description in driver_nbody.nbody_gufunc
    :returns: *return type of func* |w--w| The data.
    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| data and wavefunction with energy/gradient/hessian set appropriately when **return_wfn** specified.
    """
    from psi4.driver.driver_nbody import _print_nbody_energy

    ptype = self.driver 
    natoms = self.molecule.natom()
    supersystem = {k: v for k, v in self.task_list.items() if k.startswith('supersystem')}

    # Initialize with zeros
    energy_result, gradient_result, hessian_result = 0, None, None
    energy_body_contribution = {b: {} for b in self.bsse_type}
    energy_body_dict = {b: {} for b in self.bsse_type}
    if ptype in ['gradient', 'hessian']:
        gradient_result = np.zeros((natoms, 3))
    if ptype == 'hessian':
        hessian_result = np.zeros((natoms * 3, natoms * 3))

    levels = {int(i.split('_')[0]) for i in self.task_list if 'supersystem' not in i}

    for n in sorted(levels)[::-1]:
        energy_bsse_dict = {b: 0 for b in self.bsse_type}
        self.max_nbody = n
        results = {k: v for k, v in self.task_list.items() if k.startswith(str(n))}
        results = self.prepare_results(results=results)

        for m in range(n - 1, n + 1):
            if m == 0: continue
            # Subtract the (n-1)-body contribution from the n-body contribution to get the n-body effect
            sign = (-1)**(1 - m // n)
            for b in self.bsse_type:
                energy_bsse_dict[b] += sign * results['%s_energy_body_dict' %b.lower()]['%i%s' %(m, b.lower())]

            if ptype == 'hessian':
                hessian_result += sign * results['ptype_body_dict'][m]
                gradient_result += sign * results['gradient_body_dict'][m]
                if n == 1:
                    hessian1 = results['ptype_body_dict'][n]
                    gradient1 = results['gradient_body_dict'][n]

            elif ptype == 'gradient':
                gradient_result += sign * results['ptype_body_dict'][m]
                # Keep 1-body contribution to compute interaction data
                if n == 1:
                    gradient1 = results['ptype_body_dict'][n]

        energy_result += energy_bsse_dict[self.bsse_type[0]]
        for b in self.bsse_type:
            energy_body_contribution[b][n] = energy_bsse_dict[b]

    if supersystem:
        # Super system recovers higher order effects at a lower level
        supersystem_result = supersystem.pop('supersystem_' + str(self.max_frag)).get_results()
        self.max_nbody = max(levels)
        components = self.prepare_results(results=supersystem)

        energy_result += supersystem_result['properties']['return_energy'] - components['energy_body_dict'][self.max_nbody]
        for b in self.bsse_type:
            energy_body_contribution[b][self.molecule.nfragments()] = (supersystem_result['properties']['return_energy'] -
                                                                       components['energy_body_dict'][self.max_nbody])

        if ptype == 'hessian':
            gradient_result += supersystem_result['extras']['qcvars']['CURRENT GRADIENT'] - components['gradient_body_dict'][self.max_nbody]
            hessian_result += supersystem_result['return_result'] - components['ptype_body_dict'][self.max_nbody]

        elif ptype == 'gradient':
            gradient_result += np.array(supersystem_result['return_result']).reshape((-1, 3)) - components['ptype_body_dict'][self.max_nbody]


    for b in self.bsse_type:
        for n in energy_body_contribution[b]:
            energy_body_dict[b][n] = sum(
                [energy_body_contribution[b][i] for i in range(1, n + 1) if i in energy_body_contribution[b]])

    is_embedded = self.embedding_charges
    for b in self.bsse_type:
        _print_nbody_energy(energy_body_dict[b], '%s-corrected multilevel many-body expansion' % b.upper(),
                            is_embedded)

    if not self.return_total_data:
        # Remove monomer cotribution for interaction data
        energy_result -= energy_body_dict[self.bsse_type[0]][1]
        if ptype in ['gradient', 'hessian']:
            gradient_result -= gradient1
        if ptype == 'hessian':
            hessian_result -= hessian1


    energy_body_dict = {str(k) + b: v for b in energy_body_dict for k, v in energy_body_dict[b].items()}

    return {'ret_energy': energy_result, 'ret_ptype': locals()[ptype + '_result'], 'energy_body_dict': energy_body_dict, 'ptype_body_dict': {}}
