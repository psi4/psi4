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

import pydantic

from typing import Dict, List, Any, Union

from psi4.driver import p4util
from psi4.driver.driver_nbody import NBodyComputer
from psi4.driver.task_base import BaseTask


class NBodyComputerMultiLevel(BaseTask):
    driver: str
    molecule: Any

    # nbody multi-level kwargs
    levels: Dict[Any, str] = {}
    return_total_data: bool = False
    return_wfn: bool = False

    task_list: Dict[str, Any] = {}

    def __init__(self, **data):
        BaseTask.__init__(self, **data)

    def build_tasks(self, packets, **kwargs):
        supersystem = packets.pop('supersystem', None)

        # Create NBodyComputers for each requested truncation order
        for n, val in packets.items():
            self.task_list[n] = NBodyComputer(**val['packet'], **kwargs, max_nbody=n)
            val['packet'].pop('molecule', None)
            self.task_list[n].build_tasks(val['computer'], **val['packet'])

        # Create computers for calculations on the whole cluster if requested
        if supersystem is not None:
            self.task_list['supersystem'] = supersystem['computer'](**supersystem['packet'])
            # Subtract MBE result of the highest requested MBE truncation order
            n = max(packets)
            self.task_list[str(n) + '_supersystem'] = NBodyComputer(**supersystem['packet'], **kwargs, max_nbody=n)
            supersystem['packet'].pop('molecule', None)
            self.task_list[str(n) + '_supersystem'].build_tasks(supersystem['computer'], **supersystem['packet'])


    def plan(self):
        ret = []
        for k, v in self.task_list.items():
            ret.append(v.plan())
        return ret

    def compute(self):
        with p4util.hold_options_state():
            # gof = core.get_output_file()
            # core.close_outfile()
            for k, v in self.task_list.items():
                v.compute()

            # core.set_output_file(gof, True)

    def prepare_results(self):

        from psi4.driver.driver_nbody import _print_nbody_energy

        ptype = self.driver
        natoms = self.molecule.natom()
        bsse_type = list(self.task_list.values())[0].bsse_type
        max_nbody = max([i for i in self.task_list if isinstance(i, int)])
        supersystem = self.task_list.pop('supersystem', None)
        n_supersystem = self.task_list.pop(str(max_nbody) + '_supersystem', None)

        # Initialize with zeros
        energy_result, gradient_result, hessian_result = 0, None, None
        energy_body_contribution = {b: {} for b in bsse_type}
        energy_body_dict = {b: {} for b in bsse_type}
        if ptype in ['gradient', 'hessian']:
            gradient_result = np.zeros((natoms, 3))
        if ptype == 'hessian':
            hessian_result = np.zeros((natoms * 3, natoms * 3))

        for n in sorted(self.task_list)[::-1]:
            result_n = self.task_list[n]
            energy_bsse_dict = {b: 0 for b in result_n.bsse_type}
            results = result_n.prepare_results()

            for m in range(n - 1, n + 1):
                if m == 0: continue
                # Subtract the (n-1)-body contribution from the n-body contribution to get the n-body effect
                sign = (-1)**(1 - m // n)
                for b in result_n.bsse_type:
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

            energy_result += energy_bsse_dict[bsse_type[0]]
            for b in bsse_type:
                energy_body_contribution[b][n] = energy_bsse_dict[b]

        if supersystem:
            # Super system recovers higher order effects at a lower level
            supersystem_result = supersystem.get_results()
            components = n_supersystem.prepare_results()

            energy_result += supersystem_result['properties']['return_energy'] - components['energy_body_dict'][max_nbody]
            for b in bsse_type:
                energy_body_contribution[b][self.molecule.nfragments()] = (supersystem_result['properties']['return_energy'] -
                                                                           components['energy_body_dict'][max_nbody])

            if ptype == 'hessian':
                gradient_result += supersystem_result['psi4:qcvars']['CURRENT GRADIENT'] - components['gradient_body_dict'][max_nbody]
                hessian_result += supersystem_result['return_result'] - components['ptype_body_dict'][max_nbody]

            elif ptype == 'gradient':
                gradient_result += supersystem_result['return_result'] - components['ptype_body_dict'][max_nbody]


        for b in bsse_type:
            for n in energy_body_contribution[b]:
                energy_body_dict[b][n] = sum(
                    [energy_body_contribution[b][i] for i in range(1, n + 1) if i in energy_body_contribution[b]])

        is_embedded = list(self.task_list.values())[0].embedding_charges
        for b in bsse_type:
            _print_nbody_energy(energy_body_dict[b], '%s-corrected multilevel many-body expansion' % b.upper(),
                                is_embedded)

        if not self.return_total_data:
            # Remove monomer cotribution for interaction data
            energy_result -= energy_body_dict[bsse_type[0]][1]
            if ptype in ['gradient', 'hessian']:
                gradient_result -= gradient1
            if ptype == 'hessian':
                hessian_result -= hessian1


        energy_body_dict = {str(k) + b: v for b in energy_body_dict for k, v in energy_body_dict[b].items()}

        results = {'ret_energy': energy_result, 'ret_ptype': eval(ptype + '_result'),
                   'energy_body_dict': energy_body_dict, 'ptype_body_dict': {}}

        return results


    def get_results(self):
        # Use NBodyComputer functions as it has the same data structure
        return NBodyComputer.get_results(self)

    def get_psi_results(self, return_wfn=False):
        # Use NBodyComputer functions as it has the same data structure
       return NBodyComputer.get_psi_results(self, return_wfn)
