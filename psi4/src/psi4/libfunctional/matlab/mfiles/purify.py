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

#!/usr/bin/env python
import os, re, sys

target = sys.argv[1]    

partials = [
    'v',
    'v_rho_a',
    'v_rho_b',
    'v_gamma_aa',
    'v_gamma_ab',
    'v_gamma_bb',
    'v_tau_a',
    'v_tau_b',
    'v_rho_a_rho_a',
    'v_rho_a_rho_b',
    'v_rho_b_rho_b',
    'v_gamma_aa_gamma_aa',
    'v_gamma_aa_gamma_ab',
    'v_gamma_aa_gamma_bb',
    'v_gamma_ab_gamma_ab',
    'v_gamma_ab_gamma_bb',
    'v_gamma_bb_gamma_bb',
    'v_tau_a_tau_a',
    'v_tau_a_tau_b',
    'v_tau_b_tau_b',
    'v_rho_a_gamma_aa',
    'v_rho_a_gamma_ab',
    'v_rho_a_gamma_bb',
    'v_rho_b_gamma_aa',
    'v_rho_b_gamma_ab',
    'v_rho_b_gamma_bb',
    'v_rho_a_tau_a',
    'v_rho_a_tau_b',
    'v_rho_b_tau_a',
    'v_rho_b_tau_b',
    'v_gamma_aa_tau_a',
    'v_gamma_aa_tau_b',
    'v_gamma_ab_tau_a',
    'v_gamma_ab_tau_b',
    'v_gamma_bb_tau_a',
    'v_gamma_bb_tau_b',
    ];
deriv_table = {
    'v': 0,
    'v_rho_a': 1,
    'v_rho_b': 1,
    'v_gamma_aa': 1,
    'v_gamma_ab': 1,
    'v_gamma_bb': 1,
    'v_tau_a': 1,
    'v_tau_b': 1,
    'v_rho_a_rho_a': 2,
    'v_rho_a_rho_b': 2,
    'v_rho_b_rho_b': 2,
    'v_gamma_aa_gamma_aa': 2,
    'v_gamma_aa_gamma_ab': 2,
    'v_gamma_aa_gamma_bb': 2,
    'v_gamma_ab_gamma_ab': 2,
    'v_gamma_ab_gamma_bb': 2,
    'v_gamma_bb_gamma_bb': 2,
    'v_tau_a_tau_a': 2,
    'v_tau_a_tau_b': 2,
    'v_tau_b_tau_b': 2,
    'v_rho_a_gamma_aa': 2,
    'v_rho_a_gamma_ab': 2,
    'v_rho_a_gamma_bb': 2,
    'v_rho_b_gamma_aa': 2,
    'v_rho_b_gamma_ab': 2,
    'v_rho_b_gamma_bb': 2,
    'v_rho_a_tau_a': 2,
    'v_rho_a_tau_b': 2,
    'v_rho_b_tau_a': 2,
    'v_rho_b_tau_b': 2,
    'v_gamma_aa_tau_a': 2,
    'v_gamma_aa_tau_b': 2,
    'v_gamma_ab_tau_a': 2,
    'v_gamma_ab_tau_b': 2,
    'v_gamma_bb_tau_a': 2,
    'v_gamma_bb_tau_b': 2,
    };

null_re = re.compile(r'^  t0 = 0.0;')

definition = [];

for partial in partials:
    if (not os.path.exists(partial)):
        continue;

    deriv = deriv_table[partial]
    
    fh = open(partial);
    lines = fh.readlines();
    fh.close(); 

    os.system('rm ' + partial)

    if (re.match(null_re, lines[0])):
        continue;

    for k in range(len(lines) - 1):
        lines[k] = 'double ' + lines[k][2:];
   
    lines[-1] = partial + '[Q] += scale * (' + lines[-1][7:-2] + ');\n'; 

    definition.append('// ' + partial + '\n') 
    definition.append('if (deriv >= %d) {\n' % (deriv))
    for line in lines:
        definition.append('    ' + line)
    definition.append('}\n')
    definition.append('\n')

fh = open(target, 'w')
for line in definition:
    fh.write(line)
fh.close() 