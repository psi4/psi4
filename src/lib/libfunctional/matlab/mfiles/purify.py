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

null_re = re.compile(r'^  t0 = 0.0;')

definition = [];

for partial in partials:
    if (not os.path.exists(partial)):
        continue;
    
    fh = open(partial);
    lines = fh.readlines();
    fh.close(); 

    os.system('rm ' + partial)

    if (re.match(null_re, lines[0])):
        continue;

    for k in range(len(lines) - 1):
        lines[k] = 'double ' + lines[k][2:];
   
    lines[-1] = partial + '[Q] += scale * ' + lines[-1][7:]; 

    definition.append('// ' + partial + '\n') 
    definition.append('{\n')
    for line in lines:
        definition.append(line)
    definition.append('}\n')
    definition.append('\n')

fh = open(target, 'w')
for line in definition:
    fh.write(line)
fh.close() 
