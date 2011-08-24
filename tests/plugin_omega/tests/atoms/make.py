#!/usr/bin/env python
import os;

Atoms = ["He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"]; 
S = [0,1,0,1,2,3,2,1,0,1,0,1,2,3,2,1,0,1,0,1,2,3,6,5,4,3,2,1,0,1,2,3,2,1,0];

for k in range(0,len(Atoms)):
    print ' Building input for %2s: %d' % (Atoms[k],S[k]+1);

    fh = open('%s.dat' % (Atoms[k]), 'w'); 

    fh.write('molecule {\n')
    fh.write('0 %d\n' %(S[k] + 1)) 
    fh.write('%s\n' %(Atoms[k]))
    fh.write('symmetry c1\n')
    fh.write('}\n')
    fh.write('\n')
    fh.write('plugin_load("../../plugin_omega.so")\n')
    fh.write('\n')
    fh.write('set globals {\n')
    fh.write('  scf_type direct\n')
    fh.write('  basis cc-pvdz\n')
    fh.write('  ri_basis_scf cc-pvdz-ri\n')
    fh.write('  reference uks\n')
    fh.write('  dft_functional wB97\n')
    fh.write('  dft_order_spherical 31 \n')
    fh.write('  dft_n_radial 50\n')
    fh.write('  omega_procedure ip\n')
    fh.write('  d_converge 5\n')
    fh.write('  e_converge 7\n')
    fh.write('  #debug 3\n')
    fh.write('}\n')
    fh.write('\n')
    fh.write('\n')
    fh.write("energy('scf')\n")
    fh.write('plugin("../../plugin_omega.so")\n')

    fh.close();

    os.system('$PSIPLUGIN %s.dat %s.out' %(Atoms[k], Atoms[k]))
