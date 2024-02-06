#!/usr/bin/env python

#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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

import os
import sys

import psi4

path = psi4.core.get_datadir()
sys.path.append('%s/fsapt' % path)
from fsapt import *  # isort:skip


# => Driver Code <= #

if __name__ == '__main__':

    # > Working Dirname < #

    if len(sys.argv) == 3:
        dirA = sys.argv[1]
        dirB = sys.argv[2]
        dirD = '.'
    elif len(sys.argv) == 4:
        dirA = sys.argv[1]
        dirB = sys.argv[2]
        dirD = sys.argv[3]
    else:
        raise Exception('Usage: fsapt.py dirnameA dirnameB [dirnameD]')

    # Make dirD if needed
    if not os.path.exists(dirD):
        os.makedirs(dirD)

    # > Order-2 Analysis < #

    fh = open('%s/fsapt.dat' % dirA, 'w')
    fh, sys.stdout = sys.stdout, fh
    print('  ==> F-ISAPT: Links by Charge <==\n')
    stuffA = compute_fsapt(dirA, False)
    print('   => Full Analysis <=\n')
    print_order2(stuffA['order2'], stuffA['fragkeys']) 
    print('   => Reduced Analysis <=\n')
    print_order2(stuffA['order2r'], stuffA['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    fh = open('%s/fsapt.dat' % dirB, 'w')
    fh, sys.stdout = sys.stdout, fh
    print('  ==> F-ISAPT: Links by Charge <==\n')
    stuffB = compute_fsapt(dirB, False)
    print('   => Full Analysis <=\n')
    print_order2(stuffB['order2'], stuffB['fragkeys']) 
    print('   => Reduced Analysis <=\n')
    print_order2(stuffB['order2r'], stuffB['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    fh = open('%s/fsapt.dat' % dirD, 'w')
    fh, sys.stdout = sys.stdout, fh
    print('  ==> F-ISAPT: Links by Charge <==\n')
    order2D = diff_order2(stuffA['order2r'], stuffB['order2r'])
    print('   => Reduced Analysis <=\n')
    print_order2(order2D, stuffB['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    fh = open('%s/fsapt.dat' % dirA, 'a')
    fh, sys.stdout = sys.stdout, fh
    print('  ==> F-ISAPT: Links 50-50 <==\n')
    stuffA = compute_fsapt(dirA, True)
    print('   => Full Analysis <=\n')
    print_order2(stuffA['order2'], stuffA['fragkeys']) 
    print('   => Reduced Analysis <=\n')
    print_order2(stuffA['order2r'], stuffA['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    fh = open('%s/fsapt.dat' % dirB, 'a')
    fh, sys.stdout = sys.stdout, fh
    print('  ==> F-ISAPT: Links 50-50 <==\n')
    stuffB = compute_fsapt(dirB, True)
    print('   => Full Analysis <=\n')
    print_order2(stuffB['order2'], stuffB['fragkeys']) 
    print('   => Reduced Analysis <=\n')
    print_order2(stuffB['order2r'], stuffB['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    fh = open('%s/fsapt.dat' % dirD, 'a')
    fh, sys.stdout = sys.stdout, fh
    print('  ==> F-ISAPT: Links 50-50 <==\n')
    order2D = diff_order2(stuffA['order2r'], stuffB['order2r'])
    print('   => Reduced Analysis <=\n')
    print_order2(order2D, stuffB['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    # > Order-1 PBD Files < #

    pdbA = PDB.from_geom(stuffA['geom'])
    print_order1(dirA, stuffA['order2r'], pdbA, stuffA['frags'])

    pdbB = PDB.from_geom(stuffB['geom'])
    print_order1(dirB, stuffB['order2r'], pdbB, stuffB['frags'])

    # Using A geometry
    print_order1(dirD, order2D, pdbA, stuffA['frags'])



