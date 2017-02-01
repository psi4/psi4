#!/usr/bin/env python

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

import sys, os

sys.path.append('%s/fsapt' % os.environ['PSIDATADIR'])
from fsapt import *

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
    stuffA = computeFsapt(dirA, False)
    print('   => Full Analysis <=\n')
    printOrder2(stuffA['order2'], stuffA['fragkeys']) 
    print('   => Reduced Analysis <=\n')
    printOrder2(stuffA['order2r'], stuffA['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    fh = open('%s/fsapt.dat' % dirB, 'w')
    fh, sys.stdout = sys.stdout, fh
    print('  ==> F-ISAPT: Links by Charge <==\n')
    stuffB = computeFsapt(dirB, False)
    print('   => Full Analysis <=\n')
    printOrder2(stuffB['order2'], stuffB['fragkeys']) 
    print('   => Reduced Analysis <=\n')
    printOrder2(stuffB['order2r'], stuffB['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    fh = open('%s/fsapt.dat' % dirD, 'w')
    fh, sys.stdout = sys.stdout, fh
    print('  ==> F-ISAPT: Links by Charge <==\n')
    order2D = diffOrder2(stuffA['order2r'], stuffB['order2r'])
    print('   => Reduced Analysis <=\n')
    printOrder2(order2D, stuffB['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    fh = open('%s/fsapt.dat' % dirA, 'a')
    fh, sys.stdout = sys.stdout, fh
    print('  ==> F-ISAPT: Links 50-50 <==\n')
    stuffA = computeFsapt(dirA, True)
    print('   => Full Analysis <=\n')
    printOrder2(stuffA['order2'], stuffA['fragkeys']) 
    print('   => Reduced Analysis <=\n')
    printOrder2(stuffA['order2r'], stuffA['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    fh = open('%s/fsapt.dat' % dirB, 'a')
    fh, sys.stdout = sys.stdout, fh
    print('  ==> F-ISAPT: Links 50-50 <==\n')
    stuffB = computeFsapt(dirB, True)
    print('   => Full Analysis <=\n')
    printOrder2(stuffB['order2'], stuffB['fragkeys']) 
    print('   => Reduced Analysis <=\n')
    printOrder2(stuffB['order2r'], stuffB['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    fh = open('%s/fsapt.dat' % dirD, 'a')
    fh, sys.stdout = sys.stdout, fh
    print('  ==> F-ISAPT: Links 50-50 <==\n')
    order2D = diffOrder2(stuffA['order2r'], stuffB['order2r'])
    print('   => Reduced Analysis <=\n')
    printOrder2(order2D, stuffB['fragkeysr']) 
    fh, sys.stdout = sys.stdout, fh
    fh.close()

    # > Order-1 PBD Files < #

    pdbA = PDB.fromGeom(stuffA['geom'])
    printOrder1(dirA, stuffA['order2r'], pdbA, stuffA['frags'])

    pdbB = PDB.fromGeom(stuffB['geom'])
    printOrder1(dirB, stuffB['order2r'], pdbB, stuffB['frags'])

    # Using A geometry
    printOrder1(dirD, order2D, pdbA, stuffA['frags'])



