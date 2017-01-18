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

#!/usr/bin/python
# -*- python -*-
# -*- coding: utf-8 -*-
# vim:filetype=python:
# Create CMakeLists.txt template for leaf directories 
# (c) Roberto Di Remigio  <roberto.d.remigio@uit.no>
# licensed under the GNU Lesser General Public License

import os
import sys

sys.path.append('../cmake')
import argparse

parser = argparse.ArgumentParser(description='Create CMakeLists.txt template')
parser.add_argument('--labels',
        action='store',
        default=None,
        help='Labels for the test',
        metavar='STRING')

args = parser.parse_args() 
testname = os.path.basename(os.getcwd())
labels   = args.labels

for root, dirs, filenames in os.walk(os.getcwd()):
   for f in filenames:
       f = open('CMakeLists.txt', 'w')
       f.write('include(TestingMacros)\n\n')
       f.write('add_regression_test(' + testname + ' \"' + labels + '\")\n')
       f.close()

print('Template for {} created'.format(testname))

# vim:et:ts=4:sw=4
