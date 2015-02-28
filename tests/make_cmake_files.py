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
