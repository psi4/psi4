#!/usr/bin/python
# -*- python -*-
# -*- coding: utf-8 -*-
# vim:filetype=python:
# Create CMakeLists.txt template for leaf directories 
# (c) Roberto Di Remigio  <roberto.d.remigio@uit.no>
# licensed under the GNU Lesser General Public License

import sys

sys.path.append('../cmake')
import argparse

parser = argparse.ArgumentParser(description='Create CMakeLists.txt template')
parser.add_argument('testname',
                     action='store',
                     metavar='TESTNAME',
                     nargs='?',
                     help='Name of the test to be created')
parser.add_argument('--labels',
        action='store',
        default=None,
        help='Labels for the test',
        metavar='STRING')

args = parser.parse_args()

args = parser.parse_args() 
testname = args.testname
labels   = args.labels

f = open('CMakeLists.txt', 'w')
f.write('include(TestingMacros)\n\n')
f.write('add_regression_test(' + testname + ' \"' + labels + '\")\n')
f.close()

print('Template for {} created'.format(args.testname))

# vim:et:ts=4:sw=4
