# -*- python -*-
# -*- coding: utf-8 -*-
# vim:filetype=python:
#
# Roberto Di Remigio <roberto.d.remigio@uit.no>
# University of Tromso 2015
#

"""
This is a placeholder for the real pcmsolver.py script.
The location of the real pcmsolver.py is configured by CMake
to point to the proper install prefix.
In this way we avoid to transform inputparser.py into a file
that has to be configured by CMake
With conda, this starts to get complicated. Bottom option
works for build-in-place and build-psi-w-prebuilt-pcmsolver.
Upper option necessary when psi4metapackage is conda build
dependency. At least I think that's what's going on.

"""

if '' == 'ON':
    PCMSolver_PARSE_DIR = '/opt/anaconda1anaconda2anaconda3/bin'
else:
    PCMSolver_PARSE_DIR = '/theoryfs2/ds/richard/SrcFiles/psi4/build/stage/theoryfs2/ds/richard/SrcFiles/psi4/Install/external/pcmsolver/bin'

