#!/usr/bin/python

#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
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
# Extracts last git changelog item into html for psicode

gitlast = 'feed/latest_trac_changeset.txt'  # path to output file with latest changeset
githist = 'feed/history_trac_changeset.txt'  # path to output file with last Nhist changesets
Nhist = 100

# So that this can be called from anywhere
writedir = sys.argv[1] if len(sys.argv) == 2 else './'
thisdir = os.path.abspath(os.path.dirname(__file__))
os.chdir(thisdir)

log_cmd = 'git log --branches=master -%s --no-merges --pretty=format:\'<title> %%cr [%%h] by %%cn: <br/> %%s<br/></title>\' > %s'\
          % (Nhist, writedir + os.sep + githist)
os.system(log_cmd)

log_cmd = 'git log --branches=master -1 --no-merges --pretty=format:\'<title> %%cr by %%cn: <br/> %%s<br/></title>\' > %s'\
          % (writedir + os.sep + gitlast)
os.system(log_cmd)
