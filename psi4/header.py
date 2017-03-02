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

import datetime
import os

from . import core
from .metadata import __version__, version_formatter

time_string = datetime.datetime.now().strftime('%A, %d %B %Y %I:%M%p')
pid = os.getpid()

def sizeof_fmt(num, suffix='B'):
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)

def print_header():
    driver_info = version_formatter("""{version} {release}""")
    git_info = version_formatter("""{{{branch}}} {githash} {clean}""")
    datadir = core.get_environment("PSIDATADIR")
    memory = sizeof_fmt(core.get_memory())
    threads = str(core.get_num_threads())

    header = """
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 %s

                         Git: Rev %s


    R. M. Parrish, L. A. Burns, D. G. A. Smith, A. C. Simmonett,
    A. E. DePrince III, E. G. Hohenstein, U. Bozkaya, A. Yu. Sokolov,
    R. Di Remigio, R. M. Richard, J. F. Gonthier, A. M. James,
    H. R. McAlexander, A. Kumar, M. Saitow, X. Wang, B. P. Pritchard,
    P. Verma, H. F. Schaefer III, K. Patkowski, R. A. King, E. F. Valeev,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, and C. D. Sherrill,
    submitted.

    -----------------------------------------------------------------------


    Psi4 started on: %s

    Process ID: %6d
    PSIDATADIR: %s
    Memory:     %s
    Threads:    %s
    """ % (driver_info, git_info, time_string, pid, datadir, memory, threads)
    core.print_out(header)
