#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
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
from .metadata import __version__, version_formatter  # noqa: E401

time_string = datetime.datetime.now().strftime('%A, %d %B %Y %I:%M%p')
pid = os.getpid()


def sizeof_fmt(num, suffix='B'):
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)


def print_header():
    driver_info = version_formatter("""{version} {release}""")
    git_info = version_formatter("""{{{branch}}} {githash} {clean}""")
    datadir = core.get_environment("PSIDATADIR")
    memory = sizeof_fmt(core.get_memory())
    threads = str(core.nthread())

    header = """
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 %s

                         Git: Rev %s


    J. M. Turney, A. C. Simmonett, R. M. Parrish, E. G. Hohenstein,
    F. A. Evangelista, J. T. Fermann, B. J. Mintz, L. A. Burns, J. J. Wilke,
    M. L. Abrams, N. J. Russ, M. L. Leininger, C. L. Janssen, E. T. Seidl,
    W. D. Allen, H. F. Schaefer, R. A. King, E. F. Valeev, C. D. Sherrill,
    and T. D. Crawford, WIREs Comput. Mol. Sci. 2, 556-565 (2012)
    (doi: 10.1002/wcms.93)


                         Additional Contributions by
    A. E. DePrince, U. Bozkaya, A. Yu. Sokolov, D. G. A. Smith, R. Di Remigio,
    R. M. Richard, J. F. Gonthier, H. R. McAlexander, M. Saitow, and
    B. P. Pritchard
    -----------------------------------------------------------------------


    Psi4 started on: %s

    Process ID: %6d
    PSIDATADIR: %s
    Memory:     %s
    Threads:    %s
    """ % (driver_info, git_info, time_string, pid, datadir, memory, threads)
    core.print_out(header)
