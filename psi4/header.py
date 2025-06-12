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

import datetime
import os
import socket

from . import core
from .metadata import __version__, version_formatter

if "undef" in __version__:
    raise TypeError(
        """Using custom build without tags. Please pull git tags with `git pull origin master --tags`. If building from source, `git fetch upstream "refs/tags/*:refs/tags/*"` and re-make."""
    )

time_string = datetime.datetime.now().strftime('%A, %d %B %Y %I:%M%p')
pid = os.getpid()

__version__ = '1.10a1.dev86'
__version_branch_name = 'releaseproc'
__version_cmake = '1.9.0.999'
__version_is_clean = 'False'
__version_last_release = '1.9'
__version_long = '1.10a1.dev86+7934858'
__version_prerelease = 'False'
__version_release = 'False'

__citation_reference = "J. Chem. Phys. 152(18) 184108 (2020)"
__citation_doi = "10.1063/5.0006002"
__citation_doilink = "https://doi.org/" + __citation_doi
__citation_title = "Psi4 1.4: Open-Source Software for High-Throughput Quantum Chemistry"
__citation_authors = "D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish, M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio, A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer, R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni, J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein, B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov, K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King, F. A. Evangelista, J. M. Turney, T. D. Crawford, and C. D. Sherrill"


def citation_formatter(formatstring='oneline'):
    """Return citation information string when supplied with *formatstring* suitable
    for ``formatstring.format()``.  Use plaintext and any placeholders among: doi,
    doilink, reference, title, authors. For example '{doi} {title}' returns something like
    '10.1063/5.0006002 Psi4 1.4: Open-Source Software for High-Throughput Quantum Chemistry'.

    """
    if formatstring == 'all':
        formatstring = '{authors}, "{title}", {reference}, {doilink}'
    elif formatstring == 'oneline':
        formatstring = f'Psi4 v{version_formatter()}' + ', {reference} ({doi})'

    ans = formatstring.format(doi=__citation_doi,
                              doilink=__citation_doilink,
                              reference=__citation_reference,
                              title=__citation_title,
                              authors=__citation_authors)
    return ans


def sizeof_fmt(num, suffix='B'):
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)

def print_header():
    driver_info = version_formatter("""{version} {release}""")
    git_info = version_formatter("""{{{branch}}} {githash} {clean}""")
    journal_doi = citation_formatter("""{reference}. {doilink}""")
    datadir = core.get_datadir()
    memory = sizeof_fmt(core.get_memory())
    hostname = socket.gethostname()
    threads = str(core.get_num_threads())

    header = f"""
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 {driver_info}

                         Git: Rev {git_info}


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    {journal_doi}

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, and D. L. Poole

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: {time_string}

    Process ID: {pid}
    Host:       {hostname}
    PSIDATADIR: {datadir}
    Memory:     {memory}
    Threads:    {threads}
    """
    core.print_out(header)
