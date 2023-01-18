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

import re

yes = re.compile(r'^(yes|true|on|1)', re.IGNORECASE)
no = re.compile(r'^(no|false|off|0)', re.IGNORECASE)
der0th = re.compile(r'^(0|none|energy)', re.IGNORECASE)
der1st = re.compile(r'^(1|first|gradient)', re.IGNORECASE)
der2nd = re.compile(r'^(2|second|hessian)', re.IGNORECASE)
der3rd = re.compile(r'^(3|third)', re.IGNORECASE)
der4th = re.compile(r'^(4|fourth)', re.IGNORECASE)
der5th = re.compile(r'^(5|fifth)', re.IGNORECASE)


def parse_dertype(dertype, max_derivative=2):
    """Apply generous regex to `dertype` to return regularized integer and driver values for derivative level.

    Parameters
    ----------
    dertype : int or str
        Interpretable as a derivative level, regardless of case or type.
    max_derivative : int, optional
        Derivative level above which should throw FeatureNotImplemented error.

    Returns
    -------
    (int, {'energy', 'gradient', 'hessian'})
        Returns dertype as an integer and a driver-valid string.

    """
    derdriver = dict(enumerate(['energy', 'gradient', 'hessian', 'third', 'fourth', 'fifth']))

    if der0th.match(str(dertype)):
        derint = 0
    elif der1st.match(str(dertype)):
        derint = 1
    elif der2nd.match(str(dertype)):
        derint = 2
    elif der3rd.match(str(dertype)):
        derint = 3
    elif der4th.match(str(dertype)):
        derint = 4
    elif der5th.match(str(dertype)):
        derint = 5
    else:
        raise ValidationError("""Requested derivative level ({}) not recognized.""".format(dertype))

    if derint > max_derivative:
        raise FeatureNotImplemented("""derivative level ({})""".format(derint))

    return (derint, derdriver[derint])
