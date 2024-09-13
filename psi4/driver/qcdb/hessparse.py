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

import numpy as np

import qcelemental as qcel


def load_hessian(shess, dtype):

    # list o'lines w/o comments or blanks
    shess = qcel.util.filter_comments(shess)
    lhess = list(filter(None, map(str.strip, shess.splitlines())))

    if dtype in ['fcmfinal', 'cfour']:
        nat = int(lhess[0].split()[0])
        ndof = 3 * nat
        datastr = '\n'.join(lhess[1:])
        nhess = np.fromstring(datastr, sep=' ')
        nhess = nhess.reshape(ndof, ndof)
    else:
        raise ValidationError('Unknown dtype: {}'.format(dtype))

    return nhess


def to_string(hess, handle, dtype='psi4'):
    """Writes Hessian in various formats.

    Parameters
    ----------
    hess : ndarray
        (3 * nat, 3 * nat) Hessian array.
    handle : filename or file handle
        If file handle, it must have been opened in binary, `wb`.
    dtype : {'fcmfinal', 'cfour', 'psi4', 'intder'}, optional
        Format to write Hessian.

    Returns
    -------
    None

    """
    nat = hess.shape[0] // 3
    assert hess.shape == (3 * nat, 3 * nat)

    if dtype in ['fcmfinal', 'cfour', 'psi4', 'intder']:
        second_number = (6 * nat) if dtype == 'intder' else (3 * nat)
        header = '{:5}{:5}'.format(nat, second_number)

        np.savetxt(handle, hess.reshape((-1, 3)), fmt='%20.10f', delimiter='', newline='\n', header=header, comments='')

        # Bounty! a Psi4 mug or similar gear to anyone who trace the `6 * nat` above to a pre-PSI/CCQC source.
        #   See discussion starting https://github.com/psi4/psi4/pull/953#issuecomment-381447849
    else:
        raise ValidationError('Unknown dtype: {}'.format(dtype))
