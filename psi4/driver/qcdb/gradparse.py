#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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


# TODO: Include gradient loading feature, analagously to hessian loading.


def to_string(grad, handle, dtype='file11', mol=None, energy=None):
    """Writes gradient in various formats.

    Parameters
    ----------
    grad : ndarray
        (nat, 3) gradient array.
    handle : filename or file handle
        If file handle, it must have been opened in binary, `wb`.
    dtype : {'file11', 'GRD', 'minimal'}, optional
        Format to write Hessian.
            file11 is an old file format with the number of atoms, last energy and geometry, as well as gradient.
            GRD is an ACESII file format with the number of atoms and geometry, as well as gradient.
            minimal simply gives the number of atoms and the gradient.
    mol : psi4.molecule or qcdb.molecule
        The molecule the gradient is from. Needed only in file11 formatting.
    energy : psi4.molecule or qcdb.molecule
        The electronic energy of the molecule the gradient is evaluated from. Needed only in file11 format.

    Returns
    -------
    None

    """
    # Validate input.
    try:
        grad = grad.reshape(-1, 3)
    except ValueError:
        raise ValidationError("Number of gradient entries must be divisible by 3. Received {} entries.".format(
            grad.size))

    if mol and grad.shape[0] != mol.natom():
        raise ValidationError("Number of gradient entries must be number of atoms * 3. Received {} entries.".format(
            grad.size))

    if dtype in ['file11', 'GRD']:
        if mol is None:
            raise ValidationError("molecule must both be defined to print a gradient in file11 format.")

    # Now actually write out the gradients.
    if dtype in ['file11']:
        # We can forgive a missing energy, but not a missing molecule.
        head = "{:59} {:10}{:9}\n".format(mol.name(), "(wfn)", "(dertype)")
        head += "{:5d}{:20.10f}".format(mol.natom(), energy)
        for atom in range(mol.natom()):
            head += ("\n" + 4 * "{:20.10f}").format(mol.Z(atom), mol.x(atom), mol.y(atom), mol.z(atom))
        np.savetxt(
            handle, grad, fmt=" " * 20 + "%20.10f%20.10f%20.10f", delimiter='', newline='\n', header=head, comments='')

    elif dtype in ['GRD']:
        head = '{:5}{:20.10f}'.format(mol.natom(), 0)
        for atom in range(mol.natom()):
            head += ("\n" + 4 * "{:20.10f}").format(mol.Z(atom), mol.x(atom), mol.y(atom), mol.z(atom))
        # The need to supply the atomic charge makes writing an actual array inconvenient.
        # For consistency, save a dummy array.
        charges = np.array([mol.Z(i) for i in range(mol.natom())])
        modified_grad = np.insert(grad, 0, charges, axis=1)
        np.savetxt(handle, modified_grad, fmt='%20.10f', delimiter='', newline='\n', header=head, comments='')

    elif dtype in ['minimal']:
        head = '{:5}{:5}'.format(grad.shape[0], grad.size)
        np.savetxt(handle, grad, fmt='%20.10f', delimiter='', newline='\n', header=head, comments='')

    else:
        raise ValidationError('Unknown dtype: {}'.format(dtype))
