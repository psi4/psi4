#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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

from __future__ import division
from __future__ import print_function

import numpy as np


def random_rotation_matrix(deflection=1.0, randnums=None):
    """Generates a random 3D rotation matrix.

    Parameters
    ----------
    deflection : float, optional
        Magnitude of the rotation. For 0, no rotation; for 1, competely random
        rotation. Small deflection => small perturbation.
    randnums : array, optional
        3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.

    Returns
    -------
    3 x 3 ndarray
        Rotation matrix on random distribution.

    Sources
    -------
    improved: http://demonstrations.wolfram.com/sourcecode.html?demoname=SamplingAUniformlyRandomRotation&demodisplayname=Sampling%20a%20Uniformly%20Random%20Rotation
    py code: from http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
    c code: from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
    orig: James Arvo, Graphics Gems III (IBM Version) 1992, Pages 117-120, https://doi.org/10.1016/B978-0-08-050755-2.50034-8

    """
    if randnums is None:
        randnums = np.random.uniform(size=(3, ))

    theta, phi, z = randnums

    # rotation about the pole (Z)
    #   from Wolfram, improved by subtr half so rotation unbiased
    theta = (theta - 1 / 2) * deflection * 2 * np.pi
    # direction of pole deflection
    phi = phi * 2 * np.pi
    # magnitude of pole deflection
    z = z * 2 * deflection

    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
    r = np.sqrt(z)
    V = (np.sin(phi) * r, np.cos(phi) * r, np.sqrt(2.0 - z))

    st = np.sin(theta)
    ct = np.cos(theta)

    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
    R_z_pi = np.diag([-1., -1., 1])

    # Construct the rotation matrix  (V Transpose(V) - I) R. * R_z(pi)
    #   which is equivalent to V.S - R.
    # From Wolfram, Arno's code is missing the multiplication by R_z(pi),
    #   which is unnoticable for random rotations but causes problems
    #   for random perturbations.
    M = (np.outer(V, V) - np.eye(3)).dot(R).dot(R_z_pi)

    return M
