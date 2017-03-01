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


import numpy as np
from psi4 import core

class DIIS_helper(object):
    """
    Simple DIIS class.
    """

    def __init__(self, max_vec=6):
        self.error = []
        self.vector = []
        self.max_vec = max_vec

    def add(self, matrix, error):
        #if len(self.error) > 1:
        #    if self.error[-1].shape[0] != error.size:
        #        raise Exception("Error vector size does not match previous vector.")
        #    if self.vector[-1].shape != matrix.shape:
        #        raise Exception("Vector shape does not match previous vector.")

        self.error.append(error.clone())
        self.vector.append(matrix.clone())

    def extrapolate(self):
        # Limit size of DIIS vector
        diis_count = len(self.vector)

        if diis_count == 0:
            raise Exception("DIIS: No previous vectors.")
        if diis_count == 1:
            return self.vector[0]

        if diis_count > self.max_vec:
            # Remove oldest vector
            del self.vector[0]
            del self.error[0]
            diis_count -= 1

        # Build error matrix B
        B = np.empty((diis_count + 1, diis_count + 1))
        B[-1, :] = 1
        B[:, -1] = 1
        B[-1, -1] = 0
        for num1, e1 in enumerate(self.error):
            B[num1, num1] = e1.vector_dot(e1)
            for num2, e2 in enumerate(self.error):
                if num2 >= num1: continue
                val = e1.vector_dot(e2)
                B[num1, num2] = B[num2, num1] = val

        # Build residual vector
        resid = np.zeros(diis_count + 1)
        resid[-1] = 1

        # Solve pulay equations

        # Yea, yea this is unstable make it stable
        iszero = np.any(np.diag(B)[:-1] <= 0.0)
        if iszero:
            S = np.ones((diis_count + 1))
        else:
            S = np.ones((diis_count + 1))
            S[:-1] = np.diag(B)[:-1]
            S = S ** -0.5
            S[-1] = 1

        # Then we gotta do a custom inverse
        B *= S[:, None] * S

        invB = core.Matrix.from_array(B)
        invB.power(-1.0, 1.e-12)

        ci = np.dot(invB, resid) * S

        # combination of previous fock matrices
        V = core.Matrix("DIIS result", self.vector[0].rowdim(), self.vector[1].coldim())
        for num, c in enumerate(ci[:-1]):
            V.axpy(c, self.vector[num])

        return V
