#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2018 The Psi4 Developers.
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

from psi4 import core
from psi4.driver.p4util.exceptions import *


class SCFProducts:
    def __init__(self, wfn):

        if len(tuple(wfn.doccpi())) > 1:
            raise ValueError("Can only handle C1 wavefunctions currently.")

        self.restricted = wfn.same_a_b_dens()
        if self.restricted is False:
            raise ValueError("Can only handle restricted wavefunctions.")

        self.wfn = wfn

        # Grab sizes
        self.nmo = wfn.nmopi().sum()
        self.nalpha = wfn.nalphapi().sum()
        self.navir = self.nmo - self.nalpha
        self.narot = self.nalpha * self.navir

        self.nbeta = wfn.nbetapi().sum()
        self.nbvir = self.nmo - self.nbeta
        self.nbrot = self.nbeta * self.nbvir

        self.alpha_shape = (self.nalpha, self.navir)
        self.beta_shape = (self.nbeta, self.nbvir)

        # Grab orbitals
        self.Ca_occ = wfn.Ca_subset("AO", "OCC")
        self.Ca_vir = wfn.Ca_subset("AO", "VIR")
        self.Cb_occ = wfn.Cb_subset("AO", "OCC")
        self.Cb_vir = wfn.Cb_subset("AO", "VIR")

    def H2_product(self, vector):
        """
        H2 = H2_ai,bj * vector_bj
                          K         K.T
        [(E_a - E_i) + (ab|ij) - (aj|bi)] x_bj
        Note that we can only take a single vector, be good to take many
        """
        vector = core.Matrix.from_array(vector.reshape(self.alpha_shape))
        ret = self.wfn.onel_Hx([vector])[0]

        # Special case for DFT
        if self.wfn.functional().needs_xc():
            Jx, Kx, XCx = self.wfn.twoel_Hx([vector], False, "SO")
            Kx.axpy(-1.0, Kx.transpose())
            Kx.axpy(1.0, XCx)
        else:
            Jx, Kx = self.wfn.twoel_Hx([vector], False, "SO")
            Kx.axpy(-1.0, Kx.transpose())
        twoel_mo = core.Matrix.triplet(self.Ca_occ, Kx, self.Ca_vir, True, False, False)

        # Some signs get mixed around, but this is correct
        ret.axpy(-1.0, twoel_mo)
        ret = np.array(ret).ravel()
        return ret

    def H1_product(self, vector):
        """
        H1 (Orbital Hessian)
        """
        matvec = core.Matrix.from_array(vector.reshape(self.alpha_shape))
        ret = np.array(wfn.cphf_Hx([matvec])[0]).ravel()

        return ret

    def H2H1_product(self, vector, omega=None):
        """
        H2 * H1 * x + omega ** 2
        """
        matvec = core.Matrix.from_array(vector.reshape(self.alpha_shape))
        tmp = np.array(self.wfn.cphf_Hx([matvec])[0]).ravel()

        ret = self.H2_product(tmp)
        if omega is not None:
            ret += vector * (omega**2)

        return np.array(ret).ravel()