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
        ret = np.array(self.wfn.cphf_Hx([matvec])[0]).ravel()

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


class TDRSCFEngine(object):
    def __init__(self, wfn, ptype='rpa', triplet=False, add_tol=1.0e-8):
        self.wfn = wfn
        self.ptype = ptype
        self.triplet = triplet
        self.jk = self.wfn.jk()
        self.Co = wfn.Ca_subset("SO", "OCC")
        self.Cv = wfn.Ca_subset("SO", "VIR")
        self.ei = wfn.epsilon_a_subset("SO", "OCC")
        self.ea = wfn.epsilon_b_subset("SO", "VIR")
        # orbital energy differences
        self.Dia = None
        # symmetry of vector
        self.Gx = None
        self.occpi = self.wfn.nalphapi()
        self.virpi = self.wfn.nmopi() - self.occpi
        self.nsopi = self.wfn.nsopi()
        self.func = wfn.functional()
        self.reset_symmetry(0)
        self.add_tol = add_tol

    def new_vector(self, name=""):
        """Obtain a blank matrix object with the correct symmetry"""
        return core.Matrix(name, self.occpi, self.virpi, self.Gx)

    def build_Dia(self):
        """Rebuild the (Eps_a - Eps_i) matrix for the currently set Guess vector symmetry"""
        self.Dia = self.new_vector("ei-ea")
        for h_i in range(self.wfn.nirrep()):
            h_a = self.Gx ^ h_i
            for i in range(self.occpi[h_i]):
                for a in range(self.virpi[h_a]):
                    self.Dia.nph[h_i][i, a] = self.ea.nph[h_a][a] - self.ei.nph[h_i][i]

    def reset_symmetry(self, symmetry):
        """Reset pre-allocated storage will matrices with the correct symemtry"""
        self.Gx = symmetry
        self.build_Dia()

    def combine_a_plus_b(self, twoel_parts, Fx):
        """Combine Ax, Jx, Kx (if hybrid dft) to make (A+B)x products"""
        nvec = len(Fx)
        Px = []
        for i in range(nvec):
            Px.append(Fx[i].clone())
            if self.func.is_x_hybrid():
                Jxi = twoel_parts[2 * i]
                Kxi = twoel_parts[2 * i + 1]
                Kxi_t = Kxi.transpose()
                Kxi_t.scale(-1.0)
                Kxi_t.axpy(-1.0, Kxi)
                if self.triplet:
                    # (A+B)x = Fx + Co^T(-Kx-Kx^T)Cv
                    Px[i].add(core.Matrix.triplet(self.Co, Kxi_t, self.Cv, True, False, False))
                else:
                    # (A+B)x = Fx + Co^T(4.0J -Kx-Kx^T)Cv
                    Kxi_t.axpy(4.0, Jxi)
                    Px[i].add(core.Matrix.triplet(self.Co, Kxi_t, self.Cv, True, False, False))
            else:
                Jxi = twoel_parts[i]
                if self.triplet:
                    #(A+B)x = Fx
                    pass
                else:
                    #(A+B)x = Fx + Co^T(4.0J)Cv
                    Jxi.scale(4.0)
                    Px[i].add(core.Matrix.triplet(self.Co, Jxi, self.Cv, True, False, False))
        return Px

    def combine_a_minus_b(self, twoel_parts, Fx):
        """Combine Fx, Jx, Kx(if hybrid dft) to make (A-B)x products"""
        nvec = len(Fx)
        Mx = []
        for i in range(nvec):
            Mx.append(Fx[i].clone())
            if self.func.is_x_hybrid():
                Jxi = twoel_parts[2 * i]
                Kxi = twoel_parts[2 * i + 1]
                Kxi_t = Kxi.transpose()
                # (A-B)x = Fx + Co^T(Kx^T-Kx)Cv
                Kxi_t.axpy(-1.0, Kxi)
                Mx[i].add(core.Matrix.triplet(self.Co, Kxi_t, self.Cv, True, False, False))
            # else: (A-B) = Fx
        return Mx

    def combine_a(self, twoel_parts, Fx):
        """Combine Fx, Jx, Kx(if hybrid DFT) to make Ax products"""
        nvec = len(Fx)
        Ax = []
        for i in range(nvec):
            Ax.append(Fx[i].clone())
            if self.func.is_x_hybrid():
                Jxi = twoel_parts[2 * i]
                Kxi = twoel_parts[2 * i + 1]
                if self.triplet:
                    # Ax = Fx - Co^T(K)Cv
                    Ax[i].subtract(core.Matrix.triplet(self.Co, Kxi, self.Cv, True, False, False))
                else:
                    # Ax = Fx + Co^T(2.0Jx - K)Cv
                    Jxi.scale(2.0)
                    Jxi.axpy(-1.0, Kxi)
                    Ax[i].add(core.Matrix.triplet(self.Co, Jxi, self.Cv, True, False, False))
            else:
                Jxi = twoel_parts[i]
                if self.triplet:
                    #Ax = Fx
                    pass
                else:
                    #Ax = Fx + (2.0Jx)
                    Jxi.scale(2.0)
                    Ax[i].add(core.Matrix.triplet(self.Co, Jxi, self.Cv, True, False, False))
        return Ax

    def compute_products(self, X):
        """Given a set of vectors X Compute products
        if ptype == rpa:
           Returns pair (A+B)X, (A-B)X
        if ptype == tda:
           Returns AX
        """
        Fx = self.wfn.onel_Hx(X)
        twoel_parts = self.wfn.twoel_Hx(X, False, "SO")
        if self.ptype == 'rpa':
            Px = self.combine_a_plus_b(twoel_parts, Fx)
            Mx = self.combine_a_minus_b(twoel_parts, Fx)
            # onel/twoel returns - of what we want
            for Pxi in Px:
                Pxi.scale(-1.0)
            for Mxi in Mx:
                Mxi.scale(-1.0)
            return Px, Mx
        else:
            # onel/twoel returns - of what we want
            Ax = self.combine_a(twoel_parts, Fx)
            for Axi in Ax:
                Axi.scale(-1.0)
            return Ax

    def precondition_residual(self, rvec, shift):
        denom = self.new_vector("preconditioner")
        denom.set(shift)
        denom.axpy(-1.0, self.Dia)
        rvec.apply_denominator(denom)
        return rvec

    def approximate_eigvector(self, ax, ss_vector):
        """Approximate the eigenvector of the real system:

        Parameters:
        -----------
        ax: list, of :py:class:`~psi4.core.Matrix(x_dim, x_sym)` length == ss_dim
          The matrix times guess vector products.
        ss_vector: :py:class:`numpy.ndarray` {ss_dim}
          Subspace eigenvector

        Returns:
        --------
        :py:class:`psi4.core.Matrix`
           The approximate eigenvector: $V_i \approx \sum_j A_{ij} x_j \tilde{V}_{ji}$
        """
        ret = self.new_vector("eigenvector")
        for j in range(len(ax)):
            ret.axpy(ss_vector[j], ax[j])
        return ret

    def error(self, x, ss_vector, ss_value):
        """Compute the error piece of the residual expression:
        $\sum_{j} X_j \tilde{V}_{ji}\omega_i$
        """
        ret = self.new_vector("error")
        c = ss_vector * ss_value
        for j in range(len(x)):
            ret.axpy(c[j], x[j])
        return ret

    def vector_norm(self, v):
        return np.sqrt(v.vector_dot(v))

    def orthogonalize_subspace(self, old_X, new_X):
        X = old_X
        while len(new_X) != 0:
            new = new_X.pop()
            for o in X:
                dot = new.vector_dot(o)
                new.axpy(-dot, o)
            norm = np.sqrt(new.vector_dot(new))
            if norm >= self.add_tol:
                new.scale(1.0 / norm)
                X.append(new)
        return X, len(X)

    def residual(self, V, x, ss_vector, ss_value):
        ret = self.error(x, ss_vector, -ss_value)
        ret.axpy(1.0, V)
        return ret

    def subspace_project(self, X, Ax):
        A = np.zeros((len(X), len(X)))
        for i, xi in enumerate(X):
            for j, Axj in enumerate(Ax):
                A[i, j] = xi.vector_dot(Axj)
        return A


class TDUSCFEngine(object):
    def __init__(self, wfn, ptype='rpa', add_tol=1.0e-8):
        self.wfn = wfn
        self.ptype = ptype
        self.jk = self.wfn.jk()
        self.Co = [wfn.Ca_subset("SO", "OCC"), wfn.Cb_subset("SO", "OCC")]
        self.Cv = [wfn.Ca_subset("SO", "VIR"), wfn.Cb_subset("SO", "VIR")]
        self.ei = [wfn.epsilon_a_subset("SO", "OCC"), wfn.epsilon_b_subset("SO", "OCC")]
        self.ea = [wfn.epsilon_a_subset("SO", "VIR"), wfn.epsilon_b_subset("SO", "VIR")]
        # orbital energy differences
        self.Dia_a = None
        # symmetry of vector
        self.Gx = None
        self.occpi = [self.wfn.nalphapi(), self.wfn.nbetapi()]
        self.virpi = [self.wfn.nmopi() - self.occpi[0], self.wfn.nmopi() - self.occpi[1]]
        self.nsopi = self.wfn.nsopi()
        self.func = wfn.functional()
        self.reset_symmetry(0)
        self.add_tol = add_tol

    def new_vector(self, name=""):
        return [
            core.Matrix(name + 'a', self.occpi[0], self.virpi[0], self.Gx),
            core.Matrix(name + 'b', self.occpi[1], self.virpi[1], self.Gx)
        ]

    def build_Dia(self):
        """Rebuild the (Eps_a - Eps_i) matrix for the currently set Guess vector symmetry"""
        self.Dia = self.new_vector("ei-ea ")
        for h_i in range(self.wfn.nirrep()):
            h_a = self.Gx ^ h_i
            for i in range(self.occpi[0][h_i]):
                for a in range(self.virpi[0][h_a]):
                    self.Dia[0].nph[h_i][i, a] = self.ea[0].nph[h_a][a] - self.ei[0].nph[h_i][i]
            for i in range(self.occpi[1][h_i]):
                for a in range(self.virpi[1][h_a]):
                    self.Dia[1].nph[h_i][i, a] = self.ea[1].nph[h_a][a] - self.ei[1].nph[h_i][i]

    def reset_symmetry(self, symmetry):
        self.Gx = symmetry
        self.build_Dia()

    def combine_a_plus_b(self, twoel_parts, Fx):
        """Combine Fx, Jx, Kx (if hybrid dft) to make (A+B)x products"""
        nvec = len(Fx) // 2
        Px = []
        for i in range(nvec):
            Px_a = Fx[2 * i].clone()
            Px_b = Fx[2 * i + 1].clone()
            if self.func.is_x_hybrid():
                Jxi_a = twoel_parts[4 * i]
                Jxi_b = twoel_parts[4 * i + 1]
                Kxi_a = twoel_parts[4 * i + 2]
                Kxi_b = twoel_parts[4 * i + 3]

                # (A+B)x[s] = Fx[s] + Co[s]^T(2.0Jx[s]+2.0Jx[s']-Kx[s]-Kx[s]^T)Cv[s]
                Kxi_a_t = Kxi_a.transpose()
                Kxi_a_t.scale(-1.0)
                Kxi_a_t.axpy(-1.0, Kxi_a)
                Kxi_a_t.axpy(2.0, Jxi_a)
                Kxi_a_t.axpy(2.0, Jxi_b)

                Kxi_b_t = Kxi_b.transpose()
                Kxi_b_t.scale(-1.0)
                Kxi_b_t.axpy(-1.0, Kxi_b)
                Kxi_b_t.axpy(2.0, Jxi_a)
                Kxi_b_t.axpy(2.0, Jxi_b)

                Px_a.add(core.Matrix.triplet(self.Co[0], Kxi_a_t, self.Cv[0], True, False, False))
                Px_b.add(core.Matrix.triplet(self.Co[1], Kxi_b_t, self.Cv[1], True, False, False))

            else:
                #(A+B)x[s] = Fx[s] + Co[s]^T(Jx[s] + Jx[s'])Cv[s]
                Jxi_a = twoel_parts[2 * i]
                Jxi_b = twoel_parts[2 * i + 1]
                Jxi_a.scale(2.0)
                Jxi_a.axpy(2.0, Jxi_b)
                Px_a.add(core.Matrix.triplet(self.Co[0], Jxi_a, self.Cv[0], True, False, False))
                Px_b.add(core.Matrix.triplet(self.Co[1], Jxi_a, self.Cv[1], True, False, False))
            Px.append([Px_a, Px_b])
        return Px

    def combine_a_minus_b(self, twoel_parts, Fx):
        """Combine Fx, Jx, Kx (if hybrid dft) to make (A-B)x products"""
        nvec = len(Fx) // 2
        Mx = []
        for i in range(nvec):
            Mx_a = Fx[2 * i].clone()
            Mx_b = Fx[2 * i + 1].clone()
            if self.func.is_x_hybrid():
                # Jxi_a = twoel_parts[4*i]
                # Jxi_b = twoel_parts[4*i+1]
                Kxi_a = twoel_parts[4 * i + 2]
                Kxi_b = twoel_parts[4 * i + 3]

                Kxi_a_t = Kxi_a.transpose()
                Kxi_a_t.axpy(-1.0, Kxi_a)

                Kxi_b_t = Kxi_b.transpose()
                Kxi_b_t.axpy(-1.0, Kxi_b)

                Mx_a.add(core.Matrix.triplet(self.Co[0], Kxi_a_t, self.Cv[0], True, False, False))
                Mx_b.add(core.Matrix.triplet(self.Co[1], Kxi_b_t, self.Cv[1], True, False, False))
                # (A-B)x[s] = Fx[s] + Co[s]^T(Kx[s]^T-Kx[s])Cv[s]

            else:
                # (A-B)x[s] = Fx[s] (A-B is diagonal for non-hybrid DFT)
                pass
            Mx.append([Mx_a, Mx_b])
        return Mx

    def combine_a(self, twoel_parts, Fx):
        """Combine Fx, Jx, Kx (if hybrid dft) to make Ax products"""
        nvec = len(Fx) // 2
        Ax = []
        for i in range(nvec):
            Ax_a = Fx[2 * i].clone()
            Ax_b = Fx[2 * i + 1].clone()
            if self.func.is_x_hybrid():
                Jxi_a = twoel_parts[4 * i]
                Jxi_b = twoel_parts[4 * i + 1]
                Kxi_a = twoel_parts[4 * i + 2]
                Kxi_b = twoel_parts[4 * i + 3]

                # Ax[s] = Fx[s] + Co[s]^T(1.0Jx[s]+1.0Jx[s']-Kx[s])Cv[s]
                Jxi_a.axpy(1.0, Jxi_b)
                Jxi_b.copy(Jxi_a)

                Jxi_a.axpy(-1.0, Kxi_a)
                Jxi_b.axpy(-1.0, Kxi_b)

                Ax_a.add(core.Matrix.triplet(self.Co[0], Jxi_a, self.Cv[0], True, False, False))
                Ax_b.add(core.Matrix.triplet(self.Co[1], Jxi_b, self.Cv[1], True, False, False))

            else:
                Jxi_a = twoel_parts[2 * i]
                Jxi_b = twoel_parts[2 * i + 1]
                Jxi_a.axpy(1.0, Jxi_b)
                Jxi_b.copy(Jxi_a)

                # Ax[s] = Fx[s] + Co[s]^T(Jx[s]+Jx[s'])Cv[s]
                Ax_a.add(core.Matrix.triplet(self.Co[0], Jxi_a, self.Cv[0], True, False, False))
                Ax_b.add(core.Matrix.triplet(self.Co[1], Jxi_b, self.Cv[1], True, False, False))

            Ax.append([Ax_a, Ax_b])
        return Ax

    def compute_products(self, X):
        """Compute Products for a list of guess vectors (X).
        if ptype == 'rpa':
           returns pair (A+B)X, (A-B)X products
        if ptype == 'tda':
           returns Ax products.
        each element of return is a pair Product(X_a), Product(X_b)
        """
        X_flat = []
        for Xa, Xb in X:
            X_flat.append(Xa)
            X_flat.append(Xb)
        Fx = self.wfn.onel_Hx(X_flat)
        twoel_parts = self.wfn.twoel_Hx(X_flat, False, "SO")
        if self.ptype == 'rpa':
            Px = self.combine_a_plus_b(twoel_parts, Fx)
            Mx = self.combine_a_minus_b(twoel_parts, Fx)
            # onel/twoel return - of what we want
            for Pxi_a, Pxi_b in Px:
                Pxi_a.scale(-1.0)
                Pxi_b.scale(-1.0)
            for Mxi_a, Mxi_b in Mx:
                Mxi_a.scale(-1.0)
                Mxi_b.scale(-1.0)
            return Px, Mx
        else:
            Ax = self.combine_a(twoel_parts, Fx)
            # onel/twoel return - of what we want
            for Axi_a, Axi_b in Ax:
                Axi_a.scale(-1.0)
                Axi_b.scale(-1.0)
            return Ax

    def precondition_residual(self, rvec, shift):
        rva, rvb = rvec
        denom_a, denom_b = self.new_vector("preconditioner")
        denom_a.set(shift)
        denom_b.set(shift)
        denom_a.axpy(-1.0, self.Dia[0])
        denom_b.axpy(-1.0, self.Dia[1])
        rva.apply_denominator(denom_a)
        rvb.apply_denominator(denom_b)
        return rva, rvb

    def approximate_eigvector(self, ax, ss_vector):
        """Approximate the eigenvector of the real system:

        Parameters:
        -----------
        ax: list, of :py:class:`~psi4.core.Matrix(x_dim, x_sym)` length == 2*ss_dim (a/b pairs)
          The matrix times guess vector products.
        ss_vector: :py:class:`numpy.ndarray` {ss_dim}
          Eigenvector of ss problem

        Returns:
        --------
        pair of :py:class:`psi4.core.Matrix` (a/b pair)
           The approximate eigenvector: $V_i \approx \sum_j A_{ij} x_j \tilde{V}_{ji}$
        """
        ret = self.new_vector("eigenvector")
        for j in range(len(ss_vector)):
            ax_a, ax_b = ax[j]
            ret[0].axpy(ss_vector[j], ax_a)
            ret[1].axpy(ss_vector[j], ax_b)
        return ret

    def error(self, x, ss_vector, ss_value):
        """Compute the error piece of the residual expression:
        $\Delta_{i} = \sum_{j} X_j \tilde{V}_{ji}\omega_i$
        """
        ret = self.new_vector("error")
        c = ss_vector * ss_value
        for j,(xa, xb) in enumerate(x):
            ret[0].axpy(c[j], xa)
            ret[1].axpy(c[j], xb)
        return ret

    def orthogonalize_subspace(self, old_X, new_X):
        X = old_X
        while len(new_X) != 0:
            new_a, new_b = new_X.pop()
            for (oa, ob) in X:
                dot = new_a.vector_dot(oa)
                dot += new_b.vector_dot(ob)
                new_a.axpy(-dot, oa)
                new_b.axpy(-dot, ob)
            norm = np.sqrt(new_a.vector_dot(new_a))
            norm += np.sqrt(new_b.vector_dot(new_b))
            if norm >= self.add_tol:
                factor = 1.0 / norm
                new_a.scale(factor)
                new_b.scale(factor)
                X.append([new_a, new_b])
        return X, len(X)

    def vector_norm(self, v):
        norm = np.sqrt(v[0].vector_dot(v[0]))
        norm += np.sqrt(v[1].vector_dot(v[1]))
        return norm

    def residual(self, V, x, ss_vector, ss_value):
        ret = self.error(x, ss_vector, -ss_value)
        ret[0].axpy(1.0, V[0])
        ret[1].axpy(1.0, V[1])
        return ret

    def subspace_project(self, X, Ax):
        A = np.zeros((len(X), len(X)))
        for i, xi in enumerate(X):
            for j, Axj in enumerate(Ax):
                A[i, j] = xi[0].vector_dot(Axj[0])
                A[i, j] += xi[1].vector_dot(Axj[1])
        return A
