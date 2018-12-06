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
    def __init__(self, wfn, ptype='rpa', triplet=False):
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
        self.V_pot = wfn.V_potential()
        self.func = wfn.functional()
        self.Vx = []
        self.Dx = []
        # so that C1 wfns don't have to do this
        self.reset_symmetry(0)

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

    def Fx_product(self, v):
        Fxi = self.new_vector("Fx")
        for h in range(self.wfn.nirrep()):
            Fxi.nph[h][:, :] = v.nph[h] * self.Dia.nph[h]
        return Fxi

    def add_Dx(self, D):
        self.Dx.append(D)
        self.Vx.append(self.new_vector("Vx"))

    def compute_products(self, X):
        if self.ptype == 'rpa':
            Px = []
            Mx = []
            for xi in X:
                Pxi, Mxi = self.compute_one(xi)
                Px.append(Pxi)
                Mx.append(Mxi)
            return Px, Mx
        else:
            Ax = []
            for xi in X:
                Ax.append(self.compute_one(xi))
            return Ax

    def compute_one(self, vector):
        self.jk.C_clear()
        self.Dx.clear()
        self.Vx.clear()
        cl = self.Co
        cr = core.Matrix.doublet(self.Cv, vector, False, True)
        self.jk.C_left_add(cl)
        self.jk.C_right_add(cr)
        self.jk.compute()
        if self.func.needs_xc():
            self.add_Dx(core.Matrix.doublet(cl, cr, False, True))
            self.V_pot.compute_Vx(self.Dx, self.Vx)

        Fx = self.Fx_product(vector)

        if self.ptype == 'rpa':
            if self.triplet:
                P = core.Matrix("Px temp", self.nsopi, self.nsopi, self.Gx)
                M = core.Matrix("Mx temp", self.nsopi, self.nsopi, self.Gx)
            else:
                P = self.jk.J()[0]
                M = core.Matrix("Mx temp", self.nsopi, self.nsopi, self.Gx)
                P.scale(4.0)
                if self.func.needs_xc():
                    P.axpy(4.0, self.Vx[0])

            if self.func.is_x_hybrid():
                alpha = self.func.x_alpha()
                K = self.jk.K()[0]
                Kt = K.transpose()
                P.axpy(-alpha, K)
                P.axpy(-alpha, Kt)
                M.print_out()
                Kt.print_out()
                M.axpy(alpha, Kt)
                M.axpy(-alpha, K)

            if self.func.is_x_lrc():
                beta = self.func.x_beta()
                wK = self.jk.wK()[0]
                wKt = wK.transpose()
                P.axpy(-beta, wK)
                P.axpy(-beta, wKt)
                M.axpy(beta, wKt)
                M.axpy(-beta, wK)

            Pret = core.Matrix.triplet(self.Co, P, self.Cv, True, False, False)
            Pret.axpy(1.0, Fx)
            Mret = core.Matrix.triplet(self.Co, M, self.Cv, True, False, False)
            Mret.axpy(1.0, Fx)
            return Pret, Mret

        elif self.ptype == 'tda':
            if self.triplet:
                A = core.Matrix("Ax temp", self.nsopi, self.nsopi, self.Gx)
            else:
                A = self.jk.J()[0]
                A.scale(2.0)
            if self.func.needs_xc():
                A.axpy(2.0, self.Vx[0])
            if self.func.is_x_hybrid():
                alpha = self.func.x_alpha()
                K = self.jk.K()[0]
                A.axpy(-alpha, K)
            if self.func.is_x_lrc():
                beta = self.func.x_beta()
                wK = self.jk.wK()[0]
                A.axpy(-beta, wK)

            Aret = core.Matrix.triplet(self.Co, A, self.Cv, True, False, False)
            Aret.axpy(1.0, Fx)
            return Aret

    def precondition_residual(self, rvec, shift):
        denom = self.new_vector("preconditioner")
        denom.set(shift)
        denom.axpy(-1.0, self.Dia)
        rvec.apply_denomenator(denom)
        return rvec

    def approximate_eigvectors(self, ax, ss_vector):
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
                dot = new.dot(o)
                new.axpy(-dot, o)
            norm = np.sqrt(new.dot(new))
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
    def __init__(self, wfn, ptype='rpa'):
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
        self.V_pot = wfn.V_potential()
        self.func = wfn.functional()
        self.Vx = []
        self.Dx = []
        self.reset_symmetry(0)

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

    def Fx_product(self, vector):
        Fxi = self.new_vector("Fx")
        for h in range(self.wfn.nirrep()):
            Fxi[0].nph[h][:, :] = vector[0].nph[h] * self.Dia[0].nph[h]
            Fxi[1].nph[h][:, :] = vector[1].nph[h] * self.Dia[1].nph[h]
        return Fxi

    def add_Dx(self, D):
        self.Dx.extend(D)
        self.Vx.extend(self.new_vector("Vx"))

    def compute_products(self, X):
        if self.ptype == 'rpa':
            Px = []
            Mx = []
            for xi in X:
                Pxi, Mxi = self.compute_one(xi)
                Px.append(Pxi)
                Mx.append(Mxi)
            return Px, Mx
        else:
            Ax = []
            for xi in X:
                Ax.append(self.compute_one(xi))
            return Ax

    def compute_one(self, vector):
        vector_a, vector_b = vector
        self.jk.C_clear()
        self.Dx.clear()
        self.Vx.clear()

        cl_a, cl_b = self.Co
        Cv_a, Cv_b = self.Cv
        cr_a = core.Matrix.doublet(Cv_a, vector_a, False, True)
        cr_b = core.Matrix.doublet(Cv_b, vector_b, False, True)

        # push alpha occ onto left
        self.jk.C_left_add(cl_a)
        # push alpha vir x vector alpha^T onto right
        self.jk.C_right_add(cr_a)
        # push beta occ onto left
        self.jk.C_left_add(cl_b)
        # push beta vir x vector beta ^T onto right
        self.jk.C_right_add(cr_b)
        self.jk.compute()
        if self.func.needs_xc():
            self.add_Dx([core.Matrix.doublet(cl_a, cr_a, False, True), core.Matrix.doublet(cl_b, cr_b, False, True)])
            self.V_pot.compute_Vx(self.Dx, self.Vx)

        Fxa, Fxb = self.Fx_product(vector)

        if self.ptype == 'rpa':
            P = self.jk.J()
            P[0].scale(2.0)
            P[0].axpy(2.0, P[1])
            P[1].copy(P[0])
            M = [
                core.Matrix("Mx a", self.nsopi, self.nsopi, self.Gx),
                core.Matrix("Mx b", self.nsopi, self.nsopi, self.Gx)
            ]
            if self.func.needs_xc():
                P[0].axpy(2.0, self.Vx[0])
                P[1].axpy(2.0, self.Vx[1])

            if self.func.is_x_hybrid():
                alpha = self.func.x_alpha()
                K = self.jk.K()
                P[0].axpy(-alpha, K[0])
                P[1].axpy(-alpha, K[1])
                M[0].axpy(-alpha, K[0])
                M[1].axpy(-alpha, K[1])
                Kt = [K[0].transpose(), K[1].transpose()]
                P[0].axpy(-alpha, Kt[0])
                P[1].axpy(-alpha, Kt[1])
                M[0].axpy(alpha, Kt[0])
                M[1].axpy(alpha, Kt[1])

            if self.func.is_x_lrc():
                beta = self.func.x_beta()
                wK = self.jk.wK()
                P[0].axpy(-alpha, wK[0])
                P[1].axpy(-alpha, wK[1])
                M[0].axpy(-alpha, wK[0])
                M[1].axpy(-alpha, wK[1])
                wKt = [wK[0].transpose(), wK[1].transpose()]
                P[0].axpy(-alpha, wKt[0])
                P[1].axpy(-alpha, wKt[1])
                M[0].axpy(alpha, wKt[0])
                M[1].axpy(alpha, wKt[1])

            Pret = [
                core.Matrix.triplet(self.Co[a_or_b], P[a_or_b], self.Cv[a_or_b], True, False, False)
                for a_or_b in (0, 1)
            ]
            Pret[0].axpy(1.0, Fxa)
            Pret[1].axpy(1.0, Fxb)

            Mret = [
                core.Matrix.triplet(self.Co[a_or_b], M[a_or_b], self.Cv[a_or_b], True, False, False)
                for a_or_b in (0, 1)
            ]
            Mret[0].axpy(1.0, Fxa)
            Mret[1].axpy(1.0, Fxb)
            return Pret, Mret

        elif self.ptype == 'tda':
            J = self.jk.J()
            J[0].add(J[1])
            # AO basis w/ JK
            Ax = [
                core.Matrix("Ax a", self.nsopi, self.nsopi, self.Gx),
                core.Matrix("Ax b", self.nsopi, self.nsopi, self.Gx)
            ]
            Ax[0].copy(J[0])
            Ax[1].copy(J[0])
            if self.func.needs_xc():
                Ax[0].axpy(1.0, self.Vx[0])
                Ax[1].axpy(1.0, self.Vx[1])
            if self.func.is_x_hybrid():
                alpha = self.func.x_alpha()
                K = self.jk.K()
                Ax[0].axpy(-alpha, K[0])
                Ax[1].axpy(-alpha, K[1])
            if self.func.is_x_lrc():
                beta = self.func.x_beta()
                wK = self.jk.wK()
                Ax[0].axpy(-beta, wK[0])
                Ax[1].axpy(-beta, wK[1])

            Aret = [
                core.Matrix.triplet(self.Co[a_or_b], Ax[a_or_b], self.Cv[a_or_b], True, False, False)
                for a_or_b in (0, 1)
            ]
            Aret[0].axpy(1.0, Fxa)
            Aret[1].axpy(1.0, Fxb)
            return Aret

    def precondition_residual(self, rvec, shift):
        rva, rvb = rvec
        denom_a, denom_b = self.new_vector("preconditioner")
        denom_a.set(shift)
        denom_b.set(shift)
        denom_a.axpy(-1.0, self.Dia[0])
        denom_b.axpy(-1.0, self.Dia[1])
        rva.apply_denomenator(denom_a)
        rvb.apply_denomenator(denom_b)
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
        for ax_a, ax_b in ax:
            ret[0].axpy(ss_vector[j], ax_a)
            ret[1].axpy(ss_vector[j], ax_b)
        return ret

    def error(self, x, ss_vector, ss_value):
        """Compute the error piece of the residual expression:
        $\Delta_{i} = \sum_{j} X_j \tilde{V}_{ji}\omega_i$
        """
        ret = self.new_vector("error")
        c = ss_vector * ss_value
        for xa, xb in X:
            ret[0].axpy(c[j], xa)
            ret[1].axpy(c[j], xb)
        return ret

    def orthogonalize_subspace(self, old_X, new_X):
        X = old_X
        while len(new_X) != 0:
            new_a, new_b = new_X
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
