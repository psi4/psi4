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

import numpy as np
import time

from psi4 import core
from psi4.driver.p4util.exceptions import *

from .sapt_util import print_sapt_var


def _compute_fxc(PQrho, half_Saux, halfp_Saux):
    """
    Computes the gridless (P|fxc|Q) ALDA tensor.
    """

    naux = PQrho.shape[0]

    # Level it out
    PQrho_lvl = core.Matrix.triplet(half_Saux, PQrho, half_Saux, False, False, False)
    # PQrho_lvl = np.dot(half_Saux, PQrho).dot(half_Saux)

    # Rotate into a diagonal basis
    rho = core.Vector("rho", naux)
    U = core.Matrix("rho eigenvectors", naux, naux)
    PQrho_lvl.diagonalize(U, rho, core.DiagonalizeOrder.Descending)

    # "Gridless DFT"
    mask = rho.np < 1.e-6 # Values too small cause singularities
    rho.np[mask] = 0

    dft_size = rho.shape[0]

    inp = {"RHO_A": rho}
    out = {"V": core.Vector(dft_size),
           "V_RHO_A" : core.Vector(dft_size),
           "V_RHO_A_RHO_A": core.Vector(dft_size)}

    func_x = core.LibXCFunctional('XC_LDA_X', True)
    func_x.compute_functional(inp, out, dft_size, 2)

    func_c = core.LibXCFunctional('XC_LDA_C_VWN', True)
    func_c.compute_functional(inp, out, dft_size, 2)

    out["V_RHO_A_RHO_A"].np[mask] = 0

    # Original basis
    Ul = U.clone()
    Ul.np[:] *= out["V_RHO_A_RHO_A"].np
    tmp = core.Matrix.doublet(Ul, U, False, True)

    # Undo the leveling
    return core.Matrix.triplet(halfp_Saux, tmp, halfp_Saux, False, False, False)

def df_fdds_dispersion(primary, auxiliary, cache, leg_points=10, leg_lambda=0.3, do_print=True):

    if do_print:
        core.print_out("\n  ==> E20 Dispersion (CHF FDDS) <== \n\n")
        core.print_out("   Legendre Points: % 8d\n" % leg_points)
        core.print_out("   Lambda Shift:    % 8.3f\n" % leg_lambda)
        core.print_out("   Fxc Kernal:      % 8s\n\n" % "ALDA")

    # Build object
    df_matrix_keys = ["Cocc_A", "Cvir_A", "Cocc_B", "Cvir_B"]
    fdds_matrix_cache = {key : cache[key] for key in df_matrix_keys}

    df_vector_keys = ["eps_occ_A", "eps_vir_A", "eps_occ_B", "eps_vir_B"]
    fdds_vector_cache = {key : cache[key] for key in df_vector_keys}

    fdds_obj = core.FDDS_Dispersion(primary, auxiliary, fdds_matrix_cache, fdds_vector_cache)

    # Aux Densities
    D = fdds_obj.project_densities([cache["D_A"], cache["D_B"]])

    # Temps
    half_Saux = fdds_obj.aux_overlap().clone()
    half_Saux.power(-0.5, 1.e-12)

    halfp_Saux = fdds_obj.aux_overlap().clone()
    halfp_Saux.power(0.5, 1.e-12)

    # Builds potentials
    W_A = fdds_obj.metric().clone()
    W_A.axpy(1.0, _compute_fxc(D[0], half_Saux, halfp_Saux))

    W_B = fdds_obj.metric().clone()
    W_B.axpy(1.0, _compute_fxc(D[1], half_Saux, halfp_Saux))

    # Nuke the densities
    del D

    metric = fdds_obj.metric()
    metric_inv = fdds_obj.metric_inv()

    total_uc = 0
    total_c = 0

    core.print_out("\n   => Time Integration <= \n\n")

    val_pack = ("Omega", "Weight", "Disp20,u", "Disp20", "time [s]")
    core.print_out("% 12s % 12s % 14s % 14s % 10s\n" % val_pack)
    for point, weight in zip(*np.polynomial.legendre.leggauss(leg_points)):
        start_time = time.time()

        omega = leg_lambda * (1.0 - point) / (1.0 + point)
        lambda_scale = ( (2 * leg_lambda) / (point + 1) ** 2)

        # Monomer A
        X_A = fdds_obj.form_unc_amplitude("A", omega)

        X_A_coupled = X_A.clone()
        XSW_A = core.Matrix.triplet(X_A, metric_inv, W_A, False, False, False)

        amplitude = metric.clone()
        amplitude.axpy(1.0, XSW_A)
        amplitude.power(-1.0, 1.e-12)
        X_A_coupled.axpy(-1.0, core.Matrix.triplet(XSW_A, amplitude, X_A, False, False, False))

        del XSW_A, amplitude

        # Monomer B
        X_B = fdds_obj.form_unc_amplitude("B", omega)

        X_B_coupled = X_B.clone()
        XSW_B = core.Matrix.triplet(X_B, metric_inv, W_B, False, False, False)

        amplitude = metric.clone()
        amplitude.axpy(1.0, XSW_B)
        amplitude.power(-1.0, 1.e-12)
        X_B_coupled.axpy(-1.0, core.Matrix.triplet(XSW_B, amplitude, X_B, False, False, False))
        del XSW_B, amplitude

        tmp_uc = core.Matrix.triplet(metric_inv, X_A, metric_inv, False, False, False)
        value_uc = tmp_uc.vector_dot(X_B)
        del tmp_uc

        tmp_c = core.Matrix.triplet(metric_inv, X_A_coupled, metric_inv, False, False, False)
        value_c = tmp_c.vector_dot(X_B_coupled)
        del tmp_c

        total_uc += value_uc * weight * lambda_scale
        total_c += value_c * weight * lambda_scale

        if do_print:
            tmp_disp_unc = value_uc * weight * lambda_scale
            tmp_disp = value_c * weight * lambda_scale
            fdds_time = time.time() - start_time

            val_pack = (omega, weight, tmp_disp_unc, tmp_disp, fdds_time)
            core.print_out("% 12.3e % 12.3e % 14.3e % 14.3e %10d\n" % val_pack)

    Disp20_uc = -1.0 / (2.0 * np.pi) * total_uc
    Disp20_c = -1.0 / (2.0  * np.pi) * total_c

    core.print_out("\n")
    core.print_out(print_sapt_var("Disp20,u", Disp20_uc, short=True) + "\n")
    core.print_out(print_sapt_var("Disp20", Disp20_c, short=True) + "\n")

    return {"Disp20,FDDS (unc)" : Disp20_uc, "Disp20,FDDS" : Disp20_c}



