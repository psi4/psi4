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

import time

import numpy as np

from psi4 import core
from psi4.driver.p4util.exceptions import *
from psi4.driver import p4util
from psi4.driver import psifiles as psif

from .sapt_util import print_sapt_var


def _symmetrize(mat):
    tmp = 0.5 * (mat + mat.transpose())
    return tmp


def _compute_fxc(PQrho, half_Saux, halfp_Saux, x_alpha, rho_thresh=1.e-8):
    """
    Computes the gridless (P|fxc|Q) ALDA tensor.
    """

    naux = PQrho.shape[0]

    # Level it out
    PQrho_lvl = core.triplet(half_Saux, PQrho, half_Saux, False, False, False)

    # Rotate into a diagonal basis
    rho = core.Vector("rho eigenvalues", naux)
    U = core.Matrix("rho eigenvectors", naux, naux)
    PQrho_lvl.diagonalize(U, rho, core.DiagonalizeOrder.Ascending)

    # "Gridless DFT"
    mask = rho.np < rho_thresh  # Values too small cause singularities
    rho.np[mask] = rho_thresh

    dft_size = rho.shape[0]

    inp = {"RHO_A": rho}
    out = {"V": core.Vector(dft_size), "V_RHO_A": core.Vector(dft_size), "V_RHO_A_RHO_A": core.Vector(dft_size)}

    func_x = core.LibXCFunctional('XC_LDA_X', True)
    func_x.compute_functional(inp, out, dft_size, 2)
    out["V_RHO_A_RHO_A"].scale(1.0 - x_alpha)

    func_c = core.LibXCFunctional('XC_LDA_C_VWN', True)
    func_c.compute_functional(inp, out, dft_size, 2)

    out["V_RHO_A_RHO_A"].np[mask] = 0

    # Rotate back
    Ul = U.clone()
    Ul.np[:] *= out["V_RHO_A_RHO_A"].np
    tmp = core.doublet(Ul, U, False, True)

    # Undo the leveling
    return core.triplet(halfp_Saux, tmp, halfp_Saux, False, False, False)


def df_fdds_dispersion(primary, auxiliary, cache, is_hybrid, x_alpha, leg_points=10, leg_lambda=0.3, do_print=True):

    rho_thresh = core.get_option("SAPT", "SAPT_FDDS_V2_RHO_CUTOFF")
    if do_print:
        core.print_out("\n  ==> E20 Dispersion (CHF FDDS) <== \n\n")
        core.print_out("   Legendre Points:  % 10d\n" % leg_points)
        core.print_out("   Lambda Shift:     % 10.3f\n" % leg_lambda)
        core.print_out("   Fxc Kernal:       % 10s\n" % "ALDA")
        core.print_out("   (P|Fxc|Q) Thresh: % 8.3e\n" % rho_thresh)

    # Build object
    core.timer_on("Form FDDS object")
    df_matrix_keys = ["Cocc_A", "Cvir_A", "Cocc_B", "Cvir_B"]
    fdds_matrix_cache = {key: cache[key] for key in df_matrix_keys}

    df_vector_keys = ["eps_occ_A", "eps_vir_A", "eps_occ_B", "eps_vir_B"]
    fdds_vector_cache = {key: cache[key] for key in df_vector_keys}

    fdds_obj = core.FDDS_Dispersion(primary, auxiliary, fdds_matrix_cache, fdds_vector_cache, is_hybrid)
    core.timer_off("Form FDDS object")

    # Aux Densities
    core.timer_on("Form xc kernel")
    D = fdds_obj.project_densities([cache["D_A"], cache["D_B"]])

    # Temps
    half_Saux = fdds_obj.aux_overlap().clone()
    half_Saux.power(-0.5, 1.e-12)

    halfp_Saux = fdds_obj.aux_overlap().clone()
    halfp_Saux.power(0.5, 1.e-12)

    # Builds potentials
    W_A = fdds_obj.metric().clone()
    W_A.axpy(1.0, _compute_fxc(D[0], half_Saux, halfp_Saux, x_alpha, rho_thresh=rho_thresh))
    W_A = W_A.to_array()

    W_B = fdds_obj.metric().clone()
    W_B.axpy(1.0, _compute_fxc(D[1], half_Saux, halfp_Saux, x_alpha, rho_thresh=rho_thresh))
    W_B = W_B.to_array()

    # Nuke the densities
    del D
    core.timer_off("Form xc kernel")

    metric = fdds_obj.metric().clone().to_array()
    metric_inv = fdds_obj.metric_inv().clone().to_array()

    # Integrate
    core.timer_on("Time Integration")
    core.print_out("\n   => Time Integration <= \n\n")

    val_pack = ("Omega", "Weight", "Disp20,u", "Disp20", "time [s]")
    core.print_out("% 12s % 12s % 14s % 14s % 10s\n" % val_pack)
    start_time = time.time()

    total_uc = 0
    total_c = 0

    # Read R
    if is_hybrid:
        R_A = fdds_obj.R_A().to_array()
        R_B = fdds_obj.R_B().to_array()
        Rtinv_A = np.linalg.pinv(R_A, rcond=1.e-13).transpose()
        Rtinv_B = np.linalg.pinv(R_B, rcond=1.e-13).transpose()

    for point, weight in zip(*np.polynomial.legendre.leggauss(leg_points)):

        omega = leg_lambda * (1.0 - point) / (1.0 + point)
        lambda_scale = ((2.0 * leg_lambda) / (point + 1.0)**2)

        # Monomer A
        if is_hybrid:
            aux_dict = fdds_obj.form_aux_matrices("A", omega)
            aux_dict = {k: v.to_array() for k, v in aux_dict.items()}
            X_A_uc = aux_dict["amp"].copy()
            X_A = X_A_uc - x_alpha * aux_dict["K2L"]

            # K matrices
            K_A = -x_alpha * aux_dict["K1LD"] - x_alpha * aux_dict["K2LD"] + x_alpha * x_alpha * aux_dict["K21L"]
            KRS_A = K_A.dot(Rtinv_A).dot(metric)
        else:
            X_A = fdds_obj.form_unc_amplitude("A", omega)
            X_A.scale(-1.0)
            X_A = X_A.to_array()
            X_A_uc = X_A.copy()

        # Coupled A
        XSW_A = X_A.dot(metric_inv).dot(W_A)
        if is_hybrid:
            XSW_A += 0.25 * KRS_A

        amplitude = np.linalg.pinv(metric - XSW_A, rcond=1.e-13)
        X_A_coupled = X_A + XSW_A.dot(amplitude).dot(X_A)

        del X_A, XSW_A, amplitude
        if is_hybrid:
            del K_A, KRS_A, aux_dict

        # Monomer B
        if is_hybrid:
            aux_dict = fdds_obj.form_aux_matrices("B", omega)
            aux_dict = {k: v.to_array() for k, v in aux_dict.items()}
            X_B_uc = aux_dict["amp"].copy()
            X_B = X_B_uc - x_alpha * aux_dict["K2L"]

            # K matrices
            K_B = -x_alpha * aux_dict["K1LD"] - x_alpha * aux_dict["K2LD"] + x_alpha * x_alpha * aux_dict["K21L"]
            KRS_B = K_B.dot(Rtinv_B).dot(metric)
        else:
            X_B = fdds_obj.form_unc_amplitude("B", omega)
            X_B.scale(-1.0)
            X_B = X_B.to_array()
            X_B_uc = X_B.copy()

        # Coupled B
        XSW_B = X_B.dot(metric_inv).dot(W_B)
        if is_hybrid:
            XSW_B += 0.25 * KRS_B

        amplitude = np.linalg.pinv(metric - XSW_B, rcond=1.e-13)
        X_B_coupled = X_B + XSW_B.dot(amplitude).dot(X_B)

        del X_B, XSW_B, amplitude
        if is_hybrid:
            del K_B, KRS_B, aux_dict

        # Make sure the results are symmetrized
        X_A_uc = _symmetrize(X_A_uc)
        X_B_uc = _symmetrize(X_B_uc)
        X_A_coupled = _symmetrize(X_A_coupled)
        X_B_coupled = _symmetrize(X_B_coupled)

        # Combine
        tmp_uc = metric_inv.dot(X_A_uc).dot(metric_inv)
        value_uc = np.dot(tmp_uc.flatten(), X_B_uc.flatten())
        del tmp_uc

        tmp_c = metric_inv.dot(X_A_coupled).dot(metric_inv)
        value_c = np.dot(tmp_c.flatten(), X_B_coupled.flatten())

        # Tally
        total_uc += value_uc * weight * lambda_scale
        total_c += value_c * weight * lambda_scale

        if do_print:
            tmp_disp_unc = value_uc * weight * lambda_scale
            tmp_disp = value_c * weight * lambda_scale
            fdds_time = time.time() - start_time

            val_pack = (omega, weight, tmp_disp_unc, tmp_disp, fdds_time)
            core.print_out("% 12.3e % 12.3e % 14.3e % 14.3e %10d\n" % val_pack)

    Disp20_uc = -1.0 / (2.0 * np.pi) * total_uc
    Disp20_c = -1.0 / (2.0 * np.pi) * total_c
    core.timer_off("Time Integration")

    core.print_out("\n")
    core.print_out(print_sapt_var("Disp20,u", Disp20_uc, short=True) + "\n")
    core.print_out(print_sapt_var("Disp20", Disp20_c, short=True) + "\n")

    return {"Disp20,FDDS (unc)": Disp20_uc, "Disp20": Disp20_c}


def df_mp2_fisapt_dispersion(wfn, primary, auxiliary, cache, do_print=True):

    if do_print:
        core.print_out("\n  ==> E20 Dispersion (MP2) <== \n\n")

    # Build object
    df_matrix_keys = ["Cocc_A", "Cvir_A", "Cocc_B", "Cvir_B"]
    df_mfisapt_keys = ["Caocc0A", "Cvir0A", "Caocc0B", "Cvir0B"]
    matrix_cache = {fkey: cache[ckey] for ckey, fkey in zip(df_matrix_keys, df_mfisapt_keys)}

    other_keys = ["S", "D_A", "P_A", "V_A", "J_A", "K_A", "D_B", "P_B", "V_B", "J_B", "K_B", "K_O"]
    for key in other_keys:
        matrix_cache[key] = cache[key]

    # matrix_cache["K_O"] = matrix_cache["K_O"].transpose()

    df_vector_keys = ["eps_occ_A", "eps_vir_A", "eps_occ_B", "eps_vir_B"]
    df_vfisapt_keys = ["eps_aocc0A", "eps_vir0A", "eps_aocc0B", "eps_vir0B"]
    vector_cache = {fkey: cache[ckey] for ckey, fkey in zip(df_vector_keys, df_vfisapt_keys)}

    wfn.set_basisset("DF_BASIS_SAPT", auxiliary)
    fisapt = core.FISAPT(wfn)

    # Compute!
    fisapt.disp(matrix_cache, vector_cache, False)
    scalars = fisapt.scalars()

    core.print_out("\n")
    core.print_out(print_sapt_var("Disp20 (MP2)", scalars["Disp20"], short=True) + "\n")
    core.print_out(print_sapt_var("Exch-Disp20,u", scalars["Exch-Disp20"], short=True) + "\n")

    ret = {}
    ret["Exch-Disp20,u"] = scalars["Exch-Disp20"]
    ret["Disp20,u"] = scalars["Disp20"]

    if core.get_option("SAPT", "DO_DISP_EXCH_SINF"):
        fisapt.sinf_disp(matrix_cache, vector_cache, True)
        scalars = fisapt.scalars()

    return ret


def df_mp2_sapt_dispersion(dimer_wfn, wfn_A, wfn_B, primary_basis, aux_basis, cache, do_print=True):

    if do_print:
        core.print_out("\n  ==> E20 Dispersion (MP2) <== \n\n")

    optstash = p4util.OptionsState(['SAPT', 'SAPT0_E10'], ['SAPT', 'SAPT0_E20IND'], ['SAPT', 'SAPT0_E20DISP'],
                                   ['SAPT', 'SAPT_QUIET'])

    core.set_local_option("SAPT", "SAPT0_E10", False)
    core.set_local_option("SAPT", "SAPT0_E20IND", False)
    core.set_local_option("SAPT", "SAPT0_E20DISP", True)
    core.set_local_option("SAPT", "SAPT_QUIET", True)

    if core.get_option('SCF', 'REFERENCE') == 'RHF':
        core.IO.change_file_namespace(psif.PSIF_SAPT_MONOMERA, 'monomerA', 'dimer')
        core.IO.change_file_namespace(psif.PSIF_SAPT_MONOMERB, 'monomerB', 'dimer')

    core.IO.set_default_namespace('dimer')

    dimer_wfn.set_basisset("DF_BASIS_SAPT", aux_basis)
    dimer_wfn.set_basisset("DF_BASIS_ELST", aux_basis)
    e_sapt = core.sapt(dimer_wfn, wfn_A, wfn_B)

    optstash.restore()

    svars = dimer_wfn.variables()

    core.print_out("\n")
    core.print_out(print_sapt_var("Disp20 (MP2)", svars["E DISP20"], short=True) + "\n")
    core.print_out(print_sapt_var("Exch-Disp20,u", svars["E EXCH-DISP20"], short=True) + "\n")

    ret = {}
    ret["Exch-Disp20,u"] = svars["E EXCH-DISP20"]
    ret["Disp20,u"] = svars["E DISP20"]
    return ret
