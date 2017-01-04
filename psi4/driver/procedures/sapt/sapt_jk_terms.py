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

import psi4
import numpy as np
import time
from psi4 import core

from psi4 import extras
from psi4.driver import p4util
from psi4.driver import qcdb
from psi4.driver.p4util.exceptions import *
from psi4.driver.molutil import *
from psi4.driver.procedures.proc import scf_helper

from .sapt_util import print_sapt_var

def build_sapt_jk_cache(wfn_A, wfn_B, jk, do_print=True):
    """
    Constructs the DCBS cache data required to compute ELST/EXCH/IND
    """

    core.print_out("\n   ==> Preparing SAPT Data Cache <== \n\n")
    cache = {}
    cache["wfn_A"] = wfn_A
    cache["wfn_B"] = wfn_B

    # First grab the orbitals
    cache["Cocc_A"] = wfn_A.Ca_subset("SO", "OCC")
    cache["Cvir_A"] = wfn_A.Ca_subset("SO", "VIR")

    cache["Cocc_B"] = wfn_B.Ca_subset("SO", "OCC")
    cache["Cvir_B"] = wfn_B.Ca_subset("SO", "VIR")

    cache["eps_occ_A"] = wfn_A.epsilon_a_subset("SO", "OCC")
    cache["eps_vir_A"] = wfn_A.epsilon_a_subset("SO", "VIR")

    cache["eps_occ_B"] = wfn_B.epsilon_a_subset("SO", "OCC")
    cache["eps_vir_B"] = wfn_B.epsilon_a_subset("SO", "VIR")

    # Build the densities as HF takes an extra "step"
    cache["D_A"] = core.Matrix.doublet(cache["Cocc_A"], cache["Cocc_A"], False, True)
    cache["D_B"] = core.Matrix.doublet(cache["Cocc_B"], cache["Cocc_B"], False, True)

    # Potential ints
    mints = core.MintsHelper(wfn_A.basisset())
    cache["V_A"] = mints.ao_potential()

    mints = core.MintsHelper(wfn_B.basisset())
    cache["V_B"] = mints.ao_potential()

    # Anything else we might need
    cache["S"] = wfn_A.S().clone()

    # J and K matrices
    jk.C_clear()
    jk.C_left_add(wfn_A.Ca_subset("SO", "OCC"))
    jk.C_left_add(wfn_B.Ca_subset("SO", "OCC"))
    jk.compute()

    # Clone them as the JK object will overwrite.
    cache["J_A"] = jk.J()[0].clone()
    cache["K_A"] = jk.K()[0].clone()

    cache["J_B"] = jk.J()[1].clone()
    cache["K_B"] = jk.K()[1].clone()

    monA_nr = wfn_A.molecule().nuclear_repulsion_energy()
    monB_nr = wfn_B.molecule().nuclear_repulsion_energy()
    dimer_nr = wfn_A.molecule().extract_subsets([1, 2]).nuclear_repulsion_energy()

    cache["nuclear_repulsion_energy"] = dimer_nr - monA_nr - monB_nr

    return cache

def electrostatics(cache, do_print=True):
    """
    Computes the E10 electrostatics from a build_sapt_jk_cache datacache.
    """

    core.print_out("\n   ==> E10 Electostatics <== \n\n")
    # ELST
    Elst10  = 4.0 * cache["D_B"].vector_dot(cache["J_A"])
    Elst10 += 2.0 * cache["D_A"].vector_dot(cache["V_B"])
    Elst10 += 2.0 * cache["D_B"].vector_dot(cache["V_A"])
    Elst10 += cache["nuclear_repulsion_energy"]


    core.print_out(print_sapt_var("Elst10,r ", Elst10, short=True))
    core.print_out("\n");
    return {"Elst10,r": Elst10}

def exchange(cache, jk, do_print=True):
    """
    Computes the E10 exchange (S^2 and S^inf) from a build_sapt_jk_cache datacache.
    """

    core.print_out("\n   ==> E10 Exchange <== \n\n")

    # Build potenitals
    h_A = cache["V_A"].clone()
    h_A.axpy(2.0, cache["J_A"])
    h_A.axpy(-1.0, cache["K_A"])

    h_B = cache["V_B"].clone()
    h_B.axpy(2.0, cache["J_B"])
    h_B.axpy(-1.0, cache["K_B"])

    w_A = cache["V_A"].clone()
    w_A.axpy(2.0, cache["J_A"])

    w_B = cache["V_B"].clone()
    w_B.axpy(2.0, cache["J_B"])

    # Build inverse exchange metric
    nocc_A = cache["Cocc_A"].shape[1]
    nocc_B = cache["Cocc_B"].shape[1]
    SAB = core.Matrix.triplet(cache["Cocc_A"], cache["S"], cache["Cocc_B"], True, False, False)
    num_occ = nocc_A + nocc_B

    Sab = core.Matrix(num_occ, num_occ)
    Sab.np[:nocc_A, nocc_A:] = SAB.np
    Sab.np[nocc_A:, :nocc_A] = SAB.np.T
    Sab.np[np.diag_indices_from(Sab.np)] += 1
    Sab.power(-1.0, 1.e-12)
    Sab.np[np.diag_indices_from(Sab.np)] -= 1.0

    Tmo_AA = core.Matrix.from_array(Sab.np[:nocc_A, :nocc_A])
    Tmo_BB = core.Matrix.from_array(Sab.np[nocc_A:, nocc_A:])
    Tmo_AB = core.Matrix.from_array(Sab.np[:nocc_A, nocc_A:])

    T_A  = np.dot(cache["Cocc_A"], Tmo_AA).dot(cache["Cocc_A"].np.T)
    T_B  = np.dot(cache["Cocc_B"], Tmo_BB).dot(cache["Cocc_B"].np.T)
    T_AB = np.dot(cache["Cocc_A"], Tmo_AB).dot(cache["Cocc_B"].np.T)

    P_A = core.Matrix.doublet(cache["Cvir_A"], cache["Cvir_A"], False, True)
    P_B = core.Matrix.doublet(cache["Cvir_B"], cache["Cvir_B"], False, True)

    # Compute the J and K matrices
    jk.C_clear()

    jk.C_left_add(cache["Cocc_A"])
    jk.C_right_add(core.Matrix.doublet(cache["Cocc_A"], Tmo_AA, False, False))

    jk.C_left_add(cache["Cocc_B"])
    jk.C_right_add(core.Matrix.doublet(cache["Cocc_A"], Tmo_AB, False, False))

    jk.C_left_add(psi4.core.Matrix.from_array( np.dot(P_B, cache["S"]).dot(cache["Cocc_A"])))
    jk.C_right_add(cache["Cocc_A"])

    jk.compute()

    JT_A, JT_AB, Jij = jk.J()
    KT_A, KT_AB, Kij = jk.K()

    # Start S^2
    Exch_s2 = 0.0
    Exch_s2 -= 2.0 * np.vdot(np.dot(cache["D_A"], cache["S"]).dot(cache["D_B"]).dot(cache["S"]).dot(P_A), w_B)
    Exch_s2 -= 2.0 * np.vdot(np.dot(cache["D_B"], cache["S"]).dot(cache["D_A"]).dot(cache["S"]).dot(P_B), w_A)
    Exch_s2 -= 2.0 * np.vdot(np.dot(P_A, cache["S"]).dot(cache["D_B"]), Kij.np.T)
    core.print_out(print_sapt_var("Exch10(S^2) ", Exch_s2, short=True))
    core.print_out("\n");


    # Start Sinf
    Exch10  = 0.0
    Exch10 -= 2.0 * np.vdot(cache["D_A"], cache["K_B"])
    Exch10 += 2.0 * np.vdot(T_A, h_B.np)
    Exch10 += 2.0 * np.vdot(T_B, h_A.np)
    Exch10 += 2.0 * np.vdot(T_AB, h_A.np + h_B.np)
    Exch10 += 4.0 * np.vdot(T_B, JT_AB.np - 0.5 * KT_AB.np)
    Exch10 += 4.0 * np.vdot(T_A, JT_AB.np - 0.5 * KT_AB.np)
    Exch10 += 4.0 * np.vdot(T_B, JT_A.np - 0.5 *  KT_A.np)
    Exch10 += 4.0 * np.vdot(T_AB, JT_AB.np - 0.5 * KT_AB.np.T)
    core.print_out(print_sapt_var("Exch10", Exch10, short=True))
    core.print_out("\n");

    return {"Exch10(S^2)": Exch_s2, "Exch10": Exch10}

def induction(cache, jk, do_print=True, maxiter=20, conv=1.e-8, do_response=True):

    core.print_out("\n   ==> E20 Induction <== \n\n")

    # Build Induction and Exchange-Induction potentials
    S = cache["S"].np

    D_A = cache["D_A"].np
    V_A = cache["V_A"].np
    J_A = cache["J_A"].np
    K_A = cache["K_A"].np

    D_B = cache["D_B"].np
    V_B = cache["V_B"].np
    J_B = cache["J_B"].np
    K_B = cache["K_B"].np

    jk.C_clear()

    C_O_A = psi4.core.Matrix.from_array(D_B.dot(S).dot(cache["Cocc_A"]))
    C_P_A = psi4.core.Matrix.from_array(D_B.dot(S).dot(D_A).dot(S).dot(cache["Cocc_B"]))
    C_P_B = psi4.core.Matrix.from_array(D_A.dot(S).dot(D_B).dot(S).dot(cache["Cocc_A"]))

    jk.C_left_add(C_O_A)
    jk.C_right_add(cache["Cocc_A"])

    jk.C_left_add(C_P_A)
    jk.C_right_add(cache["Cocc_B"])

    jk.C_left_add(C_P_B)
    jk.C_right_add(cache["Cocc_A"])

    jk.compute()

    J_O, J_P_B, J_P_A = jk.J()
    K_O, K_P_B, K_P_A = jk.K()

    K_O = psi4.core.Matrix.from_array(K_O.np.T)
    W_A  = -1.0 * K_B.copy()
    W_A -= 2.0 * np.dot(S, D_B).dot(J_A)
    W_A += 1.0 * K_O.np
    W_A -= 2.0 * J_O.np

    W_A += 1.0 * np.dot(S, D_B).dot(K_A)
    W_A -= 2.0 * np.dot(J_B, D_B).dot(S)
    W_A += 1.0 * np.dot(K_B, D_B).dot(S)

    W_A += 2.0 * np.dot(S, D_B).dot(J_A).dot(D_B).dot(S)
    W_A += 2.0 * np.dot(J_B, D_A).dot(S).dot(D_B).dot(S)
    W_A -= 1.0 * np.dot(K_O, D_B).dot(S)
    W_A += 2.0 * J_P_B.np

    W_A += 2.0 * np.dot(S, D_B).dot(S).dot(D_A).dot(J_B)
    W_A -= 1.0 * np.dot(S, D_B).dot(K_O.np.T)
    W_A -= 1.0 * np.dot(S, D_B).dot(V_A)
    W_A -= 1.0 * np.dot(V_B, D_B).dot(S)
    W_A += 1.0 * np.dot(S, D_B).dot(V_A).dot(D_B).dot(S)
    W_A += 1.0 * np.dot(V_B, D_A).dot(S).dot(D_B).dot(S)
    W_A += 1.0 * np.dot(S, D_B).dot(S).dot(D_A).dot(V_B)
    W_A = np.dot(cache["Cocc_A"].np.T, W_A).dot(cache["Cvir_A"])

    K_O = psi4.core.Matrix.from_array(K_O.np.T)
    W_B  = -1.0 * K_A.copy()
    W_B -= 2.0 * np.dot(S, D_A).dot(J_B)
    W_B += 1.0 * K_O.np
    W_B -= 2.0 * J_O.np
    W_B += 1.0 * np.dot(S, D_A).dot(K_B)
    W_B -= 2.0 * np.dot(J_A, D_A).dot(S)
    W_B += 1.0 * np.dot(K_A, D_A).dot(S)
    W_B += 2.0 * np.dot(S, D_A).dot(J_B).dot(D_A).dot(S)
    W_B += 2.0 * np.dot(J_A, D_B).dot(S).dot(D_A).dot(S)
    W_B -= 1.0 * np.dot(K_O, D_A).dot(S)
    W_B += 2.0 * J_P_A.np
    W_B += 2.0 * np.dot(S, D_A).dot(S).dot(D_B).dot(J_A)
    W_B -= 1.0 * np.dot(S, D_A).dot(K_O.np.T)
    W_B -= 1.0 * np.dot(S, D_A).dot(V_B)
    W_B -= 1.0 * np.dot(V_A, D_A).dot(S)
    W_B += 1.0 * np.dot(S, D_A).dot(V_B).dot(D_A).dot(S)
    W_B += 1.0 * np.dot(V_A, D_B).dot(S).dot(D_A).dot(S)
    W_B += 1.0 * np.dot(S, D_A).dot(S).dot(D_B).dot(V_A)
    W_B = np.dot(cache["Cocc_B"].np.T, W_B).dot(cache["Cvir_B"])

    # Build electrostatic potenital
    w_A = cache["V_A"].clone()
    w_A.axpy(2.0, cache["J_A"])

    w_B = cache["V_B"].clone()
    w_B.axpy(2.0, cache["J_B"])

    w_B_MOA = core.Matrix.triplet(cache["Cocc_A"], w_B, cache["Cvir_A"], True, False, False)
    w_A_MOB = core.Matrix.triplet(cache["Cocc_B"], w_A, cache["Cvir_B"], True, False, False)

    # Do uncoupled
    core.print_out("    => Uncoupled Induction <= \n\n")
    unc_x_B_MOA = w_B_MOA.np / (cache["eps_occ_A"].np.reshape(-1, 1) - cache["eps_vir_A"].np)
    unc_x_A_MOB = w_A_MOB.np / (cache["eps_occ_B"].np.reshape(-1, 1) - cache["eps_vir_B"].np)

    unc_ind_ab = 2.0 * np.vdot(unc_x_B_MOA, w_B_MOA)
    unc_ind_ba = 2.0 * np.vdot(unc_x_A_MOB, w_A_MOB)
    unc_indexch_ba = 2.0 * np.vdot(W_B, unc_x_A_MOB)
    unc_indexch_ab = 2.0 * np.vdot(W_A, unc_x_B_MOA)

    ret = {}
    ret["Ind20,u (A<-B)"] = unc_ind_ab
    ret["Ind20,u (A->B)"] = unc_ind_ba
    ret["Ind20,u"] = unc_ind_ab + unc_ind_ba
    ret["Ind-Exch20,u (A<-B)"] = unc_indexch_ab
    ret["Ind-Exch20,u (A->B)"] = unc_indexch_ba
    ret["Ind-Exch20,u"] = unc_indexch_ba + unc_indexch_ab

    plist = ["Ind20,u (A<-B)", "Ind20,u (A->B)", "Ind20,u", "Ind-Exch20,u (A<-B)",
             "Ind-Exch20,u (A->B)", "Ind-Exch20,u"]

    for name in plist:
        core.print_out(print_sapt_var(name, ret[name], short=True))
        core.print_out("\n");

    # Do coupled
    if do_response:
        core.print_out("\n    => Coupled Induction <= \n\n")
        core.print_out("    => CPHF Monomer A <= \n")
        x_B_MOA = cache["wfn_A"].cphf_solve([w_B_MOA], conv, maxiter, 2)[0]
        core.print_out("    => CPHF Monomer B <= \n")
        x_A_MOB = cache["wfn_B"].cphf_solve([w_A_MOB], conv, maxiter, 2)[0]

        x_B = np.dot(cache["Cocc_A"], x_B_MOA).dot(cache["Cvir_A"].np.T)
        x_A = np.dot(cache["Cocc_B"], x_A_MOB).dot(cache["Cvir_B"].np.T)

        ind_ab = 2.0 * x_B_MOA.vector_dot(w_B_MOA)
        ind_ba = 2.0 * x_A_MOB.vector_dot(w_A_MOB)
        indexch_ba = 2.0 * np.vdot(W_B, x_A_MOB)
        indexch_ab = 2.0 * np.vdot(W_A, x_B_MOA)

        ret["Ind20,r (A<-B)"] = ind_ab
        ret["Ind20,r (A->B)"] = ind_ba
        ret["Ind20,r"] = ind_ab + ind_ba
        ret["Ind-Exch20,r (A<-B)"] = indexch_ab
        ret["Ind-Exch20,r (A->B)"] = indexch_ba
        ret["Ind-Exch20,r"] = indexch_ba + indexch_ab

        for name in plist:
            name = name.replace(",u", ",r")
            core.print_out(print_sapt_var(name, ret[name], short=True))
            core.print_out("\n");

    return ret
