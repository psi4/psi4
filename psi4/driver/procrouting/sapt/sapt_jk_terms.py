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

import time

import numpy as np

from psi4 import core

from ...p4util import solvers
from ...p4util.exceptions import *
from .sapt_util import print_sapt_var
from pprint import pprint as pp


def build_sapt_jk_cache(
        wfn_dimer: core.Wavefunction,
        wfn_A: core.Wavefunction,
        wfn_B: core.Wavefunction,
        jk: core.JK,
        do_print=True,
):
    """
    Constructs the DCBS cache data required to compute ELST/EXCH/IND
    """

    core.print_out("\n  ==> Preparing SAPT Data Cache <== \n\n")
    jk.print_header()

    cache = {}
    cache["wfn_A"] = wfn_A
    cache["wfn_B"] = wfn_B

    # First grab the orbitals
    cache["Cocc_A"] = wfn_A.Ca_subset("AO", "OCC")
    cache["Cvir_A"] = wfn_A.Ca_subset("AO", "VIR")

    cache["Cocc_B"] = wfn_B.Ca_subset("AO", "OCC")
    cache["Cvir_B"] = wfn_B.Ca_subset("AO", "VIR")

    cache["eps_occ_A"] = wfn_A.epsilon_a_subset("AO", "OCC")
    cache["eps_vir_A"] = wfn_A.epsilon_a_subset("AO", "VIR")

    cache["eps_occ_B"] = wfn_B.epsilon_a_subset("AO", "OCC")
    cache["eps_vir_B"] = wfn_B.epsilon_a_subset("AO", "VIR")

    # Build the densities as HF takes an extra "step"
    cache["D_A"] = core.doublet(cache["Cocc_A"], cache["Cocc_A"], False, True)
    cache["D_B"] = core.doublet(cache["Cocc_B"], cache["Cocc_B"], False, True)

    cache["P_A"] = core.doublet(cache["Cvir_A"], cache["Cvir_A"], False, True)
    cache["P_B"] = core.doublet(cache["Cvir_B"], cache["Cvir_B"], False, True)

    # External Potential
    # nA = wfn_A.molecule().natom()
    nbf = wfn_dimer.basisset().nbf()
    # print(nA, nB)
    # TODO: ensure Enucs is the right dimensionality for different sized
    # systems..., water dimer is 3x3
    ext_matrices = {
        "Enucs": core.Matrix.from_array(np.zeros((3, 3)), name="Enucs"),
        "VA": core.Matrix.from_array(np.zeros((nbf, nbf)), name="VA"),
        "VB": core.Matrix.from_array(np.zeros((nbf, nbf)), name="VB"),
        "VA_extern": core.Matrix.from_array(np.zeros((nbf, nbf)), name="VA_extern"),
        "VB_extern": core.Matrix.from_array(np.zeros((nbf, nbf)), name="VB_extern"),
        # "VA_extern": wfn_dimer.potential_variable("A").computePotentialMatrix(wfn_dimer.basisset()),
        # "VB_extern": wfn_dimer.potential_variable("B").computePotentialMatrix(wfn_dimer.basisset()),
        # "VA": wfn_dimer.potential_variable("A").computePotentialMatrix(wfn_dimer.basisset()),
        # "VB": wfn_dimer.potential_variable("B").computePotentialMatrix(wfn_dimer.basisset()),
        # "VA_extern": wfn_dimer.potential_variable("A").computePotentialMatrix(wfn_dimer.basisset()),
        # "VB_extern": wfn_dimer.potential_variable("B").computePotentialMatrix(wfn_dimer.basisset()),
    }
    print("pre sapt_nuclear_external_potential_python")
    Etot = core.sapt_nuclear_external_potential_python(
        wfn_dimer,
        ext_matrices,
        core.get_options()
    )
    print(f"Etot: {Etot}")
    print(f"Enucs:\n{ext_matrices['Enucs'].np}")
    for key, value in ext_matrices.items():
        print(f"{key}: {value.name}")
        cache[key] = value

    # Potential ints
    mints = core.MintsHelper(wfn_A.basisset())
    cache["V_A"] = mints.ao_potential()
    # check if VA is all zeros
    if np.allclose(ext_matrices["VA"].np, 0):
        mints = core.MintsHelper(wfn_A.basisset())
        cache["V_A"] = mints.ao_potential()
    else:
        mints = core.MintsHelper(wfn_A.basisset())
        cache["V_A"] = mints.ao_potential()
        # cache["V_A"] = ext_matrices['VA']
    if np.allclose(ext_matrices["VB"].np, 0):
        mints = core.MintsHelper(wfn_B.basisset())
        cache["V_B"] = mints.ao_potential()
    else:
        mints = core.MintsHelper(wfn_B.basisset())
        cache["V_B"] = mints.ao_potential()
        # cache["V_B"] = ext_matrices['VB']
    # if "VA" not in ext_matrices.keys():
    #     mints = core.MintsHelper(wfn_A.basisset())
    #     cache["V_A"] = mints.ao_potential()
    # else:
    #     cache["V_A"] = ext_matrices['VA']
    # TODO: add external potential
    # cache["V_A"].axpy(1.0, wfn_A.Va())

    # mints = core.MintsHelper(wfn_B.basisset())
    # cache["V_B"] = mints.ao_potential()
    # if "VB" not in ext_matrices.keys():
    #     mints = core.MintsHelper(wfn_B.basisset())
    #     cache["V_B"] = mints.ao_potential()
    # else:
    #     cache["V_B"] = ext_matrices['VB']
    # mints = core.MintsHelper(wfn_B.basisset())
    # cache["V_B"] = mints.ao_potential()
    # cache["V_B"].add(ext_matrices['VB'])
    # cache["V_B"].axpy(1.0, wfn_B.Va())

    # Anything else we might need
    cache["S"] = wfn_A.S().clone()

    # J and K matrices
    jk.C_clear()

    # Normal J/K for Monomer A
    jk.C_left_add(wfn_A.Ca_subset("SO", "OCC"))
    jk.C_right_add(wfn_A.Ca_subset("SO", "OCC"))

    # Normal J/K for Monomer B
    jk.C_left_add(wfn_B.Ca_subset("SO", "OCC"))
    jk.C_right_add(wfn_B.Ca_subset("SO", "OCC"))

    # K_O J/K
    C_O_A = core.triplet(cache["D_B"], cache["S"], cache["Cocc_A"], False, False, False)
    jk.C_left_add(C_O_A)
    jk.C_right_add(cache["Cocc_A"])

    jk.compute()

    # Clone them as the JK object will overwrite.
    cache["J_A"] = jk.J()[0].clone()
    cache["K_A"] = jk.K()[0].clone()

    cache["J_B"] = jk.J()[1].clone()
    cache["K_B"] = jk.K()[1].clone()

    cache["J_O"] = jk.J()[2].clone()
    cache["K_O"] = jk.K()[2].clone()
    cache["K_O"].transpose_this()

    monA_nr = wfn_A.molecule().nuclear_repulsion_energy()
    monB_nr = wfn_B.molecule().nuclear_repulsion_energy()
    dimer_nr = wfn_A.molecule().extract_subsets([1, 2]).nuclear_repulsion_energy()

    cache["nuclear_repulsion_energy"] = dimer_nr - monA_nr - monB_nr
    return cache


def electrostatics(cache, do_print=True):
    """
    Computes the E10 electrostatics from a build_sapt_jk_cache datacache.
    """
    print("CACHE")
    pp(cache)

    if do_print:
        core.print_out("\n  ==> E10 Electostatics <== \n\n")

    # ELST
    Elst10 = 4.0 * cache["D_B"].vector_dot(cache["J_A"])
    Elst10 += 2.0 * cache["D_A"].vector_dot(cache["V_B"])
    Elst10 += 2.0 * cache["D_B"].vector_dot(cache["V_A"])
    Elst10 += cache["nuclear_repulsion_energy"]

    if do_print:
        core.print_out(print_sapt_var("Elst10,r ", Elst10, short=True))
        core.print_out("\n")
    extern_extern_ie = 0
    if cache.get('extern_extern_IE'):
        extern_extern_ie = cache['extern_extern_IE'].get(0, 1) * 2.0
        core.print_out(f"    Extern-Extern               {extern_extern_ie*1000:16.8f} [mEh]\n")

    return {"Elst10,r": Elst10}, extern_extern_ie


def exchange(cache, jk, do_print=True):
    """
    Computes the E10 exchange (S^2 and S^inf) from a build_sapt_jk_cache datacache.
    """

    if do_print:
        core.print_out("\n  ==> E10 Exchange <== \n\n")

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
    SAB = core.triplet(cache["Cocc_A"], cache["S"], cache["Cocc_B"], True, False, False)
    num_occ = nocc_A + nocc_B

    Sab = core.Matrix(num_occ, num_occ)
    Sab.np[:nocc_A, nocc_A:] = SAB.np
    Sab.np[nocc_A:, :nocc_A] = SAB.np.T
    Sab.np[np.diag_indices_from(Sab.np)] += 1
    Sab.power(-1.0, 1.e-14)
    Sab.np[np.diag_indices_from(Sab.np)] -= 1.0

    Tmo_AA = core.Matrix.from_array(Sab.np[:nocc_A, :nocc_A])
    Tmo_BB = core.Matrix.from_array(Sab.np[nocc_A:, nocc_A:])
    Tmo_AB = core.Matrix.from_array(Sab.np[:nocc_A, nocc_A:])

    T_A = np.dot(cache["Cocc_A"], Tmo_AA).dot(cache["Cocc_A"].np.T)
    T_B = np.dot(cache["Cocc_B"], Tmo_BB).dot(cache["Cocc_B"].np.T)
    T_AB = np.dot(cache["Cocc_A"], Tmo_AB).dot(cache["Cocc_B"].np.T)

    S = cache["S"]

    D_A = cache["D_A"]
    P_A = cache["P_A"]

    D_B = cache["D_B"]
    P_B = cache["P_B"]

    # Compute the J and K matrices
    jk.C_clear()

    jk.C_left_add(cache["Cocc_A"])
    jk.C_right_add(core.doublet(cache["Cocc_A"], Tmo_AA, False, False))

    jk.C_left_add(cache["Cocc_B"])
    jk.C_right_add(core.doublet(cache["Cocc_A"], Tmo_AB, False, False))

    jk.C_left_add(cache["Cocc_A"])
    jk.C_right_add(core.Matrix.chain_dot(P_B, S, cache["Cocc_A"]))

    jk.compute()

    JT_A, JT_AB, Jij = jk.J()
    KT_A, KT_AB, Kij = jk.K()

    # Start S^2
    Exch_s2 = 0.0

    tmp = core.Matrix.chain_dot(D_A, S, D_B, S, P_A)
    Exch_s2 -= 2.0 * w_B.vector_dot(tmp)

    tmp = core.Matrix.chain_dot(D_B, S, D_A, S, P_B)
    Exch_s2 -= 2.0 * w_A.vector_dot(tmp)

    tmp = core.Matrix.chain_dot(P_A, S, D_B)
    Exch_s2 -= 2.0 * Kij.vector_dot(tmp)

    if do_print:
        core.print_out(print_sapt_var("Exch10(S^2) ", Exch_s2, short=True))
        core.print_out("\n")

    # Start Sinf
    Exch10 = 0.0
    Exch10 -= 2.0 * np.vdot(cache["D_A"], cache["K_B"])
    Exch10 += 2.0 * np.vdot(T_A, h_B.np)
    Exch10 += 2.0 * np.vdot(T_B, h_A.np)
    Exch10 += 2.0 * np.vdot(T_AB, h_A.np + h_B.np)
    Exch10 += 4.0 * np.vdot(T_B, JT_AB.np - 0.5 * KT_AB.np)
    Exch10 += 4.0 * np.vdot(T_A, JT_AB.np - 0.5 * KT_AB.np)
    Exch10 += 4.0 * np.vdot(T_B, JT_A.np - 0.5 * KT_A.np)
    Exch10 += 4.0 * np.vdot(T_AB, JT_AB.np - 0.5 * KT_AB.np.T)

    if do_print:
        core.set_variable("Exch10", Exch10)
        core.print_out(print_sapt_var("Exch10", Exch10, short=True))
        core.print_out("\n")

    return {"Exch10(S^2)": Exch_s2, "Exch10": Exch10}


def induction(cache, jk, do_print=True, maxiter=12, conv=1.e-8, do_response=True, Sinf=False, sapt_jk_B=None):
    """
    Compute Ind20 and Exch-Ind20 quantities from a SAPT cache and JK object.
    """

    if do_print:
        core.print_out("\n  ==> E20 Induction <== \n\n")

    # Build Induction and Exchange-Induction potentials
    S = cache["S"]

    D_A = cache["D_A"]
    V_A = cache["V_A"]
    J_A = cache["J_A"]
    K_A = cache["K_A"]

    D_B = cache["D_B"]
    V_B = cache["V_B"]
    J_B = cache["J_B"]
    K_B = cache["K_B"]

    K_O = cache["K_O"]
    J_O = cache["J_O"]

    jk.C_clear()

    jk.C_left_add(core.Matrix.chain_dot(D_B, S, cache["Cocc_A"]))
    jk.C_right_add(cache["Cocc_A"])

    jk.C_left_add(core.Matrix.chain_dot(D_B, S, D_A, S, cache["Cocc_B"]))
    jk.C_right_add(cache["Cocc_B"])

    jk.C_left_add(core.Matrix.chain_dot(D_A, S, D_B, S, cache["Cocc_A"]))
    jk.C_right_add(cache["Cocc_A"])

    jk.compute()

    J_Ot, J_P_B, J_P_A = jk.J()
    K_Ot, K_P_B, K_P_A = jk.K()

    # Exch-Ind Potential A
    EX_A = K_B.clone()
    EX_A.scale(-1.0)
    EX_A.axpy(-2.0, J_O)
    EX_A.axpy(1.0, K_O)
    EX_A.axpy(2.0, J_P_B)

    EX_A.axpy(-1.0, core.Matrix.chain_dot(S, D_B, V_A))
    EX_A.axpy(-2.0, core.Matrix.chain_dot(S, D_B, J_A))
    EX_A.axpy(1.0, core.Matrix.chain_dot(S, D_B, K_A))
    EX_A.axpy(1.0, core.Matrix.chain_dot(S, D_B, S, D_A, V_B))
    EX_A.axpy(2.0, core.Matrix.chain_dot(S, D_B, S, D_A, J_B))
    EX_A.axpy(1.0, core.Matrix.chain_dot(S, D_B, V_A, D_B, S))
    EX_A.axpy(2.0, core.Matrix.chain_dot(S, D_B, J_A, D_B, S))
    EX_A.axpy(-1.0, core.Matrix.chain_dot(S, D_B, K_O, trans=[False, False, True]))

    EX_A.axpy(-1.0, core.Matrix.chain_dot(V_B, D_B, S))
    EX_A.axpy(-2.0, core.Matrix.chain_dot(J_B, D_B, S))
    EX_A.axpy(1.0, core.Matrix.chain_dot(K_B, D_B, S))
    EX_A.axpy(1.0, core.Matrix.chain_dot(V_B, D_A, S, D_B, S))
    EX_A.axpy(2.0, core.Matrix.chain_dot(J_B, D_A, S, D_B, S))
    EX_A.axpy(-1.0, core.Matrix.chain_dot(K_O, D_B, S))

    EX_A = core.Matrix.chain_dot(cache["Cocc_A"], EX_A, cache["Cvir_A"], trans=[True, False, False])

    # Exch-Ind Potential B
    EX_B = K_A.clone()
    EX_B.scale(-1.0)
    EX_B.axpy(-2.0, J_O)
    EX_B.axpy(1.0, K_O.transpose())
    EX_B.axpy(2.0, J_P_A)

    EX_B.axpy(-1.0, core.Matrix.chain_dot(S, D_A, V_B))
    EX_B.axpy(-2.0, core.Matrix.chain_dot(S, D_A, J_B))
    EX_B.axpy(1.0, core.Matrix.chain_dot(S, D_A, K_B))
    EX_B.axpy(1.0, core.Matrix.chain_dot(S, D_A, S, D_B, V_A))
    EX_B.axpy(2.0, core.Matrix.chain_dot(S, D_A, S, D_B, J_A))
    EX_B.axpy(1.0, core.Matrix.chain_dot(S, D_A, V_B, D_A, S))
    EX_B.axpy(2.0, core.Matrix.chain_dot(S, D_A, J_B, D_A, S))
    EX_B.axpy(-1.0, core.Matrix.chain_dot(S, D_A, K_O))

    EX_B.axpy(-1.0, core.Matrix.chain_dot(V_A, D_A, S))
    EX_B.axpy(-2.0, core.Matrix.chain_dot(J_A, D_A, S))
    EX_B.axpy(1.0, core.Matrix.chain_dot(K_A, D_A, S))
    EX_B.axpy(1.0, core.Matrix.chain_dot(V_A, D_B, S, D_A, S))
    EX_B.axpy(2.0, core.Matrix.chain_dot(J_A, D_B, S, D_A, S))
    EX_B.axpy(-1.0, core.Matrix.chain_dot(K_O, D_A, S, trans=[True, False, False]))

    EX_B = core.Matrix.chain_dot(cache["Cocc_B"], EX_B, cache["Cvir_B"], trans=[True, False, False])

    # Build electrostatic potenital
    w_A = cache["V_A"].clone()
    w_A.axpy(2.0, cache["J_A"])

    w_B = cache["V_B"].clone()
    w_B.axpy(2.0, cache["J_B"])

    w_B_MOA = core.triplet(cache["Cocc_A"], w_B, cache["Cvir_A"], True, False, False)
    w_A_MOB = core.triplet(cache["Cocc_B"], w_A, cache["Cvir_B"], True, False, False)

    # Do uncoupled
    core.print_out("   => Uncoupled Induction <= \n\n")
    unc_x_B_MOA = w_B_MOA.clone()
    unc_x_B_MOA.np[:] /= (cache["eps_occ_A"].np.reshape(-1, 1) - cache["eps_vir_A"].np)
    unc_x_A_MOB = w_A_MOB.clone()
    unc_x_A_MOB.np[:] /= (cache["eps_occ_B"].np.reshape(-1, 1) - cache["eps_vir_B"].np)

    unc_ind_ab = 2.0 * unc_x_B_MOA.vector_dot(w_B_MOA)
    unc_ind_ba = 2.0 * unc_x_A_MOB.vector_dot(w_A_MOB)
    unc_indexch_ab = 2.0 * unc_x_B_MOA.vector_dot(EX_A)
    unc_indexch_ba = 2.0 * unc_x_A_MOB.vector_dot(EX_B)

    ret = {}
    ret["Ind20,u (A<-B)"] = unc_ind_ab
    ret["Ind20,u (A->B)"] = unc_ind_ba
    ret["Ind20,u"] = unc_ind_ab + unc_ind_ba
    ret["Exch-Ind20,u (A<-B)"] = unc_indexch_ab
    ret["Exch-Ind20,u (A->B)"] = unc_indexch_ba
    ret["Exch-Ind20,u"] = unc_indexch_ba + unc_indexch_ab

    plist = [
        "Ind20,u (A<-B)", "Ind20,u (A->B)", "Ind20,u", "Exch-Ind20,u (A<-B)", "Exch-Ind20,u (A->B)", "Exch-Ind20,u"
    ]

    if do_print:
        for name in plist:
            # core.set_variable(name, ret[name])
            core.print_out(print_sapt_var(name, ret[name], short=True))
            core.print_out("\n")

    # Exch-Ind without S^2
    if Sinf:
        nocc_A = cache["Cocc_A"].shape[1]
        nocc_B = cache["Cocc_B"].shape[1]
        SAB = core.triplet(cache["Cocc_A"], cache["S"], cache["Cocc_B"], True, False, False)
        num_occ = nocc_A + nocc_B

        Sab = core.Matrix(num_occ, num_occ)
        Sab.np[:nocc_A, nocc_A:] = SAB.np
        Sab.np[nocc_A:, :nocc_A] = SAB.np.T
        Sab.np[np.diag_indices_from(Sab.np)] += 1
        Sab.power(-1.0, 1.e-14)

        Tmo_AA = core.Matrix.from_array(Sab.np[:nocc_A, :nocc_A])
        Tmo_BB = core.Matrix.from_array(Sab.np[nocc_A:, nocc_A:])
        Tmo_AB = core.Matrix.from_array(Sab.np[:nocc_A, nocc_A:])

        T_A = core.triplet(cache["Cocc_A"], Tmo_AA, cache["Cocc_A"], False, False, True)
        T_B = core.triplet(cache["Cocc_B"], Tmo_BB, cache["Cocc_B"], False, False, True)
        T_AB = core.triplet(cache["Cocc_A"], Tmo_AB, cache["Cocc_B"], False, False, True)

        sT_A = core.Matrix.chain_dot(cache["Cvir_A"],
                                     unc_x_B_MOA,
                                     Tmo_AA,
                                     cache["Cocc_A"],
                                     trans=[False, True, False, True])
        sT_B = core.Matrix.chain_dot(cache["Cvir_B"],
                                     unc_x_A_MOB,
                                     Tmo_BB,
                                     cache["Cocc_B"],
                                     trans=[False, True, False, True])
        sT_AB = core.Matrix.chain_dot(cache["Cvir_A"],
                                      unc_x_B_MOA,
                                      Tmo_AB,
                                      cache["Cocc_B"],
                                      trans=[False, True, False, True])
        sT_BA = core.Matrix.chain_dot(cache["Cvir_B"],
                                      unc_x_A_MOB,
                                      Tmo_AB,
                                      cache["Cocc_A"],
                                      trans=[False, True, True, True])

        jk.C_clear()

        jk.C_left_add(core.Matrix.chain_dot(cache["Cocc_A"], Tmo_AA))
        jk.C_right_add(cache["Cocc_A"])

        jk.C_left_add(core.Matrix.chain_dot(cache["Cocc_B"], Tmo_BB))
        jk.C_right_add(cache["Cocc_B"])

        jk.C_left_add(core.Matrix.chain_dot(cache["Cocc_A"], Tmo_AB))
        jk.C_right_add(cache["Cocc_B"])

        jk.compute()

        J_AA_inf, J_BB_inf, J_AB_inf = jk.J()
        K_AA_inf, K_BB_inf, K_AB_inf = jk.K()

        # A <- B
        EX_AA_inf = V_B.clone()
        EX_AA_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_AB, V_B, trans=[False, True, False]))
        EX_AA_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_B, V_B))
        EX_AA_inf.axpy(2.00, J_AB_inf)
        EX_AA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_AB_inf, trans=[False, True, False]))
        EX_AA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_B, J_AB_inf))
        EX_AA_inf.axpy(2.00, J_BB_inf)
        EX_AA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_BB_inf, trans=[False, True, False]))
        EX_AA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_B, J_BB_inf))
        EX_AA_inf.axpy(-1.00, K_AB_inf.transpose())
        EX_AA_inf.axpy(1.00, core.Matrix.chain_dot(S, T_AB, K_AB_inf, trans=[False, True, True]))
        EX_AA_inf.axpy(1.00, core.Matrix.chain_dot(S, T_B, K_AB_inf, trans=[False, False, True]))
        EX_AA_inf.axpy(-1.00, K_BB_inf)
        EX_AA_inf.axpy(1.00, core.Matrix.chain_dot(S, T_AB, K_BB_inf, trans=[False, True, False]))
        EX_AA_inf.axpy(1.00, core.Matrix.chain_dot(S, T_B, K_BB_inf))

        EX_AB_inf = V_A.clone()
        EX_AB_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_AB, V_A, trans=[False, True, False]))
        EX_AB_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_B, V_A))
        EX_AB_inf.axpy(2.00, J_AA_inf)
        EX_AB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_AA_inf, trans=[False, True, False]))
        EX_AB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_B, J_AA_inf))
        EX_AB_inf.axpy(2.00, J_AB_inf)
        EX_AB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_AB_inf, trans=[False, True, False]))
        EX_AB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_B, J_AB_inf))
        EX_AB_inf.axpy(-1.00, K_AA_inf)
        EX_AB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_AB, K_AA_inf, trans=[False, True, False]))
        EX_AB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_B, K_AA_inf))
        EX_AB_inf.axpy(-1.00, K_AB_inf)
        EX_AB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_AB, K_AB_inf, trans=[False, True, False]))
        EX_AB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_B, K_AB_inf))

        # B <- A
        EX_BB_inf = V_A.clone()
        EX_BB_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_AB, V_A))
        EX_BB_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_A, V_A))
        EX_BB_inf.axpy(2.00, J_AB_inf)
        EX_BB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_AB_inf))
        EX_BB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_A, J_AB_inf))
        EX_BB_inf.axpy(2.00, J_AA_inf)
        EX_BB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_AA_inf))
        EX_BB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_A, J_AA_inf))
        EX_BB_inf.axpy(-1.00, K_AB_inf)
        EX_BB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_AB, K_AB_inf))
        EX_BB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_A, K_AB_inf))
        EX_BB_inf.axpy(-1.00, K_AA_inf)
        EX_BB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_AB, K_AA_inf))
        EX_BB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_A, K_AA_inf))

        EX_BA_inf = V_B.clone()
        EX_BA_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_AB, V_B))
        EX_BA_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_A, V_B))
        EX_BA_inf.axpy(2.00, J_BB_inf)
        EX_BA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_BB_inf))
        EX_BA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_A, J_BB_inf))
        EX_BA_inf.axpy(2.00, J_AB_inf)
        EX_BA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_AB_inf))
        EX_BA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_A, J_AB_inf))
        EX_BA_inf.axpy(-1.00, K_BB_inf)
        EX_BA_inf.axpy(1.00, core.Matrix.chain_dot(S, T_AB, K_BB_inf))
        EX_BA_inf.axpy(1.00, core.Matrix.chain_dot(S, T_A, K_BB_inf))
        EX_BA_inf.axpy(-1.00, K_AB_inf.transpose())
        EX_BA_inf.axpy(1.00, core.Matrix.chain_dot(S, T_AB, K_AB_inf, trans=[False, False, True]))
        EX_BA_inf.axpy(1.00, core.Matrix.chain_dot(S, T_A, K_AB_inf, trans=[False, False, True]))

        unc_ind_ab_total = 2.0 * (sT_A.vector_dot(EX_AA_inf) + sT_AB.vector_dot(EX_AB_inf))
        unc_ind_ba_total = 2.0 * (sT_B.vector_dot(EX_BB_inf) + sT_BA.vector_dot(EX_BA_inf))
        unc_indexch_ab_inf = unc_ind_ab_total - unc_ind_ab
        unc_indexch_ba_inf = unc_ind_ba_total - unc_ind_ba

        ret["Exch-Ind20,u (A<-B) (S^inf)"] = unc_indexch_ab_inf
        ret["Exch-Ind20,u (A->B) (S^inf)"] = unc_indexch_ba_inf
        ret["Exch-Ind20,u (S^inf)"] = unc_indexch_ba_inf + unc_indexch_ab_inf

        if do_print:
            for name in plist[3:]:
                name = name + ' (S^inf)'

                core.print_out(print_sapt_var(name, ret[name], short=True))
                core.print_out("\n")

    # Do coupled
    if do_response:
        core.print_out("\n   => Coupled Induction <= \n\n")

        cphf_r_convergence = core.get_option("SAPT", "CPHF_R_CONVERGENCE")

        x_B_MOA, x_A_MOB = _sapt_cpscf_solve(cache, jk, w_B_MOA, w_A_MOB, 20, cphf_r_convergence, sapt_jk_B=sapt_jk_B)

        ind_ab = 2.0 * x_B_MOA.vector_dot(w_B_MOA)
        ind_ba = 2.0 * x_A_MOB.vector_dot(w_A_MOB)
        indexch_ab = 2.0 * x_B_MOA.vector_dot(EX_A)
        indexch_ba = 2.0 * x_A_MOB.vector_dot(EX_B)

        ret["Ind20,r (A<-B)"] = ind_ab
        ret["Ind20,r (A->B)"] = ind_ba
        ret["Ind20,r"] = ind_ab + ind_ba
        ret["Exch-Ind20,r (A<-B)"] = indexch_ab
        ret["Exch-Ind20,r (A->B)"] = indexch_ba
        ret["Exch-Ind20,r"] = indexch_ba + indexch_ab

        if do_print:
            core.print_out("\n")
            for name in plist:
                name = name.replace(",u", ",r")

                # core.set_variable(name, ret[name])
                core.print_out(print_sapt_var(name, ret[name], short=True))
                core.print_out("\n")

        # Exch-Ind without S^2
        if Sinf:
            cT_A = core.Matrix.chain_dot(cache["Cvir_A"],
                                         x_B_MOA,
                                         Tmo_AA,
                                         cache["Cocc_A"],
                                         trans=[False, True, False, True])
            cT_B = core.Matrix.chain_dot(cache["Cvir_B"],
                                         x_A_MOB,
                                         Tmo_BB,
                                         cache["Cocc_B"],
                                         trans=[False, True, False, True])
            cT_AB = core.Matrix.chain_dot(cache["Cvir_A"],
                                          x_B_MOA,
                                          Tmo_AB,
                                          cache["Cocc_B"],
                                          trans=[False, True, False, True])
            cT_BA = core.Matrix.chain_dot(cache["Cvir_B"],
                                          x_A_MOB,
                                          Tmo_AB,
                                          cache["Cocc_A"],
                                          trans=[False, True, True, True])

            ind_ab_total = 2.0 * (cT_A.vector_dot(EX_AA_inf) + cT_AB.vector_dot(EX_AB_inf))
            ind_ba_total = 2.0 * (cT_B.vector_dot(EX_BB_inf) + cT_BA.vector_dot(EX_BA_inf))
            indexch_ab_inf = ind_ab_total - ind_ab
            indexch_ba_inf = ind_ba_total - ind_ba

            ret["Exch-Ind20,r (A<-B) (S^inf)"] = indexch_ab_inf
            ret["Exch-Ind20,r (A->B) (S^inf)"] = indexch_ba_inf
            ret["Exch-Ind20,r (S^inf)"] = indexch_ba_inf + indexch_ab_inf

            if do_print:
                for name in plist[3:]:
                    name = name.replace(",u", ",r") + ' (S^inf)'

                    core.print_out(print_sapt_var(name, ret[name], short=True))
                    core.print_out("\n")

    return ret


def _sapt_cpscf_solve(cache, jk, rhsA, rhsB, maxiter, conv, sapt_jk_B=None):
    """
    Solve the SAPT CPHF (or CPKS) equations.
    """

    cache["wfn_A"].set_jk(jk)
    if sapt_jk_B:
        cache["wfn_B"].set_jk(sapt_jk_B)
    else:
        cache["wfn_B"].set_jk(jk)

    # Make a preconditioner function
    P_A = core.Matrix(cache["eps_occ_A"].shape[0], cache["eps_vir_A"].shape[0])
    P_A.np[:] = (cache["eps_occ_A"].np.reshape(-1, 1) - cache["eps_vir_A"].np)

    P_B = core.Matrix(cache["eps_occ_B"].shape[0], cache["eps_vir_B"].shape[0])
    P_B.np[:] = (cache["eps_occ_B"].np.reshape(-1, 1) - cache["eps_vir_B"].np)

    # Preconditioner function
    def apply_precon(x_vec, act_mask):
        if act_mask[0]:
            pA = x_vec[0].clone()
            pA.apply_denominator(P_A)
        else:
            pA = False

        if act_mask[1]:
            pB = x_vec[1].clone()
            pB.apply_denominator(P_B)
        else:
            pB = False

        return [pA, pB]

    # Hx function
    def hessian_vec(x_vec, act_mask):
        if act_mask[0]:
            xA = cache["wfn_A"].cphf_Hx([x_vec[0]])[0]
        else:
            xA = False

        if act_mask[1]:
            xB = cache["wfn_B"].cphf_Hx([x_vec[1]])[0]
        else:
            xB = False

        return [xA, xB]

    # Manipulate the printing
    sep_size = 51
    core.print_out("   " + ("-" * sep_size) + "\n")
    core.print_out("   " + "SAPT Coupled Induction Solver".center(sep_size) + "\n")
    core.print_out("   " + ("-" * sep_size) + "\n")
    core.print_out("    Maxiter             = %11d\n" % maxiter)
    core.print_out("    Convergence         = %11.3E\n" % conv)
    core.print_out("   " + ("-" * sep_size) + "\n")

    tstart = time.time()
    core.print_out("     %4s %12s     %12s     %9s\n" % ("Iter", "(A<-B)", "(B->A)", "Time [s]"))
    core.print_out("   " + ("-" * sep_size) + "\n")

    start_resid = [rhsA.sum_of_squares(), rhsB.sum_of_squares()]

    # print function
    def pfunc(niter, x_vec, r_vec):
        if niter == 0:
            niter = "Guess"
        else:
            niter = ("%5d" % niter)

        # Compute IndAB
        valA = (r_vec[0].sum_of_squares() / start_resid[0])**0.5
        if valA < conv:
            cA = "*"
        else:
            cA = " "

        # Compute IndBA
        valB = (r_vec[1].sum_of_squares() / start_resid[1])**0.5
        if valB < conv:
            cB = "*"
        else:
            cB = " "

        core.print_out("    %5s %15.6e%1s %15.6e%1s %9d\n" % (niter, valA, cA, valB, cB, time.time() - tstart))
        return [valA, valB]

    # Compute the solver
    vecs, resid = solvers.cg_solver([rhsA, rhsB],
                                    hessian_vec,
                                    apply_precon,
                                    maxiter=maxiter,
                                    rcond=conv,
                                    printlvl=0,
                                    printer=pfunc)
    core.print_out("   " + ("-" * sep_size) + "\n")

    return vecs
