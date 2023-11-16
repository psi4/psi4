#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
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

from collections import OrderedDict
from typing import Dict

import numpy as np

from psi4 import core

from ...p4util import solvers
from ...p4util.exceptions import *
from .sapt_util import print_sapt_var

__all__ = ["compute_sapt_sf"]


def _sf_compute_JK(jk, Cleft, Cright, rotation=None):
    """
    A specialized JK computer class for terms that arrise from SF-SAPT.

    The density is computed as (Cl_mu,i rotation_ij Cr_nu,j) where the rotation
    is an arbitrary perturbation on the density.
    """

    # Handle both list and single value input
    return_single = False
    if not isinstance(Cleft, (list, tuple)):
        Cleft = [Cleft]
        return_single = True
    if not isinstance(Cright, (list, tuple)):
        Cright = [Cright]
        return_single = True
    if (not isinstance(rotation, (list, tuple))) and (rotation is not None):
        rotation = [rotation]
        return_single = True

    if len(Cleft) != len(Cright):
        raise ValidationError("Cleft list is not the same length as Cright list")

    jk.C_clear()

    zero_append = []
    num_compute = 0

    for num in range(len(Cleft)):
        Cl = Cleft[num]
        Cr = Cright[num]

        if (Cr.shape[1] == 0) or (Cl.shape[1] == 0):
            zero_append.append(num)
            continue

        if (rotation is not None) and (rotation[num] is not None):
            mol = Cl.shape[1]
            mor = Cr.shape[1]

            if (rotation[num].shape[0] != mol) or (rotation[num].shape[1] != mor):
                raise ValidationError("_sf_compute_JK: Tensor size does not match Cl (%d) /Cr (%d) : %s" %
                                (mol, mor, str(rotation[num].shape)))

            # Figure out the small MO index to contract to
            if mol < mor:
                Cl = np.dot(Cl, rotation[num])
            else:
                Cr = np.dot(Cr, rotation[num].T)

        Cl = core.Matrix.from_array(Cl)
        Cr = core.Matrix.from_array(Cr)

        jk.C_left_add(Cl)
        jk.C_right_add(Cr)
        num_compute += 1

    jk.compute()

    J_list = []
    K_list = []
    for num in range(num_compute):
        J_list.append(np.array(jk.J()[num]))
        K_list.append(np.array(jk.K()[num]))

    jk.C_clear()

    nbf = J_list[0].shape[0]
    zero_mat = np.zeros((nbf, nbf))
    for num in zero_append:
        J_list.insert(num, zero_mat)
        K_list.insert(num, zero_mat)

    if return_single:
        return J_list[0], K_list[0]
    else:
        return J_list, K_list


def _chain_dot(*dot_list):
    """
    A simple chain dot function unpacked from *args.
    """
    result = dot_list[0]
    for x in range(len(dot_list) - 1):
        result = np.dot(result, dot_list[x + 1])
    return result


def compute_sapt_sf(
    dimer: core.Molecule,
    jk: core.JK,
    wfn_A: core.Wavefunction,
    wfn_B: core.Wavefunction,
    do_print: bool = True,
) -> Dict[str, float]:
    """
    Computes Elst and Spin-Flip SAPT0 for ROHF wavefunctions
    """

    if do_print:
        core.print_out("\n  ==> Preparing SF-SAPT Data Cache <== \n\n")
        jk.print_header()

    ### Build intermediates

    # Pull out Wavefunction A quantities
    ndocc_A = wfn_A.doccpi().sum()
    nsocc_A = wfn_A.soccpi().sum()

    Cocc_A = np.asarray(wfn_A.Ca_subset("AO", "OCC"))
    Ci = Cocc_A[:, :ndocc_A]
    Ca = Cocc_A[:, ndocc_A:]
    Pi = np.dot(Ci, Ci.T)
    Pa = np.dot(Ca, Ca.T)

    mints = core.MintsHelper(wfn_A.basisset())
    V_A = mints.ao_potential()

    # Pull out Wavefunction B quantities
    ndocc_B = wfn_B.doccpi().sum()
    nsocc_B = wfn_B.soccpi().sum()

    Cocc_B = np.asarray(wfn_B.Ca_subset("AO", "OCC"))
    Cj = Cocc_B[:, :ndocc_B]
    Cb = Cocc_B[:, ndocc_B:]
    Pj = np.dot(Cj, Cj.T)
    Pb = np.dot(Cb, Cb.T)

    mints = core.MintsHelper(wfn_B.basisset())
    V_B = mints.ao_potential()

    # Pull out generic quantities
    S = np.asarray(wfn_A.S())

    intermonomer_nuclear_repulsion = dimer.nuclear_repulsion_energy()
    intermonomer_nuclear_repulsion -= wfn_A.molecule().nuclear_repulsion_energy()
    intermonomer_nuclear_repulsion -= wfn_B.molecule().nuclear_repulsion_energy()

    num_el_A = (2 * ndocc_A + nsocc_A)
    num_el_B = (2 * ndocc_B + nsocc_B)

    ### Build JK Terms
    if do_print:
        core.print_out("\n  ==> Computing required JK matrices <== \n\n")

    # Writen so that we can reorganize order to save on DF-JK cost.
    pairs = [("ii", Ci, None, Ci),
             ("ij", Ci, _chain_dot(Ci.T, S, Cj), Cj),
             ("jj", Cj, None, Cj),
             ("aa", Ca, None, Ca),
             ("aj", Ca, _chain_dot(Ca.T, S, Cj), Cj),
             ("ib", Ci, _chain_dot(Ci.T, S, Cb), Cb),
             ("bb", Cb, None, Cb),
             ("ab", Ca, _chain_dot(Ca.T, S, Cb), Cb)]

    # Reorganize
    names = [x[0] for x in pairs]
    Cleft = [x[1] for x in pairs]
    rotations = [x[2] for x in pairs]
    Cright = [x[3] for x in pairs]

    tmp_J, tmp_K = _sf_compute_JK(jk, Cleft, Cright, rotations)

    J = {key: val for key, val in zip(names, tmp_J)}
    K = {key: val for key, val in zip(names, tmp_K)}

    ### Compute Terms
    if do_print:
        core.print_out("\n  ==> Computing Spin-Flip Exchange and Electrostatics <== \n\n")

    w_A = V_A + 2 * J["ii"] + J["aa"]
    w_B = V_B + 2 * J["jj"] + J["bb"]

    h_Aa = V_A + 2 * J["ii"] + J["aa"] - K["ii"] - K["aa"]
    h_Ab = V_A + 2 * J["ii"] + J["aa"] - K["ii"]

    h_Ba = V_B + 2 * J["jj"] + J["bb"] - K["jj"]
    h_Bb = V_B + 2 * J["jj"] + J["bb"] - K["jj"] - K["bb"]

    ### Build electrostatics

    # socc/socc term
    two_el_repulsion = np.vdot(Pa, J["bb"])
    attractive_a = np.vdot(V_A, Pb) * nsocc_A / num_el_A
    attractive_b = np.vdot(V_B, Pa) * nsocc_B / num_el_B
    nuclear_repulsion = intermonomer_nuclear_repulsion * nsocc_A * nsocc_B / (num_el_A * num_el_B)
    elst_abab = two_el_repulsion + attractive_a + attractive_b + nuclear_repulsion

    # docc/socc term
    two_el_repulsion = np.vdot(Pi, J["bb"])
    attractive_a = np.vdot(V_A, Pb) * ndocc_A / num_el_A
    attractive_b = np.vdot(V_B, Pi) * nsocc_B / num_el_B
    nuclear_repulsion = intermonomer_nuclear_repulsion * ndocc_A * nsocc_B / (num_el_A * num_el_B)
    elst_ibib = 2 * (two_el_repulsion + attractive_a + attractive_b + nuclear_repulsion)

    # socc/docc term
    two_el_repulsion = np.vdot(Pa, J["jj"])
    attractive_a = np.vdot(V_A, Pj) * nsocc_A / num_el_A
    attractive_b = np.vdot(V_B, Pa) * ndocc_B / num_el_B
    nuclear_repulsion = intermonomer_nuclear_repulsion * nsocc_A * ndocc_B / (num_el_A * num_el_B)
    elst_jaja = 2 * (two_el_repulsion + attractive_a + attractive_b + nuclear_repulsion)

    # docc/docc term
    two_el_repulsion = np.vdot(Pi, J["jj"])
    attractive_a = np.vdot(V_A, Pj) * ndocc_A / num_el_A
    attractive_b = np.vdot(V_B, Pi) * ndocc_B / num_el_B
    nuclear_repulsion = intermonomer_nuclear_repulsion * ndocc_A * ndocc_B / (num_el_A * num_el_B)
    elst_ijij = 4 * (two_el_repulsion + attractive_a + attractive_b + nuclear_repulsion)

    elst = elst_abab + elst_ibib + elst_jaja + elst_ijij
    # print(print_sapt_var("Elst,10", elst))

    ### Start diagonal exchange

    exch_diag = 0.0
    exch_diag -= np.vdot(Pj, 2 * K["ii"] + K["aa"])
    exch_diag -= np.vdot(Pb, K["ii"])
    exch_diag -= np.vdot(_chain_dot(Pi, S, Pj), (h_Aa + h_Ab + h_Ba + h_Bb))
    exch_diag -= np.vdot(_chain_dot(Pa, S, Pj), (h_Aa + h_Ba))
    exch_diag -= np.vdot(_chain_dot(Pi, S, Pb), (h_Ab + h_Bb))

    exch_diag += 2.0 * np.vdot(_chain_dot(Pj, S, Pi, S, Pb), w_A)
    exch_diag += 2.0 * np.vdot(_chain_dot(Pj, S, Pi, S, Pj), w_A)
    exch_diag += np.vdot(_chain_dot(Pb, S, Pi, S, Pb), w_A)
    exch_diag += np.vdot(_chain_dot(Pj, S, Pa, S, Pj), w_A)

    exch_diag += 2.0 * np.vdot(_chain_dot(Pi, S, Pj, S, Pi), w_B)
    exch_diag += 2.0 * np.vdot(_chain_dot(Pi, S, Pj, S, Pa), w_B)
    exch_diag += np.vdot(_chain_dot(Pi, S, Pb, S, Pi), w_B)
    exch_diag += np.vdot(_chain_dot(Pa, S, Pj, S, Pa), w_B)

    exch_diag -= 2.0 * np.vdot(_chain_dot(Pi, S, Pj), K["ij"])
    exch_diag -= 2.0 * np.vdot(_chain_dot(Pa, S, Pj), K["ij"])
    exch_diag -= 2.0 * np.vdot(_chain_dot(Pi, S, Pb), K["ij"])

    exch_diag -= np.vdot(_chain_dot(Pa, S, Pj), K["aj"])
    exch_diag -= np.vdot(_chain_dot(Pi, S, Pb), K["ib"])
    # print(print_sapt_var("Exch10,offdiagonal", exch_diag))

    ### Start off-diagonal exchange

    exch_offdiag = 0.0
    exch_offdiag -= np.vdot(Pb, K["aa"])
    exch_offdiag -= np.vdot(_chain_dot(Pa, S, Pb), (h_Aa + h_Bb))
    exch_offdiag += np.vdot(_chain_dot(Pa, S, Pj), K["bb"])
    exch_offdiag += np.vdot(_chain_dot(Pi, S, Pb), K["aa"])

    exch_offdiag += 2.0 * np.vdot(_chain_dot(Pj, S, Pa, S, Pb), w_A)
    exch_offdiag += np.vdot(_chain_dot(Pb, S, Pa, S, Pb), w_A)

    exch_offdiag += 2.0 * np.vdot(_chain_dot(Pi, S, Pb, S, Pa), w_B)
    exch_offdiag += np.vdot(_chain_dot(Pa, S, Pb, S, Pa), w_B)

    exch_offdiag -= 2.0 * np.vdot(_chain_dot(Pa, S, Pb), K["ij"])
    exch_offdiag -= 2.0 * np.vdot(_chain_dot(Pa, S, Pb), K["ib"])
    exch_offdiag -= 2.0 * np.vdot(_chain_dot(Pa, S, Pj), K["ab"])
    exch_offdiag -= 2.0 * np.vdot(_chain_dot(Pa, S, Pj), K["ib"])

    exch_offdiag -= np.vdot(_chain_dot(Pa, S, Pb), K["ab"])
    # print(print_sapt_var("Exch10,off-diagonal", exch_offdiag))
    # print(print_sapt_var("Exch10(S^2)", exch_offdiag + exch_diag))

    ret_values = OrderedDict({
        "Elst10": elst,
        "Exch10(S^2) [diagonal]": exch_diag,
        "Exch10(S^2) [off-diagonal]": exch_offdiag,
        "Exch10(S^2) [highspin]": exch_offdiag + exch_diag,
    })

    return ret_values
