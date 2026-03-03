#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
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

from ... import p4util
from ...constants import constants
from ...p4util.exceptions import ValidationError
from .dftd4_r4r2 import r4r2_dftd4


def dftd4_c6_intermolecular_dispersion(
    atomic_numbers: np.ndarray,  # shape (n_atoms,)
    geometry: np.ndarray,  # shape (n_atoms, 3)
    C6s: np.ndarray,  # shape (n_atoms, n_atoms)
    monAs: np.ndarray,  # shape (nA,)
    monBs: np.ndarray,  # shape (nB,)
    params: list,  # [s6, s8, a1, a2]
) -> dict:
    s6, s8, a1, a2 = params
    energy = 0.0
    pairwise_energies = np.zeros_like(C6s)
    for A in monAs:
        el1 = atomic_numbers[A]
        Q_A = np.sqrt(0.5 * np.sqrt(el1) * r4r2_dftd4[el1 - 1])

        for B in monBs:
            el2 = atomic_numbers[B]
            Q_B = np.sqrt(0.5 * np.sqrt(el2) * r4r2_dftd4[el2 - 1])

            rrij = 3.0 * Q_A * Q_B
            r0ij = a1 * np.sqrt(rrij) + a2

            rij_vec = geometry[A] - geometry[B]
            dis2 = np.dot(rij_vec, rij_vec)

            t6 = 1.0 / (dis2**3 + r0ij**6)
            t8 = 1.0 / (dis2**4 + r0ij**8)

            edisp = s6 * t6 + s8 * rrij * t8

            de = -C6s[A, B] * edisp
            energy += de
            pairwise_energies[A, B] = de
    return {
        "D4 IE": energy,
        "FSAPT_EMPIRICAL_DISP": pairwise_energies,
    }


def sapt_dft_d4_interaction_energy(
    sapt_dimer, monomerA, monomerB, dftd4_functional_name, d4_type, data
):
    if d4_type == "supermolecular":
        core.print_out(
            "         "
            + "Supermolecular D4 Interaction Energy E_IE = E_IJ - E_I - E_J".center(58)
            + "\n"
        )
        params = {
            "s6": 1.00000000e00,
            "s8": 1.20417708e00,
            "a1": 9.09018333e-01,
            "a2": 3.24886637e-10,
            "s9": 0.00000000e00,
        }
        dimer_d4, _ = sapt_dimer.run_dftd4(dashparam=params, dashlvl="d4bjeeqatm")
        data["D4 DIMER"] = dimer_d4
        monA_d4, _ = monomerA.run_dftd4(dashparam=params, dashlvl="d4bjeeqatm")
        data["D4 MONOMER A"] = monA_d4
        monB_d4, _ = monomerB.run_dftd4(dashparam=params, dashlvl="d4bjeeqatm")
        data["D4 MONOMER B"] = monB_d4
        data["D4 IE"] = dimer_d4 - monA_d4 - monB_d4
    elif d4_type == "gd4_supermolecular":
        # This uses the default Grimme -D4 parameters for the given
        # functional with BJ damping and ATM three-body terms. Only use
        # this option in compbination with another baseline form of
        # dispersion like the delta_DFT correction:
        # "SAPT_DFT_DO_DDFT = True"
        core.print_out(
            "         "
            + "Supermolecular GD4(BJ)+ATM Interaction Energy E_IE = E_IJ - E_I - E_J".center(
                58
            )
            + "\n"
        )
        dimer_d4, _ = sapt_dimer.run_dftd4(
            func=dftd4_functional_name, dashlvl="d4bjeeqatm"
        )
        data["D4 DIMER"] = dimer_d4
        monA_d4, _ = monomerA.run_dftd4(func=dftd4_functional_name, dashlvl="d4bjeeqatm")
        data["D4 MONOMER A"] = monA_d4
        monB_d4, _ = monomerB.run_dftd4(func=dftd4_functional_name, dashlvl="d4bjeeqatm")
        data["D4 MONOMER B"] = monB_d4
        data["D4 IE"] = dimer_d4 - monA_d4 - monB_d4
    elif d4_type == "intermolecular":
        dimer_d4, _ = sapt_dimer.run_dftd4(dftd4_functional_name, property=True)
        geom, _, _, elez, _ = sapt_dimer.to_arrays()
        elez = np.array(elez, dtype=np.int32)
        monAs = np.array([i for i in range(monomerA.natom()) if monomerA.Z(i) > 0])
        monBs = np.array([i for i in range(monomerB.natom()) if monomerB.Z(i) > 0])
        if dftd4_functional_name.lower() == "pbe0":
            # SAPT_DFT_pbe0_adz_3_IE_supra
            params = [1.0, 0.89529649, -0.82043591, 0.03264695]
        elif dftd4_functional_name.lower() == "hf":
            # SAPT0_adz_3_IE_2B_BJ_inter
            params = [1.0, 0.56063068, 0.65540802, 1.06422537]
            # sadz (supermolecular damping...)
            # params = [1.0, 0.83055196, 0.70628586, 1.12379695]
        else:
            raise ValueError(
                "SAPT(DFT): DFTD4 parameters not available for functional %s"
                % dftd4_functional_name
            )
        v = dftd4_c6_intermolecular_dispersion(
            elez,
            geom,
            C6s=core.variable("DFTD4 C6 COEFFICIENTS").np,
            monAs=monAs,
            monBs=monBs,
            params=params,
        )
        data["D4 IE"] = v["D4 IE"]
        data["FSAPT_EMPIRICAL_DISP"] = v["FSAPT_EMPIRICAL_DISP"]
    else:
        raise ValueError(
            "SAPT(DFT): d4_type must be 'supermolecular' or 'intermolecular'."
        )
    return data


def sapt_dft_d3_interaction_energy(
    sapt_dimer, monomerA, monomerB, dftd3_functional_name, d3_type, data
):
    if d3_type == "supermolecular":
        core.print_out(
            "         "
            + "Supermolecular D3 Interaction Energy E_IE = E_IJ - E_I - E_J".center(58)
            + "\n"
        )
        params = {
            "s6": 1.00000000e00,
            "s8": 1.20317708e00,
            "a1": 9.09018333e-01,
            "a2": 3.23886637e-10,
            "s9": 0.00000000e00,
        }
        dimer_d3, _ = sapt_dimer.run_sdftd3(func=dftd3_functional_name, dashlvl="d3mbj", verbose=True)
        data["D3 DIMER"] = dimer_d3
        monA_d3, _ = monomerA.run_sdftd3(func=dftd3_functional_name, dashlvl="d3mbj", verbose=True)
        data["D3 MONOMER A"] = monA_d3
        monB_d3, _ = monomerB.run_sdftd3(func=dftd3_functional_name, dashlvl="d3mbj", verbose=True)
        data["D3 MONOMER B"] = monB_d3
        data["D3 IE"] = dimer_d3 - monA_d3 - monB_d3
    elif d3_type == "gd3_supermolecular":
        # This uses the default Grimme -D3 parameters for the given
        # functional with BJ damping and ATM three-body terms. Only use
        # this option in compbination with another baseline form of
        # dispersion like the delta_DFT correction:
        # "SAPT_DFT_DO_DDFT = True"
        core.print_out(
            "         "
            + "Supermolecular GD3(BJ)+ATM Interaction Energy E_IE = E_IJ - E_I - E_J".center(
                58
            )
            + "\n"
        )
        dimer_d3, _ = sapt_dimer.run_sdftd3(
            func=dftd3_functional_name, dashlvl="d3bjeeqatm"
        )
        data["D3 DIMER"] = dimer_d3
        monA_d3, _ = monomerA.run_sdftd3(func=dftd3_functional_name, dashlvl="d3bj2b")
        data["D3 MONOMER A"] = monA_d3
        monB_d3, _ = monomerB.run_sdftd3(func=dftd3_functional_name, dashlvl="d3bj2b")
        data["D3 MONOMER B"] = monB_d3
        data["D3 IE"] = dimer_d3 - monA_d3 - monB_d3
    elif d3_type == "intermolecular":
        print(dftd3_functional_name)
        dimer_d3, _ = sapt_dimer.run_sdftd3(dftd3_functional_name, dashlvl="d3bj2b", property=True)
        geom, _, _, elez, _ = sapt_dimer.to_arrays()
        elez = np.array(elez, dtype=np.int32)
        monAs = np.array([i for i in range(monomerA.natom()) if monomerA.Z(i) > 0])
        monBs = np.array([i for i in range(monomerB.natom()) if monomerB.Z(i) > 0])
        E_disp = 0.0
        FSAPT_EMPIRICAL_DISP = np.zeros_like(core.variable("DFTD3 ADDITIVE PAIRWISE ENERGY").np)
        for A in monAs:
            for B in monBs:
                E_disp += core.variable("DFTD3 ADDITIVE PAIRWISE ENERGY").np[A, B]
                FSAPT_EMPIRICAL_DISP[A, B] = core.variable("DFTD3 ADDITIVE PAIRWISE ENERGY").np[A, B]
        data["D3 IE"] = E_disp
        data["FSAPT_EMPIRICAL_DISP"] = FSAPT_EMPIRICAL_DISP
    else:
        raise ValueError(
            "SAPT(DFT): d3_type must be 'supermolecular' or 'intermolecular'."
        )
    return data
