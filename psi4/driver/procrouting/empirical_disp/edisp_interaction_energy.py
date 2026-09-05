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

from . import empirical_dispersion


def _sapt_dft_dispersion_interaction_energy(
    sapt_dimer: core.Molecule,
    monomerA: core.Molecule,
    monomerB: core.Molecule,
    dimer_wfn: core.Wavefunction,
    functional_name: str,
    d_type: str,
    data: dict[str, object],
    *,
    disp_label: str,
    gd_type: str,
    supermolecular_text: str,
    gd_supermolecular_text: str,
    two_body_level_hint: str,
    atm_level_hint: str,
) -> dict[str, object]:
    """Compute SAPT(DFT) empirical dispersion interaction terms.

    Parameters
    ----------
    sapt_dimer : core.Molecule
        Dimer molecule used for the SAPT(DFT) dispersion evaluation.
    monomerA : core.Molecule
        Monomer A molecule.
    monomerB : core.Molecule
        Monomer B molecule.
    dimer_wfn : core.Wavefunction
        Dimer wavefunction object used to store/read pairwise variables.
    functional_name : str
        Dispersion parameter lookup key (usually functional name).
    d_type : str
        Dispersion mode selector (e.g., ``supermolecular``).
    data : dict[str, object]
        Result dictionary updated in place.
    disp_label : str
        Label prefix for output keys (e.g., ``D3`` or ``D4``).
    gd_type : str
        Selector value for Grimme ATM-enabled supermolecular mode.
    supermolecular_text : str
        Header text for the regular supermolecular branch.
    gd_supermolecular_text : str
        Header text for the ATM-enabled supermolecular branch.
    two_body_level_hint : str
        Dispersion level hint for two-body-only calculations.
    atm_level_hint : str
        Dispersion level hint for ATM-enabled calculations.

    Returns
    -------
    dict[str, object]
        Updated ``data`` dictionary containing dispersion terms.
    """

    if core.has_option_changed("SCF", "DFT_DISPERSION_PARAMETERS"):
        modified_disp_params = core.get_option("SCF", "DFT_DISPERSION_PARAMETERS")
    else:
        modified_disp_params = None

    if d_type == "supermolecular":
        core.print_out(
            "         "
            + supermolecular_text.center(58)
            + "\n\n"
        )
        _disp_functor = empirical_dispersion.EmpiricalDispersion(
            name_hint=functional_name,
            level_hint=two_body_level_hint,
            param_tweaks=modified_disp_params,
            save_pairwise_disp=True,
        )
        dimer_disp = _disp_functor.compute_energy(sapt_dimer, dimer_wfn)
        monA_disp = _disp_functor.compute_energy(monomerA)
        monB_disp = _disp_functor.compute_energy(monomerB)
        data[f"{disp_label} DIMER"] = dimer_disp
        data[f"{disp_label} MONOMER A"] = monA_disp
        data[f"{disp_label} MONOMER B"] = monB_disp
        data[f"{disp_label} IE"] = dimer_disp - monA_disp - monB_disp
        data['FSAPT_EMPIRICAL_DISP'] = dimer_wfn.variable("PAIRWISE DISPERSION CORRECTION ANALYSIS")
        _disp_functor.print_out()
    elif d_type == "intermolecular":
        monAs = np.array([i for i in range(monomerA.natom()) if monomerA.Z(i) > 0])
        monBs = np.array([i for i in range(monomerB.natom()) if monomerB.Z(i) > 0])
        _disp_functor = empirical_dispersion.EmpiricalDispersion(
            name_hint=functional_name,
            level_hint=two_body_level_hint,
            param_tweaks=modified_disp_params,
            save_pairwise_disp=True,
        )
        E_disp = 0.0
        _ = _disp_functor.compute_energy(sapt_dimer, dimer_wfn)
        pairwise_energies = dimer_wfn.variable("PAIRWISE DISPERSION CORRECTION ANALYSIS").np
        FSAPT_EMPIRICAL_DISP = np.zeros_like(pairwise_energies)
        for A in monAs:
            for B in monBs:
                # account for A-B and B-A with factor of 2
                E_disp += 2 * pairwise_energies[A, B]
                FSAPT_EMPIRICAL_DISP[A, B] = pairwise_energies[A, B]
                FSAPT_EMPIRICAL_DISP[B, A] = pairwise_energies[A, B]
        data[f"{disp_label} IE"] = E_disp
        data['FSAPT_EMPIRICAL_DISP'] = FSAPT_EMPIRICAL_DISP
        _disp_functor.print_out()
    elif d_type == gd_type:
        # This uses the default Grimme parameters for the given
        # functional with BJ damping and ATM three-body terms. Only use
        # this option in compbination with another baseline form of
        # dispersion like the delta_DFT correction:
        # "SAPT_DFT_DO_DDFT = True"
        core.print_out(
            "         "
            + gd_supermolecular_text.center(58)
            + "\n\n"
        )
        _disp_functor = empirical_dispersion.EmpiricalDispersion(
            name_hint=functional_name,
            level_hint=atm_level_hint,
            param_tweaks=modified_disp_params,
            save_pairwise_disp=True,
        )

        dimer_disp = _disp_functor.compute_energy(sapt_dimer, dimer_wfn)
        monA_disp = _disp_functor.compute_energy(monomerA)
        monB_disp = _disp_functor.compute_energy(monomerB)
        data[f"{disp_label} DIMER"] = dimer_disp
        data[f"{disp_label} MONOMER A"] = monA_disp
        data[f"{disp_label} MONOMER B"] = monB_disp
        data[f"{disp_label} IE"] = dimer_disp - monA_disp - monB_disp
        data['FSAPT_EMPIRICAL_DISP'] = dimer_wfn.variable("PAIRWISE DISPERSION CORRECTION ANALYSIS")
        _disp_functor.print_out()
    else:
        raise ValueError(
            "SAPT(DFT): d_type must be 'supermolecular' or 'intermolecular'."
        )
    return data


def sapt_dft_d4_interaction_energy(
    sapt_dimer: core.Molecule,
    monomerA: core.Molecule,
    monomerB: core.Molecule,
    dimer_wfn: core.Wavefunction,
    dftd4_functional_name: str,
    d4_type: str,
    data: dict[str, object],
) -> dict[str, object]:
    """Compute SAPT(DFT) D4 interaction energy contributions.

    Parameters
    ----------
    sapt_dimer : core.Molecule
        Dimer molecule used for the SAPT(DFT) dispersion evaluation.
    monomerA : core.Molecule
        Monomer A molecule.
    monomerB : core.Molecule
        Monomer B molecule.
    dimer_wfn : core.Wavefunction
        Dimer wavefunction object used to store/read pairwise variables.
    dftd4_functional_name : str
        Dispersion parameter lookup key for D4.
    d4_type : str
        D4 mode selector (e.g., ``supermolecular``, ``intermolecular``, or ``gd4_supermolecular``).
    data : dict[str, object]
        Result dictionary updated in place.

    Returns
    -------
    dict[str, object]
        Updated ``data`` dictionary containing D4 terms.
    """

    return _sapt_dft_dispersion_interaction_energy(
        sapt_dimer,
        monomerA,
        monomerB,
        dimer_wfn,
        dftd4_functional_name,
        d4_type,
        data,
        disp_label="D4",
        gd_type="gd4_supermolecular",
        supermolecular_text="Supermolecular D4 Interaction Energy E_IE = E_IJ - E_I - E_J",
        gd_supermolecular_text="Supermolecular GD4(BJ)+ATM Interaction Energy E_IE = E_IJ - E_I - E_J",
        two_body_level_hint='d4bj2b',
        atm_level_hint='d4bjeeqatm',
    )


def sapt_dft_d3_interaction_energy(
    sapt_dimer: core.Molecule,
    monomerA: core.Molecule,
    monomerB: core.Molecule,
    dimer_wfn: core.Wavefunction,
    dftd3_functional_name: str,
    d3_type: str,
    data: dict[str, object],
) -> dict[str, object]:
    """Compute SAPT(DFT) D3 interaction energy contributions.

    Parameters
    ----------
    sapt_dimer : core.Molecule
        Dimer molecule used for the SAPT(DFT) dispersion evaluation.
    monomerA : core.Molecule
        Monomer A molecule.
    monomerB : core.Molecule
        Monomer B molecule.
    dimer_wfn : core.Wavefunction
        Dimer wavefunction object used to store/read pairwise variables.
    dftd3_functional_name : str
        Dispersion parameter lookup key for D3.
    d3_type : str
        D3 mode selector (e.g., ``supermolecular``, ``intermolecular``, or ``gd3_supermolecular``).
    data : dict[str, object]
        Result dictionary updated in place.

    Returns
    -------
    dict[str, object]
        Updated ``data`` dictionary containing D3 terms.
    """

    return _sapt_dft_dispersion_interaction_energy(
        sapt_dimer,
        monomerA,
        monomerB,
        dimer_wfn,
        dftd3_functional_name,
        d3_type,
        data,
        disp_label="D3",
        gd_type="gd3_supermolecular",
        supermolecular_text="Supermolecular -D3MBJ Interaction Energy E_IE = E_IJ - E_I - E_J",
        gd_supermolecular_text="Supermolecular GD3(BJ)+ATM Interaction Energy E_IE = E_IJ - E_I - E_J",
        two_body_level_hint='d3mbj2b',
        atm_level_hint='d3mbjatm',
    )
