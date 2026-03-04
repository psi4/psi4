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

# Need to import FISAPT to set fdrop, plot, save_fsapt_variables methods
from . import fisapt_proc
from pprint import pprint as pp


def to_matrix(obj):
    """Convert object to psi4.core.Matrix, handling both numpy arrays and Matrix objects."""
    if isinstance(obj, core.Matrix):
        return obj.clone()
    else:
        return core.Matrix.from_array(obj)


def to_vector(obj):
    """Convert object to psi4.core.Vector, handling both numpy arrays and Vector objects."""
    if isinstance(obj, core.Vector):
        return obj.clone()
    else:
        return core.Vector.from_array(obj)


def setup_fisapt_object(
    wfn, wfn_A, wfn_B, cache, scalars, basis_set=None, do_flocalize=False
):
    """
    Setup FISAPT object for F-SAPT calculations.

    Parameters
    ----------
    wfn : psi4.core.Wavefunction
        Dimer wavefunction
    wfn_A : psi4.core.Wavefunction
        Monomer A wavefunction
    wfn_B : psi4.core.Wavefunction
        Monomer B wavefunction
    cache : dict
        SAPT(DFT) cache containing orbital data
    scalars : dict
        SAPT energy components dictionary
    basis_set : psi4.core.BasisSet, optional
        Auxiliary basis set for DF. If None, uses primary basis.
    do_flocalize : bool, optional
        If True, call flocalize() to localize monomer orbitals using C++ IBOLocalizer2.
        The localized orbitals will be stored in the FISAPT object's internal matrices.
        Default is False.

    Returns
    -------
    psi4.core.FISAPT
        FISAPT object ready for felst(), fexch(), find()
    """
    # Setup FISAPT object
    if basis_set is None:
        basis_set = wfn.basisset()
    wfn.set_basisset("DF_BASIS_SAPT", basis_set)

    # Set MINAO basis for IBOLocalizer2 (needed if do_flocalize)
    if do_flocalize:
        minao = core.BasisSet.build(
            wfn.molecule(), "BASIS", core.get_global_option("MINAO_BASIS")
        )
        wfn.set_basisset("MINAO", minao)

    fisapt = core.FISAPT(wfn)

    # Used to slice arrays later if frozen core is requested
    nfrozen_A = wfn_A.basisset().n_frozen_core(
        core.get_global_option("FREEZE_CORE"), wfn_A.molecule()
    )
    nfrozen_B = wfn_B.basisset().n_frozen_core(
        core.get_global_option("FREEZE_CORE"), wfn_B.molecule()
    )

    # Basic matrix keys always needed from cache
    basic_matrix_keys = {
        "Cocc_A": "Cocc0A",
        "Cvir_A": "Cvir0A",
        "Cocc_B": "Cocc0B",
        "Cvir_B": "Cvir0B",
    }

    matrix_cache = {
        fisapt_key: to_matrix(cache[sdft_key])
        for sdft_key, fisapt_key in basic_matrix_keys.items()
    }

    other_keys = [
        "S",
        "D_A",
        "J_A",
        "K_A",
        "V_A",
        "D_B",
        "J_B",
        "K_B",
        "V_B",
        "K_O",
        "J_O",
        # "J_P_A",
        # "J_P_B",
        # fdisp
        "P_A",
        "P_B",
    ]
    for key in other_keys:
        matrix_cache[key] = to_matrix(cache[key])
    # J_P_A and J_P_B have flipped terminology in FISAPT...
    matrix_cache["J_P_A"] = to_matrix(cache["J_P_B"])
    matrix_cache["J_P_B"] = to_matrix(cache["J_P_A"])

    # Vector keys for eigenvalues
    vector_keys = {
        "eps_occ_A": "eps_occ0A",
        "eps_vir_A": "eps_vir0A",
        "eps_occ_B": "eps_occ0B",
        "eps_vir_B": "eps_vir0B",
    }
    vector_cache = {
        fisapt_key: to_vector(cache[sdft_key])
        for sdft_key, fisapt_key in vector_keys.items()
    }
    other_vector_keys = [
        "ZA",
        "ZA_orig",
        "ZB",
        "ZB_orig",
        "ZC",
        "ZC_orig",
    ]
    # When not using einsums, we are using localization from the FISAPT object
    if "ZA" not in cache:
        mol = wfn.molecule()
        natoms = mol.natom()
        fragments = mol.get_fragments()
        indA = np.arange(*fragments[0], dtype=int)
        indB = np.arange(*fragments[1], dtype=int)
        indC = np.arange(*fragments[2], dtype=int) if len(fragments) == 3 else np.array([], dtype=int)
        # --- Z per fragment originals ---
        ZA = core.Vector(natoms)
        ZB = core.Vector(natoms)
        ZC = core.Vector(natoms)
        ZA.np[:] = 0.0
        ZB.np[:] = 0.0
        ZC.np[:] = 0.0

        Z_all = np.array([mol.Z(i) for i in range(natoms)], dtype=float)
        ZA.np[indA] = Z_all[indA]
        ZB.np[indB] = Z_all[indB]
        if indC.size:
            ZC.np[indC] = Z_all[indC]
        cache["ZA"] = ZA
        cache["ZB"] = ZB
        cache["ZC"] = ZC
        cache["ZA_orig"] = core.Vector.from_array(ZA.np.copy())
        cache["ZB_orig"] = core.Vector.from_array(ZB.np.copy())
        cache["ZC_orig"] = core.Vector.from_array(ZC.np.copy())
    for key in other_vector_keys:
        vector_cache[key] = to_vector(cache[key])

    # Set up frozen/active occupied orbitals for flocalize()
    Cocc_A_np = np.asarray(matrix_cache["Cocc0A"])
    Cocc_B_np = np.asarray(matrix_cache["Cocc0B"])

    # Monomer A: frozen and active
    if nfrozen_A > 0:
        matrix_cache["Cfocc0A"] = core.Matrix.from_array(Cocc_A_np[:, :nfrozen_A])
        matrix_cache["Caocc0A"] = core.Matrix.from_array(Cocc_A_np[:, nfrozen_A:])
    else:
        matrix_cache["Cfocc0A"] = core.Matrix.from_array(
            np.zeros((Cocc_A_np.shape[0], 0))
        )
        matrix_cache["Caocc0A"] = core.Matrix.from_array(Cocc_A_np.copy())

    # Monomer B: frozen and active
    if nfrozen_B > 0:
        matrix_cache["Cfocc0B"] = core.Matrix.from_array(Cocc_B_np[:, :nfrozen_B])
        matrix_cache["Caocc0B"] = core.Matrix.from_array(Cocc_B_np[:, nfrozen_B:])
    else:
        matrix_cache["Cfocc0B"] = core.Matrix.from_array(
            np.zeros((Cocc_B_np.shape[0], 0))
        )
        matrix_cache["Caocc0B"] = core.Matrix.from_array(Cocc_B_np.copy())

    # Set initial matrices/vectors before flocalize
    fisapt.set_matrix(matrix_cache)
    fisapt.set_vector(vector_cache)

    if do_flocalize:
        # Call C++ flocalize() to localize monomer orbitals
        # This populates the FISAPT object's internal matrices with:
        # Locc0A, Locc0B, Uocc0A, Uocc0B, Qocc0A, Qocc0B,
        # Lfocc0A, Lfocc0B, Laocc0A, Laocc0B, Ufocc0A, Ufocc0B, Uaocc0A, Uaocc0B
        fisapt.flocalize()

        # Get the localized matrices back from FISAPT object
        matrices = fisapt.matrices()

        # Update matrix_cache with localized orbitals for subsequent processing
        matrix_cache["Locc0A"] = matrices["Locc0A"]
        matrix_cache["Locc0B"] = matrices["Locc0B"]
        matrix_cache["Uocc0A"] = matrices["Uocc0A"]
        matrix_cache["Uocc0B"] = matrices["Uocc0B"]
        matrix_cache["Qocc0A"] = matrices["Qocc0A"]
        matrix_cache["Qocc0B"] = matrices["Qocc0B"]
        matrix_cache["Laocc0A"] = matrices["Laocc0A"]
        matrix_cache["Laocc0B"] = matrices["Laocc0B"]
        matrix_cache["Lfocc0A"] = matrices["Lfocc0A"]
        matrix_cache["Lfocc0B"] = matrices["Lfocc0B"]
        matrix_cache["Uaocc0A"] = matrices["Uaocc0A"]
        matrix_cache["Uaocc0B"] = matrices["Uaocc0B"]
        matrix_cache["Ufocc0A"] = matrices["Ufocc0A"]
        matrix_cache["Ufocc0B"] = matrices["Ufocc0B"]
    else:
        # Localized orbitals already in cache - get them from there
        localized_matrix_keys = {
            "Locc_A": "Locc0A",
            "Locc_B": "Locc0B",
            "Uocc_A": "Uocc0A",
            "Uocc_B": "Uocc0B",
            "Qocc0A": "Qocc0A",
            "Qocc0B": "Qocc0B",
            "Laocc0A": "Laocc0A",
            "Laocc0B": "Laocc0B",
            "Lfocc0A": "Lfocc0A",
            "Lfocc0B": "Lfocc0B",
            "Uaocc0A": "Uaocc0A",
            "Uaocc0B": "Uaocc0B",
        }
        for sdft_key, fisapt_key in localized_matrix_keys.items():
            matrix_cache[fisapt_key] = to_matrix(cache[sdft_key])

        # Update Uaocc matrices based on frozen core slicing
        if nfrozen_A > 0:
            matrix_cache["Uaocc0A"] = core.Matrix.from_array(
                np.asarray(matrix_cache["Uocc0A"])[:, nfrozen_A:]
            )
        else:
            matrix_cache["Uaocc0A"] = core.Matrix.from_array(
                np.asarray(matrix_cache["Uocc0A"]).copy()
            )
        if nfrozen_B > 0:
            matrix_cache["Uaocc0B"] = core.Matrix.from_array(
                np.asarray(matrix_cache["Uocc0B"])[:, nfrozen_B:]
            )
        else:
            matrix_cache["Uaocc0B"] = core.Matrix.from_array(
                np.asarray(matrix_cache["Uocc0B"]).copy()
            )

    # Set eps_aocc vectors
    if nfrozen_A > 0:
        vector_cache["eps_aocc0A"] = core.Vector.from_array(
            np.asarray(vector_cache["eps_occ0A"])[nfrozen_A:]
        )
    else:
        vector_cache["eps_aocc0A"] = core.Vector.from_array(
            np.asarray(vector_cache["eps_occ0A"]).copy()
        )
    if nfrozen_B > 0:
        vector_cache["eps_aocc0B"] = core.Vector.from_array(
            np.asarray(vector_cache["eps_occ0B"])[nfrozen_B:]
        )
    else:
        vector_cache["eps_aocc0B"] = core.Vector.from_array(
            np.asarray(vector_cache["eps_occ0B"]).copy()
        )

    # Set all matrices and vectors on FISAPT object
    fisapt.set_matrix(matrix_cache)
    fisapt.set_vector(vector_cache)

    scalar_keys = {
        "Ind20,r (A<-B)": "Ind20,r (A<-B)",
        "Ind20,r (A->B)": "Ind20,r (B<-A)",
        "Ind20,u (A<-B)": "Ind20,u (A<-B)",
        "Ind20,u (A->B)": "Ind20,u (B<-A)",
        "Exch-Ind20,r (A<-B)": "Exch-Ind20,r (A<-B)",
        "Exch-Ind20,r (A->B)": "Exch-Ind20,r (B<-A)",
        "Exch-Ind20,u (A<-B)": "Exch-Ind20,u (A<-B)",
        "Exch-Ind20,u (A->B)": "Exch-Ind20,u (B<-A)",
        "Exch10": "Exch10",
        "Exch10(S^2)": "Exch10(S^2)",
        "Elst10,r": "Elst10,r",
        "Ind20,r": "Ind20,r",
        "Exch-Ind20,r": "Exch-Ind20,r",
    }
    # "DHF VALUE": "HF",
    if core.get_option("SAPT", "SAPT_DFT_DO_DHF"):
        scalar_keys["DHF VALUE"] = "HF"
    scalar_cache = {
        fisapt_key: scalars[sdft_key] for sdft_key, fisapt_key in scalar_keys.items()
    }
    fisapt.set_scalar(scalar_cache)
    return fisapt


def drop_saptdft_variables(wfn, wfn_A, wfn_B, cache, scalars):
    """
    Setup FISAPT object to call fisapt_fdrop for dropping SAPT(DFT) variables.

    Parameters
    ----------
    wfn : psi4.core.Wavefunction
        Dimer wavefunction
    wfn_A : psi4.core.Wavefunction
        Monomer A wavefunction
    wfn_B : psi4.core.Wavefunction
        Monomer B wavefunction
    cache : dict
        SAPT(DFT) cache containing orbital data
    scalars : dict
        SAPT energy components dictionary
    """
    fisapt = core.FISAPT(wfn)
    # iterate through cache and scalars to set these labels for fisapt_fdrop:
    cache_keys = {
        "Qocc0A": "Qocc0A",
        "Qocc0B": "Qocc0B",
        "Elst_AB": "Elst_AB",
        "Exch_AB": "Exch_AB",
        "IndAB_AB": "IndAB_AB",
        "IndBA_AB": "IndBA_AB",
    }
    # Set whether to drop dispersion matrix... fisapt has specific option for
    # this...
    core.set_local_option("FISAPT", "FISAPT_DO_FSAPT_DISP", core.get_option("SAPT", "SAPT_DFT_DO_DISP"))
    if core.get_option("SAPT", "SAPT_DFT_DO_DISP"):
        cache_keys["Disp_AB"] = "Disp_AB"
    matrix_cache = {
        fisapt_key: cache[sdft_key] for sdft_key, fisapt_key in cache_keys.items()
    }
    vector_cache = {
        "ZA": to_vector(cache["ZA"]),
        "ZA_orig": to_vector(cache["ZA_orig"]),
        "ZB": to_vector(cache["ZB"]),
        "ZB_orig": to_vector(cache["ZB_orig"]),
        "ZC": to_vector(cache["ZC"]),
        "ZC_orig": to_vector(cache["ZC_orig"]),
    }
    for key in vector_cache.keys():
        vector_cache[key].name = key.upper()
    fisapt.set_matrix(matrix_cache)
    fisapt.set_vector(vector_cache)
    fisapt.fdrop()
    fisapt.save_variables_to_wfn(wfn, sapt_type='SAPT(DFT)')
    # Now drop empirical dispersion if computed
    if core.get_option("SAPT", "SAPT_DFT_D4_IE"):
        pw_disp = cache["FSAPT_EMPIRICAL_DISP"]
        pw_disp.name = "Empirical_Disp"
        filepath = core.get_option("FISAPT", "FISAPT_FSAPT_FILEPATH")
        if filepath.lower() != "none":
            fisapt_proc._drop(pw_disp, filepath)
        core.set_variable("FSAPT_" + pw_disp.name.upper(), pw_disp)
    return wfn
