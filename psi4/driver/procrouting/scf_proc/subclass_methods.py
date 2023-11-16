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
"""
Extensions to :class:`psi4.core.RHF`, :class:`psi4.core.UHF`,
:class:`psi4.core.CUHF`, and :class:`psi4.core.ROHF` for ``diis`` and
``compute_orbital_gradient`` methods.

"""

__all__ = []

import math
from typing import Set

import numpy as np

from psi4 import core

from ...p4util import solvers
from ..diis import DIIS, RemovalPolicy, StoragePolicy
from ..response.scf_products import TDUSCFEngine


def diis_engine_helper(self):
    engines = set()
    if core.get_option('SCF', 'DIIS'):
        engines.add('diis')
    restricted_open = self.same_a_b_orbs() and not self.same_a_b_dens()
    if not restricted_open:
        aediis = core.get_option('SCF', 'SCF_INITIAL_ACCELERATOR')
        if aediis != "NONE":
            engines.add(aediis.lower())
    return engines

def _RHF_orbital_gradient(self, save_fock: bool, max_diis_vectors: int) -> float:
    """Form :math:`| FDS - SDF |` quantity and return either its
    root-mean-square when :term:`DIIS_RMS_ERROR <DIIS_RMS_ERROR (SCF)>` is True
    or its absolute maximum element otherwise. Used for second-order SCF.

    Parameters
    ----------
    save_fock
        Whether to include step in DIIS.
    max_diis_vectors
        When `save_fock` is True and a new DIIS object is needed, the maximum
        number of vectors to initialize it for.

    """
    gradient = self.form_FDSmSDF(self.Fa(), self.Da())

    if save_fock:
        if not self.initialized_diis_manager_:
            storage_policy = StoragePolicy.InCore if self.scf_type() == "DIRECT" else StoragePolicy.OnDisk
            self.diis_manager_ = DIIS(max_diis_vectors, "HF DIIS vector", RemovalPolicy.LargestError, storage_policy, engines=diis_engine_helper(self))
            self.initialized_diis_manager_ = True

        entry = {"target": [self.Fa()]}
        if core.get_option('SCF', 'DIIS'):
            entry["error"] = [gradient]
        aediis = core.get_option('SCF', 'SCF_INITIAL_ACCELERATOR')
        if aediis != "NONE":
            entry["densities"] = [self.Da()]
            if aediis == "EDIIS":
                entry["energy"] = [self.compute_E()]
        self.diis_manager_.add_entry(entry)

    if self.options().get_bool("DIIS_RMS_ERROR"):
        return gradient.rms()
    else:
        return gradient.absmax()

def _UHF_orbital_gradient(self, save_fock: bool, max_diis_vectors: int) -> float:
    """Form :math:`| FDS - SDF |` quantity and return either its
    root-mean-square when :term:`DIIS_RMS_ERROR <DIIS_RMS_ERROR (SCF)>` is True
    or its absolute maximum element otherwise. Used for second-order SCF.

    Parameters
    ----------
    save_fock
        Whether to include step in DIIS.
    max_diis_vectors
        When `save_fock` is True and a new DIIS object is needed, the maximum
        number of vectors to initialize it for.

    """
    gradient_a = self.form_FDSmSDF(self.Fa(), self.Da())
    gradient_b = self.form_FDSmSDF(self.Fb(), self.Db())

    if save_fock:
        if not self.initialized_diis_manager_:
            self.diis_manager_ = DIIS(max_diis_vectors, "HF DIIS vector", RemovalPolicy.LargestError,
                                                          StoragePolicy.OnDisk, False, engines=diis_engine_helper(self))
            self.initialized_diis_manager_ = True

        entry = {"target": [self.Fa(), self.Fb()]}
        if core.get_option('SCF', 'DIIS'):
            entry["error"] = [gradient_a, gradient_b]
        aediis = core.get_option('SCF', 'SCF_INITIAL_ACCELERATOR')
        if aediis != "NONE":
            entry["densities"] = [self.Da(), self.Db()]
            if aediis == "EDIIS":
                entry["energy"] = [self.compute_E()]
        self.diis_manager_.add_entry(entry)

    if self.options().get_bool("DIIS_RMS_ERROR"):
        return math.sqrt(0.5 * (gradient_a.rms() ** 2 + gradient_b.rms() ** 2))
    else:
        return max(gradient_a.absmax(), gradient_b.absmax())

def _ROHF_orbital_gradient(self, save_fock: bool, max_diis_vectors: int) -> float:
    # Only the inact-act, inact-vir, and act-vir rotations are non-redundant
    dim_zero = core.Dimension(self.nirrep(), "Zero Dim")
    noccpi = self.doccpi() + self.soccpi()
    row_slice = core.Slice(dim_zero, noccpi)
    col_slice = core.Slice(self.doccpi(), self.nmopi())
    MOgradient = self.moFeff().get_block(row_slice, col_slice)

    # Zero the active-active block
    for h in range(MOgradient.nirrep()):
        socc = self.soccpi()[h]
        docc = self.doccpi()[h]

        MOgradient.nph[h][docc:docc+socc, 0:socc] = 0

    # Grab inact-act and act-vir orbs
    # Ct is (nmo x nmo), not the (nso x nmo) you would expect
    row_slice = core.Slice(dim_zero, self.nmopi())
    col_slice = core.Slice(dim_zero, noccpi)
    Cia = self.Ct().get_block(row_slice, col_slice)
    col_slice = core.Slice(self.doccpi(), self.nmopi())
    Cav = self.Ct().get_block(row_slice, col_slice)

    # Back transform MOgradient
    gradient = core.triplet(Cia, MOgradient, Cav, False, False, True)

    if save_fock:
        if not self.initialized_diis_manager_:
            self.diis_manager_ = DIIS(max_diis_vectors, "HF DIIS vector", RemovalPolicy.LargestError, StoragePolicy.OnDisk, engines=diis_engine_helper(self))
            self.diis_manager_.set_error_vector_size(gradient)
            self.diis_manager_.set_vector_size(self.soFeff())
            self.initialized_diis_manager_ = True

        self.diis_manager_.add_entry({"error": [gradient], "target": [self.soFeff()]})

    if self.options().get_bool("DIIS_RMS_ERROR"):
        return gradient.rms()
    else:
        return gradient.absmax()

core.RHF.compute_orbital_gradient = _RHF_orbital_gradient
core.UHF.compute_orbital_gradient = core.CUHF.compute_orbital_gradient = _UHF_orbital_gradient
core.ROHF.compute_orbital_gradient = _ROHF_orbital_gradient

_diis_docstring = """Perform DIIS extrapolation.

    Parameters
    ----------
    Dnorm
        Error metric used to blend certain algorithms like ADIIS and EDIIS.

    Returns
    -------
    ~typing.Set[str]
        All DIIS algorithms performed.

    """

def _RHF_diis(self, Dnorm: float) -> Set[str]:
    return self.diis_manager_.extrapolate(self.Fa(), Dnorm=Dnorm)

def _UHF_diis(self, Dnorm: float) -> Set[str]:
    return self.diis_manager_.extrapolate(self.Fa(), self.Fb(), Dnorm=Dnorm)

def _ROHF_diis(self, Dnorm: float) -> Set[str]:
    return self.diis_manager_.extrapolate(self.soFeff(), Dnorm=Dnorm)

_RHF_diis.__doc__ = _diis_docstring
_UHF_diis.__doc__ = _diis_docstring
_ROHF_diis.__doc__ = _diis_docstring

core.RHF.diis = _RHF_diis
core.UHF.diis = core.CUHF.diis = _UHF_diis
core.ROHF.diis = _ROHF_diis

def _UHF_stability_analysis(self):
    # => Validate options <=
    # TODO: Stability analysis is supported for any functional UKS functional where its one-
    # and two-body Hamiltonian matrix-vector products are implemented. This is true for LDA
    # at least, but probably not other functionals. This restriction exists for now because
    # we've always had it, but we should lift it as much as we can and implement more matrix-
    # vector products so we can lift it further.
    # TODO: It should be up to the SolverEngine to validate whether it can do Hx products for the input wfn.
    if self.functional().is_meta() or self.functional().needs_vv10():
        raise ValidationError("Stability Analysis: Unrestricted Kohn-Sham Vx kernel does not support meta or VV10 functionals.")

    # => Prep options for eigenvector solver <=
    if not core.has_option_changed("SCF", "SOLVER_ROOTS_PER_IRREP"):
        roots = [core.get_option("SCF", "SOLVER_N_ROOT")] * self.nirrep()
    else:
        roots = core.get_option("SCF", "SOLVER_ROOTS_PER_IRREP")
        if len(roots) != wfn.nirrep():
            raise ValidationError(f"SOLVER_ROOTS_PER_IRREP specified {wfn.nirrep()} irreps, but there are {len(roots)} irreps.")
    r_convergence = core.get_option("SCF", "SOLVER_CONVERGENCE")
    # Below formula borrowed from TDSCF code.
    max_vecs_per_root = int(-np.log10(r_convergence) * 50)
    engine = TDUSCFEngine(self, ptype="hess")

    # => Compute eigenvectors and do trivial data processing <=
    eval_sym = core.Matrix("SCF STABILITY EIGENVALUES", core.Dimension(roots), core.Dimension([1] * self.nirrep()))
    current_eigenvalue = None
    unstable = False
    for h, nroot in enumerate(roots):
        if not nroot: continue
        if not core.has_option_changed("SCF", "SOLVER_N_GUESS"):
            nguess = nroot * 4
        else:
            nguess = core.get_option("SCF", "SOLVER_N_GUESS")
        # The below line changes the guess the engine generates, which controls the final states.
        # This selects for eigenvectors of irrep h.
        engine.reset_for_state_symm(engine.G_gs ^ h)
        ret = solvers.davidson_solver(engine=engine,
                                      nroot=nroot,
                                      guess=engine.generate_guess(nguess),
                                      r_convergence=r_convergence,
                                      max_ss_size=max_vecs_per_root * nroot,
                                      verbose=0)
        if not ret["stats"][-1]["done"]:
            raise SCFConvergenceError(maxiter, self, f"hessian eigenvectors in irrep {irrep_ES}", ret["stats"][-1])
        if h == 0:
            current_eigenvalue = ret["eigvals"][0]
            # Distinction between left and right eigenvectors is a formality for TDA-type solvers but forces the extra [0].
            current_eigenvector = ret["eigvecs"][0][0]
            unstable = current_eigenvalue < 0
        for i, eigval in enumerate(ret["eigvals"]):
            eval_sym.set(h, i, 0, eigval)

    # => Print out whether unstable <=
    if unstable:
        core.print_out(f"    Negative totally symmetric eigenvalue detected: {current_eigenvalue:.6f} \n")
        core.print_out(f"    Wavefunction unstable!\n")
    else:
        core.print_out(f"    Wavefunction stable under totally symmetric rotations.\n")
        core.print_out(f"    Lowest totally symmetric eigenvalue: {current_eigenvalue:.6f} \n")

    # => Print out and save stability eigenvalues <=
    core.print_out("    Lowest UHF->UHF stability eigenvalues: \n");
    eval_sym_pairs = []
    for h in range(eval_sym.nirrep()):
        for i in range(eval_sym.rows(h)):
            eval_sym_pairs.append((eval_sym.get(h, i, 0), h))
    self.print_stability_analysis(eval_sym_pairs)
    self.set_variable("SCF STABILITY EIGENVALUES", eval_sym)


    # => Follow instability or print out that there's nothing left to do <=
    # Legacy instability took orbital steps based on the following algorithm:
    # * Normalize the orbital eigenvector X to 1
    # * Apply exp(t(X-X^)) for t = lambda pi / 2 (defaults to lambda = 0.5)
    # * If that step returns to the same minimum, increment lambda (defaulting to 0.2) and repeat previous step
    # Rigorous mathematical analysis on the true minimum is hard to come by: the rotated orbitals need not even be periodic in t.
    # (See DOI 10.1063/1.467504 eq. 8 for explicit formulas. You can show non-periodicity in general in the simple case that P^1/2 is 2-by-2 diagonal.)
    # As such, this algorithm is best regarded as a first attempt open to improvements.
    # Example improvement: if the orbital rotation increases the energy, take a smaller step, not a larger.
    if unstable and core.get_option("SCF", "STABILITY_ANALYSIS") == "FOLLOW":
        # ==> Increment step_scale_ if necessary <==
        if hasattr(self, "last_hess_eigval") and abs(self.last_hess_eigval - current_eigenvalue) < 1e-4:
            core.print_out("    Negative eigenvalue similar to previous one, wavefunction\n")
            core.print_out("    likely to be in the same minimum.\n")
            self.step_scale += core.get_option("SCF", "FOLLOW_STEP_INCREMENT")
            core.print_out(f"    Modifying FOLLOW_STEP_SCALE to {step_scale}.\n")
        else:
            self.step_scale = core.get_option("SCF", "FOLLOW_STEP_SCALE")
            self.last_hess_eigval = current_eigenvalue
        # ==> Perform the orbital rotation! <==
        # The current eigenvector is normalized to 1/2.
        core.print_out(f"    Rotating orbitals by {self.step_scale} * pi / 2 radians along unstable eigenvector.\n");
        current_eigenvector[0].scale(self.step_scale * np.pi)
        self.rotate_orbitals(self.Ca(), current_eigenvector[0])
        current_eigenvector[1].scale(self.step_scale * np.pi)
        self.rotate_orbitals(self.Cb(), current_eigenvector[1])
        return True
    else:
        core.print_out("    Stability analysis over.\n")
        return False


core.UHF.stability_analysis = _UHF_stability_analysis
