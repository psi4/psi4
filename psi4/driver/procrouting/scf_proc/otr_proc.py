
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2026 The Psi4 Developers.
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
Psi4's side of the interface to OpenTrustRegion, the external second-order SCF
optimizer. The bare names ``opentrustregion`` and ``otr`` belong to that package's own
Python interface, so this module does not take either.

It takes over the iteration loop once :py:func:`~psi4.driver.procrouting.scf_proc.
scf_iterator.scf_iterate` reaches the second-order handoff, runs to convergence or to
its own floating-point precision floor, and hands the orbitals back. This module holds
that run and everything around it: the occupation-restart loop, what to do at the
precision floor, and the first-order polish that may follow.
"""
import math

from psi4 import core

from ...p4util.exceptions import PsiException, SCFConvergenceError


#: OpenTrustRegion's error code for running out of macro-iterations which is subject to
#: FAIL_ON_MAXITER. Every other code is a real failure: an invalid setting, a numerical
#: breakdown inside the solver, or a callback into psi4 that raised (from 1100 up) whose
#: message opentrustregion_callback_error() then holds.
_OTR_MAX_ITER_ERROR = 102


#: Returned by :py:func:`run` when OpenTrustRegion stopped at its precision
#: floor short of both convergence criteria and first-order iterations should finish the SCF.
NEEDS_POLISH = object()


def shortfall(self, e_conv, d_conv):
    """How far OpenTrustRegion finished from this SCF's thresholds, for the output."""
    _, grad, dE, _ = self.otr_outcome_
    # A run of a single macro-iteration had no earlier point to measure the energy change
    # against, and reports NaN rather than a number for it.
    dE_text = "unmeasured" if math.isnan(dE) else f"{dE:.2e}"
    return (f"orbital gradient {grad:.2e} against D_CONVERGENCE {d_conv:.2e}, last energy change "
            f"{dE_text} against E_CONVERGENCE {e_conv:.2e}")


#: How many times to restart OpenTrustRegion at a new occupation when the aufbau rule
#: disagrees with the one it optimized. Diagonalizing preserves occupation, but the rule
#: applied to the resulting orbital energies may not agree with it. Most calculations
#: never restart: the first-order iterations settle the occupation before the handoff, so the
#: loop exits on its first pass. Of the exits that do follow a restart, the energy test
#: usually comes first -- an occupation that keeps moving without lowering the energy is a
#: relabeling of degenerate orbitals rather than a better solution -- and this count is the backstop.
_OTR_MAX_OCCUPATION_RESTARTS = 5


def _occupation(wfn):
    """Snapshot the alpha and beta occupations per irrep as plain tuples."""
    na, nb = wfn.nalphapi(), wfn.nbetapi()
    return (tuple(na[h] for h in range(na.n())), tuple(nb[h] for h in range(nb.n())))


def _orbital_state(wfn):
    """Snapshot what defines the orbitals: C, orbital energies and occupation, and for ROHF the
    orthogonalized coefficients Ct, which it keeps alongside Ca and builds its Fock matrix from."""
    state = {"Ca": wfn.Ca().clone(), "Cb": wfn.Cb().clone(),
             "epsilon_a": wfn.epsilon_a().clone(), "epsilon_b": wfn.epsilon_b().clone(),
             "occupation": _occupation(wfn)}
    if isinstance(wfn, core.ROHF):
        state["Ct"] = wfn.Ct().clone()
    return state


def _restore_orbital_state(wfn, state):
    """Put back an :py:func:`_orbital_state` snapshot. The caller rebuilds D, G and F from it."""
    na, nb = state["occupation"]
    wfn.force_occpi(core.Dimension(list(nb)), core.Dimension([a - b for a, b in zip(na, nb)]))
    wfn.Ca().copy(state["Ca"])
    wfn.Cb().copy(state["Cb"])
    wfn.epsilon_a().copy(state["epsilon_a"])
    wfn.epsilon_b().copy(state["epsilon_b"])
    if "Ct" in state:
        wfn.Ct().copy(state["Ct"])


def _aufbau_occupation(wfn):
    """What the next find_occupation would leave the occupation as, without assigning it.

    HF::next_occupation answers this where find_occupation's own aufbau branch lives, so the
    two cannot drift. A given DOCC/SOCC and MOM each keep the occupation already set; only the
    aufbau branch can move it.
    """
    na, nb = wfn.next_occupation()
    return (tuple(na[h] for h in range(na.n())), tuple(nb[h] for h in range(nb.n())))


def run(self, e_conv, d_conv):
    """Hand the orbitals to OpenTrustRegion and finish the SCF with them.

    Called at the second-order handoff point, once the first-order package has brought the
    orbital gradient inside |scf__soscf_start_convergence| at or past iteration
    |scf__soscf_start|. Leaves the wavefunction canonical and the iteration bookkeeping
    consistent with the internal optimizer.
    """
    # OpenTrustRegion optimizes at fixed occupation, where the internal solver re-runs
    # find_occupation every iteration. So the aufbau rule, applied to what OpenTrustRegion
    # converged, can disagree with the occupation it optimized at, and this loop reconverges
    # at the rule's answer. The default |scf__soscf_start| usually leaves nothing to move, but
    # earlier soscf can reach the branch.
    SCFE_prev_attempt = None
    prev_attempt = None
    for attempt in range(_OTR_MAX_OCCUPATION_RESTARTS + 1):
        occupation = _occupation(self)
        if attempt:
            # OpenTrustRegion's budget takes its first macro-iteration as already counted, as the
            # handoff is. A restart's first evaluation, at the new occupation, is a J/K build no
            # iteration has covered, so charge it here.
            self.iteration_ += 1
        otr_error = self.opentrustregion_scf(e_conv, d_conv, bool(attempt))
        if otr_error and otr_error != _OTR_MAX_ITER_ERROR:
            # A real failure, so it must not be swallowed by FAIL_ON_MAXITER, which callers
            # such as run_qchf set False to tolerate a deliberately truncated SCF. Raise before
            # rebuilding from orbitals the failure may have left half-updated.
            core.print_out(f"    OpenTrustRegion solver failed with error {otr_error}.\n")
            callback_error = self.opentrustregion_callback_error()
            if callback_error:
                raise PsiException(f"OpenTrustRegion solver failed with error {otr_error}, raised "
                                   f"in a callback into Psi4:\n{callback_error}")
            raise PsiException(f"OpenTrustRegion solver failed with error {otr_error}. See the "
                               f"output file for the solver's own message.")
        # OTR's last call to Psi4 is often a Hessian-vector product (cphf_Hx reuses JK)
        # so rebuild from converged orbs lest compute_E return garbage.
        self.form_D()
        self.form_G()
        self.form_F()

        # What this run finished on. The orbital gradient has to be measured here, where F
        # and D are each other's and the measure means what it does in the first-order loop.
        # After the diagonalization below, D is the projector onto eigenvectors of this F,
        # which sends the commutator to roundoff whether or not the orbitals are converged.
        self.otr_outcome_ = (self.opentrustregion_stopped_short(),
                             self.opentrustregion_final_grad_rms(),
                             self.opentrustregion_final_energy_delta(),
                             self.compute_orbital_gradient(
                                 False, core.get_option("SCF", "DIIS_MAX_VECS")))

        SCFE = self.compute_E()
        self.set_energies("Total Energy", SCFE)
        self.set_variable("SCF ITERATION ENERGY", SCFE)
        # The first-order loop refreshes core's copies every iteration and does not run here,
        # so they would otherwise still describe the pre-handoff point. SCF D NORM is read
        # back by DirectJK/CompositeJK to decide whether an incremental Fock build is safe.
        core.set_variable("SCF ITERATION ENERGY", SCFE)
        core.set_variable("SCF D NORM", self.otr_outcome_[3])

        # Report one entry per OTR macro-iteration so SCF ITERATIONS and SCF TOTAL ENERGIES
        # mean the same thing as for internal solver. OTR's macro-iter=0 is already counted, and
        # after the handoff also recorded; a restart's was only charged above, so record it here.
        otr_energies = list(self.opentrustregion_iteration_energies())
        otr_energies[-1:] = [SCFE]
        if attempt:
            self.iteration_energies.append(otr_energies[0])
        otr_energies = otr_energies[1:]
        self.iteration_energies.extend(otr_energies)
        self.iteration_ += len(otr_energies)

        # Leave the orbitals canonical for post-SCF methods, with the occupation where OTR used it.
        # Ask what the aufbau rule would do, and commit to that answer only when another attempt is going to use it.
        self.diagonalize_orbitals()
        # That diagonalization is only block diagonal to w/i convergence, so rebuild the density.
        self.form_D()

        if otr_error:
            # Out of macro-iterations. Reported below; nothing here can improve on it.
            break
        if _aufbau_occupation(self) == occupation:
            # The occupation OpenTrustRegion optimized is the aufbau one. Done.
            break
        if SCFE_prev_attempt is not None and SCFE >= SCFE_prev_attempt - e_conv:
            # Reconverging didn't lower the energy. Degenerate orbitals in different irreps --
            # singlet O2's pi*, say -- let the aufbau assignment flip back and forth without
            # changing the solution, so stop rather than oscillate. Where the energy actually
            # rose, the previous attempt is the better solution, so return that one.
            if SCFE > SCFE_prev_attempt + e_conv:
                state, self.otr_outcome_, n_energies = prev_attempt
                _restore_orbital_state(self, state)
                # The rejected attempt's energies describe orbitals the wavefunction no longer
                # holds, and SCF TOTAL ENERGIES would end on one of them.
                del self.iteration_energies[n_energies:]
                self.form_D()
                self.form_G()
                self.form_F()
                SCFE = self.compute_E()
                self.set_energies("Total Energy", SCFE)
                self.set_variable("SCF ITERATION ENERGY", SCFE)
                core.set_variable("SCF ITERATION ENERGY", SCFE)
                core.set_variable("SCF D NORM", self.otr_outcome_[3])
                core.print_out("    Note: reconverging at the aufbau occupation raised the energy; "
                               "keeping the previous solution.\n")
            break
        if self.iteration_ >= core.get_option("SCF", "MAXITER"):
            core.print_out("    Note: aufbau occupation differs from the one OpenTrustRegion "
                           "optimized, but MAXITER leaves no iteration to reconverge; accepting "
                           "this solution.\n")
            break
        if attempt == _OTR_MAX_OCCUPATION_RESTARTS:
            core.print_out(f"    Note: occupation still moving after "
                           f"{_OTR_MAX_OCCUPATION_RESTARTS} OpenTrustRegion restarts; "
                           f"accepting the last solution.\n")
            break
        # Committing to the aufbau assignment, so rebuild the density at it before handing
        # the orbitals back to OpenTrustRegion for another attempt. Keep this attempt's
        # solution first, in case the next one ends higher.
        prev_attempt = (_orbital_state(self), self.otr_outcome_, len(self.iteration_energies))
        self.find_occupation()
        self.form_D()
        SCFE_prev_attempt = SCFE
        core.print_out("    Note: aufbau occupation differs from the one OpenTrustRegion "
                       "optimized. Reconverging.\n")

    # Write orbitals after OTR call which is otherwise only done at the end of a polish.
    if core.has_option_changed("SCF", "ORBITALS_WRITE"):
        self.to_file(core.get_option("SCF", "ORBITALS_WRITE"))

    if otr_error:
        # Only running out of macro-iterations gets this far, so subject to FAIL_ON_MAXITER
        # like any other.
        core.print_out(f"    OpenTrustRegion solver returned error {otr_error}.\n")
        raise SCFConvergenceError("""SCF iterations""", self.iteration_, self, 0.0, 0.0)

    if self.otr_outcome_[0]:
        # OpenTrustRegion stopped where the objective function stops resolving, short of both
        # criteria. Its step acceptance compares the measured against the predicted energy
        # change, so it cannot verify progress it cannot measure.
        handling = core.get_option("SCF", "OTR_MAX_PRECISION_REACHED_HANDLING")
        if handling == "DIIS":
            # First-order iterations converge on the orbital gradient instead, which is still
            # measurable here, so let the caller finish the SCF with them.
            return NEEDS_POLISH
        if handling == "FAIL":
            text = shortfall(self, e_conv, d_conv)
            raise SCFConvergenceError(
                "SCF to the requested precision", self.iteration_, self, 0.0, 0.0,
                "OpenTrustRegion stopped at its floating-point precision floor: a further\n"
                "trust-region step changed the objective function by less than can be resolved\n"
                "at this energy, so the step could not be verified and was rejected. The\n"
                "orbitals are the best it found -- " + text + ".\n\n"
                "More iterations will not help. Loosen D_CONVERGENCE and E_CONVERGENCE to what\n"
                "this system supports, or set OTR_MAX_PRECISION_REACHED_HANDLING to DIIS to\n"
                "finish with first-order iterations, or to WARN to accept this result.")
        # WARN: the note in the output is the whole report; accept the orbitals as they are.



class Polish:
    """What a first-order polish after a precision-floor stop has to remember.

    Built once the solver hands back :py:data:`NEEDS_POLISH`, and consulted again only if
    the iterations it covers run out of MAXITER without converging.
    """

    def __init__(self, wfn, e_conv, d_conv):
        self.shortfall = shortfall(wfn, e_conv, d_conv)
        self.start_iteration = wfn.iteration_
        # Keep what OpenTrustRegion converged. First-order iterations usually improve on it,
        # but from a solution this close they can also wander off -- a pinned-occupation ROHF
        # will -- and the better orbitals should be the ones left behind. run() measured its
        # gradient the way the first-order loop does, so this comparison is like for like.
        self.entry_dnorm = wfn.otr_outcome_[3]
        self.entry_state = _orbital_state(wfn)

    def give_up(self, wfn, Ediff, Dnorm, e_conv, d_conv):
        """Keep the better of the two solutions, then report the SCF as unconverged."""
        if Dnorm > self.entry_dnorm:
            _restore_orbital_state(wfn, self.entry_state)
            wfn.form_D()
            wfn.form_G()
            wfn.form_F()
            SCFE = wfn.compute_E()
            wfn.set_energies("Total Energy", SCFE)
            wfn.set_variable("SCF ITERATION ENERGY", SCFE)
            core.set_variable("SCF ITERATION ENERGY", SCFE)
            core.set_variable("SCF D NORM", self.entry_dnorm)
            core.print_out(
                f"    The first-order iterations left the orbitals worse than\n"
                f"    OpenTrustRegion did ({Dnorm:.2e} against "
                f"{self.entry_dnorm:.2e}); keeping its solution.\n")
            Dnorm = self.entry_dnorm
        raise SCFConvergenceError(
            "SCF to the requested precision", wfn.iteration_, wfn, Ediff, Dnorm,
            "OpenTrustRegion stopped at its floating-point precision floor and the\n"
            "first-order iterations that followed did not reach the thresholds either.\n\n"
            "Where OpenTrustRegion stopped: " + self.shortfall + ".\n"
            f"After {wfn.iteration_ - self.start_iteration} first-order iterations: "
            f"last energy change {abs(Ediff):.2e} against E_CONVERGENCE {e_conv:.2e}, "
            f"Dnorm {Dnorm:.2e} against D_CONVERGENCE {d_conv:.2e}.\n\n"
            "A second-order and a first-order method have both failed to reach the\n"
            "thresholds within MAXITER, so this request is past what this system and\n"
            "basis support. Loosen D_CONVERGENCE and E_CONVERGENCE.")

def selected(wfn, soscf_enabled, potentials):
    """Whether OpenTrustRegion drives this SCF's second-order iterations.

    ``potentials`` is :py:func:`~scf_iterator._python_applied_potentials`, the per-iteration
    Fock contributions shared with the OpenOrbitalOptimizer guard.

    False where it is not the package asked for, and where it is but the calculation has
    something it cannot carry. It takes over the iteration loop at the handoff, so anything
    applied per-iteration from Python (EFP/PCM/DDX/PE) or through find_occupation (MOM) would
    stop being applied thereafter. The GRAC shift splices its potential into V_xc outside the
    kernel cphf_Hx differentiates, so the Hessian never sees it. Second-order SCF in general
    is refused for CUHF, meta/VV10 and FRAC by _validate_soscf, and for sn-LinK and IncFock
    in proc.scf_helper.
    """
    if not (soscf_enabled
            and core.get_option("SCF", "SECOND_ORDER_ORBITAL_OPTIMIZER_PACKAGE")
            in ["OTR", "OPENTRUSTREGION"]):
        return False

    # Excited-state MOM is MOM_START with MOM_OCC; MOM_START alone still reorders the
    # occupation by overlap, which is neither what OpenTrustRegion optimizes at nor what
    # _aufbau_occupation predicts, so the restart loop would reconverge at an occupation
    # find_occupation then overrides.
    mom_enabled = wfn.MOM_excited_ or core.get_option("SCF", "MOM_START") != 0
    grac_enabled = (core.get_option("SCF", "DFT_GRAC_SHIFT") != 0.0
                    or core.get_option("SAPT", "SAPT_DFT_GRAC_COMPUTE") != "NONE")
    if mom_enabled or grac_enabled or any(potentials.values()):
        core.print_out("    Note: OpenTrustRegion not compatible with at least one of the following. Falling back to second_order_orbital_optimizer_package=internal\n")
        core.print_out(f"          mom={mom_enabled}, grac={grac_enabled},\n")
        core.print_out("          "
                       + ", ".join(f"{name}={value}" for name, value in potentials.items())
                       + "\n")
        return False

    return True
