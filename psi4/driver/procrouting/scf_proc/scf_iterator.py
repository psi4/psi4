
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
"""
The SCF iteration functions
"""
import numpy as np

from psi4 import core

from ... import p4util
from ...constants import constants
from ...p4util.exceptions import SCFConvergenceError, ValidationError
from ..solvent.efp import get_qm_atoms_opts, modify_Fock_induced, modify_Fock_permanent

#import logging
#logger = logging.getLogger("scf.scf_iterator")
#logger.setLevel(logging.DEBUG)

# Q: I expect more local settings of options for part of SCF.
#    For convcrit, do we want:
#   (A) easy to grep
#    with p4util.OptionsStateCM(['SCF', 'E_CONVERGENCE'], ['SCF', 'D_CONVERGENCE']):
#        core.set_local_option('SCF', 'E_CONVERGENCE', 1.e-5)
#        core.set_local_option('SCF', 'D_CONVERGENCE', 1.e-4)
#        self.iterations()
#
#   or (B) functional. options never touched
#    self.iterations(e_conv=1.e-5, d_conv=1.e-4)


def scf_compute_energy(self):
    """Base class Wavefunction requires this function. Here it is
    simply a wrapper around initialize(), iterations(), finalize_energy(). It
    returns the SCF energy computed by finalize_energy().

    """
    if core.get_option('SCF', 'DF_SCF_GUESS') and (core.get_global_option('SCF_TYPE') == 'DIRECT'):
        # speed up DIRECT algorithm (recomputes full (non-DF) integrals
        #   each iter) by first converging via fast DF iterations, then
        #   fully converging in fewer slow DIRECT iterations. aka Andy trick 2.0
        core.print_out("  Starting with a DF guess...\n\n")
        with p4util.OptionsStateCM(['SCF_TYPE']):
            core.set_global_option('SCF_TYPE', 'DF')
            self.initialize()
            try:
                self.iterations()
            except SCFConvergenceError:
                self.finalize()
                raise SCFConvergenceError("""SCF DF preiterations""", self.iteration_, self, 0, 0)
        core.print_out("\n  DF guess converged.\n\n")

        # reset the DIIS & JK objects in prep for DIRECT
        if self.initialized_diis_manager_:
            self.diis_manager_.reset_subspace()
        self.initialize_jk(self.memory_jk_)
    else:
        self.initialize()
    self.iteration_energies = []

    try:
        self.iterations()
    except SCFConvergenceError as e:
        if core.get_option("SCF", "FAIL_ON_MAXITER"):
            core.print_out("  Failed to converge.\n")
            # energy = 0.0
            # A P::e fn to either throw or protest upon nonconvergence
            # die_if_not_converged()
            raise e
        else:
            core.print_out("  Energy and/or wave function did not converge, but proceeding anyway.\n\n")
    else:
        core.print_out("  Energy and wave function converged.\n\n")

    scf_energy = self.finalize_energy()
    return scf_energy


def _build_jk(wfn, memory):
    jk = core.JK.build(wfn.get_basisset("ORBITAL"),
                       aux=wfn.get_basisset("DF_BASIS_SCF"),
                       do_wK=wfn.functional().is_x_lrc(),
                       memory=memory)
    return jk


def initialize_jk(self, memory, jk=None):

    functional = self.functional()
    if jk is None:
        jk = _build_jk(self, memory)

    self.set_jk(jk)

    jk.set_print(self.get_print())
    jk.set_memory(memory)
    jk.set_do_K(functional.is_x_hybrid())
    jk.set_do_wK(functional.is_x_lrc())
    jk.set_omega(functional.x_omega())

    jk.set_omega_alpha(functional.x_alpha())
    jk.set_omega_beta(functional.x_beta())

    jk.initialize()
    jk.print_header()


def scf_initialize(self):
    """Specialized initialization, compute integrals and does everything to prepare for iterations"""

    # Figure out memory distributions

    # Get memory in terms of doubles
    total_memory = (core.get_memory() / 8) * core.get_global_option("SCF_MEM_SAFETY_FACTOR")

    # Figure out how large the DFT collocation matrices are
    vbase = self.V_potential()
    if vbase:
        collocation_size = vbase.grid().collocation_size()
        if vbase.functional().ansatz() == 1:
            collocation_size *= 4  # First derivs
        elif vbase.functional().ansatz() == 2:
            collocation_size *= 10  # Second derivs
    else:
        collocation_size = 0

    # Change allocation for collocation matrices based on DFT type
    initialize_jk_obj = False
    if isinstance(self.jk(), core.JK):
        core.print_out("\nRe-using passed JK object instead of rebuilding\n")
        jk = self.jk()
    else:
        initialize_jk_obj = True
        jk = _build_jk(self, total_memory)
    jk_size = jk.memory_estimate()

    # Give remaining to collocation
    if total_memory > jk_size:
        collocation_memory = total_memory - jk_size
    # Give up to 10% to collocation
    elif (total_memory * 0.1) > collocation_size:
        collocation_memory = collocation_size
    else:
        collocation_memory = total_memory * 0.1

    if collocation_memory > collocation_size:
        collocation_memory = collocation_size

    # Set constants
    self.iteration_ = 0
    self.memory_jk_ = int(total_memory - collocation_memory)
    self.memory_collocation_ = int(collocation_memory)

    if self.get_print():
        core.print_out("  ==> Integral Setup <==\n\n")

    # Initialize EFP
    efp_enabled = hasattr(self.molecule(), 'EFP')
    if efp_enabled:
        # EFP: Set QM system, options, and callback. Display efp geom in [A]
        efpobj = self.molecule().EFP
        core.print_out(efpobj.banner())
        core.print_out(efpobj.geometry_summary(units_to_bohr=constants.bohr2angstroms))

        efpptc, efpcoords, efpopts = get_qm_atoms_opts(self.molecule())
        efpobj.set_point_charges(efpptc, efpcoords)
        efpobj.set_opts(efpopts, label='psi', append='psi')

        efpobj.set_electron_density_field_fn(efp_field_fn)

    # Initialize all integrals and perform the first guess
    if self.attempt_number_ == 1:
        mints = core.MintsHelper(self.basisset())

        if initialize_jk_obj:
            self.initialize_jk(self.memory_jk_, jk=jk)
        if self.V_potential():
            self.V_potential().build_collocation_cache(self.memory_collocation_)
        core.timer_on("HF: Form core H")
        self.form_H()
        core.timer_off("HF: Form core H")

        if efp_enabled:
            # EFP: Add in permanent moment contribution and cache
            core.timer_on("HF: Form Vefp")
            verbose = core.get_option('SCF', "PRINT")
            Vefp = modify_Fock_permanent(self.molecule(), mints, verbose=verbose - 1)
            Vefp = core.Matrix.from_array(Vefp)
            self.H().add(Vefp)
            Horig = self.H().clone()
            self.Horig = Horig
            core.print_out("  QM/EFP: iterating Total Energy including QM/EFP Induction\n")
            core.timer_off("HF: Form Vefp")

        core.timer_on("HF: Form S/X")
        self.form_Shalf()
        core.timer_off("HF: Form S/X")

        core.print_out("\n  ==> Pre-Iterations <==\n\n")

        # force SCF_SUBTYPE to AUTO during SCF guess
        optstash = p4util.OptionsState(["SCF", "SCF_SUBTYPE"])
        core.set_local_option("SCF", "SCF_SUBTYPE", "AUTO")

        core.timer_on("HF: Guess")
        self.guess()
        core.timer_off("HF: Guess")

        optstash.restore()

        # Print out initial docc/socc/etc data
        if self.get_print():
            lack_occupancy = core.get_local_option('SCF', 'GUESS') in ['SAD']
            if core.get_global_option('GUESS') in ['SAD']:
                lack_occupancy = core.get_local_option('SCF', 'GUESS') in ['AUTO']
                self.print_preiterations(small=lack_occupancy)
            else:
                self.print_preiterations(small=lack_occupancy)

    else:
        # We're reading the orbitals from the previous set of iterations.
        self.form_D()
        self.set_energies("Total Energy", self.compute_initial_E())

    # turn off VV10 for iterations
    if core.get_option('SCF', "DFT_VV10_POSTSCF") and self.functional().vv10_b() > 0.0:
        core.print_out("  VV10: post-SCF option active \n \n")
        self.functional().set_lock(False)
        self.functional().set_do_vv10(False)
        self.functional().set_lock(True)

    # Print iteration header
    is_dfjk = core.get_global_option('SCF_TYPE').endswith('DF')
    diis_rms = core.get_option('SCF', 'DIIS_RMS_ERROR')
    core.print_out("  ==> Iterations <==\n\n")
    core.print_out("%s                        Total Energy        Delta E     %s |[F,P]|\n\n" %
                   ("   " if is_dfjk else "", "RMS" if diis_rms else "MAX"))


def scf_iterate(self, e_conv=None, d_conv=None):

    is_dfjk = core.get_global_option('SCF_TYPE').endswith('DF')
    verbose = core.get_option('SCF', "PRINT")
    reference = core.get_option('SCF', "REFERENCE")

    # self.member_data_ signals are non-local, used internally by c-side fns
    self.diis_enabled_ = self.validate_diis()
    self.MOM_excited_ = _validate_MOM()
    self.diis_start_ = core.get_option('SCF', 'DIIS_START')
    damping_enabled = _validate_damping()
    soscf_enabled = _validate_soscf()
    frac_enabled = _validate_frac()
    efp_enabled = hasattr(self.molecule(), 'EFP')
    cosx_enabled = "COSX" in core.get_option('SCF', 'SCF_TYPE')
    ooo_scf = core.get_global_option('OOO_SCF')
    if ooo_scf:
        # SAD needs some special work since the guess doesn't actually make the orbitals in Psi4
        if self.sad_ and self.iteration_ <= 0:
            self.form_G()
            self.form_initial_F()
            self.form_initial_C()
            self.reset_occupation()
            self.find_occupation()
        self.openorbital_scf()
        self.set_energies("Total Energy", self.compute_E())
        return

    # does the JK algorithm use severe screening approximations for early SCF iterations?
    early_screening = False
    if cosx_enabled:
        early_screening = True
        self.jk().set_COSX_grid("Initial")

    # maximum number of scf iterations to run after early screening is disabled
    scf_maxiter_post_screening = core.get_option('SCF', 'COSX_MAXITER_FINAL')

    if scf_maxiter_post_screening < -1:
        raise ValidationError('COSX_MAXITER_FINAL ({}) must be -1 or above. If you wish to attempt full SCF converge on the final COSX grid, set COSX_MAXITER_FINAL to -1.'.format(scf_maxiter_post_screening))

    # has early_screening changed from True to False?
    early_screening_disabled = False

    # SCF iterations!
    SCFE_old = 0.0
    Dnorm = 0.0
    scf_iter_post_screening = 0
    while True:
        self.iteration_ += 1

        diis_performed = False
        soscf_performed = False
        self.frac_performed_ = False
        #self.MOM_performed_ = False  # redundant from common_init()

        self.save_density_and_energy()

        if efp_enabled:
            # EFP: Add efp contribution to Fock matrix
            self.H().copy(self.Horig)
            global mints_psi4_yo
            mints_psi4_yo = core.MintsHelper(self.basisset())
            Vefp = modify_Fock_induced(self.molecule().EFP, mints_psi4_yo, verbose=verbose - 1)
            Vefp = core.Matrix.from_array(Vefp)
            self.H().add(Vefp)

        SCFE = 0.0
        self.clear_external_potentials()

        # Two-electron contribution to Fock matrix from self.jk()
        core.timer_on("HF: Form G")
        self.form_G()
        core.timer_off("HF: Form G")

        # Check if special J/K construction algorithms were used
        incfock_performed = hasattr(self.jk(), "do_incfock_iter") and self.jk().do_incfock_iter()
        upcm = 0.0
        if core.get_option('SCF', 'PCM'):
            calc_type = core.PCM.CalcType.Total
            if core.get_option("PCM", "PCM_SCF_TYPE") == "SEPARATE":
                calc_type = core.PCM.CalcType.NucAndEle
            Dt = self.Da().clone()
            Dt.add(self.Db())
            upcm, Vpcm = self.get_PCM().compute_PCM_terms(Dt, calc_type)
            SCFE += upcm
            self.push_back_external_potential(Vpcm)
        self.set_variable("PCM POLARIZATION ENERGY", upcm)  # P::e PCM
        self.set_energies("PCM Polarization", upcm)

        uddx = 0.0
        if core.get_option('SCF', 'DDX'):
            Dt = self.Da().clone()
            Dt.add(self.Db())
            uddx, Vddx, self.ddx_state = self.ddx.get_solvation_contributions(Dt, self.ddx_state)
            SCFE += uddx
            self.push_back_external_potential(Vddx)
        self.set_variable("DD SOLVATION ENERGY", uddx)  # P::e DDX
        self.set_energies("DD Solvation Energy", uddx)

        upe = 0.0
        if core.get_option('SCF', 'PE'):
            Dt = self.Da().clone()
            Dt.add(self.Db())
            upe, Vpe = self.pe_state.get_pe_contribution(
                Dt, elec_only=False
            )
            SCFE += upe
            self.push_back_external_potential(Vpe)
        self.set_variable("PE ENERGY", upe)  # P::e PE
        self.set_energies("PE Energy", upe)

        core.timer_on("HF: Form F")
        # SAD: since we don't have orbitals yet, we might not be able
        # to form the real Fock matrix. Instead, build an initial one
        if (self.iteration_ == 0) and self.sad_:
            self.form_initial_F()
        else:
            self.form_F()
        core.timer_off("HF: Form F")

        if verbose > 3:
            self.Fa().print_out()
            self.Fb().print_out()

        SCFE += self.compute_E()
        if efp_enabled:
            global efp_Dt_psi4_yo

            # EFP: Add efp contribution to energy
            efp_Dt_psi4_yo = self.Da().clone()
            efp_Dt_psi4_yo.add(self.Db())
            SCFE += self.molecule().EFP.get_wavefunction_dependent_energy()

        self.set_energies("Total Energy", SCFE)
        core.set_variable("SCF ITERATION ENERGY", SCFE)
        self.iteration_energies.append(SCFE)

        Ediff = SCFE - SCFE_old
        SCFE_old = SCFE

        status = []

        # Check if we are doing SOSCF
        if (soscf_enabled and (self.iteration_ >= 3) and (Dnorm < core.get_option('SCF', 'SOSCF_START_CONVERGENCE'))):
            Dnorm = self.compute_orbital_gradient(False, core.get_option('SCF', 'DIIS_MAX_VECS'))
            diis_performed = False
            if self.functional().needs_xc():
                base_name = "SOKS, nmicro="
            else:
                base_name = "SOSCF, nmicro="

            if not _converged(Ediff, Dnorm, e_conv=e_conv, d_conv=d_conv):
                nmicro = self.soscf_update(core.get_option('SCF', 'SOSCF_CONV'),
                                           core.get_option('SCF', 'SOSCF_MIN_ITER'),
                                           core.get_option('SCF', 'SOSCF_MAX_ITER'),
                                           core.get_option('SCF', 'SOSCF_PRINT'))
                # if zero, the soscf call bounced for some reason
                soscf_performed = (nmicro > 0)

                if soscf_performed:
                    self.find_occupation()
                    status.append(base_name + str(nmicro))
                else:
                    if verbose > 0:
                        core.print_out("Did not take a SOSCF step, using normal convergence methods\n")

            else:
                # need to ensure orthogonal orbitals and set epsilon
                status.append(base_name + "conv")
                core.timer_on("HF: Form C")
                self.form_C()
                core.timer_off("HF: Form C")
                soscf_performed = True  # Stops DIIS

        if not soscf_performed:
            # Normal convergence procedures if we do not do SOSCF

            # SAD: form initial orbitals from the initial Fock matrix, and
            # reset the occupations. The reset is necessary because SAD
            # nalpha_ and nbeta_ are not guaranteed physical.
            # From here on, the density matrices are correct.
            if (self.iteration_ == 0) and self.sad_:
                self.form_initial_C()
                self.reset_occupation()
                self.find_occupation()

            else:
                # Run DIIS
                core.timer_on("HF: DIIS")
                diis_performed = False
                add_to_diis_subspace = self.diis_enabled_ and self.iteration_ >= self.diis_start_

                Dnorm = self.compute_orbital_gradient(add_to_diis_subspace, core.get_option('SCF', 'DIIS_MAX_VECS'))

                if add_to_diis_subspace:
                    for engine_used in self.diis(Dnorm):
                        status.append(engine_used)

                core.timer_off("HF: DIIS")

                if verbose > 4 and diis_performed:
                    core.print_out("  After DIIS:\n")
                    self.Fa().print_out()
                    self.Fb().print_out()

                # frac, MOM invoked here from Wfn::HF::find_occupation
                core.timer_on("HF: Form C")
                level_shift = core.get_option("SCF", "LEVEL_SHIFT")
                if level_shift > 0 and Dnorm > core.get_option('SCF', 'LEVEL_SHIFT_CUTOFF'):
                    status.append("SHIFT")
                    self.form_C(level_shift)
                else:
                    self.form_C()
                core.timer_off("HF: Form C")

                if self.MOM_performed_:
                    status.append("MOM")

                if self.frac_performed_:
                    status.append("FRAC")

                if incfock_performed:
                    status.append("INCFOCK")

                # Reset occupations if necessary
                if (self.iteration_ == 0) and self.reset_occ_:
                    self.reset_occupation()
                    self.find_occupation()

        # Form new density matrix
        core.timer_on("HF: Form D")
        self.form_D()
        core.timer_off("HF: Form D")

        self.set_variable("SCF ITERATION ENERGY", SCFE)
        core.set_variable("SCF D NORM", Dnorm)

        # After we've built the new D, damp the update
        if (damping_enabled and self.iteration_ > 1 and Dnorm > core.get_option('SCF', 'DAMPING_CONVERGENCE')):
            damping_percentage = core.get_option('SCF', "DAMPING_PERCENTAGE")
            self.damping_update(damping_percentage * 0.01)
            status.append("DAMP={}%".format(round(damping_percentage)))

        if core.has_option_changed("SCF", "ORBITALS_WRITE"):
            filename = core.get_option("SCF", "ORBITALS_WRITE")
            self.to_file(filename)

        if verbose > 3:
            self.Ca().print_out()
            self.Cb().print_out()
            self.Da().print_out()
            self.Db().print_out()

        # Print out the iteration
        core.print_out(
            "   @%s%s iter %3s: %20.14f   %12.5e   %-11.5e %s\n" %
            ("DF-" if is_dfjk else "", reference, "SAD" if
             ((self.iteration_ == 0) and self.sad_) else self.iteration_, SCFE, Ediff, Dnorm, '/'.join(status)))

        # if a an excited MOM is requested but not started, don't stop yet
        # Note that MOM_performed_ just checks initialization, and our convergence measures used the pre-MOM orbitals
        if self.MOM_excited_ and ((not self.MOM_performed_) or self.iteration_ == core.get_option('SCF', "MOM_START")):
            continue

        # if a fractional occupation is requested but not started, don't stop yet
        if frac_enabled and not self.frac_performed_:
            continue

        # have we completed our post-early screening SCF iterations?
        if early_screening_disabled:
            scf_iter_post_screening += 1
            if scf_iter_post_screening >= scf_maxiter_post_screening and scf_maxiter_post_screening > 0:
                break

        # Call any postiteration callbacks
        if not ((self.iteration_ == 0) and self.sad_) and _converged(Ediff, Dnorm, e_conv=e_conv, d_conv=d_conv):

            if early_screening:

                # we've reached convergence with early screning enabled; disable it
                early_screening = False

                # make note of the change to early screening; next SCF iteration(s) will be the last
                early_screening_disabled = True

                # cosx uses the largest grid for its final SCF iteration(s)
                if cosx_enabled:
                    self.jk().set_COSX_grid("Final")

                # clear any cached matrices associated with incremental fock construction
                # the change in the screening spoils the linearity in the density matrix
                if hasattr(self.jk(), 'clear_D_prev'):
                    self.jk().clear_D_prev()

                if scf_maxiter_post_screening == 0:
                    break
                else:
                    core.print_out("  Energy and wave function converged with early screening.\n")
                    core.print_out("  Continuing SCF iterations with tighter screening.\n\n")
            else:
                break

        if self.iteration_ >= core.get_option('SCF', 'MAXITER'):
            raise SCFConvergenceError("""SCF iterations""", self.iteration_, self, Ediff, Dnorm)


def scf_finalize_energy(self):
    """Performs stability analysis and calls back SCF with new guess
    if needed, Returns the SCF energy. This function should be called
    once orbitals are ready for energy/property computations, usually
    after iterations() is called.

    """

    # post-scf vv10 correlation
    if core.get_option('SCF', "DFT_VV10_POSTSCF") and self.functional().vv10_b() > 0.0:
        self.functional().set_lock(False)
        self.functional().set_do_vv10(True)
        self.functional().set_lock(True)
        core.print_out("  ==> Computing Non-Self-Consistent VV10 Energy Correction <==\n\n")
        SCFE = 0.0
        self.form_V()
        SCFE += self.compute_E()
        self.set_energies("Total Energy", SCFE)

    # Perform wavefunction stability analysis before doing
    # anything on a wavefunction that may not be truly converged.
    if core.get_option('SCF', 'STABILITY_ANALYSIS') != "NONE":

        # We need the integral file, make sure it is written and
        # compute it if needed
        if core.get_option('SCF', 'REFERENCE') not in {"UHF", "UKS"}:
            # Don't bother computing needed integrals if we can't do anything with them.
            if self.functional().needs_xc():
                raise ValidationError("Stability analysis not yet supported for XC functionals.")

            #psio = core.IO.shared_object()
            #psio.open(constants.PSIF_SO_TEI, 1)  # PSIO_OPEN_OLD
            #try:
            #    psio.tocscan(constants.PSIF_SO_TEI, "IWL Buffers")
            #except TypeError:
            #    # "IWL Buffers" actually found but psio_tocentry can't be returned to Py
            #    psio.close(constants.PSIF_SO_TEI, 1)
            #else:
            #    # tocscan returned None
            #    psio.close(constants.PSIF_SO_TEI, 1)

            # logic above foiled by psio_tocentry not returning None<--nullptr in pb11 2.2.1
            #   so forcibly recomputing for now until stability revamp
            core.print_out("    SO Integrals not on disk. Computing...")
            mints = core.MintsHelper(self.basisset())

            mints.integrals()
            core.print_out("done.\n")

            # Q: Not worth exporting all the layers of psio, right?

        follow = self.stability_analysis()

        while follow and self.attempt_number_ <= core.get_option('SCF', 'MAX_ATTEMPTS'):
            self.attempt_number_ += 1
            core.print_out("    Running SCF again with the rotated orbitals.\n")

            if self.initialized_diis_manager_:
                self.diis_manager_.reset_subspace()
            # reading the rotated orbitals in before starting iterations
            self.form_D()
            self.set_energies("Total Energy", self.compute_initial_E())
            self.iterations()
            follow = self.stability_analysis()

        if follow and self.attempt_number_ > core.get_option('SCF', 'MAX_ATTEMPTS'):
            core.print_out("    There's still a negative eigenvalue. Try modifying FOLLOW_STEP_SCALE\n")
            core.print_out("    or increasing MAX_ATTEMPTS (not available for PK integrals).\n")

    # At this point, we are not doing any more SCF cycles
    #   and we can compute and print final quantities.

    if hasattr(self.molecule(), 'EFP'):
        efpobj = self.molecule().EFP

        efpobj.compute()  # do_gradient=do_gradient)
        efpene = efpobj.get_energy(label='psi')
        efp_wfn_independent_energy = efpene['total'] - efpene['ind']
        self.set_energies("EFP", efpene['total'])

        SCFE = self.get_energies("Total Energy")
        SCFE += efp_wfn_independent_energy
        self.set_energies("Total Energy", SCFE)
        core.print_out(efpobj.energy_summary(scfefp=SCFE, label='psi'))

        self.set_variable("EFP ELST ENERGY", efpene['electrostatic'] + efpene['charge_penetration'] + efpene['electrostatic_point_charges'])  # P::e EFP
        self.set_variable("EFP IND ENERGY", efpene['polarization'])  # P::e EFP
        self.set_variable("EFP DISP ENERGY", efpene['dispersion'])  # P::e EFP
        self.set_variable("EFP EXCH ENERGY", efpene['exchange_repulsion'])  # P::e EFP
        self.set_variable("EFP TOTAL ENERGY", efpene['total'])  # P::e EFP
        self.set_variable("CURRENT ENERGY", efpene['total'])  # P::e EFP

    core.print_out("\n  ==> Post-Iterations <==\n\n")

    if self.V_potential():
        quad = self.V_potential().quadrature_values()
        rho_a = quad['RHO_A']/2 if self.same_a_b_dens() else quad['RHO_A']
        rho_b = quad['RHO_B']/2 if self.same_a_b_dens() else quad['RHO_B']
        rho_ab = (rho_a + rho_b)
        self.set_variable("GRID ELECTRONS TOTAL",rho_ab)  # P::e SCF
        self.set_variable("GRID ELECTRONS ALPHA",rho_a)  # P::e SCF
        self.set_variable("GRID ELECTRONS BETA",rho_b)  # P::e SCF
        dev_a = rho_a - self.nalpha()
        dev_b = rho_b - self.nbeta()
        core.print_out(f"   Electrons on quadrature grid:\n")
        if self.same_a_b_dens():
            core.print_out(f"      Ntotal   = {rho_ab:15.10f} ; deviation = {dev_b+dev_a:.3e} \n\n")
        else:
            core.print_out(f"      Nalpha   = {rho_a:15.10f} ; deviation = {dev_a:.3e}\n")
            core.print_out(f"      Nbeta    = {rho_b:15.10f} ; deviation = {dev_b:.3e}\n")
            core.print_out(f"      Ntotal   = {rho_ab:15.10f} ; deviation = {dev_b+dev_a:.3e} \n\n")
        if ((dev_b+dev_a) > 0.1):
            core.print_out("   WARNING: large deviation in the electron count on grid detected. Check grid size!")
    self.check_phases()
    self.compute_spin_contamination()
    self.frac_renormalize()
    reference = core.get_option("SCF", "REFERENCE")

    energy = self.get_energies("Total Energy")

    #    fail_on_maxiter = core.get_option("SCF", "FAIL_ON_MAXITER")
    #    if converged or not fail_on_maxiter:
    #
    #        if print_lvl > 0:
    #            self.print_orbitals()
    #
    #        if converged:
    #            core.print_out("  Energy converged.\n\n")
    #        else:
    #            core.print_out("  Energy did not converge, but proceeding anyway.\n\n")

    if core.get_option('SCF', 'PRINT') > 0:
        self.print_orbitals()

    is_dfjk = core.get_global_option('SCF_TYPE').endswith('DF')
    core.print_out("  @%s%s Final Energy: %20.14f" % ('DF-' if is_dfjk else '', reference, energy))
    # if (perturb_h_) {
    #     core.print_out(" with %f %f %f perturbation" %
    #                    (dipole_field_strength_[0], dipole_field_strength_[1], dipole_field_strength_[2]))
    # }
    core.print_out("\n\n")
    self.print_energies()

    # force list into Matrix for storage
    iteration_energies = np.array(self.iteration_energies).reshape(-1, 1)
    iteration_energies = core.Matrix.from_array(iteration_energies)
    core.set_variable("SCF TOTAL ENERGIES", core.Matrix.from_array(iteration_energies))
    self.set_variable("SCF TOTAL ENERGIES", core.Matrix.from_array(iteration_energies))

    self.clear_external_potentials()
    if core.get_option('SCF', 'PCM'):
        calc_type = core.PCM.CalcType.Total
        if core.get_option("PCM", "PCM_SCF_TYPE") == "SEPARATE":
            calc_type = core.PCM.CalcType.NucAndEle
        Dt = self.Da().clone()
        Dt.add(self.Db())
        _, Vpcm = self.get_PCM().compute_PCM_terms(Dt, calc_type)
        self.push_back_external_potential(Vpcm)
        # Set callback function for CPSCF
        self.set_external_cpscf_perturbation("PCM", lambda pert_dm : self.get_PCM().compute_V(pert_dm))

    if core.get_option('SCF', 'PE'):
        Dt = self.Da().clone()
        Dt.add(self.Db())
        _, Vpe = self.pe_state.get_pe_contribution(
            Dt, elec_only=False
        )
        self.push_back_external_potential(Vpe)
        # Set callback function for CPSCF
        self.set_external_cpscf_perturbation("PE", lambda pert_dm : self.pe_state.get_pe_contribution(pert_dm, elec_only=True)[1])

    if core.get_option('SCF', 'DDX'):
        Dt = self.Da().clone()
        Dt.add(self.Db())
        Vddx = self.ddx.get_solvation_contributions(Dt)[1]
        self.push_back_external_potential(Vddx)
        # Set callback function for CPSCF
        self.set_external_cpscf_perturbation(
            "DDX", lambda pert_dm : self.ddx.get_solvation_contributions(pert_dm, elec_only=True, nonequilibrium=True)[1])

    # Orbitals are always saved, in case an MO guess is requested later
    # save_orbitals()

    # Shove variables into global space
    for k, v in self.variables().items():
        core.set_variable(k, v)

    # TODO re-enable
    self.finalize()
    if self.V_potential():
        self.V_potential().clear_collocation_cache()

    core.print_out("\nComputation Completed\n")
    core.del_variable("SCF D NORM")

    return energy


def scf_print_energies(self):
    enuc = self.get_energies('Nuclear')
    e1 = self.get_energies('One-Electron')
    e2 = self.get_energies('Two-Electron')
    exc = self.get_energies('XC')
    ed = self.get_energies('-D')
    self.del_variable('-D Energy')
    evv10 = self.get_energies('VV10')
    eefp = self.get_energies('EFP')
    epcm = self.get_energies('PCM Polarization')
    edd = self.get_energies('DD Solvation Energy')
    epe = self.get_energies('PE Energy')
    ke = self.get_energies('Kinetic')

    hf_energy = enuc + e1 + e2
    dft_energy = hf_energy + exc + ed + evv10
    total_energy = dft_energy + eefp + epcm + edd + epe
    full_qm = (not core.get_option('SCF', 'PCM') and not core.get_option('SCF', 'DDX') and not core.get_option('SCF', 'PE')
               and not hasattr(self.molecule(), 'EFP'))

    core.print_out("   => Energetics <=\n\n")
    core.print_out("    Nuclear Repulsion Energy =        {:24.16f}\n".format(enuc))
    core.print_out("    One-Electron Energy =             {:24.16f}\n".format(e1))
    core.print_out("    Two-Electron Energy =             {:24.16f}\n".format(e2))
    if self.functional().needs_xc():
        core.print_out("    DFT Exchange-Correlation Energy = {:24.16f}\n".format(exc))
        core.print_out("    Empirical Dispersion Energy =     {:24.16f}\n".format(ed))
        core.print_out("    VV10 Nonlocal Energy =            {:24.16f}\n".format(evv10))
    if core.get_option('SCF', 'PCM'):
        core.print_out("    PCM Polarization Energy =         {:24.16f}\n".format(epcm))
    if core.get_option('SCF', 'DDX'):
        core.print_out("    DD Solvation Energy =            {:24.16f}\n".format(edd))
    if core.get_option('SCF', 'PE'):
        core.print_out("    PE Energy =                       {:24.16f}\n".format(epe))
    if hasattr(self.molecule(), 'EFP'):
        core.print_out("    EFP Energy =                      {:24.16f}\n".format(eefp))
    core.print_out("    Total Energy =                    {:24.16f}\n".format(total_energy))

    if core.get_option('SCF', 'PE'):
        core.print_out(self.pe_state.cppe_state.summary_string)

    self.set_variable("NUCLEAR REPULSION ENERGY", enuc)  # P::e SCF
    self.set_variable("ONE-ELECTRON ENERGY", e1)  # P::e SCF
    self.set_variable("TWO-ELECTRON ENERGY", e2)  # P::e SCF
    if self.functional().needs_xc():
        self.set_variable("DFT XC ENERGY", exc)  # P::e SCF
        self.set_variable("DFT VV10 ENERGY", evv10)  # P::e SCF
        self.set_variable("DFT FUNCTIONAL TOTAL ENERGY", hf_energy + exc + evv10)  # P::e SCF
        #self.set_variable(self.functional().name() + ' FUNCTIONAL TOTAL ENERGY', hf_energy + exc + evv10)
        self.set_variable("DFT TOTAL ENERGY", dft_energy)  # overwritten later for DH  # P::e SCF
    else:
        potential = total_energy - ke
        self.set_variable("HF KINETIC ENERGY", ke)  # P::e SCF
        self.set_variable("HF POTENTIAL ENERGY", potential)  # P::e SCF
        if full_qm:
            self.set_variable("HF VIRIAL RATIO", - potential / ke)  # P::e SCF
        self.set_variable("HF TOTAL ENERGY", hf_energy)  # P::e SCF
    if hasattr(self, "_disp_functor"):
        self.set_variable("DISPERSION CORRECTION ENERGY", ed)  # P::e SCF
    #if abs(ed) > 1.0e-14:
    #    for pv, pvv in self.variables().items():
    #        if abs(pvv - ed) < 1.0e-14:
    #            if pv.endswith('DISPERSION CORRECTION ENERGY') and pv.startswith(self.functional().name()):
    #                fctl_plus_disp_name = pv.split()[0]
    #                self.set_variable(fctl_plus_disp_name + ' TOTAL ENERGY', dft_energy)  # overwritten later for DH
    #else:
    #    self.set_variable(self.functional().name() + ' TOTAL ENERGY', dft_energy)  # overwritten later for DH

    self.set_variable("SCF ITERATIONS", self.iteration_)  # P::e SCF


def scf_print_preiterations(self,small=False):
    # small version does not print Nalpha,Nbeta,Ndocc,Nsocc, e.g. for SAD guess where they are not
    # available
    ct = self.molecule().point_group().char_table()

    if not small:
        core.print_out("   -------------------------------------------------------\n")
        core.print_out("    Irrep   Nso     Nmo     Nalpha   Nbeta   Ndocc  Nsocc\n")
        core.print_out("   -------------------------------------------------------\n")

        for h in range(self.nirrep()):
            core.print_out(
                f"     {ct.gamma(h).symbol():<3s}   {self.nsopi()[h]:6d}  {self.nmopi()[h]:6d}  {self.nalphapi()[h]:6d}  {self.nbetapi()[h]:6d}  {self.doccpi()[h]:6d}  {self.soccpi()[h]:6d}\n"
            )

        core.print_out("   -------------------------------------------------------\n")
        core.print_out(
            f"    Total  {self.nso():6d}  {self.nmo():6d}  {self.nalpha():6d}  {self.nbeta():6d}  {self.nbeta():6d}  {self.nalpha() - self.nbeta():6d}\n"
        )
        core.print_out("   -------------------------------------------------------\n\n")
    else:
        core.print_out("   -------------------------\n")
        core.print_out("    Irrep   Nso     Nmo    \n")
        core.print_out("   -------------------------\n")

        for h in range(self.nirrep()):
            core.print_out(
                f"     {ct.gamma(h).symbol():<3s}   {self.nsopi()[h]:6d}  {self.nmopi()[h]:6d} \n"
            )

        core.print_out("   -------------------------\n")
        core.print_out(
            f"    Total  {self.nso():6d}  {self.nmo():6d}\n"
        )
        core.print_out("   -------------------------\n\n")


# Bind functions to core.HF class
core.HF.initialize = scf_initialize
core.HF.initialize_jk = initialize_jk
core.HF.iterations = scf_iterate
core.HF.compute_energy = scf_compute_energy
core.HF.finalize_energy = scf_finalize_energy
core.HF.print_energies = scf_print_energies
core.HF.print_preiterations = scf_print_preiterations
core.HF.iteration_energies = []


def _converged(e_delta, d_rms, e_conv=None, d_conv=None):
    if e_conv is None:
        e_conv = core.get_option("SCF", "E_CONVERGENCE")
    if d_conv is None:
        d_conv = core.get_option("SCF", "D_CONVERGENCE")

    return (abs(e_delta) < e_conv and d_rms < d_conv)


def _validate_damping():
    """Sanity-checks DAMPING control options

    Raises
    ------
    ValidationError
        If any of |scf__damping_percentage|, |scf__damping_convergence|
        don't play well together.

    Returns
    -------
    bool
        Whether DAMPING is enabled during scf.

    """
    # Q: I changed the enabled criterion get_option <-- has_option_changed
    enabled = (core.get_option('SCF', 'DAMPING_PERCENTAGE') > 0.0)
    if enabled:
        parameter = core.get_option('SCF', "DAMPING_PERCENTAGE")
        if parameter < 0.0 or parameter > 100.0:
            raise ValidationError('SCF DAMPING_PERCENTAGE ({}) must be between 0 and 100'.format(parameter))

        stop = core.get_option('SCF', 'DAMPING_CONVERGENCE')
        if stop < 0.0:
            raise ValidationError('SCF DAMPING_CONVERGENCE ({}) must be > 0'.format(stop))

    return enabled


def _validate_diis(self):
    """Sanity-checks DIIS control options

    Raises
    ------
    psi4.driver.p4util.exceptions.ValidationError
        If any of DIIS options don't play well together.

    Returns
    -------
    bool
        Whether some form of DIIS is enabled during SCF.

    """

    restricted_open = self.same_a_b_orbs() and not self.same_a_b_dens()
    aediis_active = core.get_option('SCF', 'SCF_INITIAL_ACCELERATOR') != "NONE" and not restricted_open

    if aediis_active:
        start = core.get_option('SCF', 'SCF_INITIAL_START_DIIS_TRANSITION')
        stop = core.get_option('SCF', 'SCF_INITIAL_FINISH_DIIS_TRANSITION')
        if start < stop:
            raise ValidationError('SCF_INITIAL_START_DIIS_TRANSITION error magnitude cannot be less than SCF_INITIAL_FINISH_DIIS_TRANSITION.')
        elif start < 0:
            raise ValidationError('SCF_INITIAL_START_DIIS_TRANSITION cannot be negative.')
        elif stop < 0:
            raise ValidationError('SCF_INITIAL_FINISH_DIIS_TRANSITION cannot be negative.')

    enabled = bool(core.get_option('SCF', 'DIIS')) or aediis_active
    if enabled:
        start = core.get_option('SCF', 'DIIS_START')
        if start < 1:
            raise ValidationError('SCF DIIS_START ({}) must be at least 1'.format(start))

    return enabled


def _validate_frac():
    """Sanity-checks FRAC control options

    Raises
    ------
    ValidationError
        If any of |scf__frac_start| don't play well together.

    Returns
    -------
    bool
        Whether FRAC is enabled during scf.

    """
    enabled = (core.get_option('SCF', 'FRAC_START') != 0)
    if enabled:
        if enabled < 0:
            raise ValidationError('SCF FRAC_START ({}) must be at least 1'.format(enabled))

    return enabled


def _validate_MOM():
    """Sanity-checks MOM control options

    Raises
    ------
    ValidationError
        If any of |scf__mom_start|, |scf__mom_occ| don't play well together.

    Returns
    -------
    bool
        Whether excited-state MOM (not just the plain stabilizing MOM) is enabled during scf.

    """
    enabled = (core.get_option('SCF', "MOM_START") != 0 and len(core.get_option('SCF', "MOM_OCC")) > 0)
    if enabled:
        start = core.get_option('SCF', "MOM_START")
        if enabled < 0:
            raise ValidationError('SCF MOM_START ({}) must be at least 1'.format(start))

    return enabled


def _validate_soscf():
    """Sanity-checks SOSCF control options

    Raises
    ------
    ValidationError
        If any of |scf__soscf|, |scf__soscf_start_convergence|,
        |scf__soscf_min_iter|, |scf__soscf_max_iter| don't play well together.

    Returns
    -------
    bool
        Whether SOSCF is enabled during scf.

    """
    enabled = core.get_option('SCF', 'SOSCF')
    if enabled:
        start = core.get_option('SCF', 'SOSCF_START_CONVERGENCE')
        if start < 0.0:
            raise ValidationError('SCF SOSCF_START_CONVERGENCE ({}) must be positive'.format(start))

        miniter = core.get_option('SCF', 'SOSCF_MIN_ITER')
        if miniter < 1:
            raise ValidationError('SCF SOSCF_MIN_ITER ({}) must be at least 1'.format(miniter))

        maxiter = core.get_option('SCF', 'SOSCF_MAX_ITER')
        if maxiter < miniter:
            raise ValidationError('SCF SOSCF_MAX_ITER ({}) must be at least SOSCF_MIN_ITER ({})'.format(
                maxiter, miniter))

        conv = core.get_option('SCF', 'SOSCF_CONV')
        if conv < 1.e-10:
            raise ValidationError('SCF SOSCF_CONV ({}) must be achievable'.format(conv))

    return enabled

core.HF.validate_diis = _validate_diis

def efp_field_fn(xyz):
    """Callback function for PylibEFP to compute electric field from electrons
    in ab initio part for libefp polarization calculation.

    Parameters
    ----------
    xyz : list
        (3 * npt, ) flat array of points at which to compute electric field

    Returns
    -------
    list
        (3 * npt, ) flat array of electric field at points in `xyz`.

    Notes
    -----
    Function signature defined by libefp, so function uses number of
    basis functions and integrals factory `mints_psi4_yo` and total density
    matrix `efp_Dt_psi4_yo` from global namespace.

    """
    points = core.Matrix.from_array(np.array(xyz).reshape(-1, 3))
    field = mints_psi4_yo.electric_field_value(points, efp_Dt_psi4_yo).np.flatten()
    return field
