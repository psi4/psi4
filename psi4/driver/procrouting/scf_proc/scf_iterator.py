"""
The SCF iteration functions
"""

from psi4.driver import p4util
from psi4.driver import constants
from psi4.driver.p4util.exceptions import ConvergenceError, ValidationError
from psi4 import core

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
    if core.get_option('SCF', 'DF_SCF_GUESS') and (core.get_option('SCF', 'SCF_TYPE') == 'DIRECT'):
        # speed up DIRECT algorithm (recomputes full (non-DF) integrals
        #   each iter) by first converging via fast DF iterations, then
        #   fully converging in fewer slow DIRECT iterations. aka Andy trick 2.0
        core.print_out("  Starting with a DF guess...\n\n")
        with p4util.OptionsStateCM(['SCF', 'SCF_TYPE']):
            core.set_local_option('SCF', 'SCF_TYPE', 'DF')
            self.initialize()
            try:
                self.iterations()
            except ConvergenceError:
                raise ConvergenceError("""SCF DF preiterations""", self.iteration_)
        core.print_out("\n  DF guess converged.\n\n")

        # reset the DIIS & JK objects in prep for DIRECT
        if self.initialized_diis_manager_:
            self.diis_manager().reset_subspace()
        self.integrals()
    else:
        self.initialize()

    try:
        self.iterations()
    except ConvergenceError as e:
        if core.get_option("SCF", "FAIL_ON_MAXITER"):
            core.print_out("  Failed to converge.\n")
            # energy = 0.0
            # A P::e fn to either throw or protest upon nonconvergence
            # die_if_not_converged()
            raise e
        else:
            core.print_out("  Energy did not converge, but proceeding anyway.\n\n")
    else:
        core.print_out("  Energy converged.\n\n")

    scf_energy = self.finalize_energy()
    return scf_energy


def scf_initialize(self):
    """Specialized initialization, compute integrals and does everything to prepare for iterations"""

    self.iteration_ = 0

    if core.get_option('SCF', "PRINT") > 0:
        core.print_out("  ==> Pre-Iterations <==\n\n")
        self.print_preiterations()

    if self.attempt_number_ == 1:
        mints = core.MintsHelper(self.basisset())
        if core.get_global_option('RELATIVISTIC') in ['X2C', 'DKH']:
            mints.set_rel_basisset(self.get_basisset('BASIS_RELATIVISTIC'))

        mints.one_electron_integrals()
        self.integrals()

        core.timer_on("HF: Form core H")
        self.form_H()
        core.timer_off("HF: Form core H")

        # #ifdef USING_libefp
        # // EFP: Add in permanent moment contribution and cache
        # if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
        #     std::shared_ptr<Matrix> Vefp = Process::environment.get_efp()->modify_Fock_permanent();
        #     H_->add(Vefp);
        #     Horig_ = SharedMatrix(new Matrix("H orig Matrix", basisset_->nbf(), basisset_->nbf()));
        #     Horig_->copy(H_);
        #     outfile->Printf( "  QM/EFP: iterating Total Energy including QM/EFP Induction\n");
        # }
        # #endif

        core.timer_on("HF: Form S/X")
        self.form_Shalf()
        core.timer_off("HF: Form S/X")

        core.timer_on("HF: Guess")
        self.guess()
        core.timer_off("HF: Guess")

    else:
        # We're reading the orbitals from the previous set of iterations.
        self.form_D()
        self.set_energies("Total Energy", self.compute_initial_E())

    # #ifdef USING_libefp
    # if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
    #     Process::environment.get_efp()->set_qm_atoms();
    # }
    # #endif


def scf_iterate(self, e_conv=None, d_conv=None):

    df = core.get_option('SCF', "SCF_TYPE") == "DF"
    verbose = core.get_option('SCF', "PRINT")
    reference = core.get_option('SCF', "REFERENCE")

    # self.member_data_ signals are non-local, used internally by c-side fns
    self.diis_enabled_ = _validate_diis()
    self.MOM_excited_ = _validate_MOM()
    self.diis_start_ = core.get_option('SCF', 'DIIS_START')
    self.pcm_enabled_ = core.get_option('SCF', 'PCM')
    damping_enabled = _validate_damping()
    soscf_enabled = _validate_soscf()
    frac_enabled = _validate_frac()

    if self.iteration_ < 2:
        core.print_out("  ==> Iterations <==\n\n")
        core.print_out("%s                        Total Energy        Delta E     RMS |[F,P]|\n\n" % ("   "
                                                                                                      if df else ""))

    # SCF iterations!
    SCFE_old = 0.0
    SCFE = 0.0
    Drms = 0.0
    while True:
        self.iteration_ += 1

        diis_performed = False
        soscf_performed = False
        self.frac_performed_ = False
        #self.MOM_performed_ = False  # redundant from common_init()

        self.save_density_and_energy()

        # #ifdef USING_libefp
        # Horig = self.H().clone()
        # self.H().copy(Horig)
        # self.H().axpy(1.0, Vefp)

        #         # add efp contribution to Fock matrix
        #         if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
        #             H_->copy(Horig_)
        #             std::shared_ptr<Matrix> Vefp = Process::environment.get_efp()->modify_Fock_induced()
        #             H_->add(Vefp)
        #         }
        # #endif

        SCFE = 0.0
        self.clear_external_potentials()

        core.timer_on("HF: Form G")
        self.form_G()
        core.timer_off("HF: Form G")

        # reset fractional SAD occupation
        if (self.iteration_ == 0) and self.reset_occ_:
            self.reset_occupation()

        upcm = 0.0
        if self.pcm_enabled_:
            calc_type = core.PCM.CalcType.Total
            if core.get_option("PCM", "PCM_SCF_TYPE") is "SEPARATE":
                calc_type = core.PCM.CalcType.NucAndEle
            Dt = self.Da().clone()
            Dt.add(self.Db())
            upcm, Vpcm = self.get_PCM().compute_PCM_terms(Dt, calc_type)
            SCFE += upcm
            self.set_energies("PCM Polarization", upcm)
            self.set_variable("PCM POLARIZATION ENERGY", upcm)
            self.push_back_external_potential(Vpcm)
        else:
            self.set_energies("PCM polarization", upcm)
            self.set_variable("PCM POLARIZATION ENERGY", upcm)

        core.timer_on("HF: Form F")
        self.form_F()
        core.timer_off("HF: Form F")

        if verbose > 3:
            self.Fa().print_out()
            self.Fb().print_out()

        SCFE += self.compute_E()

        # #ifdef USING_libefp
        #        # add efp contribution to energy
        #        if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
        #            double efp_wfn_dependent_energy = Process::environment.get_efp()->scf_energy_update()
        #            E_ += efp_wfn_dependent_energy
        #        }
        # #endif

        self.set_energies("Total Energy", SCFE)
        Ediff = SCFE - SCFE_old
        SCFE_old = SCFE

        status = []

        # We either do SOSCF or DIIS
        if (soscf_enabled and (self.iteration_ > 3) and (Drms < core.get_option('SCF', 'SOSCF_START_CONVERGENCE'))):

            Drms = self.compute_orbital_gradient(False, core.get_option('SCF', 'DIIS_MAX_VECS'))
            diis_performed = False
            if self.functional().needs_xc():
                base_name = "SOKS, nmicro="
            else:
                base_name = "SOSCF, nmicro="

            if not _converged(Ediff, Drms, e_conv=e_conv, d_conv=d_conv):
                nmicro = self.soscf_update(
                    core.get_option('SCF', 'SOSCF_CONV'),
                    core.get_option('SCF', 'SOSCF_MIN_ITER'),
                    core.get_option('SCF', 'SOSCF_MAX_ITER'), core.get_option('SCF', 'SOSCF_PRINT'))
                if nmicro > 0:
                    # if zero, the soscf call bounced for some reason
                    self.find_occupation()
                    status.append(base_name + str(nmicro))
                    soscf_performed = True  # Stops DIIS
                else:
                    if verbose > 0:
                        core.print_out("Did not take a SOSCF step, using normal convergence methods\n")
                    soscf_performed = False  # Back to DIIS

            else:
                # need to ensure orthogonal orbitals and set epsilon
                status.append(base_name + "conv")
                core.timer_on("HF: Form C")
                self.form_C()
                core.timer_off("HF: Form C")
                soscf_performed = True  # Stops DIIS

        if not soscf_performed:
            # Normal convergence procedures if we do not do SOSCF

            core.timer_on("HF: DIIS")
            diis_performed = False
            add_to_diis_subspace = False

            if self.diis_enabled_ and self.iteration_ >= self.diis_start_:
                add_to_diis_subspace = True

            Drms = self.compute_orbital_gradient(add_to_diis_subspace, core.get_option('SCF', 'DIIS_MAX_VECS'))

            if (self.diis_enabled_
                    and self.iteration_ >= self.diis_start_ + core.get_option('SCF', 'DIIS_MIN_VECS') - 1):
                diis_performed = self.diis()

            if diis_performed:
                status.append("DIIS")

            core.timer_off("HF: DIIS")

            if verbose > 4 and diis_performed:
                core.print_out("  After DIIS:\n")
                self.Fa().print_out()
                self.Fb().print_out()

            # frac, MOM invoked here from Wfn::HF::find_occupation
            core.timer_on("HF: Form C")
            self.form_C()
            core.timer_off("HF: Form C")

        if self.MOM_performed_:
            status.append("MOM")

        if self.frac_performed_:
            status.append("FRAC")

        core.timer_on("HF: Form D")
        self.form_D()
        core.timer_off("HF: Form D")

        core.set_variable("SCF ITERATION ENERGY", SCFE)

        # After we've built the new D, damp the update
        if (damping_enabled and self.iteration_ > 1 and Drms > core.get_option('SCF', 'DAMPING_CONVERGENCE')):
            damping_percentage = core.get_option('SCF', "DAMPING_PERCENTAGE")
            self.damping_update(damping_percentage * 0.01)
            status.append("DAMP={}%".format(round(damping_percentage)))

        if verbose > 3:
            self.Ca().print_out()
            self.Cb().print_out()
            self.Da().print_out()
            self.Db().print_out()

        # Print out the iteration
        core.print_out("   @%s%s iter %3d: %20.14f   %12.5e   %-11.5e %s\n" %
                       ("DF-" if df else "", reference, self.iteration_, SCFE, Ediff, Drms, '/'.join(status)))

        # if a an excited MOM is requested but not started, don't stop yet
        if self.MOM_excited_ and not self.MOM_performed_:
            continue

        # if a fractional occupation is requested but not started, don't stop yet
        if frac_enabled and not self.frac_performed_:
            continue

        # Call any postiteration callbacks

        if _converged(Ediff, Drms, e_conv=e_conv, d_conv=d_conv):
            break
        if self.iteration_ >= core.get_option('SCF', 'MAXITER'):
            raise ConvergenceError("""SCF iterations""", self.iteration_)


def scf_finalize_energy(self):
    """Performs stability analysis and calls back SCF with new guess
    if needed, Returns the SCF energy. This function should be called
    once orbitals are ready for energy/property computations, usually
    after iterations() is called.

    """
    # Perform wavefunction stability analysis before doing
    # anything on a wavefunction that may not be truly converged.
    if core.get_option('SCF', 'STABILITY_ANALYSIS') != "NONE":
        # We need the integral file, make sure it is written and
        # compute it if needed
        if core.get_option('SCF', 'REFERENCE') != "UHF":
            psio = core.IO.shared_object()
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
                self.diis_manager().reset_subspace()
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

    # #ifdef USING_libefp
    #     if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
    #         Process::environment.get_efp()->compute()

    #         double efp_wfn_independent_energy = Process::environment.globals["EFP TOTAL ENERGY"] -
    #                                             Process::environment.globals["EFP IND ENERGY"]
    #         energies_["EFP"] = Process::environment.globals["EFP TOTAL ENERGY"]

    #         core.print_out("    EFP excluding EFP Induction   %20.12f [Eh]\n", efp_wfn_independent_energy)
    #         core.print_out("    SCF including EFP Induction   %20.12f [Eh]\n", E_)

    #         E_ += efp_wfn_independent_energy

    #         core.print_out("    Total SCF including Total EFP %20.12f [Eh]\n", E_)
    #     }
    # #endif

    core.print_out("\n  ==> Post-Iterations <==\n\n")

    self.check_phases()
    self.compute_spin_contamination()
    self.frac_renormalize()
    reference = core.get_option("SCF", "REFERENCE")

    energy = self.get_energies("Total Energy")

    if core.get_option('SCF', 'PRINT') > 0:
        self.print_orbitals()

    prefix = ""
    if core.get_option("SCF", "SCF_TYPE") == "DF":
        prefix = "DF-"

    core.print_out("  @%s%s Final Energy: %20.14f" % (prefix, reference, energy))
    # if (perturb_h_) {
    #     core.print_out(" with %f %f %f perturbation" %
    #                    (dipole_field_strength_[0], dipole_field_strength_[1], dipole_field_strength_[2]))
    # }
    core.print_out("\n\n")
    self.print_energies()

    self.clear_external_potentials()
    if self.pcm_enabled_:
        calc_type = core.PCM.CalcType.Total
        if core.get_option("PCM", "PCM_SCF_TYPE") == "SEPARATE":
            calc_type = core.PCM.CalcType.NucAndEle
        Dt = self.Da().clone()
        Dt.add(self.Db())
        _, Vpcm = self.get_PCM().compute_PCM_terms(Dt, calc_type)
        self.push_back_external_potential(Vpcm)

    # recompute the Fock matrices as they are modified during the SCF
    #   iteration and might need to be dumped to checkpoint later
    self.form_F()

    # Properties
    #  Comments so that autodoc utility will find these PSI variables
    #  Process::environment.globals["SCF DIPOLE X"] =
    #  Process::environment.globals["SCF DIPOLE Y"] =
    #  Process::environment.globals["SCF DIPOLE Z"] =
    #  Process::environment.globals["SCF QUADRUPOLE XX"] =
    #  Process::environment.globals["SCF QUADRUPOLE XY"] =
    #  Process::environment.globals["SCF QUADRUPOLE XZ"] =
    #  Process::environment.globals["SCF QUADRUPOLE YY"] =
    #  Process::environment.globals["SCF QUADRUPOLE YZ"] =
    #  Process::environment.globals["SCF QUADRUPOLE ZZ"] =

    # Orbitals are always saved, in case an MO guess is requested later
    # save_orbitals()

    # Shove variables into global space
    for k, v in self.variables().items():
        core.set_variable(k, v)

    self.finalize()

    core.print_out("\nComputation Completed\n")

    return energy


def scf_print_energies(self):
    enuc = self.get_energies('Nuclear')
    e1 = self.get_energies('One-Electron')
    e2 = self.get_energies('Two-Electron')
    exc = self.get_energies('XC')
    ed = self.get_energies('-D')
    evv10 = self.get_energies('VV10')
    eefp = self.get_energies('EFP')
    epcm = self.get_energies('PCM Polarization')

    hf_energy = enuc + e1 + e2
    dft_energy = hf_energy + exc + ed + evv10
    total_energy = dft_energy + eefp + epcm

    core.print_out("   => Energetics <=\n\n")
    core.print_out("    Nuclear Repulsion Energy =        {:24.16f}\n".format(enuc))
    core.print_out("    One-Electron Energy =             {:24.16f}\n".format(e1))
    core.print_out("    Two-Electron Energy =             {:24.16f}\n".format(e2))
    if self.functional().needs_xc():
        core.print_out("    DFT Exchange-Correlation Energy = {:24.16f}\n".format(exc))
        core.print_out("    Empirical Dispersion Energy =     {:24.16f}\n".format(ed))
        core.print_out("    VV10 Nonlocal Energy =            {:24.16f}\n".format(evv10))
    if self.pcm_enabled_:
        core.print_out("    PCM Polarization Energy =         {:24.16f}\n".format(epcm))
    #if (Process::environment.get_efp()->get_frag_count() > 0) {
    #    outfile->Printf("    EFP Energy =                      %24.16f\n", energies_["EFP"]);
    #}
    core.print_out("    Total Energy =                    {:24.16f}\n".format(total_energy))

    self.set_variable('NUCLEAR REPULSION ENERGY', enuc)
    self.set_variable('ONE-ELECTRON ENERGY', e1)
    self.set_variable('TWO-ELECTRON ENERGY', e2)
    if abs(exc) > 1.0e-14:
        self.set_variable('DFT XC ENERGY', exc)
        self.set_variable('DFT VV10 ENERGY', evv10)
        self.set_variable('DFT FUNCTIONAL TOTAL ENERGY', hf_energy + exc + evv10)
        self.set_variable('DFT TOTAL ENERGY', dft_energy)
    else:
        self.set_variable('HF TOTAL ENERGY', hf_energy)
    if abs(ed) > 1.0e-14:
        self.set_variable('DISPERSION CORRECTION ENERGY', ed)

    self.set_variable('SCF ITERATIONS', self.iteration_)
    self.set_variable('SCF N ITERS', self.iteration_)


core.HF.initialize = scf_initialize
core.HF.iterations = scf_iterate
core.HF.compute_energy = scf_compute_energy
core.HF.finalize_energy = scf_finalize_energy
core.HF.print_energies = scf_print_energies


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


def _validate_diis():
    """Sanity-checks DIIS control options

    Raises
    ------
    ValidationError
        If any of |scf__diis|, |scf__diis_start|,
        |scf__diis_min_vecs|, |scf__diis_max_vecs| don't play well together.

    Returns
    -------
    bool
        Whether DIIS is enabled during scf.

    """
    enabled = bool(core.get_option('SCF', 'DIIS'))
    if enabled:
        start = core.get_option('SCF', 'DIIS_START')
        if start < 1:
            raise ValidationError('SCF DIIS_START ({}) must be at least 1'.format(start))

        minvecs = core.get_option('SCF', 'DIIS_MIN_VECS')
        if minvecs < 1:
            raise ValidationError('SCF DIIS_MIN_VECS ({}) must be at least 1'.format(minvecs))

        maxvecs = core.get_option('SCF', 'DIIS_MAX_VECS')
        if maxvecs < minvecs:
            raise ValidationError(
                'SCF DIIS_MAX_VECS ({}) must be at least DIIS_MIN_VECS ({})'.format(maxvecs, minvecs))

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
            raise ValidationError(
                'SCF SOSCF_MAX_ITER ({}) must be at least SOSCF_MIN_ITER ({})'.format(maxiter, miniter))

        conv = core.get_option('SCF', 'SOSCF_CONV')
        if conv < 1.e-10:
            raise ValidationError('SCF SOSCF_CONV ({}) must be achievable'.format(conv))

    return enabled
