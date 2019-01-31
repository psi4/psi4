#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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

from psi4.driver import p4util
from psi4.driver import constants
from psi4.driver.p4util.exceptions import SCFConvergenceError, ValidationError
from psi4 import core

from .efp import get_qm_atoms_opts, modify_Fock_permanent, modify_Fock_induced

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
            self.diis_manager().reset_subspace()
        self.initialize_jk(self.memory_jk_)
    else:
        self.initialize()

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
    jk = core.JK.build(
        wfn.get_basisset("ORBITAL"),
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
            collocation_size *= 4 # First derivs
        elif vbase.functional().ansatz() == 2:
            collocation_size *= 10 # Second derivs
    else:
        collocation_size = 0

    # Change allocation for collocation matrices based on DFT type
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

    # Print out initial docc/socc/etc data
    if self.get_print():
        core.print_out("  ==> Pre-Iterations <==\n\n")
        self.print_preiterations()

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

        efpobj.set_electron_density_field_fn(field_fn)

    # Initilize all integratals and perform the first guess
    if self.attempt_number_ == 1:
        mints = core.MintsHelper(self.basisset())
        if core.get_global_option('RELATIVISTIC') in ['X2C', 'DKH']:
            mints.set_rel_basisset(self.get_basisset('BASIS_RELATIVISTIC'))

        mints.one_electron_integrals()
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
            Vefp = modify_Fock_permanent(self.molecule(), mints, verbose=verbose-1)
            Vefp = core.Matrix.from_array(Vefp)
            self.H().add(Vefp)
            Horig = self.H().clone()
            self.Horig = Horig
            core.print_out("  QM/EFP: iterating Total Energy including QM/EFP Induction\n")
            core.timer_off("HF: Form Vefp")

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


    # turn off VV10 for iterations
    if core.get_option('SCF', "DFT_VV10_POSTSCF") and self.functional().vv10_b() > 0.0:
        core.print_out("  VV10: post-SCF option active \n \n")
        self.functional().set_lock(False)
        self.functional().set_do_vv10(False)
        self.functional().set_lock(True)



def scf_iterate(self, e_conv=None, d_conv=None):

    is_dfjk = core.get_global_option('SCF_TYPE').endswith('DF')
    verbose = core.get_option('SCF', "PRINT")
    reference = core.get_option('SCF', "REFERENCE")

    # self.member_data_ signals are non-local, used internally by c-side fns
    self.diis_enabled_ = _validate_diis()
    self.MOM_excited_ = _validate_MOM()
    self.diis_start_ = core.get_option('SCF', 'DIIS_START')
    damping_enabled = _validate_damping()
    soscf_enabled = _validate_soscf()
    frac_enabled = _validate_frac()
    efp_enabled = hasattr(self.molecule(), 'EFP')
    diis_rms = core.get_option('SCF', 'DIIS_RMS_ERROR')

    if self.iteration_ < 2:
        core.print_out("  ==> Iterations <==\n\n")
        core.print_out("%s                        Total Energy        Delta E     %s |[F,P]|\n\n" % ("   "
                                                                                                     if is_dfjk else "", "RMS" if diis_rms else "MAX"))

    # SCF iterations!
    SCFE_old = 0.0
    SCFE = 0.0
    Dnorm = 0.0
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
            Vefp = modify_Fock_induced(self.molecule().EFP, mints_psi4_yo, verbose=verbose-1)
            Vefp = core.Matrix.from_array(Vefp)
            self.H().add(Vefp)

        SCFE = 0.0
        self.clear_external_potentials()

        core.timer_on("HF: Form G")
        self.form_G()
        core.timer_off("HF: Form G")

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
        self.set_variable("PCM POLARIZATION ENERGY", upcm)
        self.set_energies("PCM Polarization", upcm)

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
                nmicro = self.soscf_update(
                    core.get_option('SCF', 'SOSCF_CONV'),
                    core.get_option('SCF', 'SOSCF_MIN_ITER'),
                    core.get_option('SCF', 'SOSCF_MAX_ITER'), core.get_option('SCF', 'SOSCF_PRINT'))
                # if zero, the soscf call bounced for some reason
                soscf_performed = (nmicro>0)

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
            # reset the occupations. From here on, the density matrices
            # are correct.
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

                if (add_to_diis_subspace and core.get_option('SCF', 'DIIS_MIN_VECS') - 1):
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

                # Reset occupations if necessary
                if (self.iteration_ == 0) and self.reset_occ_:
                    self.reset_occupation()
                    self.find_occupation()

        # Form new density matrix
        core.timer_on("HF: Form D")
        self.form_D()
        core.timer_off("HF: Form D")

        # After we've built the new D, damp the update
        if (damping_enabled and self.iteration_ > 1 and Dnorm > core.get_option('SCF', 'DAMPING_CONVERGENCE')):
            damping_percentage = core.get_option('SCF', "DAMPING_PERCENTAGE")
            self.damping_update(damping_percentage * 0.01)
            status.append("DAMP={}%".format(round(damping_percentage)))

        if verbose > 3:
            self.Ca().print_out()
            self.Cb().print_out()
            self.Da().print_out()
            self.Db().print_out()

        # Print out the iteration
        core.print_out("   @%s%s iter %3s: %20.14f   %12.5e   %-11.5e %s\n" %
                       ("DF-" if is_dfjk else "", reference, "SAD" if ((self.iteration_ == 0) and self.sad_) else self.iteration_, SCFE, Ediff, Dnorm, '/'.join(status)))

        # if a an excited MOM is requested but not started, don't stop yet
        if self.MOM_excited_ and not self.MOM_performed_:
            continue

        # if a fractional occupation is requested but not started, don't stop yet
        if frac_enabled and not self.frac_performed_:
            continue

        # Call any postiteration callbacks
        if not ((self.iteration_ == 0) and self.sad_) and _converged(Ediff, Dnorm, e_conv=e_conv, d_conv=d_conv):
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

        # Don't bother computing needed integrals if we can't do anything with them.
        if self.functional().needs_xc():
            raise ValidationError("Stability analysis not yet supported for XC functionals.")

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
            #next 2 lines fix a bug that prohibits relativistic stability analysis
            if core.get_global_option('RELATIVISTIC') in ['X2C', 'DKH']:
                mints.set_rel_basisset(self.get_basisset('BASIS_RELATIVISTIC'))
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

        core.set_variable('EFP ELST ENERGY', efpene['electrostatic'] + efpene['charge_penetration'] + efpene['electrostatic_point_charges'])
        core.set_variable('EFP IND ENERGY', efpene['polarization'])
        core.set_variable('EFP DISP ENERGY', efpene['dispersion'])
        core.set_variable('EFP EXCH ENERGY', efpene['exchange_repulsion'])
        core.set_variable('EFP TOTAL ENERGY', efpene['total'])
        core.set_variable('CURRENT ENERGY', efpene['total'])

    core.print_out("\n  ==> Post-Iterations <==\n\n")

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

    self.clear_external_potentials()
    if core.get_option('SCF', 'PCM'):
        calc_type = core.PCM.CalcType.Total
        if core.get_option("PCM", "PCM_SCF_TYPE") == "SEPARATE":
            calc_type = core.PCM.CalcType.NucAndEle
        Dt = self.Da().clone()
        Dt.add(self.Db())
        _, Vpcm = self.get_PCM().compute_PCM_terms(Dt, calc_type)
        self.push_back_external_potential(Vpcm)

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

    # TODO re-enable
    self.finalize()
    if self.V_potential():
        self.V_potential().clear_collocation_cache()

    core.print_out("\nComputation Completed\n")

    return energy


def scf_print_energies(self):
    enuc = self.get_energies('Nuclear')
    e1 = self.get_energies('One-Electron')
    e2 = self.get_energies('Two-Electron')
    exc = self.get_energies('XC')
    ed = self.get_energies('-D')
    #self.del_variable('-D Energy')
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
    if core.get_option('SCF', 'PCM'):
        core.print_out("    PCM Polarization Energy =         {:24.16f}\n".format(epcm))
    if hasattr(self.molecule(), 'EFP'):
        core.print_out("    EFP Energy =                      {:24.16f}\n".format(eefp))
    core.print_out("    Total Energy =                    {:24.16f}\n".format(total_energy))

    self.set_variable('NUCLEAR REPULSION ENERGY', enuc)
    self.set_variable('ONE-ELECTRON ENERGY', e1)
    self.set_variable('TWO-ELECTRON ENERGY', e2)
    if self.functional().needs_xc():
        self.set_variable('DFT XC ENERGY', exc)
        self.set_variable('DFT VV10 ENERGY', evv10)
        self.set_variable('DFT FUNCTIONAL TOTAL ENERGY', hf_energy + exc + evv10)
        #self.set_variable(self.functional().name() + ' FUNCTIONAL TOTAL ENERGY', hf_energy + exc + evv10)
        self.set_variable('DFT TOTAL ENERGY', dft_energy)  # overwritten later for DH
    else:
        self.set_variable('HF TOTAL ENERGY', hf_energy)
    if hasattr(self, "_disp_functor"):
        self.set_variable('DISPERSION CORRECTION ENERGY', ed)
    #if abs(ed) > 1.0e-14:
    #    for pv, pvv in self.variables().items():
    #        if abs(pvv - ed) < 1.0e-14:
    #            if pv.endswith('DISPERSION CORRECTION ENERGY') and pv.startswith(self.functional().name()):
    #                fctl_plus_disp_name = pv.split()[0]
    #                self.set_variable(fctl_plus_disp_name + ' TOTAL ENERGY', dft_energy)  # overwritten later for DH
    #else:
    #    self.set_variable(self.functional().name() + ' TOTAL ENERGY', dft_energy)  # overwritten later for DH

    self.set_variable('SCF ITERATIONS', self.iteration_)

# Bind functions to core.HF class
core.HF.initialize = scf_initialize
core.HF.initialize_jk = initialize_jk
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


def field_fn(xyz):
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
    points = np.array(xyz).reshape(-1, 3)
    npt = len(points)

    # Cartesian basis one-electron EFP perturbation
    nbf = mints_psi4_yo.basisset().nbf()
    field_ints = np.zeros((3, nbf, nbf))

    # Electric field at points
    field = np.zeros((npt, 3))

    for ipt in range(npt):
        # get electric field integrals from Psi4
        p4_field_ints = mints_psi4_yo.electric_field(origin=points[ipt])

        field[ipt] = [np.vdot(efp_Dt_psi4_yo, np.asarray(p4_field_ints[0])),  # Ex
                      np.vdot(efp_Dt_psi4_yo, np.asarray(p4_field_ints[1])),  # Ey
                      np.vdot(efp_Dt_psi4_yo, np.asarray(p4_field_ints[2]))]  # Ez

    field = np.reshape(field, 3 * npt)

    return field
