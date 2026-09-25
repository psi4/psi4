
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
The SCF iteration functions
"""
import os

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


def _release_collocation_cache(self):
    """Free the DFT collocation cache held by this wavefunction's V_potential.

    The cache is sized from the *whole* memory budget (see scf_initialize), so it
    must be released as soon as the SCF is done with it -- including when the SCF
    fails, since a non-converged SCF never reaches finalize_energy().
    """
    if self.V_potential():
        self.V_potential().clear_collocation_cache()


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
                _release_collocation_cache(self)
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
            _release_collocation_cache(self)
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


def _resident_doubles():
    """What this process is actually holding right now, in doubles, or None off Linux.

    ``memory_committed()`` only knows about the stores that declare themselves -- an
    in-core (Q|mn), a collocation cache -- and that is less than the process is holding.
    The rest is the holes glibc leaves behind: a released cache is tens of thousands of
    small matrices, so ``malloc_trim`` can hand back the free tops of the arenas but not
    a hole with a live allocation above it.  An SCF that divides up
    ``setting - committed`` therefore hands out memory the process has already spent;
    measured at 1.5 GiB by the last SCF of a SAPT(DFT) job, enough to put the peak over
    the declaration even though every printed grant fit.  VmRSS is what a cgroup kills
    on, so budget against that when we can read it.
    """
    try:
        with open("/proc/self/statm") as fh:
            resident_pages = int(fh.read().split()[1])
    except (OSError, IndexError, ValueError):
        return None
    return resident_pages * os.sysconf("SC_PAGE_SIZE") / 8


def _scf_memory_reserve(wfn):
    """Estimate, in doubles, what this SCF will spend *outside* the two budgets it hands out.

    ``scf_initialize`` divides the memory setting between the JK and the DFT collocation
    cache, but those are not the only two consumers.  The SCF matrices, the grid itself,
    the per-thread point-function scratch and the JK's own per-call transients are all
    paid for out of whatever the split happens to leave behind, and none of them is
    declared to the memory ledger.  That residue is what ``SCF_MEM_SAFETY_FACTOR`` has
    really been covering, and a fixed fraction of a user-chosen number is the wrong shape
    for it: the transient term grows as naux*nocc*nbf while the fraction stays put, so
    the factor is wasteful on small jobs and too small on exactly the large ones where
    running out of memory is most expensive.

    Every input below is known before a single integral is computed, so this can run
    before the budget is divided.  The coefficients come from per-stage high-water
    measurements on eight systems (nanotube and peptide dimers, cc-pVDZ through
    aug-cc-pVTZ, 4-32 threads); see ~/docs/saptdft-memory.html for the derivations.
    """
    nbf = wfn.basisset().nbf()
    nthread = core.get_num_threads()

    # Columns of C that the JK is handed: one set of occupieds if the densities are
    # constrained equal, alpha and beta separately otherwise.
    nocc = max(1, wfn.nalpha())
    if not wfn.same_a_b_dens():
        nocc += wfn.nbeta()

    # Interpreter, library statics and per-thread scratch that no module declares.  Flat
    # at 0.21 GiB across every probe, plus ~0.9 MB per thread of OpenMP/BLAS stacks.
    base = (0.21 * 1024**3 + 0.9e6 * nthread) / 8

    # The nbf^2 matrices an SCF holds: H, S, X, Fa/Fb, Ca/Cb, Da/Db, Va/Vb and friends,
    # plus two DIIS vectors (error and target) per subspace entry.
    matrices = float(nbf**2 * (12 + 2 * core.get_option("SCF", "DIIS_MAX_VECS")))

    # Grid storage and the point-function scratch.  Both are zero for a grid-less SCF.
    grid_mem = 0.0
    points_mem = 0.0
    vbase = wfn.V_potential()
    grid = vbase.grid() if vbase else None
    if grid:
        nblocks = len(grid.blocks())
        max_points = grid.max_points()
        max_functions = grid.max_functions()
        # BlockOPoints: x/y/z/w per point (32 B), plus a local-function map per block
        # (26 B per function held), plus ~4 MB/thread of libxc workspace.
        grid_mem = (32.0 * grid.npoints() + 26.0 * nblocks * max_functions + 4.0e6 * nthread) / 8
        # PointFunctions basis/point values, per thread.  The fit is for a GGA (phi and
        # three gradient components); an LDA needs a quarter of it, a meta-GGA 2.5x.
        ansatz_scale = {0: 0.25, 1: 1.0, 2: 2.5}.get(vbase.functional().ansatz(), 1.0)
        points_mem = ansatz_scale * nthread * (3.38 * max_points * max_functions +
                                               1.887 * max_functions**2)

    # compute_JK's per-call buffers.  DFHelper does size these against its own grant, but
    # the BLAS/libint working set that sits on top of the accounted T1/T2/C_buffers is
    # not: the measured ratio of real high-water to accounted buffers is 1.36-1.69 over
    # six MemDFJK configurations, hence the 1.5.  The naux*nocc*nbf term is the largest
    # Q-block case and so an upper bound -- deliberately, since this is a guard.
    scf_type = core.get_global_option("SCF_TYPE").upper()
    naux = 0
    if "DF" in scf_type and "+" not in scf_type:
        try:
            naux = wfn.get_basisset("DF_BASIS_SCF").nbf()
        except Exception:
            naux = 0
    jk_mem = 1.5 * (naux * nocc * nbf + nthread * nbf * max(nocc, nbf) + nbf**2)

    return {
        "base": base,
        "matrices": matrices,
        "grid": grid_mem,
        "points": points_mem,
        "jk": jk_mem,
    }


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

    # Get memory in terms of doubles.  The memory setting describes an empty process, but
    # an SCF is often started with earlier wavefunctions still alive -- SAPT(DFT) holds the
    # dimer and both monomers, and a GRAC shift holds the neutral while the cation runs --
    # whose JK integrals and collocation caches are already spending it.  Divide up what is
    # actually left instead of handing out the whole setting again.  Our own cache, if this
    # wavefunction has one from a previous SCF, would otherwise be double-counted against
    # us; drop it only on a first attempt, which is the only case that rebuilds it below.
    if self.attempt_number_ == 1:
        _release_collocation_cache(self)
    # Return whatever the caches we just dropped were still holding before measuring, so
    # that recoverable arena space is not mistaken for memory this SCF cannot have.
    core.release_freed_memory()
    committed_memory = core.memory_committed()
    resident_memory = _resident_doubles()
    unclaimed_memory = max(0.0, (core.get_memory() / 8) - committed_memory)

    # Set aside the undeclared-but-mandatory part before dividing up the rest.  Without
    # this the JK and the cache are handed a fraction of the whole setting and everything
    # else -- SCF matrices, grid, JK transients -- has to fit in the rounding error, which
    # is why a job whose printed JK estimate fits its allocation can still be killed.  The
    # reserve is capped at half the budget so that an under-declared run degrades into an
    # out-of-core JK (which DFHelper handles) rather than a zero grant.
    reserve_terms = _scf_memory_reserve(self)
    # The interpreter-and-libraries term is resident before psi4 allocates anything, so it
    # is not part of what the memory keyword describes -- that keyword has always sized the
    # big arrays, not the process.  Subtracting it from the budget would silently take a
    # fifth of a gigabyte away from every job that declares a small one, which is enough to
    # push a Cholesky or in-core DF JK under the floor it needs to start.  It is still real
    # memory, so it is still counted against the cache ceiling below and in the advice about
    # how much this SCF actually needs; it just is not taken out of the split.
    process_reserve = reserve_terms["base"]
    reserve_total = sum(reserve_terms.values()) - process_reserve
    if reserve_total >= unclaimed_memory:
        # The declaration is smaller than the parts of this SCF that are not negotiable, so
        # it is not a budget at all -- a test that asks for two megabytes to force a disk
        # algorithm, say.  Holding a reserve back cannot make such a job fit; it only pushes
        # the JK under the floor its out-of-core algorithm needs to start.  Hand over the
        # whole (fictional) budget and let the JK pick the smallest algorithm it has.
        reserve_memory = 0.0
    else:
        reserve_memory = min(reserve_total, 0.5 * unclaimed_memory)

    # What the process is resident for beyond both the stores that declared themselves and
    # the interpreter we just reserved for: the holes glibc leaves where a released cache
    # used to be, and anything else in this process that never announced itself.  It is not
    # subtracted from the budget -- the memory keyword has always described the big arrays
    # rather than the process, and a job that declares a megabyte to force an out-of-core
    # algorithm would be handed nothing at all -- it only caps the cache below.
    if resident_memory is None:
        overhead_memory = 0.0
    else:
        overhead_memory = max(0.0, resident_memory - committed_memory - process_reserve)
    held_memory = committed_memory + overhead_memory

    # With an absolute reserve subtracted, the factor is only residual slop, so the 0.75
    # default is far too conservative; honour it if the user set it, otherwise use 0.95.
    # When the reserve was not taken, nothing else is held back, and 0.95 hands an
    # out-of-core DiskDFJK -- which really does fill its grant -- the whole setting: the
    # SAPT(DFT) dHF dimer got 60.8 GiB of a 64 GiB job and was killed where the 0.75
    # default had given it 48 GiB and finished.  Keep the default there.
    if core.has_option_changed("SCF", "SCF_MEM_SAFETY_FACTOR"):
        safety_factor = core.get_option("SCF", "SCF_MEM_SAFETY_FACTOR")
    elif reserve_memory == 0.0 and reserve_total > 0.0:
        safety_factor = core.get_option("SCF", "SCF_MEM_SAFETY_FACTOR")
    else:
        safety_factor = 0.95
    total_memory = max(0.0, unclaimed_memory - reserve_memory) * safety_factor

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

    if initialize_jk_obj:
        # What this SCF is about to allocate for a JK of its own.  DiskJK, DirectJK and
        # CompositeJK all return 0 here because they cannot predict their footprint; that
        # means "unknown", not "nothing", so the remainder is not the cache's to take --
        # treating it that way is what lets an out-of-core SCF claim a *larger* collocation
        # cache than the DF one it was meant to be cheaper than.
        jk_size = jk.memory_estimate()
        jk_size_known = jk_size > 0
    else:
        # A re-used JK's integrals are already counted in committed_memory, but the buffers it
        # allocates on every build are not.  A MemDFJK on its disk algorithm (SCF_SUBTYPE
        # OUT_OF_CORE, or AUTO when the AOs do not fit) holds nothing in core and sizes its
        # (Q|mn) blocks and M/T/C buffers from its whole grant on every build, so treating it
        # as free hands that grant to the cache a second time: the second SCF of a GRAC pair,
        # and monomer B after monomer A, then cached the full grid on top of a 70 GiB JK and
        # were killed in their first iteration.  Charge what the JK may still allocate.  JKs
        # that cannot predict their footprint report 0 and keep the old treatment.
        reused_jk_predictable = jk.memory_estimate() > 0
        if reused_jk_predictable:
            jk_size = max(0, jk.memory() - jk.memory_held())
        else:
            jk_size = 0
        jk_size_known = True

    # Give remaining to collocation
    if jk_size_known and total_memory > jk_size:
        collocation_memory = total_memory - jk_size
    # Give up to 10% to collocation
    elif (total_memory * 0.1) > collocation_size:
        collocation_memory = collocation_size
    else:
        collocation_memory = total_memory * 0.1

    if collocation_memory > collocation_size:
        collocation_memory = collocation_size

    # Keep the cache under what is actually left of the declaration once the undeclared
    # residue is counted.  Everything else in this SCF is sized by the problem; the cache is
    # the only part that can be asked to take less, and it is what pushed the peak past the
    # declaration when fragmentation reached a GiB and a half by the last SCF of a SAPT(DFT)
    # job.  This is a ceiling rather than a subtraction on purpose: a job with headroom to
    # spare keeps its whole cache, and with no residue to account for the safety factor
    # already holds the split below this line, so nothing changes.  The JK keeps its budget.
    cache_ceiling = safety_factor * max(
        0.0, (core.get_memory() / 8) - held_memory - process_reserve - reserve_memory - jk_size)
    withheld_memory = max(0.0, collocation_memory - cache_ceiling)
    collocation_memory -= withheld_memory

    # Set constants
    self.iteration_ = 0
    self.memory_jk_ = int(total_memory - collocation_memory - withheld_memory)
    # When the cache takes everything above the JK's own estimate, the line above is
    # total_memory - (total_memory - jk_size), which is jk_size only up to the rounding
    # of two float subtractions -- and int() then truncates, so the JK can be handed one
    # double less than the size it just asked for.  DFHelper and CDJK both treat that as
    # "not enough memory to do in-core" and throw, which turns a grant that is exactly
    # right into a fatal error.  Never hand the JK less than it reported it needs -- but
    # only when that need fits the budget.  An out-of-core DiskDFJK reports its in-core
    # size, which can be larger than the whole setting; flooring at it then gave the JK
    # 108 GiB of a 96 GiB job, which it sized its blocks from and was killed.
    if jk_size_known and 0 < jk_size <= total_memory:
        self.memory_jk_ = max(self.memory_jk_, int(jk_size))
    self.memory_collocation_ = int(collocation_memory)

    if self.get_print():
        gib = lambda doubles: doubles * 8 / 1024**3

        core.print_out("  ==> SCF Memory <==\n\n")
        core.print_out("    Memory setting                  {:11.3f} [GiB]\n".format(core.get_memory() / 1024**3))
        if committed_memory:
            core.print_out("    Held by live JK / grid caches   {:11.3f} [GiB]\n".format(gib(committed_memory)))
        if overhead_memory > 0.01 * 1024**3 / 8:
            core.print_out("    Resident but unaccounted        {:11.3f} [GiB]".format(gib(overhead_memory)))
            if withheld_memory:
                core.print_out("  {:.3f} [GiB] withheld from the cache".format(gib(withheld_memory)))
            core.print_out("\n")
        core.print_out("    Reserve estimate                {:11.3f} [GiB]\n".format(gib(reserve_memory)))
        for label, key in (("SCF matrices", "matrices"), ("DFT grid", "grid"),
                           ("DFT point functions", "points"), ("JK per-call transients", "jk")):
            if reserve_terms[key]:
                core.print_out("      {:30s}{:11.3f}\n".format(label, gib(reserve_terms[key])))
        core.print_out("    Interpreter and libraries       {:11.3f} [GiB]  (not taken from the budget)\n".format(
            gib(process_reserve)))
        if reserve_total > reserve_memory:
            core.print_out("      (reserve {}; this job is under-declared)\n".format(
                "not taken" if reserve_memory == 0.0 else "capped at half the remaining budget"))
        if initialize_jk_obj:
            core.print_out("    JK allocation                   {:11.3f} [GiB]\n".format(gib(self.memory_jk_)))
        # A re-used JK keeps the grant it was built with and never sees memory_jk_, so show
        # what it may still allocate per build, which is what the split above charged for it.
        elif reused_jk_predictable:
            core.print_out("    Re-used JK working budget       {:11.3f} [GiB]  of {:.3f} [GiB] granted\n".format(
                gib(jk_size), gib(jk.memory())))
        else:
            core.print_out("    Re-used JK working budget           unknown  {} cannot predict its footprint\n".format(
                jk.name()))
        if collocation_size:
            core.print_out("    Collocation cache               {:11.3f} [GiB]  of {:.3f} [GiB] full\n".format(
                gib(self.memory_collocation_), gib(collocation_size)))

        if jk_size_known:
            # jk_size is what a JK built here will allocate; for a re-used one it is what
            # it may still allocate beyond the integrals already inside committed_memory.
            required = held_memory + process_reserve + reserve_memory + jk_size
            peak = required + self.memory_collocation_
            core.print_out("    Estimated peak                  {:11.3f} [GiB]\n".format(gib(peak)))
            core.print_out("    Minimum memory for this SCF     {:11.3f} [GiB]".format(gib(1.05 * required)))
            if collocation_size:
                core.print_out("  {:.3f} [GiB] to also cache the grid".format(
                    gib(1.05 * required + collocation_size)))
            core.print_out("\n")
        else:
            core.print_out("    Estimated peak                      unknown  {} cannot predict its footprint\n".format(
                core.get_global_option("SCF_TYPE")))
        core.print_out("\n")
        if held_memory:
            core.print_out("  The memory already held above belongs to wavefunctions that are still alive\n"
                           "  (a SAPT dimer and its monomers, or the neutral behind a GRAC shift); this SCF\n"
                           "  divides up the remainder rather than the whole setting.\n\n")

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
    ooo_scf = core.get_option("SCF", "ORBITAL_OPTIMIZER_PACKAGE") in ["OOO", "OPENORBITALOPTIMIZER"]
    if ooo_scf:
        pcm_enabled = core.get_option('SCF', 'PCM')
        ddx_enabled = core.get_option('SCF', 'DDX')
        pe_enabled = core.get_option('SCF', 'PE')
        level_shift_enabled = core.get_option("SCF", "LEVEL_SHIFT") != 0.0
        autograc_enabled = core.get_option("SAPT", "SAPT_DFT_GRAC_COMPUTE") != "NONE"
        guessmix_enabled = core.get_option("SCF", "GUESS_MIX")
        if (reference in ["ROHF", "CUHF"] or soscf_enabled or self.MOM_excited_ or frac_enabled or
            efp_enabled or pcm_enabled or ddx_enabled or pe_enabled or autograc_enabled or
            level_shift_enabled or guessmix_enabled):
            core.print_out(f"    Note: OpenOrbitalOptimizer not compatible with at least one of the following. Falling back to orbital_optimizer_package=internal\n")
            core.print_out(f"          {reference=}, soscf={soscf_enabled}, mom={self.MOM_excited_}, frac={frac_enabled}, efp={efp_enabled},\n")
            core.print_out(f"          pcm={pcm_enabled}, ddx={ddx_enabled}, pe={pe_enabled}, autograc={autograc_enabled}, level_shift={level_shift_enabled},\n")
            core.print_out(f"          guess_mix={guessmix_enabled}\n")
        else:
            # SAD needs some special work since the guess doesn't actually make the orbitals in Psi4
            if self.sad_ and self.iteration_ <= 0:
                self.iteration_ += 1
                self.form_G()
                self.form_initial_F()
                self.form_initial_C()
                self.reset_occupation()
                self.find_occupation()
                ene_sad = self.compute_E()
                core.print_out(
                    "   @%s%s iter %3s: %20.14f   %12.5e   %-11.5e %s\n" %
                    ("DF-" if is_dfjk else "", reference, "SAD", ene_sad, ene_sad, 0.0, ""))
            if core.get_option("SCF", "GUESS") == "READ" and self.iteration_ <= 0:
                self.form_G()
                self.form_initial_F()
                self.form_initial_C()
                self.reset_occupation()
                self.find_occupation()
                ene_sad = self.compute_E()

            try:
                self.openorbital_scf()
            except RuntimeError as ex:
                if "openorbital_scf is virtual; it has not been implemented for your class" in str(ex):
                    core.print_out(f"    Note: OpenOrbitalOptimizer NYI for {reference}. Falling back to Internal.\n")
                else:
                    raise ex
            else:
                SCFE = self.compute_E()
                self.set_energies("Total Energy", SCFE)
                self.set_variable("SCF ITERATION ENERGY", SCFE)
                self.iteration_energies.append(SCFE)  # note 1-len array, not niter-len array like INTERNAL

                self.form_G()
                self.form_F()
                self.form_C()
                self.form_D()
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
    _release_collocation_cache(self)

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
