 /*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"

#include "psi4/libfock/jk.h"
#include "psi4/libfock/soscf.h"

#include "psi4/detci/ciwave.h"
#include "psi4/detci/civect.h"

#include "psi4/libscf_solver/hf.h"
#include "psi4/libscf_solver/rhf.h"
#include "psi4/libscf_solver/uhf.h"
#include "psi4/libscf_solver/rohf.h"
#include "psi4/libscf_solver/cuhf.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libfock/v.h"

#include "psi4/dfep2/dfep2.h"

#include "psi4/fisapt/fisapt.h"

#include "psi4/libdiis/diismanager.h"

#include <string>

using namespace psi;
void export_wavefunction(py::module& m) {

    py::class_<scf::HF, std::shared_ptr<scf::HF>, Wavefunction>(m, "HF", "docstring")
        .def("form_C", &scf::HF::form_C, "Forms the Orbital Matrices from the current Fock Matrices.")
        .def("form_D", &scf::HF::form_D, "Forms the Density Matrices from the current Orbitals Matrices")
        .def("form_V", &scf::HF::form_V, "Form the Kohn-Sham Potential Matrices from the current Density Matrices")
        .def("form_G", &scf::HF::form_G, "Forms the G matrix.")
        .def("form_F", &scf::HF::form_F, "Forms the F matrix.")
        .def("form_H", &scf::HF::form_H, "Forms the core Hamiltonian")
        .def("form_Shalf", &scf::HF::form_Shalf, "Forms the S^1/2 matrix")
        .def("guess", &scf::HF::guess, "Forms the guess (guarantees C, D, and E)")
        .def("integrals", &scf::HF::integrals, "Sets up the JK object")
        .def("onel_Hx", &scf::HF::onel_Hx, "One-electron Hessian-vector products.")
        .def("twoel_Hx", &scf::HF::twoel_Hx, "Two-electron Hessian-vector products")
        .def("cphf_Hx", &scf::HF::cphf_Hx, "CPHF Hessian-vector prodcuts (4 * J - K - K.T).")
        .def("cphf_solve", &scf::HF::cphf_solve, py::arg("x_vec"), py::arg("conv_tol"), py::arg("max_iter"),
             py::arg("print_lvl") = 2, "Solves the CPHF equations for a given set of x vectors.")
        .def("cphf_converged", &scf::HF::cphf_converged, "Adds occupied guess alpha orbitals.")
        .def("guess_Ca", &scf::HF::guess_Ca, "Sets the guess Alpha Orbital Matrix")
        .def("guess_Cb", &scf::HF::guess_Cb, "Sets the guess Beta Orbital Matrix")
        .def_property("reset_occ_", &scf::HF::reset_occ, &scf::HF::set_reset_occ,
                      "Do reset the occupation after the guess to the inital occupation.")
        .def("set_sad_basissets", &scf::HF::set_sad_basissets, "Sets the Superposition of Atomic Densities basisset.")
        .def("set_sad_fitting_basissets", &scf::HF::set_sad_fitting_basissets,
             "Sets the Superposition of Atomic Densities density-fitted basisset.")
        .def("Va", &scf::HF::Va, "Returns the Alpha Kohn-Sham Potential Matrix.")
        .def("Vb", &scf::HF::Vb, "Returns the Beta Kohn-Sham Potential Matrix.")
        .def("jk", &scf::HF::jk, "Returns the internal JK object.")
        .def("set_jk", &scf::HF::set_jk, "Sets the internal JK object !expert.")
        .def("functional", &scf::HF::functional, "Returns the internal DFT Superfunctional.")
        .def("V_potential", &scf::HF::V_potential, "Returns the internal DFT V object.")
        .def("finalize", &scf::HF::finalize, "Cleans up the the Wavefunction's temporary data.")
        .def("soscf_update", &scf::HF::soscf_update, "Computes a second-order SCF update.")
        .def("occupation_a", &scf::HF::occupation_a, "Returns the Alpha occupation numbers.")
        .def("occupation_b", &scf::HF::occupation_b, "Returns the Beta occupation numbers.")
        .def("reset_occupation", &scf::HF::reset_occupation, "docstring")
        .def("compute_E", &scf::HF::compute_E, "docstring")
        .def("compute_initial_E", &scf::HF::compute_initial_E, "docstring")
        .def("save_density_and_energy", &scf::HF::save_density_and_energy, "docstring")
        .def("compute_orbital_gradient", &scf::HF::compute_orbital_gradient, "docstring")
        .def("find_occupation", &scf::HF::find_occupation, "docstring")
        .def("diis", &scf::HF::diis, "docstring")
        .def("diis_manager", &scf::HF::diis_manager, "docstring")
        .def_property("initialized_diis_manager_", &scf::HF::initialized_diis_manager,
                      &scf::HF::set_initialized_diis_manager, "docstring")
        .def("damping_update", &scf::HF::damping_update, "docstring")
        .def("check_phases", &scf::HF::check_phases, "docstring")
        .def("print_orbitals", &scf::HF::print_orbitals, "docstring")
        .def("print_header", &scf::HF::print_header, "docstring")
        .def("get_energies", &scf::HF::get_energies, "docstring")
        .def("set_energies", &scf::HF::set_energies, "docstring")
        .def("clear_external_potentials", &scf::HF::clear_external_potentials, "Clear external_potentials_ list")
        .def("push_back_external_potential", &scf::HF::push_back_external_potential,
             "Add an external potential to the external_potentials_ list", py::arg("V"))
        .def("print_preiterations", &scf::HF::print_preiterations, "docstring")
        .def_property("iteration_", &scf::HF::iteration, &scf::HF::set_iteration, "docstring")
        .def_property("diis_enabled_", &scf::HF::diis_enabled, &scf::HF::set_diis_enabled, "docstring")
        .def_property("diis_start_", &scf::HF::diis_start, &scf::HF::set_diis_start, "docstring")
        .def_property("frac_performed_", &scf::HF::frac_performed, &scf::HF::set_frac_performed,
                      "Frac performed current iteration?")
        .def_property("MOM_excited_", &scf::HF::MOM_excited, &scf::HF::set_MOM_excited,
                      "Are we to do excited-state MOM?")
        .def_property("MOM_performed_", &scf::HF::MOM_performed, &scf::HF::set_MOM_performed,
                      "MOM performed current iteration?")
        .def_property("attempt_number_", &scf::HF::attempt_number, &scf::HF::set_attempt_number,
                      "Current macroiteration (1-indexed) for stability analysis")
        .def("stability_analysis", &scf::HF::stability_analysis, "Assess wfn stability and correct if requested")
        .def("frac_renormalize", &scf::HF::frac_renormalize, "docstring")
        .def("compute_spin_contamination", &scf::HF::compute_spin_contamination, "docstring")
        .def("semicanonicalize", &scf::HF::semicanonicalize, "Semicanonicalizes the orbitals for ROHF.");

    /// HF Functions
    py::class_<scf::RHF, std::shared_ptr<scf::RHF>, scf::HF>(m, "RHF", "docstring")
        .def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>())
        .def("c1_deep_copy", &scf::RHF::c1_deep_copy,
             "Returns a new wavefunction with internal data converted to C_1 symmetry, using pre-c1-constructed "
             "BasisSet *basis*",
             py::arg("basis"));

    py::class_<scf::ROHF, std::shared_ptr<scf::ROHF>, scf::HF>(m, "ROHF", "docstring")
        .def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>())
        .def("moFeff", &scf::ROHF::moFeff, "docstring")
        .def("moFa", &scf::ROHF::moFa, "docstring")
        .def("moFb", &scf::ROHF::moFb, "docstring")
        .def("c1_deep_copy", &scf::ROHF::c1_deep_copy,
             "Returns a new wavefunction with internal data converted to C_1 symmetry, using pre-c1-constructed "
             "BasisSet *basis*",
             py::arg("basis"));

    py::class_<scf::UHF, std::shared_ptr<scf::UHF>, scf::HF>(m, "UHF", "docstring")
        .def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>())
        .def("c1_deep_copy", &scf::UHF::c1_deep_copy,
             "Returns a new wavefunction with internal data converted to C_1 symmetry, using pre-c1-constructed "
             "BasisSet *basis*",
             py::arg("basis"));

    py::class_<scf::CUHF, std::shared_ptr<scf::CUHF>, scf::HF>(m, "CUHF", "docstring")
        .def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>())
        .def("c1_deep_copy", &scf::CUHF::c1_deep_copy,
             "Returns a new wavefunction with internal data converted to C_1 symmetry, using pre-c1-constructed "
             "BasisSet *basis*",
             py::arg("basis"));

    /// EP2 functions
    py::class_<dfep2::DFEP2Wavefunction, std::shared_ptr<dfep2::DFEP2Wavefunction>, Wavefunction>(
        m, "DFEP2Wavefunction", "A density-fitted second-order Electron Propagator Wavefunction.")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def("compute", &dfep2::DFEP2Wavefunction::compute,
             "Computes the density-fitted EP2 energy for the input orbitals");

    /// FISAPT functions
    py::class_<fisapt::FISAPT, std::shared_ptr<fisapt::FISAPT>>(m, "FISAPT", "A Fragment-SAPT Wavefunction")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def("molecule", &fisapt::FISAPT::molecule, "Returns the FISAPT's molecule.")
        .def("scalars", &fisapt::FISAPT::scalars, "Return the interally computed scalars (not copied).")
        .def("vectors", &fisapt::FISAPT::vectors, "Return the interally computed vectors (not copied).")
        .def("matrices", &fisapt::FISAPT::matrices, "Return the interally computed matrices (not copied).")
        .def("print_header", &fisapt::FISAPT::print_header, "Print header, bases, sizes, etc.")
        .def("localize", &fisapt::FISAPT::localize, "Localize the active occupied orbitals via IBO2.")
        .def("partition", &fisapt::FISAPT::partition, "Partition the nuclei and electrons.")
        .def("overlap", &fisapt::FISAPT::overlap, "Build the overlap integrals S.")
        .def("kinetic", &fisapt::FISAPT::kinetic, "Build the kinetic integrals T.")
        .def("nuclear", &fisapt::FISAPT::nuclear, "Build the nuclear potentials V and interaction energies.")
        .def("coulomb", &fisapt::FISAPT::coulomb, "Build the J/K potentials for C, D, and E.")
        .def("scf", &fisapt::FISAPT::scf, "Solve the relaxed SCF equations for A0 and B0.")
        .def("freeze_core", &fisapt::FISAPT::freeze_core, "Freeze the core orbitals.")
        .def("unify", &fisapt::FISAPT::unify, "Produce unified matrices for A', B', and C'.")
        .def("dHF", &fisapt::FISAPT::dHF, "SAPT0 delta HF.")
        .def("elst", &fisapt::FISAPT::elst, "SAPT0 electrostatics.")
        .def("exch", &fisapt::FISAPT::exch, "SAPT0 exchange.")
        .def("ind", &fisapt::FISAPT::ind, "SAPT0 induction.")
        .def("disp", &fisapt::FISAPT::disp,
             "Computes the MP2-based DispE20 and Exch-DispE20 energy.")
        .def("flocalize", &fisapt::FISAPT::flocalize, "F-SAPT0 localize.")
        .def("felst", &fisapt::FISAPT::felst, "F-SAPT0 electrostatics.")
        .def("fexch", &fisapt::FISAPT::fexch, "F-SAPT0 exchange.")
        .def("find", &fisapt::FISAPT::find, "F-SAPT0 induction.")
        .def("fdisp", &fisapt::FISAPT::fdisp, "F-SAPT0 dispersion.")
        .def("raw_plot", &fisapt::FISAPT::raw_plot, "Plot some analysis files.")
        .def("print_trailer", &fisapt::FISAPT::print_trailer, "Print SAPT results.");

    /// CIWavefunction functions
    void (detci::CIWavefunction::*py_ci_sigma)(std::shared_ptr<psi::detci::CIvect>, std::shared_ptr<psi::detci::CIvect>,
                                               int, int) = &detci::CIWavefunction::sigma;
    void (detci::CIWavefunction::*py_ci_int_sigma)(std::shared_ptr<psi::detci::CIvect>,
                                                   std::shared_ptr<psi::detci::CIvect>, int, int, SharedVector,
                                                   SharedVector) = &detci::CIWavefunction::sigma;

    typedef std::vector<SharedMatrix> (detci::CIWavefunction::*form_density_sig)(
        std::shared_ptr<psi::detci::CIvect>, std::shared_ptr<psi::detci::CIvect>, int, int);

    py::class_<detci::CIWavefunction, std::shared_ptr<detci::CIWavefunction>, Wavefunction>(m, "CIWavefunction",
                                                                                            "docstring")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def("get_dimension", &detci::CIWavefunction::get_dimension, "docstring")
        .def("diag_h", &detci::CIWavefunction::diag_h, "docstring")
        .def("ndet", &detci::CIWavefunction::ndet, "docstring")
        .def("transform_ci_integrals", &detci::CIWavefunction::transform_ci_integrals, "docstring")
        .def("transform_mcscf_integrals", &detci::CIWavefunction::transform_mcscf_integrals, "docstring")
        .def("rotate_mcscf_integrals", &detci::CIWavefunction::rotate_mcscf_integrals, "docstring")
        .def("pitzer_to_ci_order_onel", &detci::CIWavefunction::pitzer_to_ci_order_onel, "docstring")
        .def("pitzer_to_ci_order_twoel", &detci::CIWavefunction::pitzer_to_ci_order_twoel, "docstring")
        .def("get_orbitals", &detci::CIWavefunction::get_orbitals, "docstring")
        .def("set_orbitals", &detci::CIWavefunction::set_orbitals, "docstring")
        .def("form_opdm", &detci::CIWavefunction::form_opdm, "docstring")
        .def("form_tpdm", &detci::CIWavefunction::form_tpdm, "docstring")
        .def("get_opdm", &detci::CIWavefunction::get_opdm, "docstring")
        .def("get_tpdm", &detci::CIWavefunction::get_tpdm, "docstring")
        .def("opdm", form_density_sig(&detci::CIWavefunction::opdm), "docstring")
        .def("tpdm", form_density_sig(&detci::CIWavefunction::tpdm), "docstring")
        .def("ci_nat_orbs", &detci::CIWavefunction::ci_nat_orbs, "docstring")
        .def("semicanonical_orbs", &detci::CIWavefunction::semicanonical_orbs, "docstring")
        .def("hamiltonian", &detci::CIWavefunction::hamiltonian, "docstring")
        .def("new_civector", &detci::CIWavefunction::new_civector, "docstring")
        .def("print_vector", &detci::CIWavefunction::print_vector, "docstring")
        .def("Hd_vector", &detci::CIWavefunction::Hd_vector, "docstring")
        .def("D_vector", &detci::CIWavefunction::D_vector, "docstring")
        .def("mcscf_object", &detci::CIWavefunction::mcscf_object, "docstring")
        .def("compute_state_transfer", &detci::CIWavefunction::compute_state_transfer, "docstring")
        .def("sigma", py_ci_sigma, "docstring")
        .def("sigma", py_ci_int_sigma, "docstring")
        .def("cleanup_ci", &detci::CIWavefunction::cleanup_ci, "docstring")
        .def("cleanup_dpd", &detci::CIWavefunction::cleanup_dpd, "docstring")
        .def("set_ci_guess", &detci::CIWavefunction::set_ci_guess, "docstring");

    void (detci::CIvect::*py_civ_copy)(std::shared_ptr<psi::detci::CIvect>, int, int) = &detci::CIvect::copy;
    void (detci::CIvect::*py_civ_scale)(double, int) = &detci::CIvect::scale;

    py::class_<detci::CIvect, std::shared_ptr<detci::CIvect>>(m, "CIVector", py::buffer_protocol(), "docstring")
        .def("vdot", &detci::CIvect::vdot, "docstring")
        .def("axpy", &detci::CIvect::axpy, "docstring")
        .def("vector_multiply", &detci::CIvect::vector_multiply, "docstring")
        .def("copy", py_civ_copy, "docstring")
        .def("zero", &detci::CIvect::zero, "docstring")
        .def("divide", &detci::CIvect::divide, "docstring")
        .def("scale", py_civ_scale, "docstring")
        .def("norm", &detci::CIvect::norm, "docstring")
        .def("shift", &detci::CIvect::shift, "docstring")
        .def("dcalc", &detci::CIvect::dcalc3, "docstring")
        .def("symnormalize", &detci::CIvect::symnormalize, "docstring")
        .def("read", &detci::CIvect::read, "docstring")
        .def("write", &detci::CIvect::write, "docstring")
        .def("init_io_files", &detci::CIvect::init_io_files, "docstring")
        .def("close_io_files", &detci::CIvect::close_io_files, "docstring")
        .def("set_nvec", &detci::CIvect::set_nvect, "docstring")
        .def_buffer([](detci::CIvect& vec) { return vec.array_interface(); });
}
