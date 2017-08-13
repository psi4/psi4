 /*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/orbitalspace.h"
#include "psi4/libmints/extern.h"

#include "psi4/libfock/jk.h"
#include "psi4/libfock/soscf.h"

#include "psi4/detci/ciwave.h"
#include "psi4/detci/civect.h"

#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"

#include "psi4/libscf_solver/hf.h"
#include "psi4/libscf_solver/rhf.h"
#include "psi4/libscf_solver/uhf.h"
#include "psi4/libscf_solver/rohf.h"
#include "psi4/libscf_solver/cuhf.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libfock/v.h"

#include "psi4/dfep2/dfep2.h"

#include "psi4/fisapt/fisapt.h"

#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

#include <string>

using namespace psi;
void export_wavefunction(py::module& m) {
    typedef void (Wavefunction::*take_sharedwfn)(SharedWavefunction);
    py::class_<Wavefunction, std::shared_ptr<Wavefunction>>(m, "Wavefunction", "docstring",
                                                            py::dynamic_attr())
        .def(py::init<std::shared_ptr<Molecule>, std::shared_ptr<BasisSet>, Options&>())
        .def(py::init<std::shared_ptr<Molecule>, std::shared_ptr<BasisSet>>())
        .def("reference_wavefunction", &Wavefunction::reference_wavefunction, "Returns the reference wavefunction.")
        .def("set_reference_wavefunction", &Wavefunction::set_reference_wavefunction, "docstring")
        .def("shallow_copy", take_sharedwfn(&Wavefunction::shallow_copy),
             "Copies the pointers to the internal data.")
        .def("deep_copy", take_sharedwfn(&Wavefunction::deep_copy),
             "Deep copies the internal data.")
        .def("same_a_b_orbs", &Wavefunction::same_a_b_orbs,
             "Returns true if the alpha and beta orbitals are the same.")
        .def("same_a_b_dens", &Wavefunction::same_a_b_dens,
             "Returns true if the alpha and beta densities are the same.")
        .def("nfrzc", &Wavefunction::nfrzc, "Number of frozen core electrons.")
        .def("nalpha", &Wavefunction::nalpha, "Number of Alpha electrons.")
        .def("nbeta", &Wavefunction::nbeta, "Number of Beta electrons.")
        .def("nso", &Wavefunction::nso, "Number of symmetry orbitals.")
        .def("nmo", &Wavefunction::nmo, "Number of molecule orbitals.")
        .def("nirrep", &Wavefunction::nirrep, "Number of irreps in the system.")
        .def("Ca", &Wavefunction::Ca, "Returns the Alpha Orbitals.")
        .def("Cb", &Wavefunction::Cb, "Returns the Beta Orbitals.")
        .def("Ca_subset", &Wavefunction::Ca_subset, py::return_value_policy::take_ownership,
             "Returns the requested Alpha Orbital subset.")
        .def("Cb_subset", &Wavefunction::Cb_subset, py::return_value_policy::take_ownership,
             "Returns the requested Beta Orbital subset.")
        .def("Fa", &Wavefunction::Fa, "Returns the Alpha Fock Matrix.")
        .def("Fb", &Wavefunction::Fb, "Returns the Beta Fock Matrix.")
        .def("Da", &Wavefunction::Da, "Returns the Alpha Density Matrix.")
        .def("Db", &Wavefunction::Db, "Returns the Beta Density Matrix.")
        .def("Da_subset", &Wavefunction::Da_subset, py::return_value_policy::take_ownership,
             "Returns the requested Alpha Density subset.")
        .def("Db_subset", &Wavefunction::Db_subset, py::return_value_policy::take_ownership,
             "Returns the requested Beta Density subset.")
        .def("epsilon_a", &Wavefunction::epsilon_a, "Returns the Alpha Eigenvalues.")
        .def("epsilon_b", &Wavefunction::epsilon_b, "Returns the Beta Eigenvalues.")
        .def("epsilon_a_subset", &Wavefunction::epsilon_a_subset,
             "Returns the requested Alpha Eigenvalues subset.")
        .def("epsilon_b_subset", &Wavefunction::epsilon_b_subset,
             "Returns the requested Beta Eigenvalues subset.")
        .def("X", &Wavefunction::X, "Returns the Lagrangian Matrix.")
        .def("basis_projection", &Wavefunction::basis_projection,
             "Projects a orbital matrix from one basis to another.")
        .def("H", &Wavefunction::H, "Returns the 'Core' Matrix (Potential + Kinetic) Integrals.")
        .def("S", &Wavefunction::S, "Returns the One-electron Overlap Matrix.")
        .def("aotoso", &Wavefunction::aotoso,
             "Returns the Atomic Orbital to Symmetry Orbital transformer.")
        .def("basisset", &Wavefunction::basisset, "Returns the current orbital basis.")
        .def("sobasisset", &Wavefunction::sobasisset, "Returns the symmetry orbitals basis.")
        .def("get_basisset", &Wavefunction::get_basisset, "Returns the requested auxiliary basis.")
        .def("set_basisset", &Wavefunction::set_basisset, "Sets the requested auxiliary basis.")
        .def("energy", &Wavefunction::reference_energy, "Returns the Wavefunctions energy.")
        .def("gradient", &Wavefunction::gradient, "Returns the Wavefunctions gradient.")
        .def("set_gradient", &Wavefunction::set_gradient, "Sets the Wavefunctions gradient.")
        .def("hessian", &Wavefunction::hessian, "Returns the Wavefunctions Hessian.")
        .def("set_hessian", &Wavefunction::set_hessian, "Sets the Wavefunctions Hessian.")
        .def("frequencies", &Wavefunction::frequencies, "Returns the frequencies of the Hessian.")
        .def("set_frequencies", &Wavefunction::set_frequencies,
             "Sets the frequencies of the Hessian.")
        .def("esp_at_nuclei", &Wavefunction::get_esp_at_nuclei,
             "returns electrostatic potentials at nuclei")
        .def("mo_extents", &Wavefunction::get_mo_extents,
             "returns the wavefunction's electronic orbital extents.")
        .def("atomic_point_charges", &Wavefunction::get_atomic_point_charges,
             "Returns the set atomic point charges.")
        .def("no_occupations", &Wavefunction::get_no_occupations,
             "returns the natural orbital occupations on the wavefunction.")
        .def("normalmodes", &Wavefunction::get_normalmodes,
             "Returns the normal modes of the Wavefunction.")
        .def("normalmode_displacements", &Wavefunction::get_normalmodes_displacements,
             "Returns the displacements for normalmodes")
        .def("set_name", &Wavefunction::set_name,
             "Sets the level of theory this wavefunction corresponds to.")
        .def("name", &Wavefunction::name, py::return_value_policy::copy,
             "The level of theory this wavefunction corresponds to.")
        .def("alpha_orbital_space", &Wavefunction::alpha_orbital_space, "docstring")
        .def("beta_orbital_space", &Wavefunction::beta_orbital_space, "docstring")
        .def("molecule", &Wavefunction::molecule, "Returns the Wavefunctions molecule.")
        .def("doccpi", &Wavefunction::doccpi, py::return_value_policy::copy,
             "Returns the number of doubly occupied orbitals per irrep.")
        .def("soccpi", &Wavefunction::soccpi, py::return_value_policy::copy,
             "Returns the number of singly occupied orbitals per irrep.")
        .def("nsopi", &Wavefunction::nsopi, py::return_value_policy::copy,
             "Returns the number of symmetry orbitals per irrep.")
        .def("nmopi", &Wavefunction::nmopi, py::return_value_policy::copy,
             "Returns the number of molecular orbitals per irrep.")
        .def("nalphapi", &Wavefunction::nalphapi, py::return_value_policy::copy,
             "Returns the number of alpha orbitals per irrep.")
        .def("nbetapi", &Wavefunction::nbetapi, py::return_value_policy::copy,
             "Returns the number of beta orbitals per irrep.")
        .def("frzcpi", &Wavefunction::frzcpi, py::return_value_policy::copy,
             "Returns the number of frozen core orbitals per irrep.")
        .def("frzvpi", &Wavefunction::frzvpi, py::return_value_policy::copy,
             "Returns the number of frozen virtual orbitals per irrep.")
        .def("set_print", &Wavefunction::set_print, "Sets the print level of the Wavefunction.")
        .def("get_print", &Wavefunction::get_print, "Get the print level of the Wavefunction.")
        .def("compute_energy", &Wavefunction::compute_energy,
             "Computes the energy of the Wavefunction.")
        .def("compute_gradient", &Wavefunction::compute_gradient,
             "Computes the gradient of the Wavefunction")
        .def("compute_hessian", &Wavefunction::compute_hessian,
             "Computes the Hessian of the Wavefunction.")
        .def("set_external_potential", &Wavefunction::set_external_potential, "Sets the requested external potential.")
        .def("set_variable", &Wavefunction::set_variable, "Sets the requested internal variable.")
        .def("get_variable", &Wavefunction::get_variable,
             "Returns the requested internal variable.")
        .def("variables", &Wavefunction::variables, "Returns the map of all internal variables.")
        .def("get_array", &Wavefunction::get_array, "Sets the requested internal array.")
        .def("set_array", &Wavefunction::set_array, "Returns the requested internal array.")
        .def("arrays", &Wavefunction::arrays, "Returns the map of all internal arrays.");

    py::class_<scf::HF, std::shared_ptr<scf::HF>, Wavefunction>(m, "HF", "docstring")
        .def("form_C", &scf::HF::form_C,
             "Forms the Orbital Matrices from the current Fock Matrices.")
        .def("form_D", &scf::HF::form_D,
             "Forms the Density Matrices from the current Orbitals Matrices")
        .def("form_V", &scf::HF::form_V,
             "Form the Kohn-Sham Potential Matrices from the current Density Matrices")
        .def("form_G", &scf::HF::form_G, "Forms the G matrix.")
        .def("form_F", &scf::HF::form_V, "Forms the F matrix.")
        .def("onel_Hx", &scf::HF::onel_Hx, "One-electron Hessian-vector products.")
        .def("twoel_Hx", &scf::HF::twoel_Hx, "Two-electron Hessian-vector products")
        .def("cphf_Hx", &scf::HF::cphf_Hx, "CPHF Hessian-vector prodcuts (4 * J - K - K.T).")
        .def("cphf_solve", &scf::HF::cphf_solve, py::arg("x_vec"), py::arg("conv_tol"),
             py::arg("max_iter"), py::arg("print_lvl") = 2,
             "Solves the CPHF equations for a given set of x vectors.")
        .def("cphf_converged", &scf::HF::cphf_converged, "Adds occupied guess alpha orbitals.")
        .def("guess_Ca", &scf::HF::guess_Ca, "Sets the guess Alpha Orbital Matrix")
        .def("guess_Cb", &scf::HF::guess_Cb, "Sets the guess Beta Orbital Matrix")
        .def("reset_occ", &scf::HF::reset_occ,
             "If True, the occupation will be reset after the guess to the inital occupation.")
        .def("set_sad_basissets", &scf::HF::set_sad_basissets,
             "Sets the Superposition of Atomic Densities basisset.")
        .def("set_sad_fitting_basissets", &scf::HF::set_sad_fitting_basissets,
             "Sets the Superposition of Atomic Densities density-fitted basisset.")
        .def("Va", &scf::HF::Va, "Returns the Alpha Kohn-Shame Potential Matrix.")
        .def("Vb", &scf::HF::Vb, "Returns the Alpha Kohn-Shame Potential Matrix.")
        .def("jk", &scf::HF::jk, "Returns the internal JK object.")
        .def("functional", &scf::HF::functional, "Returns the internal DFT Superfunctional.")
        .def("V_potential", &scf::HF::V_potential, "Returns the internal DFT V object.")
        .def("initialize", &scf::HF::initialize, "Initializes the Wavefunction.")
        .def("iterations", &scf::HF::iterations,
             "Iterates the Wavefunction until convergence criteria have been met.")
        .def("finalize_E", &scf::HF::finalize_E, "Computes the final SCF energy.")
        .def("occupation_a", &scf::HF::occupation_a, "Returns the Alpha occupation numbers.")
        .def("occupation_b", &scf::HF::occupation_b, "Returns the Beta occupation numbers.")
        .def("semicanonicalize", &scf::HF::semicanonicalize,
             "Semicanonicalizes the orbitals for ROHF.");
    
    py::class_<findif::VIBRATION, std::shared_ptr<findif::VIBRATION>>(m, "VIBRATION" , "docstring")
        .def(py::init<int, int, double>());
    
    py::class_<scf::RHF, std::shared_ptr<scf::RHF>, scf::HF>(m, "RHF", "docstring")
        .def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>());

    py::class_<scf::ROHF, std::shared_ptr<scf::ROHF>, scf::HF>(m, "ROHF", "docstring")
        .def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>())
        .def("moFeff", &scf::ROHF::moFeff, "docstring")
        .def("moFa", &scf::ROHF::moFa, "docstring")
        .def("moFb", &scf::ROHF::moFb, "docstring");

    py::class_<scf::UHF, std::shared_ptr<scf::UHF>, scf::HF>(m, "UHF", "docstring")
        .def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>());

    py::class_<scf::CUHF, std::shared_ptr<scf::CUHF>, scf::HF>(m, "CUHF", "docstring")
        .def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>());

    py::class_<dfep2::DFEP2Wavefunction, std::shared_ptr<dfep2::DFEP2Wavefunction>, Wavefunction>(
        m, "DFEP2Wavefunction", "A density-fitted second-order Electron Propagator Wavefunction.")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def("compute", &dfep2::DFEP2Wavefunction::compute,
             "Computes the density-fitted EP2 energy for the input orbitals");

    py::class_<fisapt::FISAPT, std::shared_ptr<fisapt::FISAPT>>(m, "FISAPT",
                                                                "A Fragment-SAPT Wavefunction")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def("compute_energy", &fisapt::FISAPT::compute_energy, "Computes the FSAPT energy.")
        .def("scalars", &fisapt::FISAPT::scalars, "Return the interally computed scalars.")
        .def("disp", &fisapt::FISAPT::disp,
             "Computes the MP2-based DispE20 and Exch-DispE20 energy.");

    /// CIWavefunction data
    void (detci::CIWavefunction::*py_ci_sigma)(std::shared_ptr<psi::detci::CIvect>,
                                               std::shared_ptr<psi::detci::CIvect>, int, int) =
        &detci::CIWavefunction::sigma;
    void (detci::CIWavefunction::*py_ci_int_sigma)(
        std::shared_ptr<psi::detci::CIvect>, std::shared_ptr<psi::detci::CIvect>, int, int,
        SharedVector, SharedVector) = &detci::CIWavefunction::sigma;

    typedef std::vector<SharedMatrix>(detci::CIWavefunction::*form_density_sig)(
        std::shared_ptr<psi::detci::CIvect>, std::shared_ptr<psi::detci::CIvect>, int, int);

    py::class_<detci::CIWavefunction, std::shared_ptr<detci::CIWavefunction>, Wavefunction>(
        m, "CIWavefunction", "docstring")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def("get_dimension", &detci::CIWavefunction::get_dimension, "docstring")
        .def("diag_h", &detci::CIWavefunction::diag_h, "docstring")
        .def("ndet", &detci::CIWavefunction::ndet, "docstring")
        .def("transform_ci_integrals", &detci::CIWavefunction::transform_ci_integrals, "docstring")
        .def("transform_mcscf_integrals", &detci::CIWavefunction::transform_mcscf_integrals,
             "docstring")
        .def("rotate_mcscf_integrals", &detci::CIWavefunction::rotate_mcscf_integrals, "docstring")
        .def("pitzer_to_ci_order_onel", &detci::CIWavefunction::pitzer_to_ci_order_onel,
             "docstring")
        .def("pitzer_to_ci_order_twoel", &detci::CIWavefunction::pitzer_to_ci_order_twoel,
             "docstring")
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

    void (detci::CIvect::*py_civ_copy)(std::shared_ptr<psi::detci::CIvect>, int, int) =
        &detci::CIvect::copy;
    void (detci::CIvect::*py_civ_scale)(double, int) = &detci::CIvect::scale;

    py::class_<detci::CIvect, std::shared_ptr<detci::CIvect>>(m, "CIVector", py::buffer_protocol(),
                                                              "docstring")
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
