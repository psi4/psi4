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
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/orbitalspace.h"

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

#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

#include <string>

using namespace psi;
void export_wavefunction(py::module& m)
{

    typedef void (Wavefunction::*take_sharedwfn)(SharedWavefunction);
    py::class_<Wavefunction, std::shared_ptr<Wavefunction>>(m, "Wavefunction", "docstring", py::dynamic_attr()).
            def(py::init<std::shared_ptr<Molecule>, std::shared_ptr<BasisSet>, Options&>()).
            def(py::init<std::shared_ptr<Molecule>, std::shared_ptr<BasisSet>>()).
            def("reference_wavefunction", &Wavefunction::reference_wavefunction, "docstring").
            def("set_reference_wavefunction", &Wavefunction::set_reference_wavefunction, "docstring").
            def("shallow_copy", take_sharedwfn(&Wavefunction::shallow_copy), "docstring").
            def("deep_copy", take_sharedwfn(&Wavefunction::deep_copy), "docstring").
            def("same_a_b_orbs", &Wavefunction::same_a_b_orbs, "docstring").
            def("same_a_b_dens", &Wavefunction::same_a_b_dens, "docstring").
            def("nfrzc", &Wavefunction::nfrzc, "docstring").
            def("nalpha", &Wavefunction::nalpha, "docstring").
            def("nbeta", &Wavefunction::nbeta, "docstring").
            def("nso", &Wavefunction::nso, "docstring").
            def("nmo", &Wavefunction::nmo, "docstring").
            def("nirrep", &Wavefunction::nirrep, "docstring").
            def("Ca_subset", &Wavefunction::Ca_subset, py::return_value_policy::take_ownership, "docstring").
            def("Cb_subset", &Wavefunction::Cb_subset, py::return_value_policy::take_ownership, "docstring").
            def("epsilon_a_subset", &Wavefunction::epsilon_a_subset, "docstring").
            def("epsilon_b_subset", &Wavefunction::epsilon_b_subset, "docstring").
            def("Ca", &Wavefunction::Ca, "docstring").
            def("Cb", &Wavefunction::Cb, "docstring").
            def("Fa", &Wavefunction::Fa, "docstring").
            def("Fb", &Wavefunction::Fb, "docstring").
            def("Da", &Wavefunction::Da, "docstring").
            def("Db", &Wavefunction::Db, "docstring").
            def("X", &Wavefunction::X, "docstring").
            def("basis_projection", &Wavefunction::basis_projection, "docstring").
            def("H", &Wavefunction::H, "docstring").
            def("S", &Wavefunction::S, "docstring").
            def("aotoso", &Wavefunction::aotoso, "docstring").
            def("epsilon_a", &Wavefunction::epsilon_a, "docstring").
            def("epsilon_b", &Wavefunction::epsilon_b, "docstring").
            def("basisset", &Wavefunction::basisset, "docstring").
            def("get_basisset", &Wavefunction::get_basisset, "docstring").
            def("set_basisset", &Wavefunction::set_basisset, "docstring").
            def("sobasisset", &Wavefunction::sobasisset, "docstring").
            def("energy", &Wavefunction::reference_energy, "docstring").
            def("gradient", &Wavefunction::gradient, "docstring").
            def("set_gradient", &Wavefunction::set_gradient, "docstring").
            def("hessian", &Wavefunction::hessian, "docstring").
            def("set_hessian", &Wavefunction::set_hessian, "docstring").
            def("frequencies", &Wavefunction::frequencies, "docstring").
            def("set_frequencies", &Wavefunction::set_frequencies, "docstring").
            def("atomic_point_charges", &Wavefunction::get_atomic_point_charges, "docstring").
            def("normalmodes", &Wavefunction::normalmodes, "docstring").
            def("name", &Wavefunction::name, py::return_value_policy::copy, "The level of theory this wavefunction corresponds to.").
            def("alpha_orbital_space", &Wavefunction::alpha_orbital_space, "docstring").
            def("beta_orbital_space", &Wavefunction::beta_orbital_space, "docstring").
            def("molecule", &Wavefunction::molecule, "docstring").
            def("doccpi", &Wavefunction::doccpi, py::return_value_policy::copy, "docstring").
            def("soccpi", &Wavefunction::soccpi, py::return_value_policy::copy, "docstring").
            def("nsopi", &Wavefunction::nsopi, py::return_value_policy::copy, "docstring").
            def("nmopi", &Wavefunction::nmopi, py::return_value_policy::copy, "docstring").
            def("nalphapi", &Wavefunction::nalphapi, py::return_value_policy::copy, "docstring").
            def("nbetapi", &Wavefunction::nbetapi, py::return_value_policy::copy, "docstring").
            def("frzcpi", &Wavefunction::frzcpi, py::return_value_policy::copy, "docstring").
            def("frzvpi", &Wavefunction::frzvpi, py::return_value_policy::copy, "docstring").
            def("nalpha", &Wavefunction::nalpha, "docstring").
            def("nbeta", &Wavefunction::nbeta, "docstring").
            def("set_oeprop", &Wavefunction::set_oeprop, "Associate an OEProp object with this wavefunction").
            def("oeprop", &Wavefunction::get_oeprop, "Get the OEProp object associated with this wavefunction").
            def("set_print", &Wavefunction::set_print, "docstring").
            def("compute_energy", &Wavefunction::compute_energy, "docstring").
            def("compute_gradient", &Wavefunction::compute_gradient, "docstring").
            def("get_variable", &Wavefunction::get_variable, "docstring").
            def("set_variable", &Wavefunction::set_variable, "docstring").
            def("variables", &Wavefunction::variables, "docstring").
            def("get_array", &Wavefunction::get_array, "docstring").
            def("set_array", &Wavefunction::set_array, "docstring").
            def("arrays", &Wavefunction::arrays, "docstring");

    py::class_<scf::HF, std::shared_ptr<scf::HF>, Wavefunction>(m, "HF", "docstring").
            def("form_C", &scf::HF::form_C, "docstring").
            def("form_D", &scf::HF::form_D, "docstring").
            def("form_V", &scf::HF::form_V, "docstring").
            def("guess_Ca", &scf::HF::guess_Ca, "docstring").
            def("guess_Cb", &scf::HF::guess_Cb, "docstring").
            def("reset_occ", &scf::HF::reset_occ, "docstring").
            def("set_sad_basissets", &scf::HF::set_sad_basissets, "docstring").
            def("set_sad_fitting_basissets", &scf::HF::set_sad_fitting_basissets, "docstring").
            def("Va", &scf::HF::Va, "docstring").
            def("Vb", &scf::HF::Vb, "docstring").
            def("jk", &scf::HF::jk, "docstring").
            def("functional", &scf::HF::functional, "docstring").
            def("V_potential", &scf::HF::V_potential, "docstring").
            def("initialize", &scf::HF::initialize, "docstring").
            def("iterations", &scf::HF::iterations, "docstring").
            def("finalize_E", &scf::HF::finalize_E, "docstring").
            def("occupation_a", &scf::HF::occupation_a, "docstring").
            def("occupation_b", &scf::HF::occupation_b, "docstring").
            def("semicanonicalize", &scf::HF::semicanonicalize, "docstring");

    py::class_<scf::RHF, std::shared_ptr<scf::RHF>, scf::HF>(m, "RHF", "docstring").
            def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>());

    py::class_<scf::ROHF, std::shared_ptr<scf::ROHF>, scf::HF>(m, "ROHF", "docstring").
            def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>()).
            def("moFeff", &scf::ROHF::moFeff, "docstring").
            def("moFa", &scf::ROHF::moFa, "docstring").
            def("moFb", &scf::ROHF::moFb, "docstring");

    py::class_<scf::UHF, std::shared_ptr<scf::UHF>, scf::HF>(m, "UHF", "docstring").
            def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>());

    py::class_<scf::CUHF, std::shared_ptr<scf::CUHF>, scf::HF>(m, "CUHF", "docstring").
            def(py::init<std::shared_ptr<Wavefunction>, std::shared_ptr<SuperFunctional>>());

    /// CIWavefunction data
    void (detci::CIWavefunction::*py_ci_sigma)(std::shared_ptr<psi::detci::CIvect>,
                                    std::shared_ptr<psi::detci::CIvect>, int, int) =
                                    &detci::CIWavefunction::sigma;
    void (detci::CIWavefunction::*py_ci_int_sigma)(std::shared_ptr<psi::detci::CIvect>,
                                    std::shared_ptr<psi::detci::CIvect>, int, int,
                                    SharedVector, SharedVector) =
                                    &detci::CIWavefunction::sigma;

    typedef std::vector<SharedMatrix> (detci::CIWavefunction::*form_density_sig)(
                                          std::shared_ptr<psi::detci::CIvect>,
                                          std::shared_ptr<psi::detci::CIvect>,
                                          int, int);

    py::class_<detci::CIWavefunction, std::shared_ptr<detci::CIWavefunction>, Wavefunction>(m, "CIWavefunction", "docstring")
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

    py::class_<detci::CIvect, std::shared_ptr<detci::CIvect> >(m, "CIVector", py::buffer_protocol(), "docstring")
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
        .def_buffer([](detci::CIvect &vec){
            return vec.array_interface();
            });

}
