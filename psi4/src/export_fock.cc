/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "psi4/pybind11.h"

#include "psi4/libfock/jk.h"
#include "psi4/libfock/soscf.h"
#include "psi4/lib3index/denominator.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/lib3index/dfhelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libscf_solver/sad.h"

using namespace psi;
namespace py = pybind11;
using namespace pybind11::literals;

void export_fock(py::module &m) {
    py::class_<JK, std::shared_ptr<JK>>(m, "JK", "docstring")
        .def_static("build_JK",
                    [](std::shared_ptr<BasisSet> basis, std::shared_ptr<BasisSet> aux) {
                        return JK::build_JK(basis, aux, Process::environment.options);
                    })
        .def_static("build_JK",
                    [](std::shared_ptr<BasisSet> basis, std::shared_ptr<BasisSet> aux, bool do_wK, size_t doubles) {
                        return JK::build_JK(basis, aux, Process::environment.options, do_wK, doubles);
                    })
        .def("name", &JK::name)
        .def("memory_estimate", &JK::memory_estimate)
        .def("initialize", &JK::initialize)
        .def("basisset", &JK::basisset)
        .def("set_print", &JK::set_print)
        .def("set_cutoff", &JK::set_cutoff)
        .def("set_memory", &JK::set_memory)
        .def("set_omp_nthread", &JK::set_omp_nthread)
        .def("set_do_J", &JK::set_do_J)
        .def("set_do_K", &JK::set_do_K)
        .def("set_do_wK", &JK::set_do_wK)
        .def("set_omega", &JK::set_omega, "Dampening term for range separated DFT", "omega"_a)
        .def("get_omega", &JK::get_omega, "Dampening term for range separated DFT")
        .def("set_wcombine", &JK::set_wcombine, "Are Exchange terms in one Matrix", "wcombine"_a )
        .def("get_wcombine", &JK::get_wcombine, "Are Exchange terms in one Matrix", "wcombine")
        .def("set_omega_alpha", &JK::set_omega_alpha, "Weight for HF exchange term in range-separated DFT", "alpha"_a)
        .def("get_omega_alpha", &JK::get_omega_alpha, "Weight for HF exchange term in range-separated DFT")
        .def("set_omega_beta", &JK::set_omega_beta, "Weight for dampened exchange term in range-separated DFT", "beta"_a)
        .def("get_omega_beta", &JK::get_omega_beta, "Weight for dampened exchange term in range-separated DFT")
        .def("set_early_screening", &JK::set_early_screening, "Use severe screening techniques? Useful in early SCF iterations.", "early_screening"_a)
        .def("get_early_screening", &JK::get_early_screening, "Use severe screening techniques? Useful in early SCF iterations.")
        .def("compute", &JK::compute)
        .def("finalize", &JK::finalize)
        .def("C_clear",
             [](JK &jk) {
                 jk.C_left().clear();
                 jk.C_right().clear();
             })
        .def("C_add",
             [](JK &jk, SharedMatrix Cl) {
                 jk.C_left().push_back(Cl);
                 jk.C_right().push_back(Cl);
             })
        .def("C_left_add", [](JK &jk, SharedMatrix Cl) { jk.C_left().push_back(Cl); })
        .def("C_right_add", [](JK &jk, SharedMatrix Cr) { jk.C_right().push_back(Cr); })
        .def("J", &JK::J, py::return_value_policy::reference_internal)
        .def("K", &JK::K, py::return_value_policy::reference_internal)
        .def("wK", &JK::wK, py::return_value_policy::reference_internal)
        .def("D", &JK::D, py::return_value_policy::reference_internal)
        .def("computed_shells_per_iter", &JK::computed_shells_per_iter, "Array containing the number of ERI shell quartets computed (not screened out) during each compute call.")
        .def("print_header", &JK::print_header, "docstring");

    py::class_<LaplaceDenominator, std::shared_ptr<LaplaceDenominator>>(m, "LaplaceDenominator", "Computer class for a Laplace factorization of the four-index energy denominator in MP2 and coupled-cluster")
        .def(py::init<std::shared_ptr<Vector>, std::shared_ptr<Vector>, double>())
        .def("denominator_occ", &LaplaceDenominator::denominator_occ, "Returns the occupied orbital Laplace weights of the factorized doubles denominator (nweights * nocc)")
        .def("denominator_vir", &LaplaceDenominator::denominator_vir, "Returns the virtual orbital Laplace weights of the factorized doubles denominator (nweights * nvirt)");

    py::class_<TLaplaceDenominator, std::shared_ptr<TLaplaceDenominator>>(m, "TLaplaceDenominator", "Computer class for a Laplace factorization of the six-index energy denominator in coupled-cluster theory")
        .def(py::init<std::shared_ptr<Vector>, std::shared_ptr<Vector>, double>())
        .def("denominator_occ", &TLaplaceDenominator::denominator_occ, "Returns the occupied orbital Laplace weights of the factorized triples denominator (nweights * nocc)")
        .def("denominator_vir", &TLaplaceDenominator::denominator_vir, "Returns the virtual orbital Laplace weights of the factorized triples denominator (nweights * nvirt)");

    py::class_<DFTensor, std::shared_ptr<DFTensor>>(m, "DFTensor", "docstring")
        .def(py::init<std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, std::shared_ptr<Matrix>, int, int>())
        .def("Qso", &DFTensor::Qso, "doctsring")
        .def("Qmo", &DFTensor::Qmo, "doctsring")
        .def("Qoo", &DFTensor::Qoo, "doctsring")
        .def("Qov", &DFTensor::Qov, "doctsring")
        .def("Qvv", &DFTensor::Qvv, "doctsring")
        .def("Imo", &DFTensor::Imo, "doctsring")
        .def("Idfmo", &DFTensor::Idfmo, "doctsring");

    py::class_<FittingMetric, std::shared_ptr<FittingMetric>>(m, "FittingMetric", "docstring")
        .def(py::init<std::shared_ptr<BasisSet>, bool>())
        .def("get_algorithm", &FittingMetric::get_algorithm, "docstring")
        .def("is_poisson", &FittingMetric::is_poisson, "docstring")
        .def("is_inverted", &FittingMetric::is_inverted, "docstring")
        .def("get_metric", &FittingMetric::get_metric, "docstring")
        .def("get_pivots", &FittingMetric::get_pivots, "docstring")
        .def("get_reverse_pivots", &FittingMetric::get_reverse_pivots, "docstring")
        .def("form_fitting_metric", &FittingMetric::form_fitting_metric, "docstring")
        .def("form_cholesky_inverse", &FittingMetric::form_cholesky_inverse, "docstring")
        .def("form_QR_inverse", &FittingMetric::form_QR_inverse, "docstring")
        .def("form_eig_inverse", &FittingMetric::form_eig_inverse, "docstring")
        .def("form_full_inverse", &FittingMetric::form_full_inverse, "docstring");

    py::class_<SOMCSCF, std::shared_ptr<SOMCSCF>>(m, "SOMCSCF", "docstring")
        // .def(init<std::shared_ptr<JK>, SharedMatrix, SharedMatrix >())
        .def("Ck", &SOMCSCF::Ck)
        .def("form_rotation_matrix", &SOMCSCF::form_rotation_matrix, "x"_a, "order"_a = 2)
        .def("rhf_energy", &SOMCSCF::rhf_energy)
        .def("update", &SOMCSCF::update)
        .def("approx_solve", &SOMCSCF::approx_solve)
        .def("solve", &SOMCSCF::solve)
        .def("H_approx_diag", &SOMCSCF::H_approx_diag)
        .def("compute_Hk", &SOMCSCF::Hk)
        .def("compute_Q", &SOMCSCF::compute_Q)
        .def("compute_Qk", &SOMCSCF::compute_Qk)
        .def("compute_AFock", &SOMCSCF::compute_AFock)
        .def("current_total_energy", &SOMCSCF::current_total_energy)
        .def("current_docc_energy", &SOMCSCF::current_docc_energy)
        .def("current_ci_energy", &SOMCSCF::current_ci_energy)
        .def("current_AFock", &SOMCSCF::current_AFock)
        .def("current_IFock", &SOMCSCF::current_IFock)
        .def("zero_redundant", &SOMCSCF::zero_redundant)
        .def("gradient", &SOMCSCF::gradient)
        .def("gradient_rms", &SOMCSCF::gradient_rms);

    py::class_<DFSOMCSCF, std::shared_ptr<DFSOMCSCF>, SOMCSCF>(m, "DFSOMCSCF", "docstring");
    py::class_<DiskSOMCSCF, std::shared_ptr<DiskSOMCSCF>, SOMCSCF>(m, "DiskSOMCSCF", "docstring");

    // DF Helper
    typedef SharedMatrix (DFHelper::*take_string)(std::string);
    typedef SharedMatrix (DFHelper::*tensor_access3)(std::string, std::vector<size_t>, std::vector<size_t>,
                                                     std::vector<size_t>);

    py::class_<DFHelper, std::shared_ptr<DFHelper>>(m, "DFHelper", "docstring")
        .def(py::init<std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet> >())
        .def("set_memory", &DFHelper::set_memory)
        .def("get_memory", &DFHelper::get_memory)
        .def("set_method", &DFHelper::set_method)
        .def("get_method", &DFHelper::get_method)
        .def("set_subalgo", &DFHelper::set_subalgo)
        .def("get_AO_size", &DFHelper::get_AO_size)
        .def("set_nthreads", &DFHelper::set_nthreads)
        .def("hold_met", &DFHelper::hold_met)
        .def("set_schwarz_cutoff", &DFHelper::set_schwarz_cutoff)
        .def("get_schwarz_cutoff", &DFHelper::get_schwarz_cutoff)
        .def("set_AO_core", &DFHelper::set_AO_core)
        .def("get_AO_core", &DFHelper::get_AO_core)
        .def("set_MO_core", &DFHelper::set_MO_core)
        .def("get_MO_core", &DFHelper::get_MO_core)
        .def("add_space", &DFHelper::add_space)
        .def("initialize", &DFHelper::initialize)
        .def("print_header", &DFHelper::print_header)
        .def("add_transformation", &DFHelper::add_transformation, "name"_a, "key1"_a, "key2"_a, "order"_a = "Qpq")
        .def("transform", &DFHelper::transform)
        .def("clear_spaces", &DFHelper::clear_spaces)
        .def("clear_all", &DFHelper::clear_all)
        .def("transpose", &DFHelper::transpose)
        .def("get_space_size", &DFHelper::get_space_size)
        .def("get_tensor_size", &DFHelper::get_tensor_size)
        .def("get_tensor_shape", &DFHelper::get_tensor_shape)
        .def("get_tensor", take_string(&DFHelper::get_tensor))
        .def("get_tensor", tensor_access3(&DFHelper::get_tensor));

    py::class_<MemDFJK, std::shared_ptr<MemDFJK>, JK>(m, "MemDFJK", "docstring")
        .def("dfh", &MemDFJK::dfh, "Return the DFHelper object.");

    py::class_<DirectJK, std::shared_ptr<DirectJK>, JK>(m, "DirectJK", "docstring")
        .def("do_incfock_iter", &DirectJK::do_incfock_iter, "Was the last Fock build incremental?");

    py::class_<CompositeJK, std::shared_ptr<CompositeJK>, JK>(m, "CompositeJK", "docstring")
        .def("do_incfock_iter", &CompositeJK::do_incfock_iter, "Was the last Fock build incremental?")
        .def("clear_D_prev", &CompositeJK::clear_D_prev, "Clear previous D matrices.");

    py::class_<scf::SADGuess, std::shared_ptr<scf::SADGuess>>(m, "SADGuess", "docstring")
        .def_static("build_SAD",
                    [](std::shared_ptr<BasisSet> basis, std::vector<std::shared_ptr<BasisSet>> atomic_bases) { return scf::SADGuess(basis, atomic_bases, Process::environment.options); })
        .def("compute_guess", &scf::SADGuess::compute_guess)
        .def("set_print", &scf::SADGuess::set_print)
        .def("set_debug", &scf::SADGuess::set_debug)
        .def("set_atomic_fit_bases", &scf::SADGuess::set_atomic_fit_bases)
        .def("Da", &scf::SADGuess::Da)
        .def("Db", &scf::SADGuess::Db)
        .def("Ca", &scf::SADGuess::Ca)
        .def("Cb", &scf::SADGuess::Cb);
}
