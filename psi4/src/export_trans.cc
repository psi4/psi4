/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include <vector>

#include <pybind11/stl.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>

#include "psi4/libtrans/mospace.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libtrans/fcidump_helper.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsio/psio.hpp"

using namespace psi;
namespace py = pybind11;
using namespace pybind11::literals;

void export_trans(py::module& m) {
    py::class_<MOSpace, std::shared_ptr<MOSpace>>(m, "MOSpace",
                                                  "Defines orbital spaces in which to transform integrals")
        .def(py::init<const char>())
        .def(py::init<const char, const std::vector<int>, const std::vector<int>, const std::vector<int>,
                      const std::vector<int>>())
        .def(py::init<const char, const std::vector<int>, const std::vector<int>>())
        .def_static("fzc", []() { return MOSpace::fzc; })
        .def_static("occ", []() { return MOSpace::occ; })
        .def_static("fzv", []() { return MOSpace::fzv; })
        .def_static("vir", []() { return MOSpace::vir; })
        .def_static("all", []() { return MOSpace::all; })
        .def_static("nil", []() { return MOSpace::nil; })
        .def_static("dum", []() { return MOSpace::dum; })
        .def(py::self == py::self)
        .def(py::self >= py::self)
        .def(py::self != py::self)
        .def(py::self <= py::self)
        .def(py::self < py::self)
        .def(py::self > py::self)
        .def("label", &MOSpace::label, "Get the unique identifier for this space")
        .def("aOrbs", &MOSpace::aOrbs, "Get the alpha orbitals")
        .def("bOrbs", &MOSpace::bOrbs, "Get the beta orbitals")
        .def("aIndex", &MOSpace::aIndex, "Get the alpha orbital indexing array")
        .def("bIndex", &MOSpace::bIndex, "Get the beta orbital indexing array");

    py::class_<IntegralTransform, std::shared_ptr<IntegralTransform>> int_trans_bind(
        m, "IntegralTransform", "IntegralTransform transforms one- and two-electron integrals within general spaces",
        py::dynamic_attr());

    py::enum_<IntegralTransform::HalfTrans>(int_trans_bind, "HalfTrans")
        .value("MakeAndKeep", IntegralTransform::HalfTrans::MakeAndKeep)
        .value("ReadAndKeep", IntegralTransform::HalfTrans::ReadAndKeep)
        .value("MakeAndNuke", IntegralTransform::HalfTrans::MakeAndNuke)
        .value("ReadAndNuke", IntegralTransform::HalfTrans::ReadAndNuke);
    py::enum_<IntegralTransform::TransformationType>(int_trans_bind, "TransformationType")
        .value("Restricted", IntegralTransform::TransformationType::Restricted)
        .value("Unrestricted", IntegralTransform::TransformationType::Unrestricted)
        .value("SemiCanonical", IntegralTransform::TransformationType::SemiCanonical);
    py::enum_<IntegralTransform::MOOrdering>(int_trans_bind, "MOOrdering")
        .value("QTOrder", IntegralTransform::MOOrdering::QTOrder)
        .value("PitzerOrder", IntegralTransform::MOOrdering::PitzerOrder);
    py::enum_<IntegralTransform::OutputType>(int_trans_bind, "OutputType")
        .value("DPDOnly", IntegralTransform::OutputType::DPDOnly)
        .value("IWLOnly", IntegralTransform::OutputType::IWLOnly)
        .value("IWLAndDPD", IntegralTransform::OutputType::IWLAndDPD);
    py::enum_<IntegralTransform::FrozenOrbitals>(int_trans_bind, "FrozenOrbitals")
        .value("None", IntegralTransform::FrozenOrbitals::None)
        .value("OccOnly", IntegralTransform::FrozenOrbitals::OccOnly)
        .value("VirOnly", IntegralTransform::FrozenOrbitals::VirOnly)
        .value("OccAndVir", IntegralTransform::FrozenOrbitals::OccAndVir);
    py::enum_<IntegralTransform::SpinType>(int_trans_bind, "SpinType")
        .value("Alpha", IntegralTransform::SpinType::Alpha)
        .value("Beta", IntegralTransform::SpinType::Beta);

    int_trans_bind.def(py::init<std::shared_ptr<Wavefunction>, std::vector<std::shared_ptr<MOSpace>>,
                                IntegralTransform::TransformationType, IntegralTransform::OutputType,
                                IntegralTransform::MOOrdering, IntegralTransform::FrozenOrbitals, bool>(),
                       "wfn"_a, "spaces"_a, "transformationType"_a = IntegralTransform::TransformationType::Restricted,
                       "outputType"_a = IntegralTransform::OutputType::DPDOnly,
                       "moOrdering"_a = IntegralTransform::MOOrdering::QTOrder,
                       "FrozenOrbitals"_a = IntegralTransform::FrozenOrbitals::OccAndVir, "initialize"_a = true);

    int_trans_bind.def(py::init<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>, std::shared_ptr<Matrix>,
                                std::shared_ptr<Matrix>, std::shared_ptr<Matrix>, std::vector<std::shared_ptr<MOSpace>>,
                                IntegralTransform::TransformationType, IntegralTransform::OutputType,
                                IntegralTransform::MOOrdering, IntegralTransform::FrozenOrbitals, bool>(),
                       "H"_a, "c"_a, "i"_a, "a"_a, "v"_a, "spaces"_a,
                       "transformationType"_a = IntegralTransform::TransformationType::Restricted,
                       "outputType"_a = IntegralTransform::OutputType::DPDOnly,
                       "moOrdering"_a = IntegralTransform::MOOrdering::QTOrder,
                       "FrozenOrbitals"_a = IntegralTransform::FrozenOrbitals::OccAndVir, "initialize"_a = true);

    typedef int (IntegralTransform::*DPD_ID_1)(const std::string&);
    typedef int (IntegralTransform::*DPD_ID_2)(const char);
    typedef int (IntegralTransform::*DPD_ID_3)(const std::shared_ptr<MOSpace>, const std::shared_ptr<MOSpace>,
                                               IntegralTransform::SpinType, bool);

    int_trans_bind.def("initialize", &IntegralTransform::initialize, "Initialize an IntegralTransform")
        .def("presort_so_tei", &IntegralTransform::presort_so_tei, "docstring")
        .def("update_orbitals", &IntegralTransform::update_orbitals, "docstring")
        .def("transform_tei", &IntegralTransform::transform_tei, "Transform two-electron integrals", "s1"_a, "s2"_a,
             "s3"_a, "s4"_a, "half_trans"_a = IntegralTransform::HalfTrans::MakeAndNuke)
        .def("transform_tei_first_half", &IntegralTransform::transform_tei_first_half,
             "First half-transform two-electron integrals", "s1"_a, "s2"_a)
        .def("transform_tei_second_half", &IntegralTransform::transform_tei_second_half,
             "Second half-transform two-electron integrals", "s1"_a, "s2"_a, "s3"_a, "s4"_a)
        .def("backtransform_density", &IntegralTransform::backtransform_density)
        .def("backtransform_tpdm_restricted", &IntegralTransform::backtransform_tpdm_restricted)
        .def("backtransform_tpdm_unrestricted", &IntegralTransform::backtransform_tpdm_unrestricted)
        .def("print_dpd_lookup", &IntegralTransform::print_dpd_lookup)
        .def("compute_fock_like_matrices", &IntegralTransform::compute_fock_like_matrices, "Hcore"_a, "Cmats"_a)
        .def("DPD_ID", DPD_ID_1(&IntegralTransform::DPD_ID), "docstring", "c"_a)
        .def("DPD_ID", DPD_ID_2(&IntegralTransform::DPD_ID), "docstring", "str"_a)
        .def("DPD_ID", DPD_ID_3(&IntegralTransform::DPD_ID), "docstring", "s1"_a, "s2"_a, "spin"_a, "pack"_a);

    int_trans_bind.def("set_so_tei_file", &IntegralTransform::set_so_tei_file)
        .def("set_write_dpd_so_tpdm", &IntegralTransform::set_write_dpd_so_tpdm)
        .def("set_print", &IntegralTransform::set_print)
        .def("set_orbitals", &IntegralTransform::set_orbitals)
        .def("get_print", &IntegralTransform::get_print)
        .def("get_frozen_core_energy", &IntegralTransform::get_frozen_core_energy)
        .def("set_keep_ht_ints", &IntegralTransform::set_keep_ht_ints)
        .def("get_keep_ht_ints", &IntegralTransform::get_keep_ht_ints)
        .def("set_keep_dpd_so_ints", &IntegralTransform::set_keep_dpd_so_ints)
        .def("get_keep_dpd_so_ints", &IntegralTransform::get_keep_dpd_so_ints)
        .def("set_keep_iwl_so_ints", &IntegralTransform::set_keep_iwl_so_ints)
        .def("get_keep_iwl_so_ints", &IntegralTransform::get_keep_iwl_so_ints)
        .def("set_tpdm_already_presorted", &IntegralTransform::set_tpdm_already_presorted)
        .def("get_tei_already_presorted", &IntegralTransform::get_tei_already_presorted)
        .def("set_tei_already_presorted", &IntegralTransform::set_tei_already_presorted)
        .def("set_memory", &IntegralTransform::set_memory)
        .def("get_memory", &IntegralTransform::get_memory)
        .def("set_dpd_id", &IntegralTransform::set_dpd_id)
        .def("get_dpd_id", &IntegralTransform::get_dpd_id)
        .def("get_psio", &IntegralTransform::get_psio)
        .def("set_psio", &IntegralTransform::set_psio)
        .def("set_dpd_int_file", &IntegralTransform::set_dpd_int_file)
        .def("set_aa_int_name", &IntegralTransform::set_aa_int_name)
        .def("set_ab_int_name", &IntegralTransform::set_ab_int_name)
        .def("set_bb_int_name", &IntegralTransform::set_bb_int_name)
        .def("alpha_corr_to_pitzer", &IntegralTransform::alpha_corr_to_pitzer)
        .def("beta_corr_to_pitzer", &IntegralTransform::beta_corr_to_pitzer)
        .def("nirrep", &IntegralTransform::nirrep)
        .def("reset_so_int", &IntegralTransform::reset_so_int);

    m.def("fcidump_tei_helper", &fcidump::fcidump_tei_helper, "Write integrals to file in FCIDUMP format", "nirrep"_a,
          "restricted"_a, "DPD_info"_a, "ints_tolerance"_a, "fname"_a = "INTDUMP");
}
