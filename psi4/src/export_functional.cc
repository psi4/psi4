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

#include "psi4/pybind11.h"
#include "psi4/libmints/vector.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libfunctional/functional.h"
#include "psi4/libfunctional/LibXCfunctional.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libdisp/dispersion.h"
#include "psi4/libfock/v.h"
#include "psi4/libfock/points.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libqt/qt.h"

using namespace psi;

void export_functional(py::module &m)
{
    py::class_<SuperFunctional, std::shared_ptr<SuperFunctional>>(m, "SuperFunctional", "docstring").
        // TODO add init
        def(py::init<>()).
        def_static("blank", &SuperFunctional::blank, "Initialize a blank SuperFunctional.").
        def_static("XC_build", &SuperFunctional::XC_build, "Builds a SuperFunctional from a XC string.").
        def("allocate", &SuperFunctional::allocate, "Allocates the vectors, should be called after ansatz or npoint changes.").
        def("compute_functional", &SuperFunctional::compute_functional, "Computes the SuperFunctional.").
        def("x_functional", &SuperFunctional::x_functional, "Returns the desired X Functional.").
        def("c_functional", &SuperFunctional::c_functional, "Returns the desired C Functional").
        def("add_x_functional", &SuperFunctional::add_x_functional, "Add a exchange Functional.").
        def("add_c_functional", &SuperFunctional::add_c_functional, "Add a correlation Functional.").
        def("test_functional", &SuperFunctional::test_functional, "Quick testing capabilities.").
        def("values", &SuperFunctional::values, "Return all internal values.").
        def("value", &SuperFunctional::value, "Returns a given internal value.").
        def("name", &SuperFunctional::name, "The name of the SuperFunctional.").
        def("description", &SuperFunctional::description, "The description of the SuperFunctional").
        def("citation", &SuperFunctional::citation, "SuperFunctional citation.").
        def("ansatz", &SuperFunctional::ansatz, "SuperFunctional rung.").
        def("max_points", &SuperFunctional::max_points, "Maximum number of grid points per block.").
        def("deriv", &SuperFunctional::deriv, "Maximum derivative to compute.").
        def("x_omega", &SuperFunctional::x_omega, "Range-seperated exchange parameter.").
        def("c_omega", &SuperFunctional::c_omega, "Range-seperated correlation parameter.").
        def("x_alpha", &SuperFunctional::x_alpha, "Amount of exact HF exchange.").
        def("x_beta", &SuperFunctional::x_beta, "Amount of exact HF exchange.").
        def("c_alpha", &SuperFunctional::c_alpha, "Amount of MP2 correlation.").
        def("vv10_b", &SuperFunctional::vv10_b, "The VV10 b parameter.").
        def("vv10_c", &SuperFunctional::vv10_c, "The VV10 c parameter.").
        def("grac_shift", &SuperFunctional::grac_shift, "Shift of the bulk potenital.").
        def("grac_alpha", &SuperFunctional::grac_alpha, "GRAC Alpha.").
        def("grac_beta", &SuperFunctional::grac_beta, "GRAC Beta.").
        def("is_gga", &SuperFunctional::is_gga, "Is this a GGA?").
        def("is_meta", &SuperFunctional::is_meta, "Is this a MGGA?").
        def("is_x_lrc", &SuperFunctional::is_x_lrc, "Contains range-seperated exchange?").
        def("is_c_lrc", &SuperFunctional::is_c_lrc, "Contains range-seperated correlation?").
        def("is_x_hybrid", &SuperFunctional::is_x_hybrid, "Requires exact exchange?").
        def("is_c_hybrid", &SuperFunctional::is_c_hybrid, "Requires MP2 correlation?").
        def("set_name", &SuperFunctional::set_name, "Sets the SuperFunctional name.").
        def("set_description", &SuperFunctional::set_description, "Sets the SuperFunctional description.").
        def("set_citation", &SuperFunctional::set_citation, "Sets the SuperFunctional citation.").
        def("set_max_points", &SuperFunctional::set_max_points, "Sets the maximum number of points.").
        def("set_deriv", &SuperFunctional::set_deriv, "Sets the derivative level.").
        def("set_lock", &SuperFunctional::set_lock, "Locks the functional to prevent changes.").
        def("set_x_omega", &SuperFunctional::set_x_omega, "Sets the range-seperation exchange parameter.").
        def("set_c_omega", &SuperFunctional::set_c_omega, "Sets the range-seperation correlation parameter.").
        def("set_x_alpha", &SuperFunctional::set_x_alpha, "Sets the amount of exact HF exchange.").
        def("set_c_alpha", &SuperFunctional::set_c_alpha, "Sets the amount of MP2 correlation.").
        def("set_vv10_b", &SuperFunctional::set_vv10_b, "Sets the VV10 b parameter.").
        def("set_vv10_c", &SuperFunctional::set_vv10_c, "Sets the VV10 c parameter.").
        def("set_grac_shift", &SuperFunctional::set_grac_shift, "Sets the GRAC bulk shift value.").
        def("set_grac_alpha", &SuperFunctional::set_grac_alpha, "Sets the GRAC alpha parameter.").
        def("set_grac_beta", &SuperFunctional::set_grac_beta, "Sets the GRAC beta parameter.").
        def("needs_xc", &SuperFunctional::needs_xc, "Does this functional need XC quantities.").
        def("needs_vv10", &SuperFunctional::needs_vv10, "Does this functional need VV10 dispersion.").
        def("needs_grac", &SuperFunctional::needs_grac, "Does this functional need GRAC.").
        def("print_out",&SuperFunctional::py_print, "Prints out functional details.").
        def("print_detail",&SuperFunctional::py_print_detail, "Prints all SuperFunctional information.");

    py::class_<Functional, std::shared_ptr<Functional> >(m, "Functional", "docstring").
        // TODO need init
        def_static("build_base", &Functional::build_base,
            py::arg("alias"), "docstring").
        def("compute_functional", &Functional::compute_functional, "docstring").
        def("name", &Functional::name, "docstring").
        def("description", &Functional::description, "docstring").
        def("citation", &Functional::citation, "docstring").
        def("alpha", &Functional::alpha, "docstring").
        def("omega", &Functional::omega, "docstring").
        def("lsda_cutoff", &Functional::lsda_cutoff, "docstring").
        def("meta_cutoff", &Functional::meta_cutoff, "docstring").
        def("is_gga", &Functional::is_gga, "docstring").
        def("is_meta", &Functional::is_meta, "docstring").
        def("is_lrc", &Functional::is_lrc, "docstring").
        def("set_name", &Functional::set_name, "docstring").
        def("set_description", &Functional::set_description, "docstring").
        def("set_citation", &Functional::set_citation, "docstring").
        def("set_gga", &Functional::set_gga, "docstring").
        def("set_meta", &Functional::set_meta, "docstring").
        def("set_alpha", &Functional::set_alpha, "docstring").
        def("set_omega", &Functional::set_omega, "docstring").
        def("set_lsda_cutoff", &Functional::set_lsda_cutoff, "docstring").
        def("set_meta_cutoff", &Functional::set_meta_cutoff, "docstring").
        def("set_parameter", &Functional::set_parameter, "docstring").
        def("print_out", &Functional::py_print, "docstring").
        def("print_detail",&Functional::py_print_detail, "docstring");

    py::class_<LibXCFunctional, std::shared_ptr<LibXCFunctional>, Functional>(m, "LibXCFunctional", "docstring").
        def(py::init<std::string, bool>()).
        def("get_mix_data", &LibXCFunctional::get_mix_data, "docstring").
        def("set_omega", &LibXCFunctional::set_omega, "docstring");

    py::class_<VBase, std::shared_ptr<VBase> >(m, "VBase", "docstring").
        def_static("build", [](std::shared_ptr<BasisSet> &basis, std::shared_ptr<SuperFunctional> &func, std::string type){
            return VBase::build_V(basis, func, Process::environment.options, type);
        }).
        def("initialize", &VBase::initialize, "doctsring").
        def("finalize", &VBase::finalize, "doctsring").
        def("set_D", &VBase::set_D, "Sets the internal density.").
        def("Dao", &VBase::set_D, "Returns internal AO density.").
        def("compute_V", &VBase::compute_V, "doctsring").
        def("compute_Vx", &VBase::compute_Vx, "doctsring").
        def("compute_gradient", &VBase::compute_gradient, "Compute the DFT nuclear gradient contribution.").
        def("compute_hessain", &VBase::compute_hessian, "Compute the DFT nuclear Hessian contribution.").
        def("initialize", &VBase::initialize, "Initializes the V object.").
        def("finalize", &VBase::finalize, "Finalizes the V object.").

        def("basis", &VBase::basis, "Returns the internal basis set.").
        def("functional", &VBase::functional, "Returns the interal superfunctional.").
        def("properties", &VBase::properties, "Returns the properties computer.").
        def("get_block", &VBase::get_block, "Returns the requested BlockOPoints.").
        def("nblocks", &VBase::nblocks, "Total number of blocks.").
        def("quadrature_values", &VBase::quadrature_values, "Returns the quadrature values.");

    py::class_<BasisFunctions, std::shared_ptr<BasisFunctions> >(m, "BasisFunctions", "docstring").
        def("max_functions", &BasisFunctions::max_functions, "docstring").
        def("max_points", &BasisFunctions::max_points, "docstring").
        def("deriv", &BasisFunctions::deriv, "docstring").
        def("set_deriv", &BasisFunctions::set_deriv, "docstring").
        def("compute_functions", &BasisFunctions::compute_functions, "docstring").
        def("basis_values", &BasisFunctions::basis_values, "docstring");

    py::class_<PointFunctions, std::shared_ptr<PointFunctions>, BasisFunctions>(m, "PointFunctions", "docstring").
        def("print_out", &PointFunctions::print, py::arg("OutFileRMR")="outfile", py::arg("print") = 2, "docstring").
        def("compute_points", &PointFunctions::compute_points, "docstring").
        def("point_values", &PointFunctions::point_values, "docstring").
        def("orbital_values", &PointFunctions::orbital_values, "docstring");

    py::class_<BlockOPoints, std::shared_ptr<BlockOPoints> >(m, "BlockOPoints", "docstring").
        def("x", [](BlockOPoints &grid){
            SharedVector ret = std::shared_ptr<Vector>(new Vector("X Grid points", grid.npoints()));
            C_DCOPY(grid.npoints(), grid.x(), 1, ret->pointer(), 1);
            return ret;
        }).
        def("y", [](BlockOPoints &grid){
            SharedVector ret = std::shared_ptr<Vector>(new Vector("Y Grid points", grid.npoints()));
            C_DCOPY(grid.npoints(), grid.y(), 1, ret->pointer(), 1);
            return ret;
        }).
        def("z", [](BlockOPoints &grid){
            SharedVector ret = std::shared_ptr<Vector>(new Vector("Z Grid points", grid.npoints()));
            C_DCOPY(grid.npoints(), grid.z(), 1, ret->pointer(), 1);
            return ret;
        }).
        def("w", [](BlockOPoints &grid){
            SharedVector ret = std::shared_ptr<Vector>(new Vector("Grid Weights", grid.npoints()));
            C_DCOPY(grid.npoints(), grid.w(), 1, ret->pointer(), 1);
            return ret;
        }).
        def("refresh", &BlockOPoints::refresh, "docstring").
        def("npoints", &BlockOPoints::npoints, "docstring").
        def("print_out", &BlockOPoints::print, py::arg("OutFileRMR")="outfile", py::arg("print") = 2, "docstring").
        def("shells_local_to_global", &BlockOPoints::shells_local_to_global, "docstring").
        def("functions_local_to_global", &BlockOPoints::functions_local_to_global, "docstring");


    py::class_<Dispersion, std::shared_ptr<Dispersion> >(m, "Dispersion", "docstring").
        // TODO need init
        def_static("build", &Dispersion::build,
            py::arg("type"), py::arg("s6")=0.0, py::arg("p1")=0.0, py::arg("p2")=0.0, py::arg("p3")=0.0, "docstring").
        def("name", &Dispersion::name, "docstring").
        def("description", &Dispersion::description, "docstring").
        def("citation", &Dispersion::citation, "docstring").
        def("bibtex", &Dispersion::bibtex, "Get the BibTeX key for the literature reference.").
        def("set_name", &Dispersion::set_name, "docstring").
        def("set_description", &Dispersion::set_description, "docstring").
        def("set_citation", &Dispersion::set_citation, "docstring").
        def("set_bibtex", &Dispersion::set_bibtex, "Set the BibTeX key for the literature reference.").
        def("print_energy", &Dispersion::print_energy, "docstring").
        def("print_gradient", &Dispersion::print_gradient, "docstring").
        def("print_hessian", &Dispersion::print_hessian, "docstring").
        def("compute_energy", &Dispersion::compute_energy, "docstring").
        def("compute_gradient", &Dispersion::compute_gradient, "docstring").
        def("compute_hessian", &Dispersion::compute_hessian, "docstring").
        def("d", &Dispersion::get_d, "docstring").
        def("s6", &Dispersion::get_s6, "docstring").
        def("sr6", &Dispersion::get_sr6, "docstring").
        def("s8", &Dispersion::get_s8, "docstring").
        def("a1", &Dispersion::get_a1, "docstring").
        def("a2", &Dispersion::get_a2, "docstring").
        def("print_out",&Dispersion::py_print, "docstring");

}
