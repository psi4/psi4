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

#include "psi4/pybind11.h"
#include "psi4/libmints/vector.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libfunctional/functional.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libdisp/dispersion.h"
#include "psi4/libfock/v.h"
#include "psi4/libmints/basisset.h"

using namespace psi;

void export_functional(py::module &m)
{
    py::class_<SuperFunctional, std::shared_ptr<SuperFunctional>>(m, "SuperFunctional", "docstring").
        // TODO add init
        def(py::init<>()).
        def_static("blank", &SuperFunctional::blank, "docstring").
        def("allocate", &SuperFunctional::allocate, "docstring").
        def("x_functional", &SuperFunctional::x_functional, "docstring").
        def("c_functional", &SuperFunctional::c_functional, "docstring").
        def("add_x_functional", &SuperFunctional::add_x_functional, "docstring").
        def("add_c_functional", &SuperFunctional::add_c_functional, "docstring").
        def("test_functional", &SuperFunctional::test_functional, "docstring").
        def("value", &SuperFunctional::value, "docstring").
        def("name", &SuperFunctional::name, "docstring").
        def("description", &SuperFunctional::description, "docstring").
        def("citation", &SuperFunctional::citation, "docstring").
        def("ansatz", &SuperFunctional::ansatz, "docstring").
        def("max_points", &SuperFunctional::max_points, "docstring").
        def("deriv", &SuperFunctional::deriv, "docstring").
        def("x_omega", &SuperFunctional::x_omega, "docstring").
        def("c_omega", &SuperFunctional::c_omega, "docstring").
        def("x_alpha", &SuperFunctional::x_alpha, "docstring").
        def("c_alpha", &SuperFunctional::c_alpha, "docstring").
        def("c_ss_alpha", &SuperFunctional::c_ss_alpha, "docstring").
        def("c_os_alpha", &SuperFunctional::c_os_alpha, "docstring").
        def("is_gga", &SuperFunctional::is_gga, "docstring").
        def("is_meta", &SuperFunctional::is_meta, "docstring").
        def("is_x_lrc", &SuperFunctional::is_x_lrc, "docstring").
        def("is_c_lrc", &SuperFunctional::is_c_lrc, "docstring").
        def("is_x_hybrid", &SuperFunctional::is_x_hybrid, "docstring").
        def("is_c_hybrid", &SuperFunctional::is_c_hybrid, "docstring").
        def("is_c_scs_hybrid", &SuperFunctional::is_c_scs_hybrid, "docstring").
        def("set_name", &SuperFunctional::set_name, "docstring").
        def("set_description", &SuperFunctional::set_description, "docstring").
        def("set_citation", &SuperFunctional::set_citation, "docstring").
        def("set_max_points", &SuperFunctional::set_max_points, "docstring").
        def("set_deriv", &SuperFunctional::set_deriv, "docstring").
        def("set_x_omega", &SuperFunctional::set_x_omega, "docstring").
        def("set_c_omega", &SuperFunctional::set_c_omega, "docstring").
        def("set_x_alpha", &SuperFunctional::set_x_alpha, "docstring").
        def("set_c_alpha", &SuperFunctional::set_c_alpha, "docstring").
        def("set_c_ss_alpha", &SuperFunctional::set_c_ss_alpha, "docstring").
        def("set_c_os_alpha", &SuperFunctional::set_c_os_alpha, "docstring").
        def("print_out",&SuperFunctional::py_print, "docstring").
        def("print_detail",&SuperFunctional::py_print_detail, "docstring");

    py::class_<Functional, std::shared_ptr<Functional> >(m, "Functional", "docstring").
        // TODO need init
        def_static("build_base", &Functional::build_base,
            py::arg("alias"), "docstring").
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
        def("print_detail",&SuperFunctional::py_print_detail, "docstring");

    py::class_<VBase, std::shared_ptr<VBase> >(m, "VBase", "docstring").
        def_static("build", [](std::shared_ptr<BasisSet> &basis, std::shared_ptr<SuperFunctional> &func, std::string type){
            return VBase::build_V(basis, func, Process::environment.options, type);
        }).
        def("initialize", &VBase::initialize, "doctsring").
        def("finalize", &VBase::finalize, "doctsring").
        def("compute", &VBase::compute, "doctsring").
        def("compute_gradient", &VBase::compute_gradient, "doctsring").

        def("basis", &VBase::basis, "doctsring").
        def("functional", &VBase::functional, "doctsring").
        // def("properties", &VBase::properties, "doctsring").
        // def("grid", &VBase::grid, "doctsring").
        def("quadrature_values", &VBase::quadrature_values, "doctsring").

        def("C", &VBase::C, "doctsring").
        def("V", &VBase::V, "doctsring").
        def("D", &VBase::D, "doctsring").
        def("C_clear", [](VBase &v){
                v.C().clear();
            }).
        def("C_add", [](VBase &v, SharedMatrix C){
                v.C().push_back(C);
            });

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
