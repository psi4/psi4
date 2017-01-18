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

#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"

using namespace psi;

void export_oeprop(py::module &m)
{
    py::class_<Prop, std::shared_ptr<Prop> >(m, "Prop", "docstring").
            def("set_title", &OEProp::set_title, "docstring");

    //     def(init<std::shared_ptr<Wavefunction> >()).
    //     def("print_header", pure_virtual(&Prop::print_header)).
    //     def("compute", pure_virtual(&Prop::compute)).
    //     def("set_Da_ao", &Prop::set_Da_ao, "docstring").
    //     def("set_Db_ao", &Prop::set_Db_ao, "docstring").
    //     def("set_Da_so", &Prop::set_Da_so, "docstring").
    //     def("set_Db_so", &Prop::set_Db_so, "docstring").
    //     def("set_Da_mo", &Prop::set_Da_mo, "docstring").
    //     def("set_Db_mo", &Prop::set_Db_mo, "docstring");

    py::class_<OEProp, std::shared_ptr<OEProp>, Prop>(m, "OEProp", "docstring").
            // TODO had no_init but init member present
            def(py::init<std::shared_ptr<Wavefunction> >()).
            def("add", &OEProp::oepy_add, "docstring").
            def("compute", &OEProp::oepy_compute, "docstring").
//        def("set_title", &OEProp::set_title, "docstring").
            def("clear", &OEProp::clear, "docstring").
            def("set_Da_ao", &OEProp::set_Da_ao, "docstring", py::arg("Da"), py::arg("symmetry") = 0).
            def("set_Db_ao", &OEProp::set_Db_ao, "docstring", py::arg("Db"), py::arg("symmetry") = 0).
            def("set_Da_so", &OEProp::set_Da_so, "docstring").
            def("set_Db_so", &OEProp::set_Db_so, "docstring").
            def("set_Da_mo", &OEProp::set_Da_mo, "docstring").
            def("set_Db_mo", &OEProp::set_Db_mo, "docstring").
            def("Vvals", &OEProp::Vvals, "The electrostatic potential (in a.u.) at each grid point").
            def("Exvals", &OEProp::Exvals, "The x component of the field (in a.u.) at each grid point").
            def("Eyvals", &OEProp::Eyvals, "The y component of the field (in a.u.) at each grid point").
            def("Ezvals", &OEProp::Ezvals, "The z component of the field (in a.u.) at each grid point");


    //class_<GridProp, std::shared_ptr<GridProp> >("GridProp", "docstring").
    //    def("add", &GridProp::gridpy_add, "docstring").
    //    def("set_filename", &GridProp::set_filename, "docstring").
    //    def("add_alpha_mo", &GridProp::add_alpha_mo, "docstring").
    //    def("add_beta_mo", &GridProp::add_beta_mo, "docstring").
    //    def("add_basis_fun", &GridProp::add_basis_fun, "docstring").
    //    def("build_grid_overages", &GridProp::build_grid_overages, "docstring").
    //    def("set_n", &GridProp::set_n, "docstring").
    //    def("set_o", &GridProp::set_o, "docstring").
    //    def("set_l", &GridProp::set_l, "docstring").
    //    def("get_n", &GridProp::get_n, "docstring").
    //    def("get_o", &GridProp::get_o, "docstring").
    //    def("get_l", &GridProp::get_l, "docstring").
    //    def("set_caxis", &GridProp::set_caxis, "docstring").
    //    def("set_format", &GridProp::set_format, "docstring").
    //    def("compute", &GridProp::gridpy_compute, "docstring");
}
