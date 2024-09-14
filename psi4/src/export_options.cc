/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "psi4/liboptions/liboptions.h"
#include "psi4/liboptions/options.h"
#include "psi4/pybind11.h"

#include <string>

using namespace psi;
using namespace pybind11::literals;

void export_foptions(py::module& m) {
    py::class_<FOptions>(m, "FOptions", "docstring", py::dynamic_attr())
        .def(py::init<>())
        .def("get_group", &FOptions::get_group, "Get the options group.")
        .def("set_group", &FOptions::set_group, "Set the options group.")
        .def("add_bool", &FOptions::add_bool, "Add a Boolean option.", "label"_a, "default_value"_a, "description"_a = "")
        .def("add_int", &FOptions::add_int, "Add an int option.", "label"_a, "default_value"_a, "description"_a = "")
        .def("add_double", &FOptions::add_double, "Add a double option.", "label"_a, "default_value"_a, "description"_a = "")
        .def("add_str",
             (void(FOptions::*)(const std::string&, py::object, const std::string&)) &
                 FOptions::add_str,
             "Add a string option", "label"_a, "default_value"_a, "description"_a = "")
        .def("add_str",
             (void(FOptions::*)(const std::string&, py::object, const std::vector<std::string>&,
                                    const std::string&)) &
                 FOptions::add_str,
             "Add a string option", "label"_a, "default_value"_a, "allowed_values"_a, "description"_a = "")
        .def("add_int_list", &FOptions::add_int_array, "Add a list of integers option", "label"_a, "description"_a = "")
        .def("add_double_list", &FOptions::add_double_array, "Add a list of doubles option", "label"_a, "description"_a = "")
        .def("add_list", &FOptions::add_array, "Add an array option for general elements", "label"_a, "description"_a = "")
        .def("get_bool", &FOptions::get_bool, "Get a Boolean option.")
        .def("get_int", &FOptions::get_int, "Get an int option.")
        .def("get_double", &FOptions::get_double, "Get a double option.")
        .def("get_str", &FOptions::get_str, "Get a string option.")
        .def("get_int_list", &FOptions::get_int_list, "Get a list of integers option")
        .def("get_double_list", &FOptions::get_double_list,
             "Get a vector of doubles (py::float) option")
        .def("get_list", &FOptions::get_gen_list, "Get a general list")
        .def("set_bool", &FOptions::set_bool, "Set a Boolean option.")
        .def("set_int", &FOptions::set_int, "Set an int option.")
        .def("set_double", &FOptions::set_double, "Set a double option.")
        .def("set_str", &FOptions::set_str, "Set a string option.")
        .def("set_int_list", &FOptions::set_int_list, "Set a vector of integers option")
        .def("set_double_list", &FOptions::set_double_list,
             "Set a vector of doubles (py::float) option")
        .def("set_list", &FOptions::set_gen_list,
             "Set a vector of python objects (py::object) option");
}
void export_options(py::module& m) {
    py::class_<Options>(m, "Options", "docstring", py::dynamic_attr())
        .def("add_bool", &Options::add_bool, "add bool option")
        .def("add_int", &Options::add_int, "add int option")
        .def("add_str", &Options::add_str, "add string option")
        .def("add_str_i", &Options::add_str_i, "add string option")
        .def("add_array", &Options::add_array, "add array option")
        .def("get_bool", &Options::get_bool, "get boolean option")
        .def("get_int", &Options::get_int, "get integer option")
        .def("get_double", &Options::get_double, "get double option")
        .def("get_str", &Options::get_str, "get string option")
        .def("get_str", &Options::get_str, "get string option")
        .def("get_int_vector", &Options::get_int_vector, "get int vector option")
        .def("set_bool", &Options::set_bool, "set bool option")
        .def("set_int", &Options::set_int, "set int option")
        .def("set_double", &Options::set_double, "set double option")
        .def("set_str", &Options::set_str, "set string option")
        .def("set_str_i", &Options::set_str_i, "set string option")
        .def("set_array", &Options::set_array, "set array option")
        .def("read_globals", &Options::read_globals, "expert")
        .def("set_read_globals", &Options::set_read_globals, "expert")
        .def("set_current_module", &Options::set_current_module, "sets *arg0* (all CAPS) as current module")
        .def("get_current_module", &Options::get_current_module, "gets current module")
        .def("validate_options", &Options::validate_options, "validate options for *arg0* module")
        .def("print_module_options", &Options::print, "print global and local options prepared for current module")
        .def("print_global_options", &Options::print_globals, "print the global, cross-module options");
}
