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

#include "psi4/liboptions/liboptions.h"

#include <string>

using namespace psi;


void export_options(py::module& m)
{
    py::class_<Options>(m, "Options", "docstring")
    .def("get_bool", &Options::get_bool, "get boolean option")
    .def("get_int", &Options::get_int, "get integer option");
    .def("get_double", &Options::get_double, "get double option")
    .def("get_str", &Options::get_str, "get string option")
    .def("get_str", &Options::get_str, "get string option")
    .def("get_int_vector", &Options::get_int_vector, "get int vector")
    .def("add_bool", &Options::add_bool, "add bool option")
    .def("add_bool", &Options::add_bool, "add bool option");
}

/*
 *
 *   void add_bool(std::string key, bool b);
    void add_int(std::string key, int i);
    void add_double(std::string key, double d);
    void add_str(std::string key, std::string s, std::string c = "");
    void add_str_i(std::string key, std::string s, std::string c = "");
    void add_array(std::string key);
    void set_bool(const std::string &module, const std::string &key, bool b);
    void set_int(const std::string &module, const std::string &key, int i);
    void set_double(const std::string & module, const std::string &key, double d);
    void set_str(const std::string & module, const std::string &key, std::string s);
    void set_str_i(const std::string & module, const std::string &key, std::string s);
    void set_python(const std::string &module, const std::string& key, const py::object &p);
    void set_array(const std::string &module, const std::string& key);

    bool get_bool(std::string key);
    int get_int(std::string key);
    double get_double(std::string key);
    std::string get_str(std::string key);
    int* get_int_array(std::string key);
    void fill_int_array(std::string key, int* empty_array);
    std::vector<int> get_int_vector(std::string key);
    double* get_double_array(std::string key);

*/
