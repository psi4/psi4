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

#include "psi4/libpsi4util/libpsi4util.h"
#include "options.h"

#include <pybind11/stl.h>

using namespace pybind11::literals;

namespace psi {

// Utility functions; internal use only
py::dict make_dict_entry(const std::string& type, const std::string& group,
                         py::object default_value, const std::string& description) {
    return py::dict("type"_a = type, "group"_a = py::str(group), "value"_a = default_value,
                    "default_value"_a = default_value, "description"_a = description.c_str(),
                    "validator"_a = py::none());
}

py::dict make_dict_entry(const std::string& type, const std::string& group,
                         py::object default_value, py::list allowed_values,
                         const std::string& description) {
    if (type != "str") {
        throw std::runtime_error("allowed_values are only sypported for str-type.");
    }
    return py::dict("type"_a = type, "group"_a = py::str(group), "value"_a = default_value,
                    "default_value"_a = default_value, "allowed_values"_a = allowed_values,
                    "description"_a = description.c_str(), "validator"_a = py::none());
}

void throw_options_get_error(const std::string& correct_type, const std::string& label,
                             const std::string& type) {
    std::string msg = "Called FOptions::get_" + correct_type + "(" + label +
                      ") but the type for this option is " + type;
    throw std::runtime_error(msg);
}

void check_options_none(py::object obj, const std::string& type, const std::string& label) {
    if (obj.is_none()) {
        std::string msg = 
            "Called FOptions::get_" + type + "(" + label + ") but the value is set to None";
        throw std::runtime_error(msg);
    }   
}  

void FOptions::set_group(const std::string& group) {
    group_ = to_upper_copy(group);
}   

bool FOptions::exists(const std::string& group, const std::string& label) const {
    const auto group_uc = to_upper_copy(group);
    const auto label_uc = to_upper_copy(label);
    auto tuple = py::make_tuple(group_uc, label_uc);
    return dict_.contains(tuple);
}   

void FOptions::add(const std::string& group, const std::string& label, const std::string& type, py::object default_value,
                   const std::string& description) {
    const auto group_uc = to_upper_copy(group);
    const auto label_uc = to_upper_copy(label);
    if (exists(group_uc, label_uc)) {
        std::string msg = "Called FOptions::add_" +type + "(" + label +
                          ") but this option name is already used.";
        throw std::runtime_error(msg);
    }
    auto tuple = py::make_tuple(group_uc, label_uc);
    dict_[tuple] = make_dict_entry(type, group_uc, default_value, description);
}

void FOptions::add(const std::string& group, const std::string& label, const std::string& type, py::object default_value,
                   py::list allowed_values, const std::string& description) {
    const auto group_uc = to_upper_copy(group);
    const auto label_uc = to_upper_copy(label);
    if (exists(group_uc, label_uc)) {
        std::string msg = "Called FOptions::add_" +type + "(" + label +
                          ") but this option name is already used.";
        throw std::runtime_error(msg);
    }
    auto tuple = py::make_tuple(group_uc, label_uc);
    dict_[tuple] =
        make_dict_entry(type, group_uc, default_value, allowed_values, description);
}

void FOptions::add_bool(const std::string& group, const std::string& label, py::object default_value,
                        const std::string& description) {
    if (default_value.is_none()) {
        add(group, label, "bool", py::none(), description);
    } else {
        add(group, label, "bool", py::bool_(default_value), description);
    }
}


void FOptions::add_int(const std::string& group, const std::string& label, py::object default_value,
                       const std::string& description) {
    if (default_value.is_none()) {
        add(group, label, "int", py::none(), description);
    } else {
        add(group, label, "int", py::int_(default_value), description);
    }                        
}

void FOptions::add_double(const std::string& group, const std::string& label, py::object default_value,
                          const std::string& description) {
    if (default_value.is_none()) {
        add(group, label, "float", py::none(), description);
    } else {
        add(group, label, "float", py::float_(default_value), description);
    }
} 

void FOptions::add_str(const std::string& group, const std::string& label, py::object default_value,
                       const std::string& description) {
    if (default_value.is_none()) {
        add(group, label, "str", py::none(), description);
    } else {
        auto upper_ver = to_upper_copy(py::cast<std::string>(default_value));
        add(group, label, "str", py::str(upper_ver), description);
    }           
}
                
void FOptions::add_str(const std::string& group, const std::string& label, py::object default_value,
                       const std::vector<std::string>& allowed_values,
                       const std::string& description) {
    auto allowed_values_list = py::list();
    for (const auto& s : allowed_values) {
        allowed_values_list.append(py::str(s)); 
    }               
    if (default_value.is_none()) {
        add(group, label, "str", py::none(), allowed_values_list, description);
    } else {    
        add(group, label, "str", py::str(default_value), allowed_values_list, description);
    }   
}      

void FOptions::add_array(const std::string& group, const std::string& label, const std::string& description) {
    add(group, label, "gen_list", py::list(), description);
}

void FOptions::add_int_array(const std::string& group, const std::string& label, const std::string& description) {
    add(group, label, "int_list", py::list(), description);
}

void FOptions::add_double_array(const std::string& group, const std::string& label, const std::string& description) {
    add(group, label, "float_list", py::list(), description);
}


std::pair<py::object, std::string> FOptions::get(const std::string& group, const std::string& label) const {
    auto label_uc = to_upper_copy(label);
    auto group_uc = to_upper_copy(group);
    auto result = std::make_pair(py::cast<py::object>(Py_None), std::string("None"));
    auto tuple = py::make_tuple(group_uc, label_uc);
    if (dict_.contains(tuple)) {
        auto dict_entry = dict_[tuple];
        result = std::make_pair(dict_entry["value"], py::cast<std::string>(dict_entry["type"]));
    } else {
        std::string msg = "Called FOptions::get(" + label + ") this option is not registered.";
        throw std::runtime_error(msg);
    }
    return result;
}

bool FOptions::get_bool(const std::string& group, const std::string& label) const {
    auto value_type = get(group, label);
    check_options_none(value_type.first, "bool", label);
    if (value_type.second == "bool") {
        return py::cast<bool>(value_type.first);
    }           
    throw_options_get_error("bool", label, value_type.second);
    return false;   
}   

int FOptions::get_int(const std::string& group, const std::string& label) const {
    auto value_type = get(group, label);
    check_options_none(value_type.first, "int", label);
    if (value_type.second == "int") {
        return py::cast<int>(value_type.first);
    }
    throw_options_get_error("int", label, value_type.second);
    return 0;
}

double FOptions::get_double(const std::string& group, const std::string& label) const {
    auto value_type = get(group, label);
    check_options_none(value_type.first, "double", label);
    if (value_type.second == "float") { 
        return py::cast<double>(value_type.first);
    }   
    throw_options_get_error("double", label, value_type.second);
    return 0.0; 
}

std::string FOptions::get_str(const std::string& group, const std::string& label) const {
    auto value_type = get(group, label);
    check_options_none(value_type.first, "str", label);
    if (value_type.second == "str") {
        return py::cast<std::string>(value_type.first);
    }   
    throw_options_get_error("str", label, value_type.second);
    return std::string();
}

py::list FOptions::get_gen_list(const std::string& group, const std::string& label) const {
    auto value_type = get(group, label);
    check_options_none(value_type.first, "gen_list", label);
    if (value_type.second == "gen_list") {
        return value_type.first;
    }
    throw_options_get_error("gen_list", label, value_type.second);
    return py::list();
}

std::vector<int> FOptions::get_int_list(const std::string& group, const std::string& label) const {
    std::vector<int> result;
    auto value_type = get(group, label);
    check_options_none(value_type.first, "int_vec", label);
    if (value_type.second == "int_list") {
        for (const auto& s : value_type.first) {
            result.push_back(py::cast<int>(s));
        }
        return result;
    }
    throw_options_get_error("int_vec", label, value_type.second);
    return result;
}

std::vector<double> FOptions::get_double_list(const std::string& group, const std::string& label) const {
    std::vector<double> result;
    auto value_type = get(group, label);
    check_options_none(value_type.first, "double_vec", label);
    if (value_type.second == "float_list") {
        for (const auto& s : value_type.first) {
            result.push_back(py::cast<double>(s));
        }
        return result;
    }
    throw_options_get_error("double_vec", label, value_type.second);

    return result;
}

void FOptions::set(const std::string& group, const std::string& label, const py::object val) {
    auto label_uc = to_upper_copy(label);
    auto group_uc = to_upper_copy(group);
    auto tuple = py::make_tuple(group_uc, label_uc);
    if (dict_.contains(tuple)) {
        auto dict_entry = dict_[tuple];
        // Validate allowed_values.
        if (dict_entry.contains("allowed_values")) {
             const auto type = py::cast<std::string>(dict_entry["type"]);
             if (type != "str") {
                 std::string msg = label + " has allowed_options, but it's not a str. This is not supported.";
                 throw std::runtime_error(msg);
             }
             const auto val_str = py::cast<std::string>(val);
             std::vector<std::string> allowed_values_vec;
             for (const auto& s : dict_entry["allowed_values"]) {
                 allowed_values_vec.push_back(py::str(s));
             }
             if (std::find(allowed_values_vec.begin(), allowed_values_vec.end(), val_str) == allowed_values_vec.end()) {
                 std::string choices_joined;
                 for (const auto& s : allowed_values_vec) {
                     if (!choices_joined.empty()) {
                         choices_joined.append(" ");
                     }
                     choices_joined.append(s);
                 };
                 std::string msg = "value " + val_str + " not in allowed_values.\n Please select one of "
                                   "the following values [" + choices_joined + "]";
                 throw std::runtime_error(msg);
             }
        }
        dict_entry["value"] = val;
    } else {              
        std::string msg = "Called FOptions::set(" + label + ") this option is not registered.";
        throw std::runtime_error(msg);
    }
}

void FOptions::set_bool(const std::string& group, const std::string& label, bool val) {
    auto value_type = get(group, label);
    if (value_type.second == "bool") {
        set(group, label, py::bool_(val));
    } else {
        std::string msg = "Called FOptions::set_bool(" + label +
                          ") but the type for this option is " + value_type.second;
        throw std::runtime_error(msg);
    }
}

void FOptions::set_int(const std::string& group, const std::string& label, int val) {
    auto value_type = get(group, label);
    if (value_type.second == "int") {
        set(group, label, py::int_(val));
    } else {
        std::string msg = "Called FOptions::set_int(" + label +
                          ") but the type for this option is " + value_type.second;
        throw std::runtime_error(msg);
    }
}   

void FOptions::set_double(const std::string& group, const std::string& label, double val) {
    auto value_type = get(group, label);
    if (value_type.second == "float") {
        set(group, label, py::float_(val));
    } else {
        std::string msg = "Called FOptions::set_double(" + label +
                          ") but the type for this option is " + value_type.second;
        throw std::runtime_error(msg);
    }
}

void FOptions::set_str(const std::string& group, const std::string& label, const std::string& val) {
    auto value_type = get(group, label);
    if (value_type.second == "str") {
        auto ucase_val = to_upper_copy(val);
        set(group, label, py::str(ucase_val));
    } else {
        std::string msg = "Called FOptions::set_str(" + label +
                          ") but the type for this option is " + value_type.second;
        throw std::runtime_error(msg);
    }
}  

void FOptions::set_gen_list(const std::string& group, const std::string& label, py::list val) {
    auto value_type = get(group, label);
    if (value_type.second == "gen_list") {
        set(group, label, val);
    } else {
        std::string msg = "Called ForteOptions::set_gen_list(" + label +
                          ") but the type for this option is " + value_type.second;
        throw std::runtime_error(msg);
    }
}

void FOptions::set_int_list(const std::string& group, const std::string& label, const std::vector<int>& val) {
    std::vector<int> result;
    auto value_type = get(group, label);
    if (value_type.second == "int_list") {
        set(group, label, py::cast(val));
    } else {
        std::string msg = "Called ForteOptions::get_int_list(" + label +
                          ") but the type for this option is " + value_type.second;
        throw std::runtime_error(msg);
    }
}

void FOptions::set_double_list(const std::string& group, const std::string& label, const std::vector<double>& val) {
    std::vector<double> result;
    auto value_type = get(group, label);
    if (value_type.second == "float_list") {
        set(group, label, py::cast(val));
    } else {
        std::string msg = "Called ForteOptions::set_double_list(" + label +
                          ") but the type for this option is " + value_type.second;
        throw std::runtime_error(msg);
    }
}

}  // namespace psi

