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

#ifndef _psi_src_lib_liboptions_options_hpp
#define _psi_src_lib_liboptions_options_hpp

#include <pybind11/pybind11.h>

#include <string>
#include <vector>

namespace py = pybind11;

namespace psi {

class FOptions {
  public:
    /**
     * @brief Set the group to which options are added
     * @param group a string with the group name (default = "")
     */
    void set_group(const std::string& group = "");

    /**
     * @brief Get the group to which options are added
     * @return the group name
     */
    const std::string& get_group() const { return group_; };

    /**
     * @brief Check if an options exists
     * @param group Group label
     * @param label Option label
     * @return does this option exist?
     */
    bool exists(const std::string& group, const std::string& label) const;

    /**
     * @brief Add a python object option
     * @param group Group label
     * @param label Option label
     * @param type the option type
     * @param default_value default value of the option
     * @param description description of the option
     */
    void add(const std::string& group, const std::string& label, const std::string& type, pybind11::object default_value,
             const std::string& description);

    /**
     * @brief Add a python object option
     * @param label Option label
     * @param type the option type
     * @param default_value default value of the option
     * @param description description of the option
     */
    void add(const std::string& label, const std::string& type, pybind11::object default_value,
             const std::string& description) { add(group_, type, default_value, description) ;};

    /**
     * @brief Add a python object option
     * @param label Option label
     * @param type the option type
     * @param default_value default value of the option
     * @param description description of the option
     */
    void add(const std::string& group, const std::string& label, const std::string& type, pybind11::object default_value,
             pybind11::list allowed_values, const std::string& description);

    /**
     * @brief Add a python object option
     * @param group Group label
     * @param label Option label
     * @param type the option type
     * @param default_value default value of the option
     * @param description description of the option
     */
    void add(const std::string& label, const std::string& type, pybind11::object default_value,
             pybind11::list allowed_values, const std::string& description) { add(group_, label, type, default_value, allowed_values, description); };

    /**
     * @brief Add a boolean option
     * @param group Group label
     * @param label Option label
     * @param value Default value of the option
     * @param description Description of the option
     */
    void add_bool(const std::string& group, const std::string& label, py::object default_value,
                  const std::string& description = "");
    /**
     * @brief Add a boolean option
     * @param label Option label
     * @param value Default value of the option
     * @param description Description of the option
     */
    void add_bool(const std::string& label, py::object default_value,
                  const std::string& description = "") { add_bool(group_, label, default_value, description);  };

    /**
     * @brief Add an integer option
     * @param group Group label
     * @param label Option label
     * @param value Default value of the option
     * @param description Description of the option
     */
    void add_int(const std::string& group, const std::string& label, py::object default_value,
                 const std::string& description = "");
    /**
     * @brief Add an integer option
     * @param label Option label
     * @param value Default value of the option
     * @param description Description of the option
     */
    void add_int(const std::string& label, py::object default_value,
                 const std::string& description = "") { add_int(group_, label, default_value, description); };

    /**
     * @brief Add a double option
     * @param group Group label
     * @param label Option label
     * @param value Default value of the option
     * @param description Description of the option
     */
    void add_double(const std::string& group, const std::string& label, py::object default_value,
                    const std::string& description = "");
    /**
     * @brief Add a double option
     * @param label Option label
     * @param value Default value of the option
     * @param description Description of the option
     */
    void add_double(const std::string& label, py::object default_value,
                    const std::string& description = "") { add_double(group_, label, default_value, description); };

    /**
     * @brief Add a string option
     * @param group Group label
     * @param label Option label
     * @param value Default value of the option
     * @param description Description of the option
     */
    void add_str(const std::string& group, const std::string& label, py::object default_value,
                 const std::string& description = "");
    /**
     * @brief Add a string option
     * @param label Option label
     * @param value Default value of the option
     * @param description Description of the option
     */
    void add_str(const std::string& label, py::object default_value,
                 const std::string& description = "") { add_str(group_, label, default_value, description); };
    /**
     * @brief Add a string option and provide a list of allowed option values
     * @param group Group label
     * @param label Option label
     * @param value Default value of the option
     * @param description Description of the option
     * @param allowed_values An array of allowed option values
     */
    // TODO: Can I have allowed_values default to empty?
    void add_str(const std::string& group, const std::string& label, py::object default_value,
                 const std::vector<std::string>& allowed_values,
                 const std::string& description = "");
    /**
     * @brief Add a string option and provide a list of allowed option values
     * @param label Option label
     * @param value Default value of the option
     * @param description Description of the option
     * @param allowed_values An array of allowed option values
     */
    // TODO: Can I have allowed_values default to empty?
    void add_str(const std::string& label, py::object default_value,
                 const std::vector<std::string>& allowed_values,
                 const std::string& description = "") { add_str(group_, label, default_value, allowed_values, description); };

    /**
     * @brief Add a general array option
     * @param group Group label
     * @param label Option label
     * @param description Description of the option
     */
    void add_array(const std::string& group, const std::string& label, const std::string& description = "");
    /**
     * @brief Add a general array option
     * @param label Option label
     * @param description Description of the option
     */
    void add_array(const std::string& label, const std::string& description = "") { add_array(group_, label, description); };

    /**
     * @brief Add a integer array option
     * @param group Group label
     * @param label Option label
     * @param description Description of the option
     */
    void add_int_array(const std::string& group, const std::string& label, const std::string& description = "");
    /**
     * @brief Add a integer array option
     * @param label Option label
     * @param description Description of the option
     */
    void add_int_array(const std::string& label, const std::string& description = "") { add_int_array(group_, label, description) ; };

    /**
     * @brief Add a double array option
     * @param group Group label
     * @param label Option label
     * @param description Description of the option
     */
    void add_double_array(const std::string& group, const std::string& label, const std::string& description = "");
    /**
     * @brief Add a double array option
     * @param label Option label
     * @param description Description of the option
     */
    void add_double_array(const std::string& label, const std::string& description = "") { add_double_array(group_, label, description); };

    /**
     * @brief Get a python object option
     * @param group Group label
     * @param label Option label
     * @return a py::object containing the result
     */
    std::pair<py::object, std::string> get(const std::string& group, const std::string& label) const;

    /**
     * @brief Get a boolean option
     * @param group Group label
     * @param label Option label
     */
    bool get_bool(const std::string& group, const std::string& label) const;
    /**
     * @brief Get a boolean option
     * @param label Option label
     */
    bool get_bool(const std::string& label) const { return get_bool(group_, label); };

    /**          
     * @brief Get a integer option
     * @param group Group label
     * @param label Option label
     */
    int get_int(const std::string& group, const std::string& label) const;
    /**          
     * @brief Get a integer option
     * @param label Option label
     */
    int get_int(const std::string& label) const { return get_int(group_, label); };

    /**
     * @brief Get a double option
     * @param group Group label
     * @param label Option label
     */
    double get_double(const std::string& group, const std::string& label) const;
    /**
     * @brief Get a double option
     * @param label Option label
     */
    double get_double(const std::string& label) const { return get_double(group_, label); };

    /**
     * @brief Get a string option
     * @param group Group label
     * @param label Option label
     */
    std::string get_str(const std::string& group, const std::string& label) const;
    /**
     * @brief Get a string option
     * @param label Option label
     */
    std::string get_str(const std::string& label) const { return get_str(group_, label); };

    /**
     * @brief Get a general python list
     * @param group Group label
     * @param label
     * @return a py list
     */
    py::list get_gen_list(const std::string& group, const std::string& label) const;
    /**
     * @brief Get a general python list
     * @param label
     * @return a py list
     */
    py::list get_gen_list(const std::string& label) const { return get_gen_list(group_, label); };

    /**
     * @brief Get a vector of int option
     * @param group Group label
     * @param label Option label
     */
    std::vector<int> get_int_list(const std::string& group, const std::string& label) const;
    /**
     * @brief Get a vector of int option
     * @param label Option label
     */
    std::vector<int> get_int_list(const std::string& label) const { return get_int_list(group_, label); };

    /**
     * @brief Get a vector of int option
     * @param group Group label
     * @param label Option label
     */
    std::vector<double> get_double_list(const std::string& group, const std::string& label) const;
    /**
     * @brief Get a vector of int option
     * @param label Option label
     */
    std::vector<double> get_double_list(const std::string& label) const { return get_double_list(group_, label); };

    /**
     * @brief Set a python object option
     * @param group Group label
     * @param label Option label
     * @param a py::object containing the value to be set
     */
    void set(const std::string& group, const std::string& label, const py::object val);

    /**
     * @brief Set a boolean option
     * @param group Group label
     * @param label Option label
     * @param val Option value
     */
    void set_bool(const std::string& group, const std::string& label, bool val);
    /**
     * @brief Set a boolean option
     * @param label Option label
     * @param val Option value
     */
    void set_bool(const std::string& label, bool val) { set_bool(group_, label, val); };

    /**
     * @brief Set a integer option
     * @param group Group label
     * @param label Option label
     * @param val Option value
     */
    void set_int(const std::string& group, const std::string& label, int val);
    /**
     * @brief Set a integer option
     * @param label Option label
     * @param val Option value
     */
    void set_int(const std::string& label, int val) { set_int(group_, label, val); };

    /**
     * @brief Set a double option
     * @param group Group label
     * @param label Option label
     * @param val Option value
     */
    void set_double(const std::string& group, const std::string& label, double val);
    /**
     * @brief Set a double option
     * @param label Option label
     * @param val Option value
     */
    void set_double(const std::string& label, double val) { set_double(group_, label, val); } ;

    /**
     * @brief Set a string option
     * @param group Group label
     * @param label Option label
     * @param val Option value
     */
    void set_str(const std::string& group, const std::string& label, const std::string& val);
    /**
     * @brief Set a string option
     * @param label Option label
     * @param val Option value
     */
    void set_str(const std::string& label, const std::string& val) { set_str(group_, label, val); };

    /**
     * @brief Set a general python list
     * @param group Group label
     * @param label Option label
     * @param val Option value (a python list)
     */
    void set_gen_list(const std::string& group, const std::string& label, py::list val);
    /**
     * @brief Set a general python list
     * @param label Option label
     * @param val Option value (a python list)
     */
    void set_gen_list(const std::string& label, py::list val) { set_gen_list(group_, label, val); };

    /**
     * @brief Set a vector of int option
     * @param group Group label
     * @param label Option label
     * @param val Option value
     */
    void set_int_list(const std::string& group, const std::string& label, const std::vector<int>& val);
    /**
     * @brief Set a vector of int option
     * @param label Option label
     * @param val Option value
     */
    void set_int_list(const std::string& label, const std::vector<int>& val) { set_int_list(group_, label, val); } ;

    /**
     * @brief Set a vector of int option
     * @param group Group label
     * @param label Option label
     * @param val Option value
     */
    void set_double_list(const std::string& group, const std::string& label, const std::vector<double>& val);
    /**
     * @brief Set a vector of int option
     * @param label Option label
     * @param val Option value
     */
    void set_double_list(const std::string& label, const std::vector<double>& val) { set_double_list(group_, label, val); } ;

  private:
    /// a python dictionary object
    pybind11::dict dict_;
    /// the current option group
    std::string group_ = "";
};
}  // namespace psi
#endif

