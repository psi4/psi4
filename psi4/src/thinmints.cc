/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#include <cstdio>
#include <sstream>
#include <map>
#include <iomanip>
#include <sys/stat.h>
#include "psi4/pybind11.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/matrix.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/liboptions/liboptions_python.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/process.h"

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/writer_file_prefix.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

// definitions for Fortran Runtime library init/finalize
extern "C" {
void for_rtl_init_(int*, char**);
int for_rtl_finish_();
int for_set_reentrancy(int*);
}

using namespace psi;

// Python helper wrappers
void export_blas_lapack(py::module&);
void export_mints(py::module&);
void export_oeprop(py::module&);
void export_wavefunction(py::module&);

// Forward declare /src/bin/ methods
namespace psi {

// Declare some globals
char* psi_file_prefix;
std::string outfile_name;
std::string restart_id;
std::shared_ptr<PsiOutStream> outfile;

extern int read_options(const std::string& name, Options& options, bool suppress_printing = false);
// extern void print_version(std::string);
}

std::string to_upper(const std::string& key) {
    std::string nonconst_key = key;
    std::transform(nonconst_key.begin(), nonconst_key.end(), nonconst_key.begin(), ::toupper);
    return nonconst_key;
}

void py_flush_outfile() {}

void py_close_outfile() {
    if (outfile) {
        outfile = std::shared_ptr<PsiOutStream>();
    }
}

void py_reopen_outfile() {
    if (outfile_name == "stdout") {
        // outfile = stdout;
    } else {
        outfile = std::make_shared<PsiOutStream>(outfile_name, std::ostream::app);
        if (!outfile) throw PSIEXCEPTION("Psi4: Unable to reopen output file.");
    }
}

void py_be_quiet() {
    py_close_outfile();
    outfile = std::make_shared<PsiOutStream>("/dev/null", std::ostream::app);
    if (!outfile) throw PSIEXCEPTION("Psi4: Unable to redirect output to /dev/null.");
}

std::string py_get_outfile_name() { return outfile_name; }

void py_psi_prepare_options_for_module(std::string const& name) {
    // Tell the options object which module is about to run
    Process::environment.options.set_current_module(name);
    // Figure out the defaults for any options that have not been specified
    read_options(name, Process::environment.options, false);
    // Now we've read in the defaults, make sure that user-specified options are recognized by the current module
    Process::environment.options.validate_options();
}

void py_psi_clean() { PSIOManager::shared_object()->psiclean(); }

void py_psi_print_options() { Process::environment.options.print(); }

void py_psi_print_global_options() { Process::environment.options.print_globals(); }

std::vector<std::string> py_psi_get_global_option_list() { return Process::environment.options.list_globals(); }

void py_psi_clean_options() {
    Process::environment.options.clear();
    Process::environment.options.set_read_globals(true);
    read_options("", Process::environment.options, true);
    Process::environment.options.set_read_globals(false);
}

void py_psi_print_out(std::string s) { (*outfile->stream()) << s; }

/**
 * @return whether key describes a convergence threshold or not
 */
bool specifies_convergence(std::string const& key) {
    return ((key.find("CONV") != key.npos) || (key.find("TOL") != key.npos));
}

Options& py_psi_get_options() { return Process::environment.options; }

bool py_psi_set_local_option_string(std::string const& module, std::string const& key, std::string const& value) {
    std::string nonconst_key = to_upper(key);
    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "string") {
        Process::environment.options.set_str(module, nonconst_key, value);
    } else if (data.type() == "istring") {
        Process::environment.options.set_str_i(module, nonconst_key, value);
    } else if (data.type() == "boolean") {
        if (to_upper(value) == "TRUE" || to_upper(value) == "YES" || to_upper(value) == "ON")
            Process::environment.options.set_bool(module, nonconst_key, true);
        else if (to_upper(value) == "FALSE" || to_upper(value) == "NO" || to_upper(value) == "OFF")
            Process::environment.options.set_bool(module, nonconst_key, false);
        else
            throw std::domain_error("Required option type is boolean, no boolean specified");
    }
    return true;
}

bool py_psi_set_local_option_int(std::string const& module, std::string const& key, int value) {
    std::string nonconst_key = to_upper(key);
    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "double" && specifies_convergence(nonconst_key)) {
        double val = pow(10.0, -value);
        Process::environment.options.set_double(module, nonconst_key, val);
    } else if (data.type() == "boolean") {
        Process::environment.options.set_bool(module, nonconst_key, value ? true : false);
    } else if (data.type() == "string" || data.type() == "istring") {
        Process::environment.options.set_str(module, nonconst_key, std::to_string(value));
    } else {
        Process::environment.options.set_int(module, nonconst_key, value);
    }
    return true;
}

bool py_psi_set_local_option_double(std::string const& module, std::string const& key, double value) {
    std::string nonconst_key = to_upper(key);

    Process::environment.options.set_double(module, nonconst_key, value);
    return true;
}

bool py_psi_set_global_option_string(std::string const& key, std::string const& value) {
    std::string nonconst_key = to_upper(key);
    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "string" || data.type() == "istring") {
        Process::environment.options.set_global_str(nonconst_key, value);
    } else if (data.type() == "boolean") {
        if (to_upper(value) == "TRUE" || to_upper(value) == "YES" || to_upper(value) == "ON")
            Process::environment.options.set_global_bool(nonconst_key, true);
        else if (to_upper(value) == "FALSE" || to_upper(value) == "NO" || to_upper(value) == "OFF")
            Process::environment.options.set_global_bool(nonconst_key, false);
        else
            throw std::domain_error("Required option type is boolean, no boolean specified");
    }
    return true;
}

bool py_psi_set_global_option_int(std::string const& key, int value) {
    std::string nonconst_key = to_upper(key);
    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "double" && specifies_convergence(nonconst_key)) {
        double val = pow(10.0, -value);
        Process::environment.options.set_global_double(nonconst_key, val);
    } else if (data.type() == "boolean") {
        Process::environment.options.set_global_bool(nonconst_key, value ? true : false);
    } else if (data.type() == "string" || data.type() == "istring") {
        Process::environment.options.set_global_str(nonconst_key, std::to_string(value));
    } else {
        Process::environment.options.set_global_int(nonconst_key, value);
    }
    return true;
}

bool py_psi_set_global_option_double(std::string const& key, double value) {
    std::string nonconst_key = to_upper(key);

    Process::environment.options.set_global_double(nonconst_key, value);
    return true;
}

bool py_psi_set_local_option_array(std::string const& module, std::string const& key, const py::list& values,
                                   DataType* entry = nullptr) {
    std::string nonconst_key = to_upper(key);
    // Assign a new head entry on the first time around only
    if (entry == nullptr) {
        // We just do a cheesy "get" to make sure keyword is valid.  This get will throw if not.
        Data& data = Process::environment.options[nonconst_key];
        // This "if" statement is really just here to make sure the compiler doesn't optimize out the get, above.
        if (data.type() == "array") Process::environment.options.set_array(module, nonconst_key);
    }
    size_t size = len(values);
    for (int n = 0; n < size; ++n) {
        if (py::isinstance<py::list>(values[n])) {
            py::list l = values[n].cast<py::list>();
            DataType* newentry = Process::environment.options.set_local_array_array(module, nonconst_key, entry);
            // Now we need to recurse, to fill in the data
            py_psi_set_local_option_array(module, key, l, newentry);
        } else {
            // This is not a list; try to cast to a string
            try {
                std::string s = values[n].cast<std::string>();
                Process::environment.options.set_local_array_string(module, nonconst_key, s, entry);
            } catch (py::cast_error e) {
                try {
                    // This is not a list or string; try to cast to an integer
                    int i = values[n].cast<int>();
                    Process::environment.options.set_local_array_int(module, nonconst_key, i, entry);
                } catch (py::cast_error e) {
                    // This had better be castable to a float.  We don't catch the exception here
                    // because if we encounter one, something bad has happened
                    double f = values[n].cast<double>();
                    Process::environment.options.set_local_array_double(module, nonconst_key, f, entry);
                }
            }
        }
    }
    return true;
}

bool py_psi_set_local_option_array_wrapper(std::string const& module, std::string const& key, py::list values) {
    // A wrapper to help pybind11 handle default values
    return py_psi_set_local_option_array(module, key, values);
}

bool py_psi_set_global_option_array(std::string const& key, py::list values, DataType* entry = nullptr) {
    std::string nonconst_key = to_upper(key);
    // Assign a new head entry on the first time around only
    if (entry == nullptr) {
        // We just do a cheesy "get" to make sure keyword is valid.  This get will throw if not.
        Data& data = Process::environment.options[nonconst_key];
        // This "if" statement is really just here to make sure the compiler doesn't optimize out the get, above.
        if (data.type() == "array") Process::environment.options.set_global_array(nonconst_key);
    }
    size_t size = len(values);
    for (int n = 0; n < size; ++n) {
        if (py::isinstance<py::list>(values[n])) {
            py::list l = values[n].cast<py::list>();
            DataType* newentry = Process::environment.options.set_global_array_array(nonconst_key, entry);
            // Now we need to recurse, to fill in the data
            py_psi_set_global_option_array(key, l, newentry);
        } else {
            // This is not a list; try to cast to a string
            try {
                std::string s = values[n].cast<std::string>();
                Process::environment.options.set_global_array_string(nonconst_key, s, entry);
            } catch (py::cast_error e) {
                try {
                    // This is not a list or string; try to cast to an integer
                    int i = values[n].cast<int>();
                    Process::environment.options.set_global_array_int(nonconst_key, i, entry);
                } catch (py::cast_error e) {
                    // This had better be castable to a float.  We don't catch the exception here
                    // because if we encounter one, something bad has happened
                    double f = values[n].cast<double>();
                    Process::environment.options.set_global_array_double(nonconst_key, f, entry);
                }
            }
        }
    }
    return true;
}

bool py_psi_set_global_option_array_wrapper(std::string const& key, py::list values) {
    // A wrapper to help pybind11 handle default values
    return py_psi_set_global_option_array(key, values);
}

void py_psi_set_local_option_python(const std::string& key, py::object& obj) {
    std::string nonconst_key = to_upper(key);
    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "python")
        dynamic_cast<PythonDataType*>(data.get())->assign(obj);
    else
        throw PSIEXCEPTION("Unable to set option to a Python object.");
}

bool py_psi_has_local_option_changed(std::string const& module, std::string const& key) {
    std::string nonconst_key = to_upper(key);
    Process::environment.options.set_current_module(module);
    py_psi_prepare_options_for_module(module);
    Data& data = Process::environment.options.get_local(nonconst_key);

    return data.has_changed();
}

bool py_psi_has_global_option_changed(std::string const& key) {
    std::string nonconst_key = to_upper(key);
    Data& data = Process::environment.options.get_global(nonconst_key);

    return data.has_changed();
}

bool py_psi_has_option_changed(std::string const& module, std::string const& key) {
    std::string nonconst_key = to_upper(key);
    Process::environment.options.set_current_module(module);
    py_psi_prepare_options_for_module(module);
    Data& data = Process::environment.options.use_local(nonconst_key);

    return data.has_changed();
}

bool py_psi_option_exists_in_module(std::string const& module, std::string const& key) {
    std::string nonconst_key = to_upper(key);
    Process::environment.options.set_current_module(module);
    py_psi_prepare_options_for_module(module);
    bool in_module = Process::environment.options.exists_in_active(nonconst_key);

    return in_module;
}

void py_psi_revoke_global_option_changed(std::string const& key) {
    std::string nonconst_key = to_upper(key);
    Data& data = Process::environment.options.get_global(nonconst_key);
    data.dechanged();
}

void py_psi_revoke_local_option_changed(std::string const& module, std::string const& key) {
    std::string nonconst_key = to_upper(key);
    Process::environment.options.set_current_module(module);
    py_psi_prepare_options_for_module(module);
    Data& data = Process::environment.options.get_local(nonconst_key);
    data.dechanged();
}

// Quick function that unpacks a data type
py::list data_to_list(py::list l, Data d) {
    if (d.is_array()) {
        // Recurse
        py::list row;
        for (int i = 0; i < d.size(); ++i) {
            data_to_list(row, d[i]);
        }
        l.append(row);
    } else if (d.type() == "double") {
        l.append(py::float_(d.to_double()));
    } else if (d.type() == "string") {
        l.append(py::str(d.to_string()));
    } else if (d.type() == "boolean") {
        l.append(py::bool_(d.to_integer()));
    } else if (d.type() == "int") {
        l.append(py::int_(d.to_integer()));
    } else {
        throw PSIEXCEPTION("Unknown data type in fill_list");
    }
    return l;
}

py::object py_psi_get_local_option(std::string const& module, std::string const& key) {
    std::string nonconst_key = to_upper(key);
    Process::environment.options.set_current_module(module);
    py_psi_prepare_options_for_module(module);
    Data& data = Process::environment.options.get_local(nonconst_key);

    if (data.type() == "string" || data.type() == "istring")
        return py::cast(data.to_string());
    else if (data.type() == "boolean" || data.type() == "int")
        return py::cast(data.to_integer());
    else if (data.type() == "double")
        return py::cast(data.to_double());
    else if (data.type() == "array") {
        py::list l;
        for (size_t i = 0; i < data.size(); i++) {
            data_to_list(l, data[i]);
        }
        return l;
    }

    return py::object();
}

py::object py_psi_get_global_option(std::string const& key) {
    std::string nonconst_key = to_upper(key);
    Data& data = Process::environment.options.get_global(nonconst_key);

    if (data.type() == "string" || data.type() == "istring")
        return py::cast(data.to_string());
    else if (data.type() == "boolean" || data.type() == "int")
        return py::cast(data.to_integer());
    else if (data.type() == "double")
        return py::cast(data.to_double());
    else if (data.type() == "array") {
        py::list l;
        for (size_t i = 0; i < data.size(); i++) {
            data_to_list(l, data[i]);
        }
        return l;
    }

    return py::object();
}

py::object py_psi_get_option(std::string const& module, std::string const& key) {
    std::string nonconst_key = to_upper(key);
    Process::environment.options.set_current_module(module);
    py_psi_prepare_options_for_module(module);
    Data& data = Process::environment.options.use_local(nonconst_key);

    if (data.type() == "string" || data.type() == "istring")
        return py::cast(data.to_string());
    else if (data.type() == "boolean" || data.type() == "int")
        return py::cast(data.to_integer());
    else if (data.type() == "double")
        return py::cast(data.to_double());
    else if (data.type() == "array") {
        py::list l;
        for (size_t i = 0; i < data.size(); i++) {
            data_to_list(l, data[i]);
        }
        return l;
    }

    return py::object();
}

void py_psi_set_active_molecule(std::shared_ptr<Molecule> molecule) { Process::environment.set_molecule(molecule); }
std::shared_ptr<Molecule> py_psi_get_active_molecule() { return Process::environment.molecule(); }

void py_psi_set_gradient(SharedMatrix grad) { Process::environment.set_gradient(grad); }

SharedMatrix py_psi_get_gradient() { return Process::environment.gradient(); }

std::shared_ptr<Vector> py_psi_get_atomic_point_charges() {
    auto empty = std::make_shared<psi::Vector>();
    return empty;  // charges not added to process.h for environment - yet(?)
}

bool py_psi_has_variable(const std::string& key) { return Process::environment.globals.count(to_upper(key)); }

double py_psi_get_variable(const std::string& key) { return Process::environment.globals[to_upper(key)]; }

SharedMatrix py_psi_get_array_variable(const std::string& key) { return Process::environment.arrays[to_upper(key)]; }

void py_psi_set_variable(const std::string& key, double val) { Process::environment.globals[to_upper(key)] = val; }

void py_psi_set_array_variable(const std::string& key, SharedMatrix val) {
    Process::environment.arrays[to_upper(key)] = val;
}

void py_psi_clean_variable_map() {
    Process::environment.globals.clear();
    Process::environment.arrays.clear();
}

void py_psi_set_memory(size_t mem, bool quiet) {
    Process::environment.set_memory(mem);
    if (!quiet) {
        outfile->Printf("\n  Memory set to %7.3f %s by Python driver.\n",
                        (mem > 1073741824 ? mem / 1073741824.0 : mem / 1048576.0), (mem > 1073741824 ? "GiB" : "MiB"));
    }
}

size_t py_psi_get_memory() { return Process::environment.get_memory(); }

void py_psi_set_n_threads(size_t nthread, bool quiet) {
#ifdef _OPENMP
    Process::environment.set_n_threads(nthread);
    if (!quiet) {
        outfile->Printf("  Threads set to %d by Python driver.\n", nthread);
    }
#else
    Process::environment.set_n_threads(1);
    if (!quiet) {
        outfile->Printf(
            "  Python driver attempted to set threads to %d.\n"
            "  Psi4 was compiled without OpenMP, setting threads to 1.\n",
            nthread);
    }
#endif
}

int py_psi_get_n_threads() { return Process::environment.get_n_threads(); }

void py_psi_print_variable_map() {
    int largest_key = 0;
    for (std::map<std::string, double>::iterator it = Process::environment.globals.begin();
         it != Process::environment.globals.end(); ++it) {
        if (it->first.size() > largest_key) largest_key = it->first.size();
    }
    largest_key += 2;  // for quotation marks

    std::stringstream line;
    std::string first_tmp;
    for (std::map<std::string, double>::iterator it = Process::environment.globals.begin();
         it != Process::environment.globals.end(); ++it) {
        first_tmp = "\"" + it->first + "\"";
        line << "  " << std::left << std::setw(largest_key) << first_tmp << " => " << std::setw(20) << std::right
             << std::fixed << std::setprecision(12) << it->second << std::endl;
    }

    outfile->Printf("\n\n  Variable Map:");
    outfile->Printf("\n  ----------------------------------------------------------------------------\n");
    outfile->Printf("%s\n\n", line.str().c_str());
}

std::map<std::string, double> py_psi_return_variable_map() { return Process::environment.globals; }

std::map<std::string, SharedMatrix> py_psi_return_array_variable_map() { return Process::environment.arrays; }

std::string py_psi_top_srcdir() { return TOSTRING(PSI_TOP_SRCDIR); }

bool psi4_python_module_initialize() {
    static bool initialized = false;

    if (initialized) {
        printf("Psi4 already initialized.\n");
        return true;
    }

    // Setup the environment
    Process::environment.initialize();  // Defaults to obtaining the environment from the global environ variable
    Process::environment.set_memory(524288000);

    // There should only be one of these in Psi4
    Wavefunction::initialize_singletons();

    outfile = std::make_shared<PsiOutStream>();
    outfile_name = "stdout";
    std::string fprefix = PSI_DEFAULT_FILE_PREFIX;
    psi_file_prefix = strdup(fprefix.c_str());

    // There is only one timer:
    timer_init();

    // Initialize the I/O library
    psio_init();

    // Setup globals options
    Process::environment.options.set_read_globals(true);
    read_options("", Process::environment.options, true);
    Process::environment.options.set_read_globals(false);

#ifdef INTEL_Fortran_ENABLED
    static int argc = 1;
    static char* argv = (char*)"";
    for_rtl_init_(&argc, &argv);
#endif

    initialized = true;

    return true;
}

void psi4_python_module_finalize() {
#ifdef INTEL_Fortran_ENABLED
    for_rtl_finish_();
#endif

    // Shut things down:
    // There is only one timer:
    timer_done();

    outfile = std::shared_ptr<PsiOutStream>();
    psi_file_prefix = nullptr;
}

PYBIND11_MODULE(thinmints, core) {
    core.doc() = "Libmints Innards of C++ Innards of Psi4: Open-Source Quantum Chemistry";
    //    py::module core("core", R"pbdoc(
    //
    //        Psi4: An Open-Source Ab Initio Electronic Structure Package
    //        -----------------------------------------------------------
    //
    //        .. currentmodule:: core
    //
    //        .. autosummary::
    //           :toctree: _generate
    //
    //           version
    //           clean
    //           set_local_option
    //)pbdoc");

    core.def("initialize", &psi4_python_module_initialize);
    core.def("finalize", &psi4_python_module_finalize);

    core.def("clean", py_psi_clean, "Function to remove scratch files. Call between independent jobs.");
    core.def("clean_options", py_psi_clean_options, "Function to reset options to clean state.");

    // BLAS/LAPACK Static Wrappers
    export_blas_lapack(core);

    // OEProp/GridProp
    export_oeprop(core);

    // Options
    core.def("prepare_options_for_module", py_psi_prepare_options_for_module,
             "Sets the options module up to return options pertaining to the named argument (e.g. SCF).");
    core.def("set_active_molecule", py_psi_set_active_molecule,
             "Activates a previously defined (in the input) molecule, by name.");
    core.def("get_active_molecule", &py_psi_get_active_molecule, "Returns the currently active molecule object.");
    core.def("get_gradient", py_psi_get_gradient,
             "Returns the most recently computed gradient, as a N by 3 :py:class:`~psi4.core.Matrix` object.");
    core.def("set_gradient", py_psi_set_gradient,
             "Assigns the global gradient to the values stored in the N by 3 Matrix argument.");
    core.def("get_atomic_point_charges", py_psi_get_atomic_point_charges,
             "Returns the most recently computed atomic point charges, as a double * object.");
    core.def("set_memory_bytes", py_psi_set_memory, py::arg("memory"), py::arg("quiet") = false,
             "Sets the memory available to Psi (in bytes).");
    core.def("get_memory", py_psi_get_memory, "Returns the amount of memory available to Psi (in bytes).");
    core.def("set_datadir", [](const std::string& pdd) { Process::environment.set_datadir(pdd); },
             "Returns the amount of memory available to Psi (in bytes).");
    core.def("get_datadir", []() { return Process::environment.get_datadir(); },
             "Sets the path to shared text resources, PSIDATADIR");
    core.def("set_num_threads", py_psi_set_n_threads, py::arg("nthread"), py::arg("quiet") = false,
             "Sets the number of threads to use in SMP parallel computations.");
    core.def("get_num_threads", py_psi_get_n_threads,
             "Returns the number of threads to use in SMP parallel computations.");
    //    core.def("mol_from_file",&LibBabel::ParseFile,"Reads a molecule from another input file");
    core.def("print_options", py_psi_print_options,
             "Prints the currently set options (to the output file) for the current module.");
    core.def("print_global_options", py_psi_print_global_options,
             "Prints the currently set global (all modules) options to the output file.");
    core.def("print_out", py_psi_print_out, "Prints a string (using sprintf-like notation) to the output file.");

    // Set the different local option types
    core.def("set_local_option", py_psi_set_local_option_array_wrapper,
             "Sets value *arg3* to array keyword *arg2* scoped only to a specific module *arg1*.");
    core.def("set_local_option", py_psi_set_local_option_int,
             "Sets value *arg3* to integer keyword *arg2* scoped only to a specific module *arg1*.");
    core.def("set_local_option", py_psi_set_local_option_double,
             "Sets value *arg3* to double keyword *arg2* scoped only to a specific module *arg1*.");
    core.def("set_local_option", py_psi_set_local_option_string,
             "Sets value *arg3* to string keyword *arg2* scoped only to a specific module *arg1*.");
    core.def("set_local_option_python", py_psi_set_local_option_python,
             "Sets an option to a Python object, but scoped only to a single module.");

    // Set the different global option types
    core.def("set_global_option", py_psi_set_global_option_array_wrapper,
             "Sets value *arg2* to array keyword *arg1* for all modules.");
    core.def("set_global_option", py_psi_set_global_option_int,
             "Sets value *arg2* to integer keyword *arg1* for all modules.");
    core.def("set_global_option", py_psi_set_global_option_double,
             "Sets value *arg2* to double keyword *arg1* for all modules.");
    core.def("set_global_option", py_psi_set_global_option_string,
             "Sets value *arg2* to string keyword *arg1* for all modules.");

    // Print options list
    core.def("get_global_option_list", py_psi_get_global_option_list, "Returns a list of all global options.");

    // Get the option; either global or local or let liboptions decide whether to use global or local
    core.def("get_global_option", py_psi_get_global_option,
             "Given a string of a keyword name *arg1*, returns the value associated with the keyword from the global "
             "options. Returns error if keyword is not recognized.");
    core.def("get_local_option", py_psi_get_local_option,
             "Given a string of a keyword name *arg2* and a particular module *arg1*, returns the value associated "
             "with the keyword in the module options scope. Returns error if keyword is not recognized for the "
             "module.");
    core.def("get_option", py_psi_get_option,
             "Given a string of a keyword name *arg2* and a particular module *arg1*, returns the local value "
             "associated with the keyword if it's been set, else the global value if it's been set, else the local "
             "core.default value. Returns error if keyword is not recognized globally or if keyword is not recognized "
             "for the module.");

    // Returns whether the option has changed/revoke has changed for silent resets
    core.def("has_global_option_changed", py_psi_has_global_option_changed,
             "Returns boolean for whether the keyword *arg1* has been touched in the global scope, by either user or "
             "code. Notwithstanding, code is written such that in practice, this returns whether the option has been "
             "touched in the global scope by the user.");
    core.def("has_local_option_changed", py_psi_has_local_option_changed,
             "Returns boolean for whether the keyword *arg2* has been touched in the scope of the specified module "
             "*arg1*, by either user or code. Notwithstanding, code is written such that in practice, this returns "
             "whether the option has been touched in the module scope by the user.");
    core.def("has_option_changed", py_psi_has_option_changed,
             "Returns boolean for whether the option *arg2* has been touched either locally to the specified module "
             "*arg1* or globally, by either user or code. Notwithstanding, code is written such that in practice, this "
             "returns whether the option has been touched by the user.");
    core.def("revoke_global_option_changed", py_psi_revoke_global_option_changed,
             "Given a string of a keyword name *arg1*, sets the has_changed attribute in the global options scope to "
             "false. Used in python driver when a function sets the value of an option. Before the function exits, "
             "this command is called on the option so that has_changed reflects whether the user (not the program) has "
             "touched the option.");
    core.def("revoke_local_option_changed", py_psi_revoke_local_option_changed,
             "Given a string of a keyword name *arg2* and a particular module *arg1*, sets the has_changed attribute "
             "in the module options scope to false. Used in python driver when a function sets the value of an option. "
             "Before the function exits, this command is called on the option so that has_changed reflects whether the "
             "user (not the program) has touched the option.");
    core.def("option_exists_in_module", py_psi_option_exists_in_module,
             "Given a string of a keyword name *arg1* and a particular module *arg0*, returns whether *arg1* is a "
             "valid option for *arg0*.");

    // These return/set/print PSI variables found in Process::environment.globals
    core.def("has_variable", py_psi_has_variable, "Returns true if the PSI variable exists/is set.");
    core.def("get_variable", py_psi_get_variable,
             "Returns one of the PSI variables set internally by the modules or python driver (see manual for full "
             "listing of variables available).");
    core.def("set_variable", py_psi_set_variable, "Sets a PSI variable, by name.");
    core.def("print_variables", py_psi_print_variable_map, "Prints all PSI variables that have been set internally.");
    core.def("clean_variables", py_psi_clean_variable_map, "Empties all PSI variables that have set internally.");
    core.def("get_variables", py_psi_return_variable_map,
             "Returns dictionary of the PSI variables set internally by the modules or python driver.");
    core.def("get_array_variable", py_psi_get_array_variable,
             "Returns one of the PSI variables set internally by the modules or python driver (see manual for full "
             "listing of variables available).");
    core.def("set_array_variable", py_psi_set_array_variable, "Sets a PSI variable, by name.");
    core.def("get_array_variables", py_psi_return_array_variable_map,
             "Returns dictionary of the PSI variables set internally by the modules or python driver.");

    // Returns the location where the Psi4 source is located.
    core.def("psi_top_srcdir", py_psi_top_srcdir, "Returns the location of the source code.");

    core.def("flush_outfile", py_flush_outfile, "Flushes the output file.");
    core.def("close_outfile", py_close_outfile, "Closes the output file.");
    core.def("reopen_outfile", py_reopen_outfile, "Reopens the output file.");
    core.def("outfile_name", py_get_outfile_name, "Returns the name of the output file.");
    core.def("be_quiet", py_be_quiet,
             "Redirects output to /dev/null.  To switch back to regular output mode, use reopen_outfile()");

    core.def("get_options", py_psi_get_options, py::return_value_policy::reference, "Get options");
    core.def("set_output_file", [](const std::string ofname) {
        outfile = std::make_shared<PsiOutStream>(ofname, std::ostream::trunc);
        outfile_name = ofname;
    });
    core.def("set_output_file", [](const std::string ofname, bool append) {
        outfile = std::make_shared<PsiOutStream>(ofname, (append ? std::ostream::app : std::ostream::trunc));
        outfile_name = ofname;
    });
    core.def("get_output_file", []() { return outfile_name; });
    core.def("set_psi_file_prefix", [](std::string fprefix) { psi_file_prefix = strdup(fprefix.c_str()); });

    // Define library classes
    export_mints(core);
}
