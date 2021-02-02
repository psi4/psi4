/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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
#include <iomanip>
#include <map>
#include <sstream>
#include <sys/stat.h>

#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
//#include "psi4/pybind11.h"

#include "psi4/cc/cclambda/cclambda.h"
#include "psi4/cc/ccwave.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/writer_file_prefix.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libplugin/plugin.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"

//#include "python_data_type.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

// definitions for Fortran Runtime library init/finalize
extern "C" {
void for_rtl_init_(int*, char**);
int for_rtl_finish_();
int for_set_reentrancy(int*);
}

using namespace psi;
//namespace py = pybind11;
//using namespace pybind11::literals;

#ifdef USING_BrianQC
#include <use_brian_wrapper.h>
#include <brian_module.h>
#include <brian_macros.h>

bool brianEnableEnvFound = false;
bool brianEnableEnvValue = false;

BrianCookie brianCookie = 0;

bool brianEnable = false;
bool brianEnableDFT = true;
brianInt brianRestrictionType = 0;
bool brianCPHFFlag = false;
bool brianCPHFLeftSideFlag = false;

void checkBrian() {
    brianBool err = brianAPIGetError(&brianCookie);
    if (err) {
        throw PSIEXCEPTION("BrianQC error detected");
    }
}

void brianInit() {
    if (brianCookie != 0) {
        throw PSIEXCEPTION("Attempting to reinitialize BrianQC without releasing it first");
    }
    
    brianInt brianAPIVersionOfHost = BRIAN_API_VERSION;
    brianInt hostID = BRIAN_HOST_PSI4;
    brianCookie = brianAPIInit(&brianAPIVersionOfHost, &hostID);
    checkBrian();
    outfile->Printf("BrianQC initialization successful\n");
}

void brianRelease()
{
    if (brianCookie == 0) {
        throw PSIEXCEPTION("Attempting to release the BrianQC module when it hasn't been initialized\n");
    }
    
    outfile->Printf("Releasing the BrianQC module\n");
    brianAPIRelease(&brianCookie);
    brianCookie = 0;
}

void handleBrianOption(bool value) {
    if (brianEnableEnvFound) {
        outfile->Printf("BRIANQC_ENABLE option found, but overridden by BRIANQC_ENABLE environment variable\n");
    } else {
        outfile->Printf("BRIANQC_ENABLE option found, checking value\n");
        brianEnable = value;
        if (value && (brianCookie == 0)) {
            outfile->Printf("BRIANQC_ENABLE option set to true, initializing BrianQC\n");
            brianInit();
        } else if (!value && (brianCookie != 0)) {
            outfile->Printf("BRIANQC_ENABLE option set to false, releasing BrianQC\n");
            brianRelease();
        }
    }
}
#endif

// In export_plugins.cc
void py_psi_plugin_close_all();

extern std::map<std::string, plugin_info> plugins;

namespace opt {
psi::PsiReturnType optking(psi::Options&);
void opt_clean();
}  // namespace opt
// Forward declare /src/bin/ methods
namespace psi {

// Declare some globals
char* psi_file_prefix;
std::string outfile_name;
std::string restart_id;
std::shared_ptr<PsiOutStream> outfile;

// Wavefunction returns
namespace adc {
SharedWavefunction adc(SharedWavefunction, Options&);
}
namespace dct {
SharedWavefunction dct(SharedWavefunction, Options&);
}
namespace detci {
SharedWavefunction detci(SharedWavefunction, Options&);
}
namespace dfmp2 {
SharedWavefunction dfmp2(SharedWavefunction, Options&);
}
namespace dfoccwave {
SharedWavefunction dfoccwave(SharedWavefunction, Options&);
}
namespace libfock {
SharedWavefunction libfock(SharedWavefunction, Options&);
}
namespace fnocc {
SharedWavefunction fnocc(SharedWavefunction, Options&);
}
namespace occwave {
SharedWavefunction occwave(SharedWavefunction, Options&);
}
namespace mcscf {
SharedWavefunction mcscf(SharedWavefunction, Options&);
}
namespace psimrcc {
SharedWavefunction psimrcc(SharedWavefunction, Options&);
}

#ifdef USING_gdma
namespace gdma_interface {
SharedWavefunction gdma_interface(SharedWavefunction, Options&, const std::string& datfilename);
}
#endif

// Matrix returns
namespace scfgrad {
SharedMatrix scfgrad(SharedWavefunction, Options&);
}
namespace scfgrad {
SharedMatrix scfhess(SharedWavefunction, Options&);
}

// Does not create a wavefunction
// namespace fisapt { PsiReturnType fisapt(SharedWavefunction, Options&); }
namespace sapt {
PsiReturnType sapt(SharedWavefunction, SharedWavefunction, SharedWavefunction, Options&);
}

#ifdef USING_CheMPS2
namespace dmrg {
SharedWavefunction dmrg(SharedWavefunction, Options&);
}
#endif

// CC functions
namespace cctransort {
PsiReturnType cctransort(SharedWavefunction, Options&);
}
namespace cctriples {
PsiReturnType cctriples(SharedWavefunction, Options&);
}
namespace cchbar {
PsiReturnType cchbar(SharedWavefunction, Options&);
}
namespace cclambda {
PsiReturnType cclambda(SharedWavefunction, Options&);
}
namespace ccdensity {
PsiReturnType ccdensity(SharedWavefunction, Options&);
}
namespace ccresponse {
PsiReturnType ccresponse(SharedWavefunction, Options&);
void scatter(std::shared_ptr<Molecule> molecule, Options&, double step, std::vector<SharedMatrix> dip,
             std::vector<SharedMatrix> rot, std::vector<SharedMatrix> quad);
}  // namespace ccresponse
namespace cceom {
PsiReturnType cceom(SharedWavefunction, Options&);
}

extern int read_options(const std::string& name, Options& options, bool suppress_printing = false);
// extern void print_version(std::string);
}  // namespace psi

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
        // Default constructor corresponds to stdout
        outfile = std::make_shared<PsiOutStream>();
    } else {
        auto mode = std::ostream::app;
        outfile = std::make_shared<PsiOutStream>(outfile_name, mode);
        if (!outfile) throw PSIEXCEPTION("Psi4: Unable to reopen output file.");
    }
}

void py_be_quiet() {
    py_close_outfile();
    auto mode = std::ostream::app;
    outfile = std::make_shared<PsiOutStream>("/dev/null", mode);
    if (!outfile) throw PSIEXCEPTION("Psi4: Unable to redirect output to /dev/null.");
}

std::string py_get_outfile_name() { return outfile_name; }

void py_psi_prepare_options_for_module(std::string const& name) {
    // Tell the options object which module is about to run
    Process::environment.options.set_current_module(name);
    // Figure out the defaults for any options that have not been specified
    read_options(name, Process::environment.options, false);
    if (plugins.count(name)) {
        // Easy reference
        plugin_info& info = plugins[name];

        // Tell the plugin to load in its options into the current environment.
        info.read_options(info.name, Process::environment.options);
    }
    // Now we've read in the defaults, make sure that user-specified options are recognized by the current module
    Process::environment.options.validate_options();
}

int py_psi_optking() {
    py_psi_prepare_options_for_module("OPTKING");
    return opt::optking(Process::environment.options);
}

void py_psi_opt_clean(void) { opt::opt_clean(); }

// int py_psi_mints()
// {
//     py_psi_prepare_options_for_module("MINTS");
//     return mints::mints(Process::environment.options);
// }

SharedMatrix py_psi_scfgrad(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("SCF");
    return scfgrad::scfgrad(ref_wfn, Process::environment.options);
}

SharedMatrix py_psi_scfhess(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("SCF");
    return scfgrad::scfhess(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_occ(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("OCC");
    return occwave::occwave(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_dfocc(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("DFOCC");
    return dfoccwave::dfoccwave(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_mcscf(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("MCSCF");
    return mcscf::mcscf(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_dct(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("DCT");
    return dct::dct(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_dfmp2(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("DFMP2");
    return dfmp2::dfmp2(ref_wfn, Process::environment.options);
}

double py_psi_sapt(SharedWavefunction Dimer, SharedWavefunction MonomerA, SharedWavefunction MonomerB) {
    py_psi_prepare_options_for_module("SAPT");
    if (sapt::sapt(Dimer, MonomerA, MonomerB, Process::environment.options) == Success) {
        return Process::environment.globals["SAPT ENERGY"];
    } else
        return 0.0;
}

void py_psi_cctransort(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("CCTRANSORT");
    cctransort::cctransort(ref_wfn, Process::environment.options);
    // return 0.0;
}
SharedWavefunction py_psi_ccenergy(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("CCENERGY");
    auto ccwave = std::make_shared<ccenergy::CCEnergyWavefunction>(ref_wfn, Process::environment.options);

    double energy = ccwave->compute_energy();
    return ccwave;
}

double py_psi_cctriples(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("CCTRIPLES");
    if (cctriples::cctriples(ref_wfn, Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    } else
        return 0.0;
}

SharedWavefunction py_psi_fnocc(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("FNOCC");
    return fnocc::fnocc(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_detci(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("DETCI");
    return detci::detci(ref_wfn, Process::environment.options);
}

#ifdef USING_gdma
double py_psi_gdma(SharedWavefunction ref_wfn, const std::string& datfilename) {
    py_psi_prepare_options_for_module("GDMA");
    gdma_interface::gdma_interface(ref_wfn, Process::environment.options, datfilename);
    return 0.0;
}
#else
double py_psi_gdma(SharedWavefunction ref_wfn, const std::string& datfilename) {
    throw PSIEXCEPTION("GDMA not enabled. Recompile with -DENABLE_gdma.");
}
#endif

#ifdef USING_CheMPS2
SharedWavefunction py_psi_dmrg(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("DMRG");
    return dmrg::dmrg(ref_wfn, Process::environment.options);
}
#else
double py_psi_dmrg(SharedWavefunction ref_wfn) {
    throw PSIEXCEPTION("DMRG not enabled. Recompile with -DENABLE_CheMPS2");
}
#endif

void py_psi_cchbar(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("CCHBAR");
    cchbar::cchbar(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_cclambda(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("CCLAMBDA");
    std::shared_ptr<Wavefunction> cclambda =
        std::make_shared<cclambda::CCLambdaWavefunction>(ref_wfn, Process::environment.options);

    double energy = cclambda->compute_energy();
    return cclambda;
}

double py_psi_ccdensity(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("CCDENSITY");
    ccdensity::ccdensity(ref_wfn, Process::environment.options);
    return 0.0;
}

double py_psi_ccresponse(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("CCRESPONSE");
    ccresponse::ccresponse(ref_wfn, Process::environment.options);
    return 0.0;
}

//void py_psi_scatter(std::shared_ptr<Molecule> molecule, double step, py::list dip_polar_list, py::list opt_rot_list,
//                    py::list dip_quad_polar_list) {
//    py_psi_prepare_options_for_module("CCRESPONSE");
//
//    // Convert python tensor lists into vectors of sharedmatrices
//    std::vector<SharedMatrix> dip_polar_tensors;
//    std::vector<SharedMatrix> opt_rot_tensors;
//    std::vector<SharedMatrix> dip_quad_polar_tensors;
//
//    int list_len = len(dip_polar_list);
//    for (int i = 0; i < list_len; ++i) {
//        py::list dip_list = dip_polar_list[i].cast<py::list>();
//        py::list rot_list = opt_rot_list[i].cast<py::list>();
//        py::list quad_list = dip_quad_polar_list[i].cast<py::list>();
//        auto dip_mat = std::make_shared<Matrix>(3, 3);
//        auto rot_mat = std::make_shared<Matrix>(3, 3);
//        auto quad_mat = std::make_shared<Matrix>(9, 3);
//        for (int row = 0, j = 0; row < 3; ++row) {
//            for (int col = 0; col < 3; ++col, ++j) {
//                dip_mat->set(row, col, dip_list[j].cast<double>());
//                rot_mat->set(row, col, rot_list[j].cast<double>());
//            }
//        }
//        for (int row = 0, j = 0; row < 9; ++row) {
//            for (int col = 0; col < 3; ++col, ++j) {
//                quad_mat->set(row, col, quad_list[j].cast<double>());
//            }
//        }
//        dip_polar_tensors.push_back(dip_mat);
//        opt_rot_tensors.push_back(rot_mat);
//        dip_quad_polar_tensors.push_back(quad_mat);
//    }
//
//    //    for(std::vector<SharedMatrix>::iterator i=dip_polar_tensors.begin(); i != dip_polar_tensors.end(); ++i)
//    //        (*i)->print(stdout);
//    //    for(std::vector<SharedMatrix>::iterator i=opt_rot_tensors.begin(); i != opt_rot_tensors.end(); ++i)
//    //        (*i)->print(stdout);
//    //    for(std::vector<SharedMatrix>::iterator i=dip_quad_polar_tensors.begin(); i != dip_quad_polar_tensors.end();
//    //    ++i)
//    //        (*i)->print(stdout);
//
//    ccresponse::scatter(molecule, Process::environment.options, step, dip_polar_tensors, opt_rot_tensors,
//                        dip_quad_polar_tensors);
//}

double py_psi_cceom(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("CCEOM");
    if (cceom::cceom(ref_wfn, Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    } else
        return 0.0;
}

SharedWavefunction py_psi_psimrcc(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("PSIMRCC");
    return psimrcc::psimrcc(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_adc(SharedWavefunction ref_wfn) {
    py_psi_prepare_options_for_module("ADC");
    SharedWavefunction adc_wfn = adc::adc(ref_wfn, Process::environment.options);
    return adc_wfn;
}

char const* py_psi_version() {
#ifdef PSI_VERSION
    return PSI_VERSION;
#else
    return "";
#endif
}

char const* py_psi_git_version() {
#ifdef GIT_VERSION
    return GIT_VERSION;
#else
    return "";
#endif
}

void py_psi_clean() { PSIOManager::shared_object()->psiclean(); }

void py_psi_print_options() { Process::environment.options.print(); }

void py_psi_print_global_options() { Process::environment.options.print_globals(); }

std::vector<std::string> py_psi_get_global_option_list() { return Process::environment.options.list_globals(); }

void py_psi_clean_options() {
    Process::environment.options.clear();
    Process::environment.options.set_read_globals(true);
    read_options("", Process::environment.options, true);
    for (std::map<std::string, plugin_info>::iterator it = plugins.begin(); it != plugins.end(); ++it) {
        // Get the plugin options back into the global space
        it->second.read_options(it->second.name, Process::environment.options);
    }
    Process::environment.options.set_read_globals(false);
}

void py_psi_print_out(std::string s) { (*outfile->stream()) << s << std::flush; }

/**
 * @return whether key describes a convergence threshold or not
 */
bool specifies_convergence(std::string const& key) {
    return ((key.find("CONV") != key.npos) || (key.find("TOL") != key.npos));
}

// DCFT deprecation errors first added in 1.4. Feel free to retire after "enough" time.
void throw_deprecation_errors(std::string const& key, std::string const& module = "") {
    if (module == "DCFT") {
        throw PsiException(
            "Rename local options block. All instances of 'dcft' should be replaced with 'dct'. The method was renamed "
            "in v1.4.",
            __FILE__, __LINE__);
    }
    if (key.find("DCFT") != std::string::npos) {
        throw PsiException(
            "Rename keyword " + key +
                ". All instances of 'dcft' should be replaced with 'dct'. The method was renamed in v1.4.",
            __FILE__, __LINE__);
    }
}

Options& py_psi_get_options() { return Process::environment.options; }

bool py_psi_set_local_option_string(std::string const& module, std::string const& key, std::string const& value) {
    std::string nonconst_key = to_upper(key);

    throw_deprecation_errors(key, module);

    std::string module_temp = Process::environment.options.get_current_module();
    Process::environment.options.set_current_module(module);
    Data& data = Process::environment.options[nonconst_key];
    Process::environment.options.set_current_module(module_temp);

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

    throw_deprecation_errors(key, module);

    std::string module_temp = Process::environment.options.get_current_module();
    Process::environment.options.set_current_module(module);
    Data& data = Process::environment.options[nonconst_key];
    Process::environment.options.set_current_module(module_temp);

    if (data.type() == "double") {
        double val = (specifies_convergence(nonconst_key)) ? pow(10.0, -value) : double(value);
        Process::environment.options.set_double(module, nonconst_key, val);
    } else if (data.type() == "boolean") {
        Process::environment.options.set_bool(module, nonconst_key, value ? true : false);
    } else if (data.type() == "string" || data.type() == "istring") {
        Process::environment.options.set_str(module, nonconst_key, std::to_string(value));
    } else if (data.type() == "array") {
        Process::environment.options.set_local_array_int(module, nonconst_key, value, nullptr);
    } else {
        Process::environment.options.set_int(module, nonconst_key, value);
    }
    return true;
}

bool py_psi_set_local_option_double(std::string const& module, std::string const& key, double value) {
    std::string nonconst_key = to_upper(key);

    throw_deprecation_errors(key, module);

    Process::environment.options.set_double(module, nonconst_key, value);
    return true;
}

bool py_psi_set_global_option_string(std::string const& key, std::string const& value) {
    std::string nonconst_key = to_upper(key);

    throw_deprecation_errors(key);

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
    
#ifdef USING_BrianQC
    if (nonconst_key == "BRIANQC_ENABLE") {
        if (to_upper(value) == "TRUE" || to_upper(value) == "YES" || to_upper(value) == "ON")
            handleBrianOption(true);
        else if (to_upper(value) == "FALSE" || to_upper(value) == "NO" || to_upper(value) == "OFF")
            handleBrianOption(false);
        else
            throw std::domain_error("Required option type is boolean, no boolean specified");
    }
#endif

    return true;
}

bool py_psi_set_global_option_int(std::string const& key, int value) {
    std::string nonconst_key = to_upper(key);

    throw_deprecation_errors(key);

    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "double") {
        double val = (specifies_convergence(nonconst_key)) ? pow(10.0, -value) : double(value);
        Process::environment.options.set_global_double(nonconst_key, val);
    } else if (data.type() == "boolean") {
        Process::environment.options.set_global_bool(nonconst_key, value ? true : false);
    } else if (data.type() == "string" || data.type() == "istring") {
        Process::environment.options.set_global_str(nonconst_key, std::to_string(value));
    } else if (data.type() == "array") {
        Process::environment.options.set_global_array_int(nonconst_key, value, nullptr);
    } else {
        Process::environment.options.set_global_int(nonconst_key, value);
    }
    
#ifdef USING_BrianQC
    if (nonconst_key == "BRIANQC_ENABLE") {
        handleBrianOption(value);
    }
#endif
    
    return true;
}

bool py_psi_set_global_option_double(std::string const& key, double value) {
    std::string nonconst_key = to_upper(key);

    throw_deprecation_errors(key);

    Process::environment.options.set_global_double(nonconst_key, value);
    return true;
}

//bool py_psi_set_local_option_array(std::string const& module, std::string const& key, const py::list& values,
//                                   DataType* entry = nullptr) {
//    std::string nonconst_key = to_upper(key);
//
//    throw_deprecation_errors(key, module);
//
//    // Assign a new head entry on the first time around only
//    if (entry == nullptr) {
//        // We just do a cheesy "get" to make sure keyword is valid.  This get will throw if not.
//        std::string module_temp = Process::environment.options.get_current_module();
//        Process::environment.options.set_current_module(module);
//        Data& data = Process::environment.options[nonconst_key];
//        Process::environment.options.set_current_module(module_temp);
//        // This "if" statement is really just here to make sure the compiler doesn't optimize out the get, above.
//        if (data.type() == "array") Process::environment.options.set_array(module, nonconst_key);
//    }
//    size_t size = len(values);
//    for (int n = 0; n < size; ++n) {
//        if (py::isinstance<py::list>(values[n])) {
//            py::list l = values[n].cast<py::list>();
//            DataType* newentry = Process::environment.options.set_local_array_array(module, nonconst_key, entry);
//            // Now we need to recurse, to fill in the data
//            py_psi_set_local_option_array(module, key, l, newentry);
//        } else {
//            // This is not a list; try to cast to a string
//            try {
//                std::string s = values[n].cast<std::string>();
//                Process::environment.options.set_local_array_string(module, nonconst_key, s, entry);
//            } catch (const py::cast_error& e) {
//                try {
//                    // This is not a list or string; try to cast to an integer
//                    int i = values[n].cast<int>();
//                    Process::environment.options.set_local_array_int(module, nonconst_key, i, entry);
//                } catch (const py::cast_error& e) {
//                    // This had better be castable to a float.  We don't catch the exception here
//                    // because if we encounter one, something bad has happened
//                    double f = values[n].cast<double>();
//                    Process::environment.options.set_local_array_double(module, nonconst_key, f, entry);
//                }
//            }
//        }
//    }
//    return true;
//}
//
//bool py_psi_set_local_option_array_wrapper(std::string const& module, std::string const& key, py::list values) {
//    // A wrapper to help pybind11 handle default values
//    return py_psi_set_local_option_array(module, key, values);
//}
//
//bool py_psi_set_global_option_array(std::string const& key, py::list values, DataType* entry = nullptr) {
//    std::string nonconst_key = to_upper(key);
//
//    throw_deprecation_errors(key);
//
//    // Assign a new head entry on the first time around only
//    if (entry == nullptr) {
//        // We just do a cheesy "get" to make sure keyword is valid.  This get will throw if not.
//        Data& data = Process::environment.options[nonconst_key];
//        // This "if" statement is really just here to make sure the compiler doesn't optimize out the get, above.
//        if (data.type() == "array") Process::environment.options.set_global_array(nonconst_key);
//    }
//    size_t size = len(values);
//    for (int n = 0; n < size; ++n) {
//        if (py::isinstance<py::list>(values[n])) {
//            py::list l = values[n].cast<py::list>();
//            DataType* newentry = Process::environment.options.set_global_array_array(nonconst_key, entry);
//            // Now we need to recurse, to fill in the data
//            py_psi_set_global_option_array(key, l, newentry);
//        } else {
//            // This is not a list; try to cast to a string
//            try {
//                std::string s = values[n].cast<std::string>();
//                Process::environment.options.set_global_array_string(nonconst_key, s, entry);
//            } catch (const py::cast_error& e) {
//                try {
//                    // This is not a list or string; try to cast to an integer
//                    int i = values[n].cast<int>();
//                    Process::environment.options.set_global_array_int(nonconst_key, i, entry);
//                } catch (const py::cast_error& e) {
//                    // This had better be castable to a float.  We don't catch the exception here
//                    // because if we encounter one, something bad has happened
//                    double f = values[n].cast<double>();
//                    Process::environment.options.set_global_array_double(nonconst_key, f, entry);
//                }
//            }
//        }
//    }
//    return true;
//}
//
//bool py_psi_set_global_option_array_wrapper(std::string const& key, py::list values) {
//    // A wrapper to help pybind11 handle default values
//    return py_psi_set_global_option_array(key, values);
//}
//
//void py_psi_set_local_option_python(const std::string& key, py::object& obj) {
//    std::string nonconst_key = to_upper(key);
//    Data& data = Process::environment.options[nonconst_key];
//
//    if (data.type() == "python")
//        dynamic_cast<PythonDataType*>(data.get())->assign(obj);
//    else
//        throw PSIEXCEPTION("Unable to set option to a Python object.");
//}

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

//py::dict py_psi_options_to_python(std::string const& module) {
//    Process::environment.options.set_current_module(module);
//    py_psi_prepare_options_for_module(module);
//    std::vector<std::string> all_options = Process::environment.options.list_globals();
//
//    auto mopt = py::dict();
//    for (size_t i = 0; i < all_options.size(); i++) {
//        std::string nonconst_key = all_options[i];
//        bool in_module = Process::environment.options.exists_in_active(nonconst_key);
//        if (in_module) {
//            Data& ldata = Process::environment.options.get_local(nonconst_key);
//            bool lhoc = ldata.has_changed();
//            Data& odata = Process::environment.options.use_local(nonconst_key);
//            bool ohoc = odata.has_changed();
//            mopt[py::str(nonconst_key)] = py::make_tuple(lhoc, ohoc);
//        }
//    }
//    return mopt;
//}

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

//// Quick function that unpacks a data type
//py::list data_to_list(py::list l, Data d) {
//    if (d.is_array()) {
//        // Recurse
//        py::list row;
//        for (int i = 0; i < d.size(); ++i) {
//            data_to_list(row, d[i]);
//        }
//        l.append(row);
//    } else if (d.type() == "double") {
//        l.append(py::float_(d.to_double()));
//    } else if (d.type() == "string") {
//        l.append(py::str(d.to_string()));
//    } else if (d.type() == "boolean") {
//        l.append(py::bool_(d.to_integer()));
//    } else if (d.type() == "int") {
//        l.append(py::int_(d.to_integer()));
//    } else {
//        throw PSIEXCEPTION("Unknown data type in fill_list");
//    }
//    return l;
//}
//
//py::object py_psi_get_local_option(std::string const& module, std::string const& key) {
//    std::string nonconst_key = to_upper(key);
//    Process::environment.options.set_current_module(module);
//    py_psi_prepare_options_for_module(module);
//    Data& data = Process::environment.options.get_local(nonconst_key);
//
//    if (data.type() == "string" || data.type() == "istring")
//        return py::cast(data.to_string());
//    else if (data.type() == "boolean" || data.type() == "int")
//        return py::cast(data.to_integer());
//    else if (data.type() == "double")
//        return py::cast(data.to_double());
//    else if (data.type() == "array") {
//        py::list l;
//        for (size_t i = 0; i < data.size(); i++) {
//            data_to_list(l, data[i]);
//        }
//        return l;
//    }
//
//    return py::object();
//}
//
//py::object py_psi_get_global_option(std::string const& key) {
//    std::string nonconst_key = to_upper(key);
//    Data& data = Process::environment.options.get_global(nonconst_key);
//
//    if (data.type() == "string" || data.type() == "istring")
//        return py::cast(data.to_string());
//    else if (data.type() == "boolean" || data.type() == "int")
//        return py::cast(data.to_integer());
//    else if (data.type() == "double")
//        return py::cast(data.to_double());
//    else if (data.type() == "array") {
//        py::list l;
//        for (size_t i = 0; i < data.size(); i++) {
//            data_to_list(l, data[i]);
//        }
//        return l;
//    }
//
//    return py::object();
//}
//
//py::object py_psi_get_option(std::string const& module, std::string const& key) {
//    std::string nonconst_key = to_upper(key);
//    Process::environment.options.set_current_module(module);
//    py_psi_prepare_options_for_module(module);
//    Data& data = Process::environment.options.use_local(nonconst_key);
//
//    if (data.type() == "string" || data.type() == "istring")
//        return py::cast(data.to_string());
//    else if (data.type() == "boolean" || data.type() == "int")
//        return py::cast(data.to_integer());
//    else if (data.type() == "double")
//        return py::cast(data.to_double());
//    else if (data.type() == "array") {
//        py::list l;
//        for (size_t i = 0; i < data.size(); i++) {
//            data_to_list(l, data[i]);
//        }
//        return l;
//    }
//
//    return py::object();
//}

void py_psi_set_active_molecule(std::shared_ptr<Molecule> molecule) { Process::environment.set_molecule(molecule); }
void py_psi_set_legacy_molecule(std::shared_ptr<Molecule> legacy_molecule) {
    Process::environment.set_legacy_molecule(legacy_molecule);
}

std::shared_ptr<Molecule> py_psi_get_active_molecule() { return Process::environment.molecule(); }
std::shared_ptr<Molecule> py_psi_get_legacy_molecule() { return Process::environment.legacy_molecule(); }

void py_psi_set_gradient(SharedMatrix grad) { Process::environment.set_gradient(grad); }

SharedMatrix py_psi_get_gradient() { return Process::environment.gradient(); }

std::shared_ptr<Vector> py_psi_get_atomic_point_charges() {
    auto empty = std::make_shared<psi::Vector>();
    return empty;  // charges not added to process.h for environment - yet(?)
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
        outfile->Printf("  Threads set to %zu by Python driver.\n", nthread);
    }
#else
    Process::environment.set_n_threads(1);
    if (!quiet) {
        outfile->Printf(
            "  Python driver attempted to set threads to %zu.\n"
            "  Psi4 was compiled without OpenMP, setting threads to 1.\n",
            nthread);
    }
#endif
}

int py_psi_get_n_threads() { return Process::environment.get_n_threads(); }

std::shared_ptr<Wavefunction> py_psi_legacy_wavefunction() { return Process::environment.legacy_wavefunction(); }
void py_psi_set_legacy_wavefunction(SharedWavefunction wfn) { Process::environment.set_legacy_wavefunction(wfn); }

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

std::string py_psi_top_srcdir() { return TOSTRING(PSI_TOP_SRCDIR); }

bool psi4_python_module_initialize() {
    static bool initialized = false;

    if (initialized) {
        printf("Psi4 already initialized.\n");
        return true;
    }

    // There should only be one of these in Psi4
    Wavefunction::initialize_singletons();

    outfile = std::make_shared<PsiOutStream>();
    outfile_name = "stdout";
    std::string fprefix = PSI_DEFAULT_FILE_PREFIX;
    psi_file_prefix = strdup(fprefix.c_str());

    // There is only one timer:
    timer_init();

    // Initialize the I/O library
    // Must be done before initializing Process::environment as that needs
    // to access some globals from psio
    psio_init();

    // Setup the environment
    Process::environment.initialize();  // Defaults to obtaining the environment from the global environ variable
    Process::environment.set_memory(524288000);

    // Setup globals options
    Process::environment.options.set_read_globals(true);
    read_options("", Process::environment.options, true);
    Process::environment.options.set_read_globals(false);

#ifdef INTEL_Fortran_ENABLED
    static int argc = 1;
    static char* argv = (char*)"";
    for_rtl_init_(&argc, &argv);
#endif
    
#ifdef USING_BrianQC
    const char* brianEnableEnv = getenv("BRIANQC_ENABLE");
    brianEnableEnvFound = (bool)brianEnableEnv;
    if (brianEnableEnvFound) {
        outfile->Printf("BRIANQC_ENABLE environment variable found, checking value\n");
        brianEnableEnvValue = (bool)atoi(brianEnableEnv);
        brianEnable = brianEnableEnvValue;
        if (brianEnableEnvValue) {
            outfile->Printf("BRIANQC_ENABLE is true, attempting to initialize BrianQC\n");
            brianInit();
        }
    }
    
    const char* brianEnableDFTEnv = getenv("BRIANQC_ENABLE_DFT");
    brianEnableDFT = brianEnableDFTEnv ? (bool)atoi(brianEnableDFTEnv) : true;
#endif

    initialized = true;

    return true;
}

void psi4_python_module_finalize() {
#ifdef USING_BrianQC
    if (brianCookie != 0) {
        brianRelease();
    }
#endif
    
#ifdef INTEL_Fortran_ENABLED
    for_rtl_finish_();
#endif

    py_psi_plugin_close_all();

    // Shut things down:
    // There is only one timer:
    timer_done();

    outfile = std::shared_ptr<PsiOutStream>();
    psi_file_prefix = nullptr;
}

