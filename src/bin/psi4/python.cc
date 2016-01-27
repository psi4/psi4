/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

// This tells us the Python version number
#include <boost/python.hpp>
#include <boost/python/detail/wrap_python.hpp>
#include <boost/python/module.hpp>

#include <boost/tokenizer.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <cstdio>
#include <sstream>
#include <map>
#include <iomanip>

#include <libefp_solver/efp_solver.h>
#include <libmints/mints.h>
#include <libmints/matrix.h>
#include <libplugin/plugin.h>
#include "libparallel/mpi_wrapper.h"
#include "libparallel/local.h"
#include <liboptions/liboptions.h>
#include <liboptions/liboptions_python.h>
#include <libpsi4util/libpsi4util.h>
#include <psiconfig.h>

#include <psi4-dec.h>
#include "script.h"
#include "psi4.h"
#include "gitversion.h"
#include "libparallel/ParallelPrinter.h"
#include "../ccenergy/ccwave.h"
#include "../cclambda/cclambda.h"
//#include "../../lib/libbabel/LibBabel.h"
//#include "../mp2/mp2wave.h"

#if defined(MAKE_PYTHON_MODULE)
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include <libmints/wavefunction.h>
#include <psifiles.h>
#include "libparallel2/ParallelEnvironment.h"
namespace psi {
    int psi_start(int argc, char *argv[]);
    int psi_stop(FILE* infile, std::string, char* psi_file_prefix);
}
#endif

using namespace psi;
using namespace boost;
using namespace boost::python;
using namespace std;

// Python helper wrappers
void export_benchmarks();
void export_blas_lapack();
void export_plugins();
void export_psio();
void export_chkpt();
void export_mints();
void export_functional();
void export_oeprop();
void export_cubefile();
void export_libparallel();
void export_efp();
void export_cubeprop();

// In export_plugins.cc
void py_psi_plugin_close_all();

extern std::map<std::string, plugin_info> plugins;

#define PY_TRY(ptr, command)  \
     if(!(ptr = command)){    \
         PyErr_Print();       \
         exit(1);             \
     }

namespace opt {
psi::PsiReturnType optking(psi::Options&);
void opt_clean(void);
}

// Forward declare /src/bin/ methods
namespace psi {

// Wavefunction returns
namespace adc { SharedWavefunction     adc(SharedWavefunction, Options&); }
namespace dcft { SharedWavefunction   dcft(SharedWavefunction, Options&); }
namespace detci { SharedWavefunction detci(SharedWavefunction, Options&); }
namespace dfmp2 { SharedWavefunction dfmp2(SharedWavefunction, Options&); }
namespace dfoccwave { SharedWavefunction dfoccwave(SharedWavefunction, Options&); }
namespace libfock { SharedWavefunction libfock(SharedWavefunction, Options&); }
namespace fnocc { SharedWavefunction fnocc(SharedWavefunction, Options&); }
namespace occwave { SharedWavefunction occwave(SharedWavefunction, Options&); }
namespace mcscf { SharedWavefunction mcscf(SharedWavefunction, Options&); }
namespace scf { SharedWavefunction     scf(SharedWavefunction, Options&, PyObject *pre, PyObject *post); }

// Matrix returns
namespace deriv   { SharedMatrix     deriv(SharedWavefunction, Options&); }
namespace scfgrad { SharedMatrix   scfgrad(SharedWavefunction, Options&); }
namespace scfgrad { SharedMatrix   scfhess(SharedWavefunction, Options&); }

// Does not create a wavefunction
namespace fisapt { PsiReturnType fisapt(SharedWavefunction, Options&); }
namespace psimrcc { PsiReturnType psimrcc(SharedWavefunction, Options&); }
namespace sapt { PsiReturnType sapt(SharedWavefunction, SharedWavefunction, SharedWavefunction, Options&); }
namespace thermo { PsiReturnType thermo(SharedWavefunction, SharedVector, Options&); }

#ifdef ENABLE_CHEMPS2
namespace dmrg       { PsiReturnType dmrg(SharedWavefunction, Options&);     }
#endif

namespace mrcc {
PsiReturnType mrcc_generate_input(SharedWavefunction, Options&, const boost::python::dict&);
PsiReturnType mrcc_load_ccdensities(SharedWavefunction, Options&, const boost::python::dict&);
}

// Should die soon
namespace transqt2 { PsiReturnType transqt2(SharedWavefunction, Options&); }

// Finite difference functions
namespace findif {
std::vector<SharedMatrix> fd_geoms_1_0(boost::shared_ptr<Molecule>, Options&);
std::vector<SharedMatrix> fd_geoms_freq_0(boost::shared_ptr<Molecule>, Options&,
                                          int irrep = -1);
std::vector<SharedMatrix> fd_geoms_freq_1(boost::shared_ptr<Molecule>, Options&,
                                          int irrep = -1);
std::vector<SharedMatrix> atomic_displacements(boost::shared_ptr<Molecule>, Options&);

SharedMatrix fd_1_0(boost::shared_ptr<Molecule>, Options&, const boost::python::list&);
SharedMatrix fd_freq_0(boost::shared_ptr<Molecule>, Options&, const boost::python::list&, int irrep = -1);
SharedMatrix fd_freq_1(boost::shared_ptr<Molecule>, Options&, const boost::python::list&, int irrep = -1);
SharedMatrix displace_atom(SharedMatrix geom, const int atom,
                           const int coord, const int sign,
                           const double disp_size);
}

// CC functions
namespace ccsort { PsiReturnType ccsort(Options&); }
namespace cctransort { PsiReturnType cctransort(SharedWavefunction, Options&); }
namespace cctriples { PsiReturnType cctriples(SharedWavefunction, Options&); }
namespace cchbar { PsiReturnType cchbar(SharedWavefunction, Options&); }
namespace cclambda { PsiReturnType cclambda(SharedWavefunction, Options&); }
namespace ccdensity { PsiReturnType ccdensity(SharedWavefunction, Options&); }
namespace ccresponse {
PsiReturnType ccresponse(SharedWavefunction, Options&);
void scatter(Options&, double step, std::vector<SharedMatrix> dip, std::vector<SharedMatrix> rot,
             std::vector<SharedMatrix> quad);
}
namespace cceom { PsiReturnType cceom(SharedWavefunction, Options&); }

// No idea what to do with these yet
namespace efp { PsiReturnType efp_init(Options&); }
namespace efp { PsiReturnType efp_set_options(); }


extern int read_options(const std::string& name, Options& options, bool suppress_printing = false);
extern void print_version(std::string);

}

void py_flush_outfile()
{

}

void py_close_outfile()
{
    if (outfile) {
        outfile = boost::shared_ptr<OutFile>();
    }
}

void py_reopen_outfile()
{
    if (outfile_name == "stdout") {
        //outfile = stdout;
    }
    else {
        outfile = boost::shared_ptr<OutFile>(new OutFile(outfile_name, APPEND));
        if (!outfile)
            throw PSIEXCEPTION("PSI4: Unable to reopen output file.");
    }
}

void py_be_quiet()
{
    py_close_outfile();
    outfile = boost::shared_ptr<OutFile>(new OutFile("/dev/null", APPEND));
    if (!outfile)
        throw PSIEXCEPTION("PSI4: Unable to redirect output to /dev/null.");
}

std::string py_get_outfile_name()
{
    return outfile_name;
}

void py_psi_prepare_options_for_module(std::string const& name)
{
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

int py_psi_optking()
{
    py_psi_prepare_options_for_module("OPTKING");
    return opt::optking(Process::environment.options);
}

void py_psi_opt_clean(void)
{
    opt::opt_clean();
}

// int py_psi_mints()
// {
//     py_psi_prepare_options_for_module("MINTS");
//     return mints::mints(Process::environment.options);
// }

SharedMatrix py_psi_scfgrad(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("SCF");
    return scfgrad::scfgrad(ref_wfn, Process::environment.options);
}

SharedMatrix py_psi_scfhess(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("SCF");
    return scfgrad::scfhess(ref_wfn, Process::environment.options);
}

SharedMatrix py_psi_deriv(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("DERIV");
    return deriv::deriv(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_occ(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("OCC");
    return occwave::occwave(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_dfocc(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("DFOCC");
    return dfoccwave::dfoccwave(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_libfock(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("CPHF");
    return libfock::libfock(ref_wfn, Process::environment.options);
}

// double py_psi_scf_callbacks(PyObject *precallback, PyObject *postcallback)
// {
//     py_psi_prepare_options_for_module("SCF");
//     if (scf::scf(Process::environment.options, precallback, postcallback) == Success) {
//         return Process::environment.globals["CURRENT ENERGY"];
//     }
//     else
//         return 0.0;
// }

SharedWavefunction py_psi_mcscf(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("MCSCF");
    return mcscf::mcscf(ref_wfn, Process::environment.options);
}

PsiReturnType py_psi_mrcc_generate_input(SharedWavefunction ref_wfn, const boost::python::dict& level)
{
    py_psi_prepare_options_for_module("MRCC");
    return mrcc::mrcc_generate_input(ref_wfn, Process::environment.options, level);
}

PsiReturnType py_psi_mrcc_load_densities(SharedWavefunction ref_wfn, const boost::python::dict& level)
{
    py_psi_prepare_options_for_module("MRCC");
    return mrcc::mrcc_load_ccdensities(ref_wfn, Process::environment.options, level);
}

std::vector<SharedMatrix> py_psi_fd_geoms_1_0(boost::shared_ptr<Molecule> mol)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_geoms_1_0(mol, Process::environment.options);
}

std::vector<SharedMatrix> py_psi_fd_geoms_freq_0(boost::shared_ptr<Molecule> mol, int irrep)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_geoms_freq_0(mol, Process::environment.options, irrep);
}

std::vector<SharedMatrix> py_psi_fd_geoms_freq_1(boost::shared_ptr<Molecule> mol, int irrep)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_geoms_freq_1(mol, Process::environment.options, irrep);
}

std::vector<SharedMatrix> py_psi_atomic_displacements(boost::shared_ptr<Molecule> mol)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::atomic_displacements(mol, Process::environment.options);
}

SharedMatrix py_psi_fd_1_0(boost::shared_ptr<Molecule> mol, const boost::python::list& energies)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_1_0(mol, Process::environment.options, energies);
}

SharedMatrix py_psi_fd_freq_0(boost::shared_ptr<Molecule> mol, const boost::python::list& energies, int irrep)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_freq_0(mol, Process::environment.options, energies, irrep);
}

SharedMatrix py_psi_fd_freq_1(boost::shared_ptr<Molecule> mol, const boost::python::list& grads, int irrep)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_freq_1(mol, Process::environment.options, grads, irrep);
}

SharedMatrix py_psi_displace_atom(SharedMatrix geom, const int atom,
                                  const int coord, const int sign,
                                  const double disp_size)
{
    return findif::displace_atom(geom, atom, coord, sign, disp_size);
}

SharedWavefunction py_psi_scf(SharedWavefunction ref_wfn, PyObject *precallback,
                              PyObject *postcallback)
{
    py_psi_prepare_options_for_module("SCF");
    return scf::scf(ref_wfn, Process::environment.options, precallback, postcallback);
}

// double py_psi_scf_dummy()
// {
//     py_psi_prepare_options_for_module("SCF");
//     return scf::scf_dummy(Process::environment.options);
// }

SharedWavefunction py_psi_dcft(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("DCFT");
    return dcft::dcft(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_dfmp2(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("DFMP2");
    return dfmp2::dfmp2(ref_wfn, Process::environment.options);
}

// double py_psi_dfmp2grad()
// {
//     py_psi_prepare_options_for_module("DFMP2");
//     if (dfmp2::dfmp2grad(Process::environment.options) == Success) {
//         return Process::environment.globals["CURRENT ENERGY"];
//     }
//     else
//         return 0.0;
// }

double py_psi_sapt(SharedWavefunction Dimer, SharedWavefunction MonomerA,
                   SharedWavefunction MonomerB)
{
    py_psi_prepare_options_for_module("SAPT");
    if (sapt::sapt(Dimer, MonomerA, MonomerB, Process::environment.options) == Success) {
        return Process::environment.globals["SAPT ENERGY"];
    }
    else
        return 0.0;
}

double py_psi_fisapt(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("FISAPT");
    if (fisapt::fisapt(ref_wfn, Process::environment.options) == Success) {
        return Process::environment.globals["SAPT ENERGY"];
    }
    else
        return 0.0;
}

double py_psi_transqt2(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("TRANSQT2");
    transqt2::transqt2(ref_wfn, Process::environment.options);
    return 0.0;
}

double py_psi_ccsort()
{
    py_psi_prepare_options_for_module("CCSORT");
    ccsort::ccsort(Process::environment.options);
    return 0.0;
}

void py_psi_cctransort(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("CCTRANSORT");
    cctransort::cctransort(ref_wfn, Process::environment.options);
    // return 0.0;
}
SharedWavefunction py_psi_ccenergy(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("CCENERGY");
    SharedWavefunction ccwave(new ccenergy::CCEnergyWavefunction(
            ref_wfn,
            Process::environment.options)
    );

    std::string name = Process::environment.wavefunction()->name();
    std::string wfn_name = Process::environment.options.get_str("WFN");
    if (wfn_name != name) {
        ccwave->set_name(wfn_name);
        Process::environment.set_wavefunction(ccwave);
    }

    double energy = ccwave->compute_energy();
    return ccwave;

}

double py_psi_cctriples(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("CCTRIPLES");
    if (cctriples::cctriples(ref_wfn, Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

boost::shared_ptr<psi::efp::EFP> py_psi_efp_init()
{
    py_psi_prepare_options_for_module("EFP");
    if (psi::efp::efp_init(Process::environment.options) == Success) {
        return Process::environment.get_efp();
    }
    else
        throw PSIEXCEPTION("Unable to initialize EFP library.");
}

void py_psi_efp_set_options()
{
    py_psi_prepare_options_for_module("EFP");
    Process::environment.get_efp()->set_options();
}

SharedWavefunction py_psi_fnocc(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("FNOCC");
    return fnocc::fnocc(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_detci(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("DETCI");
    return detci::detci(ref_wfn, Process::environment.options);
}

#ifdef ENABLE_CHEMPS2
double py_psi_dmrg(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("DMRG");
    if (dmrg::dmrg(ref_wfn, Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}
#else
double py_psi_dmrg(SharedWavefunction ref_wfn)
{
    throw PSIEXCEPTION("DMRG not enabled.");
}
#endif

void py_psi_cchbar(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("CCHBAR");
    cchbar::cchbar(ref_wfn, Process::environment.options);
}

SharedWavefunction py_psi_cclambda(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("CCLAMBDA");
    boost::shared_ptr<Wavefunction> cclambda(new cclambda::CCLambdaWavefunction(
            ref_wfn,
            Process::environment.options)
    );
    Process::environment.set_wavefunction(cclambda);

    double energy = cclambda->compute_energy();
    return cclambda;
}


double py_psi_ccdensity(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("CCDENSITY");
    ccdensity::ccdensity(ref_wfn, Process::environment.options);
    return 0.0;
}

double py_psi_ccresponse(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("CCRESPONSE");
    ccresponse::ccresponse(ref_wfn, Process::environment.options);
    return 0.0;
}

void py_psi_print_list(python::list py_list)
{
    return;
}

void py_psi_scatter(double step, python::list dip_polar_list, python::list opt_rot_list,
                    python::list dip_quad_polar_list)
{
    py_psi_prepare_options_for_module("CCRESPONSE");

    // Convert python tensor lists into vectors of sharedmatrices
    std::vector<SharedMatrix> dip_polar_tensors;
    std::vector<SharedMatrix> opt_rot_tensors;
    std::vector<SharedMatrix> dip_quad_polar_tensors;

    int list_len = len(dip_polar_list);
    for (int i = 0; i < list_len; ++i) {
        python::list dip_list = extract<python::list>(dip_polar_list[i]);
        python::list rot_list = extract<python::list>(opt_rot_list[i]);
        python::list quad_list = extract<python::list>(dip_quad_polar_list[i]);
        SharedMatrix dip_mat(new Matrix(3, 3));
        SharedMatrix rot_mat(new Matrix(3, 3));
        SharedMatrix quad_mat(new Matrix(9, 3));
        for (int row = 0, j = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col, ++j) {
                dip_mat->set(row, col, extract<double>(dip_list[j]));
                rot_mat->set(row, col, extract<double>(rot_list[j]));
            }
        }
        for (int row = 0, j = 0; row < 9; ++row) {
            for (int col = 0; col < 3; ++col, ++j) {
                quad_mat->set(row, col, extract<double>(quad_list[j]));
            }
        }
        dip_polar_tensors.push_back(dip_mat);
        opt_rot_tensors.push_back(rot_mat);
        dip_quad_polar_tensors.push_back(quad_mat);
    }

//    for(std::vector<SharedMatrix>::iterator i=dip_polar_tensors.begin(); i != dip_polar_tensors.end(); ++i)
//        (*i)->print(stdout);
//    for(std::vector<SharedMatrix>::iterator i=opt_rot_tensors.begin(); i != opt_rot_tensors.end(); ++i)
//        (*i)->print(stdout);
//    for(std::vector<SharedMatrix>::iterator i=dip_quad_polar_tensors.begin(); i != dip_quad_polar_tensors.end(); ++i)
//        (*i)->print(stdout);

    ccresponse::scatter(Process::environment.options, step, dip_polar_tensors, opt_rot_tensors, dip_quad_polar_tensors);
}

double py_psi_cceom(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("CCEOM");
    if (cceom::cceom(ref_wfn, Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

double py_psi_psimrcc(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("PSIMRCC");
    psimrcc::psimrcc(ref_wfn, Process::environment.options);
    return 0.0;
}

SharedWavefunction py_psi_adc(SharedWavefunction ref_wfn)
{
    py_psi_prepare_options_for_module("ADC");
    SharedWavefunction adc_wfn = adc::adc(ref_wfn, Process::environment.options);
    return adc_wfn;
}

double py_psi_thermo(SharedWavefunction ref_wfn, SharedVector vib_freqs)
{
    py_psi_prepare_options_for_module("THERMO");
    thermo::thermo(ref_wfn, vib_freqs, Process::environment.options);
    return 0.0;
}

char const *py_psi_version()
{
#ifdef PSI_VERSION
    return PSI_VERSION;
#else
    return "";
#endif
}

char const *py_psi_git_version()
{
#ifdef GIT_VERSION
    return GIT_VERSION;
#else
    return "";
#endif
}

void py_psi_clean()
{
    PSIOManager::shared_object()->psiclean();
}

void py_psi_print_options()
{
    Process::environment.options.print();
}

void py_psi_print_global_options()
{
    Process::environment.options.print_globals();
}

boost::python::list py_psi_get_global_option_list()
{
    std::vector<std::string> options_list = Process::environment.options.list_globals();

    boost::python::list options_list_py;
    BOOST_FOREACH(const string& s, options_list) options_list_py.append(s);

    return options_list_py;
}

void py_psi_print_out(std::string s)
{
    (*outfile) << s;
}

/**
 * @return whether key describes a convergence threshold or not
 */
bool specifies_convergence(std::string const& key)
{
    return ((key.find("CONV") != key.npos) || (key.find("TOL") != key.npos));
}

bool py_psi_set_local_option_string(std::string const& module, std::string const& key, std::string const& value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "string") {
        Process::environment.options.set_str(module, nonconst_key, value);
    } else if (data.type() == "istring") {
        Process::environment.options.set_str_i(module, nonconst_key, value);
    } else if (data.type() == "boolean") {
        if (boost::to_upper_copy(value) == "TRUE" || boost::to_upper_copy(value) == "YES" || \
          boost::to_upper_copy(value) == "ON")
            Process::environment.options.set_bool(module, nonconst_key, true);
        else if (boost::to_upper_copy(value) == "FALSE" || boost::to_upper_copy(value) == "NO" || \
          boost::to_upper_copy(value) == "OFF")
            Process::environment.options.set_bool(module, nonconst_key, false);
        else
            throw std::domain_error("Required option type is boolean, no boolean specified");
    }
    return true;
}

bool py_psi_set_local_option_int(std::string const& module, std::string const& key, int value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "double" && specifies_convergence(nonconst_key)) {
        double val = pow(10.0, -value);
        Process::environment.options.set_double(module, nonconst_key, val);
    } else if (data.type() == "boolean") {
        Process::environment.options.set_bool(module, nonconst_key, value ? true : false);
    } else if (data.type() == "string" || data.type() == "istring") {
        Process::environment.options.set_str(module, nonconst_key, boost::lexical_cast<std::string>(value));
    } else {
        Process::environment.options.set_int(module, nonconst_key, value);
    }
    return true;

}

bool py_psi_set_local_option_double(std::string const& module, std::string const& key, double value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Process::environment.options.set_double(module, nonconst_key, value);
    return true;
}

template<class T>
bool is_int(T x)
{
    using boost::python::type_id;
    return type_id<T>() == type_id<int>();
}
template<class T>
bool is_double(T x)
{
    using boost::python::type_id;
    return type_id<T>() == type_id<double>();
}

bool py_psi_set_global_option_string(std::string const& key, std::string const& value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "string" || data.type() == "istring") {
        Process::environment.options.set_global_str(nonconst_key, value);
    } else if (data.type() == "boolean") {
        if (boost::to_upper_copy(value) == "TRUE" || boost::to_upper_copy(value) == "YES" || \
          boost::to_upper_copy(value) == "ON")
            Process::environment.options.set_global_bool(nonconst_key, true);
        else if (boost::to_upper_copy(value) == "FALSE" || boost::to_upper_copy(value) == "NO" || \
          boost::to_upper_copy(value) == "OFF")
            Process::environment.options.set_global_bool(nonconst_key, false);
        else
            throw std::domain_error("Required option type is boolean, no boolean specified");
    }
    return true;
}

bool py_psi_set_global_option_int(std::string const& key, int value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "double" && specifies_convergence(nonconst_key)) {
        double val = pow(10.0, -value);
        Process::environment.options.set_global_double(nonconst_key, val);
    } else if (data.type() == "boolean") {
        Process::environment.options.set_global_bool(nonconst_key, value ? true : false);
    } else if (data.type() == "string" || data.type() == "istring") {
        Process::environment.options.set_global_str(nonconst_key, boost::lexical_cast<std::string>(value));
    } else {
        Process::environment.options.set_global_int(nonconst_key, value);
    }
    return true;
}

bool py_psi_set_global_option_double(std::string const& key, double value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Process::environment.options.set_global_double(nonconst_key, value);
    return true;
}

bool py_psi_set_global_option_python(std::string const& key, boost::python::object& obj)
{
    string nonconst_key = boost::to_upper_copy(key);
    Process::environment.options.set_global_python(nonconst_key, obj);
    return true;
}

bool py_psi_set_local_option_array(std::string const& module, std::string const& key, const python::list& values,
                                   DataType *entry = NULL)
{
    string nonconst_key = boost::to_upper_copy(key);
    // Assign a new head entry on the first time around only
    if (entry == NULL) {
        // We just do a cheesy "get" to make sure keyword is valid.  This get will throw if not.
        Data& data = Process::environment.options[nonconst_key];
        // This "if" statement is really just here to make sure the compiler doesn't optimize out the get, above.
        if (data.type() == "array")
            Process::environment.options.set_array(module, nonconst_key);
    }
    size_t size = len(values);
    for (int n = 0; n < size; ++n) {
        extract<python::list> lval(values[n]);
        extract<std::string> sval(values[n]);
        extract<double> fval(values[n]);
        extract<int> ival(values[n]);
        if (lval.check()) {
            python::list l = extract<python::list>(values[n]);
            DataType *newentry = Process::environment.options.set_local_array_array(module, nonconst_key, entry);
            // Now we need to recurse, to fill in the data
            py_psi_set_local_option_array(module, key, l, newentry);
        } else if (sval.check()) {
            std::string s = extract<std::string>(values[n]);
            Process::environment.options.set_local_array_string(module, nonconst_key, s, entry);
        } else if (ival.check()) {
            int i = extract<int>(values[n]);
            Process::environment.options.set_local_array_int(module, nonconst_key, i, entry);
        } else if (fval.check()) {
            double f = extract<double>(values[n]);
            Process::environment.options.set_local_array_double(module, nonconst_key, f, entry);
        }
    }
    return true;
}


bool py_psi_set_global_option_array(std::string const& key, python::list values, DataType *entry = NULL)
{
    string nonconst_key = boost::to_upper_copy(key);
    // Assign a new head entry on the first time around only
    if (entry == NULL) {
        // We just do a cheesy "get" to make sure keyword is valid.  This get will throw if not.
        Data& data = Process::environment.options[nonconst_key];
        // This "if" statement is really just here to make sure the compiler doesn't optimize out the get, above.
        if (data.type() == "array")
            Process::environment.options.set_global_array(nonconst_key);
    }
    size_t size = len(values);
    for (int n = 0; n < size; ++n) {
        extract<python::list> lval(values[n]);
        extract<std::string> sval(values[n]);
        extract<double> fval(values[n]);
        extract<int> ival(values[n]);
        if (lval.check()) {
            python::list l = extract<python::list>(values[n]);
            DataType *newentry = Process::environment.options.set_global_array_array(nonconst_key, entry);
            // Now we need to recurse, to fill in the data
            py_psi_set_global_option_array(key, l, newentry);
        } else if (sval.check()) {
            std::string s = extract<std::string>(values[n]);
            Process::environment.options.set_global_array_string(nonconst_key, s, entry);
        } else if (ival.check()) {
            int i = extract<int>(values[n]);
            Process::environment.options.set_global_array_int(nonconst_key, i, entry);
        } else if (fval.check()) {
            double f = extract<double>(values[n]);
            Process::environment.options.set_global_array_double(nonconst_key, f, entry);
        }
    }
    return true;
}

void py_psi_set_local_option_python(const string& key, boost::python::object& obj)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "python")
        dynamic_cast<PythonDataType *>(data.get())->assign(obj);
    else
        throw PSIEXCEPTION("Unable to set option to a Python object.");
}

bool py_psi_has_local_option_changed(std::string const& module, std::string const& key)
{
    string nonconst_key = key;
    Process::environment.options.set_current_module(module);
    py_psi_prepare_options_for_module(module);
    Data& data = Process::environment.options.get_local(nonconst_key);

    return data.has_changed();
}

bool py_psi_has_global_option_changed(std::string const& key)
{
    string nonconst_key = key;
    Data& data = Process::environment.options.get_global(nonconst_key);

    return data.has_changed();
}

bool py_psi_has_option_changed(std::string const& module, std::string const& key)
{
    string nonconst_key = key;
    Process::environment.options.set_current_module(module);
    py_psi_prepare_options_for_module(module);
    Data& data = Process::environment.options.use_local(nonconst_key);

    return data.has_changed();
}

void py_psi_revoke_global_option_changed(std::string const& key)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options.get_global(nonconst_key);
    data.dechanged();
}

void py_psi_revoke_local_option_changed(std::string const& module, std::string const& key)
{

    string nonconst_key = boost::to_upper_copy(key);
    Process::environment.options.set_current_module(module);
    py_psi_prepare_options_for_module(module);
    Data& data = Process::environment.options.get_local(nonconst_key);
    data.dechanged();
}

object py_psi_get_local_option(std::string const& module, std::string const& key)
{
    string nonconst_key = key;
    Process::environment.options.set_current_module(module);
    py_psi_prepare_options_for_module(module);
    Data& data = Process::environment.options.get_local(nonconst_key);

    if (data.type() == "string" || data.type() == "istring")
        return str(data.to_string());
    else if (data.type() == "boolean" || data.type() == "int")
        return object(data.to_integer());
    else if (data.type() == "double")
        return object(data.to_double());
    else if (data.type() == "array")
        return object(data.to_list());

    return object();
}

object py_psi_get_global_option(std::string const& key)
{
    string nonconst_key = key;
    Data& data = Process::environment.options.get_global(nonconst_key);

    if (data.type() == "string" || data.type() == "istring")
        return str(data.to_string());
    else if (data.type() == "boolean" || data.type() == "int")
        return object(data.to_integer());
    else if (data.type() == "double")
        return object(data.to_double());
    else if (data.type() == "array")
        return object(data.to_list());

    return object();
}

object py_psi_get_option(std::string const& module, std::string const& key)
{
    string nonconst_key = key;
    Process::environment.options.set_current_module(module);
    py_psi_prepare_options_for_module(module);
    Data& data = Process::environment.options.use_local(nonconst_key);

    if (data.type() == "string" || data.type() == "istring")
        return str(data.to_string());
    else if (data.type() == "boolean" || data.type() == "int")
        return object(data.to_integer());
    else if (data.type() == "double")
        return object(data.to_double());
    else if (data.type() == "array")
        return object(data.to_list());

    return object();
}

void py_psi_set_active_molecule(boost::shared_ptr<Molecule> molecule)
{
    Process::environment.set_molecule(molecule);
}

void py_psi_set_parent_symmetry(std::string pg)
{
    boost::shared_ptr<PointGroup> group = boost::shared_ptr<PointGroup>();
    if (pg != "") {
        group = boost::shared_ptr<PointGroup>(new PointGroup(pg));
    }

    Process::environment.set_parent_symmetry(group);
}

boost::shared_ptr<Molecule> py_psi_get_active_molecule()
{
    return Process::environment.molecule();
}

boost::shared_ptr<psi::efp::EFP> py_psi_get_active_efp()
{
    return Process::environment.get_efp();
}

void py_psi_set_gradient(SharedMatrix grad)
{
    if (Process::environment.wavefunction()) {
        Process::environment.wavefunction()->set_gradient(grad);
    } else {
        Process::environment.set_gradient(grad);
    }
}

SharedMatrix py_psi_get_gradient()
{
    if (Process::environment.wavefunction()) {
        boost::shared_ptr<Wavefunction> wf = Process::environment.wavefunction();
        return wf->gradient();
    } else {
        return Process::environment.gradient();
    }
}

void py_psi_set_efp_torque(SharedMatrix torq)
{
    if (Process::environment.get_efp()->get_frag_count() > 0) {
        Process::environment.get_efp()->set_torque(torq);
    } else {
        Process::environment.set_efp_torque(torq);
    }
}

SharedMatrix py_psi_get_efp_torque()
{
    if (Process::environment.get_efp()->get_frag_count() > 0) {
        boost::shared_ptr<psi::efp::EFP> efp = Process::environment.get_efp();
        return efp->torque();
    } else {
        return Process::environment.efp_torque();
    }
}

void py_psi_set_frequencies(boost::shared_ptr<Vector> freq)
{
    if (Process::environment.wavefunction()) {
        Process::environment.wavefunction()->set_frequencies(freq);
    } else {
        Process::environment.set_frequencies(freq);
    }
}

void py_psi_set_normalmodes(boost::shared_ptr<Vector> norm)
{
    if (Process::environment.wavefunction()) {
        Process::environment.wavefunction()->set_normalmodes(norm);
    } else {
//      Process::environment.set_normalmodes(norm);
    }
}

boost::shared_ptr<Vector> py_psi_get_frequencies()
{
    if (Process::environment.wavefunction()) {
        boost::shared_ptr<Wavefunction> wf = Process::environment.wavefunction();
        return wf->frequencies();
    } else {
        return Process::environment.frequencies();
    }
}

boost::shared_ptr<Vector> py_psi_get_atomic_point_charges()
{
    if (Process::environment.wavefunction()) {
        boost::shared_ptr<Wavefunction> wf = Process::environment.wavefunction();
        return wf->get_atomic_point_charges();
    }
    else {
        boost::shared_ptr<psi::Vector> empty(new psi::Vector());
        return empty; // charges not added to process.h for environment - yet(?)
    }
}

boost::shared_ptr<Vector> py_psi_get_normalmodes()
{
    if (Process::environment.wavefunction()) {
        boost::shared_ptr<Wavefunction> wf = Process::environment.wavefunction();
        return wf->normalmodes();
    } else {
//      return Process::environment.normalmodes();
    }
}

double py_psi_get_variable(const std::string& key)
{
    string uppercase_key = key;
    transform(uppercase_key.begin(), uppercase_key.end(), uppercase_key.begin(), ::toupper);
    return Process::environment.globals[uppercase_key];
}

SharedMatrix py_psi_get_array_variable(const std::string& key)
{
    string uppercase_key = key;
    transform(uppercase_key.begin(), uppercase_key.end(), uppercase_key.begin(), ::toupper);
    return Process::environment.arrays[uppercase_key];
}

void py_psi_set_variable(const std::string& key, double val)
{
    string uppercase_key = key;
    transform(uppercase_key.begin(), uppercase_key.end(), uppercase_key.begin(), ::toupper);
    Process::environment.globals[uppercase_key] = val;
}

void py_psi_set_array_variable(const std::string& key, SharedMatrix val)
{
    string uppercase_key = key;
    transform(uppercase_key.begin(), uppercase_key.end(), uppercase_key.begin(), ::toupper);
    Process::environment.arrays[uppercase_key] = val;
}

void py_psi_clean_variable_map()
{
    Process::environment.globals.clear();
    Process::environment.arrays.clear();
}

void py_psi_set_memory(unsigned long int mem)
{
    Process::environment.set_memory(mem);
    outfile->Printf("\n  Memory set to %7.3f %s by Python script.\n", (mem > 1000000000 ? mem / 1.0E9 : mem / 1.0E6), \
        (mem > 1000000000 ? "GiB" : "MiB"));
}

unsigned long int py_psi_get_memory()
{
    return Process::environment.get_memory();
}

void py_psi_set_n_threads(int nthread)
{
    Process::environment.set_n_threads(nthread);
}

int py_psi_get_n_threads()
{
    return Process::environment.get_n_threads();
}

boost::shared_ptr<Wavefunction> py_psi_wavefunction()
{
    return Process::environment.wavefunction();
}
void py_psi_set_wavefunction(SharedWavefunction wfn)
{
    Process::environment.set_wavefunction(wfn);
}
SharedWavefunction py_psi_new_wavefunction(boost::shared_ptr<Molecule> molecule,
                                           const std::string& basis)
{
    // Ultimately options will not go here
    return SharedWavefunction(new Wavefunction(molecule, basis, Process::environment.options));
}

string py_psi_get_input_directory()
{
    return infile_directory;
}

void py_psi_print_variable_map()
{
    int largest_key = 0;
    for (std::map<string, double>::iterator it = Process::environment.globals.begin();
         it != Process::environment.globals.end(); ++it) {
        if (it->first.size() > largest_key)
            largest_key = it->first.size();
    }
    largest_key += 2;  // for quotation marks

    std::stringstream line;
    std::string first_tmp;
    for (std::map<string, double>::iterator it = Process::environment.globals.begin();
         it != Process::environment.globals.end(); ++it) {
        first_tmp = "\"" + it->first + "\"";
        line << "  " << std::left << std::setw(largest_key) << first_tmp << " => " << std::setw(20) << std::right <<
        std::fixed << std::setprecision(12) << it->second << std::endl;
    }

    outfile->Printf("\n\n  Variable Map:");
    outfile->Printf("\n  ----------------------------------------------------------------------------\n");
    outfile->Printf("%s\n\n", line.str().c_str());
}

// Converts a C++ map to a python dict
// from https://gist.github.com/octavifs/5362297
template<class K, class V>
boost::python::dict toPythonDict(std::map<K, V> map)
{
    typename std::map<K, V>::iterator iter;
    boost::python::dict dictionary;
    for (iter = map.begin(); iter != map.end(); ++iter) {
        dictionary[iter->first] = iter->second;
    }
    return dictionary;
}

boost::python::dict py_psi_return_variable_map()
{
    return toPythonDict(Process::environment.globals);
}

boost::python::dict py_psi_return_array_variable_map()
{
    return toPythonDict(Process::environment.arrays);
}

std::string py_psi_top_srcdir()
{
    return PSI_TOP_SRCDIR;
}

#if defined(MAKE_PYTHON_MODULE)
bool psi4_python_module_initialize()
{
    static bool initialized = false;

    if (initialized) {
        printf("Psi4 already initialized.\n");
        return true;
    }

    print_version("stdout");

    // Track down the location of PSI4's python script directory.
    std::string psiDataDirName = Process::environment("PSIDATADIR");
    std::string psiDataDirWithPython = psiDataDirName + "/psi4";
    boost::filesystem::path bf_path;
    bf_path = boost::filesystem::system_complete(psiDataDirWithPython);
    // printf("Python dir is at %s\n", psiDataDirName.c_str());
    if(!boost::filesystem::is_directory(bf_path)) {
        printf("Unable to read the PSI4 Python folder - check the PSIDATADIR environmental variable\n"
                "      Current value of PSIDATADIR is %s\n", psiDataDirName.c_str());
        return false;
    }

    // Add PSI library python path
    PyObject *path, *sysmod, *str;
    PY_TRY(sysmod , PyImport_ImportModule("sys"));
    PY_TRY(path   , PyObject_GetAttrString(sysmod, "path"));
#if PY_MAJOR_VERSION == 2
    PY_TRY(str    , PyString_FromString(psiDataDirWithPython.c_str()));
#else
    PY_TRY(str    , PyUnicode_FromString(psiDataDirWithPython.c_str()));
#endif
    PyList_Append(path, str);
    Py_DECREF(str);
    Py_DECREF(path);
    Py_DECREF(sysmod);

    initialized = true;

    return true;
}

void psi4_python_module_finalize()
{
    Process::environment.wavefunction().reset();
    py_psi_plugin_close_all();

    // Shut things down:
    // There is only one timer:
    timer_done();

    psi_stop(infile, "outfile", psi_file_prefix);
    Script::language->finalize();

}
#endif

void translate_psi_exception(const PsiException& e)
{
#ifdef DEBUG
#if PY_MAJOR_VERSION == 2
    PyObject *message = PyString_FromFormat("%s (%s:%d)", e.what(), e.file(), e.line());
#else
    PyObject *message = PyUnicode_FromFormat("%s (%s:%d)", e.what(), e.file(), e.line());
#endif
    PyErr_SetObject(PyExc_RuntimeError, message);
    Py_DECREF(message);
#else
    PyErr_SetString(PyExc_RuntimeError, e.what());
#endif
}

// Tell python about the default final argument to the array setting functions
BOOST_PYTHON_FUNCTION_OVERLOADS(set_global_option_overloads, py_psi_set_global_option_array, 2, 3)

BOOST_PYTHON_FUNCTION_OVERLOADS(set_local_option_overloads, py_psi_set_local_option_array, 3, 4)

BOOST_PYTHON_MODULE (psi4)
{
#if defined(MAKE_PYTHON_MODULE)
    // Initialize the world communicator
    WorldComm = boost::shared_ptr<LibParallel::ParallelEnvironment>(
          new LibParallel::ParallelEnvironment(0, 0));

    // Setup the environment
    Process::arguments.initialize(0, 0);
    Process::environment.initialize(); // Defaults to obtaining the environment from the global environ variable

    // There is only one timer:
    timer_init();

    // There should only be one of these in Psi4
    Wavefunction::initialize_singletons();

    if(psi_start(0, 0) == PSI_RETURN_FAILURE) return;

    // Initialize the I/O library
    psio_init();

    // Setup globals options
    Process::environment.options.set_read_globals(true);
    read_options("", Process::environment.options, true);
    Process::environment.options.set_read_globals(false);

    def("initialize", &psi4_python_module_initialize);
    def("finalize", &psi4_python_module_finalize);
#endif


    register_exception_translator<PsiException>(&translate_psi_exception);

//    docstring_options sphx_doc_options(true, true, false);

    enum_<PsiReturnType>("PsiReturnType", "docstring")
            .value("Success", Success)
            .value("Failure", Failure)
            .value("Balk", Balk)
            .value("EndLoop", EndLoop)
            .export_values();


    def("version", py_psi_version, "Returns the version ID of this copy of Psi.");
    def("git_version", py_psi_git_version, "Returns the git version of this copy of Psi.");
    def("clean", py_psi_clean, "Function to remove scratch files. Call between independent jobs.");

    // Benchmarks
    export_benchmarks();

    // BLAS/LAPACK Static Wrappers
    export_blas_lapack();

    // Plugins
    export_plugins();

    // OEProp/GridProp
    export_oeprop();

    // CubeFile
    export_cubefile();

    // EFP
    export_efp();

    // CubeProperties
    export_cubeprop();

    // Options
// The following line was conflct between master and roa branch (TDC, 10/29/2014)
//    def("print_list", py_psi_print_list, "Prints a python list using a C function.");
    def("prepare_options_for_module",
        py_psi_prepare_options_for_module,
        "Sets the options module up to return options pertaining to the named argument (e.g. SCF).");
    def("set_active_molecule",
        py_psi_set_active_molecule,
        "Activates a previously defined (in the input) molecule, by name.");
    def("get_active_molecule", &py_psi_get_active_molecule, "Returns the currently active molecule object.");
    def("wavefunction",
        py_psi_wavefunction,
        "Returns the current wavefunction object from the most recent computation.");
    def("set_wavefunction",
        py_psi_set_wavefunction,
        "Returns the current wavefunction object from the most recent computation.");
    def("new_wavefunction",
        py_psi_new_wavefunction,
        "Builds a new wavefunction from scratch.");
    def("get_gradient", py_psi_get_gradient, "Returns the most recently computed gradient, as a N by 3 Matrix object.");
    def("set_gradient",
        py_psi_set_gradient,
        "Assigns the global gradient to the values stored in the N by 3 Matrix argument.");
    def("get_active_efp", &py_psi_get_active_efp, "Returns the currently active EFP object.");
    def("get_efp_torque",
        py_psi_get_efp_torque,
        "Returns the most recently computed gradient for the EFP portion, as a Nefp by 6 Matrix object.");
    def("set_efp_torque",
        py_psi_set_efp_torque,
        "Assigns the global EFP gradient to the values stored in the Nefp by 6 Matrix argument.");
    def("get_frequencies",
        py_psi_get_frequencies,
        "Returns the most recently computed frequencies, as a 3N-6 Vector object.");
    def("get_atomic_point_charges",
        py_psi_get_atomic_point_charges,
        "Returns the most recently computed atomic point charges, as a double * object.");
    def("get_normalmodes", py_psi_get_normalmodes, "Returns the most recently computed normal modes Vector object.");
    def("set_frequencies",
        py_psi_set_frequencies,
        "Assigns the global frequencies to the values stored in the 3N-6 Vector argument.");
    def("set_normalmodes",
        py_psi_set_normalmodes,
        "Assigns the global normalmodes to the values stored in a Vector argument.");
    def("set_memory", py_psi_set_memory, "Sets the memory available to Psi (in bytes).");
    def("get_memory", py_psi_get_memory, "Returns the amount of memory available to Psi (in bytes).");
    def("set_nthread", &py_psi_set_n_threads, "Sets the number of threads to use in SMP parallel computations.");
    def("nthread", &py_psi_get_n_threads, "Returns the number of threads to use in SMP parallel computations.");
//    def("mol_from_file",&LibBabel::ParseFile,"Reads a molecule from another input file");
    def("set_parent_symmetry",
        py_psi_set_parent_symmetry,
        "Sets the symmetry of the 'parent' (undisplaced) geometry, by Schoenflies symbol, at the beginning of a finite difference computation.");
    def("print_options",
        py_psi_print_options,
        "Prints the currently set options (to the output file) for the current module.");
    def("print_global_options",
        py_psi_print_global_options,
        "Prints the currently set global (all modules) options to the output file.");
    def("print_out", py_psi_print_out, "Prints a string (using sprintf-like notation) to the output file.");

    // Set the different local option types
    def("set_local_option",
        py_psi_set_local_option_string,
        "Sets value *arg3* to string keyword *arg2* scoped only to a specific module *arg1*.");
    def("set_local_option",
        py_psi_set_local_option_double,
        "Sets value *arg3* to double keyword *arg2* scoped only to a specific module *arg1*.");
    def("set_local_option",
        py_psi_set_local_option_int,
        "Sets value *arg3* to integer keyword *arg2* scoped only to a specific module *arg1*.");
    def("set_local_option", py_psi_set_local_option_array, set_local_option_overloads());
    def("set_local_option_python",
        py_psi_set_local_option_python,
        "Sets an option to a Python object, but scoped only to a single module.");

    // Set the different global option types
    def("set_global_option",
        py_psi_set_global_option_string,
        "Sets value *arg2* to string keyword *arg1* for all modules.");
    def("set_global_option",
        py_psi_set_global_option_double,
        "Sets value *arg2* to double keyword *arg1* for all modules.");
    def("set_global_option",
        py_psi_set_global_option_int,
        "Sets value *arg2* to integer keyword *arg1* for all modules.");
    def("set_global_option", py_psi_set_global_option_array, set_global_option_overloads());
    def("set_global_option_python", py_psi_set_global_option_python, "Sets a global option to a Python object type.");

    // Print options list
    def("get_global_option_list", py_psi_get_global_option_list, "Returns a list of all global options.");

    // Get the option; either global or local or let liboptions decide whether to use global or local
    def("get_global_option",
        py_psi_get_global_option,
        "Given a string of a keyword name *arg1*, returns the value associated with the keyword from the global options. Returns error if keyword is not recognized.");
    def("get_local_option",
        py_psi_get_local_option,
        "Given a string of a keyword name *arg2* and a particular module *arg1*, returns the value associated with the keyword in the module options scope. Returns error if keyword is not recognized for the module.");
    def("get_option",
        py_psi_get_option,
        "Given a string of a keyword name *arg2* and a particular module *arg1*, returns the local value associated with the keyword if it's been set, else the global value if it's been set, else the local default value. Returns error if keyword is not recognized globally or if keyword is not recognized for the module.");

    // Returns whether the option has changed/revoke has changed for silent resets
    def("has_global_option_changed",
        py_psi_has_global_option_changed,
        "Returns boolean for whether the keyword *arg1* has been touched in the global scope, by either user or code. Notwithstanding, code is written such that in practice, this returns whether the option has been touched in the global scope by the user.");
    def("has_local_option_changed",
        py_psi_has_local_option_changed,
        "Returns boolean for whether the keyword *arg2* has been touched in the scope of the specified module *arg1*, by either user or code. Notwithstanding, code is written such that in practice, this returns whether the option has been touched in the module scope by the user.");
    def("has_option_changed",
        py_psi_has_option_changed,
        "Returns boolean for whether the option *arg2* has been touched either locally to the specified module *arg1* or globally, by either user or code. Notwithstanding, code is written such that in practice, this returns whether the option has been touched by the user.");
    def("revoke_global_option_changed",
        py_psi_revoke_global_option_changed,
        "Given a string of a keyword name *arg1*, sets the has_changed attribute in the global options scope to false. Used in python driver when a function sets the value of an option. Before the function exits, this command is called on the option so that has_changed reflects whether the user (not the program) has touched the option.");
    def("revoke_local_option_changed",
        py_psi_revoke_local_option_changed,
        "Given a string of a keyword name *arg2* and a particular module *arg1*, sets the has_changed attribute in the module options scope to false. Used in python driver when a function sets the value of an option. Before the function exits, this command is called on the option so that has_changed reflects whether the user (not the program) has touched the option.");

    // These return/set/print PSI variables found in Process::environment.globals
    def("get_variable",
        py_psi_get_variable,
        "Returns one of the PSI variables set internally by the modules or python driver (see manual for full listing of variables available).");
    def("set_variable", py_psi_set_variable, "Sets a PSI variable, by name.");
    def("print_variables", py_psi_print_variable_map, "Prints all PSI variables that have been set internally.");
    def("clean_variables", py_psi_clean_variable_map, "Empties all PSI variables that have set internally.");
    def("get_variables",
        py_psi_return_variable_map,
        "Returns dictionary of the PSI variables set internally by the modules or python driver.");
    def("get_array_variable",
        py_psi_get_array_variable,
        "Returns one of the PSI variables set internally by the modules or python driver (see manual for full listing of variables available).");
    def("set_array_variable", py_psi_set_array_variable, "Sets a PSI variable, by name.");
//    def("print_array_variables", py_psi_print_array_variable_map, "Prints all PSI variables that have been set internally.");
    def("get_array_variables",
        py_psi_return_array_variable_map,
        "Returns dictionary of the PSI variables set internally by the modules or python driver.");

    // Get the name of the directory where the input file is at
    def("get_input_directory", py_psi_get_input_directory, "Returns the location of the input file.");

    // Returns the location where the Psi4 source is located.
    def("psi_top_srcdir", py_psi_top_srcdir, "Returns the location of the source code.");

    def("flush_outfile", py_flush_outfile, "Flushes the output file.");
    def("close_outfile", py_close_outfile, "Closes the output file.");
    def("reopen_outfile", py_reopen_outfile, "Reopens the output file.");
    def("outfile_name", py_get_outfile_name, "Returns the name of the output file.");
    def("be_quiet",
        py_be_quiet,
        "Redirects output to /dev/null.  To switch back to regular output mode, use reopen_outfile()");

    // modules
    def("deriv",
        py_psi_deriv,
        "Runs deriv, which contracts density matrices with derivative integrals, to compute gradients.");
    def("scfgrad", py_psi_scfgrad, "Run scfgrad, which is a specialized DF-SCF gradient program.");
    def("scfhess", py_psi_scfhess, "Run scfhess, which is a specialized DF-SCF hessian program.");

    def("scf", py_psi_scf, "Runs the SCF code.");
    def("dcft", py_psi_dcft, "Runs the density cumulant functional theory code.");
    def("libfock", py_psi_libfock, "Runs a CPHF calculation, using libfock.");
    def("dfmp2", py_psi_dfmp2, "Runs the DF-MP2 code.");
    // def("dfmp2grad", py_psi_dfmp2grad, "Runs the DF-MP2 gradient.");
//    def("mp2", py_psi_mp2, "Runs the conventional (slow) MP2 code.");
    def("mcscf", py_psi_mcscf, "Runs the MCSCF code, (N.B. restricted to certain active spaces).");
    def("mrcc_generate_input", py_psi_mrcc_generate_input, "Generates an input for Kallay's MRCC code.");
    def("mrcc_load_densities", py_psi_mrcc_load_densities, "Reads in the density matrices from Kallay's MRCC code.");
    def("fd_geoms_1_0", py_psi_fd_geoms_1_0,
        "Gets list of displacements needed for a finite difference gradient computation, from energy points.");
    def("fd_geoms_freq_0", py_psi_fd_geoms_freq_0,
        "Gets list of displacements needed for a finite difference frequency computation, from energy points, for a given irrep.");
    def("fd_geoms_freq_1", py_psi_fd_geoms_freq_1,
        "Gets list of displacements needed fof a finite difference frequency computation, from gradients, for a given irrep");
    def("fd_1_0", py_psi_fd_1_0,
        "Performs a finite difference gradient computation, from energy points.");
    def("fd_freq_0", py_psi_fd_freq_0,
        "Performs a finite difference frequency computation, from energy points, for a given irrep.");
    def("fd_freq_1", py_psi_fd_freq_1,
        "Performs a finite difference frequency computation, from gradients, for a given irrep.");
    def("atomic_displacements",
        py_psi_atomic_displacements,
        "Returns list of displacements generated by displacing each atom in the +/- x, y, z directions");
    def("displace_atom", py_psi_displace_atom, "Displaces one coordinate of single atom.");
    def("sapt", py_psi_sapt, "Runs the symmetry adapted perturbation theory code.");
    def("fisapt", py_psi_fisapt, "Runs the functional-group intramolecular symmetry adapted perturbation theory code.");
    def("psimrcc", py_psi_psimrcc, "Runs the multireference coupled cluster code.");
    def("optking", py_psi_optking, "Runs the geometry optimization / frequency analysis code.");
//    def("transqt", py_psi_transqt, "Runs the (deprecated) transformation code.");
    def("transqt2", py_psi_transqt2, "Runs the (deprecated) transformation code.");
    def("ccsort", py_psi_ccsort, "Runs CCSORT, which reorders integrals for use in the coupled cluster codes.");
    def("cctransort", py_psi_cctransort, "Runs CCTRANSORT, which transforms and reorders integrals for use in the coupled cluster codes.");
    def("ccenergy", py_psi_ccenergy, "Runs the coupled cluster energy code.");
    def("cctriples", py_psi_cctriples, "Runs the coupled cluster (T) energy code.");
    def("detci", py_psi_detci, "Runs the determinant-based configuration interaction code.");
    def("dmrg", py_psi_dmrg, "Runs the DMRG code.");
// DGAS    def("detcas", py_psi_detcas, "Runs the determinant-based complete active space self consistent field.");
    def("fnocc", py_psi_fnocc, "Runs the fno-ccsd(t)/qcisd(t)/mp4/cepa energy code");
    def("efp_init", py_psi_efp_init, "Initializes the EFP library and returns an EFP object.");
    def("efp_set_options", py_psi_efp_set_options, "Set EFP options from environment options object.");
    def("cchbar", py_psi_cchbar, "Runs the code to generate the similarity transformed Hamiltonian.");
    def("cclambda", py_psi_cclambda, "Runs the coupled cluster lambda equations code.");
    def("ccdensity", py_psi_ccdensity, "Runs the code to compute coupled cluster density matrices.");
    def("ccresponse", py_psi_ccresponse, "Runs the coupled cluster response theory code.");
    def("scatter", py_psi_scatter, "New Scatter function.");
    def("cceom", py_psi_cceom, "Runs the equation of motion coupled cluster code, for excited states.");
    def("occ", py_psi_occ, "Runs the orbital optimized CC codes.");
    def("dfocc", py_psi_dfocc, "Runs the density-fitted orbital optimized CC codes.");
    def("adc", py_psi_adc, "Runs the ADC propagator code, for excited states.");
    def("thermo", py_psi_thermo, "Computes thermodynamic data.");
    def("opt_clean", py_psi_opt_clean, "Cleans up the optimizer's scratch files.");

    // Define library classes
    export_psio();
    export_chkpt();
    export_mints();
    export_functional();
    export_libparallel();

    typedef string (Process::Environment::*environmentStringFunction)(const string&);

    class_<Process::Environment>("Environment").
            def("__getitem__", environmentStringFunction(&Process::Environment::operator()), "docstring");

    class_<Process>("Process").
            add_static_property("environment", Process::get_environment, "docstring");
}

Python::Python() : Script()
{

}

Python::~Python()
{

}

void Python::initialize()
{
}

void Python::finalize()
{
//    Py_Finalize();
}

void Python::run(FILE *input)
{
    using namespace boost::python;
    char *s = 0;

    if (!Py_IsInitialized()) {
        s = strdup("psi");

#if PY_MAJOR_VERSION == 2
        if (PyImport_AppendInittab(strdup("psi4"), initpsi4) == -1) {
            outfile->Printf("Unable to register psi4 with your Python.\n");
            abort();
        }
#else
        if (PyImport_AppendInittab(strdup("psi4"), PyInit_psi4) == -1) {
            outfile->Printf( "Unable to register psi4 with your Python.\n");
            abort();
        }
#endif

        // Py_InitializeEx(0) causes sig handlers to not be installed.
        Py_InitializeEx(0);
#if PY_VERSION_HEX >= 0x03000000
        Py_SetProgramName(L"psi");
#else
        Py_SetProgramName(s);
#endif

        // Track down the location of PSI4's auxiliary directories path
        std::string psiPath = Process::environment("PSIPATH") + ":./";
        boost::char_separator<char> sep(":");
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
        tokenizer tokens(psiPath, sep);
        PyObject *path, *sysmod, *str;
        PY_TRY(sysmod, PyImport_ImportModule("sys"));
        PY_TRY(path, PyObject_GetAttrString(sysmod, "path"));
        for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
            boost::filesystem::path bf_path2;
            bf_path2 = boost::filesystem::system_complete(*tok_iter);
            if (!boost::filesystem::is_directory(bf_path2)) {
                printf("Unable to read the PSI4 Auxililary folder - check the PSIPATH environmental variable\n"
                               "      Current value of PSIPATH is %s\n", psiPath.c_str());
                exit(1);
            }
#if PY_MAJOR_VERSION == 2
            PY_TRY(str, PyString_FromString((*tok_iter).c_str()));
#else
            PY_TRY(str    , PyUnicode_FromString((*tok_iter).c_str()));
#endif
            PyList_Append(path, str);
        }
        Py_DECREF(str);
        Py_DECREF(path);
        Py_DECREF(sysmod);

        // Track down the location of PSI4's python script directory.
        std::string psiDataDirName = Process::environment("PSIDATADIR");
        std::string psiDataDirWithPython = psiDataDirName + "/python";
        boost::filesystem::path bf_path;
        bf_path = boost::filesystem::system_complete(psiDataDirWithPython);
        // printf("Python dir is at %s\n", psiDataDirName.c_str());
        if (!boost::filesystem::is_directory(bf_path)) {
            printf("Unable to read the PSI4 Python folder - check the PSIDATADIR environmental variable\n"
                           "      Current value of PSIDATADIR is %s\n", psiDataDirName.c_str());
            exit(1);
        }

        // Add PSI library python path
        PY_TRY(sysmod, PyImport_ImportModule("sys"));
        PY_TRY(path, PyObject_GetAttrString(sysmod, "path"));
#if PY_MAJOR_VERSION == 2
        PY_TRY(str, PyString_FromString(psiDataDirWithPython.c_str()));
#else
        PY_TRY(str    , PyUnicode_FromString(psiDataDirWithPython.c_str()));
#endif

        // Append to the path list
        PyList_Append(path, str);

        Py_DECREF(str);
        Py_DECREF(path);
        Py_DECREF(sysmod);
    }
    if (Py_IsInitialized()) {

        try {
            string inputfile;
            object objectMain(handle<>(borrowed(PyImport_AddModule("__main__"))));
            object objectDict = objectMain.attr("__dict__");
            s = strdup("import psi4");
            PyRun_SimpleString(s);

            if (!interactive_python) {

                // Stupid way to read in entire file.
                char line[256];
                std::stringstream file;
                while (fgets(line, sizeof(line), input)) {
                    file << line;
                }

                if (!skip_input_preprocess) {
                    // Process the input file
                    PyObject *input;
                    PY_TRY(input, PyImport_ImportModule("inputparser"));
                    PyObject *function;
                    PY_TRY(function, PyObject_GetAttrString(input, "process_input"));
                    PyObject *pargs;
                    PY_TRY(pargs, Py_BuildValue("(s)", file.str().c_str()));
                    PyObject *ret;
                    PY_TRY(ret, PyEval_CallObject(function, pargs));

                    char *val;
                    PyArg_Parse(ret, "s", &val);
                    inputfile = val;

                    Py_DECREF(ret);
                    Py_DECREF(pargs);
                    Py_DECREF(function);
                    Py_DECREF(input);
                }
                else
                    inputfile = file.str();

                if (verbose) {
                    outfile->Printf("\n Input file to run:\n%s", inputfile.c_str());

                }

                str strStartScript(inputfile);

                object objectScriptInit = exec(strStartScript, objectDict, objectDict);
            }
            else { // interactive python
                // Process the input file
                PyObject *input;

                PY_TRY(input, PyImport_ImportModule("interactive"));
                PyObject *function;
                PY_TRY(function, PyObject_GetAttrString(input, "run"));
                PyObject *ret;
                PY_TRY(ret, PyEval_CallObject(function, NULL));
            }
        }
        catch (error_already_set const& e) {
            PyErr_Print();
            exit(1);
        }
    }
    else {
        outfile->Printf("Unable to run Python input file.\n");
        return;
    }

    if (s)
        free(s);
    Process::environment.molecule().reset();
    Process::environment.wavefunction().reset();
    py_psi_plugin_close_all();
}
