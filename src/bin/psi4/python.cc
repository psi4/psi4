#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <cstdio>
#include <sstream>
#include <map>
#include <iomanip>

#include <libmints/mints.h>
#include <libplugin/plugin.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <liboptions/python.h>
#include <psiconfig.h>

#include <psi4-dec.h>
#include "script.h"
#include "psi4.h"

#include "../ccenergy/ccwave.h"
#include "../mp2/mp2wave.h"

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

// In export_plugins.cc
void py_psi_plugin_close_all();

extern std::map<std::string, plugin_info> plugins;

namespace opt {
  psi::PsiReturnType optking(psi::Options &);
  void opt_clean(void);
}

namespace psi {
    namespace mints      { PsiReturnType mints(Options &);    }
    namespace deriv      { PsiReturnType deriv(Options &);    }
    namespace scfgrad    { PsiReturnType scfgrad(Options &); }
    namespace scf        { PsiReturnType scf(Options &, PyObject* pre, PyObject* post);   }
    namespace libfock    { PsiReturnType libfock(Options &);  }
    namespace dfmp2      { PsiReturnType dfmp2(Options &);    }
    namespace dfmp2      { PsiReturnType dfmp2grad(Options &);}
    namespace dfcc       { PsiReturnType dfcc(Options &);     }
    namespace sapt       { PsiReturnType sapt(Options &);     }
    namespace dcft       { PsiReturnType dcft(Options &);     }
    namespace mcscf      { PsiReturnType mcscf(Options &);    }
    namespace psimrcc    { PsiReturnType psimrcc(Options &);  }
    namespace transqt    { PsiReturnType transqt(Options &);  }
    namespace transqt2   { PsiReturnType transqt2(Options &); }
    namespace ccsort     { PsiReturnType ccsort(Options&);    }
//    namespace lmp2       { PsiReturnType lmp2(Options&);      }
    namespace cctriples  { PsiReturnType cctriples(Options&); }
    namespace cchbar     { PsiReturnType cchbar(Options&);    }
    namespace cclambda   { PsiReturnType cclambda(Options&);  }
    namespace ccdensity  { PsiReturnType ccdensity(Options&); }
    namespace ccresponse { PsiReturnType ccresponse(Options&);}
    namespace cceom      { PsiReturnType cceom(Options&);     }
    namespace detci      { PsiReturnType detci(Options&);     }
    namespace cepa { PsiReturnType cepa(Options&);}
    namespace stable     { PsiReturnType stability(Options&); }
    namespace omp2wave   { PsiReturnType omp2wave(Options&);  }
    namespace adc        { PsiReturnType adc(Options&);       }
    namespace mrcc       {
        PsiReturnType mrcc_generate_input(Options&, const boost::python::dict&);
        PsiReturnType mrcc_load_ccdensities(Options&, const boost::python::dict&);
    }
    namespace findif    {
      std::vector< boost::shared_ptr<Matrix> > fd_geoms_1_0(Options &);
      //std::vector< boost::shared_ptr<Matrix> > fd_geoms_2_0(Options &);
      std::vector< boost::shared_ptr<Matrix> > fd_geoms_freq_0(Options &, int irrep=-1);
      std::vector< boost::shared_ptr<Matrix> > fd_geoms_freq_1(Options &, int irrep=-1);
      std::vector< boost::shared_ptr<Matrix> > fd_geoms_hessian_0(Options &);

      PsiReturnType fd_1_0(Options &, const boost::python::list&);
      //PsiReturnType fd_2_0(Options &, const boost::python::list&);
      PsiReturnType fd_freq_0(Options &, const boost::python::list&, int irrep=-1);
      PsiReturnType fd_freq_1(Options &, const boost::python::list&, int irrep=-1);
      PsiReturnType fd_hessian_0(Options &, const boost::python::list&);
    }

    extern int read_options(const std::string &name, Options & options, bool suppress_printing = false);
    extern FILE *outfile;
}

void py_flush_outfile()
{
    fflush(outfile);
}

void py_close_outfile()
{
    if (outfile != stdout) {
        fclose(outfile);
        outfile = NULL;
    }
}

void py_reopen_outfile()
{
    if (outfile_name == "stdout")
        outfile = stdout;
    else {
        outfile = fopen(outfile_name.c_str(), "a");
        if (outfile == NULL)
            throw PSIEXCEPTION("PSI4: Unable to reopen output file.");
    }
}

std::string py_get_outfile_name()
{
    return outfile_name;
}

void py_psi_prepare_options_for_module(std::string const & name)
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

int py_psi_stability()
{
    py_psi_prepare_options_for_module("STABILITY");
    return stable::stability(Process::environment.options);
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

int py_psi_mints()
{
    py_psi_prepare_options_for_module("MINTS");
    return mints::mints(Process::environment.options);
}

int py_psi_scfgrad()
{
    py_psi_prepare_options_for_module("SCF");
    return scfgrad::scfgrad(Process::environment.options);
}

int py_psi_deriv()
{
    py_psi_prepare_options_for_module("DERIV");
    return deriv::deriv(Process::environment.options);
}

int py_psi_omp2()
{
    py_psi_prepare_options_for_module("OMP2");
    return omp2wave::omp2wave(Process::environment.options);
}

int py_psi_libfock()
{
    py_psi_prepare_options_for_module("CPHF");
    return libfock::libfock(Process::environment.options);
}

double py_psi_scf_callbacks(PyObject* precallback, PyObject* postcallback)
{
    py_psi_prepare_options_for_module("SCF");
    if (scf::scf(Process::environment.options, precallback, postcallback) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

/*double py_psi_lmp2()
{
    py_psi_prepare_options_for_module("LMP2");
    if (lmp2::lmp2(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}
*/

double py_psi_mcscf()
{
    py_psi_prepare_options_for_module("MCSCF");
    if (mcscf::mcscf(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

PsiReturnType py_psi_mrcc_generate_input(const boost::python::dict& level)
{
    py_psi_prepare_options_for_module("MRCC");
    return mrcc::mrcc_generate_input(Process::environment.options, level);
}

PsiReturnType py_psi_mrcc_load_densities(const boost::python::dict& level)
{
    py_psi_prepare_options_for_module("MRCC");
    return mrcc::mrcc_load_ccdensities(Process::environment.options, level);
}

std::vector< SharedMatrix > py_psi_fd_geoms_1_0()
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_geoms_1_0(Process::environment.options);
}

//std::vector< boost::shared_ptr<Matrix> > py_psi_fd_geoms_2_0()
//{
    //py_psi_prepare_options_for_module("FINDIF");
    //return findif::fd_geoms_2_0(Process::environment.options);
//}

std::vector< SharedMatrix > py_psi_fd_geoms_freq_0(int irrep)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_geoms_freq_0(Process::environment.options, irrep);
}

std::vector< SharedMatrix > py_psi_fd_geoms_hessian_0()
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_geoms_hessian_0(Process::environment.options);
}

std::vector< SharedMatrix > py_psi_fd_geoms_freq_1(int irrep)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_geoms_freq_1(Process::environment.options, irrep);
}

PsiReturnType py_psi_fd_1_0(const boost::python::list& energies)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_1_0(Process::environment.options, energies);
}

//PsiReturnType py_psi_fd_2_0(const boost::python::list& energies)
//{
    //py_psi_prepare_options_for_module("FINDIF");
    //return findif::fd_2_0(Process::environment.options, energies);
//}

PsiReturnType py_psi_fd_freq_0(const boost::python::list& energies, int irrep)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_freq_0(Process::environment.options, energies, irrep);
}

PsiReturnType py_psi_fd_hessian_0(const boost::python::list& energies)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_hessian_0(Process::environment.options, energies);
}

PsiReturnType py_psi_fd_freq_1(const boost::python::list& grads, int irrep)
{
    py_psi_prepare_options_for_module("FINDIF");
    return findif::fd_freq_1(Process::environment.options, grads, irrep);
}

double py_psi_scf()
{
    return py_psi_scf_callbacks(Py_None, Py_None);
}

double py_psi_dcft()
{
    py_psi_prepare_options_for_module("DCFT");
    if (dcft::dcft(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

double py_psi_dfmp2()
{
    py_psi_prepare_options_for_module("DFMP2");
    if (dfmp2::dfmp2(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

double py_psi_dfmp2grad()
{
    py_psi_prepare_options_for_module("DFMP2");
    if (dfmp2::dfmp2grad(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

double py_psi_dfcc()
{
    py_psi_prepare_options_for_module("DFCC");
    if (dfcc::dfcc(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

double py_psi_sapt()
{
    py_psi_prepare_options_for_module("SAPT");
    if (sapt::sapt(Process::environment.options) == Success) {
        return Process::environment.globals["SAPT ENERGY"];
    }
    else
        return 0.0;
}

double py_psi_transqt()
{
    py_psi_prepare_options_for_module("TRANSQT");
    transqt::transqt(Process::environment.options);
    return 0.0;
}

double py_psi_transqt2()
{
    py_psi_prepare_options_for_module("TRANSQT2");
    transqt2::transqt2(Process::environment.options);
    return 0.0;
}

double py_psi_ccsort()
{
    py_psi_prepare_options_for_module("CCSORT");
    ccsort::ccsort(Process::environment.options);
    return 0.0;
}

double py_psi_ccenergy()
{
    py_psi_prepare_options_for_module("CCENERGY");
    boost::shared_ptr<Wavefunction> ccwave(new ccenergy::CCEnergyWavefunction(
                                               Process::environment.reference_wavefunction(),
                                               Process::environment.options)
                                           );
    Process::environment.set_reference_wavefunction(ccwave);

    double energy = ccwave->compute_energy();
    return energy;

//    if (ccenergy::ccenergy(Process::environment.options) == Success) {
//        return Process::environment.globals["CURRENT ENERGY"];
//    }
//    else
//        return 0.0;
}

double py_psi_mp2()
{
    py_psi_prepare_options_for_module("MP2");
    boost::shared_ptr<Wavefunction> mp2wave(new mp2::MP2Wavefunction(
                                                Process::environment.reference_wavefunction(),
                                                Process::environment.options)
                                           );
    Process::environment.set_reference_wavefunction(mp2wave);

    double energy = mp2wave->compute_energy();
    return energy;
}

double py_psi_cctriples()
{
    py_psi_prepare_options_for_module("CCTRIPLES");
    if (cctriples::cctriples(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

double py_psi_cepa()
{
    py_psi_prepare_options_for_module("CEPA");
    if (cepa::cepa(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}
double py_psi_detci()
{
    py_psi_prepare_options_for_module("DETCI");
    if (detci::detci(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

double py_psi_cchbar()
{
    py_psi_prepare_options_for_module("CCHBAR");
    cchbar::cchbar(Process::environment.options);
    return 0.0;
}

double py_psi_cclambda()
{
    py_psi_prepare_options_for_module("CCLAMBDA");
    cclambda::cclambda(Process::environment.options);
    return 0.0;
}

double py_psi_ccdensity()
{
    py_psi_prepare_options_for_module("CCDENSITY");
    ccdensity::ccdensity(Process::environment.options);
    return 0.0;
}

double py_psi_ccresponse()
{
    py_psi_prepare_options_for_module("CCRESPONSE");
    ccresponse::ccresponse(Process::environment.options);
    return 0.0;
}

double py_psi_cceom()
{
    py_psi_prepare_options_for_module("CCEOM");
    if (cceom::cceom(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

double py_psi_psimrcc()
{
    py_psi_prepare_options_for_module("PSIMRCC");
    psimrcc::psimrcc(Process::environment.options);
    return 0.0;
}

double py_psi_adc()
{
    py_psi_prepare_options_for_module("ADC");
    if (adc::adc(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
}

char const* py_psi_version()
{
    return PSI_VERSION;
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
    fprintf(outfile,"%s",s.c_str());
}

/**
 * @return whether key describes a convergence threshold or not
 */
bool specifies_convergence(std::string const & key){
    return ((key.find("CONV") != key.npos) || (key.find("TOL") != key.npos));
}

bool check_for_basis(std::string const & name, std::string const & type)
{
    if (type.find("BASIS") != type.npos) {
        // check the exceptions
        if (type == "BASIS_PATH" ||
            type == "AO_BASIS" ||
            type == "DUAL_BASIS")
            return false;

        // else set the basis for all atoms.
        Process::environment.molecule()->set_basis_all_atoms(name, type);
        return true;
    }
    return false;
}

bool py_psi_set_option_string(std::string const & module, std::string const & key, std::string const & value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options[nonconst_key];

    if (data.type() == "string") {
        Process::environment.options.set_str(module, nonconst_key, value);
        check_for_basis(value, nonconst_key);
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

bool py_psi_set_option_int(std::string const & module, std::string const & key, int value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options.use(nonconst_key);

    if(data.type() == "double" && specifies_convergence(nonconst_key)){
        double val = pow(10.0, -value);
        Process::environment.options.set_double(module, nonconst_key, val);
    }else if (data.type() == "boolean") {
        Process::environment.options.set_bool(module, nonconst_key, value ? true : false);
    }else{
        Process::environment.options.set_int(module, nonconst_key, value);
    }
    return true;

}

bool py_psi_set_option_double(std::string const & module, std::string const & key, double value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Process::environment.options.set_double(module, nonconst_key, value);
    return true;
}

template <class T>
bool is_int(T x)
{
   using boost::python::type_id;
   return type_id<T>() == type_id<int>();
}
template <class T>
bool is_double(T x)
{
   using boost::python::type_id;
   return type_id<T>() == type_id<double>();
}

bool py_psi_set_global_option_string(std::string const & key, std::string const & value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options.use(nonconst_key);

    if (data.type() == "string") {
        Process::environment.options.set_global_str(nonconst_key, value);
        check_for_basis(value, nonconst_key);
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

bool py_psi_set_global_option_int(std::string const & key, int value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options.use(nonconst_key);

    if( data.type() == "double" && specifies_convergence(nonconst_key)){
        double val = pow(10.0, -value);
        Process::environment.options.set_global_double(nonconst_key, val);
    }else if (data.type() == "boolean") {
        Process::environment.options.set_global_bool(nonconst_key, value ? true : false);
    }else{
        Process::environment.options.set_global_int(nonconst_key, value);
    }
    return true;
}

bool py_psi_set_global_option_double(std::string const & key, double value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Process::environment.options.set_global_double(nonconst_key, value);
    return true;
}

bool py_psi_set_global_option_python(std::string const & key, boost::python::object& obj)
{
    string nonconst_key = boost::to_upper_copy(key);
    Process::environment.options.set_global_python(nonconst_key, obj);
    return true;
}

bool py_psi_set_option_array(std::string const & module, std::string const & key, const python::list &values, DataType *entry = NULL)
{
    string nonconst_key = boost::to_upper_copy(key);
    // Assign a new head entry on the first time around only
    if(entry == NULL)
        Process::environment.options.set_array(module, nonconst_key);
    size_t size = len(values);
    for(int n = 0; n < size; ++n){
        extract<python::list> lval(values[n]);
        extract<std::string> sval(values[n]);
        extract<double> fval(values[n]);
        extract<int> ival(values[n]);
        if(lval.check()){
            python::list l = extract<python::list>(values[n]);
            DataType *newentry = Process::environment.options.set_local_array_array(module, nonconst_key, entry);
            // Now we need to recurse, to fill in the data
            py_psi_set_option_array(module, key, l, newentry);
        }else if(sval.check()){
            std::string s = extract<std::string>(values[n]);
            Process::environment.options.set_local_array_string(module, nonconst_key, s, entry);
        }else if(ival.check()){
            int i = extract<int>(values[n]);
            Process::environment.options.set_local_array_int(module, nonconst_key, i, entry);
        }else if(fval.check()){
            double f = extract<double>(values[n]);
            Process::environment.options.set_local_array_double(module, nonconst_key, f, entry);
        }
    }
    return true;
}

bool py_psi_set_global_option_array(std::string const & key, python::list values, DataType *entry=NULL)
{
    string nonconst_key = boost::to_upper_copy(key);
    // Assign a new head entry on the first time around only
    if(entry == NULL)
        Process::environment.options.set_global_array(nonconst_key);
    size_t size = len(values);
    for(int n = 0; n < size; ++n){
        extract<python::list> lval(values[n]);
        extract<std::string> sval(values[n]);
        extract<double> fval(values[n]);
        extract<int> ival(values[n]);
        if(lval.check()){
            python::list l = extract<python::list>(values[n]);
            DataType *newentry = Process::environment.options.set_global_array_array(nonconst_key, entry);
            // Now we need to recurse, to fill in the data
            py_psi_set_global_option_array(key, l, newentry);
        }else if(sval.check()){
            std::string s = extract<std::string>(values[n]);
            Process::environment.options.set_global_array_string(nonconst_key, s, entry);
        }else if(ival.check()){
            int i = extract<int>(values[n]);
            Process::environment.options.set_global_array_int(nonconst_key, i, entry);
        }else if(fval.check()){
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
        dynamic_cast<PythonDataType*>(data.get())->assign(obj);
    else
        throw PSIEXCEPTION("Unable to set option to a Python object.");
}

object py_psi_get_option(const string& key)
{
    string nonconst_key = key;
    Data& data = Process::environment.options.use(nonconst_key);

    if (data.type() == "string")
        return str(data.to_string());
    else if (data.type() == "boolean" || data.type() == "int")
        return object(data.to_integer());
    else if (data.type() == "double")
        return object(data.to_double());

    return object();
}

bool py_psi_has_option_changed(const string& key)
{
    string nonconst_key = key;
    Data& data = Process::environment.options.use(nonconst_key);

    return data.has_changed();
}

bool py_psi_has_global_option_changed(std::string const & key)
{
    string nonconst_key = key;
    Data& data = Process::environment.options.get_global(nonconst_key);

    return data.has_changed();
}

bool py_psi_has_local_option_changed(std::string const & module, std::string const & key)
{
    string nonconst_key = key;
    Process::environment.options.set_current_module(module);
    Data& data = Process::environment.options.use(nonconst_key);

    return data.has_changed();
}

void py_psi_revoke_option_changed(std::string const & key)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options.use(nonconst_key);
    data.dechanged();
}

void py_psi_revoke_global_option_changed(std::string const & key)
{
    string nonconst_key = boost::to_upper_copy(key);
    Data& data = Process::environment.options.get_global(nonconst_key);
    data.dechanged();
}

void py_psi_revoke_local_option_changed(std::string const & module, std::string const & key)
{
    string nonconst_key = boost::to_upper_copy(key);
    Process::environment.options.set_current_module(module);
    Data& data = Process::environment.options.use(nonconst_key);
    data.dechanged();
}

object py_psi_get_global_option(std::string const & key)
{
    string nonconst_key = key;
    Data& data = Process::environment.options.get_global(nonconst_key);

    if (data.type() == "string")
        return str(data.to_string());
    else if (data.type() == "boolean" || data.type() == "int")
        return object(data.to_integer());
    else if (data.type() == "double")
        return object(data.to_double());

    return object();
}

object py_psi_get_local_option(std::string const & module, std::string const & key)
{
    string nonconst_key = key;
    Process::environment.options.set_current_module(module);
    Data& data = Process::environment.options.use(nonconst_key);

    if (data.type() == "string")
        return str(data.to_string());
    else if (data.type() == "boolean" || data.type() == "int")
        return object(data.to_integer());
    else if (data.type() == "double")
        return object(data.to_double());

    return object();
}

void py_psi_set_active_molecule(boost::shared_ptr<Molecule> molecule)
{
    Process::environment.set_molecule(molecule);
}

void py_psi_set_parent_symmetry(std::string pg)
{
    boost::shared_ptr<PointGroup> group = boost::shared_ptr<PointGroup>();
    if(pg != ""){
        group = boost::shared_ptr<PointGroup>(new PointGroup(pg));
    }

    Process::environment.set_parent_symmetry(group);
}


boost::shared_ptr<Molecule> py_psi_get_active_molecule()
{
    return Process::environment.molecule();
}

void py_psi_set_gradient(SharedMatrix grad)
{
    if (Process::environment.reference_wavefunction()) {
        Process::environment.reference_wavefunction()->set_gradient(grad);
    } else {
        Process::environment.set_gradient(grad);
    }
}

SharedMatrix py_psi_get_gradient()
{
    if (Process::environment.reference_wavefunction()) {
        boost::shared_ptr<Wavefunction> wf = Process::environment.reference_wavefunction();
        return wf->gradient();
    } else {
        return Process::environment.gradient();
    }
}

double py_psi_get_variable(const std::string & key)
{
    string uppercase_key = key;
    transform(uppercase_key.begin(), uppercase_key.end(), uppercase_key.begin(), ::toupper);
    return Process::environment.globals[uppercase_key];
}

void py_psi_set_variable(const std::string & key, double val)
{
    string uppercase_key = key;
    transform(uppercase_key.begin(), uppercase_key.end(), uppercase_key.begin(), ::toupper);
    Process::environment.globals[uppercase_key] = val;
}

void py_psi_set_memory(unsigned long int mem)
{
    Process::environment.set_memory(mem);
    fprintf(outfile,"\n  Memory set to %7.3f %s by Python script.\n",(mem > 1000000000 ? mem/1.0E9 : mem/1.0E6), \
        (mem > 1000000000 ? "GiB" : "MiB" ));
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

int py_psi_get_nproc()
{
    return Communicator::world->nproc();
}

int py_psi_get_me()
{
    return Communicator::world->me();
}

boost::shared_ptr<Wavefunction> py_psi_reference_wavefunction()
{
    return Process::environment.reference_wavefunction();
}

void py_psi_add_user_specified_basis_file(const string& file)
{
    Process::environment.user_basis_files.push_front(file);
}

string py_psi_get_input_directory()
{
    return infile_directory;
}

void py_psi_print_variable_map()
{
    int largest_key = 0;
    for ( std::map<string,double>::iterator it = Process::environment.globals.begin() ; it != Process::environment.globals.end(); it++ ) {
        if (it->first.size() > largest_key)
            largest_key = it->first.size();
    }
    largest_key += 2;  // for quotation marks

    std::stringstream line;
    std::string first_tmp;
    for ( std::map<string,double>::iterator it = Process::environment.globals.begin() ; it != Process::environment.globals.end(); it++ ) {
        first_tmp = "\"" + it->first + "\"";
        line << "  " << std::left << std::setw(largest_key) << first_tmp << " => " << std::setw(20) << std::right <<
            std::fixed << std::setprecision(12) << it->second << std::endl;
    }

    fprintf(outfile, "\n\n  Variable Map:");
    fprintf(outfile, "\n  ----------------------------------------------------------------------------\n");
    fprintf(outfile, "%s\n\n", line.str().c_str());
}

std::string py_psi_top_srcdir()
{
    return PSI_TOP_SRCDIR;
}

// Tell python about the default final argument to the array setting functions
BOOST_PYTHON_FUNCTION_OVERLOADS(set_global_option_overloads, py_psi_set_global_option_array, 2, 3)
BOOST_PYTHON_FUNCTION_OVERLOADS(set_local_option_overloads, py_psi_set_option_array, 3, 4)

BOOST_PYTHON_MODULE(PsiMod)
{
    docstring_options sphx_doc_options(true, true, false);

    enum_<PsiReturnType>("PsiReturnType", "docstring")
            .value("Success", Success)
            .value("Failure", Failure)
            .value("Balk", Balk)
            .value("EndLoop", EndLoop)
            .export_values();

    def("version", py_psi_version, "docstring");
    def("clean", py_psi_clean, "Function to remove scratch files. Call between independent jobs.");

    // Benchmarks
    export_benchmarks();

    // BLAS/LAPACK Static Wrappers
    export_blas_lapack();

    // Plugins
    export_plugins();

    // OEProp/GridProp
    export_oeprop();


    // Options
    def("prepare_options_for_module", py_psi_prepare_options_for_module, "docstring");
    def("set_active_molecule", py_psi_set_active_molecule, "docstring");
    def("get_active_molecule", &py_psi_get_active_molecule, "docstring");
    def("reference_wavefunction", py_psi_reference_wavefunction, "docstring");
    def("get_gradient", py_psi_get_gradient, "docstring");
    def("set_gradient", py_psi_set_gradient, "docstring");
    def("set_memory", py_psi_set_memory, "docstring");
    def("get_memory", py_psi_get_memory, "docstring");
    def("set_nthread", &py_psi_set_n_threads, "docstring");
    def("nthread", &py_psi_get_n_threads, "docstring");
    def("nproc", &py_psi_get_nproc, "docstring");
    def("me", &py_psi_get_me, "docstring");

    def("set_parent_symmetry", py_psi_set_parent_symmetry, "docstring");
    def("print_options", py_psi_print_options, "docstring");
    def("print_global_options", py_psi_print_global_options, "docstring");
    def("print_out", py_psi_print_out, "docstring");

    // Set the different local option types
    def("set_local_option", py_psi_set_option_string, "docstring");
    def("set_local_option", py_psi_set_option_double, "docstring");
    def("set_local_option", py_psi_set_option_int, "docstring");
    def("set_local_option", py_psi_set_option_array, set_local_option_overloads());

    // Set the different global option types
    def("set_global_option", py_psi_set_global_option_string, "docstring");
    def("set_global_option", py_psi_set_global_option_double, "docstring");
    def("set_global_option", py_psi_set_global_option_int, "docstring");
    def("set_global_option", py_psi_set_global_option_array, set_global_option_overloads());

    def("set_global_option_python", py_psi_set_global_option_python, "docstring");
    def("set_local_option_python", py_psi_set_local_option_python, "docstring");

    def("get_global_option_list", py_psi_get_global_option_list, "docstring");

    // Get the option; letting liboptions decide whether to use global or local
    def("get_option", py_psi_get_option, "docstring");

    // Get the option; specify whether to use global or local
    def("get_global_option", py_psi_get_global_option, "docstring");
    def("get_local_option", py_psi_get_local_option, "docstring");

    // Returns whether the option has changed/revoke has changed for silent resets
    def("has_option_changed", py_psi_has_option_changed, "docstring");
    def("has_global_option_changed", py_psi_has_global_option_changed, "docstring");
    def("has_local_option_changed", py_psi_has_local_option_changed, "docstring");
    def("revoke_option_changed", py_psi_revoke_option_changed, "docstring");
    def("revoke_global_option_changed", py_psi_revoke_global_option_changed, "docstring");
    def("revoke_local_option_changed", py_psi_revoke_local_option_changed, "docstring");

    // These return/set variable value found in Process::environment.globals
    def("get_variable", py_psi_get_variable, "docstring");
    def("set_variable", py_psi_set_variable, "docstring");

    // Print the variables found in Process::environment.globals
    def("print_variables", py_psi_print_variable_map, "docstring");

    // Adds a custom user basis set file.
    def("add_user_basis_file", py_psi_add_user_specified_basis_file, "docstring");

    // Get the name of the directory where the input file is at
    def("get_input_directory", py_psi_get_input_directory, "docstring");

    // Returns the location where the Psi4 source is located.
    def("psi_top_srcdir", py_psi_top_srcdir, "docstring");

    def("flush_outfile", py_flush_outfile, "docstring");
    def("close_outfile", py_close_outfile, "docstring");
    def("reopen_outfile", py_reopen_outfile, "docstring");
    def("outfile_name", py_get_outfile_name, "docstring");

    // modules
    def("mints", py_psi_mints, "docstring");
    def("deriv", py_psi_deriv, "docstring");
    def("scfgrad", py_psi_scfgrad, "docstring");

    typedef double (*scf_module_none)();
    typedef double (*scf_module_two)(PyObject*, PyObject*);

    def("scf", py_psi_scf_callbacks, "docstring");
    def("scf", py_psi_scf, "docstring");
    def("dcft", py_psi_dcft, "docstring");
    def("libfock", py_psi_libfock, "docstring");
    def("dfmp2", py_psi_dfmp2, "docstring");
    def("dfmp2grad", py_psi_dfmp2grad, "docstring");
    def("dfcc", py_psi_dfcc, "docstring");
//    def("lmp2", py_psi_lmp2, "docstring");
    def("mp2", py_psi_mp2, "docstring");
    def("mcscf", py_psi_mcscf, "docstring");
    def("mrcc_generate_input", py_psi_mrcc_generate_input, "docstring");
    def("mrcc_load_densities", py_psi_mrcc_load_densities, "docstring");
    def("fd_geoms_1_0", py_psi_fd_geoms_1_0, "docstring");
    //def("fd_geoms_2_0", py_psi_fd_geoms_2_0);
    def("fd_geoms_freq_0", py_psi_fd_geoms_freq_0, "docstring");
    def("fd_geoms_freq_1", py_psi_fd_geoms_freq_1, "docstring");
    def("fd_geoms_hessian_0", py_psi_fd_geoms_hessian_0, "docstring");
    def("fd_1_0", py_psi_fd_1_0, "docstring");
    //def("fd_2_0", py_psi_fd_2_0);
    def("fd_freq_0", py_psi_fd_freq_0, "docstring");
    def("fd_freq_1", py_psi_fd_freq_1, "docstring");
    def("fd_hessian_0", py_psi_fd_hessian_0, "docstring");
    def("sapt", py_psi_sapt, "docstring");
    def("stability", py_psi_stability, "docstring");
    def("psimrcc", py_psi_psimrcc, "docstring");
    def("optking", py_psi_optking, "docstring");
    def("transqt", py_psi_transqt, "docstring");
    def("transqt2", py_psi_transqt2, "docstring");
    def("ccsort", py_psi_ccsort, "docstring");
    def("ccenergy", py_psi_ccenergy, "docstring");
    def("cctriples", py_psi_cctriples, "docstring");
    def("detci", py_psi_detci, "docstring");
    def("cepa", py_psi_cepa, "docstring");
    def("cchbar", py_psi_cchbar, "docstring");
    def("cclambda", py_psi_cclambda, "docstring");
    def("ccdensity", py_psi_ccdensity, "docstring");
    def("ccresponse", py_psi_ccresponse, "docstring");
    def("cceom", py_psi_cceom, "docstring");
    def("omp2", py_psi_omp2, "docstring");
    def("adc", py_psi_adc, "docstring");
    def("opt_clean", py_psi_opt_clean, "docstring");

    // Define library classes
    export_psio();
    export_chkpt();
    export_mints();
    export_functional();

    typedef string (Process::Environment::*environmentStringFunction)(const string&);

    class_<Process::Environment>("Environment").
        def("__getitem__", environmentStringFunction(&Process::Environment::operator ()), "docstring");

    typedef string (Process::Arguments::*argumentsStringFunction)(int);

    class_<Process::Arguments>("Arguments").
        def("__getitem__", argumentsStringFunction(&Process::Arguments::operator ()), "docstring");

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
    Py_Finalize();
}

#define PY_TRY(ptr, command)  \
     if(!(ptr = command)){    \
         PyErr_Print();       \
         exit(1);             \
     }


void Python::run(FILE *input)
{
    using namespace boost::python;
    char *s = 0;
    if (input == NULL)
        return;

    // Setup globals options
    Process::environment.options.set_read_globals(true);
    read_options("", Process::environment.options, true);
    Process::environment.options.set_read_globals(false);

    if (!Py_IsInitialized()) {
        s = strdup("psi");
        // Py_InitializeEx(0) causes sig handlers to not be installed.
        Py_InitializeEx(0);
        #if PY_VERSION_HEX >= 0x03000000
        Py_SetProgramName(L"psi");
        #else
        Py_SetProgramName(s);
        #endif

        // Track down the location of PSI4's python script directory.
        std::string psiDataDirName = Process::environment("PSIDATADIR");
        std::string psiDataDirWithPython = psiDataDirName + "/python";
        boost::filesystem::path bf_path;
        bf_path = boost::filesystem::system_complete(psiDataDirWithPython);
        if(!boost::filesystem::is_directory(bf_path)) {
            printf("Unable to read the PSI4 Python folder - check the PSIDATADIR environmental variable\n"
                    "      Current value of PSIDATADIR is %s\n", psiDataDirName.c_str());
            exit(1);
        }

        // Add PSI library python path
        PyObject *path, *sysmod, *str;
        PY_TRY(sysmod , PyImport_ImportModule("sys"));
        PY_TRY(path   , PyObject_GetAttrString(sysmod, "path"));
        PY_TRY(str    , PyString_FromString(psiDataDirWithPython.c_str()));
        PyList_Append(path, str);
        Py_DECREF(str);
        Py_DECREF(path);
        Py_DECREF(sysmod);
    }
    if (Py_IsInitialized()) {
        // Stupid way to read in entire file.
        char line[256];
        std::stringstream file;
        while(fgets(line, sizeof(line), input)) {
            file << line;
        }

        try {
            PyImport_AppendInittab(strdup("PsiMod"), initPsiMod);
            object objectMain(handle<>(borrowed(PyImport_AddModule("__main__"))));
            object objectDict = objectMain.attr("__dict__");
            s = strdup("import PsiMod");
            PyRun_SimpleString(s);

            // Process the input file
            PyObject *input;
            PY_TRY(input, PyImport_ImportModule("input") );
            PyObject *function;
            PY_TRY(function, PyObject_GetAttrString(input, "process_input"));
            PyObject *pargs;
            PY_TRY(pargs, Py_BuildValue("(s)", file.str().c_str()) );
            PyObject *ret;
            PY_TRY( ret, PyEval_CallObject(function, pargs) );

            char *val;
            PyArg_Parse(ret, "s", &val);
            string inputfile = val;

            Py_DECREF(ret);
            Py_DECREF(pargs);
            Py_DECREF(function);
            Py_DECREF(input);

            if (verbose) {
                fprintf(outfile, "\n Input file to run:\n%s", inputfile.c_str());
                fflush(outfile);
            }

            str strStartScript(inputfile);

            object objectScriptInit = exec( strStartScript, objectDict, objectDict );
        }
        catch (error_already_set const& e)
        {
            PyErr_Print();
            exit(1);
        }
    }
    else {
        fprintf(stderr, "Unable to run Python input file.\n");
        return;
    }

    if (s)
      free(s);
    Process::environment.reference_wavefunction().reset();
    py_psi_plugin_close_all();
}
