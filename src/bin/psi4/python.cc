#include <boost/algorithm/string.hpp>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>
#include <cstdio>
#include <sstream>
#include <map>

#include <psiconfig.h>
#include <libmints/mints.h>
#include <libplugin/plugin.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>

#include <psi4-dec.h>
#include "script.h"
#include "psi4.h"

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

namespace opt      { psi::PsiReturnType optking(psi::Options &); }

namespace psi {
    namespace mints     { PsiReturnType mints(Options &); }
    namespace deriv     { PsiReturnType deriv(Options &); }
    namespace scf       { PsiReturnType scf(Options &, PyObject* pre, PyObject* post);   }
    namespace dfmp2     { PsiReturnType dfmp2(Options &); }
    namespace dfcc      { PsiReturnType dfcc(Options &); }
    namespace sapt      { PsiReturnType sapt(Options &);  }
    namespace dcft      { PsiReturnType dcft(Options &);  }
    namespace mcscf     { PsiReturnType mcscf(Options &);  }

    namespace psimrcc   { PsiReturnType psimrcc(Options &);  }
    namespace transqt   { PsiReturnType transqt(Options &);  }
    namespace transqt2  { PsiReturnType transqt2(Options &); }
    namespace ccsort    { PsiReturnType ccsort(Options&);    }
    namespace ccenergy  { PsiReturnType ccenergy(Options&);  }
    namespace cctriples { PsiReturnType cctriples(Options&); }
    namespace cchbar    { PsiReturnType cchbar(Options&);    }
    namespace cclambda  { PsiReturnType cclambda(Options&);  }
    namespace ccdensity { PsiReturnType ccdensity(Options&); }
    // DETCI: uncomment
    //namespace detci     { PsiReturnType detci(Options&);     }

    extern int read_options(const std::string &name, Options & options, bool suppress_printing = false);
    extern FILE *outfile;
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

int py_psi_optking()
{
    py_psi_prepare_options_for_module("OPTKING");
    return opt::optking(Process::environment.options);
}

int py_psi_mints()
{
    py_psi_prepare_options_for_module("MINTS");
    return mints::mints(Process::environment.options);
}

int py_psi_deriv()
{
    py_psi_prepare_options_for_module("DERIV");
    return deriv::deriv(Process::environment.options);
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

double py_psi_mcscf()
{
    py_psi_prepare_options_for_module("MCSCF");
    if (mcscf::mcscf(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
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
    if (ccenergy::ccenergy(Process::environment.options) == Success) {
        return Process::environment.globals["CURRENT ENERGY"];
    }
    else
        return 0.0;
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

double py_psi_detci()
{
    py_psi_prepare_options_for_module("DETCI");

    // DETCI: Uncomment
    //if (detci::detci(Process::environment.options) == Success) {
    //    return Process::environment.globals["CURRENT ENERGY"];
    //}
    //else
    //    return 0.0;
    fprintf(outfile,"\n\nWorld's slowest quantum method goes here.\n\n");
    fflush(outfile);

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

double py_psi_psimrcc()
{
    py_psi_prepare_options_for_module("PSIMRCC");
    psimrcc::psimrcc(Process::environment.options);
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

void py_psi_print_out(std::string s)
{
    fprintf(outfile,"%s",s.c_str());
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

    if (data.type() == "boolean") {
        Process::environment.options.set_bool(module, nonconst_key, value ? true : false);
    }else{
        Process::environment.options.set_int(module, nonconst_key, value);
    }
    return true;

}

bool py_psi_set_option_float(std::string const & module, std::string const & key, float value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Process::environment.options.set_double(module, nonconst_key, value);
    return true;
}

// Right now this can only handle arrays of integers.
// Unable to handle strings.
bool py_psi_set_option_array(std::string const & module, std::string const & key, const python::list &values)
{
    string nonconst_key = boost::to_upper_copy(key);
    std::vector<double> vector;
    size_t n = len(values);
    for(int i = 0; i < n; ++i) vector.push_back(extract<double>(values[i]));

    Process::environment.options.set_array(module, nonconst_key, vector);
    return true;
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

    if (data.type() == "boolean") {
        Process::environment.options.set_global_bool(nonconst_key, value ? true : false);
    }else{
        Process::environment.options.set_global_int(nonconst_key, value);
    }
    return true;
}

bool py_psi_set_global_option_float(std::string const & key, float value)
{
    string nonconst_key = boost::to_upper_copy(key);
    Process::environment.options.set_global_double(nonconst_key, value);
    return true;
}

// Right now this can only handle arrays of integers.
// Unable to handle strings.
bool py_psi_set_global_option_array(std::string const & key, python::list values)
{
    string nonconst_key = boost::to_upper_copy(key);

    std::vector<double> vector;
    size_t n = len(values);
    for(int i = 0; i < n; ++i) vector.push_back(extract<double>(values[i]));

    Process::environment.options.set_global_array(nonconst_key, vector);
    return true;
}

object py_psi_get_option(const string& key)
{
    string nonconst_key = key;
    Data& data = Process::environment.options.use(nonconst_key);

    if (data.type() == "string")
        return str(data.to_string());
    else if (data.type() == "boolean" || data.type() == "integer")
        return object(data.to_integer());
    else if (data.type() == "double")
        return object(data.to_double());

    return object();
}

object py_psi_get_global_option(const string& key)
{
    string nonconst_key = key;
    Data& data = Process::environment.options.use(nonconst_key);

    if (data.type() == "string")
        return str(data.to_string());
    else if (data.type() == "boolean" || data.type() == "integer")
        return object(data.to_integer());
    else if (data.type() == "double")
        return object(data.to_double());

    return object();
}

void py_psi_set_active_molecule(boost::shared_ptr<Molecule> molecule)
{
    Process::environment.set_molecule(molecule);
}

boost::shared_ptr<Molecule> py_psi_get_active_molecule()
{
    return Process::environment.molecule();
}

double py_psi_get_variable(const std::string & key)
{
    string uppercase_key = key;
    transform(uppercase_key.begin(), uppercase_key.end(), uppercase_key.begin(), ::toupper);
    return Process::environment.globals[uppercase_key];
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
    for ( std::map<string,double>::iterator it = Process::environment.globals.begin() ; it != Process::environment.globals.end(); it++ ) {
        fprintf(outfile, "\"%s\" => %20.15f\n", it->first.c_str(), (double)it->second);
    }
}

BOOST_PYTHON_MODULE(PsiMod)
{
    enum_<PsiReturnType>("PsiReturnType")
            .value("Success", Success)
            .value("Failure", Failure)
            .value("Balk", Balk)
            .value("EndLoop", EndLoop)
            .export_values();

    def("version", py_psi_version);
    def("clean", py_psi_clean);

    // Benchmarks
    export_benchmarks();

    // BLAS/LAPACK Static Wrappers
    export_blas_lapack();

    // Plugins
    export_plugins();

    // OEProp/GridProp
    export_oeprop();

    // Options
    def("prepare_options_for_module", py_psi_prepare_options_for_module);
    def("set_active_molecule", py_psi_set_active_molecule);
    def("get_active_molecule", &py_psi_get_active_molecule);
    def("reference_wavefunction", py_psi_reference_wavefunction);
    def("set_memory", py_psi_set_memory);
    def("get_memory", py_psi_get_memory);
    def("set_n_threads", &py_psi_set_n_threads);
    def("get_n_threads", &py_psi_get_n_threads);

    def("print_options", py_psi_print_options);
    def("print_global_options", py_psi_print_global_options);
    def("print_out", py_psi_print_out);

    def("set_local_option", py_psi_set_option_string);
    def("set_local_option", py_psi_set_option_float);
    def("set_local_option", py_psi_set_option_int);
    def("set_local_option", py_psi_set_option_array);

    def("set_global_option", py_psi_set_global_option_string);
    def("set_global_option", py_psi_set_global_option_float);
    def("set_global_option", py_psi_set_global_option_int);
    def("set_global_option", py_psi_set_global_option_array);

    def("get_option", py_psi_get_option);
    def("get_global_option", py_psi_get_global_option);

    def("get_variable", py_psi_get_variable);
    def("print_variables", py_psi_print_variable_map);

    def("add_user_basis_file", py_psi_add_user_specified_basis_file);

    def("get_input_directory", py_psi_get_input_directory);

    // modules
    def("mints", py_psi_mints);
    def("deriv", py_psi_deriv);

    typedef double (*scf_module_none)();
    typedef double (*scf_module_two)(PyObject*, PyObject*);

    def("scf", py_psi_scf_callbacks);
    def("scf", py_psi_scf);
    def("dcft", py_psi_dcft);
    def("dfmp2", py_psi_dfmp2);
    def("dfcc", py_psi_dfcc);
    def("mcscf", py_psi_mcscf);
    def("sapt", py_psi_sapt);
    def("psimrcc", py_psi_psimrcc);
    def("optking", py_psi_optking);
    def("transqt", py_psi_transqt);
    def("transqt2", py_psi_transqt2);
    def("ccsort", py_psi_ccsort);
    def("ccenergy", py_psi_ccenergy);
    def("cctriples", py_psi_cctriples);
    def("detci", py_psi_detci);
    def("cchbar", py_psi_cchbar);
    def("cclambda", py_psi_cclambda);
    def("ccdensity", py_psi_ccdensity);

    // Define library classes
    export_psio();
    export_chkpt();
    export_mints();
    export_functional();

    /**
    class_<DFTensor, boost::shared_ptr<DFTensor> >( "DFTensor", no_init).
        def("bootstrap", &DFTensor::bootstrap_DFTensor).
        staticmethod("booststrap").
        def("form_fitting_metric", &DFTensor::form_fitting_metric).
        def("form_cholesky_metric", &DFTensor::form_cholesky_metric).
        def("form_qr_metric", &DFTensor::form_qr_metric).
        def("finalize", &DFTensor::finalize).
        def("print_out", &DFTensor::print_python);
    **/

    typedef string (Process::Environment::*environmentStringFunction)(const string&);

    class_<Process::Environment>("Environment").
        def("__getitem__", environmentStringFunction(&Process::Environment::operator ()));
//        def("set", &Process::Environment::set);

    typedef string (Process::Arguments::*argumentsStringFunction)(int);

    class_<Process::Arguments>("Arguments").
        def("__getitem__", argumentsStringFunction(&Process::Arguments::operator ()));

    class_<Process>("Process").
        add_static_property("environment", &Process::environment);

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
}

void Python::run(FILE *input)
{
    using namespace boost::python;
    char *s;
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
        std::string psiDataDirName = Process::environment("PSIDATADIR") + "/python";
        boost::filesystem::path bf_path;
        bf_path = boost::filesystem::system_complete(psiDataDirName);
        if(!boost::filesystem::is_directory(bf_path)) {
            throw PSIEXCEPTION("Unable to read the PSI4 Python folder - check the PSIDATADIR environmental variable\n"
                               "      Current value of PSIDATADIR is " + Process::environment("PSIDATADIR"));
        }

        // Add PSI library python path
        PyObject *path, *sysmod, *str;
        sysmod = PyImport_ImportModule("sys");
        path = PyObject_GetAttrString(sysmod, "path");
        str = PyString_FromString(psiDataDirName.c_str());
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
            PyObject *input = PyImport_ImportModule("input");
            PyObject *function = PyObject_GetAttrString(input, "process_input");
            PyObject *pargs = Py_BuildValue("(s)", file.str().c_str());
            PyObject *ret = PyEval_CallObject(function, pargs);

            char *val;
            PyArg_Parse(ret, "s", &val);
            string inputfile = val;

            Py_DECREF(ret);
            Py_DECREF(pargs);
            Py_DECREF(function);
            Py_DECREF(input);


            if (verbose) {
                fprintf(outfile, "Input file to run:\n%s", inputfile.c_str());
                fflush(outfile);
            }

            str strStartScript(inputfile);

            object objectScriptInit = exec( strStartScript, objectDict, objectDict );
            Communicator::world->sync();
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

    py_psi_plugin_close_all();
}
