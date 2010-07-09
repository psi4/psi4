// Include Python's header file.
// Use "python-config --includes" to determine this location.
#include <cstdio>
#include <boost/python.hpp>
#include <psiconfig.h>
#include <sstream>
#include "script.h"
#include "libchkpt/chkpt.hpp"
#include "liboptions/liboptions.h"
#include <libmints/molecule.h>
#include <libmints/pointgrp.h>

#include <libpsio/psio.hpp>

using namespace psi;
using namespace boost;
using namespace boost::python;
using namespace std;

namespace psi {
    extern void psiclean(void);
    extern FILE *outfile;
    extern Options options;
}

bool py_psi_input()
{
    /* Apply any options that need to be to global options object.

       If desired, set all the default values using read_options.
       read_options would need to be modified to not use ipv1.

       read_options("INPUT", options);

       Override options, if needed.
       options.set_str("UNITS", "BOHR");
     */

    /* Need to modify input to not take argc and argv
       input::input3(options);

       or

       dispatch_table["INPUT"](options);
     */
    printf("input: I did absolutely nothing.\n");
    return true;
}

char const* py_psi_version()
{
    return PSI_VERSION;
}

void py_psi_clean()
{
    psiclean();
}

bool py_psi_configure_psio(PSIO* obj)
{
    return psiopp_ipv1_config(obj);
}

BOOST_PYTHON_MODULE(PsiMod)
{
    def("version", py_psi_version);
    def("clean", py_psi_clean);
    def("configure_io", py_psi_configure_psio);

    // modules
    def("input",py_psi_input);

    // Define library classes
    class_<PSIO, shared_ptr<PSIO> >( "IO" ).
        def( "state", &PSIO::state ).
        def( "open", &PSIO::open ).
        def( "close", &PSIO::close ).
        def( "rehash", &PSIO::rehash ).
        def( "open_check", &PSIO::open_check ).
        def( "tocclean", &PSIO::tocclean ).
        def( "tocprint", &PSIO::tocprint ).
        def( "tocwrite", &PSIO::tocwrite ).
        def( "sharedObject", &PSIO::shared_object).
        staticmethod("sharedObject");

    class_<Chkpt, shared_ptr<Chkpt> >( "Checkpoint", init<PSIO*, int>() ).
        add_property( "enuc", &Chkpt::rd_enuc, &Chkpt::wt_enuc).
        add_property( "label", &Chkpt::rd_label, &Chkpt::wt_label).
        add_property( "escf", &Chkpt::rd_escf, &Chkpt::wt_escf).
        add_property( "eref", &Chkpt::rd_eref, &Chkpt::wt_eref).
        add_property( "ecorr", &Chkpt::rd_ecorr, &Chkpt::wt_ecorr).
        add_property( "efzc", &Chkpt::rd_efzc, &Chkpt::wt_efzc).
        add_property( "etot", &Chkpt::rd_etot, &Chkpt::wt_etot).
        add_property( "disp", &Chkpt::rd_disp, &Chkpt::wt_disp).
        add_property( "eccsd", &Chkpt::rd_eccsd, &Chkpt::wt_eccsd).
        add_property( "e_t", &Chkpt::rd_e_t, &Chkpt::wt_e_t).
        add_property( "emp2", &Chkpt::rd_emp2, &Chkpt::wt_emp2).
        def( "sharedObject", &Chkpt::shared_object).
        staticmethod("sharedObject");

    class_<Vector3>("Vector3").
        def(init<double>()).
        def(init<double, double, double>()).
        def(init<const Vector3&>()).
//      def(self = other<double>()).
        def(self += self).
        def(self -= self).
        def(self *= other<double>()).
        def(self + self).
        def(self - self).
        def(-self).
        def("dot", &Vector3::dot).
        def("distance", &Vector3::distance).
        def("normalize", &Vector3::normalize).
        def("norm", &Vector3::norm).
        def("cross", &Vector3::cross).
        def("__str__", &Vector3::to_string).
        def("__getitem__", &Vector3::get);

    typedef string (Process::Environment::*environmentStringFunction)(const string&);

    class_<Process::Environment>("Environment").
        def("__getitem__", environmentStringFunction(&Process::Environment::operator ()));
//        def("set", &Process::Environment::set);

    typedef string (Process::Arguments::*argumentsStringFunction)(int);

    class_<Process::Arguments>("Arguments").
        def("__getitem__", argumentsStringFunction(&Process::Arguments::operator ()));

    class_<Process>("Process").
        add_static_property("environment", &Process::environment);

    typedef void (SymmetryOperation::*intFunction)(int);
    typedef void (SymmetryOperation::*doubleFunction)(double);

    class_<SymmetryOperation>("SymmetryOperation").
        def(init<const SymmetryOperation& >()).
        def("trace", &SymmetryOperation::trace).
        def("zero", &SymmetryOperation::zero).
        def("operate", &SymmetryOperation::operate).
        def("transform", &SymmetryOperation::transform).
        def("unit", &SymmetryOperation::unit).
        def("E", &SymmetryOperation::E).
        def("i", &SymmetryOperation::i).
        def("sigma_h", &SymmetryOperation::sigma_h).
        def("sigma_xz", &SymmetryOperation::sigma_xz).
//        def("sigma_yz", &SymmetryOperation::sigma_yz).
        def("rotateN", intFunction(&SymmetryOperation::rotation)).
        def("rotateTheta", doubleFunction(&SymmetryOperation::rotation)).
        def("c2_x", &SymmetryOperation::c2_x).
        def("c2_y", &SymmetryOperation::c2_y).
        def("transpose", &SymmetryOperation::transpose);

    class_<PointGroup, shared_ptr<PointGroup> >("PointGroup").
        def(init<const char*>()).
        def("symbol", &PointGroup::symbol).
        //def("origin", &PointGroup::origin).
        def("setSymbol", &PointGroup::set_symbol);

    class_<Molecule, shared_ptr<Molecule> >("Molecule").
        def("initWithCheckpoint", &Molecule::init_with_chkpt).
        def("saveToCheckpoint", &Molecule::save_to_chkpt).
        def("initWithIO", &Molecule::init_with_psio).
        def("addAtom", &Molecule::add_atom).
        def("natom", &Molecule::natom).
        def("Z", &Molecule::Z).
        def("x", &Molecule::x).
        def("y", &Molecule::y).
        def("z", &Molecule::z).
        //def("xyz", &Molecule::xyz).
        def("centerOfMass", &Molecule::center_of_mass).
        def("translate", &Molecule::translate).
        def("moveToCOM", &Molecule::move_to_com).
        def("mass", &Molecule::mass).
        def("label", &Molecule::label).
        def("charge", &Molecule::charge).
        def("atomAtPosition", &Molecule::atom_at_position1).
        def("printToOutput", &Molecule::print).
        def("nuclearRepulsionEnergy", &Molecule::nuclear_repulsion_energy).
        def("reorient", &Molecule::reorient).
        def("findPointGroup", &Molecule::find_point_group).
        def("setPointGroup", &Molecule::set_point_group).
        def("formSymmetryInformation", &Molecule::form_symmetry_information);
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
    if (!Py_IsInitialized()) {
        s = strdup("psi");
        Py_Initialize();
        #if PY_VERSION_HEX >= 0x03000000
        Py_SetProgramName(L"psi");
        #else
        Py_SetProgramName(s);
        #endif

        // Track down the location of PSI4's python script directory.
        std::string psiDataDirName = Process::environment("PSIDATADIR") + "/python";

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
        // str is a Boost Python C++ wrapper for Python strings.
        str strStartScript(file.str().c_str());
        //printf(file.str().c_str());
        try {
            PyImport_AppendInittab(strdup("PsiMod"), initPsiMod);
            object objectMain(handle<>(borrowed(PyImport_AddModule("__main__"))));
            object objectDict = objectMain.attr("__dict__");
            s = strdup("import PsiMod");
            PyRun_SimpleString(s);

            object objectScriptInit = exec( strStartScript, objectDict, objectDict );
        }
        catch (error_already_set const& e)
        {
            PyErr_Print();
        }
    } else {
            fprintf(stderr, "Unable to run Python input file.\n");
            return;
    }
}
