// Include Python's header file.
// Use "python-config --includes" to determine this location.
#include <cstdio>
#include <boost/python.hpp>
//#include <Python.h>
#include <psiconfig.h>
#include <sstream>
#include "script.h"
#include "libchkpt/chkpt.hpp"
#include "liboptions/liboptions.h"
#include "libmints/molecule.h"

#include <libpsio/psio.hpp>

using namespace psi;
using namespace boost;
using namespace boost::python;

namespace psi {
    extern void psiclean(void);
    extern FILE *outfile;
    extern Options options;
}

bool py_psi_input()
{
    // Need to modify input to not take argc and argv
    // And Options object to be global.
    // input::input3();
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

/* Sets the content of the options keyword GEOMETRY
 *
 * This assumes the user is inputing geometry with the following list format
 * geometry([
 *    [ "H", 0.0, 0.0, 0.0],
 *    [ "H", 0.0, 0.0, 1.0]
 * ])
 */
void py_psi_geometry(list geometry)
{
    // Check to see that the geometry keyword exists.
    if (!options.exists("GEOMETRY")) {
        printf("GEOMETRY keyword does not exist in options.\nCreating it.\n");
        options.add_array("GEOMETRY");
    }

    // Rest the contents of the GEOEMTRY array
    options["GEOMETRY"].reset();
    
    size_t length = len(geometry);
    printf("length of the geometry: %d\n", length);
    for (size_t i=0; i<length; ++i) {
        object atom = geometry[i];

        size_t atom_length = len(atom);
        printf("Atom %d is of length %d\n", i, atom_length);

        if (atom_length != 4) {
            printf("Atom %d is of wrong length (expected 4, given %d)\n", i, atom_length);
        }
    }
}

BOOST_PYTHON_MODULE(psi)
{
    def("version", py_psi_version);
    def("clean", py_psi_clean);
    def("configure_io", py_psi_configure_psio);
    def("geometry", py_psi_geometry);
    
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
        def("xyz", &Molecule::xyz).
        def("centerOfMass", &Molecule::center_of_mass).
        def("translate", &Molecule::translate).
        def("moveToCOM", &Molecule::move_to_com).
        def("mass", &Molecule::mass).
        def("label", &Molecule::label).
        def("charge", &Molecule::charge).
        def("atomAtPosition", &Molecule::atom_at_position).
        def("printToOutput", &Molecule::print).
        def("nuclearRepulsionEnergy", &Molecule::nuclear_repulsion_energy).
        def("reorient", &Molecule::reorient);

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
//        PyImport_AppendInittab(s, initpsi);
        Py_Initialize();
        #if PY_VERSION_HEX >= 0x03000000
        Py_SetProgramName(L"psi");
        #else
        Py_SetProgramName(s);
        #endif
    }
    if (Py_IsInitialized()) {
        char line[256];
        std::stringstream file;
        while(fgets(line, sizeof(line), input)) {
            file << line;
        }
//        printf("Input file:\n%s", file.str().c_str());
        str strStartScript(file.str().c_str());

        try {
            PyImport_AppendInittab(s, initpsi);
            object objectMain(handle<>(borrowed(PyImport_AddModule("__main__"))));
            object objectDict = objectMain.attr("__dict__");
            s = strdup("import psi; from psi import *;");
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
