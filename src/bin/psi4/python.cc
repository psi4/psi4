// Include Python's header file.
// Use "python-config --includes" to determine this location.
#include <Python.h>
#include <psiconfig.h>
#include "script.h"

using namespace psi;
using namespace boost;

extern void psiclean();

static PyObject *
py_psi_version(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", PSI_VERSION);
}

static PyObject *
py_psi_input(PyObject *self, PyObject *args)
{
//    int result = input::input(options, 0, NULL);
    return Py_BuildValue("i", 0);
}

static PyObject *
py_psi_clean(PyObject *self, PyObject *args)
{
//    psiclean();
    return Py_BuildValue("i", 0);
}

static PyMethodDef PsiMethods[] = {
    { "version", py_psi_version, METH_NOARGS, "Obtain version information from PSI."},
    { "input",   py_psi_input,   METH_NOARGS, "Run the input module."},
    { "clean",   py_psi_clean,   METH_NOARGS, "Run the psiclean module."},
    { NULL, NULL, 0, NULL } /* Sentinel */
};

static struct PyModuleDef PsiModule = {
    PyModuleDef_HEAD_INIT,
    "psi",   /* name of the module */
    NULL,    /* something with documentation */
    -1,      /* size of the per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
    PsiMethods
};

PyMODINIT_FUNC
PyInit_psi()
{
    return PyModule_Create(&PsiModule);
}

Python::Python() : Script()
{
    
}

Python::~Python()
{
    
}

void Python::initialize()
{
    PyImport_AppendInittab("psi", PyInit_psi);
    Py_SetProgramName(L"psi");    
    Py_Initialize();
    PyImport_ImportModule("psi");
}

void Python::finalize()
{
    Py_Finalize();
}

void Python::run(FILE *input)
{
    if (input == NULL)
        return;
    if (!Py_IsInitialized()) {
        initialize();
    }
    if (Py_IsInitialized()) {
        PyRun_AnyFile(input, "User Input File");
    } else {
        fprintf(stderr, "Unable to run Python input file.\n");
        return;
    }
}