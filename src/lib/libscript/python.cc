// Include Python's header file.
// Use "python-config --includes" to determine this location.
#include <Python.h>
#include "script.h"

using namespace psi;
using namespace boost;

static PyObject *
py_psi_version(PyObject *self, PyObject *args)
{
    return Py_BuildValue("i", 400);
}

static PyMethodDef PsiMethods[] = {
    { "version", py_psi_version, 0,
      "Obtain version information from PSI."},
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
