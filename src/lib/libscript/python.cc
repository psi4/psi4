// Include Python's header file.
// Use "python-config --includes" to determine this location.
#include <Python.h>
#include "script.h"

Python::Python() : Script()
{
    
}

Python::~Python()
{
    
}

void Python::initialize()
{
    // Python C-interface call to add PSI4 functionality must
    // be done before the call to Py_Initialize.
    Py_Initialize();
}

void Python::finalize()
{
    Py_Finalize();
}
