#include <pyconfig.h>
#include <Python.h>

int main()
{
    PyObject * args = NULL; 
    PyObject * result = NULL;
    // Initiliaze the Python interpreter
    Py_Initialize();
    // Run a simple print('Hello, world!')
    PyRun_SimpleString("print('Hello, world!'");
    // Free all temporary Python objects.
    Py_XDECREF(args); 
    Py_XDECREF(result);
    // Finalize the interpreter 
    Py_Finalize();
   
    return 0;
}

