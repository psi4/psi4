#include <boost/shared_ptr.hpp>
#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <psi4-dec.h>
#include "superfunctional.h"
#include "functional.h"

#define PY_TRY(ptr, command)  \
     if(!(ptr = command)){    \
         PyErr_Print();       \
         exit(1);             \
     }

using namespace boost::python;

namespace psi {

boost::shared_ptr<SuperFunctional> SuperFunctional::build(const std::string& alias, int max_points, int deriv)
{
    boost::shared_ptr<SuperFunctional> super;

    if (Py_IsInitialized()) {
        try {
            // Grab the SuperFunctional off of the Python plane
            PyObject *functional;
            PY_TRY(functional, PyImport_ImportModule("functional") );
            PyObject *function;
            PY_TRY(function, PyObject_GetAttrString(functional, "build_superfunctional"));
            PyObject *pargs;
            PY_TRY(pargs, Py_BuildValue("(sii)", alias.c_str(), max_points, deriv));
            PyObject *ret;
            PY_TRY(ret, PyEval_CallObject(function, pargs));

            // Extract the SuperFunctional
            super = boost::python::extract<boost::shared_ptr<SuperFunctional> >(ret);

            // Decref Python env pointers
            Py_DECREF(ret);
            Py_DECREF(pargs);
            Py_DECREF(function);
            Py_DECREF(functional);
        }
        catch (error_already_set const& e)
        {
            PyErr_Print();
            exit(1);
        }
    }
    else {
        throw PSIEXCEPTION("Unable to parse superfunctional.\n");
    }

    return super;
}

}
