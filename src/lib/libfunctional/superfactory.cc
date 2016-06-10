/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/shared_ptr.hpp>
#include <liboptions/liboptions.h>
#include <liboptions/liboptions_python.h>
#include <psi4-dec.h>
#include "superfunctional.h"

#define PY_TRY(ptr, command)  \
     if(!(ptr = command)){    \
         PyErr_Print();       \
         exit(1);             \
     }

using namespace boost::python;

namespace psi {

boost::shared_ptr<SuperFunctional> SuperFunctional::current(Options& options, int npoints, int deriv)
{
    if (npoints == -1) {
        npoints = options.get_int("DFT_BLOCK_MAX_POINTS");
    }

    boost::shared_ptr<SuperFunctional> super;
    if (options.get_str("DFT_FUNCTIONAL") == "GEN" || options.get_str("DFT_FUNCTIONAL") == "") {
        boost::python::object pySuper = dynamic_cast<PythonDataType*>(options["DFT_CUSTOM_FUNCTIONAL"].get())->to_python();
        super = boost::python::extract<boost::shared_ptr<SuperFunctional> >(pySuper);
        if (!super) {
            throw PSIEXCEPTION("Custom Functional requested, but nothing provided in DFT_CUSTOM_FUNCTIONAL");
        }
    } else {
        super = SuperFunctional::build(options.get_str("DFT_FUNCTIONAL"), npoints, deriv);
        if (options["DFT_OMEGA"].has_changed() && super->is_x_lrc()) {
            super->set_x_omega(options.get_double("DFT_OMEGA"));
        }
        if (options["DFT_ALPHA"].has_changed()) {
            super->set_x_alpha(options.get_double("DFT_ALPHA"));
        }
        if (options["DFT_OMEGA_C"].has_changed() && super->is_c_lrc()) {
            super->set_c_omega(options.get_double("DFT_OMEGA_C"));
        }
        if (options["DFT_ALPHA_C"].has_changed()) {
            super->set_c_alpha(options.get_double("DFT_ALPHA_C"));
        }
    }

    if (npoints != super->max_points())
        super->set_max_points(npoints);
    if (deriv != super->deriv())
        super->set_deriv(deriv);

    return super;
}

boost::shared_ptr<SuperFunctional> SuperFunctional::build(const std::string& alias, int max_points, int deriv)
{
    boost::shared_ptr<SuperFunctional> super;

    if (Py_IsInitialized()) {
        try {
            // Grab the SuperFunctional off of the Python plane
            PyObject *functional;
            PY_TRY(functional, PyImport_ImportModule("procedures.functional") );
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
