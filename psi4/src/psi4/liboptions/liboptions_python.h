/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_liboptions_python_h
#define _psi_src_lib_liboptions_python_h

#include "liboptions.h"

#include "psi4/pybind11.h"

namespace psi {

class PythonDataType : public DataType
{
    py::object python_object_;
public:
    PythonDataType();
    PythonDataType(const py::object& p);
    virtual ~PythonDataType();

    virtual std::string type() const;

    const py::object& to_python() const;
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#endif
    void assign(const py::object& p);
#ifdef __clang__
#pragma clang diagnostic pop
#endif
};

}

#endif // _psi_src_lib_liboptions_python_h
