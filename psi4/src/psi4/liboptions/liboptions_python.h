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

#ifndef _psi_src_lib_liboptions_python_h
#define _psi_src_lib_liboptions_python_h

#include <psi4/pybind11.h>

namespace psi {

class PythonDataType : public DataType
{
    pybind11::object python_object_;
public:
    PythonDataType();
    PythonDataType(const pybind11::object& p);
    virtual ~PythonDataType();

    virtual std::string type() const;

    const pybind11::object& to_python() const;
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#endif
    void assign(const pybind11::object& p);
#ifdef __clang__
#pragma clang diagnostic pop
#endif
};

}

#endif // _psi_src_lib_liboptions_python_h
