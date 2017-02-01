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

#include "liboptions.h"
#include "liboptions_python.h"

namespace psi {

py::list fill_list(py::list l, Data d)
{
    if(d.is_array()){
        // Recurse
        py::list row;
        for(int i = 0; i < d.size(); ++i){
            fill_list(row, d[i]);
        }
        l.append(row);
    }else if(d.type() == "double"){
        l.append(py::float_(d.to_double()));
    }else if(d.type() == "string"){
        l.append(py::str(d.to_string()));
    }else if(d.type() == "boolean"){
        l.append(py::bool_(d.to_integer()));
    }else if(d.type() == "int"){
        l.append(py::int_(d.to_integer()));
    }else{
        throw PSIEXCEPTION("Unknown data type in fill_list");
    }
    return l;
}

py::list Data::to_list() const
{
    return ptr_->to_list();
}

py::list ArrayType::to_list() const
{
    py::list l;
    for(int i = 0; i < array_.size(); ++i)
        fill_list(l, array_[i]);
    return l;
}

PythonDataType::PythonDataType()
{
}

PythonDataType::PythonDataType(const py::object &p)
    : python_object_(p)
{
}

PythonDataType::~PythonDataType()
{ }

std::string PythonDataType::type() const{
    return std::string("python");
}

void PythonDataType::assign(const py::object &p)
{
    python_object_ = p;
    changed();
}

const py::object& PythonDataType::to_python() const
{
    return python_object_;
}

}
