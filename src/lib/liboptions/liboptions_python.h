/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef _psi_src_lib_liboptions_python_h
#define _psi_src_lib_liboptions_python_h

#include <boost/python.hpp>
#include <boost/python/object.hpp>

namespace psi {

class PythonDataType : public DataType
{
    boost::python::object python_object_;
public:
    PythonDataType();
    PythonDataType(const boost::python::object& p);
    virtual ~PythonDataType();

    virtual std::string type() const;

    const boost::python::object& to_python() const;

    void assign(const boost::python::object& p);
};

}

#endif // _psi_src_lib_liboptions_python_h
