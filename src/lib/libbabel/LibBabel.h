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
#ifndef LIBBABEL_H_
#define LIBBABEL_H_

#include <string>
#include <boost/python.hpp>
namespace psi{
namespace LibBabel{
/** \brief Creates a molecule from a file using OpenBabel
 *
 *  OpenBabel is capable of reading many types of file formats
 *  and this function is designed to interface with OpenBabel.
 *  It is important that the file has the correct extension,
 *  or else OpenBabel does not know how to open it.  If the file
 *  contains crystal data, this function also allows you to
 *  periodically replicate the unit cell.  To do this you need
 *  to use the dim argument.  dim=0, means do not replicate the
 *  unit cell.  dim=1, means that it will be replicated once
 *  creating a shell 1 cell thick (a total of 27 cells), for dim=2
 *  you would have 3125 cells.  In general you have m^m cells, where
 *  m=2*dim+1.
 *
 *
 *
 *  \param[in] FileName The file from which the molecule is read
 *  \param[in] dim Number of periodic replications
 */
boost::python::str ParseFile(const std::string& FileName,const int dim=0);

}}//End namespaces



#endif /* LIBBABEL_H_ */
