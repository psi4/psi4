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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_ITERATORGUTS_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_ITERATORGUTS_H_
#include <vector>
namespace psi{
namespace LibMolecule{
class Atom;
class IteratorGuts{
   public:
      ///Returns a pointer to a copy of the current object
      virtual boost::shared_ptr<IteratorGuts> Clone()const=0;
      IteratorGuts(){}
      IteratorGuts(const IteratorGuts& other){}
      virtual ~IteratorGuts(){}
      virtual bool IsEqual(const IteratorGuts& other)const=0;
      virtual boost::shared_ptr<const Atom> GetAtom()const=0;
      virtual void Next()=0;
      virtual void Previous()=0;
};

}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_ITERATORGUTS_H_ */
