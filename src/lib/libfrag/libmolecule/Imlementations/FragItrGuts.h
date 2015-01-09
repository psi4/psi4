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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_IMLEMENTATIONS_FRAGITRGUTS_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_IMLEMENTATIONS_FRAGITRGUTS_H_
#include <set>
#include <boost/shared_ptr.hpp>
#include "MolItrGuts.h"

namespace psi{
namespace LibMolecule{

class FragItrGuts: public MolItrGuts{
   private:
      const FragItrGuts* UpCast(const IteratorGuts& other)const;
   protected:
      typedef std::set<int>::iterator Mem_t;
      const Fragment* Frag_;
      Mem_t MemItr_;
      Mem_t MemItrEnd_;
      MolItrGuts MolItrBegin_;

   public:
      FragItrGuts(bool AtEnd,const Fragment* Frag);
      bool IsEqual(const IteratorGuts& other)const;
      /** \brief Function that defines what it means for two fragments to be
       *         less.
       *
       *  We agree to iterate through the atoms in members first, then the
       *  additional atoms.  We assume that other is a valid iterator to
       *  the same fragment, and that it's MemItr has not been iterated past
       *  MemItrEnd_.  Then our iterator is guaranteed less than other if
       *  its MemItr_ is less than other's, conversely our iterator is
       *  guaranteed greater than other's if its MemItr_ is greater.  If
       *  the two MemItr_'s are equal, then it depends on whether the base
       *  class's are equal, less, or greater (i.e. depends on the iterators
       *  to the two Atoms_ arrays)
       */
      bool IsLess(const IteratorGuts& other)const;
      boost::shared_ptr<const Atom> Atom()const;
      void Next();
      void Previous();
};

}}//End namesapces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_IMLEMENTATIONS_FRAGITRGUTS_H_ */
