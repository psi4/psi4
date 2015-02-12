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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_FRAGITRGUTS_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_FRAGITRGUTS_H_
#include <set>
#include <boost/shared_ptr.hpp>

#include "../Implementations/MolItrGuts.h"

namespace psi{
namespace LibMolecule{
class Fragment;
class FragItrGuts: public MolItrGuts{
   private:
      void Copy(const FragItrGuts& other);
      const FragItrGuts* UpCast(const IteratorGuts& other)const;
   protected:
      const Fragment* Frag_;
      std::set<int>::iterator MemItr_;
      std::set<int>::iterator MemItrEnd_;
      MolItrGuts MolItrBegin_;
   public:
      boost::shared_ptr<IteratorGuts> Clone()const;
      FragItrGuts(const FragItrGuts& other);
      const FragItrGuts& operator=(const FragItrGuts& other);
      FragItrGuts(bool AtEnd,const Fragment* Frag);
      bool IsEqual(const IteratorGuts& other)const;
      boost::shared_ptr<const Atom> GetAtom()const;
      void Next();
      void Previous();
};

}}//End namesapces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_FRAGITRGUTS_H_ */
