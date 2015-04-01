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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_MOLITRGUTS_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_MOLITRGUTS_H_
#include <vector>
#include <boost/shared_ptr.hpp>
#include "IteratorGuts.h"
namespace psi{
namespace LibMolecule{
class Atom;
class MolItrGuts:public IteratorGuts{
   private:
      const MolItrGuts* UpCast(const IteratorGuts& other)const;
      void Copy(const MolItrGuts& other);
   protected:
      ///Typedef to make the next typedef fit on one line
      typedef boost::shared_ptr<const Atom> SharedAtom;
      ///Typedef b/c this type is horrendously long...
      typedef std::vector<SharedAtom>::const_iterator Itr_t;
      Itr_t MyItr_;
   public:
      virtual boost::shared_ptr<IteratorGuts> Clone()const;
      const MolItrGuts& operator=(const MolItrGuts& other);
      MolItrGuts(const MolItrGuts& other);
      MolItrGuts(){}
      MolItrGuts(const Itr_t& MyItr);
      virtual ~MolItrGuts(){}
      virtual bool IsEqual(const IteratorGuts& other)const;
      boost::shared_ptr<const Atom> GetAtom()const{return *MyItr_;}
      virtual void Next(){++MyItr_;}
      virtual void Previous(){--MyItr_;}
};

}}//End namespaces

#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_MOLITRGUTS_H_ */
