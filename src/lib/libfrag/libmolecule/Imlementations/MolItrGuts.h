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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_IMLEMENTATIONS_MOLITRGUTS_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_IMLEMENTATIONS_MOLITRGUTS_H_
#include <vector>
#include <boost/shared_ptr.hpp>

namespace psi{
namespace LibMolecule{
class Atom;
class MolItrGuts:public IteratorGuts{
   private:
      const MolItrGuts* UpCast(const IteratorGuts& other)const{
         return dynamic_cast<const MolItrGuts*>(& other);
      }
   protected:
      typedef std::vector<boost::shared_ptr<const Atom> >::iterator Itr_t;
      Itr_t MyItr_;
   public:
      MolItrGuts(Itr_t& MyItr):MyItr_(MyItr){}
      virtual ~MolItrGuts(){}
      virtual bool IsEqual(const IteratorGuts& other)const{
         return this->MyItr_==UpCast(other)->MyItr_;
      }

      virtual bool IsLess(const IteratorGuts& other)const{
         return this->MyItr_<UpCast(other)->MyItr_;
      }
      boost::shared_ptr<const Atom> Atom()const{
         return *MyItr_;
      }
      virtual void Next(){
         ++MyItr_;
      }
      virtual void Previous(){
         --MyItr_;
      }
};

}}//End namespaces

#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_IMLEMENTATIONS_MOLITRGUTS_H_ */
