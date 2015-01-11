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

#include "FragItrGuts.h"
#include "../Fragment.h"
namespace psi{
namespace LibMolecule{

FragItrGuts::FragItrGuts(bool AtEnd,const Fragment* Frag):
      MolItrGuts(AtEnd?Frag->Atoms_.end():Frag->Atoms_.begin()),
      MemItr_(AtEnd?Frag->Members_.end():Frag->Members_.begin()),
      MemItrEnd_(Frag->Members_.end()),
      MolItrBegin_(Frag->Molecule::Begin()),
      Frag_(Frag){
}


bool FragItrGuts::IsEqual(const IteratorGuts& other)const{
   const FragItrGuts* Other=UpCast(other);
   return (this->MemItr_==Other->MemItr_&&MolItrGuts::IsEqual(other));
}

bool FragItrGuts::IsLess(const IteratorGuts& other)const{
   const FragItrGuts* Other=UpCast(other);
   if(!this->MemItr_<Other->MemItr_){//Either equal or greater
      if(this->MemItr_==this->MemItrEnd_)//Both done
         return MolItrGuts::IsLess(other);
      else return false;//Greater or equal on MemItr
   }
   else return true;
}

boost::shared_ptr<const Atom> FragItrGuts::Atom()const{
   return (this->MemItr_<this->MemItrEnd_?
     (*(Frag_->Mol_))[(*MemItr_)]:
      MolItrGuts::Atom()
   );
}

void FragItrGuts::Next(){
   if(MemItr_<MemItrEnd_)++MemItr_;
   else MolItrGuts::Next();
}

void FragItrGuts::Previous(){
   if(MolItrGuts::IsEqual(MolItrBegin_))--MemItr_;
   else MolItrGuts::Previous();
}

}}

