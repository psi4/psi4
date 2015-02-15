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
#include "LibFragFragment.h"

namespace psi{
namespace LibMolecule{

typedef boost::shared_ptr<IteratorGuts> SharedGuts;
SharedGuts FragItrGuts::Clone()const{
   SharedGuts temp(new FragItrGuts(*this));
   return temp;
}

FragItrGuts::FragItrGuts(const FragItrGuts& other):MolItrGuts(other){
   this->Copy(other);
}

const FragItrGuts& FragItrGuts::operator=(const FragItrGuts& other){
   if(this!=&other)this->Copy(other);
   return *this;
}

void FragItrGuts::Copy(const FragItrGuts& other){
   this->Frag_=other.Frag_;
   this->MemItr_=other.MemItr_;
   this->MemItrEnd_=other.MemItrEnd_;
   this->MolItrBegin_=other.MolItrBegin_;
}

const FragItrGuts* FragItrGuts::UpCast(const IteratorGuts& other)const{
   return dynamic_cast<const FragItrGuts*>(&other);
}

FragItrGuts::FragItrGuts(bool AtEnd,const Fragment* Frag):
      MolItrGuts(AtEnd?Frag->Atoms_.end():Frag->Atoms_.begin()),
      MemItr_(AtEnd?Frag->Members_.end():Frag->Members_.begin()),
      MemItrEnd_(Frag->Members_.end()),
      MolItrBegin_(Frag->Atoms_.begin()),
      Frag_(Frag){
}


bool FragItrGuts::IsEqual(const IteratorGuts& other)const{
   const FragItrGuts* Other=UpCast(other);
   return (this->MemItr_==Other->MemItr_&&MolItrGuts::IsEqual(other));
}

boost::shared_ptr<const Atom> FragItrGuts::GetAtom()const{
   return (this->MemItr_!=this->MemItrEnd_?
     (*(Frag_->Mol_))[(*MemItr_)]:
      MolItrGuts::GetAtom()
   );
}

void FragItrGuts::Next(){
   if(MemItr_!=MemItrEnd_)++MemItr_;
   else MolItrGuts::Next();
}

void FragItrGuts::Previous(){
   if(MolItrGuts::IsEqual(MolItrBegin_))--MemItr_;
   else MolItrGuts::Previous();
}

}}

