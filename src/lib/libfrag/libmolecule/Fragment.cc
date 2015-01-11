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

#include "Fragment.h"

namespace psi{
namespace LibMolecule{

typedef std::vector<int> SVec;
typedef SVec::iterator VecIt;

Fragment::Fragment(boost::shared_ptr<const Molecule> Mol):Mol_(Mol){}

MolItr Fragment::Begin()const{
   return MolItr(new FragItrGuts(false,this));
}

MolItr Fragment::End()const{
   return MolItr(new FragItrGuts(true,this));
}

Fragment& Fragment::operator<<(boost::shared_ptr<const Atom> NewAtom){
   Molecule::operator<<(NewAtom);
   OtherMembers_.push_back(-1);
   return *this;
}

void Fragment::AddRepAtom(boost::shared_ptr<const Atom> NewAtom, const int i){
   (*this)<<NewAtom;
   OtherMembers_.back()=i;
}


int Fragment::NAtoms()const{
   int TotalSize=Molecule::NAtoms();
   TotalSize+=Members_.size();
   return TotalSize;
}

boost::shared_ptr<const Atom> Fragment::LookUp(const int i)const{
   int size=Members_.size();
   if(i<size){
      std::set<int>::iterator It=Members_.begin();
      for(int j=0;j<i;j++)++It;
      return (*Mol_)[(*It)];
   }
   else return Atoms_[i-size];
   //Can't call Molecule::operator[], as that will call this fxn again...
   //return (i<size?(*Mol_)[Members_[i]]:Atoms_[i-size]);
}


/********** Stuff below here is for the operators*************/

const Fragment& Fragment::operator+=(const Fragment& other){
      if(other.Members_.empty())return *this;
      typedef std::set<int>::iterator It_t;
      std::pair<It_t,bool> ret;
      //For the first iteration we utilize the
      It_t Atoms=other.Members_.begin();
      ret=this->Members_.insert((*Atoms));
      ++Atoms;
      for(;Atoms!=other.Members_.end();++Atoms)
         this->Members_.insert(ret.first,(*Atoms));

   return *this;
}

Fragment Fragment::operator+(const Fragment& other)const{
   Fragment temp(*this);
   temp+=other;
   return temp;
}

const Fragment& Fragment::operator-=(const Fragment& other){
   if(other.Members_.empty())return *this;
   std::set<int> Intersec;
   std::set_intersection(this->Members_.begin(),this->Members_.end(),
         other.Members_.begin(),other.Members_.end(),
         std::inserter(Intersec,Intersec.begin()));
   this->Members_=Intersec;
   return *this;
}

Fragment Fragment::operator-(const Fragment&other)const{
   Fragment temp(*this);
   temp-=other;
   return temp;
}


}}

