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

#include "LibFragFragment.h"
#include "Implementations/FragItrGuts.h"
namespace psi{
namespace LibMolecule{

typedef std::vector<int> SVec;
typedef SVec::iterator VecIt;
typedef boost::shared_ptr<FragItrGuts> SharedItr;

void Fragment::Copy(const Fragment& other){
   this->Mol_=other.Mol_;
   this->Members_=other.Members_;
   this->OtherMembers_=other.OtherMembers_;
   this->SN_=other.SN_;
   MolItr AtomI=other.Molecule::Begin(),AtomEnd=other.Molecule::End();
   for(;AtomI!=AtomEnd;++AtomI)
      (*this)<<(*(*AtomI));
}

Fragment::Fragment(boost::shared_ptr<const Molecule> Mol):Mol_(Mol){}

Fragment::Fragment(boost::shared_ptr<const Molecule> Mol,const long int ID):
      Mol_(Mol){
   SN_.insert(ID);
}

MolItr Fragment::Begin()const{
   return MolItr(SharedItr(new FragItrGuts(false,this)));
}

MolItr Fragment::End()const{
   return MolItr(SharedItr(new FragItrGuts(true,this)));
}

Fragment& Fragment::operator<<(const Atom& NewAtom){
   Molecule::operator<<(NewAtom);
   return *this;
}

void Fragment::AddRepAtom(const Atom& NewAtom, const int i){
   (*this)<<NewAtom;
   OtherMembers_.push_back(i);
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

template<typename T>
void SetUnion(std::set<T>& lhs,const std::set<T>&rhs){
   typedef typename std::set<T>::iterator It_t;
   std::pair<It_t,bool> ret;
   //For the first iteration we utilize the
   It_t Atoms=rhs.begin();
   ret=lhs.insert((*Atoms));
   ++Atoms;
   for(;Atoms!=rhs.end();++Atoms)
      lhs.insert(ret.first,(*Atoms));
}

const Fragment& Fragment::operator+=(const Fragment& other){
      if(other.Members_.empty())return *this;
      SetUnion(this->Members_,other.Members_);
      SetUnion(this->SN_,other.SN_);
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

bool Fragment::operator==(const Fragment& other)const{
   //Have to be the same size to be equal
   if(this->Members_.size()!=other.Members_.size())return false;
   Fragment temp=(*this)-other;
   return(temp.Members_.size()==this->Members_.size());
}

bool Fragment::operator<(const Fragment& other)const{
   //In order for this to be a proper subset of other, it's size has
   //to be smaller
   if(this->Members_.size()>=other.Members_.size())return false;
   Fragment temp=(*this)-other;
   return (temp.Members_.size()==this->Members_.size());
}
}}

