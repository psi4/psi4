/*
 * Set.cc
 *
 *  Created on: May 13, 2014
 *      Author: richard
 */
#include <algorithm>
#include "Set.h"
#include "psi4-dec.h"

namespace LibFrag {

typedef std::vector<int> SVec;
typedef SVec::iterator VecIt;
typedef VecIt (*Operator)(VecIt, VecIt, VecIt, VecIt, VecIt);

void SetCombine(SVec& ThisSet, const SVec& OtherSet, Operator opr) {
   SVec temp=ThisSet,temp2=OtherSet;
   int MaxSize=ThisSet.size()+OtherSet.size();
   ThisSet.resize(MaxSize); //Make sure there is enough room in ThisSet
   VecIt it=opr(temp.begin(), temp.end(), temp2.begin(), temp2.end(),
         ThisSet.begin());
   ThisSet.resize(it-ThisSet.begin()); //Shrink ThisSet down
}

void Set::Sort() {
   std::sort(Atoms.begin(), Atoms.end());
}
void Set::Union(const Set&other) {
   SetCombine(this->Atoms, other.Atoms, &std::set_union);
}
void Set::Intersection(const Set& other) {
   SetCombine(this->Atoms, other.Atoms, &std::set_intersection);
}

bool Set::contains(const int i){
   SVec temp(1,i);
   SetCombine(temp,this->Atoms,&std::set_intersection);
   return (temp.size()==1);
}

bool Set::IsSuperset(const Set& other)const{
   bool is=false;
   if(other.size()<this->size()){
      Set temp=other/(*this);
      if(temp.size()==other.size())is=true;
   }
   return is;
}

bool Set::IsSubset(const Set& other)const{
   bool is=false;
   if(other.IsSuperset(*this))is=true;
   return is;
}

bool Set::IsEqual(const Set& other)const{
   bool is=false;
   if(other.size()==this->size()){
      Set temp=other/(*this);
      if(temp.size()==other.size())is=true;
   }
   return is;
}

void Set::print_out()const{
   for(int i=0;i<size();i++)psi::fprintf(psi::outfile,"%d ",Atoms[i]);
   psi::fprintf(psi::outfile,"\n");
}

}//End namespace LibFrag

