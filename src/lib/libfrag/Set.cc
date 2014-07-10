/*
 * Set.cc
 *
 *  Created on: May 13, 2014
 *      Author: richard
 */
#include <algorithm>
#include "Set.h"
namespace LibFrag{

typedef std::vector<int> SVec;
typedef SVec::iterator VecIt;
typedef VecIt (*Operator)(VecIt,VecIt,VecIt,VecIt,VecIt);

void SetCombine(SVec& ThisSet,const SVec& OtherSet,Operator opr){
	SVec temp=ThisSet,temp2=OtherSet;
	int MaxSize=ThisSet.size()+OtherSet.size();
	ThisSet.resize(MaxSize);//Make sure there is enough room in ThisSet
	VecIt it=opr(temp.begin(),temp.end(),temp2.begin(),temp2.end(),ThisSet.begin());
	ThisSet.resize(it-ThisSet.begin());//Shrink ThisSet down
}

void Set::Union(const Set&other){SetCombine(this->Atoms,other.Atoms,&std::set_union);}
void Set::Intersection(const Set& other){SetCombine(this->Atoms,other.Atoms,&std::set_intersection);}


}


