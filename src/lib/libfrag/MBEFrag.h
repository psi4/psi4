/*
 * MBEFrag.h
 *
 *  Created on: May 20, 2014
 *      Author: richard
 */

#ifndef MBEFRAG_H_
#define MBEFRAG_H_
#include "AtomSet.h"
#include "psi4-dec.h"
#include <vector>
namespace LibFrag{


class MBEFrag: public AtomSet{
	private:
		std::vector<int> Parents;
		int MBEOrder;
		void Copy(const MBEFrag& other){
		   this->Parents=other.Parents;
		   this->MBEOrder=other.MBEOrder;
		}
	public:
		int GetMBEOrder(){return MBEOrder;}
		int ParentI(const int I){return Parents[I];}
		void SetMBEOrder(const int N){MBEOrder=N;}

		void SetParents(const int* Ps){
		   Parents.clear();
		   for(int i=0;i<MBEOrder;i++)Parents.push_back(Ps[i]);
		}

		void PrintParents(){
		   for(int i=0;i<MBEOrder;i++)psi::fprintf(psi::outfile,"%d ",Parents[i]);
		}
		MBEFrag(const int MBEOrder_=0,const int *Parents_=NULL):
			MBEOrder(MBEOrder_){if(Parents_!=NULL)SetParents(Parents_);}

		MBEFrag(const MBEFrag& other):AtomSet(other){Copy(other);}

		MBEFrag& operator=(const MBEFrag& other){
		   if(this!=&other){
		      AtomSet::operator=(other);
		      Copy(other);
		   }
		   return *this;
		}
};

}


#endif /* MBEFRAG_H_ */
