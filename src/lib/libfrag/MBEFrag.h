/*
 * MBEFrag.h
 *
 *  Created on: May 20, 2014
 *      Author: richard
 */

#ifndef MBEFRAG_H_
#define MBEFRAG_H_
#include "Set.h"
#include <vector>
namespace LibFrag{


class MBEFrag: public Set{
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
		   for(int i=0;i<MBEOrder;i++)Parents.push_back(Ps[i]);
		}

		MBEFrag(const int MBEOrder_=0,const int *Parents_=NULL):
			MBEOrder(MBEOrder_){if(Parents_!=NULL)SetParents(Parents_);}

		MBEFrag(const MBEFrag& other):Set(other){Copy(other);}

		MBEFrag& operator=(const MBEFrag& other){
		   if(this!=&other){
		      Set::operator=(other);
		      Copy(other);
		   }
		   return *this;
		}
};

}


#endif /* MBEFRAG_H_ */
