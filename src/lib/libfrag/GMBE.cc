/*
 * GMBE.cc
 *
 *  Created on: May 15, 2014
 *      Author: richard
 */

#include "GMBE.h"
#include "MBEFrag.h"

namespace LibFrag {

void GMBE::MakeNmers(const NMerSet& Monomers,NMerSet& NMers){
	int NFrags=Monomers.size();
	std::vector<int> indices;
	for(int i=0;i<N;i++)indices.push_back(i);
	int index=0;
	bool done=false;
	while(!done){
		MBEFrag temp=(*Monomers[indices[0]]);
		temp.SetMBEOrder(N);
		temp.SetParents(&indices[0]);
		for(int i=1;i<N;i++)
			temp*=(*Monomers[indices[i]]);
		NMers.push_back(boost::shared_ptr<MBEFrag>(new MBEFrag(temp)));
		bool igood=false;
		for(int i=N-1;i>=0&&!igood;i--){
			int MaxMono=NFrags-(N-i);
			if(indices[i]<MaxMono){ //We can increase index i
				indices[i]++;
				igood=true;
				for(int j=i+1;j<N;j++)
					indices[j]=indices[i]+(j-i); //Reset indices after i
			}
			else if(i==0){
				done=true;
				igood=true;
			}
		}
	}
}

}

