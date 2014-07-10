/*
 * Fragmenter.cc
 *
 *  Created on: May 15, 2014
 *      Author: richard
 */

#include "Fragmenter.h"
#include "MBEFrag.h"
#include "../libmints/molecule.h"
#include <boost/shared_ptr.hpp>
#include <iostream>

namespace LibFrag {


void UDFragmenter::Fragment(SharedMol& AMol,NMerSet& Monomers){
	for(int frags=0,index=0;frags<AMol->nfragments();frags++){
		SharedMol Frag=AMol->py_extract_subsets_6(frags+1);
		int AtomsInFrag=Frag->natom();
		int zero=0;
		SharedFrag temp(new LibFrag::MBEFrag(1,&zero));
		Monomers.push_back(temp);
		for(int Atoms=0;Atoms<AtomsInFrag;Atoms++){
			boost::shared_ptr<Set> MonoI=Monomers[frags];
			(*MonoI)<<index++;
		}
	}
}
}

