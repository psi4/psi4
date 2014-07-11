/*
 * libfrag.cc
 *
 *  Created on: May 13, 2014
 *      Author: richard
 */

#include "libfrag.h"
#include "psi4-dec.h" //Where we get the molecule from
#include "MBEFrag.h"
#include <sstream>
#include <iomanip>


#include "Fragmenter.h"
#include "../libparallel/parallel.h"
#include "../../bin/psi4/script.h"
#include "../libmints/molecule.h"
#include "LibFragPrinter.h"
#include "InputManip.h"
#include "MBE.h"
#include "BSSEer.h"

typedef psi::Molecule Mol;
typedef boost::shared_ptr<Mol> SharedMol;
typedef std::vector<int> SVec;
typedef LibFrag::MBEFrag* pSet;
typedef std::vector<pSet> FragSet;
typedef std::string str;

namespace LibFrag{

//Handy template for switching python things to c++ things
template <typename T,typename S>
void ToC(T& cvalue,S& pyvalue){
	cvalue=boost::python::extract<T>(pyvalue);
}

boost::python::list baseGetCall(const int NMer,const int N,std::vector<NMerSet>& Systems,
      bool IsGhost){
   boost::python::list DaList;
   SharedFrag DaNMer=(Systems[NMer])[N];
   if(!IsGhost){
      for(int i=0;i<DaNMer->size();i++)DaList.append((*DaNMer)[i]);
   }
   else{
      for(int i=0;i<DaNMer->Ghosts.size();i++)
         DaList.append(DaNMer->Ghosts[i]);
   }
   return DaList;
}

boost::python::list LibFragHelper::GetNMerN(const int NMer,const int N){
   return baseGetCall(NMer,N,Systems,false);
}

boost::python::list LibFragHelper::GetGhostNMerN(const int NMer,const int N){
   return baseGetCall(NMer,N,Systems,true);
}

void LibFragHelper::Fragment_Helper(boost::python::str& BSSEMethod,
      boost::python::str& FragMethod){
	str fname,bname;
	ToC(bname,BSSEMethod);
	ToC(fname,FragMethod);
	DaOptions.SetBMethod(bname);
	DaOptions.SetFMethod(fname);
	SharedMol AMol=psi::Process::environment.molecule();
	Systems.push_back(NMerSet());

	///Fragment the system
	boost::shared_ptr<Fragmenter> FragFactory;
	if(DaOptions.FMethod==USER_DEFINED)
		FragFactory=boost::shared_ptr<Fragmenter>(new UDFragmenter);
	FragFactory->Fragment(AMol,Systems[0]);

	///Add in any BSSE corrections we may need
	if(DaOptions.BMethod==FULL)
	    BSSEFactory=boost::shared_ptr<BSSEer>(new FullBSSE(AMol->natom()));

	if(DaOptions.BMethod!=NO_BSSE)BSSEFactory->AddBSSEJobs(Systems[0]);

	//At this point need to figure out if the fragments intersect,
	//assuming they don't for now
	Expansion=boost::shared_ptr<GMBE>(new MBE);
}

void LibFragHelper::Embed_Helper(boost::python::str& EmbedMethod){
	std::string name;
	ToC(name,EmbedMethod);
}

void LibFragHelper::Cap_Helper(boost::python::str& CapMethod){
	std::string name;
	ToC(name,CapMethod);
}

double LibFragHelper::CalcEnergy(boost::python::list& Energies){
	std::vector<double*> energies;
	for (int i=0;i<len(Energies);i++){
		int size=len(boost::python::extract<boost::python::list>(Energies[i]));
		energies.push_back(new double[size]);
		for(int j=0;j<size;j++)
		energies[i][j]=boost::python::extract<double>(Energies[i][j]);
	}
	return Expansion->Energy(Systems,energies);
}

void LibFragHelper::NMer_Helper(const int N){
	DaOptions.MBEOrder=N;
    Expansion->SetN(N);
	Expansion->MakeIntersections(Systems);
	Systems.push_back(NMerSet());
	Expansion->MakeNmers(Systems[0],Systems[Systems.size()-1]);
	if(DaOptions.BMethod!=NO_BSSE){
	   for(int order=1;order<Systems.size();order++)
	   BSSEFactory->AddBSSEJobs(Systems[order]);
	}
}

}

