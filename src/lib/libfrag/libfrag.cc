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


typedef psi::Molecule Mol;
typedef boost::shared_ptr<Mol> SharedMol;
typedef std::vector<int> SVec;
typedef LibFrag::MBEFrag* pSet;
typedef std::vector<pSet> FragSet;
typedef std::string str;

static bool ShowOut=false;
namespace LibFrag{

template <typename T,typename S>
void ToC(T& cvalue,S& pyvalue){
	cvalue=boost::python::extract<T>(pyvalue);
}

boost::python::list LibFragHelper::GetNMerN(const int NMer,const int N){
	boost::python::list DaList;
	SharedFrag DaNMer=(Systems[NMer])[N];
	for(int i=0;i<DaNMer->size();i++)DaList.append((*DaNMer)[i]);
	return DaList;
}



void LibFragHelper::Fragment_Helper(boost::python::str& FragMethod){
	std::string name;
	ToC(name,FragMethod);
	SharedMol AMol=psi::Process::environment.molecule();
	Systems.push_back(NMerSet());
	boost::shared_ptr<Fragmenter> FragFactory;
	if(name=="user_defined"){
	    DaOptions.FMethod=USER_DEFINED;
		FragFactory=boost::shared_ptr<Fragmenter>(new UDFragmenter);
	}
	FragFactory->Fragment(AMol,Systems[0]);
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
	Expansion->SetN(N);
	Expansion->MakeIntersections(Systems);
	Systems.push_back(NMerSet());
	Expansion->MakeNmers(Systems[0],Systems[Systems.size()-1]);
}

}

