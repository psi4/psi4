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
#include "MBE.h"
#include "BSSEer.h"
#include "Capper.h"

typedef std::vector<int> SVec;
typedef LibFrag::MBEFrag* pSet;
typedef std::vector<pSet> FragSet;
typedef std::string str;

namespace LibFrag {

//Handy template for switching python things to c++ things
template <typename T, typename S>
void ToC(T& cvalue, S& pyvalue) {
   cvalue=boost::python::extract<T>(pyvalue);
}

boost::python::list baseGetCall(const int NMer, const int N,
      std::vector<NMerSet>& Systems, bool IsGhost) {
   boost::python::list DaList;
   SharedFrag DaNMer=(Systems[NMer])[N];
   if (!IsGhost) {
      for (int i=0; i<DaNMer->size(); i++)
         DaList.append((*DaNMer)[i]);
   }
   else {
      for (int i=0; i<DaNMer->Ghosts.size(); i++)
         DaList.append(DaNMer->Ghosts[i]);
   }
   return DaList;
}

boost::python::list LibFragHelper::GetNMerN(const int NMer, const int N) {
   return baseGetCall(NMer, N, Systems, false);
}

boost::python::list LibFragHelper::GetGhostNMerN(const int NMer, const int N) {
   return baseGetCall(NMer, N, Systems, true);
}

int LibFragHelper::RunFrags() {
   return Expansion->RunFrags();
}

void LibFragHelper::ReadMOs() {
   if (!this->IsGMBE()) {
      MOFiles.push_back(psi::MOFile());
      MOFiles.back().Init();
   }
}

void LibFragHelper::WriteMOs(const int N, const int x) {
   if (!this->IsGMBE()) {
      SharedFrag nmer=Systems[N][x];
      psi::MOFile File2Write(MOFiles[nmer->ParentI(0)]);
      for (int i=1; i<=N; i++) {
         int frag=nmer->ParentI(i);
         psi::MOFile temp=File2Write.DirectSum(MOFiles[frag]);
         File2Write=temp;
      }
      File2Write.WriteFile();
   }

}

void LibFragHelper::Fragment_Helper(boost::python::str& FragMethod, const int N,
      boost::python::str& EmbedMethod, boost::python::str& CapMethod,
      boost::python::str& BSSEMethod) {
   str fname,bname,ename,cname;
   ToC(bname, BSSEMethod);
   ToC(fname, FragMethod);
   ToC(ename, EmbedMethod);
   ToC(cname, CapMethod);
   DaOptions.SetBMethod(bname);
   DaOptions.SetFMethod(fname);
   DaOptions.SetCMethod(cname);
   DaOptions.SetEMethod(ename);
   DaOptions.MBEOrder=N;
   DaOptions.PrintOptions();
   SharedMol AMol=psi::Process::environment.molecule();
   Systems.push_back(NMerSet());

   ///Fragment the system
   boost::shared_ptr<Fragmenter> FragFactory=DaOptions.MakeFragFactory();
   FragProps FProp=FragFactory->Fragment(AMol, Systems[0]);
   Expansion=(
         FProp.disjoint ? boost::shared_ptr<GMBE>(new MBE(N)) :
               boost::shared_ptr<GMBE>(new GMBE(N)));

   ///If the user wants a one-body GMBE calc, we hack a bit here and add
   ///the intersections into the fragments
   if (Expansion->IsGMBE()&&N==1) {
      Systems.push_back(Systems[0]);
      Expansion->MakeIntersections(Systems);
      Systems[0]=Systems[1];
      Systems[0].insert(Systems[0].end(), Systems[2].begin(), Systems[2].end());
      for (int j=0; j<Systems[0].size(); j++) {
         Systems[0][j]->print_out();
      }
   }
   if (FProp.severed) {
      CapFactory=DaOptions.MakeCapFactory(AMol);
      CapFactory->MakeCaps(Systems[0]);
   }
   ///Add in any BSSE corrections we may need
   BSSEFactory=DaOptions.MakeBSSEFactory(AMol->natom());

   if (DaOptions.BMethod!=NO_BSSE) BSSEFactory->AddBSSEJobs(Systems[0]);
}

void LibFragHelper::Embed_Helper(boost::python::str& EmbedMethod) {

}

std::string LibFragHelper::Cap_Helper(int order, int frag) {
   std::stringstream caps;
   for (int i=0; i<Systems[order][frag]->Caps.size(); i++)
      caps<<Systems[order][frag]->Caps[i].print_out();
   return caps.str();
}

double LibFragHelper::CalcEnergy(boost::python::list& Energies) {
   std::vector<double*> energies;
   for (int i=0; i<len(Energies); i++) {
      int size=len(boost::python::extract<boost::python::list>(Energies[i]));
      energies.push_back(new double[size]);
      for (int j=0; j<size; j++)
         energies[i][j]=boost::python::extract<double>(Energies[i][j]);
   }
   std::string dashes;
   dashes="------------------------------------------------------------------\n";
   for(int i=0;i<Systems.size();i++){
      if(i==0)psi::fprintf(psi::outfile,"Monomer #     Energy (a.u.)\n");
      else psi::fprintf(psi::outfile,  "%2d-mer  #     Energy (a.u.)\n",i+1);
      psi::fprintf(psi::outfile,dashes.c_str());
      for(int j=0;j<Systems[i].size();j++){
         if(i!=0)Systems[i][j]->PrintParents();
         else psi::fprintf(psi::outfile,"%d ",j);
         psi::fprintf(psi::outfile," %16.15f\n",energies[i][j]);
      }
   }

   return Expansion->Energy(Systems, energies);
}

void LibFragHelper::NMer_Helper(const int N) {
   DaOptions.MBEOrder=N;
   Expansion->SetN(N);
   if (!Expansion->IsGMBE()) {
      Expansion->MakeIntersections(Systems);
      Systems.push_back(NMerSet());
      Expansion->MakeNmers(Systems[0], Systems.back());
   }
   else {	//We make the n-Mers first for the GMBE
      Systems.push_back(NMerSet());
      Expansion->MakeNmers(Systems[0], Systems.back());
      Expansion->MakeIntersections(Systems);
      for (int i=0; i<Systems.size(); i++) {
         for (int j=0; j<Systems[i].size(); j++) {
            Systems[i][j]->print_out();
         }
      }
   }
   if (DaOptions.BMethod!=NO_BSSE) {
      for (int order=1; order<Systems.size(); order++)
         BSSEFactory->AddBSSEJobs(Systems[order]);
   }
}

int LibFragHelper::IsGMBE() {
   return Expansion->IsGMBE();
}

void LibFragHelper::Synchronize(boost::python::str& Comm, const int N) {
   if (N==0) {
      //Synchronize MO coefficients
      std::string comm;
      ToC(comm,Comm);
      std::vector<psi::MOFile> tempfiles;
      int NProc=psi::WorldComm->nproc(comm);
      int remainder=Systems[0].size()%NProc;
      int batchsize=(Systems[0].size()-remainder)/NProc;
      int me=psi::WorldComm->me(comm);
      for (int i=0; i<NProc; i++) {
         for (int j=0; j<batchsize; j++) {
            if (me==i) {
               MOFiles[j].Broadcast(comm, i);
               tempfiles.push_back(MOFiles[i]);
            }
            else {
               tempfiles.push_back(psi::MOFile());
               tempfiles.back().Receive(comm, i);
            }
         }
      }
      for (int i=0; i<remainder; i++) {
         if (me==i) {
            MOFiles[batchsize].Broadcast(comm, i);
            tempfiles.push_back(MOFiles[batchsize]);
         }
         else {
            tempfiles.push_back(psi::MOFile());
            tempfiles.back().Receive(comm, i);
         }
      }
      MOFiles=tempfiles;
   }
}

}//End namespace LibFrag
