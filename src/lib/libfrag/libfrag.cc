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
#include "libmints/wavefunction.h"
#include "MBE.h"
#include "BSSEer.h"
#include "Capper.h"
#include "process.h"
#include "Embedder.h"
#include "libciomr/libciomr.h"
#include "MBEFrag.h"

typedef std::vector<int> SVec;
typedef std::string str;
typedef const int cint;
typedef std::vector<boost::shared_ptr<double[]> > GMatrix; //Ghetto matrix
typedef boost::python::list PyList;
typedef boost::python::str PyStr;

namespace psi{
typedef std::vector<MOFile> MOtype;
namespace LibFrag {

typedef MBEFrag* pSet;
typedef std::vector<pSet> FragSet;
typedef boost::shared_ptr<Embedder> ShareEmbed;

ShareEmbed LibFragHelper::EmbedFactory=ShareEmbed();

//Handy template for switching python things to c++ things
template <typename T, typename S>
void ToC(T& cvalue, S& pyvalue) {
   cvalue=boost::python::extract<T>(pyvalue);
}

LibFragHelper::LibFragHelper(){
   tstart();
}

LibFragHelper::~LibFragHelper(){
   tstop();
}

PyList baseGetCall(cint NMer, cint N, std::vector<NMerSet>& Systems,
      bool IsGhost) {
   PyList DaList;
   SharedFrag DaNMer=(Systems[NMer])[N];
   if(IsGhost){
      for (int i=0; i<DaNMer->Ghosts().size(); i++)
         DaList.append(DaNMer->Ghosts()[i]);
   }
   else{
      for (int i=0;i<DaNMer->Atoms().size();i++)
         DaList.append(DaNMer->Atoms()[i]);
   }
   return DaList;
}

PyList LibFragHelper::GetNMerN(const int NMer,const int N){
   return baseGetCall(NMer,N,Systems,false);
}

PyList LibFragHelper::GetGhostNMerN(const int NMer, const int N){
   return baseGetCall(NMer,N,Systems,true);
}

void LibFragHelper::GatherData() {
   if (!this->IsGMBE()) {
      MOFiles.push_back(psi::MOFile());
      MOFiles.back().Read();
   }
   if (EmbedFactory) {
      ChargeSets.push_back(
            psi::Process::environment.wavefunction()->atomic_point_charges());
   }

}
int LibFragHelper::Iterate(cint itr) {
   return (EmbedFactory ? EmbedFactory->Iterate(itr) : false);
}

void LibFragHelper::WriteMOs(cint N, cint x) {
   if (!this->IsGMBE()) {
      if (N>0) {
         SharedFrag nmer=Systems[N][x];
         psi::MOFile File2Write(MOFiles[nmer->ParentI(0)]);
         for (int i=1; i<=N; i++) {
            int frag=nmer->ParentI(i);
            psi::MOFile temp=File2Write.DirectSum(MOFiles[frag]);
            File2Write=temp;
         }
         //File2Write.print_out();
         File2Write.Write();
      }
      else MOFiles[x].Write();
   }

}

void LibFragHelper::Fragment_Helper(PyStr& FragMethod, cint N,
      PyStr& EmbedMethod, PyStr& CapMethod, PyStr& BSSEMethod) {
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
         FProp.Disjoint_ ? boost::shared_ptr<GMBE>(new MBE(N)) :
               boost::shared_ptr<GMBE>(new GMBE(N)));

   ///If the user wants a one-body GMBE calc, we hack a bit here and add
   ///the intersections into the fragments
   if (Expansion->IsGMBE()&&N==1) {
      Systems.push_back(Systems[0]);
      Expansion->MakeIntersections(Systems);
      Systems[0]=Systems[1];
      Systems[0].insert(Systems[0].end(), Systems[2].begin(), Systems[2].end());
   }

   ///If we broke bonds add some caps
   if (FProp.Severed_) {
      CapFactory=DaOptions.MakeCapFactory(AMol);
      CapFactory->MakeCaps(Systems[0]);
   }

   ///If the user wants embedding, anything that's not a cap or an atom is
   ///a charge
   EmbedFactory=DaOptions.MakeEmbedFactory(AMol);
   ///Actual embedding of the monomers is done in synchronize, as this is
   ///the first time we have the charges/densities

   ///Add in any BSSE corrections we may need
   BSSEFactory=DaOptions.MakeBSSEFactory(AMol->natom());

   ///Caps and real atoms are real, embedding and ghost atoms are ghosts
   ///in particular this means we have ghost atoms, literally on top of
   ///point charges, this needs tested...
   if (BSSEFactory) BSSEFactory->AddBSSEJobs(Systems[0]);
}

PyList LibFragHelper::Embed_Helper(cint N, cint x) {
   PyList DaList;
   if (EmbedFactory) {//Don't derefrence pointer if it's not set
      if (EmbedFactory->HaveCharges()) {
         SharedFrag NMer=Systems[N][x];
         for (int i=0; i<NMer->Charges().size(); i++)
            DaList.append(NMer->Charges()[i]);
      }
   }
   return DaList;
}

str LibFragHelper::Cap_Helper(int order, int frag) {
   std::stringstream caps;
   //for (int i=0; i<Systems[order][frag]->Caps.size(); i++)
   //   caps<<Systems[order][frag]->Caps[i].print_out();
   return caps.str();
}

double LibFragHelper::CalcEnergy(PyList& Energies) {
   std::vector<double*> energies;
   for (int i=0; i<len(Energies); i++) {
      int size=len(boost::python::extract<PyList>(Energies[i]));
      energies.push_back(new double[size]);
      for (int j=0; j<size; j++)
         energies[i][j]=boost::python::extract<double>(Energies[i][j]);
   }
   std::string dashes;
   dashes=
         "------------------------------------------------------------------\n";

   this->PrintBanner('*',"N-mer energies");
   outfile->Printf( "");
   for (int i=0; i<Systems.size(); i++) {
      if (i==0) outfile->Printf( "Monomer #     Energy (a.u.)\n");
      else outfile->Printf( "%2d-mer  #     Energy (a.u.)\n", i+1);
      outfile->Printf( dashes.c_str());
      for (int j=0; j<Systems[i].size(); j++) {
         if (i!=0) Systems[i][j]->PrintParents();
         else outfile->Printf( "%d ", j);
         outfile->Printf( " %16.15f\n", energies[i][j]);
      }
   }

   return Expansion->Energy(Systems, energies);
}

bool LibFragHelper::SelfIntOff() {
   return (EmbedFactory);
}

void LibFragHelper::NMer_Helper(cint N) {
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
   }

}

int LibFragHelper::IsGMBE() {
   return Expansion->IsGMBE();
}

void LibFragHelper::BroadCastWrapper(cint i, cint j, str& comm,
      MOtype& tempfiles, GMatrix& tempChargeSets, bool bcast) {
   if (bcast) {
      MOFiles[j].Broadcast(comm, i);
      tempfiles.push_back(MOFiles[j]);
   }
   else {
      tempfiles.push_back(psi::MOFile());
      tempfiles.back().Receive(comm, i);
   }
   if (EmbedFactory) {
      int natoms=Systems[0][tempfiles.size()]->Atoms().size();
      if (bcast) {
         psi::WorldComm->bcast(&(ChargeSets[j][0]), natoms, i, comm);
         tempChargeSets.push_back(ChargeSets[j]);
      }
      else {
         boost::shared_ptr<double[]>(new double[natoms]);
         psi::WorldComm->bcast(&(tempChargeSets.back()[0]), natoms, i, comm);
      }
   }
}

void LibFragHelper::Synchronize(boost::python::str& Comm, cint N, cint itr) {
   std::string comm;
   ToC(comm, Comm);
   bool UpdateCharges=this->Iterate(itr);
   bool Parallel=(psi::WorldComm->nproc(comm)>1);
   if (N==0) {
      if (UpdateCharges) {
         for (int frag=0; frag<Systems[0].size(); frag++) {
            for (int atom=0; atom<Systems[0][frag]->Atoms().size(); atom++) {
               int atomI=Systems[0][frag]->Atoms()[atom];
               EmbedFactory->SetCharge(atomI, ChargeSets[frag][atom]);
            }
         }
         EmbedFactory->print_out();
         EmbedFactory->Embed(Systems[0]);
      }
      if (itr>=1) {
         ///Update our fragment MO coefficients (we've been putting the new
         ///MO coefs behind the old ones)
         std::vector<psi::MOFile> temp(MOFiles.begin()+Systems[0].size(),
               MOFiles.end());
         MOFiles=temp;
      }
      if (Parallel) {
         //Synchronize MO coefficients
         GMatrix tempChargeSets;
         std::vector<psi::MOFile> tempfiles;
         int NProc=psi::WorldComm->nproc(comm);
         int remainder=Systems[0].size()%NProc;
         int batchsize=(Systems[0].size()-remainder)/NProc;
         int me=psi::WorldComm->me(comm);
         for (int i=0; i<NProc; i++) {
            for (int j=0; j<batchsize; j++)
               BroadCastWrapper(i, j, comm, tempfiles, tempChargeSets, (me==i));
         }
         for (int i=0; i<remainder; i++)
            BroadCastWrapper(i, batchsize, comm, tempfiles, tempChargeSets,
                  (me==i));
         MOFiles=tempfiles;
      }
   }
}

bool LibFragHelper::RunFrags(){return Expansion->RunFrags();}
}}	//End namespaces
