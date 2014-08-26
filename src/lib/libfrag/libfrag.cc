/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
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
#include "libparallel/TableSpecs.h"

typedef std::vector<int> SVec;
typedef std::string str;
typedef const int cint;
typedef std::vector<boost::shared_ptr<double[]> > GMatrix; //Ghetto matrix


namespace psi{
typedef std::vector<MOFile> MOtype;
namespace LibFrag {

typedef MBEFrag* pSet;
typedef std::vector<pSet> FragSet;
typedef boost::shared_ptr<Embedder> ShareEmbed;

ShareEmbed LibFragHelper::EmbedFactory_=ShareEmbed();

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
   return baseGetCall(NMer,N,Systems_,false);
}

PyList LibFragHelper::GetGhostNMerN(const int NMer, const int N){
   return baseGetCall(NMer,N,Systems_,true);
}

void LibFragHelper::GatherData() {
   if (!this->IsGMBE()) {
      MOFiles_.push_back(psi::MOFile());
      MOFiles_.back().Read();
   }
   if (EmbedFactory_) {
      ChargeSets_.push_back(
            psi::Process::environment.wavefunction()->atomic_point_charges());
   }

}
int LibFragHelper::Iterate(cint itr) {
   return (EmbedFactory_ ? EmbedFactory_->Iterate(itr) : false);
}

void LibFragHelper::WriteMOs(cint N, cint x) {
   if (!this->IsGMBE()) {
      if (N>0) {
         SharedFrag nmer=Systems_[N][x];
         psi::MOFile File2Write(MOFiles_[nmer->ParentI(0)]);
         for (int i=1; i<=N; i++) {
            int frag=nmer->ParentI(i);
            psi::MOFile temp=File2Write.DirectSum(MOFiles_[frag]);
            File2Write=temp;
         }
         //File2Write.print_out();
         File2Write.Write();
      }
      else MOFiles_[x].Write();
   }

}

void LibFragHelper::Fragment_Helper(PyStr& FragMethod, cint N,
      PyStr& EmbedMethod, PyStr& CapMethod, PyStr& BSSEMethod) {
   str fname,bname,ename,cname;
   ToC(bname, BSSEMethod);
   ToC(fname, FragMethod);
   ToC(ename, EmbedMethod);
   ToC(cname, CapMethod);
   DaOptions_=boost::shared_ptr<LibFragOptions>(
         new LibFragOptions(N,fname,ename,cname,bname));
   DaOptions_->PrintOptions();
   SharedMol AMol=psi::Process::environment.molecule();
   Systems_.push_back(NMerSet());

   ///Fragment the system
   boost::shared_ptr<Fragmenter> FragFactory=
         DaOptions_->FOptions().MakeFactory(AMol);
   FragProps FProp=FragFactory->Fragment(AMol, Systems_[0]);

   ///Set-up either a MBE or a GMBE based on whether frags are disjoint
   Expansion_=(
         FProp.Disjoint_ ? boost::shared_ptr<GMBE>(new MBE(N)) :
               boost::shared_ptr<GMBE>(new GMBE(N)));

   ///If the user wants a one-body GMBE calc, we hack a bit here and add
   ///the intersections into the fragments
   if (Expansion_->IsGMBE()&&N==1) {
      Systems_.push_back(Systems_[0]);
      Expansion_->MakeIntersections(Systems_);
      Systems_[0]=Systems_[1];
      Systems_[0].insert(Systems_[0].end(), Systems_[2].begin(),
            Systems_[2].end());
   }

   ///If we broke bonds add some caps
   if (FProp.Severed_) {
      boost::shared_ptr<Capper> CapFactory=
            DaOptions_->COptions().MakeFactory(AMol);
      CapFactory->MakeCaps(Systems_[0]);
   }

   ///If the user wants embedding, anything that's not a cap or an atom is
   ///a charge
   EmbedFactory_=DaOptions_->EOptions().MakeFactory(AMol);

   ///Actual embedding of the monomers is done in synchronize, as this is
   ///the first time we have the charges/densities

   ///Add in any BSSE corrections we may need
   boost::shared_ptr<BSSEer>BSSEFactory=
         DaOptions_->BOptions().MakeFactory(AMol);

   ///Caps and real atoms are real, embedding and ghost atoms are ghosts
   ///in particular this means we have ghost atoms, literally on top of
   ///point charges, this needs tested...
   if (BSSEFactory) BSSEFactory->AddBSSEJobs(Systems_[0]);
}

PyList LibFragHelper::Embed_Helper(cint N, cint x) {
   PyList DaList;
   if (EmbedFactory_) {
      if (EmbedFactory_->HaveCharges()) {
         SharedFrag NMer=Systems_[N][x];
         for (int i=0; i<NMer->Charges().size(); i++)
            DaList.append(NMer->Charges()[i]);
      }
   }
   return DaList;
}

str LibFragHelper::Cap_Helper(int order, int frag) {
   std::stringstream caps;
   SharedFrag Frag=Systems_[order][frag];
   for (int i=0; i<Frag->Caps().size(); i++){
      int ActualCap=Frag->Caps()[i];
      caps<<Frag->Caps().Object(ActualCap)->print_out();
   }
   return caps.str();
}

GMatrix ExtractEnergies(PyList& Energies){
   GMatrix energies;
   for (int i=0; i<len(Energies); i++) {
      int size=len(boost::python::extract<PyList>(Energies[i]));
      boost::shared_ptr<double[]> temp(new double[size]);
      energies.push_back(temp);
      for (int j=0; j<size; j++)
         energies[i][j]=boost::python::extract<double>(Energies[i][j]);
   }
   return energies;
}

void LibFragHelper::PrintEnergy(PyList& Energies,const int N){
   GMatrix energies=ExtractEnergies(Energies);
   std::vector<std::string> Titles;
   if(N==1) Titles.push_back("Monomer #");
   else{
      std::stringstream temp;
      temp<<N<<"-Mer #";
      Titles.push_back(temp.str());
   }
   Titles.push_back("Energy (a.u.)");
   std::vector<std::string> Parents;
   int NRows=Systems_[N-1].size();
   for(int i=0;i<NRows;i++){
      Parents.push_back(Systems_[N-1][i]->PrintParents());
   }
   TableSpecs<std::string,double> Specs(NRows);
   Specs.Init(&Parents[0],&(energies[N-1][0]));
   Specs.SetTitles(Titles);
   (*outfile)<<(Specs.Table());
   (*outfile)<<std::endl;
}


double LibFragHelper::CalcEnergy(PyList& Energies,bool IsCorr) {
   GMatrix energies=ExtractEnergies(Energies);
   return Expansion_->Energy(Systems_, energies,IsCorr);
}

//bool LibFragHelper::SelfIntOff() {
 //  return (EmbedFactory);
//}

void LibFragHelper::NMer_Helper(cint N) {
   DaOptions_->MBEOrder_=N;
   Expansion_->SetN(N);
   if (!Expansion_->IsGMBE()) {
      Expansion_->MakeIntersections(Systems_);
      Systems_.push_back(NMerSet());
      Expansion_->MakeNmers(Systems_[0], Systems_.back());
   }
   else {	//We make the n-Mers first for the GMBE
      Systems_.push_back(NMerSet());
      Expansion_->MakeNmers(Systems_[0], Systems_.back());
      Expansion_->MakeIntersections(Systems_);
   }

}

int LibFragHelper::IsGMBE() {
   return Expansion_->IsGMBE();
}

void LibFragHelper::BroadCastWrapper(cint i, cint j, str& comm,
      MOtype& tempfiles, GMatrix& tempChargeSets, bool bcast) {
   if (bcast) {
      MOFiles_[j].Broadcast(comm, i);
      tempfiles.push_back(MOFiles_[j]);
   }
   else {
      tempfiles.push_back(psi::MOFile());
      tempfiles.back().Receive(comm, i);
   }
   if (EmbedFactory_) {
      int natoms=Systems_[0][tempfiles.size()]->Atoms().size();
      if (bcast) {
         psi::WorldComm->bcast(&(ChargeSets_[j][0]), natoms, i, comm);
         tempChargeSets.push_back(ChargeSets_[j]);
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
         for (int frag=0; frag<Systems_[0].size(); frag++) {
            for (int atom=0; atom<Systems_[0][frag]->Atoms().size(); atom++) {
               int atomI=Systems_[0][frag]->Atoms()[atom];
               EmbedFactory_->SetCharge(atomI, ChargeSets_[frag][atom]);
            }
         }
         EmbedFactory_->print_out();
         EmbedFactory_->Embed(Systems_[0]);
      }
      if (itr>=1) {
         ///Update our fragment MO coefficients (we've been putting the new
         ///MO coefs behind the old ones)
         std::vector<psi::MOFile> temp(MOFiles_.begin()+Systems_[0].size(),
               MOFiles_.end());
         MOFiles_=temp;
      }
      if (Parallel) {
         //Synchronize MO coefficients
         GMatrix tempChargeSets;
         std::vector<psi::MOFile> tempfiles;
         int NProc=psi::WorldComm->nproc(comm);
         int remainder=Systems_[0].size()%NProc;
         int batchsize=(Systems_[0].size()-remainder)/NProc;
         int me=psi::WorldComm->me(comm);
         for (int i=0; i<NProc; i++) {
            for (int j=0; j<batchsize; j++)
               BroadCastWrapper(i, j, comm, tempfiles, tempChargeSets, (me==i));
         }
         for (int i=0; i<remainder; i++)
            BroadCastWrapper(i, batchsize, comm, tempfiles, tempChargeSets,
                  (me==i));
         MOFiles_=tempfiles;
      }
   }
}

bool LibFragHelper::RunFrags(){return Expansion_->RunFrags();}
}}	//End namespaces
