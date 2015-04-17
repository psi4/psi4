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
#include <set>
#include <sstream>
#include "FragSysGuts.h"
#include "FragmentedSys.h"
#include "Utils/Fragmenter.h"
#include "LibFragFragment.h"
#include "PsiMap.h"
#include "MoleculeTypes.h"
#include "Utils/GeomManipulator.h"
namespace psi {
namespace LibMolecule {

typedef boost::shared_ptr<const Molecule> SharedMol;
typedef boost::shared_ptr<Fragment> SharedFrag;
typedef std::vector<SharedFrag> vSharedFrag;
typedef vSharedFrag::const_iterator vSFragcItr;
typedef vSharedFrag::iterator vSFragItr;
typedef boost::shared_ptr<NMers> SharedNMers;
typedef boost::shared_ptr<SNOMolecule> SharedSNO;
typedef PsiMap<int, SharedSNO> SNO_t;

static void UpdateNMerI(NMers& NMers_,const SerialNumber& SN1,
                 const SerialNumber& SN2,const int N){
   NMers_.ScaleFacts_.Change(N, SN1, 1.0);
   NMers_.SNLookUp_[SN2]=SN1;
}

static boost::shared_ptr<Fragmenter> PickFragmenter(boost::shared_ptr<const Molecule> Mol){
   boost::shared_ptr<Fragmenter> DaPointer;
   std::string FragType=
         psi::Process::environment.options["FRAG_METHOD"].to_string();
   if(FragType=="BOND_BASED"){
      DaPointer=boost::shared_ptr<BondFragmenter>(
            new BondFragmenter(Mol));
   }
   else if(FragType=="MONOMER_BASED"){
      DaPointer=boost::shared_ptr<MonomerFragmenter>(
            new MonomerFragmenter(Mol));
   }
   return DaPointer;
}


FragSysGuts::FragSysGuts(boost::shared_ptr<Molecule> System2Frag, const int N) :
      NMers_(N) {
   boost::shared_ptr<Fragmenter> FragMaker=
         PickFragmenter(System2Frag);
   vSharedFrag Fragments=FragMaker->MakeFrags();
   if (N>1)MakeNMers(Fragments);
   NMers_[0]=Fragments;
   for (int i=0; i<N; i++) {
      NMers::NMerItr_t NMerI=NMers_.NMerBegin(i+1),NMerEnd=NMers_.NMerEnd(i+1);
      for(;NMerI!=NMerEnd; ++NMerI)
         UpdateNMerI(NMers_,(*NMerI)->GetSN(),(*NMerI)->GetSN(),i);
   }
}

//A 27 by 3 array of the three mirror plane rotations
const double Mirror[]={-1.0, 0.0, 0.0,0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, //MirrorY
      0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,0.0, 1.0, 0.0, 0.0, 0.0, -1.0};

static void CheckNMer(NMers& NMers_,const int N) {
   NMers::NMerItr_t NMerI=NMers_.NMerBegin(N+1),
                         NMerEnd=NMers_.NMerEnd(N+1);
   PsiMap<int, SerialNumber> Val2SN;
   vSharedFrag ReturnFrags;
   SNO_t SymmUnique;
   for (int counter=0; NMerI!=NMerEnd; ++NMerI) {
      SharedSNO temp(new SNOMolecule((*(*NMerI))));
      SNO_t::const_iterator SNOI=SymmUnique.begin(),SNOEnd=SymmUnique.end();
      bool IsGood=true;
      for (; SNOI!=SNOEnd; ++SNOI) {
         if ((*(SNOI->second))==(*temp)) {
            UpdateNMerI(NMers_,Val2SN[SNOI->first],(*NMerI)->GetSN(),N);
            IsGood=false;
            break;
         }
      }
      if (!IsGood) continue;
      for (int k=0; k<3; k++) {
         SharedSNO temp2(new SNOMolecule(*temp));
         GeomManipulator Manip(temp2.get());
         Manip.Rotate(&Mirror[k*9]);
         Manip.Set();
         for (SNOI=SymmUnique.begin(); SNOI!=SNOEnd; ++SNOI) {
            if ((*(SNOI->second))==(*temp2)) {
               UpdateNMerI(NMers_,Val2SN[SNOI->first],(*NMerI)->GetSN(),N);
               IsGood=false;
               break;
            }
         }
         if (!IsGood) break;
      }
      if (IsGood) {
         SymmUnique[counter]=temp;
         const SerialNumber& SN=(*NMerI)->GetSN();
         Val2SN[counter]=SN;
         UpdateNMerI(NMers_,SN,SN,N);
         ReturnFrags.push_back(*NMerI);
      }
   }
   NMers_[N]=ReturnFrags;
}

FragSysGuts::FragSysGuts(const boost::shared_ptr<SuperCell> System2Frag, const int N) :
      NMers_(N) {
   SharedMol FullSys=System2Frag;
   boost::shared_ptr<Fragmenter> FragMaker=
         PickFragmenter(System2Frag->GetUnitCell()),
         FragMaker2=PickFragmenter(FullSys);
   vSharedFrag UnitFrags=FragMaker->MakeFrags(),
               AllFrags=FragMaker2->MakeFrags();
   vSFragItr FragI=UnitFrags.begin(),FragEnd=UnitFrags.end();
   //Here we exploit the fact that the atoms in the unitcell are the first
   //atoms in the supercell (i.e. the numbering didn't change)
   vSharedFrag TempFrags;
   for(; FragI!=FragEnd; ++FragI) {
      TempFrags.push_back(SharedFrag(new Fragment(FullSys)));
      (*TempFrags.back())+=(*(*FragI));
   }
   //Get the n-mers with at least one monomer in unitcell
   if (N>1)MakeNMers(TempFrags, AllFrags);
   NMers_[0]=TempFrags;
   for(int i=0;i<N;i++)CheckNMer(NMers_,i);
   //Now for all orders less than N need the full set of n-mers
   ScaleFactors TempSF=NMers_.ScaleFacts_;
   std::vector<boost::shared_ptr<Fragment> > TempNMers=NMers_[N-1];
   if (N>1) {
      NMers BlankNMers(N);
      NMers_=BlankNMers;
      //Now get all n-mers
      MakeNMers(AllFrags);
      NMers_[0]=AllFrags;
      //Abuse this call to find symmetry unique n-mers up to N-1
      //(We are going to use the previously found N-mers)
      for(int i=0;i<N-1;i++)CheckNMer(NMers_,i);
      //Only want the N-Mers with at least one monomer in the UC
      NMers_[N-1]=TempNMers;
      //Want our old scalefactors
      NMers_.ScaleFacts_=TempSF;

      for(int i=0;i<N-1;i++){
         NMers::NMerItr_t NMerI=NMers_.NMerBegin(i+1),
                       NMerIEnd=NMers_.NMerEnd(i+1);
         for(;NMerI!=NMerIEnd;++NMerI){
            const SerialNumber& SN=(*NMerI)->GetSN();
            if(NMers_.ScaleFacts_[i].count(SN)!=1)
               NMers_.ScaleFacts_[i][SN]=0.0;
         }
      }
      //Set the SNLookUp for the N-mers
      NMers::NMerItr_t NMerI=NMers_.NMerBegin(N),NMerIEnd=NMers_.NMerEnd(N);
      for(;NMerI!=NMerIEnd;++NMerI){
         const SerialNumber& SN=(*NMerI)->GetSN();
         NMers_.SNLookUp_[SN]=SN;
      }
   }
}

/*NMers CleanUp(NMers& NMersIn) {
   NMers NMers2(NMersIn.N());
   for (int i=2; i<=NMersIn.N(); i++) {
      NMers::NMerItr_t NMerI=NMersIn.NMerBegin(i),NMerEnd=NMersIn.NMerEnd(i);
      for (int offset=0; NMerI!=NMerEnd; ++NMerI, ++offset) {
         if (!(*NMerI)) continue;
         bool good=true;
         int offset2=0;
         for (NMers::NMerItr_t NMerJ=NMerI+1; NMerJ!=NMerEnd;
               ++NMerJ, ++offset2) {
            //Check if we already removed NMerJ
            if (!(*NMerJ)) continue;
            if ((*(*NMerI))<=(*(*NMerJ))) {
               std::cout<<offset<<"="<<offset2<<std::endl;
               NMerJ->reset();
               good=false;
               break;
            }
         }
         if (good) NMers2.AddNMer(i, (*NMerI));
      }
   }
   return NMers2;
}*/

static void Recurse(const std::vector<double>& Distances,
              vSFragcItr FragJ, vSFragcItr& FragEnd2, const Fragment& FragI,
      int n, NMers& NMers) {
   if (n>NMers.N()) return;
   for (; FragJ!=FragEnd2; ++FragJ) {
      boost::shared_ptr<Fragment> temp(new Fragment(FragI));
      (*temp)+=(*(*FragJ));
      if((int)psi::Process::environment.options["MBE_DISTANCE_THRESHOLDS"].size()>=(n-1)){
         double DistanceThresh=
            psi::Process::environment.options["MBE_DISTANCE_THRESHOLDS"]
                        [n-2].to_double()*LibMoleculeBase().AngToBohr();
         const SerialNumber& SN=temp->GetSN();
         SerialNumber::const_iterator FragK=SN.begin(),FragL,
                                       FragEndK=SN.end();
         bool IsGood=true;
         for(;FragK!=FragEndK&&IsGood;++FragK){
            FragL=FragK;
            ++FragL;
            for(;FragL!=FragEndK&&IsGood;++FragL){
               if(Distances[(*FragK)*NMers[0].size()+(*FragL)]>DistanceThresh)IsGood=false;
            }
         }
         if(!IsGood)continue;
      }
      NMers.AddNMer(n, temp);
      Recurse(Distances,FragJ+1, FragEnd2, (*temp), n+1, NMers);
   }
}

typedef std::vector<std::vector<double> > CoM_t;

static void CalcCoM(int i,vSharedFrag::const_iterator FragI, CoM_t& CoMs){
   MolItr AtomI=(*FragI)->Begin(),AtomEnd=(*FragI)->End();
   double Total=0.0;
   for(;AtomI!=AtomEnd;++AtomI){
      Total+=(*AtomI)->Mass();
      for(int j=0;j<3;j++)CoMs[i][j]+=(*AtomI)->Mass()*(*(*AtomI))[j];
   }
   for(int j=0;j<3;j++)CoMs[i][j]/=Total;
}

static std::vector<double> GetFragDistances(const vSharedFrag& Set){
   int NFrags=Set.size();
   std::vector<double> Distances(NFrags*NFrags,0.0);
   vSharedFrag::const_iterator FragI=Set.begin(),FragEnd=Set.end(),FragJ;
   CoM_t CoMs(NFrags,std::vector<double>(3,0.0));
      for(int i=0;FragI!=FragEnd;++FragI,++i){
         if(i==0)CalcCoM(i,FragI,CoMs);
         int j=i+1;
         FragJ=FragI;
         ++FragJ;
         for(;FragJ!=FragEnd;++FragJ,++j){
            if(i==0)CalcCoM(j,FragJ,CoMs);
            for(int k=0;k<3;k++){
               double quant=CoMs[i][k]-CoMs[j][k];
               Distances[i*NFrags+j]+=quant*quant;
            }
            Distances[j*NFrags+i]=sqrt(Distances[i*NFrags+j]);
            Distances[i*NFrags+j]=Distances[j*NFrags+i];
         }
      }
      return Distances;
}


void FragSysGuts::MakeNMers(const vSharedFrag& Set1, const vSharedFrag& Set2) {
   std::vector<double> Distances= GetFragDistances(Set2);
   vSharedFrag::const_iterator FragI=Set1.begin(),FragEnd=Set1.end(),
         FragJ,FragEnd2=Set2.end();
   for(int offset=1; FragI!=FragEnd; ++FragI, ++offset) {
      FragJ=Set2.begin();
      FragJ+=offset;
      Recurse(Distances,FragJ, FragEnd2, (*(*FragI)), 2, NMers_);
   }
   //If the frags are non-disjoint we need to call this
   //NMers_=boost::shared_ptr<NMers>(new NMers(CleanUp(*NMers_)));
}

void FragSysGuts::MakeNMers(const vSharedFrag& Set) {MakeNMers(Set, Set);}

std::string NMers::PrintOut(const int Value) const {
   Base_t::const_iterator OrderI=Base_t::begin(),OrderN=Base_t::end();
   std::stringstream Result;
   for (int n=1; OrderI!=OrderN; ++OrderI, ++n) {
      NMers::cNMerItr_t NMerI=OrderI->begin(),NMerEnd=OrderI->end();
      for (int nn=0; NMerI!=NMerEnd; ++NMerI, ++nn) {
         if (n>1) Result<<n<<"-mer ";
         else Result<<"monomer ";
         Result<<(*NMerI)->GetSN().PrintOut();
         Result<<std::endl;
         //" Multiplicity :"<<ScaleFacts_(n-1,nn)<<std::endl;
         Result<<(*NMerI)->PrintOut(Value);
         Result<<"---"<<std::endl;
      }
   }
   return Result.str();

}

}}   //End namespaces

