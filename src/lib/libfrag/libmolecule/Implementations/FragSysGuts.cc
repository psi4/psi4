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
#include "FragSysGuts.h"
#include <sstream>
#include "FragmentedSys.h"
#include "Utils/Fragmenter.h"
#include "Fragment.h"
#include "PsiMap.h"
#include "MoleculeTypes.h"
#include "Utils/GeomManipulator.h"
namespace psi{
namespace LibMolecule{

typedef boost::shared_ptr<const Molecule> SharedMol;
typedef boost::shared_ptr<Fragment> SharedFrag;
typedef std::vector<SharedFrag> vSharedFrag;
typedef vSharedFrag::const_iterator vSFragcItr;
typedef vSharedFrag::iterator vSFragItr;
typedef boost::shared_ptr<NMers> SharedNMers;
typedef boost::shared_ptr<SNOMolecule> SharedSNO;
typedef PsiMap<int,SharedSNO> SNO_t;
FragSysGuts::FragSysGuts(const Molecule& System2Frag,const int N):
   NMers_(N){
   BondFragmenter FragMaker(SharedMol(new Molecule(System2Frag)));
   vSharedFrag Fragments=FragMaker.MakeFrags();
   vSFragItr FragI=Fragments.begin(),FragEnd=Fragments.end();
   for(int i=0;FragI!=FragEnd;++FragI,++i){
      NMers_.AddNMer(1,*FragI);
      NMers_.ScaleFacts_.Change(0,i,1.0);
   }
   if(N>1){
      MakeNMers(Fragments);
      for(int i=1;i<N;i++){
         NMers::NMerItr_t NMerI=NMers_.NMerBegin(i+1),
                          NMerEnd=NMers_.NMerEnd(i+1);
         for(int j=0;NMerI!=NMerEnd;++NMerI,j++)
            NMers_.ScaleFacts_.Change(i,j,1.0);
      }
   }
}

const double Mirror[]={-1.0, 0.0, 0.0,//MirrorX
                              0.0, 1.0, 0.0,
                              0.0, 0.0, 1.0,
                              1.0, 0.0, 0.0,//MirrorY
                              0.0,-1.0, 0.0,
                              0.0, 0.0, 1.0,
                              1.0, 0.0, 0.0,//MirrorZ
                              0.0, 1.0, 0.0,
                              0.0, 0.0,-1.0
};

vSharedFrag CheckNMer(const int N,vSFragItr& NMerI,
                                  vSFragItr& NMerEnd,
                            ScaleFactors& ScaleFacts_){
   vSharedFrag ReturnFrags;
   SNO_t SymmUnique;
   for(int counter=0;NMerI!=NMerEnd;++NMerI){
      SharedSNO temp(new SNOMolecule((*(*NMerI))));
      SNO_t::const_iterator SNOI=SymmUnique.begin(),
                            SNOEnd=SymmUnique.end();
      bool IsGood=true;
      for(;SNOI!=SNOEnd;++SNOI){
         if((*(SNOI->second))==(*temp)){
            ScaleFacts_.Change(N,SNOI->first,1.0);
            IsGood=false;
            break;
         }
      }
      if(!IsGood)continue;
      for(int k=0;k<3;k++){
         SharedSNO temp2(new SNOMolecule(*temp));
         GeomManipulator Manip(temp2.get());
         Manip.Rotate(&Mirror[k*9]);
         Manip.Set();
         for(SNOI=SymmUnique.begin();SNOI!=SNOEnd;++SNOI){
            if((*(SNOI->second))==(*temp2)){
               ScaleFacts_.Change(N,SNOI->first,1.0);
               IsGood=false;
               break;
             }
         }
         if(!IsGood)break;
      }
      if(IsGood){
         SymmUnique[counter]=temp;
         ScaleFacts_.Change(N,counter++,1.0);
         ReturnFrags.push_back(*NMerI);
      }
   }
   return ReturnFrags;
}

FragSysGuts::FragSysGuts(const SuperCell& System2Frag,const int N):
   NMers_(N){
   SharedMol FullSys(new Molecule(System2Frag));
   BondFragmenter FragMaker(System2Frag.GetUnitCell()),FragMaker2(FullSys);
   vSharedFrag UnitFrags=FragMaker.MakeFrags(),AllFrags=FragMaker2.MakeFrags();
   vSFragItr FragI=UnitFrags.begin(),FragEnd=UnitFrags.end();
   //Here we exploit the fact that the atoms in the unitcell are the first
   //atoms in the supercell (i.e. the numbering didn't change)
   vSharedFrag TempFrags;
   for(int i=0;FragI!=FragEnd;++FragI,++i){
      TempFrags.push_back(SharedFrag(new Fragment(FullSys)));
      (*TempFrags.back())+=(*(*FragI));
   }
   if(N>1)MakeNMers(TempFrags,AllFrags);
   vSFragItr TempI=TempFrags.begin(),TempEnd=TempFrags.end();
   NMers TempNMers(N);
   vSharedFrag Fragments=CheckNMer(0,TempI,TempEnd,TempNMers.ScaleFacts_);
   for(int i=2;i<=N;i++){
      NMers::NMerItr_t NMerI=NMers_.NMerBegin(i),NMerEnd=NMers_.NMerEnd(i);
      std::vector<SharedFrag>TempFrags2=
            CheckNMer(i-1,NMerI,NMerEnd,TempNMers.ScaleFacts_);
      std::vector<SharedFrag>::iterator FragJ,FragJEnd=TempFrags2.end();
      for(FragJ=TempFrags2.begin();FragJ!=FragJEnd;++FragJ)
         TempNMers.AddNMer(i,(*FragJ));
   }
   NMers_=TempNMers;
   FragI=Fragments.begin(),FragEnd=Fragments.end();
   for(int i=0;FragI!=FragEnd;++FragI,++i){
      NMers_.AddNMer(1,*FragI);
   }
}

NMers CleanUp(NMers& NMersIn){
   NMers NMers2(NMersIn.N());
   for(int i=2;i<=NMersIn.N();i++){
      NMers::NMerItr_t NMerI=NMersIn.NMerBegin(i),
                       NMerEnd=NMersIn.NMerEnd(i);
      for(int offset=0;NMerI!=NMerEnd;++NMerI,++offset){
         if(!(*NMerI))continue;
         bool good=true;
         int offset2=0;
         for(NMers::NMerItr_t NMerJ=NMerI+1;NMerJ!=NMerEnd;++NMerJ,++offset2){
            //Check if we already removed NMerJ
            if(!(*NMerJ))continue;
            if((*(*NMerI))<=(*(*NMerJ))){
               std::cout<<offset<<"="<<offset2<<std::endl;
               NMerJ->reset();
               good=false;
               break;
            }
         }
         if(good)NMers2.AddNMer(i,(*NMerI));
      }
   }
   return NMers2;
}

void Recurse(vSFragcItr FragJ,
             vSFragcItr& FragEnd2,
             const Fragment& FragI,int n,NMers& NMers){
   if(n>NMers.N())return;
   for(;FragJ!=FragEnd2;++FragJ){
      boost::shared_ptr<Fragment> temp(new Fragment(FragI));
      (*temp)+=(*(*FragJ));
      NMers.AddNMer(n,temp);
      Recurse(FragJ+1,FragEnd2,(*temp),n+1,NMers);
   }
}

void FragSysGuts::MakeNMers
(const vSharedFrag& Set1, const vSharedFrag& Set2){
   vSharedFrag::const_iterator FragI=Set1.begin(),FragEnd=Set1.end();
   for(int offset=1;FragI!=FragEnd;++FragI,++offset){
      vSFragcItr FragJ=Set2.begin(),FragEnd2=Set2.end();
      FragJ+=offset;
      Recurse(FragJ,FragEnd2,(*(*FragI)),2,NMers_);
   }
   //If the frags are non-disjoint we need to call this
   //NMers_=boost::shared_ptr<NMers>(new NMers(CleanUp(*NMers_)));
}

void FragSysGuts::MakeNMers(const vSharedFrag& Set){
   MakeNMers(Set,Set);
}

std::string NMers::PrintOut(const int Value)const{
   Base_t::const_iterator OrderI=Base_t::begin(),OrderN=Base_t::end();
   std::stringstream Result;
   for(int n=1;OrderI!=OrderN;++OrderI,++n){
      NMers::cNMerItr_t NMerI=OrderI->begin(),NMerEnd=OrderI->end();
      for(int nn=0;NMerI!=NMerEnd;++NMerI,++nn){
         std::set<long int>::const_iterator ID=(*NMerI)->SerialNumber().begin(),
               IDEnd=(*NMerI)->SerialNumber().end();
         if(n>1)Result<<n<<"-mer ";
         else Result<<"monomer ";
         for(;ID!=IDEnd;++ID)Result<<(*ID)<<" ";
         Result<<std::endl;
         //" Multiplicity :"<<ScaleFacts_(n-1,nn)<<std::endl;
         Result<<(*NMerI)->PrintOut(Value);
         Result<<"---"<<std::endl;
      }
   }
   return Result.str();

}

}}//End namespaces


