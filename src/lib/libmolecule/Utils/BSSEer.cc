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
#include <vector>
#include "BSSEer.h"
#include "../LibFragMolecule.h"
#include "../AtomTypes.h"
#include "../FragmentedSys.h"
#include "../MoleculeTypes.h"
typedef unsigned int uint;
namespace psi {
namespace LibMolecule {

void FullBSSEer::CalcBSSE(NMers& Sys, uint Stop, uint Start) const {
   for (uint n=Start; n<=(Stop==0 ? Sys.N() : Stop); n++) {
      NMers::NMerItr_t FragI=Sys.NMerBegin(n),FragIEnd=Sys.NMerEnd(n);
      for (; FragI!=FragIEnd; ++FragI) {
         MolItr AtomI=(*FragI)->Begin(),AtomEndI=(*FragI)->End();
         //Start by reading the atoms into an array we don't touch
         //(adding ghosts will change the iterators)
         std::vector<boost::shared_ptr<const Atom> > TempAtoms;
         for (int counter=0; AtomI!=AtomEndI; ++AtomI, ++counter)
            TempAtoms.push_back(*AtomI);
         AtomI=(*FragI)->Mol()->Begin();
         AtomEndI=(*FragI)->Mol()->End();
         for (; AtomI!=AtomEndI; ++AtomI) {
            bool IsGood=true;
            for (uint i=0; i<TempAtoms.size()&&IsGood; i++)
               if ((*(*AtomI))==(*TempAtoms[i])) IsGood=false;
            if (!IsGood) continue;
            (*(*FragI))<<GhostAtom((*(*AtomI)));
         }
      }
   }
}

void VMFCnItr::UpdateMe(uint NewGhost) {
   Mine_=Ghost_;
   //std::cout<<"Calling UpdateMe"<<std::endl;
   SerialNumber::const_iterator SNI=Real_.begin();
   for (uint counter=1; counter<=Real_.size(); counter++,++SNI)
      Mine_.insert((*SNI)*(counter==NewGhost?-1.0:1.0));
   Current_=Mine_;
   if(Real_.size()>2)
      MyItr_=boost::shared_ptr<VMFCnItr>(
         new VMFCnItr(Mine_,AllowDuplicates_));
}

VMFCnItr::VMFCnItr(const SerialNumber& SN,const bool AllowDuplicates) :
      IsDone_(false),AllowDuplicates_(AllowDuplicates), l_(0) {
   SerialNumber::const_iterator SNI=SN.begin(),SNIEnd=SN.end();
   for (; SNI!=SNIEnd; ++SNI) {
      if ((*SNI)>0) Real_.insert(*SNI);
      else Ghost_.insert(*SNI);
   }
   l_=Real_.size();
   UpdateMe(l_--);
}

void VMFCnItr::Next() {
   do{
     // std::cout<<"In Mine="<<Mine_.PrintOut()<<std::endl;
      if(Current_==Mine_){
         if(MyItr_){//(m=/=1)-mers in n-mer basis
            Current_=(*(*MyItr_));
            ++(*MyItr_);
            //std::cout<<"Current set off initalized iterator"<<std::endl;
         }
         else if(l_!=0)UpdateMe(l_--);//Monomers in n-mer basis
         else IsDone_=true;
      }
      else{//I have a subiterator and it's been running
         if(!MyItr_->Done()){//Is that iterator done?
            //std::cout<<"Current set off iterator"<<std::endl;
            Current_=(*(*MyItr_));
            ++(*MyItr_);
         }
         else if(l_!=0)UpdateMe(l_--);//Can I be updated
         else IsDone_=true;//I'm done
      }
   }while(!IsDone_&&!AllowDuplicates_&&FoundSNs_.count(Current_)==1);
   FoundSNs_.insert(Current_);
}

void VMFCn::CalcBSSE(NMers& Sys, uint Stop, uint Start) const {
   if (Start!=1&&Stop!=0)
      throw PSIEXCEPTION("VMFCn Must be performed on the full set of n-mers");
   uint N=(Stop==0 ? Sys.N() : Stop);
   //This will be the result
   NMers TempNMers(Sys);

   //Create a mapping between a Monomer's SN and it's actual self for ease
   PsiMap <long int,boost::shared_ptr<Fragment> > Monomers;
   NMers::NMerItr_t NMerI=Sys.NMerBegin(1),NMerIEnd=Sys.NMerEnd(1);
   for(;NMerI!=NMerIEnd;++NMerI){
      SerialNumber SN=(*NMerI)->GetSN();
      Monomers[(*SN.begin())]=(*NMerI);
   }

   //There are no VMFC(n) corrections for n==1
   for (uint n=Start+1; n<=N; n++) {
      NMerI=Sys.NMerBegin(n);
      NMerIEnd=Sys.NMerEnd(n);
      for (; NMerI!=NMerIEnd; ++NMerI) {
         SerialNumber SN=(*NMerI)->GetSN();
         VMFCnItr SNI(SN,false);
         for (; !SNI.Done(); ++SNI) {
            SerialNumber::const_iterator SNJ=SNI->begin(),
                                        SNJEnd=SNI->end();
            --SNJEnd;
            boost::shared_ptr<Fragment> Temp(
                  new Fragment(Monomers[(*SNJEnd)]->Mol(),(*SNI))
            );
            (*Temp)+=(*Monomers[(*SNJEnd)]);
            for(;SNJ!=SNJEnd;++SNJ){
               long int FragN=(*SNJ)*((*SNJ)>0?1:-1);
               if((*SNJ)>0)(*Temp)+=(*Monomers[FragN]);
               else{
                  MolItr AtomI=Monomers[FragN]->Begin(),
                      AtomIEnd=Monomers[FragN]->End();
                  for(;AtomI!=AtomIEnd;++AtomI)
                     (*Temp)<<GhostAtom(*(*AtomI));
               }
            }
            Temp->GetSN().PrintOut();
            TempNMers.AddNMer(n,Temp);
         }
      }
   }
   Sys=TempNMers;
}
}
}         //End namespaces

