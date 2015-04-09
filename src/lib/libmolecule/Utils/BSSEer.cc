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
#include <algorithm>
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

VMFCnItr::VMFCnItr(const SerialNumber& SN){
   SerialNumber::const_iterator It=SN.begin(),ItEnd=SN.end();
   for(;It!=ItEnd;++It){
      if((*It)<0)InitialGhosts_.insert(*It);
      else SN_.insert(*It);
   }
   It_=boost::shared_ptr<PowerSetItr<SerialNumber> >(
         new PowerSetItr<SerialNumber>(SN_));
   //Power set iterator starts on empty set
   Next();
}

void VMFCnItr::Next() {
   ++(*It_);
   if(!It_->Done()){
      SerialNumber Temp;
      std::set_difference(SN_.begin(),SN_.end(),
            (*It_)->begin(),(*It_)->end(),
            std::inserter(Temp,Temp.begin()));
      Value_=(*(*It_));
      SerialNumber::const_iterator It1=Temp.begin(),It1End=Temp.end();
      for(;It1!=It1End;++It1)
         Value_.insert(-1*(*It1));
      //Don't want the original SN
      if((*Value_.begin())>=0)++(*It_);
      //Need to add the ghosts now or else they mess w/ the last check
      Value_.insert(InitialGhosts_.begin(),InitialGhosts_.end());
   }

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
         VMFCnItr SNI(SN);
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
            TempNMers.AddNMer(n,Temp);
         }
      }
   }
   Sys=TempNMers;
}
}
}         //End namespaces

