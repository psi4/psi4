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
#include "MBEFrag.h"
#include "Embedder.h"
#include "psi4-dec.h"
#include <math.h>
namespace psi{
namespace LibFrag{

void Embedder::print_out(){
   for(int i=0;i<Charges_.size();i++){
      psi::outfile->Printf("%16.15f ",Charges_[i]);
      for(int j=0;j<3;j++){
         psi::outfile->Printf("%16.15f ",Carts_[i*3+j]);
      }
      psi::outfile->Printf("\n");
   }
   double corr=0.0;
   for(int i=0;i<Charges_.size();i++){
      for(int j=i+1;j<Charges_.size();j++){
         double dr=0.0;
         for(int k=0;k<3;k++)dr+=pow(Carts_[i*3+k]-Carts_[j*3+k],2);
         corr+=Charges_[i]*Charges_[j]/sqrt(dr);
      }
   }
   psi::outfile->Printf("Self-interaction correction is: %16.15f\n",corr);
}

bool Embedder::Iterate(const int itr){
   return ((itr==0)?true : DoIterate_);
}

void Embedder::SetCharge(int i,double q){
   NotFirstItr_=true;
   Charges_[i]=q;
}

Embedder::Embedder(SharedMol& AMol,bool Iterating):
      DoIterate_(Iterating),Charges_(AMol->natom(),0.0),
      NotFirstItr_(false){
   for(int i=0;i<AMol->natom();i++){
      Carts_.push_back(AMol->fx(i));
      Carts_.push_back(AMol->fy(i));
      Carts_.push_back(AMol->fz(i));
   }
}

void Embedder::Embed(NMerSet& Set2Embed){
   ChargeType ChargesBySet;
   EmbedImpl(ChargesBySet,Set2Embed);
   if(ChargesBySet.size()!=Set2Embed.size())
      throw PSIEXCEPTION("Each fragment must have an embedding set even"
            " if it is just empty.");
   for(int i=0;i<ChargesBySet.size();i++){
      Set2Embed[i]->Charges_=ChargesBySet[i];
   }
}

void Embedder::CommonInit(ChargeType& ChargesBySet){
   Set<SharedCharge> newcharge;
   ChargesBySet.push_back(newcharge);
   if(ChargesBySet.size()==1){
     for(int charge=0;charge<Charges_.size();charge++){
        SharedCharge temp=SharedCharge(new Charge(
                 Charges_[charge],charge,&Carts_[charge*3]));
      ChargesBySet[0].AddSet(charge,temp);
     }
   }
   else{
      ChargesBySet.back()=ChargesBySet[0];
      ChargesBySet.back().Clear();
   }
}

void APCEmbedder::EmbedImpl(ChargeType& ChargesBySet, NMerSet& Set2Embed){
   if(Set2Embed.size()==0)
      throw PSIEXCEPTION("Expected fragments to embed");
   for(int set=0;set<Set2Embed.size();set++){
      CommonInit(ChargesBySet);
      std::vector<bool> IsCharge(Charges_.size(),true);
      SharedFrag NMer=Set2Embed[set];
      for(int atom=0;atom<NMer->Atoms().size();atom++)
         IsCharge[NMer->Atoms()[atom]]=false;
      for(int cap=0;cap<NMer->Caps().size();cap++){
         int ActualCap=NMer->Caps()[cap];
         IsCharge[NMer->Caps().Object(ActualCap)->ReplacedAtom()]=false;
      }
      for(int i=0;i<IsCharge.size();i++)
         if(IsCharge[i])ChargesBySet.back()<<i;

   }
}

}}//End namespaces

