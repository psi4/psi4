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
namespace LibFrag{

void Embedder::print_out(){
   for(int i=0;i<Charges.size();i++){
      psi::outfile->Printf("%16.15f ",Charges[i]);
      for(int j=0;j<3;j++){
         psi::outfile->Printf("%16.15f ",Carts[i*3+j]);
      }
      psi::outfile->Printf("\n");
   }
   double corr=0.0;
   for(int i=0;i<Charges.size();i++){
      for(int j=i+1;j<Charges.size();j++){
         double dr=0.0;
         for(int k=0;k<3;k++)dr+=pow(Carts[i*3+k]-Carts[j*3+k],2);
         corr+=Charges[i]*Charges[j]/sqrt(dr);
      }
   }
   psi::outfile->Printf("Self-interaction correction is: %16.15f\n",corr);
}

bool Embedder::Iterate(const int itr){
   return ((itr==0)?true : DoIterate);
}

void Embedder::SetCharge(int i,double q){
   NotFirstItr=true;
   Charges[i]=q;
}

Embedder::Embedder(SharedMol& AMol,bool Iterating):
      DoIterate(Iterating),Charges(AMol->natom(),0.0),NotFirstItr(false){
   for(int i=0;i<AMol->natom();i++){
      Carts.push_back(AMol->fx(i));
      Carts.push_back(AMol->fy(i));
      Carts.push_back(AMol->fz(i));
   }
}

void APCEmbedder::Embed(NMerSet& Set2Embed){
   for(int set=0;set<Set2Embed.size();set++){
      std::vector<bool> IsCharge(Charges.size(),true);
      SharedFrag NMer=Set2Embed[set];
      NMer->Charges.clear();
      for(int atom=0;atom<NMer->size();atom++)
         IsCharge[(*NMer)[atom]]=false;
      for(int cap=0;cap<NMer->Caps.size();cap++)
         IsCharge[NMer->Caps[cap].ReplacedAtom()]=false;
      for(int i=0;i<IsCharge.size();i++){
         if(IsCharge[i]){
            NMer->Charges.push_back(Charges[i]);
            for(int j=0;j<3;j++)NMer->Charges.push_back(Carts[i*3+j]);
         }
      }
   }
}

}//End namespace LibFrag

