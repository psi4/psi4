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

#include "Capper.h"
#include "libmints/molecule.h"
#include "Connections.h"
#include "MBEFrag.h"
#include "cov_radii.h"
#include "masses.h"
namespace psi{
namespace LibFrag {

Capper::Capper(SharedMol& AMol) {
   AtomConnect_=boost::shared_ptr<Connections>(new Connections(AMol));
   for(int i=0;i<AMol->natom();i++){
      AtomCarts_.push_back(AMol->fx(i));
      AtomCarts_.push_back(AMol->fy(i));
      AtomCarts_.push_back(AMol->fz(i));
   }
}

void Capper::MakeCaps(NMerSet& Set2Cap){
   CapType Caps;
   CapImpl(Caps,Set2Cap);
   if(Caps.size()!=Set2Cap.size())
      throw PSIEXCEPTION("Each set must have a cap set even if it is just "
            "an empty set.");
   for(int i=0;i<Caps.size();i++)
      Set2Cap[i]->Caps_=Caps[i];
}

std::vector<int> Capper::Atoms2Replace(NMerSet& Set2Cap)const{
   std::vector<int> Atoms;
   for(int i=0;i<Set2Cap.size();i++){
      boost::shared_ptr<MBEFrag> NMer=Set2Cap[i];
      for(int j=0;j<NMer->Atoms().size();j++){
         int atom=NMer->Atoms()[j];
         bool severed=false;
         for(int k=0;k<AtomConnect_->GetNConnecs(atom)&&!severed;k++){
            int ConAtom=(*AtomConnect_)(atom,k)-1;
            if(!NMer->Atoms().Contains(ConAtom)){
               Atoms.push_back(i);
               Atoms.push_back(atom);
               Atoms.push_back(ConAtom);
            }
         }
      }
   }
   return Atoms;
}

void UnofficialInit(CapType& Caps,NMerSet& Set2Cap){
   CartSet<SharedCap> TempSet1;
     Caps.push_back(TempSet1);
     for(int i=1;i<Set2Cap.size();i++){
        CartSet<SharedCap> TempSet(Caps[0]);
        TempSet.Clear();
        Caps.push_back(TempSet);
     }
}

void FillCap(int* Bonds,int i, CapType& Caps,double* AtomCarts){
   int frag=Bonds[0],BA=Bonds[1],RA=Bonds[2];
   SharedCap TempCap(
               new Cap(BA,RA,1,an2masses[1],AtomCarts));
         Caps[0].AddSet(i,TempCap);
         Caps[frag]<<i;
}

void ReplaceAndCap::CapImpl(CapType& Caps,NMerSet& Set2Cap) {
   std::vector<int> Bonds=Atoms2Replace(Set2Cap);
   UnofficialInit(Caps,Set2Cap);
   for(int i=0;i<Bonds.size()/3;i++)
      FillCap(&Bonds[i*3],i,Caps,&AtomCarts_[Bonds[i*3+2]*3]);
}
ShiftAndCap::ShiftAndCap(SharedMol& AMol):Capper(AMol){
   for(int i=0;i<AMol->natom();i++)
      Zs_.push_back(AMol->Z(i));
}

void ShiftAndCap::CapImpl(CapType& Caps, NMerSet& Set2Cap){
   std::vector<int> Bonds=Atoms2Replace(Set2Cap);
   UnofficialInit(Caps,Set2Cap);
   for(int i=0;i<Bonds.size()/3;i++){
      int atom=Bonds[i*3+1], cap=Bonds[i*3+2];
      //This is the ratio of the atom-hydrogen to atom-cap bond
      double ratio=(cov_radii[Zs_[atom]]+cov_radii[1])/
            (cov_radii[Zs_[atom]]+cov_radii[Zs_[cap]]);
      std::vector<double> newcarts(3);
      for(int j=0;j<3;j++){
         newcarts[j]=AtomCarts_[atom*3+j]+
               ratio*(AtomCarts_[cap*3+j]-AtomCarts_[atom*3+j]);
      }
      FillCap(&Bonds[i*3],i,Caps,&newcarts[0]);
   }
}

}}//End namespaces

