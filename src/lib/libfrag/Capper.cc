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

namespace LibFrag {

Capper::Capper(SharedMol& AMol) {
   AtomConnect=boost::shared_ptr<Connections>(new Connections(AMol));
   for(int i=0;i<AMol->natom();i++){
      AtomCarts.push_back(AMol->fx(i));
      AtomCarts.push_back(AMol->fy(i));
      AtomCarts.push_back(AMol->fz(i));
   }
}

std::vector<int> Capper::Atoms2Replace(NMerSet& Set2Cap)const{
   std::vector<int> Atoms;
   for(int i=0;i<Set2Cap.size();i++){
      boost::shared_ptr<MBEFrag> NMer=Set2Cap[i];
      for(int j=0;j<NMer->size();j++){
         int atom=(*NMer)[j];
         bool severed=false;
         for(int k=0;k<AtomConnect->GetNConnecs(atom)&&!severed;k++){
            int ConAtom=(*AtomConnect)(atom,k)-1;
            if(!NMer->contains(ConAtom)){
               Atoms.push_back(i);
               Atoms.push_back(atom);
               Atoms.push_back(ConAtom);
            }
         }
      }
   }
   return Atoms;
}

void ReplaceAndCap::MakeCaps(NMerSet& Set2Cap) {
   std::vector<int> Bonds=Atoms2Replace(Set2Cap);
   for(int i=0;i<Bonds.size()/3;i++){
      int frag=Bonds[i*3], atom=Bonds[i*3+2];
      Set2Cap[frag]->Caps.push_back(Cap("H",AtomCarts[atom*3],
            AtomCarts[atom*3+1],AtomCarts[atom*3+2],atom));
   }
}
ShiftAndCap::ShiftAndCap(SharedMol& AMol):Capper(AMol){
   for(int i=0;i<AMol->natom();i++)
      Zs.push_back(AMol->Z(i));
}

void ShiftAndCap::MakeCaps(NMerSet& Set2Cap){
   std::vector<int> Bonds=Atoms2Replace(Set2Cap);
   for(int i=0;i<Bonds.size()/3;i++){
      int frag=Bonds[i*3], atom=Bonds[i*3+1], cap=Bonds[i*3+2];

      //This is the ratio of the atom-hydrogen to atom-cap bond
      double ratio=(cov_radii[Zs[atom]]+cov_radii[1])/
            (cov_radii[Zs[atom]]+cov_radii[Zs[cap]]);

      std::vector<double> newcarts(3);
      for(int i=0;i<3;i++){
         newcarts[i]=AtomCarts[atom*3+i]+
               ratio*(AtomCarts[cap*3+i]-AtomCarts[atom*3+i]);
      }
      Set2Cap[frag]->Caps.push_back(Cap("H",newcarts[0],newcarts[1],
            newcarts[2],cap));
   }
}

}

