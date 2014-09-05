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

#include "BSSEer.h"
#include "MBEFrag.h"

namespace psi{
namespace LibFrag{


void BSSEer::CommonInit(GhostType& Ghosts,NMerSet& NMers){
   Set<SharedGhost> TempSet;
   Ghosts.push_back(TempSet);
   if(Ghosts.size()==1){
      for(int i=0;i<NAtoms_;i++){
         SharedAtom Atom=NMers[0]->Atoms_.Object(i);
         SharedGhost temp(new Ghost(i,Atom->Z(),Atom->Carts()[0],Atom->Carts()[1],
               Atom->Carts()[2]));
      }
   }
   else{
      Ghosts.back()=Ghosts[0];
      Ghosts.back().Clear();
   }
}

void BSSEer::AddBSSEJobs(NMerSet& NMers){
   GhostType Ghosts;
   BSSEImpl(Ghosts,NMers);
   if(Ghosts.size()!=NMers.size()&&Ghosts.size()!=0)throw PSIEXCEPTION(
     "Each fragment must have a ghost atom set, or no fragment may have one");
   for(int i=0;i<Ghosts.size();i++)
      NMers[i]->Ghosts_=Ghosts[i];
}

void FullBSSE::BSSEImpl(GhostType& Ghosts,NMerSet& NMers){
   for(int NMer=0;NMer<NMers.size();NMer++){
      std::vector<bool>Is_Real(NAtoms_,false);
      SharedFrag DaNMer=NMers[NMer];
      CommonInit(Ghosts,NMers);
      for(int atoms=0;atoms<DaNMer->Atoms().size();atoms++)
         Is_Real[(*DaNMer).Atoms()[atoms]]=true;
      for(int caps=0;caps<DaNMer->Caps().size();caps++){
         int ActualCap=DaNMer->Caps()[caps];
         Is_Real[DaNMer->Caps().Object(ActualCap)->ReplacedAtom()]=true;
      }
      for(int atoms=0;atoms<NAtoms_;atoms++){
         if(!Is_Real[atoms])Ghosts.back()<<atoms;
      }
   }
}
}}//End namespaces

