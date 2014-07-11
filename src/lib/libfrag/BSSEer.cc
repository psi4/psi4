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

namespace LibFrag{

void FullBSSE::AddBSSEJobs(NMerSet& NMers){
   for(int NMer=0;NMer<NMers.size();NMer++){
      std::vector<bool>Is_Real(NAtoms,false);
      SharedFrag DaNMer=NMers[NMer];
      //If n-mer was made by copying another n-mer as a starting point, then
      //it also got it's list of ghosts
      DaNMer->Ghosts.clear();
      for(int atoms=0;atoms<DaNMer->size();atoms++)
         Is_Real[(*DaNMer)[atoms]]=true;
      for(int atoms=0;atoms<NAtoms;atoms++){
         if(!Is_Real[atoms])DaNMer->AddGhost(atoms);
      }
   }
}
}

