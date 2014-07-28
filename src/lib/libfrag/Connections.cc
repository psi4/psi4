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

#include "Connections.h"
#include "libmints/molecule.h"
#include "exception.h"
#include "psi4-dec.h"
#include "physconst.h"
#include "cov_radii.h"
#include <algorithm>

namespace LibFrag{

Connections::Connections(std::vector<AtomSet>& Groups,
      const Connections& AtomConnec):
      NConnecs(Groups.size(),0),CTable(nbonds*Groups.size(),0),
      Distance(Groups.size(),Groups.size()),
      bondthresh(AtomConnec.bondthresh){
  for(int i=0;i<Groups.size();i++){
     Distance(i,i)=0.0;
     for(int j=i+1;j<Groups.size();j++){
        Distance(i,j)=Groups[i].Distance(Groups[j]);
        Distance(j,i)=Distance(i,j);
        bool bonded=false;
        for(int k=0;k<Groups[i].size()&&!bonded;k++){
           int atom=(Groups[i])[k];
           for(int l=0;l<AtomConnec.GetNConnecs(atom)&&!bonded;l++){
              int atom2=AtomConnec(atom,l)-1;
              bonded=Groups[j].contains(atom2);
           }
        }
        if(bonded){
           CTable[i*nbonds+NConnecs[i]]=j+1;
           CTable[j*nbonds+NConnecs[j]]=i+1;
           NConnecs[i]++;NConnecs[j]++;
        }

     }
  }
}


void Connections::print_out(){
   fprintf(psi::outfile,"Connections by atom:");
   for(int i=0;i<NConnecs.size();i++){
      fprintf(psi::outfile,"\n%d 's connections: ",i+1);
      for(int j=0;j<NConnecs[i];j++){
         fprintf(psi::outfile," %d",CTable[i*nbonds+j]);
      }
   }
   fprintf(psi::outfile,"\n");
}

Connections::Connections(SharedMol& Mol):Distance(Mol->distance_matrix()),
NConnecs(Mol->natom(),0),CTable(nbonds*Mol->natom(),0),
bondthresh(2.0/pc_bohr2angstroms){
   RoughTable(Mol);
}

void Connections::RoughTable(SharedMol Mol){
   for(int i=0;i<Distance.rowdim();i++){
      for(int j=i+1;j<Distance.coldim();j++){
         //allow bonds to be up to 10% longer
         if(Mol)bondthresh=1.1*(cov_radii[(int)Mol->Z(i)] +
               cov_radii[(int)Mol->Z(j)])/pc_bohr2angstroms;
         if(Distance(i,j)<=bondthresh){
            NConnecs[i]++;
            NConnecs[j]++;
            if(NConnecs[i]>nbonds||NConnecs[j]>nbonds)
              throw psi::PSIEXCEPTION("Hard-coded number of bonds exceeded");
            CTable[i*nbonds+(NConnecs[i]-1)]=j+1;
            CTable[j*nbonds+(NConnecs[j]-1)]=i+1;
         }
      }
   }
}




}//End libfrag namespace


