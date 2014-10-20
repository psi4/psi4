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
#include "libparallel/TableSpecs.h"
namespace psi{
namespace LibFrag{

Connections::Connections(GroupType& Groups,
      const Connections& AtomConnec):
      NConnecs_(Groups.size(),0),CTable_(nbonds_*Groups.size(),0),
      Distance_(Groups.size(),Groups.size()),
      bondthresh_(AtomConnec.bondthresh_){
  for(int i=0;i<Groups.size();i++){
     Distance_(i,i)=0.0;
     for(int j=i+1;j<Groups.size();j++){
        Distance_(i,j)=Groups[i]->Distance((*Groups[j]));
        Distance_(j,i)=Distance_(i,j);
        bool bonded=false;
        for(int k=0;k<Groups[i]->size()&&!bonded;k++){
           int atom=(*Groups[i])[k];
           for(int l=0;l<AtomConnec.GetNConnecs(atom)&&!bonded;l++){
              int atom2=AtomConnec(atom,l)-1;
              bonded=Groups[j]->Contains(atom2);
           }
        }
        if(bonded){
           CTable_[i*MaxBonds()+GetNConnecs(i)]=j+1;
           CTable_[j*MaxBonds()+GetNConnecs(j)]=i+1;
           NConnecs_[i]++;NConnecs_[j]++;
        }

     }
  }
}


void Connections::print_out()const{
   std::vector<std::string> Titles;
   Titles.push_back("Atom/Group #");
   int natoms=NConnecs_.size();
   std::vector<boost::shared_ptr<int[]> > Bonds;
   for(int i=0;i<MaxBonds();i++){
      std::stringstream title;
      title<<"Bond # "<<i;
      Titles.push_back(title.str());
      boost::shared_ptr<int[]> temp(new int[natoms]);
      Bonds.push_back(temp);
   }
   std::vector<int> Atom;
   for(int i=0;i<natoms;i++){
      Atom.push_back(i);
      for(int j=0;j<MaxBonds();j++){
         Bonds[j][i]=((*this)(i,j));
      }
   }
   TableSpecs<int,int,int,int,int> Specs(natoms);
   Specs.Init(&Atom[0],&(Bonds[0][0]),&(Bonds[1][0]),&(Bonds[2][0]),
         &(Bonds[3][0]));
   Specs.SetTitles(Titles);
   (*outfile)<<Specs.Table();
   (*outfile)<<std::endl;
}

Connections::Connections(SharedMol& Mol):Distance_(Mol->distance_matrix()),
NConnecs_(Mol->natom(),0),CTable_(nbonds_*(Mol->natom()),0),
bondthresh_(2.0/pc_bohr2angstroms){
   RoughTable(Mol);
}

void Connections::RoughTable(SharedMol Mol){
   for(int i=0;i<Distance_.rowdim();i++){
      for(int j=i+1;j<Distance_.coldim();j++){
         //allow bonds to be up to 10% longer
         if(Mol)bondthresh_=1.1*(cov_radii[(int)Mol->Z(i)] +
               cov_radii[(int)Mol->Z(j)])/pc_bohr2angstroms;
         if(Distance_(i,j)<=bondthresh_){
            NConnecs_[i]++;
            NConnecs_[j]++;
            if(NConnecs_[i]>nbonds_||NConnecs_[j]>nbonds_)
              throw psi::PSIEXCEPTION("Hard-coded number of bonds exceeded");
            CTable_[i*nbonds_+(NConnecs_[i]-1)]=j+1;
            CTable_[j*nbonds_+(NConnecs_[j]-1)]=i+1;
         }
      }
   }
}




}}//End namespaces


