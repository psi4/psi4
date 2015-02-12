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
#include <set>
#include <algorithm>
#include <iterator>
#include "UnitCell.h"
#include "SuperCell.h"
#include "../Utils/GeomManipulator.h"
#include "../Utils/Geometry.h"
namespace psi{
namespace LibMolecule{
void MolRecurse(const Connections& Conns, std::set<int>& Atoms2Move,
                int Last,int index){
   Atoms2Move.insert(index);
   int NConns=Conns[index].size();
   for(int i=0;i<NConns;i++){
      int AtomJ=Conns[index][i];
      if(AtomJ==Last)continue;
      Atoms2Move.insert(AtomJ);
      MolRecurse(Conns,Atoms2Move,index,AtomJ);
  }
}

Molecule MolFromSet(const std::set<int>& Atoms,const Molecule& Other){
   Molecule TempMol;
   std::set<int>::const_iterator AtomI=Atoms.begin(),AtomIEnd=Atoms.end();
   for(;AtomI!=AtomIEnd;++AtomI)TempMol<<(*Other[(*AtomI)]);
   return TempMol;
}

std::set<int> SetDiff(const std::set<int>& One,int NAtoms){
   std::set<int> TempSet,Diff;
   for(int i=0;i<NAtoms;i++)TempSet.insert(i);
   std::set_difference(TempSet.begin(),TempSet.end(),
          One.begin(),One.end(),
          std::insert_iterator<std::set<int> >(Diff,Diff.end()));
   return Diff;
}

//1) Make a 3x3x3 supercell
//2) Grab all atoms in the unitcell, and their new bonds
//3) Make new 3x3x3 supercell
//4) Remove redundencies
void UnitCell::FixUnitCell(){
   std::vector<int> CellVecs(3,1);
   UnitCell TempUC(*this);
   SuperCell TempSC(TempUC,CellVecs);
   Geometry Geom(&TempSC);
   MolItr AtomI=this->Begin(),AtomIEnd=this->End();
   std::set<int> Atoms;
   for(int index=0;AtomI!=AtomIEnd;++AtomI,++index){
      int Last=-1;
      MolRecurse(Geom.GetConns(),Atoms,Last,index);
   }
   /*Now at this point we have some redundancies
    * Imagine:
    *   _____________________
    *  |      |      |      |
    *  | Y  X + Y  X + Y  X |
    *  |      |      |      |
    *  ----------------------
    *
    * If the central cell is our unit cell, simply replicating it will
    * make an X-Y molecule sit on top  of another X-Y molecule.
    */
   Molecule TempMol=MolFromSet(Atoms,TempSC);
   while(true){
      SuperCell TempSC2(UnitCell(TempMol,TempUC.Sides(),
            TempUC.Angles(),false,true,false)
            ,CellVecs);
      MolItr AtomK=TempSC2.Begin(),AtomKEnd=TempSC2.End();
      std::set<int> Atoms2Remove;
      for(int indexi=0;indexi<TempMol.NAtoms();++AtomK,++indexi){
         MolItr AtomL=AtomK;++AtomL;
         for(;AtomL!=AtomKEnd;++AtomL){
            double dz=0.0;
            for(int i=0;i<3;i++){
               double dq=(*(*AtomL))[i]-(*(*AtomK))[i];
               dz+=dq*dq;
            }
            if(sqrt(dz)<1e-6){
               Geometry Geom2(&TempMol);
               MolRecurse(Geom2.GetConns(),Atoms2Remove,-1,indexi);
               break;
            }
         }
         if(Atoms2Remove.size()>0)break;
      }
      if(Atoms2Remove.size()==0)break;
      std::set<int> Diff;
      Atoms=SetDiff(Atoms2Remove,TempMol.NAtoms());
      TempMol=MolFromSet(Atoms,TempMol);
   }
   UnitCell TempUC2(TempMol,TempUC.Sides(),
            TempUC.Angles(),false,true,false);
   (*this)=TempUC2;
}

}}//End namespaces

