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

#include "Grouper.h"
#include "libmints/molecule.h"
namespace psi{
namespace LibFrag{

inline void MakeSharedAtom(const SharedMol& AMol, CartSet<SharedAtom>& Atoms) {
   for (int atom=0; atom<AMol->natom(); atom++) {
      int Z=AMol->Z(atom);
      double m=AMol->fmass(atom),x=AMol->fx(atom),y=AMol->fy(atom);
      double z=AMol->fz(atom);
      SharedAtom tempAtom(new Atom(Z, m, x, y, z));
      Atoms.AddSet(atom, tempAtom);
   }
}

inline void AddAtomToGroup(std::vector<bool>& Assigned, GroupType Groups,
      int atom) {
   if (Assigned[atom])
      throw PSIEXCEPTION("Atom is already assigned to a group");
   (*(Groups.back()))<<atom;
   Assigned[atom]=true;
}

inline void CommonInit(const SharedMol& AMol, GroupType& Frags) {
   Frags.push_back(SharedAtomSet(new CartSet<SharedAtom>()));
   if (Frags.size()==1) MakeSharedAtom(AMol, *Frags[0]);
   else { //Copy existing universe
      (*(Frags.back()))=(*(Frags[0]));
      //We also got Frags[0]'s atoms, get rid of them
      Frags.back()->Clear();
   }
}

void Grouper::CondenseHydrogens(GroupType& Groups,
      std::vector<bool>& Assigned,
      const Connections& CTable,const SharedMol& Mol2Frag)const{
   for (int i=0; i<Assigned.size(); i++) {
      if (Mol2Frag->Z(i)==1&&!Assigned[i]) {
         if (CTable.GetNConnecs(i)!=1)
            throw PSIEXCEPTION(
                  "This hydrogen is not connected to one and only one atom");
         CommonInit(Mol2Frag, Groups);
         AddAtomToGroup(Assigned, Groups, i);
         int center=CTable(i, 0)-1;
         AddAtomToGroup(Assigned, Groups, center);
         for (int j=0; j<CTable.GetNConnecs(center); j++) {
            int atom=CTable(center, j)-1;
            //We already know it's attached to i
            if (atom!=i&&Mol2Frag->Z(atom)==1) {
               if (CTable.GetNConnecs(atom)!=1) throw psi::PSIEXCEPTION(
                     "This hydrogen is not connected to only one atom");
               AddAtomToGroup(Assigned, Groups, atom);
            }
         }
         (*(Groups.back())).Sort();
      }
   }
}

bool Grouper::CheckRing(const std::vector<int>& Path,GroupType& Groups,
      const int CurrGroup)const {
   bool ringnotfound=true;
   //Path[Path.size()-1] is CurrGroup so NewGroup=/=Path.back()
   //Path[Path.size()-2] was already checked
   for (int j=0; j<Path.size()-2&&ringnotfound; j++) {
      if (Path[j]==CurrGroup) {
         std::vector<bool> GoodGroups(Groups.size(), true);
         ringnotfound=false;
         //We have a Path.size()-j membered ring
         //Involving groups Path[j] to Path[Path.size()-1]
         //add the atoms in those groups to Group[j]
         for (int k=j+1; k<Path.size(); k++) {
            int Groupk=Path[k];
            for (int l=0; l<Groups[Groupk]->size(); l++) {
               (*Groups[CurrGroup])<<(*Groups[Groupk])[l];
            }
            GoodGroups[Groupk]=false;
         }
         GroupType TempGroups;
         for (int k=0; k<Groups.size(); k++) {
            if (GoodGroups[k]) TempGroups.push_back(Groups[k]);
         }
         Groups=TempGroups;
      }
   }
   return ringnotfound;
}

void Grouper::RingRecursion(const int MaxSize,const Connections& GroupTable,
      GroupType& Groups, std::vector<int>& Path, bool& ringnotfound)const {
   int depth=Path.size()-1;
   if (depth<MaxSize &&ringnotfound) {
      int CurrGroup=Path.back();
      int NConns=GroupTable.GetNConnecs(CurrGroup);
      for (int i=0; i<NConns && ringnotfound; i++) {
         int NewGroup=GroupTable(CurrGroup, i)-1;
         bool Good=(depth>1 ? NewGroup!=Path[depth-1] : true);
         if (Good) {
            if (Path.size()>=3)//Need at least 3 groups for a ring
               ringnotfound=CheckRing(Path,Groups,NewGroup);
            if(ringnotfound){
               Path.push_back(NewGroup);
               RingRecursion(MaxSize,GroupTable,Groups,Path,ringnotfound);
               Path.pop_back();
            }
         }
      }
   }
}

void Grouper::FindRings(const int MaxSize,GroupType& Groups,
      const Connections& CTable)const {
   bool done;
   Connections tempGroupTable(Groups,CTable);
   do{
      done=true;
      Connections GroupTable(Groups, CTable);
      for(int i=0;i<Groups.size()&&done;i++){
         std::vector<int> Path(1,i);
         RingRecursion(MaxSize,GroupTable,Groups,Path,done);
      }
   }while(!done);

}


GroupType Grouper::MakeGroups(const Connections& CTable,
      const SharedMol& Mol2Frag)const{
   GroupType Groups;

      std::vector<bool> Assigned(Mol2Frag->natom(), false);

      CondenseHydrogens(Groups, Assigned, CTable, Mol2Frag);

      //Here we should now do things like make carbonyls, carboxyls, etc.
      //into groups

      //Whatever is left is it's own group
      for (int i=0; i<Assigned.size(); i++) {
         if (!Assigned[i]) {
            CommonInit(Mol2Frag, Groups);
            AddAtomToGroup(Assigned, Groups, i);
         }
      }
      //With all atoms assigned to groups, our last check looks for
      //rings...The plan is if it's a six-membered (just for you benzene...)
      //or less ring make it a group.  For fused rings they should collapse
      //down to one group
      FindRings(6,Groups,CTable);

      return Groups;
}


}}//end namespaces

