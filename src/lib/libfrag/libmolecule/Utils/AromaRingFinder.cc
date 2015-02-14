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

#include "AromaRingFinder.h"
#include "FxnalGroup.h"
#include "OrganicGeom.h"
#include <iostream>
#include <set>
namespace psi {
namespace LibMolecule {

typedef boost::shared_ptr<const FxnalGroup> SharedGroup;
typedef ConnGroups PsiMap_t;
typedef const Connections Conn_t;
typedef std::vector<std::pair<int, int> > Attach_t;
typedef std::vector<SharedGroup> Group_t;
std::vector<int> UpdateAttachPoints(Group_t& Groups,
             std::set<int>& Members,Conn_t& Conns,int &Order){
   Group_t::const_iterator GroupI=Groups.begin();
   std::vector<int> AttachPoints;
   Order=0;
   std::set<int>::iterator MemI=Members.begin();
   for(;GroupI!=Groups.end();++GroupI){
      SharedGroup GI=(*GroupI);
      int OrderI=GI->Order();
      for(int i=0;i<GI->NAttachPoint();i++){
         int AttachI=GI->AttachPoint(i);
         if(AttachI==-1)continue;
         std::vector<int>::const_iterator ConnI=Conns[AttachI].begin();
         bool AllUsed=true;
         for(;ConnI!=Conns[AttachI].end();++ConnI){
            if(Members.count((*ConnI))==1){
               bool good=true;
               //Make sure it's not part of GI
               for(int j=0;j<GI->size()&&good;j++)
                  if((*ConnI)==(*GI)[j])good=false;
               if(good)OrderI--;
            }
            else AllUsed=false;
         }
         if(!AllUsed)AttachPoints.push_back(AttachI);
      }
      Order+=OrderI;
   }
   return AttachPoints;
}

bool CheckForRing(std::vector<SharedGroup>& Groups, PsiMap_t& FoundGroups,
      Attach_t& Attach, Conn_t& Conns) {
   bool Found=false;
   SharedGroup GI=Groups.back();
   /* The logic for this loop is best explained by an example, a four
    * membered ring of atoms 1-2-3-4 would look like:
    * 123412 when it comes in here.  The last pair is the pair we are
    * checking against, and it can't possibly be equal to the pair
    * before it, hence the -2
    */
   for (int i=0; i<Attach.size()-2; i++) {
      //Not a ring
      if (Attach[i]!=Attach.back()) continue;
      //Have a ring
      Found=true;
      std::vector<SharedGroup> RGroups;
      std::vector<int> RAttach;
      std::set<int> Members;
      //Don't grab the last group twice (it's also i)!!!
      for (int j=i; j<Groups.size()-1; j++) {
         SharedGroup GJ=Groups[j];
         RGroups.push_back(GJ);
         for(int k=0;k<GJ->size();k++)
            Members.insert((*GJ)[k]);
         int NAttach=GJ->NAttachPoint();
         for (int k=0; k<NAttach; k++) {
            int AttachK=GJ->AttachPoint(k);
            if (FoundGroups.count(AttachK)==1) FoundGroups.erase(AttachK);
         }
      }
      int Order=0;
      RAttach=UpdateAttachPoints(RGroups,Members,Conns,Order);
      boost::shared_ptr<AromaticRing> temp(
            new AromaticRing(RAttach, Order,RGroups));
      for (int j=0; j<temp->NAttachPoint(); j++){
         if(temp->AttachPoint(j)==-1)continue;
         FoundGroups[temp->AttachPoint(j)]=temp;
      }
      //This line accounts for the case when the ring doesn't attach
      //to anything, we just stick it somewhere
      if (temp->NAttachPoint()==0) FoundGroups[Attach[i].first]=temp;
      break;
   }
   return Found;
}

bool Recurse(PsiMap_t& FoundGroups, Conn_t& Conns,
      std::vector<SharedGroup>& Groups, Attach_t& AttachPoints) {
   bool Found=false;
   //Stop looking after 5 groups
   if (Groups.size()==5) return Found;
   int AttachPoint=AttachPoints.back().second;
   SharedGroup GI=Groups.back();
   typedef std::vector<int>::const_iterator cIntIt;
   cIntIt ConnI=Conns[AttachPoint].begin(),ConnEnd=Conns[AttachPoint].end();
   for (; ConnI!=ConnEnd&&!Found; ++ConnI) {
      //Check if atom is attach-able
      if (FoundGroups.count((*ConnI))==0) continue;
      SharedGroup GJ=FoundGroups[(*ConnI)];
      if (GJ==GI) continue;
      if (Groups.size()>=2)//Ensure we didn't go backwards okay
         if ((*ConnI)==AttachPoints[AttachPoints.size()-2].second) continue;
      //For fused rings these keep our aromaticity going...
      bool SpecialTypes=(GJ->Type()==ALKENYL3||GJ->Type()==ALKENYL2);
      bool DoubleBonds=(GJ->Type()==CCDB4||GJ->Type()==CCDB3||GJ->Type()==CCDB2);
      bool NBonds=(GJ->Type()==ALDIMINE2||GJ->Type()==KETIMINE2);
      if (DoubleBonds||SpecialTypes||NBonds||GJ->Type()==AROMATICRING) {
         //An aromatic ring must have at least two open sites
         if (GJ->Type()==AROMATICRING&&GJ->NAttachPoint()<2) continue;
         Groups.push_back(GJ);
         AttachPoints.push_back(std::pair<int, int>(*ConnI, *ConnI));
         for (int j=0; j<GJ->NAttachPoint(); j++) {
            int ConnJ=GJ->AttachPoint(j);
            if(ConnJ==-1)continue;
            if (!SpecialTypes) AttachPoints.back().second=ConnJ;
            Found=CheckForRing(Groups, FoundGroups, AttachPoints, Conns);
            if (!Found) Found=Recurse(FoundGroups, Conns, Groups, AttachPoints);
            if (Found) break;
         }
         if (!Found) {
            AttachPoints.erase(AttachPoints.end()-1, AttachPoints.end());
            Groups.erase(Groups.end()-1, Groups.end());
         }
         else break;
      }
   }
   return Found;
}

template <FxnGroup_t T1, FxnGroup_t T2>
class FindRingFunctor {
   public:
      static bool FindRing(PsiMap_t& FoundGroups, Conn_t& Conns) {
         bool Found=false;
         GroupIt GroupI=FoundGroups.begin();
         for (; GroupI!=FoundGroups.end(); ++GroupI) {
            SharedGroup GI=(*GroupI);
            //Don't bother with groups that aren't attached to things
            if (GI->NAttachPoint()==0) continue;
            int AttachPoint=GI->AttachPoint();
            if(AttachPoint==-1)continue;
            if (GI->Type()==T1||GI->Type()==T2) {
               std::vector<SharedGroup> Groups;
               Groups.push_back(GI);
               Attach_t AttachPoints;
               AttachPoints.push_back(
                     std::pair<int, int>(AttachPoint, AttachPoint));
               for (int i=1; i<GI->NAttachPoint()&&!Found; i++) {
                  AttachPoints.back().second=GI->AttachPoint(i);
                  Found=Recurse(FoundGroups, Conns, Groups, AttachPoints);
               }
            }
            if (Found) break;
         }
         return Found;
      }
};

bool AromaticRingFinder::FindRing(PsiMap_t& FoundGroups, Conn_t& Conns) {
   //Start with rings on the edges (involve a secondary 2x bond)
   bool Found=FindRingFunctor<CCDB2, ALDIMINE2>::FindRing(FoundGroups, Conns);
   if (!Found)         //Now try to find ones with ternary 2x bonds
      Found=FindRingFunctor<CCDB3, KETIMINE2>::FindRing(FoundGroups, Conns);
   if (!Found)         //The only other types that can be part of rings...
      Found=FindRingFunctor<CCDB4, AROMATICRING>::FindRing(FoundGroups, Conns);
   return Found;
}

}
}         //End namespaces

