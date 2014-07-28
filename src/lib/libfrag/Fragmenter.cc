/*
 * Fragmenter.cc
 *
 *  Created on: May 15, 2014
 *      Author: richard
 */

#include "Fragmenter.h"
#include "MBEFrag.h"
#include "../libmints/molecule.h"
#include "Connections.h"
#include <boost/shared_ptr.hpp>
#include <iostream>
#include "exception.h"

namespace LibFrag {

void AddAtom(AtomSet& temp, SharedMol& Mol2Frag, const int i) {
   temp<<i;
   temp.AddMass(i, Mol2Frag->fmass(i));
   temp.AddCarts(i, Mol2Frag->fx(i), Mol2Frag->fy(i), Mol2Frag->fz(i));
}

std::vector<AtomSet> Fragmenter::MakeGroups(Connections& CTable,
      SharedMol& Mol2Frag) {
   std::vector<AtomSet> Groups;
   //Loop and find hydrogens, assign them and all atoms connected to them
   //to a group
   std::vector<bool> Assigned(Mol2Frag->natom(), false);
   for (int i=0; i<Assigned.size(); i++) {
      if (Mol2Frag->Z(i)==1&&!Assigned[i]) {
         if (CTable.GetNConnecs(i)!=1)
            throw psi::PSIEXCEPTION(
                  "This hydrogen is not connected to one and only one atom");
         AtomSet temp;
         AddAtom(temp, Mol2Frag, i);
         Assigned[i]=true;
         int center=CTable(i, 0)-1;
         if (Assigned[center])
            throw psi::PSIEXCEPTION("This H's center is already assigned");
         AddAtom(temp, Mol2Frag, center);
         Assigned[center]=true;
         for (int j=0; j<CTable.GetNConnecs(center); j++) {
            int atom=CTable(center, j)-1;
            //We already know it's attached to i
            if (atom!=i&&Mol2Frag->Z(atom)==1) {
               if (CTable.GetNConnecs(atom)!=1)
                  throw psi::PSIEXCEPTION(
                        "This hydrogen is not connected to one and only one atom");
               AddAtom(temp, Mol2Frag, atom);
               Assigned[atom]=true;
            }
         }
         temp.Sort();
         Groups.push_back(temp);
      }
   }
   for (int i=0; i<Assigned.size(); i++) {
      if (!Assigned[i]) {
         AtomSet temp;
         AddAtom(temp, Mol2Frag, i);
         Groups.push_back(temp);
      }
   }
   return Groups;
}

void Fragmenter::Groups2Frag(std::vector<AtomSet>& Groups,
      std::vector<int>& Groups2Add, int NGroups, NMerSet& Monomers,
      bool check) {
   SharedFrag temp(new MBEFrag(1));
   for (int i=0; i<NGroups; i++) {
      (*temp)*=Groups[Groups2Add[i]];
   }
   temp->Sort();
   bool Unique=true;
   for (int i=0; i<Monomers.size()&&check&&Unique; i++) {
      if ((*Monomers[i])>=(*temp)) Unique=false;
   }
   if (Unique) Monomers.push_back(temp);
}

bool Fragmenter::SortUnique(NMerSet& Monomers, NMerSet& temp) {
   std::vector<bool> Monomers2Erase(temp.size(), false);
   bool disjoint=true;
   for (int i=0; i<temp.size(); i++) {
      for (int j=i+1; j<temp.size()&&!Monomers2Erase[i]; j++) {
         if(!Monomers2Erase[j]){
            Set temp2=(*temp[i])/(*temp[j]);
            if (temp2.size()>1) disjoint=false;
            if (temp2.size()==temp[i]->size()) Monomers2Erase[i]=true;
            else if (temp2.size()==temp[j]->size()) Monomers2Erase[j]=true;
         }
      }
   }
   for (int i=0; i<Monomers2Erase.size(); i++) {
      if (!Monomers2Erase[i]) Monomers.push_back(temp[i]);
   }
   fprintf(psi::outfile,"The unique monomers are:\n");
   for (int i=0; i<Monomers.size(); i++)
      Monomers[i]->print_out();
   fprintf(psi::outfile,"************************\n");
   return disjoint;
}

void BondFragmenter::BondRecursion(std::vector<int>& BondedGroups,
      Connections* GConnec, std::vector<AtomSet>& Groups, NMerSet& Monomers,
      int depth, bool check, bool& severed) {
   if (depth<Bonds) {
      int CurrGroup=depth;
      bool goodfound=false;   //Have we found a path to follow
      while (CurrGroup>=0&&!goodfound) {

         //This is our current center
         int Center=BondedGroups[CurrGroup];

         //This is the minimal number of bonds a group must make
         int CurrBonds=GConnec->GetNConnecs(BondedGroups[0]);
         int MinBonds=(Bonds-depth>1&&CurrBonds>1 ? 2 : CurrBonds);

         //This is the number our current group makes,last recursion cycle
         //ensured NBonds>=MinBonds
         int NBonds=GConnec->GetNConnecs(Center);

         for (int bond=0; bond<NBonds; bond++) {

            //This is our new potential group
            int group=(*GConnec)(Center, bond)-1;

            //Make sure the number of bonds is ok
            int currbonds=GConnec->GetNConnecs(group);
            bool bondsgood=(currbonds>=MinBonds);

            //Do we need to check for uniqueness
            check=(currbonds==MinBonds);

            //Is this group already included
            bool backwards=false;
            for (int i=0; i<=depth&&!backwards; i++) {
               if (BondedGroups[i]==group) backwards=true;
            }

            if (!backwards&&bondsgood) {
               depth++;
               goodfound=true;
               BondedGroups[depth]=group;
               BondRecursion(BondedGroups, GConnec, Groups, Monomers, depth,
                     check,severed);
               //When we return our depth has decreased...
               depth--;
            }
         }
         if (!goodfound) CurrGroup--;
      }
      if (!goodfound) {
         //We can't manage to get NBonds so add these groups as a fragment
         //double check that this is actually unique and not some weird
         //fringe case
         Groups2Frag(Groups, BondedGroups, depth+1, Monomers, true);
      }

   }
   else {
      //Recursion has ended normally, add the group
      Groups2Frag(Groups, BondedGroups, Bonds+1, Monomers, check);
      //Once we broke a bond we need to make caps, so don't bother
      //continuing to check
      if(!severed){
         if(GConnec->GetNConnecs(BondedGroups[0])>1||
               GConnec->GetNConnecs(BondedGroups[Bonds])>1)severed=true;
         for(int j=1;j<Bonds&&!severed;j++){
            if(GConnec->GetNConnecs(BondedGroups[j])>2)severed=true;
         }
      }
   }
}

FragProps BondFragmenter::Fragment(SharedMol& AMol, NMerSet& Monomers) {
   Connections CTable(AMol);
   std::vector<AtomSet> Groups=MakeGroups(CTable,AMol);
   Connections GConnec(Groups,CTable);
   NMerSet temp;
   bool severed=false;
   for (int nbonds=0; nbonds<=GConnec.MaxBonds(); nbonds++) {
      for (int group=0; group<Groups.size(); group++) {
         int nconnec=GConnec.GetNConnecs(group);
         if (nconnec==nbonds&&nbonds!=0) {
            //You need Bonds+1 groups to have "Bonds" bonds
            std::vector<int> BondedGroups(Bonds+1, -1);
            ///The first group in this fragment is this group
            BondedGroups[0]=group;
            BondRecursion(BondedGroups, &GConnec, Groups, temp, 0, false,
                  severed);
         }
         else if (nconnec==nbonds&&nbonds==0) {
            std::vector<int> BondedGroups(1, group);
            Groups2Frag(Groups, BondedGroups, 1, temp, false);
         }
      }
   }
   bool disjoint=SortUnique(Monomers,temp);
   return FragProps(disjoint,severed);
}

FragProps UDFragmenter::Fragment(SharedMol& AMol, NMerSet& Monomers) {
   for (int frags=0,index=0; frags<AMol->nfragments(); frags++) {
      SharedMol Frag=AMol->py_extract_subsets_6(frags+1);
      int AtomsInFrag=Frag->natom();
      int zero=0;
      SharedFrag temp(new LibFrag::MBEFrag(1, &zero));
      Monomers.push_back(temp);
      for (int Atoms=0; Atoms<AtomsInFrag; Atoms++) {
         boost::shared_ptr<AtomSet> MonoI=Monomers[frags];
         (*MonoI)<<index;
         MonoI->AddMass(index, AMol->fmass(index));
         MonoI->AddCarts(index, AMol->fx(index), AMol->fy(index),
               AMol->fz(index));
         index++;
      }
   }
   return FragProps(true,false);
}

FragProps DistFragmenter::Fragment(SharedMol& AMol, NMerSet& Monomers) {
   Connections CTable(AMol);
   std::vector<AtomSet>Groups=MakeGroups(CTable,AMol);
   NMerSet tempset;
   for (int i=0; i<Groups.size(); i++) {
      SharedFrag temp=boost::shared_ptr<MBEFrag>(new MBEFrag(1));
      (*temp)*=Groups[i];
      for (int j=0; j<Groups.size(); j++) {
         if (i!=j && Groups[i].Distance(Groups[j])<=cutoff) {
            (*temp)*=Groups[j];
         }
      }
      tempset.push_back(temp);
   }
   bool disjoint=SortUnique(Monomers,tempset);
   bool severed=true;
   return FragProps(disjoint,severed);
}

}

