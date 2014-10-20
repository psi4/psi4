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

#include "Fragmenter.h"
#include "MBEFrag.h"
#include "../libmints/molecule.h"
#include "Connections.h"
#include <boost/shared_ptr.hpp>
#include <iostream>
#include "exception.h"
#include "Grouper.h"

namespace psi {
namespace LibFrag {

inline void MakeSharedAtom(SharedMol& AMol, CartSet<SharedAtom>& Atoms) {
   for (int atom=0; atom<AMol->natom(); atom++) {
      int Z=AMol->Z(atom);
      double m=AMol->fmass(atom),x=AMol->fx(atom),y=AMol->fy(atom);
      double z=AMol->fz(atom);
      SharedAtom tempAtom(new Atom(Z, m, x, y, z));
      Atoms.AddSet(atom, tempAtom);
   }
}

inline void Fragmenter::CommonInit(SharedMol& AMol, GroupType& Frags) {
   Frags.push_back(SharedAtomSet(new CartSet<SharedAtom>()));
   if (Frags.size()==1) MakeSharedAtom(AMol, *Frags[0]);
   else { //Copy existing universe
      (*(Frags.back()))=(*(Frags[0]));
      //We also got Frags[0]'s atoms, get rid of them
      Frags.back()->Clear();
   }
}

FragProps Fragmenter::Fragment(SharedMol& AMol, NMerSet& Monomers) {
   GroupType AtomSets;
   FragProps Properties=FragmentImpl(AMol, AtomSets);
   Properties.Disjoint_=SortUnique(Monomers, AtomSets);
   if (!Properties.Disjoint_) MakeIntersections(Monomers);
   return Properties;
}

void Fragmenter::MakeNMers(const NMerSet& Monomers, const int N, NMerSet& Frags,
      bool Disjoint) {
   int NFrags=Monomers.size();
   NMerSet tempNMers;
   std::vector<int> indices;
   for (int i=0; i<N; i++)
      indices.push_back(i);
   int index=0;
   bool done=false;
   //Avoid a for loop because the ending condition increases factorially
   while (!done) {
      boost::shared_ptr<MBEFrag> temp(new MBEFrag(*Monomers[indices[0]]));
      temp->SetMBEOrder(N);
      temp->SetParents(&indices[0]);
      for (int i=1; i<N; i++)
         temp->Atoms_*=Monomers[indices[i]]->Atoms_;
      tempNMers.push_back(temp);
      bool igood=false;
      for (int i=N-1; i>=0&&!igood; i--) {
         int MaxMono=NFrags-(N-i);
         if (indices[i]<MaxMono) { //We can increase index i
            indices[i]++;
            igood=true;
            for (int j=i+1; j<N; j++)
               indices[j]=indices[i]+(j-i); //Reset indices after i
         }
         else if (i==0) {
            done=true;
            igood=true;
         }
      }
   }

   std::vector<bool> IsGood(tempNMers.size(), true);
   for (int i=0; i<tempNMers.size(); i++) {
      for (int j=i+1; j<tempNMers.size(); j++) {
         if (tempNMers[i]->Atoms()<=tempNMers[j]->Atoms()) IsGood[i]=false;
         else if (tempNMers[j]->Atoms()<=tempNMers[i]->Atoms()) IsGood[j]=false;
      }
   }
   for (int i=0,index=0; i<IsGood.size(); i++) {
      if (IsGood[i]) Frags.push_back(tempNMers[i]);
   }
   if (!Disjoint) MakeIntersections(Frags);
}

void Fragmenter::Groups2Frag(GroupType& Groups, std::vector<int>& Groups2Add,
      int NGroups, GroupType& Monomers, bool check) {
   if (NGroups>=1) {
      SharedAtomSet temp(new CartSet<SharedAtom>((*Groups[Groups2Add[0]])));
      for (int i=1; i<NGroups; i++)
         (*temp)*=(*Groups[Groups2Add[i]]);
      temp->Sort();
      bool Unique=true;
      for (int i=0; i<Monomers.size()&&check&&Unique; i++)
         if ((*Monomers[i])>=(*temp)) Unique=false;
      if (Unique) Monomers.back()=temp;
   }
}

bool Fragmenter::SortUnique(NMerSet& Monomers, GroupType& temp) {
   int NTemp=temp.size();
   std::vector<bool> Monomers2Erase(NTemp, false);
   bool disjoint=true;
   for (int i=0; i<NTemp; i++) {
      for (int j=i+1; j<NTemp&&!Monomers2Erase[i]; j++) {
         if (!Monomers2Erase[j]) {
            Set<SharedAtom> temp2=(*temp[i])/(*temp[j]);
            if (temp2.size()>1) disjoint=false;
            if (temp2.size()==temp[i]->size()) Monomers2Erase[i]=true;
            else if (temp2.size()==temp[j]->size()) Monomers2Erase[j]=true;
         }
      }
   }
   for (int i=0,index=0; i<NTemp; i++) {
      if (!Monomers2Erase[i]) {
         Monomers.push_back(SharedFrag(new MBEFrag(1, &index)));
         (Monomers.back())->Atoms_=(*temp[i]);
         index++;
      }
   }
   psi::outfile->Printf("The unique monomers are:\n");
   for (int i=0; i<Monomers.size(); i++)
      Monomers[i]->Atoms_.print_out();

   return disjoint;
}

void BondFragmenter::BondRecursion(std::vector<int>& BondedGroups,
      Connections* GConnec, GroupType& Groups, GroupType& Monomers, int depth,
      bool check, bool& severed, SharedMol& AMol) {
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
                     check, severed, AMol);
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
         CommonInit(AMol, Monomers);
         Groups2Frag(Groups, BondedGroups, depth+1, Monomers, true);
      }

   }
   else {
      //Recursion has ended normally, add the group
      CommonInit(AMol, Monomers);
      Groups2Frag(Groups, BondedGroups, Bonds+1, Monomers, check);
      //Once we broke a bond we need to make caps, so don't bother
      //continuing to check
      if (!severed) {
         if (GConnec->GetNConnecs(BondedGroups[0])>1
               ||GConnec->GetNConnecs(BondedGroups[Bonds])>1) severed=true;
         for (int j=1; j<Bonds&&!severed; j++) {
            if (GConnec->GetNConnecs(BondedGroups[j])>2) severed=true;
         }
      }
   }
}

FragProps BondFragmenter::FragmentImpl(SharedMol& AMol, GroupType& Monomers) {
   Connections CTable(AMol);   Grouper Fish;//Groupers are fish aren't they?
   GroupType Groups=Fish.MakeGroups(CTable,AMol);
   Connections GConnec(Groups, CTable);
   bool severed=false;
   for (int nbonds=0; nbonds<=GConnec.MaxBonds(); nbonds++) {
      for (int group=0; group<Groups.size(); group++) {
         int nconnec=GConnec.GetNConnecs(group);
         if (nconnec==nbonds&&nbonds!=0) {
            //You need Bonds+1 groups to have "Bonds" bonds
            std::vector<int> BondedGroups(Bonds+1, -1);
            ///The first group in this fragment is this group
            BondedGroups[0]=group;
            BondRecursion(BondedGroups, &GConnec, Groups, Monomers, 0, false,
                  severed, AMol);
         }
         else if (nconnec==nbonds&&nbonds==0) {
            std::vector<int> BondedGroups(1, group);
            CommonInit(AMol, Monomers);
            Groups2Frag(Groups, BondedGroups, 1, Monomers, false);
         }
      }
   }
   return FragProps(severed);
}

FragProps UDFragmenter::FragmentImpl(SharedMol& AMol, GroupType& Monomers) {
   int nfrags=AMol->nfragments();
   if(nfrags==1)throw PSIEXCEPTION("You didn't specify fragments in the input"
         " file.  Please do so or use a fragmentation algorithm.");

   for (int frags=0,index=0; frags<AMol->nfragments(); frags++) {
      SharedMol Frag=AMol->py_extract_subsets_6(frags+1);
      CommonInit(AMol, Monomers);
      int AtomsInFrag=Frag->natom();
      for (int Atoms=0; Atoms<AtomsInFrag; Atoms++) {
         (*Monomers.back())<<index++;
      }
   }
   return FragProps(false);
}

FragProps DistFragmenter::FragmentImpl(SharedMol& AMol, GroupType& Monomers) {
   Connections CTable(AMol);
   Grouper Fish;//Groupers are fish aren't they?
   GroupType Groups=Fish.MakeGroups(CTable, AMol);
   for (int i=0; i<Groups.size(); i++) {
      SharedAtomSet temp(new CartSet<SharedAtom>(*Groups[i]));
      for (int j=0; j<Groups.size(); j++) {
         if (i!=j&&Groups[i]->Distance((*Groups[j]))<=cutoff) {
            (*temp)*=(*Groups[j]);
         }
      }
      CommonInit(AMol, Monomers);
      Monomers.back()=temp;
   }
   bool severed=true;
   return FragProps(severed);
}

/** \brief A little "MACRO" to adjust the multiplicity correctly
 *
 *  \param[in,out] Mult The multiplicity we are adjusting
 *  \param[in] AreOpp True if the fragments we are comparing have the same
 *                    original sign
 *  \param[in] Sign   True if the fragment whose multiplicty we are adjusting
 *                     is positive
 */
inline void AdjustMult(int& Mult, const bool AreOpp, const bool Sign) {
   if (AreOpp&&Sign) Mult--;
   else if (AreOpp&&!Sign) Mult++;
   else if (!AreOpp&&Sign) Mult++;
   else if (!AreOpp&&!Sign) Mult--;
}

void Fragmenter::MakeIntersections(NMerSet& NMers) const {
   //This will be the final set of NMers, we do this by recursion so
   NMerSet TempSet;

   //In order to set the signs appropriately of each fragment
   //we need to keep track of it's original sign
   std::vector<bool> IsPositive;
   for (int i=0; i<NMers.size(); i++) {
      SharedFrag NMer=NMers[i];
      //Set for this iteration
      NMerSet ItrSet;
      for (int j=0; j<TempSet.size(); j++) {
         SharedFrag temp=boost::shared_ptr<MBEFrag>(new MBEFrag(*NMer));
         SharedFrag CompNMer=TempSet[j];
         temp->Atoms_/=CompNMer->Atoms_;
         int tempsize=temp->Atoms_.size();
         if (tempsize>0) {
            bool IsGood=true;
            temp->Mult_=(!IsPositive[j] ? 1 : -1);
            for (int k=0; k<TempSet.size()&&IsGood; k++) { //Check for uniqueness
               SharedFrag CompNMer2=TempSet[k];
               if (k==j) {   //Have a faster check
                  if (tempsize==CompNMer->Atoms().size()) {
                     IsGood=false;
                     //Always have the opposite sign
                     AdjustMult(CompNMer->Mult_, true, IsPositive[j]);
                  }
               }
               if (temp->Atoms_==CompNMer2->Atoms_) {
                  IsGood=false;
                  bool AreOpp=(IsPositive[j]&&IsPositive[k]);
                  AdjustMult(CompNMer2->Mult_, AreOpp, IsPositive[k]);
               }
            }
            if (IsGood) {
               ItrSet.push_back(temp);
               IsPositive.push_back(!IsPositive[j]);
            }

         }
      }
      if (ItrSet.size()>0)
         TempSet.insert(TempSet.end(), ItrSet.begin(), ItrSet.end());
      TempSet.push_back(NMer);
      IsPositive.push_back(true);
   }
   NMers=TempSet;
}

}
}      //End namespaces

