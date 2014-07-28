/*
 * GMBE.cc
 *
 *  Created on: May 15, 2014
 *      Author: richard
 */

#include "GMBE.h"
#include "MBEFrag.h"
#include "psi4-dec.h"
namespace LibFrag {

void GMBE::MakeNmers(const NMerSet& Monomers, NMerSet& NMers) {
   int NFrags=Monomers.size();
   std::vector<int> indices;
   for (int i=0; i<N; i++)
      indices.push_back(i);
   int index=0;
   bool done=false;
   while (!done) {
      MBEFrag temp=(*Monomers[indices[0]]);
      temp.SetMBEOrder(N);
      temp.SetParents(&indices[0]);
      for (int i=1; i<N; i++)
         temp*=(*Monomers[indices[i]]);
      NMers.push_back(boost::shared_ptr<MBEFrag>(new MBEFrag(temp)));
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

   std::vector<bool> IsGood(NMers.size(), true);
   for (int i=0; i<NMers.size(); i++) {
      for (int j=i+1; j<NMers.size(); j++) {
         if ((*NMers[i])<=(*NMers[j])) IsGood[i]=false;
         else if ((*NMers[j])<=(*NMers[i])) IsGood[j]=false;
      }
   }
   for (int i=0,index=0; i<IsGood.size(); i++) {
      if (!IsGood[i]) NMers.erase(NMers.begin()+i-(index++));
   }
   /*fprintf(psi::outfile, "The Unique N-Mers:\n");
   for (int i=0; i<NMers.size(); i++) {
      NMers[i]->print_out();
   }
   fprintf(psi::outfile, "******************\n");*/
}

inline void MakeInts(NMerSet& ExistingNMers, SharedFrag& NMer,
      NMerSet& FoundNMers, std::vector<int>& Mults) {
   for (int j=0; j<ExistingNMers.size(); j++) {
      SharedFrag temp=boost::shared_ptr<MBEFrag>(new MBEFrag((*NMer)));
      (*temp)/=(*ExistingNMers[j]);
      if (temp->size()>0) {
         if (temp->size()!=ExistingNMers[j]->size()) FoundNMers.push_back(temp);
         else //This means it is the same as the n-mer
         Mults[j]--;
      }
   }
}

inline void CheckUnique(NMerSet& SamePhase, NMerSet& OpPhase, NMerSet& Array,
      std::vector<bool>& IsGood, std::vector<int>& SameMults,
      std::vector<int>& OpMults) {
   for (int j=0; j<Array.size(); j++) {
      for (int k=0; k<SamePhase.size()&&IsGood[j]; k++) {
         if ((*SamePhase[k])==(*Array[j])) {
            IsGood[j]=false;
            SameMults[k]++;
         }
      }
      for (int k=0; k<OpPhase.size()&&IsGood[j]; k++) {
         if ((*OpPhase[k])==(*Array[j])) {
            IsGood[j]=false;
            OpMults[k]--;
         }
      }
   }
}

inline void AddInt(NMerSet& Exist, NMerSet& Curr, std::vector<int>& Mult,
      std::vector<bool>& IsGood) {
   for (int j=0; j<Curr.size(); j++) {
      if (IsGood[j]) {
         Exist.push_back(Curr[j]);
         Mult.push_back(1);
      }
   }
}

void GMBE::MakeIntersections(std::vector<NMerSet>& Systems) {
   //These will hold the intersections between odd numbers of NMers
   //and even numbers respective (that refers to the original appearance) of
   //an intersection, in the end the only thing that matters is the signs
   NMerSet PInts,NInts;
   for (int i=0; i<Systems[1].size(); i++) {
      SharedFrag NMer=Systems[1][i];
      //The equivalent P and N Ints we generate this cycle, need checked for
      //uniqueness
      NMerSet PMadeNMers,NMadeNMers;

      //Calls that make the ints
      MakeInts(PInts, NMer, NMadeNMers, PMults);
      MakeInts(NInts, NMer, PMadeNMers, NMults);

      //Check for uniqueness and then add to arrays
      std::vector<bool> PIsGood(PMadeNMers.size(), true);
      std::vector<bool> NIsGood(NMadeNMers.size(), true);
      CheckUnique(PInts, NInts, PMadeNMers, PIsGood, PMults, NMults);

      //Add the PMadeNMers in, that way the NMadeNMers will now be checked
      //against them
      AddInt(PInts, PMadeNMers, PMults, PIsGood);
      CheckUnique(NInts, PInts, NMadeNMers, NIsGood, NMults, PMults);
      AddInt(NInts, NMadeNMers, NMults, NIsGood);

      //Add our n-mer
      PInts.push_back(NMer);
      PMults.push_back(1);
   }
   Systems[1]=PInts;
   Systems.push_back(NInts);
}

double GMBE::Energy(const std::vector<NMerSet>& Systems,
      const std::vector<double*>& Energies) {
   double TEnergy=0.0;
   //If N==1 special case and energies are only in Energies[0],but
   //we get sizes from Systems[1] and Systems[2]
   int index1=(N==1?0:1);
   for (int i=0; i<Systems[1].size(); i++)
      TEnergy+=PMults[i]*Energies[index1][i];

   for (int j=0; j<Systems[2].size(); j++) {
      int index2=(N==1?Systems[1].size()+j:j);
      TEnergy-=NMults[j]*Energies[index1][index2];
   }
   fprintf(psi::outfile, "The total %d-body GMBE energy is: %16.12f (a.u.)", N,
         TEnergy);
   return TEnergy;
}
}

