/*
 * Fragmenter.h
 *
 *  Created on: May 15, 2014
 *      Author: richard
 */

#ifndef FRAGMENTER_H_
#define FRAGMENTER_H_
#include "LibFragTypes.h"
#include "FragOptions.h"
#include "AtomSet.h"
#include "physconst.h"

namespace LibFrag{
class Connections;

class FragProps{
   private:
      void Copy(const FragProps& other){
         this->disjoint=other.disjoint;
         this->severed=other.severed;
      }
   public:
      bool disjoint;
      bool severed;
      FragProps(const bool d=true,const bool s=false):
         disjoint(d),severed(s){}
      FragProps(const FragProps& other){this->Copy(other);}
      const FragProps& operator=(const FragProps& other){
         if(this!=&other)this->Copy(other);return *this;
      }
};

class Fragmenter{
   protected:
      ///Function takes a SharedMol and returns a vector of the groups in it
      std::vector<AtomSet> MakeGroups(Connections& CTable,
            SharedMol& Mol2Frag);

      /** \brief Function for converting groups to fragments
       *
       *  \param[in] Groups      The vector of groups from MakeGroups
       *  \param[in] Groups2Add  A vector of groups we are adding
       *  \param[in] NGroups     How many of the groups in Groups2Add are we
       *                         actually adding?
       *  \param[out] Monomers   The set of fragments plus our new frag
       *  \param[in]  Check      Should we check if the new frag is unique?
       */
      void Groups2Frag(std::vector<AtomSet>& Groups,
            std::vector<int>& Groups2Add, int NGroups, NMerSet& Monomers,
            bool check);

      ///Copies unique frags in temp to Monomers, returns true if disjoint
      bool SortUnique(NMerSet& Monomers,NMerSet& temp);
   public:
        ///Returns properties of resulting fragments
		virtual FragProps Fragment(SharedMol& Mol2Frag,NMerSet& Fragments)=0;
		virtual ~Fragmenter(){}
};

class UDFragmenter:public Fragmenter{
   public:
      ///Returns properties of resulting fragments
      FragProps Fragment(SharedMol& Mol2Frag,NMerSet& Fragments);
};

class BondFragmenter:public Fragmenter{
   private:
      int Bonds;
      /** \brief Function for iterating over all bonded groups that are
       *       at most BondFragmenter::Bonds apart.
       *
       *
       *  We do this by recursion.  We assume that our top-level loop is
       *  giving us groups to loop over in terms of the number of connections
       *  i.e., we first get all groups with 0 connections, then 1, etc.
       *
       *  Consider some branched molecule, with groups arranged as:
       *
       *         1-2
       *           |
       *         5-3-4
       *           |
       *           6-7
       *
       *  and the scenario where we want fragments from groups that are two
       *  bonds apart.  Clearly we need fragments 123, 234, 235, 236, 345,
       *  346, 356, and 367.  Because we go in order of the number of
       *  connections, we first obtain: 123, 234, 345, 346, 235, 356, and 367.
       *  Now when we are at two connections our first group is 2, we only
       *  need to consider groups bonded to 2 with at least two bonds
       *  (i.e. group 3) because group 1 makes one bond all fragments it is
       *  part of are already taken care of.  We thus obtain 236, our
       *  missing fragment.  Other than 6, which can go 6-3-2 no other
       *  group has a path of two groups that form at least two bonds.
       *
       *  We derive from this experience that our path is terminated when
       *  we have no avenues of groups that make the same or higher numbers
       *  of bonds.  If in the path we use another group with the same
       *  number of groups we need to take care to ensure that the resulting
       *  fragment is unique.
       *
       *  Unfortunately that scenario was too simple, consider:
       *
       *     1   8   9
       *      \ / \ /
       *       2   7
       *      /     \
       *     3       6
       *      \     /
       *       4---5
       *
       *  Our above algorithm gets everything except for 287, which is
       *  missed because the 2 and 7 never loop to the 8.  We must
       *  therefore relax our restriction to never considering lower
       *  connectivities unless there is a path through them to a higher
       *  connectivity, i.e. if we have more than one bond left in our
       *  recursion we must consider groups that make two bonds, regardless
       *  of how many bonds our seed makes.
       *
       *  \params[in] BondedGroups A preallocated, "Bonds" long vector of
       *                           the groups we have found for the current
       *                           fragment
       *  \params[in] GConnec      The Group connectivity table
       *  \params[in] Groups       The complete set of groups
       *  \params[out] Monomers    The array of monomers that we are appending
       *                           fragments to as we find them
       *  \params[in] depth        The level of recursion we are currently
       *                           at
       *  \params[in] check        Whether we have to check for uniqness
       *
       */
      void BondRecursion(std::vector<int>& BondedGroups,
            Connections* GConnec, std::vector<AtomSet>& Groups,
            NMerSet& Monomers, int depth, bool check,bool& severed);
   public:
      BondFragmenter(int nbonds=2):Bonds(nbonds){}

      ///Returns properties of resulting fragments
      FragProps Fragment(SharedMol& Mol2Frag,NMerSet& Fragments);
};

class DistFragmenter:public Fragmenter{
   private:
      double cutoff;
   public:
      DistFragmenter(double newcutoff=2.0):
         cutoff(newcutoff/pc_bohr2angstroms){}

      ///Returns properties of resulting fragments
      FragProps Fragment(SharedMol& Mol2Frag,NMerSet& Fragments);
};
}
#endif /* FRAGMENTER_H_ */
