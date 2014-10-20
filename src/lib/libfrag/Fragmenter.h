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

#ifndef FRAGMENTER_H_
#define FRAGMENTER_H_
#include "LibFragTypes.h"
#include "FragOptions.h"
#include "AtomSet.h"
#include "physconst.h"

namespace psi{
namespace LibFrag{
class Connections;

///Tiny class that tells us quickly the properties of our fragments
class FragProps{
   private:
      ///Copies our properties
      void Copy(const FragProps& other){
         this->Disjoint_=other.Disjoint_;
         this->Severed_=other.Severed_;
      }
   public:
      ///True if all fragments are disjoint
      bool Disjoint_;
      ///True if we broke any bonds
      bool Severed_;
      FragProps(const bool s=false,const bool d=true):
         Disjoint_(d),Severed_(s){}
      FragProps(const FragProps& other){this->Copy(other);}
      const FragProps& operator=(const FragProps& other){
         if(this!=&other)this->Copy(other);return *this;
      }
};
typedef boost::shared_ptr<CartSet<SharedAtom> > SharedAtomSet;

///The type of the groups
typedef std::vector<SharedAtomSet > GroupType;


/** \brief The class that is in charge of setting up the Atoms_ array of
 *  a MBEFrag that is a fragment or NMer.
 *
 *  The Fragmenter class is a friend of MBEFrag, which allows it to set
 *  the Atoms_ arrays.  The various algorithms that are implemented by
 *  deriving from this class are not.  What this means is your algorithm
 *  needs to implement the function FragmentImpl().  It should return
 *  a FragProps object, with Severed_ set correctly, and the GroupType
 *  object named Monomers should be set-up with the fragments you created.
 *  Disjoint_ of FragProps will be set when the uniqueness of your
 *  fragments is checked by Fragmenter.
 *
 *  It is the responsibility of the Fragmenter base class to set-up the
 *  MBEFrags, and in doing so will also check your Monomers for uniqueness
 *  so your algorithm need not do that, although the check will go faster
 *  if you do not generate redundant Monomers if at all possible.
 */
class Fragmenter{
   protected:
      /** \brief Function for converting copying groups to fragments
       *
       *  \param[in] Groups      The vector of groups from MakeGroups
       *  \param[in] Groups2Add  A vector of groups we are adding
       *  \param[in] NGroups     How many of the groups in Groups2Add are we
       *                         actually adding?
       *  \param[out] Monomers   The set of fragments plus our new frag
       *  \param[in]  Check      Should we check if the new frag is unique?
       */
      void Groups2Frag(GroupType& Groups,
            std::vector<int>& Groups2Add, int NGroups, GroupType& Monomers,
            bool check);

      ///Copies unique frags in temp to Monomers, returns true if disjoint
      bool SortUnique(NMerSet& Monomers,GroupType& temp);

      ///Establishes the universe for the Monomers
      void MakeUniv(const SharedMol& AMol, NMerSet& Monomers)const;

      ///Function that adds a fresh SharedAtomSet onto Frags
      inline void CommonInit(SharedMol& Mol2Frag,GroupType& Frags);

      virtual FragProps FragmentImpl(SharedMol& Mol2Frag,
            GroupType& Monomers)=0;

      ///If fragments are non-disjont this function appends the intersections
      void MakeIntersections(NMerSet& NMers)const;
   public:
        ///Returns properties of resulting fragments
		FragProps Fragment(SharedMol& Mol2Frag,NMerSet& Fragments);

		///Function that tells the Fragmenter to make NMers
		void MakeNMers(const NMerSet& Monomers,const int N, NMerSet& NMers,
		      bool Disjoint);

		virtual ~Fragmenter(){}
};

class UDFragmenter:public Fragmenter{
   protected:
      ///Returns properties of resulting fragments
      FragProps FragmentImpl(SharedMol& Mol2Frag,GroupType& Fragments);
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
            Connections* GConnec, GroupType& Groups,
            GroupType& Monomers, int depth, bool check,bool& severed,
            SharedMol& AMol);
   protected:
      ///Returns properties of resulting fragments
      FragProps FragmentImpl(SharedMol& Mol2Frag,GroupType& Fragments);
   public:
      BondFragmenter(int nbonds=2):Bonds(nbonds){}
};

class DistFragmenter:public Fragmenter{
   private:
      double cutoff;
   protected:
      ///Returns properties of resulting fragments
      FragProps FragmentImpl(SharedMol& Mol2Frag,GroupType& Fragments);
   public:
      DistFragmenter(double newcutoff=2.0):
         cutoff(newcutoff/pc_bohr2angstroms){}


};
}}//End namespaces
#endif /* FRAGMENTER_H_ */
