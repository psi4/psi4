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

#ifndef MBEFRAG_H_
#define MBEFRAG_H_
#include "AtomSet.h"
#include "psi4-dec.h"
#include <vector>
namespace psi {
namespace LibFrag {
/** \brief Fragments in the (G)MBE are comprised of atoms, caps, ghosts, and
 *   charges
 *
 *   Although, fragments are typically thought of as synonymous with the
 *   atoms that comprise them, they actually are comprised of several types
 *   of sets: atoms, caps, charges, and ghost atoms.  Each fragment needs
 *   to be able to specify a unique calculation on it's own.  Because the
 *   Caps, Charges, and Ghosts are intimately tied to the atoms in the
 *   fragment, we define superset, subset, and equality in terms of the
 *   equivalent relations on the atoms.
 *
 *   Fragments are also not easy to set up and this task is handled by
 *   four factories:
 *          Fragmenter: Makes fragments
 *          Capper:     Makes caps
 *          Embedder:   Makes charges
 *          BSSEer:     Makes ghosts
 *
 *    These four factories are friends of this class so that they can set
 *    it up.  Each factory is expected to leave the other three sets
 *    unchanged.  Furthermore they are expected to go through the factories
 *    in that order.  First, we need the bare bones of what comprises the
 *    fragment.  Second, we need to put on caps.  With that information
 *    we know which atoms should not be charges or ghost atoms.
 *
 *    The order of the other two factories is somewhat arbitrary, but I've
 *    defined it anyways to remove ambiguity.  Currently it is unclear
 *    whether one should be allowed to put basis functions on top of charges,
 *    if maybe we would only replace some atoms with charges/ghosts,
 *    or if BSSE/Embedding is an either-or thing.
 *
 */
class MBEFrag {
   private:
      ///Factory friends
      //@{
      friend class Fragmenter;
      friend class Capper;
      friend class Embedder;
      friend class BSSEer;
      //@}

      ///Performs a deep copy of the CartSets
      void Copy(const MBEFrag& other);

   protected:
      CartSet<SharedAtom> Atoms_;
      CartSet<SharedCap> Caps_;
      Set<SharedCharge> Charges_;
      Set<SharedGhost> Ghosts_;
      std::vector<int> Parents;
      int MBEOrder;

   public:
      /** \brief Accessors for the atoms,caps, ghosts, and
       *   charges.
       *
       *   If this MBEFrag doesn't have atoms, caps, charges, or
       *   ghosts the returned shared_ptr will evaluate to false.
       *   The top four take care of the boost symantics of upcasting
       *   for you.
       *
       */
      //@{
      const CartSet<SharedAtom>& Atoms() const {
         return Atoms_;
      }
      const CartSet<SharedCap>& Caps() const {
         return Caps_;
      }
      const Set<SharedGhost>& Ghosts() const {
         return Ghosts_;
      }
      const Set<SharedCharge>& Charges() const {
         return Charges_;
      }
      //@}

      /** \brief Set operations
       *
       * With the exception of the intersection and union operators,
       * these operations only worry about the atoms part of the
       * fragment.  For example: this->operator<(other) will return true
       * if the Atoms contained in this are a proper subset of those in
       * other.  The status of the other sets can be inferred from this
       * behavior and your logic should reflect this.
       *
       * For the intersection and union operations, the caps, ghosts, and
       * charges will be updated accordingly.  This means that under
       * union, the charges/ghosts will be the set difference
       * between the charges/ghosts in this and the atoms in other (taking
       * away caps/ghosts that are now atoms).  We have to check the
       * union of the caps manually
       *
       * For intersection, charges/ghosts are the union less
       * whatever atoms survive the intersection.  The search space for
       * caps is again the union of the caps, and needs checked manually.
       */
      //@{
      bool operator<(const MBEFrag& other) {
         return this->Atoms_<other.Atoms_;
      }
      bool operator>(const MBEFrag& other) {
         return this->Atoms_>other.Atoms_;
      }
      bool operator==(const MBEFrag& other) {
         return this->Atoms_==other.Atoms_;
      }
      bool operator<=(const MBEFrag& other) {
         return this->Atoms_<=other.Atoms_;
      }
      bool operator>=(const MBEFrag& other) {
         return this->Atoms_>=other.Atoms_;
      }
      bool operator!=(const MBEFrag& other) {
         return this->Atoms_!=other.Atoms_;
      }
      void operator*=(const MBEFrag& other);
      void operator/=(const MBEFrag& other);
      MBEFrag operator*(const MBEFrag& other) {
         MBEFrag temp(*this);
         temp*=other;
         return temp;
      }
      MBEFrag operator/(const MBEFrag& other) {
         MBEFrag temp(*this);
         temp/=other;
         return temp;
      }
      //@}

      int GetMBEOrder() {
         return MBEOrder;
      }
      int ParentI(const int I) {
         return Parents[I];
      }
      void SetMBEOrder(const int N) {
         MBEOrder=N;
      }

      void SetParents(const int* Ps) {
         Parents.clear();
         for (int i=0; i<MBEOrder; i++)
            Parents.push_back(Ps[i]);
      }

      void PrintParents() {
         for (int i=0; i<MBEOrder; i++)
            psi::outfile->Printf("%d ", Parents[i]);
      }

      MBEFrag(const int MBEOrder_=0, const int *Parents_=NULL) :
            MBEOrder(MBEOrder_) {
         if (Parents_!=NULL) SetParents(Parents_);
      }

      MBEFrag(const MBEFrag& other) {
         Copy(other);
      }

      MBEFrag& operator=(const MBEFrag& other) {
         if (this!=&other) Copy(other);
         return *this;
      }
};

}
}      //End namespaces

#endif /* MBEFRAG_H_ */
