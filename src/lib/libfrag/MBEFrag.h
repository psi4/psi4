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

enum MBEProps{MBE_ENERGY,MBE_DENSITY};

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
      std::vector<int> Parents_;
      int MBEOrder_;
      /** \brief The number of times this fragment appears
       *
       *  For fragments that use the MBE this will always be 1, unless
       *  the system is symmetric.  If the system is symmetric this will
       *  be the number of identical occurrences of this fragment.  In
       *  particular note that this means that it is not equal to the
       *  binomial coefficient and sign of each term in the MBE, the thought
       *  being we plan on returning the value of the property for all
       *  orders lower than or equal to the current order so the coefficient
       *  and sign changes.
       *
       *  TODO: Account for symmetry
       *  (Note: Symmetry doesn't actually work right now).
       *
       *  For the GMBE this will be the fragment/intersection's signed
       *  coefficent.  This is because I do not know a formula that allows
       *  me to recursively find the value of the property for all lower
       *  orders at the time of writing.
       */
      int Mult_;
      //std::map<MBEProps,Property> Properties_;
   public:
      /** \brief Accessors for the atoms,caps, ghosts, and
       *   charges.
       *
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

      int Mult()const{return Mult_;}
      int GetMBEOrder()const {
         return MBEOrder_;
      }
      int ParentI(const int I) const{
         return Parents_[I];
      }
      void SetMBEOrder(const int N) {
         MBEOrder_=N;
      }

      void SetParents(const int* Ps) {
         Parents_.clear();
         for (int i=0; i<MBEOrder_; i++)
            Parents_.push_back(Ps[i]);
      }

      std::string PrintParents()const {
         std::stringstream parents;
         for (int i=0; i<MBEOrder_; i++)
            parents<<Parents_[i]<<" ";
         return parents.str();
      }

      MBEFrag(const int MBEOrder=0, const int *Parents=NULL,const int Mult=1) :
            MBEOrder_(MBEOrder), Mult_(Mult) {
         if (Parents!=NULL) SetParents(Parents);
      }

      MBEFrag(const MBEFrag& other) {
         Copy(other);
      }

      const MBEFrag& operator=(const MBEFrag& other) {
         if (this!=&other) Copy(other);
         return *this;
      }
};

}
}      //End namespaces

#endif /* MBEFRAG_H_ */
