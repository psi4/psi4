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
#ifndef CAPPER_H_
#define CAPPER_H_
#include<sstream>
#include<boost/python.hpp>
#include "AtomSet.h"
#include "LibFragTypes.h"

namespace psi{
namespace LibFrag{
class Connections;

typedef std::vector<CartSet<SharedCap> > CapType;


/** \brief The factory that makes the caps for the fragments.
 *
 *  Unlike the other three factories, this one has no common initialization
 *  function, because the carts of the caps depends on the algorithm. Each
 *  algorithm is therefore responsible for making it's own CapType array.
 *  There is an unofficial initializer that is designed to set-up
 *  "number of fragments" empty cap sets.
 */
class Capper{
   protected:
      std::vector<double> AtomCarts_;
      boost::shared_ptr<Connections> AtomConnect_;
      /** \brief Returns array of bonds to cap
       *
       *  If the returned array is called bonds, nbonds=bonds.size()/3, and
       *  elements bonds[i*2+0] through bonds[i*2+2] are associated with the
       *  i-th severed bond. bonds[i*2+0] is the fragment that is getting
       *  the cap, bonds[i*2+1] is the atom in the fragment that the cap
       *  attaches to and bonds[i*2+2] is the atom *not* in the fragment that
       *  the cap is replacing.
       */
      std::vector<int> Atoms2Replace(NMerSet& Set2Cap)const;
      virtual void CapImpl(CapType& Caps,NMerSet& Set2Cap)=0;

   public:
      ///Given a molecule sets up AtomConnect for you, and saves the
      ///molecule's carts to AtomCarts
      Capper(SharedMol& AMol);
      virtual ~Capper(){}
      void MakeCaps(NMerSet& Set2Cap);
};

class ReplaceAndCap:public Capper{
   protected:
      void CapImpl(CapType& Caps, NMerSet& Set2Cap);
   public:
      ReplaceAndCap(SharedMol& AMol):Capper(AMol){}

};

class ShiftAndCap:public Capper{
   private:
      ///The atomic numbers of the system
      std::vector<int> Zs_;
      ///Taken from Collins' fragmentation method
      void CapImpl(CapType& Caps,NMerSet& Set2Cap);
   public:
      ShiftAndCap(SharedMol& AMol);

};

}}//End namespaces


#endif /* CAPPER_H_ */
