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
namespace LibFrag{
class Connections;

class Capper{
   protected:
      std::vector<double> AtomCarts;
      boost::shared_ptr<Connections> AtomConnect;
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
   public:
      ///Given a molecule sets up AtomConnect for you, and saves the
      ///molecule's carts to AtomCarts
      Capper(SharedMol& AMol);
      virtual ~Capper(){}
      virtual void MakeCaps(NMerSet& Set2Cap)=0;
};

class ReplaceAndCap:public Capper{
   public:
      ReplaceAndCap(SharedMol& AMol):Capper(AMol){}
      void MakeCaps(NMerSet& Set2Cap);
};

class ShiftAndCap:public Capper{
   private:
      ///The atomic numbers of the system
      std::vector<int> Zs;
   public:
      ShiftAndCap(SharedMol& AMol);
      ///Taken from Collins' fragmentation method
      void MakeCaps(NMerSet& Set2Cap);
};

}//End namespace


#endif /* CAPPER_H_ */
