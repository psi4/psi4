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

#ifndef CONNECTIONS_H_
#define CONNECTIONS_H_

#include<boost/shared_ptr.hpp>
#include "libmints/matrix.h"
#include "AtomSet.h"
#include "LibFragTypes.h"


namespace psi{
typedef boost::shared_ptr<Matrix> SharedMatrix;
namespace LibFrag{
typedef boost::shared_ptr<CartSet<SharedAtom> > SharedAtomSet;

///The type of the groups
typedef std::vector<SharedAtomSet> GroupType;
class Connections{
   private:
      ///The connectivity table
      std::vector<int> CTable;
      ///The number of connections each atom has
      std::vector<int> NConnecs;
      ///The distance matrix;
      Matrix Distance;
      ///Things below bondthresh are likely bonds (in a.u.)
      double bondthresh;
      ///Things can only have this many bonds
      const static int nbonds=4;
      /** Makes the connectivity table based on covalent radii*/
      void RoughTable(SharedMol Mol=SharedMol());
   public:
      ///Really should actually check, but I'm lazy...
      int MaxBonds(){return nbonds;}
      Connections(SharedMol& Mol);
      ///Given a set of groups and the atom connectivity table makes a group
      ///connectivity table
      Connections(GroupType& Groups,const Connections& AtomConnec);
      int GetNConnecs(const int Atom)const{return NConnecs[Atom];}
      /**Returns the "Conn"-th atom bonded to "Atom", a return value of 0
         means "Atom" doesn't make that many bonds or you requested a value
         greater than equal to nbonds.  Consequentially, a
         return value of 1 means the 0-th atom of the system, etc.*/
      int operator()(const int Atom,const int Conn)const{
         return (Conn<nbonds?CTable[Atom*nbonds+Conn]:0);
      }

      ///Debugging function
      void print_out();
};
}}//End namespaces


#endif /* CONNECTIONS_H_ */
