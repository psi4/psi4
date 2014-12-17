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

/** \brief A class for holding the connections of either atoms or things
 *   that are of type GroupType
 *
 *   There is a need to have place-holder data in the table (i.e. not all
 *   atoms make four bonds, but memory for four bonds was allocated). In
 *   such cases, for better or for worse, I have chosen to use 0 as my
 *   placeholder.  Thus a return of 0 from (*this)(i,j) means that atom
 *   i (counting from 0) does not form j bonds.  If the return is non zero
 *   then it is the number of the atom COUNTING FROM 1!!!!!!
 */
class Connections{
   private:
      ///The connectivity table
      std::vector<int> CTable_;
      ///The number of connections each atom has
      std::vector<int> NConnecs_;
      ///The distance matrix;
      Matrix Distance_;
      ///Things below bondthresh are likely bonds (in a.u.)
      double bondthresh_;
      ///Things can only have this many bonds
      const static int nbonds_=4;
      /** Makes the connectivity table based on covalent radii*/
      void RoughTable(SharedMol Mol=SharedMol());
   public:
      ///Right now returns the maximum number of bonds an atom may have
      int MaxBonds()const{return nbonds_;}
      Connections(SharedMol& Mol);
      ///Given a set of groups and the atom connectivity table makes a group
      ///connectivity table
      Connections(GroupType& Groups,const Connections& AtomConnec);
      int GetNConnecs(const int Atom)const{return NConnecs_[Atom];}
      /**Returns the "Conn"-th atom bonded to "Atom", a return value of 0
         means "Atom" doesn't make that many bonds or you requested a value
         greater than equal to nbonds.  Consequentially, a
         return value of 1 means the 0-th atom of the system, etc.*/
      int operator()(const int Atom,const int Conn)const{
         return (Conn<MaxBonds()?CTable_[Atom*MaxBonds()+Conn]:0);
      }

      ///Debugging function
      void print_out()const;
};
}}//End namespaces


#endif /* CONNECTIONS_H_ */
