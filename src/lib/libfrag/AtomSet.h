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
#ifndef ATOMSET_H_
#define ATOMSET_H_
#include "Set.h"
#include "LibFragTypes.h"

namespace LibFrag{
///A set of atoms, knows basic atom-y stuff like carts and masses
class AtomSet:public Set{
   private:

      ///Makes this a copy of other
      void Copy(const AtomSet& other);

      /** \brief This maps the atom number in the input file to data.
       *
       * Atoms in Elem2Atoms are specified by atom number from the input file
       * so as to prevent problems from arising from reording the elements
       * of the set, which is necessary because they must always be in
       * numeric order for the intersection and union operations to work.
       */
      std::map<int,Atom> Elem2Atoms;

      ///This set's center of mass (a.u.) with coords given to set
      std::vector<double> CoM;

      ///Calculates the CoM
      void CalcCoM();

   public:

      ///Returns the distance (a.u.) between the CoM of this set and "other"
      double Distance(AtomSet& other);

      ///Basic constructor
      AtomSet():Set(){}

      ///Default destructor
      ~AtomSet(){}

      ///Copy constructor
      AtomSet(const AtomSet& other):Set(other){Copy(other);}

      ///Assignment operator
      const AtomSet& operator=(const AtomSet&other);

      /**For BSSE purposes ghosts attached to this set, made it public for
      calls like Set.Ghosts[i] to distinguish from Set[i], which gives
      the real atom i*/
      std::vector<int> Ghosts;

      std::vector<Cap> Caps;
      ///Wrapper function to add ghost atoms
      void AddGhost(const int i){Ghosts.push_back(i);}

      ///Sets the carts of atom i, {x,y,z}
      void AddCarts(const int i, const double x,const double y,
            const double z);

      ///Sets the mass of atom i to m
      void AddMass(const int i,const double m);

      ///Returns the mass of atom i
      double Mass(const int i){return Elem2Atoms[Atoms[i]].mass;}

      ///Returns a vector of the carts
      std::vector<double> Carts(const int i){
         return Elem2Atoms[Atoms[i]].carts;
      }
};
}




#endif /* ATOMSET_H_ */
