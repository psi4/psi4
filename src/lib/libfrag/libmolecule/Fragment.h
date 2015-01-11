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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_FRAGMENT_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_FRAGMENT_H_
#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>
#include "Molecule.h"

namespace psi{
namespace LibMolecule{
class FragItrGuts;

/** \brief Part of a molecule, usable as if it was it's own molecule
 *
 *  A fragment is associated with a molecule, and the data of the molecule
 *  is only carried around via a pointer.  Operations like unions and
 *  intersections are performed on the Members_ array, with the assumption
 *  being that two different indices correspond to two different atoms.  Being
 *  a derived class of molecule it also has it's own array of atoms, this
 *  is where ghost atoms can be added (or caps, or point charges, i.e. things
 *  that are present for this molecule, but not the total molecule).
 *
 *  There are several ways to add things to this fragment, and it's important
 *  to understand the differences.  For this discussion we assume this
 *  fragment is a part of a molecule called Mol_.
 *  1) operator<<(const int i): this is how you add atoms from Mol_ to the
 *     fragment, i.e. to add atom 0 of Mol_ to the fragment the syntax is:
 *     \code
 *     (*this)<<0;
 *     \endcode
 *  2) AddRepAtom(boost::shared_ptr<const Atom> NewAtom, const int i): this
 *     is how you add an atom that replaces an atom in Mol_, e.g. this is how
 *     you replace atom i in Mol_ with a ghost atom
 *  3) operator<<(boost::shared_ptr<const Atom> NewAtom): this is how you
 *     add an atom to the fragment that doesn't have a corresponding atom in
 *     Mol_.  I'm not sure why you'd really use this...
 */
class Fragment: public Molecule{
   private:
      boost::shared_ptr<const Atom> LookUp(const int i)const;
      friend class FragItrGuts;
   protected:

      ///The molecule from which this fragment was derived
      boost::shared_ptr<const Molecule> Mol_;

      ///The atoms in Mol_ that belong to this fragment
      std::set<int> Members_;

      /** For the caps, point charges, etc. it is often useful to know which
       *  atom they replaced, that info goes here
       */
      std::vector<int> OtherMembers_;

   public:
      Fragment(boost::shared_ptr<const Molecule> Mol);

      MolItr Begin()const;
      MolItr End()const;

      Fragment& operator<<(const int i){Members_.insert(i);return *this;}

      /** \brief Adds NewAtom, that replaces atom i in Mol
       *
       *  The atoms in Mol, are never actually touched (Mol is const after
       *  all).  Instead NewAtom is added to Atoms_, and i is appended to the
       *  back of OtherMembers_.
       */
      void AddRepAtom(boost::shared_ptr<const Atom> NewAtom,const int i);

      ///Adds a new atom to this fragment
      Fragment& operator<<(boost::shared_ptr<const Atom> NewAtom);

      ///Returns the number of atoms in this fragment (members+caps etc.)
      virtual int NAtoms()const;

      ///Makes this fragment the union of itself with other
      const Fragment& operator+=(const Fragment& other);
      ///Returns the union of this fragment with other
      Fragment operator+(const Fragment& other)const;


      ///Makes this fragment the intersection of itself with other
      const Fragment& operator-=(const Fragment& other);
      ///Returns the intersection of this fragment with other
      Fragment operator-(const Fragment& other)const;

      virtual ~Fragment(){}
};

}}



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_FRAGMENT_H_ */
