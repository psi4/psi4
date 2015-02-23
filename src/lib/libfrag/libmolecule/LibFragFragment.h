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
#include "LibFragMolecule.h"

namespace psi{
namespace LibMolecule{
class FragItrGuts;

/** \brief A unique identifier for each n-mer
 *
 *  As you manipulate your fragment set you'll quickly want to start playing
 *  tricks with them.  This will inevitably disturb the ordering.  To get
 *  around this we assign to each fragment a "serial number" that is
 *  comprised of an arbitrarily assigned number for a monomer, and the
 *  union of the serial numbers of the n monomers for an n-mer.  This
 *  quantity is actually an std::set so there is no ambiguity between
 *  trimer: "1 2 3" and dimer: "1 23" for example.  When looking things up
 *  it is always done by serial number.
 */
class SerialNumber:public std::set<unsigned int>{
   public:
      SerialNumber(){}
      ///Prints out the serial number
      std::string PrintOut()const;
};


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
 *     Mol_.  The primary purpose of this is for copying the object,
 *     so you can quite likely ignore this option
 *
 *  Note that all operators other than assignment,only apply to the
 *  Members_ array, e.g. the
 *  intersection of two fragments is equal to the intersection of the
 *  two Members_ arrays.  This is because such a relation between the other
 *  two arrays (OtherMembers_ and Atoms_) is ill defined.  A priori, one may
 *  speculate that the values of OtherMembers_ and Atoms_ should just be all
 *  the atoms in Mol_, that aren't in Members_ after the operation (or in
 *  set theory terms, the complement of Members_ in the universe defined by
 *  Mol_).  However, this is not true as some atoms may have shifted from
 *  point charges to caps, or to ghost atoms, etc.  Ultimately, the result
 *  is sensitive to how the user wants to apply BSSE corrections, embedding,
 *  and caps.  Hence you are responsible for recalculating the other
 *  arrays after an operation.
 *
 *
 */
class Fragment: public Molecule{
   private:
      ///Function that returns a particular atom (slow for this class)
      boost::shared_ptr<const Atom> LookUp(const int i)const;
      //Needed for the iterator to work
      friend class FragItrGuts;
      /** \brief Makes a shallow copy of fragment
       *
       *  When copying a fragment Mol_ is shallow copied, Members_,
       *  OtherMembers_, and Molecule::Atoms_ are all deep copied.
       */
      void Copy(const Fragment& other);
   protected:
      ///The molecule from which this fragment was derived
      boost::shared_ptr<const Molecule> Mol_;

      ///The atoms in Mol_ that belong to this fragment
      std::set<int> Members_;

      ///The SN of the fragment, a unique identification of this fragment
      SerialNumber SN_;

      /** For the caps, point charges, etc. it is often useful to know which
       *  atom they replaced, that info goes here
       */
      std::vector<int> OtherMembers_;

   public:
      ///Constructor for making a monomer that is number N
      Fragment(boost::shared_ptr<const Molecule> Mol,const long int N);
      ///Constructor for making an N-Mer (use union to fill it)
      Fragment(boost::shared_ptr<const Molecule> Mol);
      Fragment(const Fragment& other){
         this->Copy(other);
      }
      const Fragment& operator=(const Fragment& other){
         if(this!=&other)this->Copy(other);
         return *this;
      }

      ///An iterator that loops over all members (not just Members_)
      MolItr Begin()const;
      ///The iterator corresponding to the end of the begin call
      MolItr End()const;

      ///Allows you to add an atom, without copying it
      Fragment& operator<<(const int i){Members_.insert(i);return *this;}

      /** \brief Adds NewAtom, that replaces atom i in Mol
       *
       *  The atoms in Mol, are never actually touched (Mol is const after
       *  all).  Instead NewAtom is added to Atoms_, and i is appended to the
       *  back of OtherMembers_.
       */
      void AddRepAtom(const Atom& NewAtom,const int i);

      ///Add a new atom to fragment, that doesn't correspond to one in Mol
      Fragment& operator<<(const Atom& NewAtom);

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

      ///Returns true if the two fragments are equal
      bool operator==(const Fragment& other)const;

      ///Returns true if this fragment is a proper subset of other
      bool operator<(const Fragment& other)const;
      ///Returns true if this fragment is a subset of other
      bool operator<=(const Fragment& other)const{
         return ((*this)<other||(*this)==other);
      }

      ///Returns true if this fragment is a proper superset of other
      bool operator>(const Fragment& other)const{
         return (other<(*this));
      }
      ///Returns true if this fragment is a superset of other
      bool operator>=(const Fragment& other)const{
         return ((*this)>other||(*this)==other);
      }

      const SerialNumber& GetSN()const{return SN_;}

      ///No memory to free, so does nothing
      virtual ~Fragment(){}
};

}}



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_FRAGMENT_H_ */
