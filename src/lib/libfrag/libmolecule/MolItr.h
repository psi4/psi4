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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_MOLITR_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_MOLITR_H_

#include <vector>
#include <boost/shared_ptr.hpp>
#include "Atom.h"

namespace psi{
namespace LibMolecule{
class Molecule;
class IteratorGuts;

/** \brief An iterator that allows one to loop over the atoms in a molecule
 *
 *   If you are not familiar with iterators, you can think of them as fancy
 *   pointers.  Basically they allow you to access members of a class
 *   quickly, and in some predetermined order.  This makes it seem like
 *   the elements of the class are laid out contingously in memory, even
 *   if they are not.  They have similar syntax to pointers; the * operator
 *   "dereferences" it, and returns the current element.  The ++ operator
 *   advances to the next element; the -- operator rewinds an element.
 *
 *   For our new molecule class the reliance on different containers other
 *   than only std::vector means that looking an element up,
 *   by it's offset (i.e. MyMol[4]) is not efficient (although we do allow
 *   it).  Thus we provide iterators to efficiently loop over the class.
 *   Our molecular iterators are what as known as bidirectional
 *   iterators, you can increment or deincrement them at will.  They are
 *   not random access (you can't ask for the element that is 5 away
 *   from the current position).  Our iterators follow the same convention
 *   as the Standard Template Library, an iterator set to the beginning
 *   can be dereferenced to give the first element.  An iterator set to the
 *   end cannot be dereferenced as it is the first unacceptable state for
 *   the iterator, i.e. forward iteration termination conditions should
 *   be not less than or equal to the ending iterator, and reverse iteration
 *   conditions should be less than a beginning iterator.
 *
 *   A couple of general notes about iterators:
 *
 *   1) Comparisons of iterators should only be performed on iterators that
 *      originate from the same instance (not class).  That is to say, if
 *      I have an iterator to a water molecule, and one to a hydrogen
 *      molecule, call them Itr1 and Itr2 respectively, Itr1<Itr2 is a
 *      valid, but nonsensical comparison.  It makes no sense to see if
 *      we have iterated further through the water molecule than
 *      the hydrogen molecule (yes I realize that one could come up with
 *      scenarios in which this does make sense--we are 90% of the way
 *      through...---but that's not how the actual iterators work).
 *
 *   2) Iterators in general are no longer valid if the instance they belong
 *      to changes.  If I add/remove an atom to a molecule the count may
 *      get screwed up (whether it does or doesn't depends on where the
 *      atom was added/removed from).  If you modify a molecule, get a new
 *      set of iterators.
 *
 *   3) Never, under any circumstances, dereference an iterator that is at
 *      its end.  This is because you will likely cause a seg fault.
 *      By analogy to pointers, dereferencing the end iterator is
 *      like the following code (which will likely produce a seg fault):
 *      /code
 *      double *p=new double[4];
 *      for(int i=0;i<=4;i++){
 *         std::cout<<p[i]<<std::endl;
 *      /endcode
 *
 */
class MolItr{
   protected:
      ///This is the object that actually implements the iterator
      boost::shared_ptr<IteratorGuts> Impl_;
   public:
      MolItr(boost::shared_ptr<IteratorGuts> Impl);

      ///Base operators
      ///@{
      bool operator==(const MolItr& other)const;
      bool operator<(const MolItr& other)const;
      const MolItr& operator++();
      const MolItr& operator--();
      boost::shared_ptr<const Atom> operator*()const;
      ///@}

      ///Derived operators (only use base operators)
      ///@{
      bool operator<=(const MolItr& other)const{
         return !((*this)<other || (*this)==other);
      }
      bool operator>=(const MolItr& other)const{
         return !((*this)<other);
      }
      bool operator>(const MolItr& other)const{
         return !((*this)<=other);
      }
      /** \brief Allows access to the current atom's members
       *
       *  It's probably best to ignore the declaration of this operator, as
       *  the rules for overloading it are odd.  What is important is
       *  what this operator does, basically assume MyMol is a Molecule,
       *  then this operator is used like:
       *  /code
       *  MolItr Itr=MyMol.begin();
       *  std::cout<<"The first atom's atomic number is "
       *           <<Itr->Z();
       *  /endcode
       *  As you can see, this operator allows us to directly use the
       *  members of the atom the iterator is currently on.
       *
       */
      boost::shared_ptr<const Atom> operator->()const{return (*(*this));}
      ///@}
};

}}//End namespaces


#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_MOLITR_H_ */
