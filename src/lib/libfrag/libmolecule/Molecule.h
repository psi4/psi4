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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_MOLECULE_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_MOLECULE_H_
#include<vector>
#include<boost/shared_ptr.hpp>

#include "MolItr.h"
namespace psi{
namespace LibMolecule{
class Atom;
class Molecule:protected LibMoleculeBase{
   private:
      int Spin_;
      int Charge_;
   protected:
      std::vector<boost::shared_ptr<const Atom> > Atoms_;

      /** \brief Derived classes will use this function to
       *         enforce a particular looping order
       *
       *   What the brief means is best explained by example.
       *   Consider the Fragment derived class, it needs
       *   to loop over its own Members array before it loops
       *   over the Atoms array.  This is the function it
       *   redefines to accomplish that.
       */
      virtual boost::shared_ptr<const Atom> LookUp(const int i)const{
         return Atoms_[i];
      }
   public:
      int Spin()const{return Spin_;}
      int Charge()const{return Charge_;}

      /** \name Iterator accessors
       *
       *  The following two functions are all that is needed to loop over
       *  the atoms in this molecule using a syntax akin to the standard
       *  template library
       */
      ///@{
      virtual MolItr Begin()const;
      virtual MolItr End()const;
      ///@}

      /** \name Standard accessors
       *
       *  The three fxns in this block are used to
       *  access atoms based on their index, i.e. give
       *  me atom 1.  When random access is not needed, i.e.
       *  you are going to loop over some contigious subset of atoms,
       *  the use of iterators are preferred       *
       */
      ///@{
      ///Returns the number of atoms
      virtual int NAtoms()const{return Atoms_.size();}
      ///Returns atom i
      boost::shared_ptr<const Atom> operator[](const int i)const{
         return LookUp(i);
      }
      ///Returns Cartesian component i of the I-th atom
      double operator()(const int AtomI, const int i)const;
      ///@}

      ///Adds an atom to this molecule
      virtual Molecule& operator<<(boost::shared_ptr<const Atom> NewAtom){
         Atoms_.push_back(NewAtom);
         return *this;
      }

      ///Prints out the molecule
      std::string PrintOut(const int DebugLevel=1)const;
      virtual ~Molecule(){}
};


}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_MOLECULE_H_ */
