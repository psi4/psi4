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
#include "Implementations/MoleculeGuts.h"

namespace psi{
namespace LibMolecule{
/** \brief The public interface to the Molecule class
 *
 *   A typical user should only need things from this class
 */
class Molecule:protected MoleculeGuts{
   public:
      ///Returns the Multiplicity of the Molecule
      int Mult()const{return Mult_;}
      ///Returns the charge of the Molecule
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
      virtual Molecule& operator<<(const Atom& NewAtom){
         Atoms_.push_back(
               boost::shared_ptr<const Atom>(new Atom(NewAtom)));
         return *this;
      }

      ///Prints out the molecule
      std::string PrintOut(const int DebugLevel=1)const;

      ///Makes this molecule the union of this and other
      virtual const Molecule& operator+=(const Molecule& other);
      ///Returns the union of this molecule and other
      Molecule operator+(const Molecule& other);

      Molecule():MoleculeGuts(){}
      Molecule(const Molecule& other):MoleculeGuts(other){}
      virtual const Molecule& operator=(const Molecule& other){
         MoleculeGuts::operator=(other);
         return *this;
      }
      virtual ~Molecule(){}
      ///True if Atom::operator==() is true for all atoms in the molecule
      virtual bool operator==(const Molecule& other)const;
      ///Returns the opposite of operator==
      bool operator!=(const Molecule& other)const{
         return !((*this)==other);
      }
};


}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_MOLECULE_H_ */
