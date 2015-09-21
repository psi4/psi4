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
#ifndef PSIAPI_GAUSSIANBASIS_H_
#define PSIAPI_GAUSSIANBASIS_H_
#include "BasisSet.h"
#include "Utils/NestedIterator.h"

///Macro for defining a clone and copy function with minimal retyping
#define CLONE(ClassName)\
   virtual ClassName* Clone()const{return new ClassName(*this);}\
   ClassName(const ClassName& Other):BasisSet(Other){}


namespace PsiAPI{
///A Gaussian Primitive
class GaussianPrim: public BasisSet{
   public:
      ///\name Constructors
      ///@{
      ///Members for constructing a Gaussian primitive

      CLONE(GaussianPrim)

      ///Makes a primitive with exponent "Alpha" and coefficient "Coef"
      GaussianPrim(double Alpha,double Coef){
        AddParam(EXPONENT,Alpha);
        AddParam(COEF,Coef);
      }
      ///@}

      ///@{ Convenience functions for accessing parameters
      ///Convenience function to return the exponent
      double Exponent()const{return GetParam(EXPONENT);}
      ///Convenience function to return the coefficient
      double Coef()const{return GetParam(COEF);}
      ///@}

      ///Prettier printing of the parameters
      operator std::string()const;
};

///A Contracted Gaussian Atomic Orbital
class GaussianAO: public BasisSet{
   public:

      CLONE(GaussianAO)

      ///Default constructor, does nothing
      GaussianAO(){}
      ///Loops over all primitive printing w/o printing params
      operator std::string()const;
};

///A set of Gaussian Atomic Orbitals all with the same angular momentum
class GaussianShell: public GaussianAO{
   public:
      ///\name Constructors
      ///@{
      ///Functions for making a GaussianShell

      //Can't use macro because BasisSet isn't a direct base class
      virtual GaussianShell* Clone()const{return new GaussianShell(*this);}
      GaussianShell(const GaussianShell& Other):GaussianAO(Other){}


      ///Makes a shell of angular momentum "L"
      GaussianShell(size_t L,bool IsPure){
         AddParam(ANGMOM,(double)L);
         AddParam(PURE,(double)IsPure);
      }
      ///@}

      ///@{Convenience access functions
      ///Convenience function for telling if the shell is spherical or not
      bool IsPure()const{return (bool)GetParam(PURE);}
      ///Convenience function returning the total angular momentum of the shell
      size_t L()const{return (size_t)GetParam(ANGMOM);}
      ///@}

      ///Overrides NBasis to be correct
      size_t NBasis()const{return IsPure()?2*L()+1:(L()+1)*(L()+2)/2;}

      ///Adds angular momentum and a flag for Cartesian/Pure
      operator std::string()const;
};

///Class to hold on shells centered on the same center
class GaussianAtomicBasis: public BasisSet{
   public:
      CLONE(GaussianAtomicBasis)

      ///Makes a GaussianAtomicBasis from contigious memory (a.u.)
      GaussianAtomicBasis(const double* Carts){
         AddParam(CENTERX,Carts[0]);
         AddParam(CENTERY,Carts[1]);
         AddParam(CENTERZ,Carts[2]);

      }
      ///Makes a GaussianAtomicBasis from separate coords (a.u.)
      GaussianAtomicBasis(double x,double y, double z){
         AddParam(CENTERX,x);
         AddParam(CENTERY,y);
         AddParam(CENTERZ,z);
      }
      ///Prettier printing
      operator std::string()const;
};
/*
///Iterator to loop over all shells in a GaussianBasis
class ShellItr: public NestedPtrItr<
   BasisSet::iterator,BasisSet::iterator,BasisSet>{
   public:
      ///Constructor.  Begin==End makes iterator just past last shell
      ShellItr(BasisSet::iterator Begin,BasisSet::iterator End):
         NestedPtrItr<BasisSet::iterator,
             BasisSet::iterator,BasisSet>(Begin,End){}

};

///Iterator to loop over all primitives in a GaussianBasis
class PrimItr: public NestedPtrItr<
   ShellItr, GaussianAO::iterator,GaussianPrim>{
   public:
      ///Constructor. Begin==End makes iterator just past last primitive
      PrimItr(ShellItr Begin,ShellItr End):
         NestedPtrItr<ShellItr,GaussianAO::iterator,
                        GaussianPrim>(Begin,End){}

};*/

///The top level interface to a Gaussian basis set.
class GaussianBasis: public BasisSet{
   public:
      ///Standard clone function, calls copy constructor on this
      CLONE(GaussianBasis)

      ///Default constructor does nothing
      GaussianBasis(){}
      ///Returns the number of centers in the basis set
      size_t NAtoms()const{return BasisSet::NPrims(0);}
      ///Returns the number of shells in the basis set
      size_t NShells()const{return BasisSet::NPrims(1);}
      ///Returns the number of orbitals in the basis set
      size_t NOrbitals()const{return BasisSet::NPrims(2);}
      ///Returns the number of primitives in the basis set
      size_t NPrims()const{return BasisSet::NPrims(3);}

      /*//@{ Iterators
      ///Returns a modifiable iterator to the first shell in the basis set
      ShellItr ShellBegin(){return ShellItr(begin(),end());}
      ///Returns a modifiable iterator just past the last shell of the basis
      ShellItr ShellEnd(){return ShellItr(end(),end());}
      ///Returns a modifiable iterator to the first primitive in the basis set
      PrimItr PrimBegin(){return PrimItr(ShellBegin(),ShellEnd());}
      ///Returns a modifiable iterator just past the last shell of the basis
      PrimItr PrimEnd(){return PrimItr(ShellEnd(),ShellEnd());}
      ///@}*/
};

}//End PsiAPI namespace
#undef CLONE


#endif /* PSIAPI_GAUSSIANBASIS_H_ */
