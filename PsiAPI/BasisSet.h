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
#ifndef PSIAPI_BASISSET_H_
#define PSIAPI_BASISSET_H_
#include <sstream>
#include <vector>
#include "boost/shared_ptr.hpp"
#include "Utils/PsiMap.h"
#include "PsiAPIBase.h"

namespace PsiAPI{
/** \page page1 PsiAPI Basis Set Philosophy
 *
 * \tableofcontents
 *
 * \section sec1 Goal
 * The goal of the PsiAPI BasisSet class is to be capable of uniformly
 * handling any type of basis set that may be encountered in quantum
 * chemistry including: Gaussian, Slater, plane wave, Dirac, etc. In
 * order to do this we distill each basis set type down to its core elements
 * and propose a common nomenclature.  The next section refreshes my (and
 * hopefully your) memory regarding basis sets.
 *
 * \section sec2 Background
 *
 * \subsection sec21 Gaussian Basis Sets
 * Perhaps the most ubiquitously found basis set in quantum chemistry is
 * the Gaussian basis.  At the most primitive level, a Gaussian basis is
 * made up of "primitives" (how original...), the \f$i\f$ -th one being
 * denoted \f$\Theta_{Prim}^i\f$.  \f$\Theta_{Prim}^i\f$ is centered on
 * \f$r_i\f$ and has the form:
   \f[
      \Theta_{GPrim}^i(r)=(x-x_i)^{l_x}(y-y_i)^{l_y}(z-z_i)^{l_z}
      e^{-\alpha {|r-r_i|}^{2}}.
   \f]
 * Here \f$\alpha\f$ is known as the "exponent" and \f$l\f$ is the angular
 * momentum of the shell (broken into its components).
 *
 * Anyways, an individual Gaussian primitive is a poor approximation
 * to an actual atomic orbital (AO) so one usually uses a linear combination
 * of Gaussian primitives to approximate a Gaussian AO.  The resulting
 * basis function is known as a contracted Gaussian AO,
 * \f$\Theta_{GAO}^i\f$, and given by:
 * \f{eqnarray*}
 *   \Theta_{GAO}^i(r)=&
 *    \sum_{j=1}^NC_j(x-x_i)^{l_x}(y-y_i)^{l_y}(z-z_i)^{l_z}
 *     e^{-\alpha_j {|r-r_i|}^{2}}\\
 *    =&
 *    (x-x_i)^{l_x}(y-y_i)^{l_y}(z-z_i)^{l_z}\sum_{j=1}^NC_j
 *    e^{-\alpha_j {|r-r_i|}^{2}}
 * \f}.
 * The \f$C_j\f$ are known as contraction coefficients.
 *
 * In this sense,
 * Gaussian primitives form a basis for a Gaussian AOs.  These
 * AOs are then grouped into sets called shells, where any particular shell
 * is made up of the AOs on the same center with the same set of \f$C_j\f$
 * and \f$\alpha_j\f$ parameters.  Hence the Gaussian AOs form a basis
 * for a Gaussian shell.  Shells come in two flavors: "Cartesian" or
 * "Pure", the distinguishing fact for our purpose being that Cartesian
 * shells contain \f$\frac{(l+1)(l+2)}{2}\f$ AOs and pure contain \f$2l+1\f$.
 * The set of shells on a common center are part of a set which
 * we call the atomic basis set and finally the set of atomic basis sets
 * form a basis for the molecule, which is the basis set we ultimately care
 * about.
 *
 * \subsection sec22 Slater Basis Sets
 *
 * Closely related to the Gaussian primitive is the Slater function,
 * \f$\Theta_{Slater}^i\f$, given by:
 *
 * \f[
 *   \Theta_{Slater}^i(r)=(x-x_i)^{l_x}(y-y_i)^{l_y}(z-z_i)^{l_z}
 *   e^{-\alpha |r-r_i|},
 * \f]
 * the only real difference being the power of the distance in the exponent
 * of \f$e\f$ and the fact that a Slater AO is a single Slater function not
 * a contraction.  Slater AOs are then ordered into shells and then atom
 * basis sets in the same way that Gaussian basis sets are.
 *
 * \section sec3 Current Implementation and History
 *
 * The previous section suggests that we really only
 * need two classes a basis function and a basis set.  Shells and all that jazz can be written
 * in terms of these two classes.  This dramatically simplifies our
 * implementation aside from a minor hiccup in that in practice there is
 * no distinction between an AO and a shell.  That is one typically stores
 * the AOs plus a flag for whether they are Pure or not and the total
 * angular momentum of the shell.  The individual AOs are not stored, i.e.
 * we have the exponents and coefficients for any "p"-like AO and also
 * the flag for whether the "p" shell is pure or not, (we also obviously
 * have the angular momentum because we know it's a p shell), but we
 * never explicitly store that say orbital 1 is a "px" AO.  Integral packages
 * expect shell pairs or shell quartets so we can't really get away from
 * this form of storage unless the interfaces to integral packages change.
 * This also means the order depends on the integral package.
 *
 * Largely for my own recollection I have documented the history of this
 * class.  Hopefully this way I will remember how I got to where I did and
 * why.  Everyone else can skip to the last paragraph of this section to
 * get a breakdown of how the class currently works.
 *
 * My original plan was to use the abstraction of the typical basis set to
 * implement the desired functionality as recursive calls to a templated
 * BasisSet class.  This solution avoids re-implementing much of the
 * inner workings of each layer of a traditional Gaussian basis set; however it means that
 * there is no common base class other than PsiAPIBase for the various basis
 * sets.  Specifically A GaussianBasis has a base class of type
 * BasisSet<GaussianAtomicBasis>, plane-wave and Dirac basis sets
 * both have a base type of BasisSet<BasisFxn>.
 * The problem is we have realized that a BasisSet can be a type of
 * BasisFxn, hence my next try was to make BasisSet derive from BasisFxn,
 * then all the classes should have a common base class of BasisFxn.
 *
 * This try likely could have been made to work, but with some slightly
 * illogical code, e.g. a basis function would need to contain an array of
 * basis functions, which actually makes it a basis set...This suggests
 * a realization that a basis function is actually a basis set.  I think
 * that what I have really stumbled upon is a sort of trivial identity,
 * e.g. a basis set for representing say a unit vector in the x-direction
 * is a unit vector in the x-direction.
 *
 *
 *
 * This leads to a usage of:
 * \code
   //Assume we have a basis (there's one in your Ref)
   BasisSet MyBasis;

   //Get its type
   BasisType Type=MyBasis.Type();
   if(Type!=GAUSSIAN)
     PSIERROR("How is our Gaussian Example not of type GAUSSIAN?");

   //Get the number of atoms it is defined for
   size_t NAtoms=MyBasis.NBasis();
   for(size_t AtomI=0;AtomI<NAtoms;++AtomI){

      //Get the number of shells on AtomI
      size_t NShells=MyBasis[AtomI].NBasis();
      for(size_t ShellI=0;ShellI<NShells;++ShellI){

        //Get the number of orbitals in ShellI
        size_t NOrbs=MyBasis[AtomI][ShellI].NBasis();

        //We could then loop over the orbitals like this
        for(size_t OrbI=0;OrbI<NOrbs;OrbI){
           //Do stuff that depends on orbitals, like contract densities
        }// End loop over orbitals

        //As mentioned what we usually pass are shell pairs/quartets
        //So get the number of primitives in current AO
        size_t NPrims=MyBasis[AtomI][ShellI].GaussianAO::NBasis();
        for(size_t PrimI=0;PrimI<NPrims;++PrimI){

           //Get the exponent of PrimI
           double alpha=MyBasis[AtomI][ShellI][PrimI].GetParam[EXPONENT];
           //Get the contraction coefficient of PrimI
           double ci=MyBasis[AtomI][ShellI][PrimI].GetParam[COEF];

        }//End Primitive Loop
      }//End Shell Loop
   }//End Atom Loop
   \endcode

   I realize most people probably don't like the idea that the number of
   atoms, shells, orbitals, and primitives all stem from a call of the same
   name so
   the GaussianBasis class includes wrappers that are appropriately named
   and return the same quantities.  Since it only makes sense to have
   such calls if the BasisSet is truly of GaussianBasis type, you will
   have to cast up using a dynamic cast, which I don't feel bad about
   because you don't need to do it.

   The above syntax is overly verbose and would need to be copy/pasted
   every time we want to access a primitive Gaussian function or a shell.
   For this reason I have made shell and primitive iterators. To use them:
   \code
   //Get a GaussianBasis
   GausianBasis MyBasis;
   ShellItr ShellI=MyBasis.begin(),ShellEnd=MyBasis.end();
   for(;ShellI!=ShellEnd;++ShellI){

      //The actual shell that is ShellI
      GaussianShell Shell=(*ShellI);

      //Is this shell pure
      bool IsPure Shell.IsPure();

      //What's its angular momentum
      size_t L=Shell->L();

      //Number of orbitals in ShellI
      size_t Norbs=ShellI->NBasis();

      //Number of primitives in ShellI
      size_t NPrims=ShellI->GaussianPrimitive::NBasis();
   }//End loop over shells

   //To iterate over primitives
   PrimItr PrimI=MyBasis.begin(),PrimEnd=MyBasis.end();
   for(;PrimI!=PrimEnd;++PrimI){

      //Get the exponent
      double alpha=PrimI->GetParam(EXPONENT);

      //Get the contraction coefficient
      double ci=PrimI->GetParam(COEF);
   }//End loop over primitives
   \endcode


 */


///The recognized basis set types
enum BasisType{SLATER,GAUSSIAN,PLANEWAVE,DIRAC};
/** \brief The types of the parameters in a basis function
 *
 *   These are all of the parameters that can be stored in a BasisFxn.
 *   Note they are all stored as doubles and cast for wrapper functions.
 *
 *   - EXPONENT The exponent of either a Gaussian or a Slater
 *   - COEF The contraction coefficient in a contracted Gaussian shell
 *   - ANGMOM The total angular momentum of the shell
 *   - PURE Whether the Gaussian is pure or Cartesian 1 (true) is pure
 *   - CENTERX The x coordinate the basis function is centered on
 *   - CENTERY The y coordinate the basis function is centered on
 *   - CENTERZ The z coordinate the basis function is centered on
 *
 */
enum ParamType{EXPONENT,COEF,ANGMOM,PURE,CENTERX,CENTERY,CENTERZ};

/** \brief The base class for a basis set.
 *
 *  You will get a pointer to this class in your reference.
 */
class BasisSet: public PsiAPIBase{
   private:
      ///Internal typedef of the array storing the basis functions
      typedef std::vector<boost::shared_ptr<BasisSet> > Fxn_t;
      ///The set of basis functions of type T
      Fxn_t Fxns_;
      ///Where the parameters relevant to this basis set are held
      PsiMap<ParamType,double> Params_;
   public:

      ///\name Constructing Members
      ///@{
      ///Members related to making a new BasisSet

      ///Standard clone function, calls copy constructor on this
      virtual BasisSet* Clone()const{return new BasisSet(*this);}
      ///Makes this BasisSet a shallow copy of other
      BasisSet(const BasisSet& Other):Params_(Other.Params_),Fxns_(Other.Fxns_){}
      ///Default constructor
      BasisSet(){}
      virtual ~BasisSet(){}
      ///@}

      ///@{ Parameter Access
      ///Returns the number of parameters
      size_t NParams()const{return Params_.size();}
      ///Returns the parameter with type "Type"
      double GetParam(ParamType Type)const{return Params_[Type];}
      ///Adds a parameter with type "Type" and value "Value"
      void AddParam(ParamType Type,double Value){Params_[Type]=Value;}
      ///@}

      ///\name Iterator Stuff for basis functions
      ///@{
      ///Types and class members related to iterating over the basis functions

      ///Typedef to a modifiable iterator to this basis
      typedef Fxn_t::iterator iterator;
      ///Typedef of a const iterator to this basis
      typedef Fxn_t::const_iterator const_iterator;
      ///Returns a modifiable iterator to the beginning of the BasisSet
      iterator begin(){return Fxns_.begin();}
      ///Returns a modifiable iterator just past the end of the BasisSet
      iterator end(){return Fxns_.end();}
      ///Returns a const iterator to the beginning of the BasisSet
      const_iterator begin()const{return Fxns_.begin();}
      ///Returns a const iterator just past end of the BasisSet
      const_iterator end()const{return Fxns_.end();}
      ///@}

      ///\name Member-wise access/insertion
      ///@{
      ///Class members to modify the basis functions

      ///Returns the number of BasisFxn
      virtual size_t NBasis()const{return Fxns_.size();}
      ///Returns a const BasisSet i
      const BasisSet& operator[](size_t i)const{return *(Fxns_[i]);}
      ///Returns a modifiable BasisSet i
      BasisSet& operator[](size_t i){return *(Fxns_[i]);}
      ///Copies object into this object via its Clone() method
      void insert(const BasisSet& object);
      ///@}

      /** \brief Function for recursively finding the number of basis fxns
       *
       *  You probably shouldn't ever use this function, but it needs to
       *  be public so that other classes can call it.  Basically, this
       *  function allows you to recursively determine the number of
       *  basis functions in a basis set at a recursion level of "Rec"
       *  below your current class.  Why is this useful?  Say you have an
       *  instance of
       *  the GaussianBasis class and you want to know the number of
       *  shells in the total basis set.  Ignoring the NShells member,
       *  because it wraps this call for you, you would have to loop over
       *  the atoms in the GaussianBasis and then the shells on each atom.
       *  This call will do that for you, if you give it "Rec=1".  The
       *  first level of recursion (Rec=0) gives you the number of basis
       *  functions in the current object (number of atoms in the previous
       *  example), the second gives you the total number of basis functions
       *  in each basis function of the current object (shells in previous
       *  example).  So on and so forth.  Again you can probably ignore
       *  this function as I believe I have made convenience functions
       *  for all useful quantities...
       */
      size_t NPrims(size_t Rec)const;


      ///Pretty printing of the BasisSet
      operator std::string()const;
};//End BasisSet declaration



}//End PsiAPI namespace



#endif /* PSIAPI_BASISSET_H_ */
