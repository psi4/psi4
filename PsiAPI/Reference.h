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
#ifndef PSIAPI_REFERENCE_H_
#define PSIAPI_REFERENCE_H_
#include <boost/shared_ptr.hpp>
#include <boost/move/unique_ptr.hpp>
#include "Utils/PsiMap.h"
namespace PsiAPI{
class Molecule;
class Options;
class BasisSet;

/** \brief The recognized reference types
 *
 *  Similar to the atom labels in "Labels.h" these are the types of
 *  references we recognize.  When you add a reference please also add
 *  documentation including what quantities are available.
 *
 *  The reference types are:
 *  - NOREF Indicates that all you have to work with is a basis set, a
 *          Molecule class, and an options array.
 *  - HF    A reference including an alpha and beta density (if the two are
 *          equal they will point to the same memory address, i.e. the
 *          address of the first element of both density matrices will be
 *          the same).  Also includes the NOREF data.
 *  - MBE   This is the reference data for a many-body expansion. It includes
 *          a series of alpha and beta densities for the monomers and the
 *          energies of each calculation that has run.  For the density
 *          matrices, the order is the same as the fragments (i.e. the
 *          first density is for the first fragment etc.)
 */
enum RefType{NOREF,HF,MBE};

/** \brief The recognized quantities
 *
 *  This is a list of the quantities that may be stored in the Reference
 *  class and their descriptions:
 *
 *  - ALPHA_DENSITY The alpha density matrix
 *  - BETA_DENSITY The beta density matrix
 */
enum DataType{ALPHA_DENSITY,BETA_DENSITY};

/** \brief The class charged with holding restart information and enabling
 *         method to method communication.
 *
 *  If you know only one thing about this class know this: it owns the
 *  memory inside it.  You do not delete this class's memory.  You
 *  allocate memory and give it to this class and do NOT delete it. This
 *  is because the memory in this class will live on after your method
 *  has returned.  The Reference class will never change your data, except
 *  when an instance goes out of scope, at which point it will free it.
 *
 *  The overall restart strategy works like this, you read and save all
 *  starting/intermediate/ and final results to this class.  That way
 *  all of the information you need is always in this class.  When your
 *  method is called you will get an instance of this class, your method
 *  should see what data is inside that instance and choose an appropriate
 *  starting point within the algorithm.  As the algorithm progresses
 *  overwrite the data inside the class such that if you were handed
 *  this object, you would be able to start at your current point again.
 *
 *  The above strategy poses several problems.  First this class needs
 *  to contain all the information to get any calculation started.  Second,
 *  it needs to be able to store any essential intermediate information
 *  along the way, from a myriad of possible methods.  Third, it needs
 *  to address memory ownership.
 *
 *  Our current philosophy is that any method should be capable of starting
 *  minimally from a Reference of type NOREF, which contains the system,
 *  the basis set, and a list of options.  In the initial phase of the
 *  implementation of this API that is the only reference anything will
 *  recognize.  As new methods are added they should define what they need
 *  to restart, and actually store that information in this class.
 *
 *  Every instance of this class contains the system that was currently
 *  being worked on, the options used to run that calculation, and the
 *  basis set.  As for the other contents of this class the main criterion
 *  is that the results be small-ish and meaningful. Anything that can
 *  easily be recomputed, should be recomputed. We provide
 *  a vector of double* for storing multiple arrays.  Exactly what is in
 *  these arrays will vary from method to method.  The label on the
 *  Reference will tell us how to read the contents.
 *
 *  To illustrate some of these ideas consider an SCF routine.  The entire
 *  SCF procedure is specified by the basis set (which in turn gives the
 *  integrals) and the density (which gives our initial guess for each
 *  iteration).  On iteration one we read the density and basis in from this
 *  class, we then perform an SCF iteration, resulting in a new density,
 *  we then save that density to this class, and repeat.
 *  If at any time we encounter an error we throw an exception.
 *  The driver will then see the exception, grab the current Reference,
 *  save it, and exit gracefully.
 *
 *
 *  Presently the reference object does not traverse different packages, that
 *  is to say we are not running a Psi4 HF and then feeding it into CFour.
 *  The latter functionality is intended at a later date as this then gives
 *  a way of passing results from one package to another; however there
 *  is a notable difficulty in that quantities are not standardized yet so
 *  the order of say the density matrices may differ from package to package.
 *
 *  When we progress to the stage of using the reference as a package
 *  traversal object it will make sense to create an object capable of
 *  reorienting a quantity to an arbitrary ordering.
 *
 *  Implementation note: in order to allow for multiple instances of the
 *  same DataType (e.g. two alpha densities) I currently use a (wrapped)
 *  std::map
 *  between the DataType and a vector of double*'s.  This is in essence
 *  a ghetto multimap.  Until C++11 is universally used std::multimap is
 *  not guaranteed to be stable (i.e. if I have two alpha density matrices
 *  the one I insert first is not guaranteed to always be the first one).
 *  The use of the vector ensures my stability.
 *
 *  Typical usage scenarios of the Reference class:
 *  \code
 *  //Assume I already have a Reference (driver will give you one)
 *  Reference ARef;
 *
 *  //Get it's type
 *  RefType ItsType=ARef.Type();
 *
 *  //Because you know the type you know what's in it, but pretending
 *  //you don't....
 *
 *  //Check if it has an ALPHA_DENSITY
 *  bool HasRhoA=(ARef.NData(ALPHA_DENSITY)!=0);
 *
 *  //If it doesn't have one, make one (nbasis=number of basis functions)
 *  ARef.AddData(ALPHA_DENSITY,new double[nbasis*nbasis]);
 *
 *  //Or (not both) if we already had an existing chunk of memory (say
 *  //from our own Matrix/Tensor class):
 *  double* ExistingRhoA;
 *  ARef.AddData(ALPHA_DENSITY,ExistingRhoA);
 *
 *  //Note:ARef is now coupled to ExistingRhoA so...
 *  ExistingRhoA[3]=3.14;
 *
 *  //...also changes ARef, hence this should be true
 *  bool IsSame=(ExistingRhoA[3]==ARef.DataI(ALPHA_DENSITY)[3]);
 *
 *  //You can also modify ARef in place:
 *  ARef.DataI(ALPHA_DENSITY)[3]=3.14;
 *
 *  \endcode
 */
class Reference{
   private:
      ///Convenient typedef of data type could be made complex
      typedef boost::move::unique_ptr<double *> Data_t;
      ///The "real" data lives here
      psi::PsiMap<DataType,std::vector<Data_t> > Data_;
      ///The molecule used to generate this data
      boost::shared_ptr<Molecule> Mol_;
      ///The basis set used to generate this data
      boost::shared_ptr<BasisSet> Basis_;
      ///The Options used to generate this data
      boost::shared_ptr<Options> Options_;
      ///What sort of reference is this?
      RefType MyType_;
   public:
      ///Creates a reference of type Ref
      Reference(RefType Ref=NOREF):MyType_(Ref){}

      ///Returns the number of data pieces of a given type
      size_t NData(DataType DType)const{return Data_[DType].size();}

      /** \brief Returns the i-th (counting from 0) data piece of a given
       *         type
       *
       *  This call is for returning an existing array only!!!
       *  To add data use the AddData member function.  Don't delete the
       *  pointer (that's what the const is there for)
       */
      const double*  DataI(DataType DType,size_t i=0)const{return Data_[DType][i].get();}

      /** \brief Returns the i-th (counting from 0) data piece of a given
       *         type
       *
       *  This is identical to the previous call except now you are
       *  allowed to modify the data.  The intent is to allow you to
       *  build inside this object not for you to delete or allocate memory.
       *  If you don't want some of the objects in the Reference instance
       *  you were given make a new Reference and let the old one go out
       *  of scope.
       */
      double* DataI(DataType DType,size_t i=0){return Data_[DType][i].get();}

      ///Creates a new array of type "DType" at memory address "Address"
      void AddData(DataType DType,const double*& Address){
         Data_[DType].push_back(Data_t(Address));
      }

      ///Returns the type of this reference
      RefType Type()const{return MyType_;}

      ///Returns the molecule
      const Molecule& Molecule()const{return *Mol_;}
      ///Returns the basis set
      const BasisSet& Basis()const{return *Basis_;}
      ///Returns the options
      const Options& Options()const{return *Options_;}
      ///Allows modification of the options, through the reference
      Options& Options(){return *Options_;}

};

}//End PsiAPI namespace

#endif /* PSIAPI_REFERENCE_H_ */
