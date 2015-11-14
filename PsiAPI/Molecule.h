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
#ifndef PSIAPI_MOLECULE_H_
#define PSIAPI_MOLECULE_H_
#include <vector>
#include <boost/shared_ptr.hpp>
#include "Labels.h"
#include "PsiAPIBase.h"
namespace PsiAPI{
class Atom;
/** \brief A container to hold a molecule
 *
 *
 *  For our purposes a molecule is any set of atoms. This means we are
 *  terming a water cluster a molecule, as well as the waters within that
 *  cluster a molecule.  Although possibly confusing in theory, this rarely
 *  presents a problem in practice, in fact the ability to treat any set
 *  of atoms the same way is often a superior code design strategy to
 *  distinguishing between various types of molecules (e.g. fragments, QM
 *  system, MM system, etc.).
 *
 *  I've made a couple of design decisions presently.  First, the coordinates
 *  sit in a buffer that is organized so that coordinates,
 *  \f$i,i+1, i+2\f$ are the \f$x,y,\f$ and \f$z\f$ components of atom
 *  \f$i\f$.  Each atom then contains a pointer to its coordinates, instead
 *  of its actual coordinates.  This is useful for rotating, translating,
 *  etc.  the entire molecule easily, as well as saving space when atoms
 *  are replicated many times.  The latter is particularly important for
 *  methods that require many copies of the molecule or parts of the
 *  molecule (e.g. fragment based methods).
 *
 *  The above brings us to the next important design decision, Atom, Fragment,
 *  UnitCell, etc. classes are associated with a molecule.  Atom in particular
 *  is only good as long as the Molecule instance it came from is still in
 *  scope.  Fragment and UnitCell instances contain shared pointers to
 *  the original Molecule and are thus good for as long as they are in scope.
 *  The idea being that Fragment and UnitCell instances are likely to be
 *  used as Molecule instances themselves, whereas Atom instances are a means to
 *  obtain data.  This means letting a Molecule instance go out of scope
 *  does not invalidate a Fragment or UnitCell instance, but does
 *  invalidate an Atom instance.  If you want to keep the data of an Atom
 *  instance, but are going to let the Molecule instance go out of scope,
 *  you need to copy the data.
 *
 *  Related to this last design decision, the carts of the molecule are
 *  kept in an std::vector and the Atom instances are linked to memory
 *  addresses in that vector.  Unfortunately, this means that any time
 *  the vector reallocates, all of the references are invalidated.  This
 *  in turn means all of the atoms need to be updated with their new
 *  locations and consequentially invalidates all of the references to
 *  the Atom instances.  This is a very technical way of saying that you
 *  should fill a Molecule before you use it.  If for some reason you
 *  absolutely can not fill a Molecule before you use it, the constructor
 *  takes a parameter for the number of atoms.  This will ensure your
 *  atom references remain valid until you insert one atom more than
 *  this parameter.  Consequentially, specification of this parameter,
 *  even when you plan to fully fill the class first, will lead to
 *  optimum performance, but is not necessary.
 *
 *
 *  As another design decision I encourage the use of the iterators, but
 *  I also realize that it is sometimes more useful to be able to pluck
 *  out a particular atom.  For this reason, I have allowed element-wise
 *  access of atoms by index, i.e. you can request the third atom.  As
 *  usual, counting starts from 0.
 *
 *
 *  The molecule class is intended to be initialized as:
    \code
    //Assume we have a list of coordinates...
    double* Carts;
    //...and a list of atomic numbers...
    size_t* Zs;
    //The number of atoms being:
    size_t NAtoms;

    //First we declare a molecule (use of NAtoms leads to better performance)
    PsiAPI::Molecule MyMol(NAtoms);
    //...or (not both)
    PsiAPI::Molecule MyMol;
    //Now we fill it
    for(size_t AtomI=0;AtomI<NAtoms;++AtomI){
              size_t Offset=3*AtomI;
       MyMol.AddAtom(Zs[AtomI],&Carts[Offset]);
       //...or (again not both):
       MyMol.AddAtom(Zs[AtomI],Carts[Offset],Carts[Offset+1],Carts[Offset+2]);
       //Say we want to make the fourth atom a ghost
       if(AtomI==3){
          Atom& LastAtom=MyMol.Last();
          LastAtom->AddLabel(GHOST);
       }
    }
    \endcode
*
*   The above code represents a basic filling of the molecule.  For more
*   advanced fillings, such as those used in automatic counterpoise
*   corrections see the Fragment class.  For typical users the process
*   of filling a molecule is irrelevant and it is only its use that matters.
*   In this case the code snippet above should be enough to illustrate how
*   to loop over the atoms.
*/
class Molecule: public PsiAPIBase{
   private:
      ///Internal typedef of pointer to an atom
      typedef boost::shared_ptr<Atom> SharedAtom_t;
      ///The carts of the system
      std::vector<double> Carts_;
      ///A vector of atoms
      std::vector<SharedAtom_t> Atoms_;
      ///The charge of the molecule
      double Charge_;
      ///The multiplicity of the molecule
      size_t Mult_;
   public:
      /* Creates a molecule with charge "Charge", multiplicity
       * "Mult" and "NAtoms" atoms.
       */
      Molecule(double Charge=0.0, size_t Mult=1,size_t NAtoms=0);

      ///Returns the multiplicity of the molecule
      size_t Mult()const{return Mult_;}
      ///Returns the charge of the molecule
      double Charge()const{return Charge_;}


      ///@{
      /** \brief Functions for adding atoms
       *
       *  For convenience I have added several overloaded functions
       *  that will add an atom.  They all do the same thing, but are
       *  set up to be more useful depending on how the coordinates
       *  are stored.
       *
       *  \param[in] Z The atomic number.  1 is hydrogen, 2 is helium, etc.
       *               Feel free to use 0 for a dummy atom/point charge or
       *               anything else where such a field is not relevant.
       *  \param[in] Carts The coordinates of the atom, in a.u. The first
       *               Carts[0] should be the x coordinate, Carts[1] y, and
       *               Carts[2] z.
       *  \param[in] Q The charge of the atom, defaults to 0, i.e. the atom
       *               is neutral.  This really only plays a role in theories
       *               where it makes sense to talk about the charge on an atom
       *               or when this atom is tagged as a point charge.  For other
       *               theories, it's the charge of the molecule that
       *               matters.
       *  \param[in] M The mass of the atom, defaults to -1.  Any negative
       *               value will cause the mass to be looked up for you. A
       *               value of 0 will be assigned to anything with Z=0.  If
       *               you choose to specify the mass it should be in Daltons,
       *               i.e. atomic mass units, i.e. the unit where carbon 12
       *               weighs exactly 12.000....
       */
      void AddAtom(size_t Z, double x, double y, double z, double Q=0, double M=-1);
      void AddAtom(size_t Z,const double* Carts,double Q=0, double M=-1){
         this->AddAtom(Z,Carts[0],Carts[1],Carts[2],Q,M);
      }
      ///@}

      ///@{
      ///Returns the most recently added atom, useful for initial filling
      Atom& Last(){return *(Atoms_.back());}
      const Atom& Last()const{return *(Atoms_.back());}
      ///@}

      ///@{
      ///Returns atom i
      Atom& operator[](size_t i){return *(Atoms_[i]);}
      ///Returns a pointer to a const atom i
      const Atom& operator[](size_t i)const {return *(Atoms_[i]);}
      ///@}

      ///Returns the number of atoms in the molecule
      size_t NAtoms()const{return Atoms_.size();}

      ///@{
      ///Typedefs of the iterators to the atoms
      typedef std::vector<SharedAtom_t>::iterator AtomItr;
      typedef std::vector<SharedAtom_t>::const_iterator cAtomItr;
      ///@}

      ///@{
      ///Functions to return (const) iterators to the first/last atom
      AtomItr AtomBegin(){return Atoms_.begin();}
      cAtomItr AtomBegin()const{return Atoms_.begin();}
      AtomItr AtomEnd(){return Atoms_.end();}
      cAtomItr AtomEnd()const{return Atoms_.end();}
      ///@}

      ///Function to turn this molecule into a string
      operator std::string()const;
};


}



#endif /* PSIAPI_MOLECULE_H_ */
