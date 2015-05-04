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
#ifndef PSIAPI_ATOM_H_
#define PSIAPI_ATOM_H_
#include <set>
#include "PsiAPIBase.h"
#include "Labels.h"
namespace PsiAPI{
/** \brief A class to wrap atomic data
 *
 *   The atom class serves as an interface to atomic data. It is intended
 *   to be, for the most part a read-only class, which is to say that
 *   once the data is set you can't change it.  The notable exception to
 *   this is the label system.
 *
 *   The label system is meant to be a flexible way of "tagging" atoms.
 *   For example it allows you to specify that an atom is to be treated
 *   as a point charge or a ghost atom.  Some typical usages of the
 *   label system:
 *   \code
     //Fake function to get a molecule
     PsiAPI::Molecule AMol=FxnReturningMolecule();
     //Loop over the atoms in the molecule
     PsiAPI::AtomItr AtomI=AMol.AtomBegin(),AtomEnd=AMol.AtomEnd();
     for(;AtomI!=AtomEnd;++AtomI){
        //Check if atom is a ghost
        bool IsGhost=AtomI->Is(GHOST);
        //Loop over labels and make a list of all of them
        Atom::LabelItr LabelI=AtomI->LabelBegin(),LabelEnd=AtomI->LabelEnd();
        std::set<AtomLabels> ItsLabels;
        for(;LabelI!=LabelEnd;++LabelI)
           ItsLabels.insert(*LabelI);
     }
     \endcode
 *
 *   The labels in the system are enumerations to ensure consistent naming
 *   across the project and a comprehensive list is stored in the "Labels.h'
 *   file.
 *
 *   On an implementation note this class does not own its coordinates so
 *   they shouldn't be freed.  The coordinates are owned, in a block, by
 *   the Molecule class.  This also means an instance of an Atom is invalidated
 *   if the Molecule it belongs to goes out of scope.
 */
class Atom: public PsiAPIBase{
   private:
      ///Atomic Number
      size_t Z_;
      ///The atomic symbol
      std::string Symbol_;
      ///Pointer to the Cartesian coordinates (in a.u.) for the atom
      const double* Carts_;
      ///The charge of the atom
      double Q_;
      ///The mass of the atom
      double Mass_;
      ///Labels the atom possess
      std::set<AtomLabels> Labels_;

   public:
      /** \brief Basic constructor for an atom
       *
       *  Parameters should largely be self-explanatory, but:
       *
       *  \param[in] Z The atomic number, 1 is hydrogen, 2 is helium, etc.
       *               Use 0 for things where such a quantity is ill defined
       *  \param[in] Carts The Cartesian coordinates of the atom, in Bohr
       *  \param[in] Q The charge of the atom, in electrons, defaults to
       *               neutral
       *  \param[in] Mass The mass of the atom, any negative number causes
       *               the mass to be looked up by atomic number.  If supplied
       *               mass should be in Daltons, i.e. atomic mass units, i.e.
       *               the units where carbon 12 has a mass of exactly 12.0000
       *
       */
      Atom(size_t Z,const double *Carts, double Q=0, double Mass=-1);
      ///Returns the atomic number
      size_t Z()const{return Z_;}
      ///Returns the atomic symbol
      std::string AtomicSymbol()const{return Symbol_;}
      ///Returns the charge of the atom
      double Q()const{return Q_;}
      ///Returns the mass of the atom
      double Mass()const{return Mass_;}
      ///Returns coordinate \f$i\f$ where \f$i=0\f$ is \f$x\f$ etc.
      double operator[](size_t i)const{return Carts_[i];}

      ///@{
      /** Functions and types relating to labels*/
      typedef std::set<AtomLabels>::const_iterator LabelItr;
      ///Returns an iterator to the first Label
      LabelItr LabelBegin()const{return Labels_.begin();}
      ///Returns an iterator just past the last Label
      LabelItr LabelEnd()const{return Labels_.end();}
      ///Returns true if atom possess a given label
      bool Is(const AtomLabels& Label)const{
         return Labels_.find(Label)!=Labels_.end();
      }
      ///Adds a label
      void AddLabel(AtomLabels NewLabel){Labels_.insert(NewLabel);}
      ///@}
      ///Converts the atom to a string
      operator std::string() const;


};

}



#endif /* PSIAPI_ATOM_H_ */
