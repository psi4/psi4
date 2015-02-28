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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_MOLECULEGUTS_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_MOLECULEGUTS_H_
#include "../LibMoleculeBase.h"
#include "MolItr.h"
namespace psi{
namespace LibMolecule{
class Atom;
/** \brief A class designed to hold the details of the Molecule class
 *
 */
class MoleculeGuts:public LibMoleculeBase{
   private:
      ///Performs a deep copy
      void Copy(const MoleculeGuts& other);
   protected:
      ///The multiplicity
      int Mult_;
      ///The charge
      int Charge_;

      ///The actual atoms contained in this Molecule
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
      virtual MolItr Begin()const;
      virtual MolItr End()const;
      void SetCharge(const int Charge){Charge_=Charge;}
      void SetMult(const int Mult){Mult_=Mult;}
      ///Go away compiler warning...
      virtual ~MoleculeGuts(){}
      MoleculeGuts(const int Charge=0,const int Mult=1):
         Charge_(Charge),Mult_(Mult){}
      MoleculeGuts(const MoleculeGuts& other){
         this->Copy(other);
      }
      const MoleculeGuts& operator=(const MoleculeGuts& other){
         if(this!=&other)this->Copy(other);
         return *this;
      }
};


}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_MOLECULEGUTS_H_ */
