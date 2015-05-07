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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_ATOMGUTS_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_ATOMGUTS_H_
#include <vector>
#include "../LibMoleculeBase.h"
namespace psi{
namespace LibMolecule{

///Allows easy editing of the type
typedef std::vector<double> Carts_t;

class AtomGuts:public LibMoleculeBase{
   protected:
      ///Deep copy of da atom
      void Copy(const AtomGuts& other);
      ///The three Cartesian coordinates (a.u.) in the order: x,y,z
      Carts_t Carts_;
      ///The atomic number
      int Z_;
      ///The chemical symbol
      std::string Label_;
      ///The number of electrons
      int NElec_;
      ///The mass in Daltons
      double Mass_;
      ///The Charge in multiples of the electron's charge
      double Charge_;
   public:
      AtomGuts(const double* Carts, const int Z=0,const bool IsBohr=true,
            const std::string Label="",const double Charge=0.0,
           const int NElec=0,const double Mass=0.0);
      AtomGuts(const AtomGuts& other){this->Copy(other);}
      const AtomGuts& operator=(const AtomGuts& other){
         if(this!=&other)this->Copy(other);
         return *this;
      }
};

}}//End namespaces
#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_ATOMGUTS_H_ */
