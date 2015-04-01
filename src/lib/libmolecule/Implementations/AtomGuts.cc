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

#include "AtomGuts.h"

namespace psi{
namespace LibMolecule{

AtomGuts::AtomGuts(const double* Carts,const int Z, const bool IsBohr,
      const std::string Label,  const double Charge, const int NElec,
      const double Mass):
           Carts_(3,0.0),Z_(Z<0?-1*Z:Z),Label_(Label),NElec_(NElec),
           Mass_(Mass),Charge_(Charge){
   for(int i=0;i<3;i++)Carts_[i]=Carts[i];
   if(!IsBohr){
      double Conv=AngToBohr();
      for(int i=0;i<3;i++)Carts_[i]*=Conv;
   }

   if(Label==""&&Z!=0)Label_=Sym(Z);
   if(NElec==0&&Z>=1)NElec_=Z;
   if(Mass==0.0&&Z>=1)Mass_=LibMoleculeBase::Mass(Z);
}

void AtomGuts::Copy(const AtomGuts& other){
   this->Carts_=other.Carts_;
   this->Mass_=other.Mass_;
   this->Charge_=other.Charge_;
   this->Z_=other.Z_;
   this->NElec_=other.NElec_;
   this->Label_=other.Label_;
}


}}//End namespaces

