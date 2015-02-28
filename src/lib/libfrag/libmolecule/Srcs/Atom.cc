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
#include <sstream>
#include "Atom.h"

namespace psi{
namespace LibMolecule{

Atom::Atom(const double* Carts,const int Z, const bool IsBohr,
      const std::string Label,  const double Charge, const int NElec,
      const double Mass):AtomGuts(Carts,Z,IsBohr,Label,Charge,NElec,Mass){
}

std::string Atom::PrintOut(const bool Bohr,
      const int DebugLevel)const{
   std::stringstream AtomStr;
   AtomStr<<this->Label();
   if(DebugLevel>=1)AtomStr<<
         (Bohr?" (a.u.):":" (Angstroms):");
   double Conv=(Bohr?1.0:1/AngToBohr());
   for(int i=0;i<3;i++)AtomStr<<" "<<Carts_[i]*Conv;
   if(DebugLevel>=2)AtomStr<<" m (Daltons)="<<Mass()
                           <<" q (a.u.)="<<Charge();
   if(DebugLevel>2)AtomStr<<" NE="<<NElec()
                          <<" Z="<<Z();
   return AtomStr.str();
}

bool Atom::operator==(const Atom& other)const{
   if(this->IsReal()&&!other.IsReal())return false;
   else if(this->IsPointCharge()&&!other.IsPointCharge())return false;
   else if(this->IsGhost()&&!other.IsGhost())return false;
   else if(this->IsDummy()&&!other.IsDummy())return false;
   if(this->Z()!=other.Z())return false;
   if(this->Mass()!=other.Mass())return false;
   for(int i=0;i<3;i++)
      if(fabs((*this)[i]-other[i])>1e-6)return false;
   return true;
}

}}//End namespaces

