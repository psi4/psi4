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
#include<sstream>
#include "ParameterType.h"
#include "exception.h"

namespace psi{
namespace LibMolecule{

ParamT::ParamT(const std::string& BaseAbbrv,
               const std::string& BaseName,
               const int Order,
               const int Priority):
    Order_(Order),Priority_(Priority){
  Base_[0]=BaseAbbrv;
  Base_[1]=BaseName;
  Atom_[0]=BaseAbbrv;
  Atom_[1]=BaseName;
  Full_[0]=Name(0);
  Full_[1]=Name(1);
}

ParamT::ParamT(const ParamT& Other,
       const std::string& BaseAbbrv,
       const std::string& BaseName,
       const int Order,
       const int Priority):
      Order_(Order),Priority_(Priority){
      Base_[0]=BaseAbbrv;
      Base_[1]=BaseName;
      Atom_[0]=Other.Atom()[0];
      Atom_[1]=Other.Atom()[1];
      Full_[0]=Name(0);
      Full_[1]=Name(1);
}

static std::string OrderLookUp(const size_t Prior);

std::string ParamT::Name(const size_t i)const{
   std::stringstream Mess,Base,Atom;
   bool HasAtom=(Atom_[i]!=Base_[i]),
        HasOrder=Order_>=0,
        HasPrior(Priority_>0);
   if(i==1)
      Base<<OrderLookUp(HasOrder?Order_:0)<<(HasOrder?" ":"");
   Base<<Base_[i];
   if (i==0&&HasOrder)
      Base<<Order_;
   if(HasAtom){
      Atom<<Atom_[i];
      if(HasPrior)Atom<<Priority_;
   }
   Mess<<(i==0?Atom.str():Base.str())
       <<(i==0?"_":" ")
       <<(i==0?Base.str():Atom.str());
   return Mess.str();
}

std::string OrderLookUp(const size_t Prior){
   std::string Result;
   switch(Prior){
      case(0):{Result="";break;}
      case(1):{Result="Primary";break;}
      case(2):{Result="Secondary";break;}
      case(3):{Result="Tertiary";break;}
      case(4):{Result="Quaternary";break;}
      case(5):{Result="Quinary";break;}
      case(6):{Result="Senary";break;}
      case(7):{Result="Septenary";break;}
      case(8):{Result="Octonary";break;}
      default:{
         throw PSIEXCEPTION("Precedence Ordinal Numbers Beyond 8 are not"
               " coded");
       break;
      }
   }
   return Result;
}

}}

