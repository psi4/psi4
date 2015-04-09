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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_PRIMITIVENODES_HH_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_PRIMITIVENODES_HH_

#include "MMAtomTypes.h"
#include "TrivialNodes.hh"
#include "RadialSearch.h"

namespace psi{
namespace LibMolecule{
#define PRIM(Name,Bonds,Abbv,Full,Centers...) \
   class Name:public RadialSearch<Bonds,Centers>{\
      private:\
         typedef RadialSearch<Bonds,Centers> Base_t;\
      public:\
         Name():Base_t(Abbv,Full){}\
};

/*********Carbon Primitives********************/
PRIM(Methane,4,"C","Methane",Carbon,Hydrogen,Hydrogen,Hydrogen,Hydrogen)
PRIM(Methyl,4,"C","Methyl",Carbon,Hydrogen,Hydrogen,Hydrogen)
PRIM(Methene,4,"C","Methene",Carbon,Hydrogen,Hydrogen)
PRIM(Methyne,4,"C","Methyne",Carbon,Hydrogen)
PRIM(C4,4,"C","Carbon",Carbon)
PRIM(Alkenyl1,3,"CDB","Alkenyl",Carbon,Hydrogen,Hydrogen)
PRIM(Alkenyl2,3,"CDB","Alkenyl",Carbon,Hydrogen)
PRIM(Alkenyl3,3,"CDB","Alkenyl",Carbon)
PRIM(Alkynyl1,2,"CTB","Alkynyl",Carbon,Hydrogen)
PRIM(Alkynyl2,2,"CTB","Alkynyl",Carbon)

class Alkenyl:public Node{
   public:
      Alkenyl():Node("CDB*","Alkenyl WildCard"){
         this->MyTypes_[0]=ParamT("CDB1","Alkenyl",1,0);
         this->MyTypes_.push_back(ParamT("CDB2","Alkenyl",2,0));
         this->MyTypes_.push_back(ParamT("CDB3","Alkenyl",3,0));
      }
};

class Alkynyl:public Node{
   public:
      Alkynyl():Node("CTB*","Alkynyl WildCard"){
         this->MyTypes_[0]=ParamT("CTB1","Alkynyl",1,0);
         this->MyTypes_.push_back(ParamT("CTB2","Alkenyl",2,0));
      }
};

/**********************Nitrogen Primitives*******************/
PRIM(Ammonium,4,"NP","Ammonium",Nitrogen,Hydrogen,Hydrogen,Hydrogen,Hydrogen)
PRIM(Ammonium1,4,"NP","Ammonium",Nitrogen,Hydrogen,Hydrogen,Hydrogen)
PRIM(Ammonium2,4,"NP","Ammonium",Nitrogen,Hydrogen,Hydrogen)
PRIM(Ammonium3,4,"NP","Ammonium",Nitrogen,Hydrogen)
PRIM(Ammonium4,4,"NP","Ammonium",Nitrogen)
PRIM(Ammonia,3,"N","Ammonia",Nitrogen,Hydrogen,Hydrogen,Hydrogen)
PRIM(Amine1,3,"N","Amine",Nitrogen,Hydrogen,Hydrogen)
PRIM(Amine2,3,"N","Amine",Nitrogen,Hydrogen)
PRIM(Amine3,3,"N","Amine",Nitrogen)
PRIM(Azo1,2,"NDB","Azo",Nitrogen,Hydrogen)
PRIM(Azo2,2,"NDB","Azo",Nitrogen)
PRIM(NTB,1,"NTB","Triple-Bonded Nitrogen",Nitrogen)

class Amine:public Node{
   public:
      Amine():Node("N*","Amine WildCard"){
         this->MyTypes_[0]=ParamT("N1","Amine",1,0);
         this->MyTypes_.push_back(ParamT("N2","Amine",2,0));
         this->MyTypes_.push_back(ParamT("N3","Amine",3,0));
      }
};


/**********************Oxygen Primitives**********************/

PRIM(Water,2,"O","Water",Oxygen,Hydrogen,Hydrogen)
PRIM(Hydroxyl,2,"O","Hydroxyl",Oxygen,Hydrogen)
PRIM(Ether,2,"O","Ether",Oxygen)
PRIM(ODB,1,"ODB","Double-Bonded Oxygen",Oxygen)

/*********************Sulfur primitives**************************/

PRIM(HydrogenSulfide,2,"S","Hydrogen Sulfide",Sulfur,Hydrogen,Hydrogen)
PRIM(Thiol,2,"S","Thiol",Sulfur,Hydrogen)
PRIM(Sulfide,2,"S","Sulfide",Sulfur)
PRIM(SDB,1,"SDB","Double-Bonded Sulfur",Sulfur)
#undef PRIM

}}//End namespaces
#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_PRIMITIVENODES_HH_ */
