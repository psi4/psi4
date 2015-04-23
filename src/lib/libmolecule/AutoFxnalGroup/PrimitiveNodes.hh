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
#define PRIM(Name,Order,Bonds,Abbv,Full,Centers...) \
   template<size_t...args>\
   class base_##Name:public RadialSearch<Bonds,Centers>{\
      private:\
         typedef RadialSearch<Bonds,Centers> Base_t;\
      protected:\
         bool IsPrim()const{return true;}\
      public:\
         base_##Name():Base_t(Abbv,Full){\
            if(sizeof...(args)!=0){\
               size_t Temp[sizeof...(args)]={args...};\
               PsiMap<size_t,size_t> Temp1;\
               for(size_t i=0;i<sizeof...(args);){\
                  size_t temp3=Temp[i++];\
                  Temp1[temp3]=Temp[i++];\
               }\
               MyTypes_[0]=ParamT(Abbv,Full,Temp1);\
            }\
         }\
};\
typedef base_##Name <> Name;\
typedef base_##Name <0,Order> Name##_t;

#define WILDCARD(Name,N,Abbv,Full)\
      class Name:public Node{\
         public:\
            Name():Node(0,Abbv,Full){\
               PsiMap<size_t,size_t> Temp;\
               Temp[0]=1;\
               this->MyTypes_[0]=ParamT(Abbv,Full,Temp,0);\
               for(int i=1;i<N;i++){\
                  ++Temp[0];\
                  this->MyTypes_.push_back(ParamT(Abbv,Full,Temp,0));\
               }\
            }\
      };


/*********Carbon Primitives********************/
PRIM(Methane,0,4,"C","Methane",Carbon,Hydrogen,Hydrogen,Hydrogen,Hydrogen)
PRIM(Methyl,1,4,"C","Methyl",Carbon,Hydrogen,Hydrogen,Hydrogen)
PRIM(Methene,2,4,"C","Methene",Carbon,Hydrogen,Hydrogen)
PRIM(Methyne,3,4,"C","Methyne",Carbon,Hydrogen)
PRIM(C4,4,4,"C","Carbon",Carbon)
PRIM(Alkenyl1,1,3,"CDB","Alkenyl",Carbon,Hydrogen,Hydrogen)
PRIM(Alkenyl2,2,3,"CDB","Alkenyl",Carbon,Hydrogen)
PRIM(Alkenyl3,3,3,"CDB","Alkenyl",Carbon)
PRIM(Alkynyl1,1,2,"CTB","Alkynyl",Carbon,Hydrogen)
PRIM(Alkynyl2,2,2,"CTB","Alkynyl",Carbon)

WILDCARD(Alkenyl,3,"CDB","Alkenyl")
WILDCARD(Alkynyl,2,"CTB","Alkynyl")

/**********************Nitrogen Primitives*******************/
PRIM(Ammonium,0,4,"N+","Ammonium",Nitrogen,Hydrogen,Hydrogen,Hydrogen,Hydrogen)
PRIM(Ammonium1,1,4,"N+","Ammonium",Nitrogen,Hydrogen,Hydrogen,Hydrogen)
PRIM(Ammonium2,2,4,"N+","Ammonium",Nitrogen,Hydrogen,Hydrogen)
PRIM(Ammonium3,3,4,"N+","Ammonium",Nitrogen,Hydrogen)
PRIM(Ammonium4,4,4,"N+","Ammonium",Nitrogen)
PRIM(Ammonia,0,3,"N","Ammonia",Nitrogen,Hydrogen,Hydrogen,Hydrogen)
PRIM(Amine1,1,3,"N","Amine",Nitrogen,Hydrogen,Hydrogen)
PRIM(Amine2,2,3,"N","Amine",Nitrogen,Hydrogen)
PRIM(Amine3,3,3,"N","Amine",Nitrogen)
PRIM(Azo1,1,2,"NDB","Azo",Nitrogen,Hydrogen)
PRIM(Azo2,2,2,"NDB","Azo",Nitrogen)
PRIM(NTB,1,1,"NTB","Triple-Bonded Nitrogen",Nitrogen)

WILDCARD(Azo,2,"NDB","Azo")
WILDCARD(Amine,3,"N","Amine")

/**********************Oxygen Primitives**********************/

PRIM(Water,0,2,"O","Water",Oxygen,Hydrogen,Hydrogen)
PRIM(Hydroxyl,1,2,"O","Hydroxyl",Oxygen,Hydrogen)
PRIM(Ether,2,2,"O","Ether",Oxygen)
PRIM(ODB,1,1,"ODB","Double-Bonded Oxygen",Oxygen)

/*********************Sulfur primitives**************************/

PRIM(HydrogenSulfide,0,2,"S","Hydrogen Sulfide",Sulfur,Hydrogen,Hydrogen)
PRIM(Thiol,1,2,"S","Thiol",Sulfur,Hydrogen)
PRIM(Sulfide,2,2,"S","Sulfide",Sulfur)
PRIM(SDB,1,1,"SDB","Double-Bonded Sulfur",Sulfur)


/*********************Conjugated-Pi Primitives*************************

class ConPiPrims:public Node{
         public:
            ConPiPrims():Node(0,"Pi*","Conjugated Pi Atom"){
               PsiMap<size_t,size_t> Temp;
               this->MyTypes_[0]=ParamT("CDB","Alkenyl",2,0);
               this->MyTypes_.push_back(ParamT("CDB","Alkenyl",3,0));
               this->MyTypes_.push_back(ParamT("NDB","Azo",2,0));
               this->MyTypes_.push_back(ParamT("N","Amine",2,0));
               this->MyTypes_.push_back(ParamT("O","Ether",2,0));
            }
      };*/

#undef PRIM
#undef WILDCARD
}}//End namespaces
#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_PRIMITIVENODES_HH_ */
