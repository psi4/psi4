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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_DERIVEDNODES_HH_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_DERIVEDNODES_HH_

#include "PrimitiveNodes.hh"
#include "LinearSearch.h"
namespace psi {
namespace LibMolecule {
#define PRIM(Name,Abbv,Full,Centers...) \
   template<size_t...args>\
   class base_##Name:public LinearSearch<Centers>{\
      private:\
         typedef LinearSearch<Centers> Base_t;\
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
typedef base_##Name <> Name;
#define DERIV(Name,Abbv,Full,Centers...) \
   template<size_t...args>\
   class base_##Name:public LinearSearch<Centers>{\
      private:\
         typedef LinearSearch<Centers> Base_t;\
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
   typedef base_##Name <> Name;

/*************C-C Multiple Bonds*************/
DERIV(Ethene, "DBCC", "Ethene", Alkenyl1_t, Alkenyl1_t)
DERIV(DBCC, "DBCC", "C=C 2X Bond", Alkenyl, Alkenyl)
DERIV(Acetylene, "TBCC", "Acetylene", Alkynyl1_t, Alkynyl1_t)
DERIV(TBCC, "TBCC", "C=C 3X Bond", Alkynyl, Alkynyl)

/*************C-N 2X Bonds*****************/
DERIV(Formaldimine, "DBCN", "Formaldimine", Alkenyl1_t, Azo1_t)
//I can't for the life of me find a name for this
DERIV(DBCN1, "DBCN", "Stupid Aldimine", Alkenyl1_t, Azo2_t)
DERIV(Aldimine1, "DBCN", "Aldimine", Alkenyl2_t, Azo1_t)
DERIV(Aldimine2, "DBCN", "Aldimine", Alkenyl2_t, Azo2_t)
DERIV(Ketimine1, "DBCN", "Ketimine", Alkenyl3_t, Azo1_t)
typedef base_Ketimine1<0,2> Ketimine1_t;
DERIV(Ketimine2, "DBCN", "Ketimine", Alkenyl3_t, Azo2_t)
typedef base_Ketimine2<0,2,1,1> Ketimine2_t;

/*****************C-N 3X Bonds*************/
DERIV(HydrogenCyanide, "TBCN", "Hydrogen Cyanide", Alkynyl1_t, NTB_t)
DERIV(Nitrile, "TBCN", "Nitrile", Alkynyl2_t, NTB_t)
typedef base_Nitrile<0,1> Nitrile_t;

/**************C-O 2X Bonds*****************/
template<size_t...args>
class base_Carboxylate:public RadialSearch<3,Alkenyl3_t,ODB_t,ODB_t>{
   private:
      typedef RadialSearch<3,Alkenyl3_t,ODB_t,ODB_t> Base_t;
   protected:
      std::vector<size_t> DeterminePrior(
            std::deque<boost::shared_ptr<Node> >&,
            bool&)const{
         std::vector<size_t> Priors(3,0);
         Priors[1]=1;
         Priors[2]=1;
         return Priors;
      }
      bool IsPrim()const{return true;}
   public:
      base_Carboxylate():Base_t("CO2-","Carboxylate"){
         if(sizeof...(args)!=0){
            size_t Temp[sizeof...(args)]={args...};
            PsiMap<size_t,size_t> Temp1;
            for(size_t i=0;i<sizeof...(args);){
               size_t temp3=Temp[i++];
               Temp1[temp3]=Temp[i++];
            }
            MyTypes_[0]=ParamT("CO2-","Carboxylate",Temp1);
         }
      }
   };
typedef base_Carboxylate <> Carboxylate;
typedef base_Carboxylate<0,1> Carboxylate_t;
DERIV(Formaldehyde, "DBCO", "Formaldehyde", Alkenyl1_t, ODB_t)
DERIV(Aldehyde, "DBCO", "Aldehyde", Alkenyl2_t, ODB_t)
typedef base_Aldehyde<0,1> Aldehyde_t;
DERIV(Carbonyl, "DBCO", "Carbonyl", Alkenyl3_t, ODB_t)
typedef base_Carbonyl<0,2> Carbonyl_t;


/**************Carbonyl Groups********************/

DERIV(Carboxyl, "CO2H", "Carboxyl", Carbonyl_t, Hydroxyl_t)
typedef base_Carboxyl<0,1> Carboxyl_t;
DERIV(Ester, "OCO", "Ester", Carbonyl_t, Ether_t)
typedef base_Ester<0,1,1,1> Ester_t;
DERIV(Carbonate, "CO3", "Carbonate", Ether_t, Carbonyl_t, Ether_t)
typedef base_Carbonate<0,1,2,1> Carbonate_t;
DERIV(Amide, "CON", "Amide", Carbonyl_t, Amine)
DERIV(Imide, "CONCO", "Imide", Carbonyl_t, Amine, Carbonyl_t)

/*************** Ether Groups ***************/
DERIV(Methanol, "HOCH3", "Methanol", Hydroxyl_t, Methyl_t)
DERIV(Methoxy, "OCH3", "Methoxy", Ether_t, Methyl_t)
typedef base_Methoxy<0,1> Methoxy_t;
DERIV(HydrogenPeroxide, "HOOH", "Hydrogen Peroxide", Hydroxyl_t, Hydroxyl_t)
DERIV(Hydroperoxide, "OOH", "Hydroperoxide", Ether_t, Hydroxyl_t)
typedef base_Hydroperoxide<0,1> Hydroperoxide_t;
DERIV(Peroxide, "OO", "Peroxide", Ether_t, Ether_t)
typedef base_Peroxide<0,1,1,1> Peroxide_t;
DERIV(Hemiacetal, "OCHOH", "Hemiacetal", Ether_t, Methyne_t, Hydroxyl_t)
DERIV(Hemiketal, "OCOH", "Hemiketal", Ether_t, C4_t, Hydroxyl_t)
DERIV(Acetal, "OCHO", "Acetal", Ether_t, Methyne_t, Ether_t)
DERIV(Ketal, "OCO", "Ketal", Ether_t, C4_t, Ether_t)
template<size_t...args>
class base_OrthoEster:
      public RadialSearch<4,C4_t, Ether_t, Ether_t, Ether_t>{
      private:
         typedef RadialSearch<4,C4_t, Ether_t, Ether_t, Ether_t> Base_t;
      public:
         base_OrthoEster():Base_t("COOO", "Orthoester"){
            if(sizeof...(args)!=0){
               std::vector<size_t> Temp{args...};
               PsiMap<size_t,size_t> Temp1;
               for(size_t i=0;i<sizeof...(args);){
                  size_t temp3=Temp[i++];
                  Temp1[temp3]=Temp[i++];
               }
               MyTypes_[0]=ParamT("COOO", "Orthoester",Temp1);
            }
         }
};
typedef base_OrthoEster <> OrthoEster;

/**************Miscellaneous Nitrogen Groups**********/
DERIV(Diazene, "DBNN", "Diazene", Azo1_t, Azo1_t)
DERIV(DBNN, "DBNN", "N=N 2x Bond", Azo, Azo)
DERIV(Azide, "NNN", "Azide", NTB_t, base_DBNN<0,1,1,1>)
DERIV(Cyanate, "OCN", "Cyanate", Nitrile_t, Ether_t)
DERIV(Isocyanate, "NCO", "Isocyanate", Azo2_t, Alkynyl2_t, ODB_t)
DERIV(Nitro, "NOO", "Nitro", ODB_t, Amine3_t, ODB_t)
typedef base_Nitro<0,1> Nitro_t;
template<size_t...args>
class base_Nitrate:public LinearSearch<Nitro_t,Ether_t>{
   private:
      typedef LinearSearch<Nitro_t,Ether_t> Base_t;
   protected:
      std::vector<size_t> DeterminePrior(
            std::deque<boost::shared_ptr<Node> >& FN,
            bool& Symm)const{
         std::vector<size_t> Priors=Base_t::DeterminePrior(FN,Symm);
         Symm=true;
         return Priors;
      }
   public:
      base_Nitrate():Base_t("NO3","Nitrate"){
         if(sizeof...(args)!=0){
            size_t Temp[sizeof...(args)]={args...};
            PsiMap<size_t,size_t> Temp1;
            for(size_t i=0;i<sizeof...(args);){
               size_t temp3=Temp[i++];
               Temp1[temp3]=Temp[i++];
            }
            MyTypes_[0]=ParamT("NO3","Nitrate",Temp1);
         }
      }
   };
typedef base_Nitrate <> Nitrate;
DERIV(Nitroso, "NO", "Nitroso", Azo2_t, ODB_t)
typedef base_Nitroso<0,1> Nitroso_t;
DERIV(Nitrite, "NO2", "Nitrite", Nitroso_t, Ether_t)

/**********Miscalaneous Carbon Groups***********/
PRIM(Isopropyl,"IPR","Isopropyl",Methyl_t, Methyne_t, Methyl_t)
typedef base_Isopropyl<1,1> Isopropyl_t;
DERIV(Butyl,"BT","2-Butyl",Methyl_t,Methyne_t,Methene_t,Methyl_t)
typedef base_Butyl<1,1> Butyl_t;
DERIV(Hydroxyethyl,"HEt","2-Hydroxyethyl",Hydroxyl_t,Methyne_t,Methyl_t)
typedef base_Hydroxyethyl<1,1> Hydroxyethyl_t;
PRIM(ArgThing,"ARGR","Argine End",Amine1_t,Alkenyl3_t,Amine1_t)
typedef base_ArgThing<1,1> ArgThing_t;

#undef DERIV
#undef PRIM
}
}

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_DERIVEDNODES_HH_ */
