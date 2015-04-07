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
namespace psi{
namespace LibMolecule{
#define DERIV(Name,Abbv,Full,Centers...) \
   class Name:public LinearSearch<Centers>{\
      private:\
         typedef LinearSearch<Centers> Base_t;\
      public:\
         Name():Base_t(Abbv,Full){}\
};

/*************C-C Multiple Bonds*************/
DERIV(DBCC,"DBCC","C=C 2X Bond",Alkenyl,Alkenyl)
DERIV(TBCC,"TBCC","C=C 3X Bond",Alkynyl,Alkynyl)


/*************C-N 2X Bonds*****************/
DERIV(Formaldimine,"DBCN","Formaldimine",Alkenyl1,Azo1)
//I can't for the life of me find a name for this
DERIV(DBCN1,"DBCN","Stupid Aldimine",Alkenyl1,Azo2)
DERIV(Aldimine1,"DBCN","Aldimine",Alkenyl2,Azo1)
DERIV(Aldimine2,"DBCN","Aldimine",Alkenyl2,Azo2)
DERIV(Ketimine1,"DBCN","Ketimine",Alkenyl3,Azo1)
DERIV(Ketimine2,"DBCN","Ketimine",Alkenyl3,Azo2)

/*****************C-N 3X Bonds*************/
DERIV(HydrogenCyanide,"TBCN","Hydrogen Cyanide",Alkynyl1,NTB)
DERIV(Nitrile,"TBCN","Nitrile",Alkynyl2,NTB)

/**************C-O 2X Bonds*****************/
DERIV(Formaldehyde,"DBC","Formaldehyde",Alkenyl1,ODB)
DERIV(Aldehyde,"DBCO","Aldehyde",Alkenyl2,ODB)
DERIV(Carbonyl,"DBCO","Carbonyl",Alkenyl3,ODB)


/**************Carbonyl Groups********************/
DERIV(Carboxylate,"CO2M","Carboxylate",Carbonyl,ODB)
DERIV(Carboxyl,"CO2H","Carboxyl",Carbonyl,Hydroxyl)
DERIV(Ester,"OCO","Ester",Carbonyl,Ether)
DERIV(Carbonate,"CO3","Carbonate",Ether,Carbonyl,Ether)
DERIV(Amide,"CON","Amide",Carbonyl,Amine)

class Imide: public LinearSearch<Carbonyl,Amide>{
   private:
      typedef LinearSearch<Carbonyl,Amide> Base_t;
      ParamT RHSType(const unsigned i){
         return (i==1?ParamT("N*","Amine",0,0):Base_t::RHSType(i));
      }
   public:
      Imide():Base_t("CONCO","Imide"){}
};

/*************** Ether Groups ***************
class Methanol: public LinearSearch<Hydroxyl,Methyl>{
   private:
      typedef LinearSearch<Hydroxyl,Methyl> Base_t;
   public:
      Methanol():Base_t(F_t::OCH4,F_t::OOCH4,F_t::H1OCH4,
                        F_t::COCH4,F_t::H2OCH4,F_t::H2OCH4,
                        F_t::H2OCH4,F_t::H2OCH4){}
};

class Methoxy: public LinearSearch<Ether,Methyl>{
   private:
      typedef LinearSearch<Ether,Methyl> Base_t;
   public:
      Methoxy():Base_t(F_t::OCH3,F_t::OOCH3,F_t::COCH3,F_t::HOCH3,
                                 F_t::HOCH3,F_t::HOCH3,F_t::HOCH3){}
};

class HydrogenPeroxide: public LinearSearch<Hydroxyl,Hydroxyl>{
   private:
      typedef LinearSearch<Hydroxyl,Hydroxyl> Base_t;
   public:
      HydrogenPeroxide():Base_t(F_t::O2H2,F_t::HO2H2,F_t::OO2H2,
                                          F_t::OO2H2,F_t::HO2H2){}
};

class Hydroperoxide: public LinearSearch<Ether,Hydroxyl>{
   private:
      typedef LinearSearch<Ether,Hydroxyl> Base_t;
   public:
      Hydroperoxide():Base_t(F_t::O2H,F_t::O1O2H,F_t::O2O2H,F_t::HO2H){}
};

class Peroxide: public LinearSearch<Ether,Ether>{
   private:
      typedef LinearSearch<Ether,Ether> Base_t;
   public:
      Peroxide():Base_t(F_t::OO,F_t::OOO,F_t::OOO){}
};

class Hemiacetal: public LinearSearch<Ether,Methyne,Hydroxyl>{
   private:
      typedef LinearSearch<Ether,Methyne,Hydroxyl> Base_t;
   public:
      Hemiacetal():Base_t(F_t::OCHOH,F_t::O1OCHOH,F_t::COCHOH,
                  F_t::H2OCHOH,F_t::O2OCHOH,F_t::H1OCHOH){}
};

class Hemiketal: public LinearSearch<Ether,C4,Hydroxyl>{
   private:
      typedef LinearSearch<Ether,C4,Hydroxyl> Base_t;
   public:
      Hemiketal():Base_t(F_t::OCOH,F_t::O1OCOH,F_t::COCOH,
                         F_t::O2OCOH,F_t::HOCOH){}
};

class Acetal: public LinearSearch<Ether,Methyne,Ether>{
   private:
      typedef LinearSearch<Ether,Methyne,Ether> Base_t;
   public:
      Acetal():Base_t(F_t::OCHO,F_t::OOCHO,F_t::COCHO,
                      F_t::HOCHO,F_t::OOCHO){}
};

class Ketal: public LinearSearch<Ether,C4,Ether>{
   private:
      typedef LinearSearch<Ether,C4,Ether> Base_t;
   public:
      Ketal():Base_t(F_t::OCO,F_t::OOCO,F_t::COCO,F_t::OOCO){}
};

class Orthoester: public RadialSearch<4,C4,Ether,Ether,Ether>{
   private:
      typedef RadialSearch<4,C4,Ether,Ether,Ether> Base_t;
   public:
      Orthoester():Base_t(F_t::COOO,F_t::CCOOO,F_t::OCOOO,F_t::OCOOO,
                          F_t::OCOOO){}
};

/**************Miscellaneous Nitrogen Groups**********
class Azide: public LinearSearch<Azo1,Azo2,Azo2>{
   private:
      typedef LinearSearch<Azo1,Azo2,Azo2> Base_t;
   public:
      Azide():Base_t(F_t::NNN,F_t::N3NNN,F_t::N2NNN,F_t::N1NNN){}
};

class Diazene:public LinearSearch<Azo1,Azo1>{
   private:
      typedef LinearSearch<Azo1,Azo1> Base_t;
   public:
      Diazene():Base_t(F_t::DBNN0,F_t::NDBNN0,F_t::HDBNN0,
                       F_t::NDBNN0,F_t::HDBNN0){}
};

class DBNN1:public LinearSearch<Azo1,Azo2>{
   private:
      typedef LinearSearch<Azo1,Azo2> Base_t;
   public:
      DBNN1():Base_t(F_t::DBNN1,F_t::N2DBNN1,F_t::HDBNN1,
                       F_t::N2DBNN1){}
};

class DBNN2:public LinearSearch<Azo2,Azo2>{
   private:
      typedef LinearSearch<Azo2,Azo2> Base_t;
   public:
      DBNN2():Base_t(F_t::DBNN2,F_t::NDBNN2,F_t::NDBNN2){}
};

class Cyanate:public LinearSearch<Nitrile,Ether>{
   private:
      typedef LinearSearch<Nitrile,Ether> Base_t;
   public:
      Cyanate():Base_t(F_t::OCN,F_t::COCN,F_t::NOCN,F_t::OOCN){}
};

class Isocyanate:public LinearSearch<Azo2,Alkynyl2,ODB>{
   private:
      typedef LinearSearch<Azo2,Alkynyl2,ODB> Base_t;
   public:
      Isocyanate():Base_t(F_t::NCO,F_t::NNCO,F_t::CNCO,F_t::ONCO){}
};

class Nitro:public LinearSearch<ODB,Amine3,ODB>{
   private:
      typedef LinearSearch<ODB,Amine3,ODB> Base_t;
   public:
      Nitro():Base_t(F_t::NOO,F_t::ONOO,F_t::NNOO,F_t::ONOO){}
};

class Nitrate:public LinearSearch<Nitro,Ether>{
   private:
      typedef LinearSearch<Nitro,Ether> Base_t;
   public:
      Nitrate():Base_t(F_t::NO3,F_t::O2NO3,F_t::NNO3,F_t::O2NO3,
                       F_t::O1NO3){}
};

class Nitroso:public LinearSearch<Amine2,ODB>{
   private:
      typedef LinearSearch<Amine2,ODB> Base_t;
   public:
      Nitroso():Base_t(F_t::NO,F_t::NNO,F_t::ONO){}
};

class Nitrite:public LinearSearch<Nitroso,Ether>{
   private:
      typedef LinearSearch<Nitroso,Ether> Base_t;
   public:
      Nitrite():Base_t(F_t::NO2,F_t::NNO2,F_t::O2NO2,F_t::O1NO2){}
};

/**********Miscalaneous Carbon Groups***********
class Isopropyl: public LinearSearch<Methyl,Methyne,Methyl>{
   private:
      typedef LinearSearch<Methyl,Methyne,Methyl> Base_t;
   public:
      Isopropyl():Base_t(F_t::IPR,F_t::C2IPR,F_t::H2IPR,F_t::H2IPR,F_t::H2IPR,
            F_t::C1IPR,F_t::H1IPR,
            F_t::C2IPR,F_t::H2IPR,F_t::H2IPR,F_t::H2IPR
            ){}
};
*/

#undef DERIV
}}



#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_DERIVEDNODES_HH_ */
