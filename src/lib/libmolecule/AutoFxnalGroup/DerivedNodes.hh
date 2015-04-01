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
typedef FxnGrpType F_t;
/*************C-C 2X Bonds*************/
class Ethene: public LinearSearch<Alkenyl1,Alkenyl1>{
   private:
      typedef LinearSearch<Alkenyl1,Alkenyl1> Base_t;
   public:
      Ethene():Base_t(F_t::DBCC0,F_t::CDBCC0,F_t::HDBCC0,F_t::HDBCC0,
                      F_t::CDBCC0,F_t::HDBCC0,F_t::HDBCC0){}
};

class DBCC1: public LinearSearch<Alkenyl1,Alkenyl2>{
   private:
      typedef LinearSearch<Alkenyl1,Alkenyl2> Base_t;
   public:
      DBCC1():Base_t(F_t::DBCC1,F_t::C2DBCC1,F_t::H2DBCC1,F_t::H2DBCC1,
                     F_t::C1DBCC1,F_t::H1DBCC1){}
};

class DBCC2: public LinearSearch<Alkenyl2,Alkenyl2>{
   private:
      typedef LinearSearch<Alkenyl2,Alkenyl2> Base_t;
   public:
      DBCC2():Base_t(F_t::DBCC2,F_t::CDBCC2,F_t::HDBCC2,
                     F_t::CDBCC2,F_t::HDBCC2){}
};

class DBCC2G: public LinearSearch<Alkenyl1,Alkenyl3>{
   private:
      typedef LinearSearch<Alkenyl1,Alkenyl3> Base_t;
   public:
      DBCC2G():Base_t(F_t::DBCC2G,F_t::C2DBCC2G,F_t::HDBCC2G,
                      F_t::HDBCC2G,F_t::C1DBCC2G){}
};

class DBCC3: public LinearSearch<Alkenyl2,Alkenyl3>{
   private:
      typedef LinearSearch<Alkenyl2,Alkenyl3> Base_t;
   public:
      DBCC3():Base_t(F_t::DBCC3,F_t::C2DBCC3,F_t::HDBCC3,
                      F_t::C1DBCC3){}
};

class DBCC4: public LinearSearch<Alkenyl3,Alkenyl3>{
   private:
      typedef LinearSearch<Alkenyl3,Alkenyl3> Base_t;
   public:
      DBCC4():Base_t(F_t::DBCC4,F_t::CDBCC4,F_t::CDBCC4){}
};

/*************C-C 3X Bonds*****************/
class Acetylene: public LinearSearch<Alkynyl1,Alkynyl1>{
   private:
      typedef LinearSearch<Alkynyl1,Alkynyl1> Base_t;
   public:
      Acetylene():Base_t(F_t::TBCC0,F_t::CTBCC0,F_t::HTBCC0,
                                  F_t::CTBCC0,F_t::HTBCC0){}
};

class Ethynyl: public LinearSearch<Alkynyl1,Alkynyl2>{
   private:
      typedef LinearSearch<Alkynyl1,Alkynyl2> Base_t;
   public:
      Ethynyl():Base_t(F_t::TBCC1,F_t::C2TBCC1,F_t::HTBCC1,
                                  F_t::C1TBCC1){}
};

class CCTB2: public LinearSearch<Alkynyl2,Alkynyl2>{
   private:
      typedef LinearSearch<Alkynyl2,Alkynyl2> Base_t;
   public:
      CCTB2():Base_t(F_t::TBCC2,F_t::CTBCC2,F_t::CTBCC2){}
};

/*************C-N 2X Bonds*****************/

class Formaldimine: public LinearSearch<Alkenyl1,Azo1>{
   private:
      typedef LinearSearch<Alkenyl1,Azo1> Base_t;
   public:
      Formaldimine():Base_t(F_t::DBCN0,F_t::CDBCN0,F_t::H2DBCN0,F_t::H2DBCN0,
                            F_t::NDBCN0,F_t::H1DBCN0){}
};

class DBCN1: public LinearSearch<Alkenyl1,Azo2>{
   private:
      typedef LinearSearch<Alkenyl1,Azo2> Base_t;
   public:
      DBCN1():Base_t(F_t::DBCN1N,F_t::CDBCN1N,F_t::HDBCN1N,F_t::HDBCN1N,
                            F_t::NDBCN1N){}
};

class Aldimine1: public LinearSearch<Alkenyl2,Azo1>{
   private:
      typedef LinearSearch<Alkenyl2,Azo1> Base_t;
   public:
      Aldimine1():Base_t(F_t::DBCN1C,F_t::CDBCN1C,F_t::H2DBCN1C,
                            F_t::NDBCN1C,F_t::H1DBCN1C){}
};

class Aldimine2: public LinearSearch<Alkenyl2,Azo2>{
   private:
      typedef LinearSearch<Alkenyl2,Azo2> Base_t;
   public:
      Aldimine2():Base_t(F_t::DBCN2,F_t::CDBCN2,F_t::HDBCN2,F_t::NDBCN2){}
};

class Ketimine1: public LinearSearch<Alkenyl3,Azo1>{
   private:
      typedef LinearSearch<Alkenyl3,Azo1> Base_t;
   public:
      Ketimine1():Base_t(F_t::DBCN2G,F_t::CDBCN2G,F_t::NDBCN2G,F_t::HDBCN2G){}
};

class Ketimine2:public LinearSearch<Alkenyl3,Azo2>{
   private:
      typedef LinearSearch<Alkenyl3,Azo2> Base_t;
   public:
      Ketimine2():Base_t(F_t::DBCN3,F_t::CDBCN3,F_t::NDBCN3){}

};

/*****************C-N 3X Bonds*************/
class HydrogenCyanide:public LinearSearch<Alkynyl1,NTB>{
   private:
      typedef LinearSearch<Alkynyl1,NTB> Base_t;
   public:
      HydrogenCyanide():Base_t(F_t::TBCN0,F_t::CTBCN0,F_t::HTBCN0,F_t::NTBCN0){}
};

class Nitrile:public LinearSearch<Alkynyl2,NTB>{
   private:
      typedef LinearSearch<Alkynyl2,NTB> Base_t;
   public:
      Nitrile():Base_t(F_t::TBCN1,F_t::CTBCN1,F_t::NTBCN1){}
};

/**************C-O 2X Bonds*****************/
class Formaldehyde: public LinearSearch<Alkenyl1,ODB>{
   private:
      typedef LinearSearch<Alkenyl1,ODB> Base_t;
   public:
      Formaldehyde():Base_t(F_t::DBCO0,F_t::CDBCO0,F_t::HDBCO0,
                            F_t::HDBCO0,F_t::ODBCO0){}
};

class Aldehyde: public LinearSearch<Alkenyl2,ODB>{
   private:
      typedef LinearSearch<Alkenyl2,ODB> Base_t;
   public:
      Aldehyde():Base_t(F_t::DBCO1,F_t::CDBCO1,F_t::HDBCO1,F_t::ODBCO1){}
};

class Carbonyl: public LinearSearch<Alkenyl3,ODB>{
   private:
      typedef LinearSearch<Alkenyl3,ODB> Base_t;
   public:
      Carbonyl():Base_t(F_t::DBCO2,F_t::CDBCO2,F_t::ODBCO2){}
};


/**************Carbonyl Groups********************/

class Carboxylate: public LinearSearch<Carbonyl,ODB>{
   private:
      typedef LinearSearch<Carbonyl,ODB> Base_t;
   public:
      Carboxylate():Base_t(F_t::CO2M,F_t::CCO2M,F_t::OCO2M,F_t::OCO2M){}
};

class Carboxyl: public LinearSearch<Carbonyl,Hydroxyl>{
   private:
      typedef LinearSearch<Carbonyl,Hydroxyl> Base_t;
   public:
      Carboxyl():Base_t(F_t::CO2H,F_t::CCO2H,F_t::O2CO2H,F_t::O1CO2H,F_t::HCO2H){}
};

class Ester:public LinearSearch<Carbonyl,Ether>{
   private:
      typedef LinearSearch<Carbonyl,Ether> Base_t;
   public:
      Ester():Base_t(F_t::CO2,F_t::CCO2,F_t::O2CO2,F_t::O1CO2){}

};

class Carbonate: public LinearSearch<Ether,Carbonyl,Ether>{
   private:
      typedef LinearSearch<Ether,Carbonyl,Ether> Base_t;
   public:
      Carbonate():Base_t(F_t::CO3,F_t::O1CO3,F_t::CCO3,F_t::O2CO3,F_t::O1CO3){}
};

class Amide1C: public LinearSearch<Carbonyl,Amine1>{
   private:
      typedef LinearSearch<Carbonyl,Amine1> Base_t;
   public:
      Amide1C():Base_t(F_t::CONH2,F_t::CCONH2,F_t::OCONH2,F_t::NCONH2,
                       F_t::HCONH2,F_t::HCONH2){}
};

class Amide1N: public LinearSearch<Aldehyde,Amine2>{
   private:
      typedef LinearSearch<Aldehyde,Amine2> Base_t;
   public:
      Amide1N():Base_t(F_t::CHONH,F_t::CCHONH,F_t::H2CHONH,F_t::OCHONH,
                       F_t::NCHONH,F_t::H1CHONH){}
};

class Amide2: public LinearSearch<Carbonyl,Amine2>{
   private:
      typedef LinearSearch<Carbonyl,Amine2> Base_t;
   public:
      Amide2():Base_t(F_t::CONH,F_t::CCONH,F_t::OCONH,F_t::NCONH,
                      F_t::HCONH){}
};

class Amide2G: public LinearSearch<Aldehyde,Amine3>{
   private:
      typedef LinearSearch<Aldehyde,Amine3> Base_t;
   public:
      Amide2G():Base_t(F_t::CHON,F_t::CCHON,F_t::HCHON,F_t::OCHON,
                      F_t::NCHON){}
};

class Amide3: public LinearSearch<Carbonyl,Amine3>{
   private:
      typedef LinearSearch<Carbonyl,Amine3> Base_t;
   public:
      Amide3():Base_t(F_t::CON,F_t::CCON,F_t::OCON,F_t::NCON){}
};

class Imide1C: public LinearSearch<Carbonyl,Amide1N>{
   private:
      typedef LinearSearch<Carbonyl,Amide1N> Base_t;
   public:
      Imide1C():Base_t(F_t::CONHCOH,F_t::C1CONHCOH,F_t::O1CONHCOH,
                       F_t::C2CONHCOH,F_t::H2CONHCOH,F_t::O2CONHCOH,
                       F_t::NCONHCOH,F_t::H1CONHCOH){}
};

class Imide1N: public LinearSearch<Aldehyde,Amide2G>{
   private:
      typedef LinearSearch<Aldehyde,Amide2G> Base_t;
   public:
      Imide1N():Base_t(F_t::CHONCOH,
                       F_t::CCHONCOH,F_t::HCHONCOH,F_t::OCHONCOH,
                       F_t::CCHONCOH,F_t::HCHONCOH,F_t::OCHONCOH,
                       F_t::NCHONCOH){}
};

class Imide2: public LinearSearch<Carbonyl,Amide2G>{
   private:
      typedef LinearSearch<Carbonyl,Amide2G> Base_t;
   public:
      Imide2():Base_t(F_t::CONCOH,
            F_t::C1CONCOH,F_t::O1CONCOH,
            F_t::C2CONCOH,F_t::O2CONCOH,F_t::HCONCOH,
            F_t::NCONCOH){}
};

class Imide2C: public LinearSearch<Carbonyl,Amide2>{
   private:
      typedef LinearSearch<Carbonyl,Amide2> Base_t;
   public:
      Imide2C():Base_t(F_t::CONHCO,
            F_t::CCONHCO,F_t::OCONHCO,
            F_t::CCONHCO,F_t::OCONHCO,
            F_t::NCONHCO,F_t::HCONHCO){}
};

class Imide3: public LinearSearch<Carbonyl,Amide3>{
   private:
      typedef LinearSearch<Carbonyl,Amide3> Base_t;
   public:
      Imide3():Base_t(F_t::CONCO,
            F_t::CCONCO,F_t::OCONCO,
            F_t::CCONCO,F_t::OCONCO,F_t::NCONCO
            ){}
};

/*************** Ether Groups ***************/
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

/**************Miscellaneous Nitrogen Groups**********/
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

/**********Miscalaneous Carbon Groups***********/
class Isopropyl: public LinearSearch<Methyl,Methyne,Methyl>{
   private:
      typedef LinearSearch<Methyl,Methyne,Methyl> Base_t;
   public:
      Isopropyl():Base_t(F_t::IPR,F_t::C2IPR,F_t::H2IPR,F_t::H2IPR,F_t::H2IPR,
            F_t::C1IPR,F_t::H1IPR,
            F_t::C2IPR,F_t::H2IPR,F_t::H2IPR,F_t::H2IPR
            ){}
};


}}



#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_DERIVEDNODES_HH_ */
