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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_DERIVEDGROUPS_HH_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_DERIVEDGROUPS_HH_
/****This file not meant to be included in any other file other than
 *   FxnalGroup.h
 */
#include "FxnalGroup.h"
namespace psi{
namespace LibMolecule{

/****************** Derived Groups*****************************/

/** \brief Helper class to reduce copy/paste code
 *
 *  \param[in] MyType the type of the resulting group
 *  \param[in] i The number of groups comprising the resulting group
 *
 *  The rest of the parameters are the groups that comprise the new
 *  group, in the order they should be looked for
 *
 */
template<FxnGroup_t MyType,int i,int order,FxnGroup_t G1,FxnGroup_t G2,
         FxnGroup_t G3=NO_GROUP,FxnGroup_t G4=NO_GROUP,
         FxnGroup_t G5=NO_GROUP,FxnGroup_t G6=NO_GROUP>
class GenDerGrp:public DerivedFxnalGrp{
   public:
      typedef GenDerGrp<MyType,i,order,G1,G2,G3,G4,G5,G6> ThisType;
      GenDerGrp<MyType,i,order,G1,G2,G3,G4,G5,G6>():
            DerivedFxnalGrp(MyType,i,order){}
      FxnGroup_t GetPrimI(const int i)const{
         FxnGroup_t RValue=NO_GROUP;
         switch(i){
            case(0):{RValue=G1;break;}
            case(1):{RValue=G2;break;}
            case(2):{RValue=G3;break;}
            case(3):{RValue=G4;break;}
            case(4):{RValue=G5;break;}
            case(5):{RValue=G6;break;}
         }
         return RValue;
      }
};

class HydrogenCyanide:public GenDerGrp<HYDROGENCYANIDE,2,0,NITROTB,ALKYNYL1>{
   public:
      HydrogenCyanide(){}
      HydrogenCyanide(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 0;}
};

class Nitrile:public GenDerGrp<NITRILE,2,1,NITROTB,ALKYNYL2>{
   public:
      Nitrile(){}
      Nitrile(const Group_t& Groups){SetGroups(Groups);}
};

class Formaldehyde:public GenDerGrp<FORMALDEHYDE,2,0,OXYDB,ALKENYL1>{
   public:
      Formaldehyde(){}
      Formaldehyde(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 0;}
};

class Aldehyde:public GenDerGrp<ALDEHYDE,2,1,OXYDB,ALKENYL2>{
   public:
      Aldehyde(){}
      Aldehyde(const Group_t& Groups){SetGroups(Groups);}
};

class Carbonyl:public GenDerGrp<CARBONYL,2,2,OXYDB,ALKENYL3>{
   public:
      Carbonyl(){}
      Carbonyl(const Group_t& Groups){SetGroups(Groups);}
};

class Carboxyl:public GenDerGrp<CARBOXYL,2,1,HYDROXYL,CARBONYL>{
   public:
      Carboxyl(){}
      Carboxyl(const Group_t& Groups){SetGroups(Groups);}
};

class HydroPeroxy:public GenDerGrp<HYDROPEROXY,2,1,HYDROXYL,OXYGEN2>{
   public:
      HydroPeroxy(){}
      HydroPeroxy(const Group_t& Groups){SetGroups(Groups);}
};

class Peroxide:public GenDerGrp<PEROXIDE,2,0,HYDROXYL,HYDROXYL>{
   public:
      Peroxide(){}
      Peroxide(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 0;}
};

class Methoxy:public GenDerGrp<METHOXY,2,1,METHYL,OXYGEN2>{
   public:
      Methoxy(){}
      Methoxy(const Group_t& Groups){SetGroups(Groups);}
};

class Methanol:public GenDerGrp<METHANOL,2,0,HYDROXYL,METHYL>{
   public:
      Methanol(){}
      Methanol(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 0;}
};

class ccdb4:public GenDerGrp<CCDB4,2,4,ALKENYL3,ALKENYL3>{
   public:
      ccdb4(){}
      ccdb4(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 2;}
};

class ccdb3:public GenDerGrp<CCDB3,2,3,ALKENYL2,ALKENYL3>{
   public:
      ccdb3(){}
      ccdb3(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 2;}
};

class ccdb2:public GenDerGrp<CCDB2,2,2,ALKENYL2,ALKENYL2>{
   public:
      ccdb2(){}
      ccdb2(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 2;}
};

class Ethenyl2:public GenDerGrp<ETHENYL2,2,2,ALKENYL1,ALKENYL3>{
   public:
      Ethenyl2(){}
      Ethenyl2(const Group_t& Groups){SetGroups(Groups);}
};

class Ethenyl1:public GenDerGrp<ETHENYL1,2,1,ALKENYL1,ALKENYL2>{
   public:
      Ethenyl1(){}
      Ethenyl1(const Group_t& Groups){SetGroups(Groups);}
};

class Ethene:public GenDerGrp<ETHENE,2,0,ALKENYL1,ALKENYL1>{
   public:
      Ethene(){}
      Ethene(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 0;}
};

class cctb:public GenDerGrp<CCTB,2,2,ALKYNYL2,ALKYNYL2>{
   public:
      cctb(){}
      cctb(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 2;}
};

class Ethynyl:public GenDerGrp<ETHYNYL,2,1,ALKYNYL1,ALKYNYL2>{
   public:
      Ethynyl(){}
      Ethynyl(const Group_t& Groups){SetGroups(Groups);}
};

class Ethyne:public GenDerGrp<ETHYNE,2,0,ALKYNYL1,ALKYNYL1>{
   public:
      Ethyne(){}
      Ethyne(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 0;}
};

class Ketimine1:public GenDerGrp<KETIMINE1,2,2,NITRODB1,ALKENYL3>{
   public:
      Ketimine1(){}
      Ketimine1(const Group_t& Groups){SetGroups(Groups);}
};

class Ketimine2:public GenDerGrp<KETIMINE2,2,3,NITRODB2,ALKENYL3>{
   public:
      Ketimine2(){}
      Ketimine2(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 2;}
};

class Aldimine1:public GenDerGrp<ALDIMINE1,2,1,NITRODB1,ALKENYL2>{
   public:
      Aldimine1(){}
      Aldimine1(const Group_t& Groups){SetGroups(Groups);}
};

class Aldimine2:public GenDerGrp<ALDIMINE2,2,2,NITRODB2,ALKENYL2>{
   public:
      Aldimine2(){}
      Aldimine2(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 2;}
};

class Methanimine:public GenDerGrp<METHANIMINE,2,0,NITRODB1,ALKENYL1>{
   public:
      Methanimine(){}
      Methanimine(const Group_t& Groups){SetGroups(Groups);}
      int NAttachPoint()const{return 0;}
};

}}
#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_DERIVEDGROUPS_HH_ */
