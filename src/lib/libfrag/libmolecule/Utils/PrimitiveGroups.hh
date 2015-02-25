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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_PRIMITIVEGROUPS_HH_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_PRIMITIVEGROUPS_HH_
/****This file not meant to be included in any other file other than
 *   FxnalGroup.h
 */
#include "FxnalGroup.h"
#include "FindDerGroup.h"
namespace psi{
namespace LibMolecule{
/*********** Primitive Groups ***************************************/

class Methane:public DerivedFxnalGrp{
   public:
      //C first
      Methane():DerivedFxnalGrp(METHANE,5,0){}
      int NAttachPoint()const{return 0;}
      typedef FindPrimGroup<4,Methane,Carbon,Hydrogen,Hydrogen,
                                             Hydrogen,Hydrogen> FindMe;
};

class Methyl:public DerivedFxnalGrp{
   public:
      //C first
      Methyl():DerivedFxnalGrp(METHYL,4,1){}
      typedef FindPrimGroup<4,Methyl,Carbon,Hydrogen,Hydrogen,
                                             Hydrogen> FindMe;
};

class Methene:public DerivedFxnalGrp{
   public:
      //C first
      Methene():DerivedFxnalGrp(METHENE,3,2){}
      typedef FindPrimGroup<4,Methene,Carbon,Hydrogen,Hydrogen> FindMe;
};

class Methyne:public DerivedFxnalGrp{
   public:
      //C first
      Methyne():DerivedFxnalGrp(METHYNE,3,3){}
      typedef FindPrimGroup<4,Methane,Carbon,Hydrogen> FindMe;
};

class Carbon4:public DerivedFxnalGrp{
   public:
      //Only has a C
      Carbon4():DerivedFxnalGrp(CARBON4,1,4){}
      typedef FindPrimGroup<4,Carbon4,Carbon> FindMe;
};

class Alkenyl1:public DerivedFxnalGrp{
   public:
      //C first
      Alkenyl1():DerivedFxnalGrp(ALKENYL1,3,1){}
      typedef FindPrimGroup<3,Alkenyl1,Carbon,Hydrogen,Hydrogen> FindMe;
};

class Alkenyl2:public DerivedFxnalGrp{
   public:
      //C first
      Alkenyl2():DerivedFxnalGrp(ALKENYL2,2,2){}
      typedef FindPrimGroup<3,Alkenyl2,Carbon,Hydrogen> FindMe;
};

class Alkenyl3:public DerivedFxnalGrp{
   public:
      //C first
      Alkenyl3():DerivedFxnalGrp(ALKENYL3,1,3){}
      typedef FindPrimGroup<3,Alkenyl3,Carbon> FindMe;
};

class Alkynyl1:public DerivedFxnalGrp{
   public:
      //C first
      Alkynyl1():DerivedFxnalGrp(ALKYNYL1,2,1){}
      typedef FindPrimGroup<2,Alkynyl1,Carbon,Hydrogen> FindMe;
};

class Alkynyl2:public DerivedFxnalGrp{
   public:
      //C first
      Alkynyl2():DerivedFxnalGrp(ALKYNYL2,1,2){}
      typedef FindPrimGroup<2,Alkynyl2,Carbon> FindMe;
};

class Water:public DerivedFxnalGrp{
   public:
      int NAttachPoint()const{return 0;}
      Water():DerivedFxnalGrp(WATER,3,0){}
      typedef FindPrimGroup<2,Water,Oxygen,Hydrogen,Hydrogen> FindMe;
};

class Hydroxyl:public DerivedFxnalGrp{
   public:
      //O first
      Hydroxyl():DerivedFxnalGrp(HYDROXYL,2,1){}
      typedef FindPrimGroup<2,Hydroxyl,Oxygen,Hydrogen> FindMe;

};

class Oxygen2:public DerivedFxnalGrp{
   public:
      //O first
      Oxygen2():DerivedFxnalGrp(OXYGEN2,1,2){}
      typedef FindPrimGroup<2,Oxygen2,Oxygen> FindMe;

};

class OxyDB:public DerivedFxnalGrp{
   public:
      //O first
      OxyDB():DerivedFxnalGrp(OXYDB,1,1){}
      typedef FindPrimGroup<1,OxyDB,Oxygen> FindMe;
};

class Ammonia:public DerivedFxnalGrp{
   public:

      int NAttachPoint()const{return 0;}
      //The N first
      Ammonia():DerivedFxnalGrp(AMMONIA,4,0){}
      typedef FindPrimGroup<3,Ammonia,Nitrogen,Hydrogen,Hydrogen,
                                                        Hydrogen> FindMe;
};

class Amine1:public DerivedFxnalGrp{
   public:
      //The N first
      Amine1():DerivedFxnalGrp(AMINE1,3,1){}
      typedef FindPrimGroup<3,Amine1,Nitrogen,Hydrogen,Hydrogen> FindMe;
};

class Amine2:public DerivedFxnalGrp{
   public:
      //The N first
      Amine2():DerivedFxnalGrp(AMINE2,2,2){}
      typedef FindPrimGroup<3,Amine2,Nitrogen,Hydrogen> FindMe;
};

class Amine3:public DerivedFxnalGrp{
   public:
      //The N first
      Amine3():DerivedFxnalGrp(AMINE3,1,3){}
      typedef FindPrimGroup<3,Amine3,Nitrogen> FindMe;
};

class NitroDB1:public DerivedFxnalGrp{
   public:
      //The N first
      NitroDB1():DerivedFxnalGrp(NITRODB1,2,1){}
      typedef FindPrimGroup<2,NitroDB1,Nitrogen,Hydrogen> FindMe;
};

class NitroDB2:public DerivedFxnalGrp{
   public:
      //The N first
      NitroDB2():DerivedFxnalGrp(NITRODB2,1,2){}
      typedef FindPrimGroup<2,Amine1,Nitrogen> FindMe;
};

class NitroTB:public DerivedFxnalGrp{
   public:
      //The N first
      NitroTB():DerivedFxnalGrp(NITROTB,1,1){}
      typedef FindPrimGroup<1,NitroTB,Nitrogen> FindMe;
};

class HydrogenFluoride: public DerivedFxnalGrp{
   public:
      int NAttachPoint()const{return 0;}
      HydrogenFluoride():DerivedFxnalGrp(HYDROGENFLUORIDE,2,0){}
      typedef FindPrimGroup<1,HydrogenFluoride,Fluorine,Hydrogen> FindMe;
};

class Fluorine1: public DerivedFxnalGrp{
   public:
      //The F first
      Fluorine1():DerivedFxnalGrp(FLUORINE,1,1){}
      typedef FindPrimGroup<1,Fluorine1,Fluorine> FindMe;
};

class HydrogenChloride: public DerivedFxnalGrp{
   public:
      //The Cl first
      HydrogenChloride():DerivedFxnalGrp(HYDROGENCHLORIDE,2,0){}
      int NAttachPoint()const{return 0;}
      typedef FindPrimGroup<1,HydrogenChloride,Chlorine,Hydrogen> FindMe;
};

class Chlorine1: public DerivedFxnalGrp{
   public:
      //The Cl first
      Chlorine1():DerivedFxnalGrp(CHLORINE1,1,1){}
      typedef FindPrimGroup<1,Chlorine1,Chlorine,Hydrogen> FindMe;
};

class HydrogenBromide: public DerivedFxnalGrp{
   public:
      //The Br first
      HydrogenBromide():DerivedFxnalGrp(HYDROGENBROMIDE,2,0){}
      int NAttachPoint()const{return 0;}
      typedef FindPrimGroup<1,HydrogenBromide,Bromine,Hydrogen> FindMe;
};

class Bromine1: public DerivedFxnalGrp{
   public:
      //The Br first
      Bromine1():DerivedFxnalGrp(BROMINE1,1,1){}
      typedef FindPrimGroup<1,Bromine1,Bromine> FindMe;
};

class HydrogenIodide: public DerivedFxnalGrp{
   public:
      //The I first
      HydrogenIodide():DerivedFxnalGrp(HYDROGENIODIDE,2,0){}
      int NAttachPoint()const{return 0;}
      typedef FindPrimGroup<1,HydrogenIodide,Iodine,Hydrogen> FindMe;
};

class Iodine1: public DerivedFxnalGrp{
   public:
      //The I first
      Iodine1():DerivedFxnalGrp(IODINE1,1,1){}
      typedef FindPrimGroup<1,Iodine1,Iodine> FindMe;
};


}}//End namespaces

#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_PRIMITIVEGROUPS_HH_ */
