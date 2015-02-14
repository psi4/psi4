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
namespace psi{
namespace LibMolecule{
/*********** Primitive Groups ***************************************/

class Carbon:public FxnalGroup{
   public:
     Carbon(const int* Members):FxnalGroup(CARBON,0){
        Members_.push_back(-1);
     }
};

class Oxygen:public FxnalGroup{
   public:
     Oxygen(const int* Members):FxnalGroup(OXYGEN,0){
        Members_.push_back(-1);
     }
};

class Nitrogen:public FxnalGroup{
   public:
     Nitrogen(const int* Members):FxnalGroup(NITROGEN,0){
        Members_.push_back(-1);
     }
};

class Hydrogen:public FxnalGroup{
   public:
     Hydrogen(const int* Members):FxnalGroup(HYDROGEN,0){
        Members_.push_back(-1);
     }
};

class Fluorine:public FxnalGroup{
   public:
     Fluorine(const int* Members):FxnalGroup(FLUORINE,0){
        Members_.push_back(-1);
     }
};

class Chlorine:public FxnalGroup{
   public:
     Chlorine(const int* Members):FxnalGroup(CHLORINE,0){
        Members_.push_back(-1);
     }
};

class Bromine:public FxnalGroup{
   public:
     Bromine(const int* Members):FxnalGroup(BROMINE,0){
        Members_.push_back(-1);
     }
};

class Iodine:public FxnalGroup{
   public:
     Iodine(const int* Members):FxnalGroup(IODINE,0){
        Members_.push_back(-1);
     }
};


class Methane:public FxnalGroup{
   public:
      //C first
      Methane(const int* Members):FxnalGroup(METHANE,0,5,Members){}
      int NAttachPoint()const{return 0;}

};

class Methyl:public FxnalGroup{
   public:
      //C first
      Methyl(const int* Members):FxnalGroup(METHYL,1,4,Members){}
};

class Methene:public FxnalGroup{
   public:
      //C first
      Methene(const int* Members):FxnalGroup(METHENE,2,3,Members){}
};

class Methyne:public FxnalGroup{
   public:
      //C first
      Methyne(const int* Members):FxnalGroup(METHYNE,3,2,Members){}
};

class Carbon4:public FxnalGroup{
   public:
      //Only has a C
      Carbon4(const int* Members):FxnalGroup(CARBON4,4,1,Members){}
};

class Alkenyl1:public FxnalGroup{
   public:
      //C first
      Alkenyl1(const int* Members):FxnalGroup(ALKENYL1,1,3,Members){}
};

class Alkenyl2:public FxnalGroup{
   public:
      //C first
      Alkenyl2(const int* Members):FxnalGroup(ALKENYL2,2,2,Members){}
};

class Alkenyl3:public FxnalGroup{
   public:
      //C first
      Alkenyl3(const int* Members):FxnalGroup(ALKENYL3,3,1,Members){}
};

class Alkynyl1:public FxnalGroup{
   public:
      //C first
      Alkynyl1(const int* Members):FxnalGroup(ALKYNYL1,1,2,Members){}
};

class Alkynyl2:public FxnalGroup{
   public:
      //C first
      Alkynyl2(const int* Members):FxnalGroup(ALKYNYL2,2,1,Members){}
};

class Water:public FxnalGroup{
   public:
      //O first
      Water(const int* Members):FxnalGroup(WATER,0,3,Members){}
      int NAttachPoint()const{return 0;}

};

class Hydroxyl:public FxnalGroup{
   public:
      //O first
      Hydroxyl(const int * Members):FxnalGroup(HYDROXYL,1,2,Members){}
};

class Oxygen2:public FxnalGroup{
   public:
      //The O
      Oxygen2(const int *Members):FxnalGroup(OXYGEN2,2,1,Members){}
};

class OxyDB:public FxnalGroup{
   public:
      OxyDB(const int* Members):FxnalGroup(OXYDB,1,1,Members){}
};

class Ammonia:public FxnalGroup{
   public:
      //The N first
      Ammonia(const int* Members):FxnalGroup(AMMONIA,0,4,Members){}
      int NAttachPoint()const{return 0;}
};

class Amine1:public FxnalGroup{
   public:
      //The N first
      Amine1(const int *Members):FxnalGroup(AMINE1,1,3,Members){}
};

class Amine2:public FxnalGroup{
   public:
      //The N first
      Amine2(const int *Members):FxnalGroup(AMINE2,2,2,Members){}
};

class Amine3:public FxnalGroup{
   public:
      //The N first
      Amine3(const int *Members):FxnalGroup(AMINE3,3,1,Members){}
};

class NitroDB1:public FxnalGroup{
   public:
      //The N first
      NitroDB1(const int* Members):FxnalGroup(NITRODB1,1,2,Members){}
};

class NitroDB2:public FxnalGroup{
   public:
      //The N first
      NitroDB2(const int* Members):FxnalGroup(NITRODB2,2,1,Members){}
};

class NitroTB:public FxnalGroup{
   public:
      //The N first
      NitroTB(const int* Members):FxnalGroup(NITROTB,1,1,Members){}
};

class HydrogenFluoride: public FxnalGroup{
   public:
      //The F first
      HydrogenFluoride(const int *Members):
         FxnalGroup(HYDROGENFLUORIDE,0,2,Members){}
      int NAttachPoint()const{return 0;}
};

class Fluorine1: public FxnalGroup{
   public:
      //The F first
      Fluorine1(const int *Members):
         FxnalGroup(FLUORINE1,1,1,Members){}
};

class HydrogenChloride: public FxnalGroup{
   public:
      //The Cl first
      HydrogenChloride(const int *Members):
         FxnalGroup(HYDROGENCHLORIDE,0,2,Members){}
      int NAttachPoint()const{return 0;}
};

class Chlorine1: public FxnalGroup{
   public:
      //The Cl first
      Chlorine1(const int *Members):
         FxnalGroup(CHLORINE1,1,1,Members){}
};

class HydrogenBromide: public FxnalGroup{
   public:
      //The Br first
      HydrogenBromide(const int *Members):
         FxnalGroup(HYDROGENBROMIDE,0,2,Members){}
      int NAttachPoint()const{return 0;}
};

class Bromine1: public FxnalGroup{
   public:
      //The Br first
      Bromine1(const int *Members):
         FxnalGroup(BROMINE1,1,1,Members){}
};

class HydrogenIodide: public FxnalGroup{
   public:
      //The I first
      HydrogenIodide(const int *Members):
         FxnalGroup(HYDROGENIODIDE,0,2,Members){}
      int NAttachPoint()const{return 0;}
};

class Iodine1: public FxnalGroup{
   public:
      //The I first
      Iodine1(const int *Members):
      FxnalGroup(IODINE1,1,1,Members){}
};


}}//End namespaces

#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_PRIMITIVEGROUPS_HH_ */
