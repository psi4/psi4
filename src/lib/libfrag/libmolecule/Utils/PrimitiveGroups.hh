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

class Methane:public GenPrimGrp<4,METHANE,Methane,0,Carbon,
	Hydrogen,Hydrogen,Hydrogen,Hydrogen>{
   public:
      int NAttachPoint()const{return 0;}
};

class Methyl:public GenPrimGrp<4,METHYL,Methyl,1,Carbon,
	Hydrogen,Hydrogen,Hydrogen>{};

class Methene:public GenPrimGrp<4,METHENE,Methene,2,Carbon,
	Hydrogen,Hydrogen>{};

class Methyne:public GenPrimGrp<4,METHYNE,Methyne,3,Carbon,Hydrogen>{};

class Carbon4:public GenPrimGrp<4,CARBON4,Carbon4,4,Carbon>{};

class Alkenyl1:public GenPrimGrp<3,ALKENYL1,Alkenyl1,1,Carbon,Hydrogen,
		Hydrogen> {};

class Alkenyl2:public GenPrimGrp<3,ALKENYL2,Alkenyl2,2,Carbon,Hydrogen>{};

class Alkenyl3:public GenPrimGrp<3,ALKENYL3,Alkenyl3,3,Carbon>{};

class Alkynyl1:public GenPrimGrp<2,ALKYNYL1,Alkynyl1,1,Carbon,Hydrogen>{};

class Alkynyl2:public GenPrimGrp<2,ALKYNYL2,Alkynyl2,2,Carbon>{};

class Water:public GenPrimGrp<2,WATER,Water,0,Oxygen,Hydrogen,Hydrogen>{
   public:
      int NAttachPoint()const{return 0;}
};

class Hydroxyl:public GenPrimGrp<2,HYDROXYL,Hydroxyl,1,Oxygen,Hydrogen>{};

class Oxygen2:public GenPrimGrp<2,OXYGEN2,Oxygen2,2,Oxygen>{};

class OxyDB:public GenPrimGrp<1,OXYDB,OxyDB,1,Oxygen> {};

class Ammonia:public GenPrimGrp<3,AMMONIA,Ammonia,0,Nitrogen,Hydrogen,Hydrogen,
Hydrogen>{
   public:
      int NAttachPoint()const{return 0;}
};

class Amine1:public GenPrimGrp<3,AMINE1,Amine1,1,Nitrogen,Hydrogen,
	Hydrogen>{};

class Amine2:public GenPrimGrp<3,AMINE2,Amine2,2,Nitrogen,Hydrogen> {};

class Amine3:public GenPrimGrp<3,AMINE3,Amine3,3,Nitrogen> {};

class NitroDB1:public GenPrimGrp<2,NITRODB1,NitroDB1,1,Nitrogen,Hydrogen>{};

class NitroDB2:public GenPrimGrp<2,NITRODB2,Amine1,2,Nitrogen>{};

class NitroTB:public GenPrimGrp<1,NITROTB,NitroTB,1,Nitrogen>{};

class HydrogenFluoride: public GenPrimGrp<1,HYDROGENFLUORIDE,HydrogenFluoride,
0,Fluorine,Hydrogen>{
   public:
      int NAttachPoint()const{return 0;}
};

class Fluorine1: public GenPrimGrp<1,FLUORINE,Fluorine1,1,Fluorine>{};

class HydrogenChloride: public GenPrimGrp<1,HYDROGENCHLORIDE,
	HydrogenChloride,0,Chlorine,Hydrogen>{
   public:
      int NAttachPoint()const{return 0;}
};

class Chlorine1: public GenPrimGrp<1,CHLORINE1,Chlorine1,1,Chlorine> {};

class HydrogenBromide: public GenPrimGrp<1,HYDROGENBROMIDE,
HydrogenBromide,0,Bromine,Hydrogen>{
   public:
      int NAttachPoint()const{return 0;}
};

class Bromine1: public GenPrimGrp<1,BROMINE1,Bromine1,1,Bromine> {};

class HydrogenIodide: public GenPrimGrp<1,HYDROGENIODIDE,
HydrogenIodide,0,Iodine,Hydrogen>{
public:
	int NAttachPoint()const{return 0;}
};

class Iodine1: public GenPrimGrp<1,IODINE1,Iodine1,1,Iodine>{};

}}//End namespaces

#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_PRIMITIVEGROUPS_HH_ */
