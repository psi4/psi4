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

class HydrogenCyanide:
		public GenDerGrp<HYDROGENCYANIDE,HydrogenCyanide,0,Alkynyl1,NitroTB>{
   public:
      int NAttachPoint()const{return 0;}
};

class Nitrile:public GenDerGrp<NITRILE,Nitrile,1,Alkynyl2,NitroTB>{};

class Formaldehyde:
		public GenDerGrp<FORMALDEHYDE,Formaldehyde,0,Alkenyl1,OxyDB>{
   public:
      int NAttachPoint()const{return 0;}
};

class Aldehyde:public GenDerGrp<ALDEHYDE,Aldehyde,1,Alkenyl2,OxyDB>{};

class Carbonyl:public GenDerGrp<CARBONYL,Carbonyl,2,Alkenyl3,OxyDB>{};


class Carboxyl:public GenDerGrp<CARBOXYL,Carboxyl,1,Carbonyl,Hydroxyl>{};

class HydroPeroxy:
		public GenDerGrp<HYDROPEROXY,HydroPeroxy,1,Oxygen2,Hydroxyl>{};

class Peroxide:public GenDerGrp<PEROXIDE,Peroxide,0,Hydroxyl,Hydroxyl>{
   public:
      int NAttachPoint()const{return 0;}
};

class Methoxy:public GenDerGrp<METHOXY,Methoxy,1,Oxygen2,Methyl>{};

class Methanol:public GenDerGrp<METHANOL,Methanol,0,Hydroxyl,Methyl>{
   public:
      int NAttachPoint()const{return 0;}
};

class ccdb4:public GenDerGrp<CCDB4,ccdb4,4,Alkenyl3,Alkenyl3>{
   public:
      int NAttachPoint()const{return 2;}
};

class ccdb3:public GenDerGrp<CCDB3,ccdb3,3,Alkenyl2,Alkenyl3>{
   public:
      int NAttachPoint()const{return 2;}
};

class ccdb2:public GenDerGrp<CCDB2,ccdb2,2,Alkenyl2,Alkenyl2>{
   public:
      int NAttachPoint()const{return 2;}
};

class Ethenyl2:public GenDerGrp<ETHENYL2,Ethenyl2,2,Alkenyl3,Alkenyl1>{};

class Ethenyl1:public GenDerGrp<ETHENYL1,Ethenyl1,1,Alkenyl2,Alkenyl1>{};

class Ethene:public GenDerGrp<ETHENE,Ethene,0,Alkenyl1,Alkenyl1>{
   public:
      int NAttachPoint()const{return 0;}
};

class cctb:public GenDerGrp<CCTB,cctb,2,Alkynyl2,Alkynyl2>{
   public:
      int NAttachPoint()const{return 2;}
};

class Ethynyl:public GenDerGrp<ETHYNYL,Ethynyl,1,Alkynyl2,Alkynyl1>{};

class Ethyne:public GenDerGrp<ETHYNE,Ethyne,0,Alkynyl1,Alkynyl1>{
   public:
      int NAttachPoint()const{return 0;}
};

class Ketimine1:
		public GenDerGrp<KETIMINE1,Ketimine1,2,Alkenyl3,NitroDB1>{};

class Ketimine2:
		public GenDerGrp<KETIMINE2,Ketimine2,3,Alkenyl3,NitroDB2>{
   public:
      int NAttachPoint()const{return 2;}
};

class Aldimine1:
		public GenDerGrp<ALDIMINE1,Aldimine1,1,Alkenyl2,NitroDB1>{};

class Aldimine2:
		public GenDerGrp<ALDIMINE2,Aldimine2,2,Alkenyl2,NitroDB2>{
   public:
      int NAttachPoint()const{return 2;}
};

class Methanimine:
		public GenDerGrp<METHANIMINE,Methanimine,0,Alkenyl1,NitroDB1>{
   public:
      int NAttachPoint()const{return 0;}
};

}}
#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_DERIVEDGROUPS_HH_ */
