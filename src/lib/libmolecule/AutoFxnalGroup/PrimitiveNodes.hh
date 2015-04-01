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
typedef FxnGrpType F_t;
/*********Carbon Primitives********************/
class Methane:public RadialSearch<4,C,H,H,H,H>{
   private:
      typedef RadialSearch<4,C,H,H,H,H> Base_t;
   public:
      Methane():Base_t(F_t::C0,F_t::CC0,F_t::HC0,F_t::HC0,F_t::HC0,F_t::HC0){}
};

class Methyl:public RadialSearch<4,C,H,H,H>{
   private:
      typedef RadialSearch<4,C,H,H,H> Base_t;
   public:
      Methyl():Base_t(F_t::C1,F_t::CC1,F_t::HC1,F_t::HC1,F_t::HC1){}
};

class Methene:public RadialSearch<4,C,H,H>{
   private:
      typedef RadialSearch<4,C,H,H> Base_t;
   public:
      Methene():Base_t(F_t::C2,F_t::CC2,F_t::HC2,F_t::HC2){}
};

class Methyne:public RadialSearch<4,C,H>{
   private:
      typedef RadialSearch<4,C,H> Base_t;
   public:
      Methyne():Base_t(F_t::C3,F_t::CC3,F_t::HC3){}
};

class C4:public RadialSearch<4,C>{
   private:
      typedef RadialSearch<4,C> Base_t;
   public:
      C4():Base_t(F_t::C4,F_t::CC4){}
};

class Alkenyl1: public RadialSearch<3,C,H,H>{
   private:
      typedef RadialSearch<3,C,H,H> Base_t;
   public:
      Alkenyl1():Base_t(F_t::CDB1,F_t::CCDB1,F_t::HCDB1,F_t::HCDB1){}
};

class Alkenyl2: public RadialSearch<3,C,H>{
   private:
      typedef RadialSearch<3,C,H> Base_t;
   public:
      Alkenyl2():Base_t(F_t::CDB2,F_t::CCDB2,F_t::HCDB2){}
};

class Alkenyl3:public RadialSearch<3,C>{
   private:
      typedef RadialSearch<3,C> Base_t;
   public:
      Alkenyl3():Base_t(F_t::CDB3,F_t::CCDB3){}
};

class Alkynyl1:public RadialSearch<2,C,H>{
   private:
      typedef RadialSearch<2,C,H> Base_t;
   public:
      Alkynyl1():Base_t(F_t::CTB1,F_t::CCTB1,F_t::HCTB1){}
};

class Alkynyl2:public RadialSearch<2,C>{
   private:
      typedef RadialSearch<2,C> Base_t;
   public:
      Alkynyl2():Base_t(F_t::CTB2,F_t::CCTB2){}
};

/**********************Nitrogen Primitives*******************/
class Ammonia:public RadialSearch<3,N,H,H,H>{
   private:
      typedef RadialSearch<3,N,H,H,H> Base_t;
   public:
      Ammonia():Base_t(F_t::N0,F_t::NN0,F_t::HN0,F_t::HN0,F_t::HN0){}
};

class Amine1:public RadialSearch<3,N,H,H>{
   private:
      typedef RadialSearch<3,N,H,H> Base_t;
   public:
      Amine1():Base_t(F_t::N1,F_t::NN1,F_t::HN1,F_t::HN1){}
};

class Amine2:public RadialSearch<3,N,H>{
   private:
      typedef RadialSearch<3,N,H> Base_t;
   public:
      Amine2():Base_t(F_t::N2,F_t::NN2,F_t::HN2){}
};

class Amine3:public RadialSearch<3,N>{
   private:
      typedef RadialSearch<3,N> Base_t;
   public:
      Amine3():Base_t(F_t::N3,F_t::NN3){}
};

class Azo1:public RadialSearch<2,N,H>{
   private:
      typedef RadialSearch<2,N,H> Base_t;
   public:
      Azo1():Base_t(F_t::NDB1,F_t::NNDB1,F_t::HNDB1){}
};

class Azo2:public RadialSearch<2,N>{
   private:
      typedef RadialSearch<2,N> Base_t;
   public:
      Azo2():Base_t(F_t::NDB2,F_t::NNDB2){}
};

class NTB:public RadialSearch<1,N>{
   private:
      typedef RadialSearch<1,N> Base_t;
   public:
      NTB():Base_t(F_t::NTB,F_t::NNTB){}
};

class Ammonium:public RadialSearch<4,N,H,H,H,H>{
   private:
      typedef RadialSearch<4,N,H,H,H,H> Base_t;
   public:
      Ammonium():Base_t(F_t::NP0,F_t::NNP0,F_t::HNP0,F_t::HNP0,
                                           F_t::HNP0,F_t::HNP0){}
};

class Ammonium1:public RadialSearch<4,N,H,H,H>{
   private:
      typedef RadialSearch<4,N,H,H,H> Base_t;
   public:
      Ammonium1():Base_t(F_t::NP1,F_t::NNP1,F_t::HNP1,F_t::HNP1,
                                           F_t::HNP1){}
};

class Ammonium2:public RadialSearch<4,N,H,H>{
   private:
      typedef RadialSearch<4,N,H,H> Base_t;
   public:
      Ammonium2():Base_t(F_t::NP2,F_t::NNP2,F_t::HNP2,F_t::HNP2){}
};

class Ammonium3:public RadialSearch<4,N,H>{
   private:
      typedef RadialSearch<4,N,H> Base_t;
   public:
      Ammonium3():Base_t(F_t::NP3,F_t::NNP3,F_t::HNP3){}
};

class Ammonium4:public RadialSearch<4,N>{
   private:
      typedef RadialSearch<4,N> Base_t;
   public:
      Ammonium4():Base_t(F_t::NP4,F_t::NNP4){}
};

/**********************Oxygen Primitives**********************/
class Water:public RadialSearch<2,O,H,H>{
   private:
      typedef RadialSearch<2,O,H,H> Base_t;
   public:
      Water():Base_t(F_t::O0,F_t::OO0,F_t::HO0,F_t::HO0){}
};

class Hydroxyl:public RadialSearch<2,O,H>{
   private:
      typedef RadialSearch<2,O,H> Base_t;
   public:
      Hydroxyl():Base_t(F_t::O1,F_t::OO1,F_t::HO1){}
};

class Ether:public RadialSearch<2,O>{
   private:
      typedef RadialSearch<2,O> Base_t;
   public:
      Ether():Base_t(F_t::O2,F_t::OO2){}
};

class ODB:public RadialSearch<1,O>{
   private:
      typedef RadialSearch<1,O> Base_t;
   public:
      ODB():Base_t(F_t::ODB,F_t::OODB){}
};

/*********************Sulfur primitives**************************/

class HydrogenSulfide:public RadialSearch<3,S,H,H>{
   private:
      typedef RadialSearch<3,S,H,H> Base_t;
   public:
      HydrogenSulfide():Base_t(F_t::S0,F_t::SS0,F_t::HS0,F_t::HS0){}
};

class Thiol:public RadialSearch<3,S,H>{
   private:
      typedef RadialSearch<3,S,H> Base_t;
   public:
      Thiol():Base_t(F_t::S1,F_t::SS1,F_t::HS1){}
};

class Sulfide:public RadialSearch<3,S>{
   private:
      typedef RadialSearch<3,S> Base_t;
   public:
      Sulfide():Base_t(F_t::S2,F_t::SS2){}
};

class SDB:public RadialSearch<2,S>{
   private:
      typedef RadialSearch<2,S> Base_t;
   public:
      SDB():Base_t(F_t::SDB,F_t::SSDB){}
};
}}//End namespaces
#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_PRIMITIVENODES_HH_ */
