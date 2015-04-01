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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_AMINOACIDNODETYPES_HH_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_AMINOACIDNODETYPES_HH_

#include "RingNodeTypes.hh"

namespace psi{
namespace LibMolecule{

class PNTerminus: public LinearSearch<Ammonium1,Methyne,Carbonyl>{
   private:
      typedef LinearSearch<Ammonium1,Methyne,Carbonyl> Base_t;
   public:
      PNTerminus(): Base_t(F_t::PNT,F_t::NPNT,F_t::HNPNT,
                           F_t::HNPNT,F_t::HNPNT,
                         F_t::CAPNT,F_t::HCAPNT,F_t::CPBPNT,F_t::OPBPNT){}
};

class NCTerminus: public LinearSearch<Carboxylate,Methyne,Amine2>{
   private:
      typedef LinearSearch<Carboxylate,Methyne,Amine2> Base_t;
   public:
      NCTerminus(): Base_t(F_t::NCT,F_t::CNCT,F_t::ONCT,F_t::ONCT,
                           F_t::CANCT,F_t::HCANCT,F_t::NPBNCT,F_t::HPBNCT){}
};

class PNTermGly: public LinearSearch<Ammonium1,Methene,Carbonyl>{
   private:
      typedef LinearSearch<Ammonium1,Methene,Carbonyl> Base_t;
   public:
      PNTermGly():Base_t(F_t::GLYPNT,
           F_t::NGLYPNT,F_t::HNGLYPNT,F_t::HNGLYPNT,F_t::HNGLYPNT,
           F_t::CAGLYPNT,F_t::HCAGLYPNT,F_t::HCAGLYPNT,
           F_t::CPBGLYPNT,F_t::OPBGLYPNT){}
};

class NCTermGly: public LinearSearch<Carboxylate,Methene,Amine2>{
   private:
      typedef LinearSearch<Carboxylate,Methene,Amine2> Base_t;
   public:
      NCTermGly():Base_t(F_t::GLYNCT,
            F_t::CGLYNCT,F_t::OGLYNCT,F_t::OGLYNCT,
            F_t::CAGLYNCT,F_t::HCAGLYNCT,F_t::HCAGLYNCT,
            F_t::NGLYNCT,F_t::HNGLYNCT
           ){}
};

class Glycine: public LinearSearch<Carbonyl,Methene,Amine2>{
   private:
      typedef LinearSearch<Carbonyl,Methene,Amine2> Base_t;
   public:
      Glycine():Base_t(F_t::GLY,F_t::CPBGLY,F_t::OPBGLY,
                       F_t::CCAGLY,F_t::HCAGLY,F_t::HCAGLY,
                       F_t::NPBGLY,F_t::HPBGLY){}
};

class AABackBone: public LinearSearch<Carbonyl,Methyne,Amine2>{
   private:
      typedef LinearSearch<Carbonyl,Methyne,Amine2> Base_t;
   public:
      AABackBone():Base_t(F_t::AABB,F_t::CPBAABB,F_t::OPBAABB,
                          F_t::CCAAABB,F_t::HCAAABB,F_t::NPBAABB,F_t::HPBAABB){}
};

class AlanineR:public LinearSearch<Methyl>{
   private:
      typedef LinearSearch<Methyl> Base_t;
   public:
      AlanineR():Base_t(F_t::ALAR,
            F_t::CALAR,F_t::HALAR,F_t::HALAR,F_t::HALAR){}
};

class ValineR: public LinearSearch<Isopropyl>{
   private:
      typedef LinearSearch<Isopropyl> Base_t;
   public:
      ValineR():Base_t(F_t::VALR,F_t::C2VALR,F_t::H2VALR,F_t::H2VALR,F_t::H2VALR,
            F_t::C1VALR,F_t::H1VALR,
            F_t::C2VALR,F_t::H2VALR,F_t::H2VALR,F_t::H2VALR
            ){}
};

class IsoleucineR: public LinearSearch<Methyl,Methyne,Methene,Methyl>{
   private:
      typedef LinearSearch<Methyl,Methyne,Methene,Methyl> Base_t;
   public:
      IsoleucineR():Base_t(F_t::ILER,
                           F_t::C3ILER,F_t::H3ILER,F_t::H3ILER,F_t::H3ILER,
                           F_t::C1ILER,F_t::H1ILER,
                           F_t::C2ILER,F_t::H2ILER,F_t::H2ILER,
                           F_t::C4ILER,F_t::H4ILER,F_t::H4ILER,F_t::H4ILER){}
};

class LeucineR: public LinearSearch<Isopropyl,Methyne>{
   private:
      typedef LinearSearch<Isopropyl,Methyne> Base_t;
   public:
      LeucineR():Base_t(F_t::LEUR,
            F_t::C3LEUR,F_t::H3LEUR,F_t::H3LEUR,F_t::H3LEUR,
            F_t::C2LEUR,F_t::H2LEUR,
            F_t::C3LEUR,F_t::H3LEUR,F_t::H3LEUR,F_t::H3LEUR,
            F_t::C1LEUR,F_t::H1LEUR
            ){}
};

class MethionineR: public LinearSearch<Methyl,Sulfide,Methene,Methene>{
   private:
      typedef LinearSearch<Methyl,Sulfide,Methene,Methene> Base_t;
   public:
      MethionineR():Base_t(F_t::METR,
            F_t::C3METR,F_t::H3METR,F_t::H3METR,F_t::H3METR,
            F_t::SMETR,
            F_t::C2METR,F_t::H2METR,F_t::H2METR,
            F_t::C1METR,F_t::H1METR,F_t::H1METR){}
};

class PhenylalanineR: public LinearSearch<Phenyl,Methene>{
   private:
      typedef LinearSearch<Phenyl,Methene> Base_t;
   public:
      PhenylalanineR():Base_t(F_t::PHER,
            F_t::C2PHER,F_t::C3PHER,F_t::H2PHER,F_t::C4PHER,F_t::H3PHER,
            F_t::C5PHER,F_t::H4PHER,F_t::C4PHER,F_t::H3PHER,
            F_t::C3PHER,F_t::H2PHER,F_t::C1PHER,F_t::H1PHER,F_t::H1PHER){}
};

class TyrosineR: public LinearSearch<Hydroxyl,ParaBenzene,Methene>{
   private:
      typedef LinearSearch<Hydroxyl,ParaBenzene,Methene> Base_t;
   public:
      TyrosineR():Base_t(F_t::TYRR,
            F_t::OTYRR,F_t::H4TYRR,F_t::C2TYRR,
            F_t::C3TYRR,F_t::H2TYRR,F_t::C4TYRR,F_t::H3TYRR,
            F_t::C5TYRR,F_t::C4TYRR,F_t::H3TYRR,F_t::H2TYRR,
            F_t::C1TYRR,F_t::H1TYRR,F_t::H1TYRR){}
};

class TryptophanR: public LinearSearch<Indolyl1_5,Methene>{
   private:
      typedef LinearSearch<Indolyl1_5,Methene> Base_t;
   public:
      TryptophanR():Base_t(F_t::TRPR,
      F_t::NTRPR,F_t::H2TRPR,F_t::C2TRPR,F_t::H3TRPR,F_t::C3TRPR,
      F_t::C4TRPR,F_t::C5TRPR,F_t::H4TRPR,F_t::C6TRPR,F_t::H5TRPR,
      F_t::C7TRPR,F_t::H6TRPR,F_t::C8TRPR,F_t::H7TRPR,F_t::C9TRPR,
      F_t::C1TRPR,F_t::H1TRPR,F_t::H1TRPR){}
};

class SerineR: public LinearSearch<Hydroxyl,Methene>{
   private:
      typedef LinearSearch<Hydroxyl,Methene> Base_t;
   public:
      SerineR():Base_t(F_t::SERR,
            F_t::OSERR,F_t::H2SERR,F_t::CSERR,F_t::H1SERR,F_t::H1SERR
            ){}
};


}}//End namespaces

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_AMINOACIDNODETYPES_HH_ */
