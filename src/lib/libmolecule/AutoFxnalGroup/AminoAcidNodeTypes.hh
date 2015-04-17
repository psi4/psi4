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
#define DERIV(Name,Abbv,Full,Centers...) \
   template<size_t...args>\
   class base_##Name:public LinearSearch<Centers>{\
      private:\
         typedef LinearSearch<Centers> Base_t;\
      public:\
         base_##Name():Base_t(Abbv,Full){\
            if(sizeof...(args)!=0){\
               std::vector<size_t> Temp{args...};\
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

class NTerm: public Node{
   public:
      NTerm():Node(0,"NT","N-Terminus N"){
         PsiMap<size_t,size_t> Temp;
         Temp[0]=1;
         this->MyTypes_[0]=ParamT("N+","Ammonium",Temp,0);
         MyTypes_.push_back(ParamT("N","Amine",Temp,0));
      }
};
class CTerm: public Node{
   public:
      CTerm():Node(0,"CT","C-Terminus C"){
         PsiMap<size_t,size_t> Temp;
         Temp[0]=1;
         this->MyTypes_[0]=ParamT("CO2-","Carboxylate",Temp,0);
         MyTypes_.push_back(ParamT("CO2H","Carboxyl",Temp,0));
      }
};

class BB: public Node{
   public:
      BB():Node(0,"BB","The non-R part of an AA"){
         PsiMap<size_t,size_t> Temp;
         Temp[0]=1;
         this->MyTypes_[0]=ParamT("AA","An amino-acid",Temp,0);
         Temp[0]++;
         MyTypes_.push_back(ParamT("CT","C-Terminus",Temp,0));
         MyTypes_.push_back(ParamT("NT","N-Terminus",Temp,0));
         Temp[0]++;
         MyTypes_.push_back(ParamT("AABB","Amino-acid back-bone",Temp,0));
      }
};

DERIV(CTerminus,"CT","C-Terminus",CTerm,Methyne_t,Amine2_t)
DERIV(NTerminus,"NT","N-Terminus",NTerm,Methyne_t,Carbonyl_t)
DERIV(NTermGly,"NTGLY","N-Terminal Gly",NTerm,Methene_t,Carbonyl_t)
DERIV(CTermGly,"CTGLY","C-Terminal Gly",CTerm,Methene_t,Amine2_t)
DERIV(Gly,"GLY","Gly",Amine2_t,Methene_t,Carbonyl_t)
DERIV(AAGly,"AGLY","A Single Gly",CTerm,Methene_t,NTerm)
DERIV(CTermPro,"CTPRO","C-Terminal Pro",CTerm,Methyne_t,Methene_t,Methene_t,Methene_t,Amine3_t)
DERIV(NPTermPro,"N+TPRO","Positive N-Terminal Pro",Carbonyl_t,Methyne_t,Methene_t,Methene_t,Methene_t,Ammonium2_t)
DERIV(NTermPro,"NTPRO","N-Terminal Pro",Carbonyl_t,Methyne_t,Methene_t,Methene_t,Methene_t,Amine2_t)
DERIV(Pro,"PRO","Pro",Carbonyl_t,Methyne_t,Methene_t,Methene_t,Methene_t,Amine3_t)
DERIV(AABB,"AABB","Amino-acid back-bone",Amine2_t,Methyne_t,Carbonyl_t)
DERIV(AnAA,"AA","An amino-acid",CTerm,Methyne_t,NTerm)

DERIV(Ala,"ALA","Alanine",BB,Methyl_t)
DERIV(Val,"VAL","Valine",BB,Isopropyl_t)
DERIV(Ile,"ILE","Isoleucine",BB,Butyl_t)
DERIV(Leu,"LEU","Leucine",BB,Methene_t,Isopropyl_t)
DERIV(Met,"MET","Methionine",BB,Methene_t,Methene_t,Sulfide_t,Methyl_t)
DERIV(Phe,"PHE","Phenylalanine",BB,Methene_t,base_Benzene<0,1>)
DERIV(Tyr,"TYR","Tyrosine",BB,Methene_t,base_Benzene<0,1,3,1>,Hydroxyl_t)
DERIV(Trp,"TRP","Tryptophan",BB,Methene_t,base_Indole<0,1>)
DERIV(Ser,"SER","Serine",BB,Methene_t,Hydroxyl_t)
DERIV(Thr,"THR","Threonine",BB,Hydroxyethyl_t)
DERIV(Asn,"ASN","Asparagine",BB,Methene_t,Carbonyl_t,Amine1_t)
DERIV(Gln,"GLN","Glutamine",BB,Methene_t,Methene_t,Carbonyl_t,Amine1_t)
DERIV(Cys,"CYS","Cysteine",BB,Methene_t,Thiol_t)
DERIV(CysSB,"SBCYS","Sulfide bond Cysteine",BB,Methene_t,Sulfide_t)
DERIV(Arg,"ARG","Arginine",BB,Methene_t,Methene_t,Methene_t,Amine2_t,ArgThing_t)
DERIV(His,"HIS","Histidine",BB,Methene_t,base_Imidazole<0,1>)
DERIV(Hip,"HIP","Protonated Histidine",BB,Methene_t,base_DBCC<0,2,1,1>,Amine2_t,Alkenyl2_t,Amine2_t)
DERIV(Lys,"LYS","Lysine",BB,Methene_t,Methene_t,Methene_t,Methene_t,Ammonium1_t)
DERIV(Asp,"ASP","Aspratic Acid",BB,Methene_t,Carboxylate_t)
DERIV(Glu,"GLU","Glutamic Acid",BB,Methene_t,Methene_t,Carboxylate_t)
#undef DERIV

/*
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

*/
}}//End namespaces

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_AMINOACIDNODETYPES_HH_ */
