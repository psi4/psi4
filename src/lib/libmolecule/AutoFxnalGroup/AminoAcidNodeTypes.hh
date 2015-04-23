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
      protected:\
         bool IsPrim()const{return true;}\
      public:\
         base_##Name():Base_t(Abbv,Full){\
            if(sizeof...(args)!=0){\
               size_t Temp[sizeof...(args)]={args...};\
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
#define PRIM(Name,Abbv,Full,Centers...) \
   template<size_t...args>\
   class base_##Name:public LinearSearch<Centers>{\
      private:\
         typedef LinearSearch<Centers> Base_t;\
      protected:\
         bool IsPrim()const{return true;}\
      public:\
         base_##Name():Base_t(Abbv,Full){\
            if(sizeof...(args)!=0){\
               size_t Temp[sizeof...(args)]={args...};\
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
         Temp[1]=1;
         MyTypes_.push_back(ParamT("CT","C-Terminus",Temp,0));
         MyTypes_.push_back(ParamT("NT","N-Terminus",Temp,0));
         Temp[2]=1;
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
DERIV(NPTermPro,"NTPRO","N-Terminal Pro",Carbonyl_t,Methyne_t,Methene_t,Methene_t,Methene_t,Ammonium2_t)
DERIV(NTermPro,"NNTPRO","Neutral N-Terminal Pro",Carbonyl_t,Methyne_t,Methene_t,Methene_t,Methene_t,Amine2_t)
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
#undef PRIM
}}//End namespaces

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_AMINOACIDNODETYPES_HH_ */
