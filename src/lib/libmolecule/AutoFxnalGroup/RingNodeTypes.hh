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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGNODETYPES_HH_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGNODETYPES_HH_
#include "RingFinder.h"

namespace psi{
namespace LibMolecule{
#define RING(Name,Abbv,Full,Centers...)\
template<size_t...args>\
class base_##Name:public RingFinder<Centers>{\
   private:\
      typedef RingFinder<Centers> Base_t;\
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

/******Fused 5 and 6 Membered Rings**********************/
RING(Benzofuran,"Ar9O","Benzofuran",Ether_t,Alkenyl3_t,Alkenyl,Alkenyl,Alkenyl,Alkenyl,Alkenyl3_t,Alkenyl,Alkenyl)
RING(Isobenzofuran,"Ar9IO","Isobenzofuran",Ether_t,Alkenyl,Alkenyl3_t,Alkenyl,Alkenyl,Alkenyl,Alkenyl,Alkenyl3_t,Alkenyl)
RING(Indole,"Ar9N","Indole",Amine,Alkenyl3_t,Alkenyl,Alkenyl,Alkenyl,Alkenyl,Alkenyl3_t,Alkenyl,Alkenyl)
RING(Isoindole,"Ar9IN","Isoindole",Amine,Alkenyl,Alkenyl3_t,Alkenyl,Alkenyl,Alkenyl,Alkenyl,Alkenyl3_t,Alkenyl)
RING(Purine,"Ar9N4","Purine",Amine,Alkenyl3_t,Azo,Alkenyl,Azo,Alkenyl,Alkenyl3_t,Azo,Alkenyl)
RING(Indazole,"Ar9N2","Indazole",Azo,Amine,Alkenyl3_t,Alkenyl,Alkenyl,Alkenyl,Alkenyl,Alkenyl3_t,Alkenyl)
RING(Benzoxazole,"Ar9ON","Benzoxazole",Ether_t,Alkenyl3_t,Alkenyl,Alkenyl,Alkenyl,Alkenyl,Alkenyl3_t,Azo,Alkenyl)
RING(Benzisoxazole,"Ar9ION","Benzisoxazole",Ether_t,Alkenyl3_t,Alkenyl,Alkenyl,Alkenyl,Alkenyl,Alkenyl3_t,Alkenyl,Azo)


/***********5 Membered Rings****************************/
RING(Furan,"Ar5O","Furan",Ether_t,Alkenyl,Alkenyl,Alkenyl,Alkenyl)
RING(Pyrrole,"Ar5N","Pyrrole",Amine,Alkenyl,Alkenyl,Alkenyl,Alkenyl)
RING(Thiophene,"Ar5S","Thiophene",Sulfide_t,Alkenyl,Alkenyl,Alkenyl,Alkenyl)
RING(Imidazole,"Ar5NCN","Imidazole",Amine,Alkenyl,Azo2_t,Alkenyl,Alkenyl)
RING(Pyrazole,"Ar5NN","Pyrazole",Amine,Azo2_t,Alkenyl,Alkenyl,Alkenyl)
RING(Oxazole,"Ar5OCN","Oxazole",Ether_t,Alkenyl,Azo2_t,Alkenyl,Alkenyl)
RING(Isoxazole,"Ar5ON","Isooxazole",Ether_t,Azo2_t,Alkenyl,Alkenyl,Alkenyl)

/*************************6 Membered Rings*****************/
RING(Benzene,"Ph","Benzene",Alkenyl,Alkenyl,Alkenyl,Alkenyl,Alkenyl,Alkenyl)
RING(Pyridine,"Ar6N","Pyridine",Azo,Alkenyl,Alkenyl,Alkenyl,Alkenyl,Alkenyl)
RING(Pyrazine,"Ar6pN2","Pyrazine",Azo,Alkenyl,Alkenyl,Azo,Alkenyl,Alkenyl)
RING(Pyrimidine,"Ar6mN2","Pyrimidine",Azo,Alkenyl,Azo,Alkenyl,Alkenyl,Alkenyl)
RING(Pyridazine,"Ar6oN2","Pyridazine",Azo,Azo,Alkenyl,Alkenyl,Alkenyl,Alkenyl)

#undef RING

}}//End namespaces

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGNODETYPES_HH_ */
