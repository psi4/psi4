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
#include "MBE.h"
namespace psi{
namespace LibFrag{

typedef LibMolecule::SerialNumber Set_t;
typedef Set_t::const_iterator cSetIt_t;

std::vector<cSetIt_t> CanonicalMBE::FillIndices(const int m,const Set_t& CurrSerial)const{
   std::vector<cSetIt_t> ReturnVector(m+1);
   for(int l=0;l<=m;l++){
      ReturnVector[l]=CurrSerial.begin();
      for(int k=0;k<l;k++)++ReturnVector[l];
   }
   return ReturnVector;
}

bool CanonicalMBE::UpdateIndices(const cSetIt_t& LastIndex, const int m,
       std::vector<Set_t::iterator>& Indices)const{
   bool goodindex=false;
   cSetIt_t tempLastIndex=LastIndex;
   for(int l=m;l>=0;l--){
      ++Indices[l];
      if(Indices[l]!=tempLastIndex){//Can increment
         goodindex=true;
         for(int k=l+1;k<=m;k++){
            Indices[k]=Indices[l];
            for(int j=0;j<k-l;j++)++Indices[k];
         }
         break;
      }
      --tempLastIndex;
   }
   return goodindex;
}

}}//End namespaces


