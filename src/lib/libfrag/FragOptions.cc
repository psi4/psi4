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

#include "FragOptions.h"

namespace LibFrag{
void FragOptions::copy(const FragOptions& other){
   this->MBEOrder=other.MBEOrder;
      this->FMethod=other.FMethod;
      this->EMethod=other.EMethod;
      this->CMethod=other.CMethod;
      this->BMethod=other.BMethod;
}

void FragOptions::SetFMethod(const std::string& Frag){
   if(Frag=="user_defined"||Frag=="ud"||Frag=="user")FMethod=USER_DEFINED;
   else if(Frag=="bond_based"||Frag=="bond")FMethod=BOND_BASED;
   else if(Frag=="distance_based"||Frag=="distance")FMethod=DISTANCE_BASED;
}

void FragOptions::SetEMethod(const std::string& Embed){
   if(Embed=="none"||Embed=="no_embed"||Embed=="no"||Embed=="false")
      EMethod=NO_EMBED;
   else if(Embed=="point_charge"||Embed=="charges")EMethod=POINT_CHARGE;
   else if(Embed=="iterative"||Embed=="iterative_point_charge")
      EMethod=ITR_POINT_CHARGE;
   else if(Embed=="density")EMethod=DENSITY;
   else if(Embed=="itr_density"||Embed=="iterative_density")
      EMethod=ITR_DENSITY;
}

void FragOptions::SetCMethod(const std::string& Cap){
   if(Cap=="none"||Cap=="no_cap"||Cap=="no"||Cap=="false")CMethod=NO_CAPS;
   else if(Cap=="h_replace")CMethod=H_REPLACE;
   else if(Cap=="h_shifted")CMethod=H_SHIFTED;
}

void FragOptions::SetBMethod(const std::string& BSSE){
   if(BSSE=="none"||BSSE=="false"||BSSE=="no_bsse"||BSSE=="no")
      BMethod=NO_BSSE;
   else if(BSSE=="full"||BSSE=="bettens")BMethod=FULL;
   else if(BSSE=="MBCPN")BMethod=MBCPN;
   else if(BSSE=="VMFCN")BMethod=VMFCN;
}

void FragOptions::DefaultOptions(){
   FMethod=USER_DEFINED;
   EMethod=NO_EMBED;
   CMethod=NO_CAPS;
   BMethod=NO_BSSE;
   MBEOrder=2;
}

}


