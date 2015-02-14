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
#include <sstream>
#include "FxnalGroup.h"
namespace psi{
namespace LibMolecule{

typedef PsiMap<FxnGroup_t,std::string> Map_t;
Map_t Init(){
   Map_t temp;
   temp[HYDROGEN]="Elemental Hydrogen";
   temp[CARBON]="Elemental Carbon";
   temp[NITROGEN]="Elemental Nitrogen";
   temp[OXYGEN]="Elemental Oxygen";
   temp[FLUORINE]="Elemental Fluorine";
   temp[CHLORINE]="Elemental Chlorine";
   temp[BROMINE]="Elemental Bromine";
   temp[IODINE]="Elemental Iodine";
   temp[METHANE]="Methane";
   temp[METHYL]="Methyl";
   temp[METHENE]="Methene";
   temp[ALKENYL1]="Primary Alkenyl";
   temp[METHYNE]="Methyne";
   temp[ALKENYL2]="Secondary Alkenyl";
   temp[ALKYNYL1]="Primary Alkynyl";
   temp[CARBON4]="Quaternary Carbon";
   temp[ALKENYL3]="Tertiary Alkenyl";
   temp[ALKYNYL2]="Secondary Alkynyl";
   temp[WATER]="Water";
   temp[HYDROXYL]="Hydroxyl";
   temp[OXYGEN2]="Secondary Oxygen";
   temp[OXYDB]="Double-Bonded Oxygen";
   temp[AMMONIA]="Ammonia";
   temp[AMINE1]="Primary Amine";
   temp[AMINE2]="Seconardy Amine";
   temp[NITRODB1]="Primary Double-Bonded Nitrogen";
   temp[AMINE3]="Tertiary Amine";
   temp[NITRODB2]="Secondary Double-Bonded Nitrogen";
   temp[NITROTB]="Triple-Bonded Nitrogen";
   temp[HYDROGENFLUORIDE]="Hydrofluoric Acid";
   temp[FLUORINE1]="Primary Fluorine";
   temp[HYDROGENCHLORIDE]="Hydrochloric Acid";
   temp[CHLORINE1]="Primary Chlorine";
   temp[HYDROGENBROMIDE]="Hydrobromic Acid";
   temp[BROMINE1]="Primary Bromine";
   temp[HYDROGENIODIDE]="Hydroiodic Acid";
   temp[IODINE1]="Primary Iodine";
   temp[NITRILE]="Nitrile";
   temp[HYDROGENCYANIDE]="Hydrogen Cyanide";
   temp[ALDEHYDE]="Aldehyde";
   temp[CARBONYL]="Carbonyl";
   temp[CARBOXYL]="Carboxyl";
   temp[FORMALDEHYDE]="Formaldehyde";
   temp[PEROXIDE]="Hydrogen Peroxide";
   temp[HYDROPEROXY]="Hydroperoxy";
   temp[METHOXY]="Methoxy";
   temp[METHANOL]="Methanol";
   temp[CCDB4]="Quaternary C-C Double Bond";
   temp[CCDB3]="Tertiary C-C Double Bond";
   temp[CCDB2]="Secondary C-C Double Bond";
   temp[ETHENYL1]="Primary Ethenyl";
   temp[ETHENYL2]="Secondardy Ethenyl";
   temp[ETHENE]="Ethene";
   temp[CCTB]="C-C Triple Bond";
   temp[ETHYNYL]="Ethynyl";
   temp[ETHYNE]="Ethyne";
   temp[KETIMINE1]="Primary Ketimine";
   temp[KETIMINE2]="Secondary Ketimine";
   temp[ALDIMINE1]="Primary Aldimine";
   temp[ALDIMINE2]="Secondary Aldimine";
   temp[METHANIMINE]="Methanimine";
   temp[AROMATICRING]="Polyaromatic Hydrocarbon";
   return temp;
}

Map_t FxnalGroup::Names_=Init();

FxnalGroup::FxnalGroup(const FxnGroup_t Type,const int Order,
      const int NMembers,const int* Members):
   Order_(Order),Type_(Type),Members_(NMembers){
   memcpy(&Members_[0],Members,sizeof(int)*NMembers);
}

std::string FxnalGroup::PrintOut(const std::string& spaces)const{
   std::stringstream result;
   result<<spaces<<Names_[Type_]<<" :";
   for(int i=0;i<Members_.size();i++)result<<Members_[i]<<" ";
   result<<std::endl;
   return result.str();
}
typedef boost::shared_ptr<const FxnalGroup> SharedGrp;

typedef std::vector<SharedGrp> GrpList;


void DerivedFxnalGrp::SetGroups(const GrpList& Groups ){
   Groups_=Groups;
   GrpList::const_iterator GrpI=Groups.begin();
      for(;GrpI!=Groups.end();++GrpI){
         int NMems=(*GrpI)->size();
         for(int i=0;i<NMems;i++)
            Members_.push_back((*(*GrpI))[i]);
      }
}


std::string DerivedFxnalGrp::PrintOut(const std::string& spaces)const{
   std::stringstream Result;
   Result<<spaces<<Names_[Type()]<<" is comprised of:"<<std::endl;
   typedef GrpList::const_iterator MyIt;
   std::string newspaces=spaces;
   newspaces+="    ";
   for(MyIt Grp=Groups_.begin();Grp!=Groups_.end();++Grp){
      Result<<(*Grp)->PrintOut(newspaces);
   }
   return Result.str();
}

AromaticRing::AromaticRing(const std::vector<int>& NAttach,
      const int Order,const Group_t& Groups):Attachment_(NAttach),
DerivedFxnalGrp(AROMATICRING,Groups.size(),Order){
   SetGroups(Groups);

}

std::string AromaticRing::PrintOut(const std::string& spaces)const{
   std::stringstream Result;
   Result<<DerivedFxnalGrp::PrintOut(spaces);
   Result<<spaces<<"Attachment Points: ";
   for(int i=0;i<NAttachPoint();i++){
      if(AttachPoint(i)==-1)continue;
      Result<<AttachPoint(i)<<" ";
   }
   Result<<std::endl;
   return Result.str();
}

}}

