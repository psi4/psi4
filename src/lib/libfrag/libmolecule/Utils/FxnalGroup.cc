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
   temp[HELIUM]="Elemental Helium";
   temp[LITHIUM]="Elemental Lithium";
   temp[BERYLLIUM]="Elemental Beryllium";
   temp[BORON]="Elemental Boron";
   temp[CARBON]="Elemental Carbon";
   temp[NITROGEN]="Elemental Nitrogen";
   temp[OXYGEN]="Elemental Oxygen";
   temp[FLUORINE]="Elemental Fluorine";
   temp[NEON]="Elemental Neon";
   temp[SODIUM]="Elemental Sodium";
   temp[MAGNESIUM]="Elemental Magnesium";
   temp[ALUMINIUM]="Elemental Aluminium";
   temp[SILICON]="Elemental Silicon";
   temp[PHOSPHORUS]="Elemental Phosphorus";
   temp[SULFUR]="Elemental Sulfur";
   temp[CHLORINE]="Elemental Chlorine";
   temp[ARGON]="Elemental Argon";
   temp[POTASSIUM]="Elemental Potassium";
   temp[CALCIUM]="Elemental Calcium";
   temp[SCANDIUM]="Elemental Scandium";
   temp[TITANIUM]="Elemental Titanium";
   temp[VANADIUM]="Elemental Vanadium";
   temp[CHROMIUM]="Elemental Chromium";
   temp[MANGANESE]="Elemental Manganese";
   temp[IRON]="Elemental Iron";
   temp[COBALT]="Elemental Cobalt";
   temp[NICKEL]="Elemental Nickel";
   temp[COPPER]="Elemental Copper";
   temp[ZINC]="Elemental Zinc";
   temp[GALLIUM]="Elemental Gallium";
   temp[GERMANIUM]="Elemental Germanium";
   temp[ARSENIC]="Elemental Arsenic";
   temp[SELENIUM]="Elemental Selenium";
   temp[BROMINE]="Elemental Bromine";
   temp[KRYPTON]="Elemental Krypton";
   temp[RUBIDIUM]="Elemental Rubidium";
   temp[STRONTIUM]="Elemental Strontium";
   temp[YTTRIUM]="Elemental Yttrium";
   temp[ZIRCONIUM]="Elemental Zirconium";
   temp[NIOBIUM]="Elemental Niobium";
   temp[MOLYBDENUM]="Elemental Molybdenum";
   temp[TECHNETIUM]="Elemental Technetium";
   temp[RUTHENIUM]="Elemental Ruthenium";
   temp[RHODIUM]="Elemental Rhodium";
   temp[PALLADIUM]="Elemental Palladium";
   temp[SILVER]="Elemental Silver";
   temp[CADMIUM]="Elemental Cadmium";
   temp[INDIUM]="Elemental Indium";
   temp[TIN]="Elemental Tin";
   temp[ANTIMONY]="Elemental Antimony";
   temp[TELLURIUM]="Elemental Tellurium";
   temp[IODINE]="Elemental Iodine";
   temp[XENON]="Elemental Xenon";
   temp[CESIUM]="Elemental Cesium";
   temp[BARIUM]="Elemental Barium";
   temp[LANTHANUM]="Elemental Lanthanum";
   temp[CERIUM]="Elemental Cerium";
   temp[PRASEODYMIUM]="Elemental Praseodymium";
   temp[NEODYMIUM]="Elemental Neodymium";
   temp[PROMETHIUM]="Elemental Promethium";
   temp[SAMARIUM]="Elemental Samarium";
   temp[EUROPIUM]="Elemental Europium";
   temp[GADOLINIUM]="Elemental Gadolinium";
   temp[TERBIUM]="Elemental Terbium";
   temp[DYSPROSIUM]="Elemental Dysprosium";
   temp[HOLMIUM]="Elemental Holmium";
   temp[ERBIUM]="Elemental Erbium";
   temp[THULIUM]="Elemental Thulium";
   temp[YTTERBIUM]="Elemental Ytterbium";
   temp[LUTETIUM]="Elemental Lutetium";
   temp[HAFNIUM]="Elemental Hafnium";
   temp[TANTALUM]="Elemental Tantalum";
   temp[TUNGSTEN]="Elemental Tungsten";
   temp[RHENIUM]="Elemental Rhenium";
   temp[OSMIUM]="Elemental Osmium";
   temp[IRIDIUM]="Elemental Iridium";
   temp[PLATINUM]="Elemental Platinum";
   temp[GOLD]="Elemental Gold";
   temp[MERCURY]="Elemental Mercury";
   temp[THALLIUM]="Elemental Thallium";
   temp[LEAD]="Elemental Lead";
   temp[BISMUTH]="Elemental Bismuth";
   temp[POLONIUM]="Elemental Polonium";
   temp[ASTATINE]="Elemental Astatine";
   temp[RADON]="Elemental Radon";
   temp[FRANCIUM]="Elemental Francium";
   temp[RADIUM]="Elemental Radium";
   temp[ACTINIUM]="Elemental Actinium";
   temp[THORIUM]="Elemental Thorium";
   temp[PROTACTINIUM]="Elemental Protactinium";
   temp[URANIUM]="Elemental Uranium";
   temp[NEPTUNIUM]="Elemental Neptunium";
   temp[PLUTONIUM]="Elemental Plutonium";
   temp[AMERICIUM]="Elemental Americium";
   temp[CURIUM]="Elemental Curium";
   temp[BERKELIUM]="Elemental Berkelium";
   temp[CALIFORNIUM]="Elemental Californium";
   temp[EINSTEINIUM]="Elemental Einsteinium";
   temp[FERMIUM]="Elemental Fermium";
   temp[MENDELEVIUM]="Elemental Mendelevium";
   temp[NOBELIUM]="Elemental Nobelium";
   temp[LAWRENCIUM]="Elemental Lawrencium";
   temp[RUTHERFORDIUM]="Elemental Rutherfordium";
   temp[DUBNIUM]="Elemental Dubnium";
   temp[SEABORGIUM]="Elemental Seaborgium";
   temp[BOHRIUM]="Elemental Bohrium";
   temp[HASSIUM]="Elemental Hassium";
   temp[MEITNERIUM]="Elemental Meitnerium";
   temp[DARMSTADTIUM]="Elemental Darmstadtium";
   temp[ROENTGENIUM]="Elemental Roentgenium";
   temp[COPERNICIUM]="Elemental Copernicium";
   temp[UNUNTRIUM]="Elemental Ununtrium";
   temp[FLEROVIUM]="Elemental Flerovium";
   temp[UNUNPENTIUM]="Elemental Ununpentium";
   temp[LIVERMORIUM]="Elemental Livermorium";
   temp[UNUNSEPTIUM]="Elemental Ununseptium";
   temp[UNUNOCTIUM]="Elemental Ununoctium";
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

