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

#include "../AtomicData.h"
#include "cov_radii.h"
#include "vdw_radii.h"
#include "Z_to_element.h"
#include "masses.h"
#include "../Units.h"
namespace psi{


void AtomData::AddIsotope(const std::string& Label,
      const double Mass){
   LookUp_[Label]=Isotopes_.size();
   Isotopes_.push_back(Isotope(Mass,Label));
}

AtomicData::AtomicData(){
   //Right now this is the (reasonable) limit
   int MaxZ=LAST_COV_RADII_INDEX+1;
   BaseUnitConverter BConv;
   double AngToBohr=BConv(ANGSTROM,BOHR);
   for(int i=0,index=0;i<MaxZ;i++){
      double Mass=(i<=LAST_ATOMIC_INDEX?
            an2masses[i]:-99.0);
      double CovRad=(i<=LAST_COV_RADII_INDEX?
                  AngToBohr*cov_radii[i]:-99.0);
      double VDWRad=(i<=LAST_VDW_RADII_INDEX?
            AngToBohr*atomic_vdw_radii[i]:-99.0);
      std::string Label=(i<=LAST_ATOMIC_INDEX?
            atomic_labels[i]:"NotNamed");
      std::string Full=(i<=103?
            Z_to_element[i]:"NotNamed");
      Data_->push_back(AtomData(Label,Full,Mass,CovRad,VDWRad));
      //Here we add isotopes
      if(i==1){//Hydrogen has fun special names
         //Don't want the bare ones
         index++;
         Data_->back().AddIsotope(mass_labels[1],atomic_masses[1]);
         index++;
         //H2 and D
         Data_->back().AddIsotope(mass_labels[2],atomic_masses[2]);
         index++;
         Data_->back().LookUp_[mass_labels[3]]=
               Data_->back().LookUp_[mass_labels[2]];
         index++;
         //H3 and T
         Data_->back().AddIsotope(mass_labels[4],atomic_masses[4]);
         index++;
         Data_->back().LookUp_[mass_labels[5]]=
               Data_->back().LookUp_[mass_labels[4]];
         index++;
      }
      //Technically we aren't adding the last isotope
      else if(i>=1&&i<LAST_ATOMIC_INDEX){
         std::string StartLabel=atomic_labels[i];
         std::string EndLabel=atomic_labels[i+1];
         //We were on our normal value
         std::string label=mass_labels[++index];
         while(label!=EndLabel){
            Data_->back().AddIsotope(mass_labels[index],atomic_masses[index]);
            label=mass_labels[++index];
         }
      }
   }
}

}//End namespace
