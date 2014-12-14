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
#include "physconst.h"
#include "Units.h"
#include "PsiMap.h"
#include<boost/shared_ptr.hpp>
namespace psi{
typedef PsiMap<Units,double> Arg2;
typedef PsiMap<Units,Arg2> Arg1;
typedef boost::shared_ptr<Arg1> SharedConversion;




void FillMap(SharedConversion& Map);

BaseUnitConverter::BaseUnitConverter(){
   if(!Conversions_)FillMap(Conversions_);
}

/************May want to stop reading here**************/
/* Again a bash script:
#!/bin/sh
#units=( "BOHR" "ANGSTROM" "METER" )
#conv=( "pc_bohr2angstroms" "pc_bohr2m" )
#size=3
#units=( "DALTON" "GRAM" "ELECTRONMASS" )
#conv=( "pc_amu2g" "1.0/pc_au2amu" )
#size=3
units=( "HARTREE" "JOULE" "CALORIE" "EV" "WAVENUMBERS" "HZ" "CALMOL" "KCALMOL" )
conv=( "pc_hartree2J" "pc_hartree2J/pc_cal2J" "pc_hartree2ev" "pc_hartree2wavenumbers" "pc_hartree2MHz*1E6" "pc_hartree2kcalmol*1E3" "pc_hartree2kcalmol")
size=7
for ((i=0;i<${size};i++)) do
   arg1=${units[$i]}
   for ((k=0;k<${size};k++)) do
       arg2=${units[${k}]}
       value="1/Conversions_[$arg2][$arg1]"
        if [ $i -eq $k ];then
            value="1.0"
        elif [ $i -eq 0 ];then
            indexk=$(($k-1))
            value="${conv[$indexk]}"
        elif  [ $i -lt $k ];then
            label="${units[0]}"
            value="Conversions_[$arg1][${label}]*Conversions_[${label}][$arg2]"
        fi
        echo "Conversions_[$arg1][$arg2]=${value};"
   done
done
*/
void FillMap(SharedConversion& Map){
   Map=SharedConversion(new Arg1());
   Arg1& Conversions_=(*Map);
   Conversions_[BOHR][BOHR]=1.0;
   Conversions_[BOHR][ANGSTROM]=pc_bohr2angstroms;
   Conversions_[BOHR][METER]=pc_bohr2m;
   Conversions_[ANGSTROM][BOHR]=1/Conversions_[BOHR][ANGSTROM];
   Conversions_[ANGSTROM][ANGSTROM]=1.0;
   Conversions_[ANGSTROM][METER]=Conversions_[ANGSTROM][BOHR]*Conversions_[BOHR][METER];
   Conversions_[METER][BOHR]=1/Conversions_[BOHR][METER];
   Conversions_[METER][ANGSTROM]=1/Conversions_[ANGSTROM][METER];
   Conversions_[METER][METER]=1.0;
   Conversions_[DALTON][DALTON]=1.0;
   Conversions_[DALTON][GRAM]=pc_amu2g;
   Conversions_[DALTON][ELECTRONMASS]=1.0/pc_au2amu;
   Conversions_[GRAM][DALTON]=1/Conversions_[DALTON][GRAM];
   Conversions_[GRAM][GRAM]=1.0;
   Conversions_[GRAM][ELECTRONMASS]=Conversions_[GRAM][DALTON]*Conversions_[DALTON][ELECTRONMASS];
   Conversions_[ELECTRONMASS][DALTON]=1/Conversions_[DALTON][ELECTRONMASS];
   Conversions_[ELECTRONMASS][GRAM]=1/Conversions_[GRAM][ELECTRONMASS];
   Conversions_[ELECTRONMASS][ELECTRONMASS]=1.0;
   Conversions_[HARTREE][HARTREE]=1.0;
   Conversions_[HARTREE][JOULE]=pc_hartree2J;
   Conversions_[HARTREE][CALORIE]=pc_hartree2J/pc_cal2J;
   Conversions_[HARTREE][EV]=pc_hartree2ev;
   Conversions_[HARTREE][WAVENUMBERS]=pc_hartree2wavenumbers;
   Conversions_[HARTREE][HZ]=pc_hartree2MHz*1E6;
   Conversions_[HARTREE][CALMOL]=pc_hartree2kcalmol*1E3;
   Conversions_[HARTREE][KCALMOL]=pc_hartree2kcalmol;
   Conversions_[JOULE][HARTREE]=1/Conversions_[HARTREE][JOULE];
   Conversions_[JOULE][JOULE]=1.0;
   Conversions_[JOULE][CALORIE]=Conversions_[JOULE][HARTREE]*Conversions_[HARTREE][CALORIE];
   Conversions_[JOULE][EV]=Conversions_[JOULE][HARTREE]*Conversions_[HARTREE][EV];
   Conversions_[JOULE][WAVENUMBERS]=Conversions_[JOULE][HARTREE]*Conversions_[HARTREE][WAVENUMBERS];
   Conversions_[JOULE][HZ]=Conversions_[JOULE][HARTREE]*Conversions_[HARTREE][HZ];
   Conversions_[JOULE][CALMOL]=Conversions_[JOULE][HARTREE]*Conversions_[HARTREE][CALMOL];
   Conversions_[JOULE][KCALMOL]=Conversions_[JOULE][HARTREE]*Conversions_[HARTREE][KCALMOL];
   Conversions_[CALORIE][HARTREE]=1/Conversions_[HARTREE][CALORIE];
   Conversions_[CALORIE][JOULE]=1/Conversions_[JOULE][CALORIE];
   Conversions_[CALORIE][CALORIE]=1.0;
   Conversions_[CALORIE][EV]=Conversions_[CALORIE][HARTREE]*Conversions_[HARTREE][EV];
   Conversions_[CALORIE][WAVENUMBERS]=Conversions_[CALORIE][HARTREE]*Conversions_[HARTREE][WAVENUMBERS];
   Conversions_[CALORIE][HZ]=Conversions_[CALORIE][HARTREE]*Conversions_[HARTREE][HZ];
   Conversions_[CALORIE][CALMOL]=Conversions_[CALORIE][HARTREE]*Conversions_[HARTREE][CALMOL];
   Conversions_[CALORIE][KCALMOL]=Conversions_[CALORIE][HARTREE]*Conversions_[HARTREE][KCALMOL];
   Conversions_[EV][HARTREE]=1/Conversions_[HARTREE][EV];
   Conversions_[EV][JOULE]=1/Conversions_[JOULE][EV];
   Conversions_[EV][CALORIE]=1/Conversions_[CALORIE][EV];
   Conversions_[EV][EV]=1.0;
   Conversions_[EV][WAVENUMBERS]=Conversions_[EV][HARTREE]*Conversions_[HARTREE][WAVENUMBERS];
   Conversions_[EV][HZ]=Conversions_[EV][HARTREE]*Conversions_[HARTREE][HZ];
   Conversions_[EV][CALMOL]=Conversions_[EV][HARTREE]*Conversions_[HARTREE][CALMOL];
   Conversions_[EV][KCALMOL]=Conversions_[EV][HARTREE]*Conversions_[HARTREE][KCALMOL];
   Conversions_[WAVENUMBERS][HARTREE]=1/Conversions_[HARTREE][WAVENUMBERS];
   Conversions_[WAVENUMBERS][JOULE]=1/Conversions_[JOULE][WAVENUMBERS];
   Conversions_[WAVENUMBERS][CALORIE]=1/Conversions_[CALORIE][WAVENUMBERS];
   Conversions_[WAVENUMBERS][EV]=1/Conversions_[EV][WAVENUMBERS];
   Conversions_[WAVENUMBERS][WAVENUMBERS]=1.0;
   Conversions_[WAVENUMBERS][HZ]=Conversions_[WAVENUMBERS][HARTREE]*Conversions_[HARTREE][HZ];
   Conversions_[WAVENUMBERS][CALMOL]=Conversions_[WAVENUMBERS][HARTREE]*Conversions_[HARTREE][CALMOL];
   Conversions_[WAVENUMBERS][KCALMOL]=Conversions_[WAVENUMBERS][HARTREE]*Conversions_[HARTREE][KCALMOL];
   Conversions_[HZ][HARTREE]=1/Conversions_[HARTREE][HZ];
   Conversions_[HZ][JOULE]=1/Conversions_[JOULE][HZ];
   Conversions_[HZ][CALORIE]=1/Conversions_[CALORIE][HZ];
   Conversions_[HZ][EV]=1/Conversions_[EV][HZ];
   Conversions_[HZ][WAVENUMBERS]=1/Conversions_[WAVENUMBERS][HZ];
   Conversions_[HZ][HZ]=1.0;
   Conversions_[HZ][CALMOL]=Conversions_[HZ][HARTREE]*Conversions_[HARTREE][CALMOL];
   Conversions_[HZ][KCALMOL]=Conversions_[HZ][HARTREE]*Conversions_[HARTREE][KCALMOL];
   Conversions_[CALMOL][HARTREE]=1/Conversions_[HARTREE][CALMOL];
   Conversions_[CALMOL][JOULE]=1/Conversions_[JOULE][CALMOL];
   Conversions_[CALMOL][CALORIE]=1/Conversions_[CALORIE][CALMOL];
   Conversions_[CALMOL][EV]=1/Conversions_[EV][CALMOL];
   Conversions_[CALMOL][WAVENUMBERS]=1/Conversions_[WAVENUMBERS][CALMOL];
   Conversions_[CALMOL][HZ]=1/Conversions_[HZ][CALMOL];
   Conversions_[CALMOL][CALMOL]=1.0;
   Conversions_[CALMOL][KCALMOL]=Conversions_[CALMOL][HARTREE]*Conversions_[HARTREE][KCALMOL];
   Conversions_[KCALMOL][HARTREE]=1/Conversions_[HARTREE][KCALMOL];
   Conversions_[KCALMOL][JOULE]=1/Conversions_[JOULE][KCALMOL];
   Conversions_[KCALMOL][CALORIE]=1/Conversions_[CALORIE][KCALMOL];
   Conversions_[KCALMOL][EV]=1/Conversions_[EV][KCALMOL];
   Conversions_[KCALMOL][WAVENUMBERS]=1/Conversions_[WAVENUMBERS][KCALMOL];
   Conversions_[KCALMOL][HZ]=1/Conversions_[HZ][KCALMOL];
   Conversions_[KCALMOL][CALMOL]=1/Conversions_[CALMOL][KCALMOL];
   Conversions_[KCALMOL][KCALMOL]=1.0;

}

}//End namespace psi

