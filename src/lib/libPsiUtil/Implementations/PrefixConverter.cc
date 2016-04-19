/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */
#include "Units.h"
#include "PsiMap.h"
#include<boost/shared_ptr.hpp>
namespace psi{
typedef PsiMap<SIPrefixes,double> Arg2;
typedef PsiMap<SIPrefixes,Arg2> Arg1;
typedef boost::shared_ptr<Arg1> SharedConversion;

void FillArray(SharedConversion& Map);

PrefixConverter::PrefixConverter(){
   if(!Conversions_)FillArray(Conversions_);
}
/********* You may want to stop reading at this point....*********/
/* Here's a bash script to generate this:
#!/bin/sh
PREFIXES=( "ATTO" "FEMTO" "PICO" "NANO" "MICRO" "MILLI" "CENTI" "BASE" "KILO" "MEGA" "GIGA" "TERA" "PETA" )
VALUES=( "-18" "-15" "-12" "-9" "-6" "-3" "-2" "0" "3" "6" "9" "12" "15" )

length=${#PREFIXES[@]}
for i in `seq 1 ${length}`;do
   indexi=$(($i-1))
   arg1=${PREFIXES[$indexi]}
   for j in `seq 1 ${length}`;do
   indexj=$(($j-1))
   arg2=${PREFIXES[$indexj]}
     if [ $i -eq $j ];then
        echo "Conversions_[$arg1][$arg2]=1.0;"
     elif [ $i -lt $j ];then
        value1=${VALUES[$indexi]}
        value2=${VALUES[$indexj]}
        diff=$((${value1}-${value2}))
        echo "Conversions_[$arg1][$arg2]=1.0E${diff};"
     else
        echo "Conversions_[$arg1][$arg2]=1/Conversions_[$arg2][$arg1];"
     fi
   done
done
 *
 */

void FillArray(SharedConversion& Map){
   Map=SharedConversion(new Arg1());
   Arg1& Conversions_=(*Map);
   Conversions_[ATTO][ATTO]=1.0;
   Conversions_[ATTO][FEMTO]=1.0E-3;
   Conversions_[ATTO][PICO]=1.0E-6;
   Conversions_[ATTO][NANO]=1.0E-9;
   Conversions_[ATTO][MICRO]=1.0E-12;
   Conversions_[ATTO][MILLI]=1.0E-15;
   Conversions_[ATTO][CENTI]=1.0E-16;
   Conversions_[ATTO][BASE]=1.0E-18;
   Conversions_[ATTO][KILO]=1.0E-21;
   Conversions_[ATTO][MEGA]=1.0E-24;
   Conversions_[ATTO][GIGA]=1.0E-27;
   Conversions_[ATTO][TERA]=1.0E-30;
   Conversions_[ATTO][PETA]=1.0E-33;
   Conversions_[FEMTO][ATTO]=1/Conversions_[ATTO][FEMTO];
   Conversions_[FEMTO][FEMTO]=1.0;
   Conversions_[FEMTO][PICO]=1.0E-3;
   Conversions_[FEMTO][NANO]=1.0E-6;
   Conversions_[FEMTO][MICRO]=1.0E-9;
   Conversions_[FEMTO][MILLI]=1.0E-12;
   Conversions_[FEMTO][CENTI]=1.0E-13;
   Conversions_[FEMTO][BASE]=1.0E-15;
   Conversions_[FEMTO][KILO]=1.0E-18;
   Conversions_[FEMTO][MEGA]=1.0E-21;
   Conversions_[FEMTO][GIGA]=1.0E-24;
   Conversions_[FEMTO][TERA]=1.0E-27;
   Conversions_[FEMTO][PETA]=1.0E-30;
   Conversions_[PICO][ATTO]=1/Conversions_[ATTO][PICO];
   Conversions_[PICO][FEMTO]=1/Conversions_[FEMTO][PICO];
   Conversions_[PICO][PICO]=1.0;
   Conversions_[PICO][NANO]=1.0E-3;
   Conversions_[PICO][MICRO]=1.0E-6;
   Conversions_[PICO][MILLI]=1.0E-9;
   Conversions_[PICO][CENTI]=1.0E-10;
   Conversions_[PICO][BASE]=1.0E-12;
   Conversions_[PICO][KILO]=1.0E-15;
   Conversions_[PICO][MEGA]=1.0E-18;
   Conversions_[PICO][GIGA]=1.0E-21;
   Conversions_[PICO][TERA]=1.0E-24;
   Conversions_[PICO][PETA]=1.0E-27;
   Conversions_[NANO][ATTO]=1/Conversions_[ATTO][NANO];
   Conversions_[NANO][FEMTO]=1/Conversions_[FEMTO][NANO];
   Conversions_[NANO][PICO]=1/Conversions_[PICO][NANO];
   Conversions_[NANO][NANO]=1.0;
   Conversions_[NANO][MICRO]=1.0E-3;
   Conversions_[NANO][MILLI]=1.0E-6;
   Conversions_[NANO][CENTI]=1.0E-7;
   Conversions_[NANO][BASE]=1.0E-9;
   Conversions_[NANO][KILO]=1.0E-12;
   Conversions_[NANO][MEGA]=1.0E-15;
   Conversions_[NANO][GIGA]=1.0E-18;
   Conversions_[NANO][TERA]=1.0E-21;
   Conversions_[NANO][PETA]=1.0E-24;
   Conversions_[MICRO][ATTO]=1/Conversions_[ATTO][MICRO];
   Conversions_[MICRO][FEMTO]=1/Conversions_[FEMTO][MICRO];
   Conversions_[MICRO][PICO]=1/Conversions_[PICO][MICRO];
   Conversions_[MICRO][NANO]=1/Conversions_[NANO][MICRO];
   Conversions_[MICRO][MICRO]=1.0;
   Conversions_[MICRO][MILLI]=1.0E-3;
   Conversions_[MICRO][CENTI]=1.0E-4;
   Conversions_[MICRO][BASE]=1.0E-6;
   Conversions_[MICRO][KILO]=1.0E-9;
   Conversions_[MICRO][MEGA]=1.0E-12;
   Conversions_[MICRO][GIGA]=1.0E-15;
   Conversions_[MICRO][TERA]=1.0E-18;
   Conversions_[MICRO][PETA]=1.0E-21;
   Conversions_[MILLI][ATTO]=1/Conversions_[ATTO][MILLI];
   Conversions_[MILLI][FEMTO]=1/Conversions_[FEMTO][MILLI];
   Conversions_[MILLI][PICO]=1/Conversions_[PICO][MILLI];
   Conversions_[MILLI][NANO]=1/Conversions_[NANO][MILLI];
   Conversions_[MILLI][MICRO]=1/Conversions_[MICRO][MILLI];
   Conversions_[MILLI][MILLI]=1.0;
   Conversions_[MILLI][CENTI]=1.0E-1;
   Conversions_[MILLI][BASE]=1.0E-3;
   Conversions_[MILLI][KILO]=1.0E-6;
   Conversions_[MILLI][MEGA]=1.0E-9;
   Conversions_[MILLI][GIGA]=1.0E-12;
   Conversions_[MILLI][TERA]=1.0E-15;
   Conversions_[MILLI][PETA]=1.0E-18;
   Conversions_[CENTI][ATTO]=1/Conversions_[ATTO][CENTI];
   Conversions_[CENTI][FEMTO]=1/Conversions_[FEMTO][CENTI];
   Conversions_[CENTI][PICO]=1/Conversions_[PICO][CENTI];
   Conversions_[CENTI][NANO]=1/Conversions_[NANO][CENTI];
   Conversions_[CENTI][MICRO]=1/Conversions_[MICRO][CENTI];
   Conversions_[CENTI][MILLI]=1/Conversions_[MILLI][CENTI];
   Conversions_[CENTI][CENTI]=1.0;
   Conversions_[CENTI][BASE]=1.0E-2;
   Conversions_[CENTI][KILO]=1.0E-5;
   Conversions_[CENTI][MEGA]=1.0E-8;
   Conversions_[CENTI][GIGA]=1.0E-11;
   Conversions_[CENTI][TERA]=1.0E-14;
   Conversions_[CENTI][PETA]=1.0E-17;
   Conversions_[BASE][ATTO]=1/Conversions_[ATTO][BASE];
   Conversions_[BASE][FEMTO]=1/Conversions_[FEMTO][BASE];
   Conversions_[BASE][PICO]=1/Conversions_[PICO][BASE];
   Conversions_[BASE][NANO]=1/Conversions_[NANO][BASE];
   Conversions_[BASE][MICRO]=1/Conversions_[MICRO][BASE];
   Conversions_[BASE][MILLI]=1/Conversions_[MILLI][BASE];
   Conversions_[BASE][CENTI]=1/Conversions_[CENTI][BASE];
   Conversions_[BASE][BASE]=1.0;
   Conversions_[BASE][KILO]=1.0E-3;
   Conversions_[BASE][MEGA]=1.0E-6;
   Conversions_[BASE][GIGA]=1.0E-9;
   Conversions_[BASE][TERA]=1.0E-12;
   Conversions_[BASE][PETA]=1.0E-15;
   Conversions_[KILO][ATTO]=1/Conversions_[ATTO][KILO];
   Conversions_[KILO][FEMTO]=1/Conversions_[FEMTO][KILO];
   Conversions_[KILO][PICO]=1/Conversions_[PICO][KILO];
   Conversions_[KILO][NANO]=1/Conversions_[NANO][KILO];
   Conversions_[KILO][MICRO]=1/Conversions_[MICRO][KILO];
   Conversions_[KILO][MILLI]=1/Conversions_[MILLI][KILO];
   Conversions_[KILO][CENTI]=1/Conversions_[CENTI][KILO];
   Conversions_[KILO][BASE]=1/Conversions_[BASE][KILO];
   Conversions_[KILO][KILO]=1.0;
   Conversions_[KILO][MEGA]=1.0E-3;
   Conversions_[KILO][GIGA]=1.0E-6;
   Conversions_[KILO][TERA]=1.0E-9;
   Conversions_[KILO][PETA]=1.0E-12;
   Conversions_[MEGA][ATTO]=1/Conversions_[ATTO][MEGA];
   Conversions_[MEGA][FEMTO]=1/Conversions_[FEMTO][MEGA];
   Conversions_[MEGA][PICO]=1/Conversions_[PICO][MEGA];
   Conversions_[MEGA][NANO]=1/Conversions_[NANO][MEGA];
   Conversions_[MEGA][MICRO]=1/Conversions_[MICRO][MEGA];
   Conversions_[MEGA][MILLI]=1/Conversions_[MILLI][MEGA];
   Conversions_[MEGA][CENTI]=1/Conversions_[CENTI][MEGA];
   Conversions_[MEGA][BASE]=1/Conversions_[BASE][MEGA];
   Conversions_[MEGA][KILO]=1/Conversions_[KILO][MEGA];
   Conversions_[MEGA][MEGA]=1.0;
   Conversions_[MEGA][GIGA]=1.0E-3;
   Conversions_[MEGA][TERA]=1.0E-6;
   Conversions_[MEGA][PETA]=1.0E-9;
   Conversions_[GIGA][ATTO]=1/Conversions_[ATTO][GIGA];
   Conversions_[GIGA][FEMTO]=1/Conversions_[FEMTO][GIGA];
   Conversions_[GIGA][PICO]=1/Conversions_[PICO][GIGA];
   Conversions_[GIGA][NANO]=1/Conversions_[NANO][GIGA];
   Conversions_[GIGA][MICRO]=1/Conversions_[MICRO][GIGA];
   Conversions_[GIGA][MILLI]=1/Conversions_[MILLI][GIGA];
   Conversions_[GIGA][CENTI]=1/Conversions_[CENTI][GIGA];
   Conversions_[GIGA][BASE]=1/Conversions_[BASE][GIGA];
   Conversions_[GIGA][KILO]=1/Conversions_[KILO][GIGA];
   Conversions_[GIGA][MEGA]=1/Conversions_[MEGA][GIGA];
   Conversions_[GIGA][GIGA]=1.0;
   Conversions_[GIGA][TERA]=1.0E-3;
   Conversions_[GIGA][PETA]=1.0E-6;
   Conversions_[TERA][ATTO]=1/Conversions_[ATTO][TERA];
   Conversions_[TERA][FEMTO]=1/Conversions_[FEMTO][TERA];
   Conversions_[TERA][PICO]=1/Conversions_[PICO][TERA];
   Conversions_[TERA][NANO]=1/Conversions_[NANO][TERA];
   Conversions_[TERA][MICRO]=1/Conversions_[MICRO][TERA];
   Conversions_[TERA][MILLI]=1/Conversions_[MILLI][TERA];
   Conversions_[TERA][CENTI]=1/Conversions_[CENTI][TERA];
   Conversions_[TERA][BASE]=1/Conversions_[BASE][TERA];
   Conversions_[TERA][KILO]=1/Conversions_[KILO][TERA];
   Conversions_[TERA][MEGA]=1/Conversions_[MEGA][TERA];
   Conversions_[TERA][GIGA]=1/Conversions_[GIGA][TERA];
   Conversions_[TERA][TERA]=1.0;
   Conversions_[TERA][PETA]=1.0E-3;
   Conversions_[PETA][ATTO]=1/Conversions_[ATTO][PETA];
   Conversions_[PETA][FEMTO]=1/Conversions_[FEMTO][PETA];
   Conversions_[PETA][PICO]=1/Conversions_[PICO][PETA];
   Conversions_[PETA][NANO]=1/Conversions_[NANO][PETA];
   Conversions_[PETA][MICRO]=1/Conversions_[MICRO][PETA];
   Conversions_[PETA][MILLI]=1/Conversions_[MILLI][PETA];
   Conversions_[PETA][CENTI]=1/Conversions_[CENTI][PETA];
   Conversions_[PETA][BASE]=1/Conversions_[BASE][PETA];
   Conversions_[PETA][KILO]=1/Conversions_[KILO][PETA];
   Conversions_[PETA][MEGA]=1/Conversions_[MEGA][PETA];
   Conversions_[PETA][GIGA]=1/Conversions_[GIGA][PETA];
   Conversions_[PETA][TERA]=1/Conversions_[TERA][PETA];
   Conversions_[PETA][PETA]=1.0;
}

}//End namespace
