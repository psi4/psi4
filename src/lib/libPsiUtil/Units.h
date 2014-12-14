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
#ifndef SRC_LIB_LIBPSIUTIL_UNITS_H_
#define SRC_LIB_LIBPSIUTIL_UNITS_H_
#include "Implementations/UnitsGuts.hh"

namespace psi{
enum Units{NOUNIT,
           SECOND,//Time
           KELVIN,//Temperature
           MOLE,//Quantity
           COULOMB,//Charge
           BOHR,ANGSTROM,METER,CENTIMETER,//Distance
           DALTON,GRAM,KILOGRAM,ELECTRONMASS,//Mass
           //Note wavenumbers=cm^-1 not m^-1
           HARTREE,JOULE,CALORIE,EV,WAVENUMBERS,HZ,MHZ,CALMOL,KCALMOL//Energy
};
enum Constants{PI,TWOPI,H,C,KB,R,E0,NA,ME};
enum SIPrefixes{ATTO,FEMTO,PICO,NANO,MICRO,MILLI,CENTI,
                BASE,KILO,MEGA,GIGA,TERA,PETA};


/** \brief A class that handles SI prefix conversions
 *
 *  This is a functor that takes two arguments the SI prefix
 *  we are converting from and the one we are converting to. Typical
 *  usage would be something like:
 *
 *  \code
 *  PrefixConverter MyConverter;
 *  std::cout<<MyConverter(KILO,BASE)<<std::endl;
    //Output should be: "1000"
 *  \endcode
 *
 */
class PrefixConverter:public Converter<SIPrefixes>{
   public:
      PrefixConverter();
};

/** \brief A class that handles base unit conversions
 *
 *  This is a functor that takes two arguments the unit
 *  we are converting from and the one we are converting to. Typical
 *  usage would be something like:
 *
 *  \code
 *  BaseUnitConverter MyConverter;
 *  std::cout<<MyConverter(BOHR,ANGSTROM)<<std::endl;
    //Output should be: "0.529177" (note C++ truncates printed doubles)
 *  \endcode
 *
 */
class BaseUnitConverter: public Converter<Units>{
   public:
      BaseUnitConverter();
};

/** \brief Class for general unit conversions involving prefixes and
 *         units
 *
 *   This class is just a wrapper around the two previous classes:
 *   PrefixConverter and BaseUnitConverter, that allows a prefix and
 *   a base to be converted simultaneously.  Typical usage:
 *
 *   \code
 *   UnitConverter MyConverter;
 *   std::cout<<MyConverter(HARTREE,CALMOL,BASE,KILO)<<std::endl;
 *   //Output should be 627.51
 *   \endcode
 *
 *   Note that the prefixes come after the units.  This is because
 *   I presume that most unit conversions will be base to base, and
 *   default arguments have to be furthest to the right.
 */
class UnitConverter{
   private:
      PrefixConverter PConv_;
      BaseUnitConverter BConv_;
   public:
      double operator()(const Units& U1, const Units& U2,
             const SIPrefixes& P1=BASE, const SIPrefixes& P2=BASE)const{
         return PConv_(P1,P2)*BConv_(U1,U2);
      }
};

}//End namespace psi




#endif /* SRC_LIB_LIBPSIUTIL_UNITS_H_ */
