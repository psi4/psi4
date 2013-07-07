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

/*! 
  \file Z_to_element.h
  \ingroup
  \brief convert element name/symbol to atomic number
*/

#ifndef _psi_include_Z_to_element_h_
#define _psi_include_Z_to_element_h_

#include <string>

namespace psi {

const std::string Z_to_element [] = {
//0-2
"GHOST", "HYDROGEN", "HELIUM",

//3-10
"LITHIUM", "BERYLLIUM", "BORON", "CARBON", "NITROGEN", "OXYGEN", "FLUORINE", "NEON",

//11-36
"SODIUM", "MAGNESIUM", "ALUMINUM", "SILICON", "PHOSPHORUS", "SULFUR", "CHLORINE", "ARGON",
"POTASSIUM", "CALCIUM", "SCANDIUM", "TITANIUM", "VANADIUM", "CHROMIUM", "MANGANESE",
"IRON", "COBALT", "NICKEL", "COPPER", "ZINC", "GALLIUM", "GERMANIUM", "ARSENIC", "SELENIUM",
"BROMINE", "KRYPTON",

//37-54
"RUBIDIUM", "STRONTIUM", "YTTRIUM", "ZIRCONIUM", "NIOBIUM", "MOLYBDENUM", "TECHNETIUM",
"RUTHENIUM", "RHODIUM", "PALLADIUM", "SILVER", "CADMIUM", "INDIUM", "TIN", "ANTIMONY",
"TELLURIUM", "IODINE", "XENON",

"CESIUM", "BARIUM", "LANTHANUM", "CERIUM", "PRASEODYMIUM", "NEODYMIUM", "PROMETHIUM",
"SAMARIUM", "EUROPIUM", "GADOLINIUM", "TERBIUM", "DYSPROSIUM", "HOLMIUM", "ERBIUM",
"THULIUM", "YTTERBIUM", "LUTETIUM", "HAFNIUM", "TANTALUM", "TUNGSTEN", "RHENIUM",
"OSMIUM", "IRIDIUM", "PLATINUM", "GOLD", "MERCURY", "THALLIUM", "LEAD", "BISMUTH",
"POLONIUM", "ASTATINE", "RADON",

"FRANCIUM", "RADIUM", "ACTINIUM", "THORIUM", "PROTACTINIUM", "URANIUM", "NEPTUNIUM",
"PLUTONIUM", "AMERICIUM", "CURIUM", "BERKELIUM", "CALIFORNIUM", "EINSTEINIUM", "FERMIUM",
"MENDELEVIUM", "NOBELIUM", "LAWRENCIUM", "UNQ", "UNP", "UNH", "UNS",

};

}

#endif

