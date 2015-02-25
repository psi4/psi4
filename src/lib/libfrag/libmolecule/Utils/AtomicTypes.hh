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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ATOMICTYPES_HH_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ATOMICTYPES_HH_


namespace psi{
namespace LibMolecule{

class Hydrogen:public FxnalGroup{
public:
Hydrogen():FxnalGroup(HYDROGEN,0){}
Hydrogen(const int* Members):FxnalGroup(HYDROGEN,0,1,Members){}
};
class Helium:public FxnalGroup{
public:
Helium():FxnalGroup(HELIUM,0){}
Helium(const int* Members):FxnalGroup(HELIUM,0,1,Members){}
};
class Lithium:public FxnalGroup{
public:
Lithium():FxnalGroup(LITHIUM,0){}
Lithium(const int* Members):FxnalGroup(LITHIUM,0,1,Members){}
};
class Beryllium:public FxnalGroup{
public:
Beryllium():FxnalGroup(BERYLLIUM,0){}
Beryllium(const int* Members):FxnalGroup(BERYLLIUM,0,1,Members){}
};
class Boron:public FxnalGroup{
public:
Boron():FxnalGroup(BORON,0){}
Boron(const int* Members):FxnalGroup(BORON,0,1,Members){}
};
class Carbon:public FxnalGroup{
public:
Carbon():FxnalGroup(CARBON,0){}
Carbon(const int* Members):FxnalGroup(CARBON,0,1,Members){}
};
class Nitrogen:public FxnalGroup{
public:
Nitrogen():FxnalGroup(NITROGEN,0){}
Nitrogen(const int* Members):FxnalGroup(NITROGEN,0,1,Members){}
};
class Oxygen:public FxnalGroup{
public:
Oxygen():FxnalGroup(OXYGEN,0){}
Oxygen(const int* Members):FxnalGroup(OXYGEN,0,1,Members){}
};
class Fluorine:public FxnalGroup{
public:
Fluorine():FxnalGroup(FLUORINE,0){}
Fluorine(const int* Members):FxnalGroup(FLUORINE,0,1,Members){}
};
class Neon:public FxnalGroup{
public:
Neon():FxnalGroup(NEON,0){}
Neon(const int* Members):FxnalGroup(NEON,0,1,Members){}
};
class Sodium:public FxnalGroup{
public:
Sodium():FxnalGroup(SODIUM,0){}
Sodium(const int* Members):FxnalGroup(SODIUM,0,1,Members){}
};
class Magnesium:public FxnalGroup{
public:
Magnesium():FxnalGroup(MAGNESIUM,0){}
Magnesium(const int* Members):FxnalGroup(MAGNESIUM,0,1,Members){}
};
class Aluminium:public FxnalGroup{
public:
Aluminium():FxnalGroup(ALUMINIUM,0){}
Aluminium(const int* Members):FxnalGroup(ALUMINIUM,0,1,Members){}
};
class Silicon:public FxnalGroup{
public:
Silicon():FxnalGroup(SILICON,0){}
Silicon(const int* Members):FxnalGroup(SILICON,0,1,Members){}
};
class Phosphorus:public FxnalGroup{
public:
Phosphorus():FxnalGroup(PHOSPHORUS,0){}
Phosphorus(const int* Members):FxnalGroup(PHOSPHORUS,0,1,Members){}
};
class Sulfur:public FxnalGroup{
public:
Sulfur():FxnalGroup(SULFUR,0){}
Sulfur(const int* Members):FxnalGroup(SULFUR,0,1,Members){}
};
class Chlorine:public FxnalGroup{
public:
Chlorine():FxnalGroup(CHLORINE,0){}
Chlorine(const int* Members):FxnalGroup(CHLORINE,0,1,Members){}
};
class Argon:public FxnalGroup{
public:
Argon():FxnalGroup(ARGON,0){}
Argon(const int* Members):FxnalGroup(ARGON,0,1,Members){}
};
class Potassium:public FxnalGroup{
public:
Potassium():FxnalGroup(POTASSIUM,0){}
Potassium(const int* Members):FxnalGroup(POTASSIUM,0,1,Members){}
};
class Calcium:public FxnalGroup{
public:
Calcium():FxnalGroup(CALCIUM,0){}
Calcium(const int* Members):FxnalGroup(CALCIUM,0,1,Members){}
};
class Scandium:public FxnalGroup{
public:
Scandium():FxnalGroup(SCANDIUM,0){}
Scandium(const int* Members):FxnalGroup(SCANDIUM,0,1,Members){}
};
class Titanium:public FxnalGroup{
public:
Titanium():FxnalGroup(TITANIUM,0){}
Titanium(const int* Members):FxnalGroup(TITANIUM,0,1,Members){}
};
class Vanadium:public FxnalGroup{
public:
Vanadium():FxnalGroup(VANADIUM,0){}
Vanadium(const int* Members):FxnalGroup(VANADIUM,0,1,Members){}
};
class Chromium:public FxnalGroup{
public:
Chromium():FxnalGroup(CHROMIUM,0){}
Chromium(const int* Members):FxnalGroup(CHROMIUM,0,1,Members){}
};
class Manganese:public FxnalGroup{
public:
Manganese():FxnalGroup(MANGANESE,0){}
Manganese(const int* Members):FxnalGroup(MANGANESE,0,1,Members){}
};
class Iron:public FxnalGroup{
public:
Iron():FxnalGroup(IRON,0){}
Iron(const int* Members):FxnalGroup(IRON,0,1,Members){}
};
class Cobalt:public FxnalGroup{
public:
Cobalt():FxnalGroup(COBALT,0){}
Cobalt(const int* Members):FxnalGroup(COBALT,0,1,Members){}
};
class Nickel:public FxnalGroup{
public:
Nickel():FxnalGroup(NICKEL,0){}
Nickel(const int* Members):FxnalGroup(NICKEL,0,1,Members){}
};
class Copper:public FxnalGroup{
public:
Copper():FxnalGroup(COPPER,0){}
Copper(const int* Members):FxnalGroup(COPPER,0,1,Members){}
};
class Zinc:public FxnalGroup{
public:
Zinc():FxnalGroup(ZINC,0){}
Zinc(const int* Members):FxnalGroup(ZINC,0,1,Members){}
};
class Gallium:public FxnalGroup{
public:
Gallium():FxnalGroup(GALLIUM,0){}
Gallium(const int* Members):FxnalGroup(GALLIUM,0,1,Members){}
};
class Germanium:public FxnalGroup{
public:
Germanium():FxnalGroup(GERMANIUM,0){}
Germanium(const int* Members):FxnalGroup(GERMANIUM,0,1,Members){}
};
class Arsenic:public FxnalGroup{
public:
Arsenic():FxnalGroup(ARSENIC,0){}
Arsenic(const int* Members):FxnalGroup(ARSENIC,0,1,Members){}
};
class Selenium:public FxnalGroup{
public:
Selenium():FxnalGroup(SELENIUM,0){}
Selenium(const int* Members):FxnalGroup(SELENIUM,0,1,Members){}
};
class Bromine:public FxnalGroup{
public:
Bromine():FxnalGroup(BROMINE,0){}
Bromine(const int* Members):FxnalGroup(BROMINE,0,1,Members){}
};
class Krypton:public FxnalGroup{
public:
Krypton():FxnalGroup(KRYPTON,0){}
Krypton(const int* Members):FxnalGroup(KRYPTON,0,1,Members){}
};
class Rubidium:public FxnalGroup{
public:
Rubidium():FxnalGroup(RUBIDIUM,0){}
Rubidium(const int* Members):FxnalGroup(RUBIDIUM,0,1,Members){}
};
class Strontium:public FxnalGroup{
public:
Strontium():FxnalGroup(STRONTIUM,0){}
Strontium(const int* Members):FxnalGroup(STRONTIUM,0,1,Members){}
};
class Yttrium:public FxnalGroup{
public:
Yttrium():FxnalGroup(YTTRIUM,0){}
Yttrium(const int* Members):FxnalGroup(YTTRIUM,0,1,Members){}
};
class Zirconium:public FxnalGroup{
public:
Zirconium():FxnalGroup(ZIRCONIUM,0){}
Zirconium(const int* Members):FxnalGroup(ZIRCONIUM,0,1,Members){}
};
class Niobium:public FxnalGroup{
public:
Niobium():FxnalGroup(NIOBIUM,0){}
Niobium(const int* Members):FxnalGroup(NIOBIUM,0,1,Members){}
};
class Molybdenum:public FxnalGroup{
public:
Molybdenum():FxnalGroup(MOLYBDENUM,0){}
Molybdenum(const int* Members):FxnalGroup(MOLYBDENUM,0,1,Members){}
};
class Technetium:public FxnalGroup{
public:
Technetium():FxnalGroup(TECHNETIUM,0){}
Technetium(const int* Members):FxnalGroup(TECHNETIUM,0,1,Members){}
};
class Ruthenium:public FxnalGroup{
public:
Ruthenium():FxnalGroup(RUTHENIUM,0){}
Ruthenium(const int* Members):FxnalGroup(RUTHENIUM,0,1,Members){}
};
class Rhodium:public FxnalGroup{
public:
Rhodium():FxnalGroup(RHODIUM,0){}
Rhodium(const int* Members):FxnalGroup(RHODIUM,0,1,Members){}
};
class Palladium:public FxnalGroup{
public:
Palladium():FxnalGroup(PALLADIUM,0){}
Palladium(const int* Members):FxnalGroup(PALLADIUM,0,1,Members){}
};
class Silver:public FxnalGroup{
public:
Silver():FxnalGroup(SILVER,0){}
Silver(const int* Members):FxnalGroup(SILVER,0,1,Members){}
};
class Cadmium:public FxnalGroup{
public:
Cadmium():FxnalGroup(CADMIUM,0){}
Cadmium(const int* Members):FxnalGroup(CADMIUM,0,1,Members){}
};
class Indium:public FxnalGroup{
public:
Indium():FxnalGroup(INDIUM,0){}
Indium(const int* Members):FxnalGroup(INDIUM,0,1,Members){}
};
class Tin:public FxnalGroup{
public:
Tin():FxnalGroup(TIN,0){}
Tin(const int* Members):FxnalGroup(TIN,0,1,Members){}
};
class Antimony:public FxnalGroup{
public:
Antimony():FxnalGroup(ANTIMONY,0){}
Antimony(const int* Members):FxnalGroup(ANTIMONY,0,1,Members){}
};
class Tellurium:public FxnalGroup{
public:
Tellurium():FxnalGroup(TELLURIUM,0){}
Tellurium(const int* Members):FxnalGroup(TELLURIUM,0,1,Members){}
};
class Iodine:public FxnalGroup{
public:
Iodine():FxnalGroup(IODINE,0){}
Iodine(const int* Members):FxnalGroup(IODINE,0,1,Members){}
};
class Xenon:public FxnalGroup{
public:
Xenon():FxnalGroup(XENON,0){}
Xenon(const int* Members):FxnalGroup(XENON,0,1,Members){}
};
class Cesium:public FxnalGroup{
public:
Cesium():FxnalGroup(CESIUM,0){}
Cesium(const int* Members):FxnalGroup(CESIUM,0,1,Members){}
};
class Barium:public FxnalGroup{
public:
Barium():FxnalGroup(BARIUM,0){}
Barium(const int* Members):FxnalGroup(BARIUM,0,1,Members){}
};
class Lanthanum:public FxnalGroup{
public:
Lanthanum():FxnalGroup(LANTHANUM,0){}
Lanthanum(const int* Members):FxnalGroup(LANTHANUM,0,1,Members){}
};
class Cerium:public FxnalGroup{
public:
Cerium():FxnalGroup(CERIUM,0){}
Cerium(const int* Members):FxnalGroup(CERIUM,0,1,Members){}
};
class Praseodymium:public FxnalGroup{
public:
Praseodymium():FxnalGroup(PRASEODYMIUM,0){}
Praseodymium(const int* Members):FxnalGroup(PRASEODYMIUM,0,1,Members){}
};
class Neodymium:public FxnalGroup{
public:
Neodymium():FxnalGroup(NEODYMIUM,0){}
Neodymium(const int* Members):FxnalGroup(NEODYMIUM,0,1,Members){}
};
class Promethium:public FxnalGroup{
public:
Promethium():FxnalGroup(PROMETHIUM,0){}
Promethium(const int* Members):FxnalGroup(PROMETHIUM,0,1,Members){}
};
class Samarium:public FxnalGroup{
public:
Samarium():FxnalGroup(SAMARIUM,0){}
Samarium(const int* Members):FxnalGroup(SAMARIUM,0,1,Members){}
};
class Europium:public FxnalGroup{
public:
Europium():FxnalGroup(EUROPIUM,0){}
Europium(const int* Members):FxnalGroup(EUROPIUM,0,1,Members){}
};
class Gadolinium:public FxnalGroup{
public:
Gadolinium():FxnalGroup(GADOLINIUM,0){}
Gadolinium(const int* Members):FxnalGroup(GADOLINIUM,0,1,Members){}
};
class Terbium:public FxnalGroup{
public:
Terbium():FxnalGroup(TERBIUM,0){}
Terbium(const int* Members):FxnalGroup(TERBIUM,0,1,Members){}
};
class Dysprosium:public FxnalGroup{
public:
Dysprosium():FxnalGroup(DYSPROSIUM,0){}
Dysprosium(const int* Members):FxnalGroup(DYSPROSIUM,0,1,Members){}
};
class Holmium:public FxnalGroup{
public:
Holmium():FxnalGroup(HOLMIUM,0){}
Holmium(const int* Members):FxnalGroup(HOLMIUM,0,1,Members){}
};
class Erbium:public FxnalGroup{
public:
Erbium():FxnalGroup(ERBIUM,0){}
Erbium(const int* Members):FxnalGroup(ERBIUM,0,1,Members){}
};
class Thulium:public FxnalGroup{
public:
Thulium():FxnalGroup(THULIUM,0){}
Thulium(const int* Members):FxnalGroup(THULIUM,0,1,Members){}
};
class Ytterbium:public FxnalGroup{
public:
Ytterbium():FxnalGroup(YTTERBIUM,0){}
Ytterbium(const int* Members):FxnalGroup(YTTERBIUM,0,1,Members){}
};
class Lutetium:public FxnalGroup{
public:
Lutetium():FxnalGroup(LUTETIUM,0){}
Lutetium(const int* Members):FxnalGroup(LUTETIUM,0,1,Members){}
};
class Hafnium:public FxnalGroup{
public:
Hafnium():FxnalGroup(HAFNIUM,0){}
Hafnium(const int* Members):FxnalGroup(HAFNIUM,0,1,Members){}
};
class Tantalum:public FxnalGroup{
public:
Tantalum():FxnalGroup(TANTALUM,0){}
Tantalum(const int* Members):FxnalGroup(TANTALUM,0,1,Members){}
};
class Tungsten:public FxnalGroup{
public:
Tungsten():FxnalGroup(TUNGSTEN,0){}
Tungsten(const int* Members):FxnalGroup(TUNGSTEN,0,1,Members){}
};
class Rhenium:public FxnalGroup{
public:
Rhenium():FxnalGroup(RHENIUM,0){}
Rhenium(const int* Members):FxnalGroup(RHENIUM,0,1,Members){}
};
class Osmium:public FxnalGroup{
public:
Osmium():FxnalGroup(OSMIUM,0){}
Osmium(const int* Members):FxnalGroup(OSMIUM,0,1,Members){}
};
class Iridium:public FxnalGroup{
public:
Iridium():FxnalGroup(IRIDIUM,0){}
Iridium(const int* Members):FxnalGroup(IRIDIUM,0,1,Members){}
};
class Platinum:public FxnalGroup{
public:
Platinum():FxnalGroup(PLATINUM,0){}
Platinum(const int* Members):FxnalGroup(PLATINUM,0,1,Members){}
};
class Gold:public FxnalGroup{
public:
Gold():FxnalGroup(GOLD,0){}
Gold(const int* Members):FxnalGroup(GOLD,0,1,Members){}
};
class Mercury:public FxnalGroup{
public:
Mercury():FxnalGroup(MERCURY,0){}
Mercury(const int* Members):FxnalGroup(MERCURY,0,1,Members){}
};
class Thallium:public FxnalGroup{
public:
Thallium():FxnalGroup(THALLIUM,0){}
Thallium(const int* Members):FxnalGroup(THALLIUM,0,1,Members){}
};
class Lead:public FxnalGroup{
public:
Lead():FxnalGroup(LEAD,0){}
Lead(const int* Members):FxnalGroup(LEAD,0,1,Members){}
};
class Bismuth:public FxnalGroup{
public:
Bismuth():FxnalGroup(BISMUTH,0){}
Bismuth(const int* Members):FxnalGroup(BISMUTH,0,1,Members){}
};
class Polonium:public FxnalGroup{
public:
Polonium():FxnalGroup(POLONIUM,0){}
Polonium(const int* Members):FxnalGroup(POLONIUM,0,1,Members){}
};
class Astatine:public FxnalGroup{
public:
Astatine():FxnalGroup(ASTATINE,0){}
Astatine(const int* Members):FxnalGroup(ASTATINE,0,1,Members){}
};
class Radon:public FxnalGroup{
public:
Radon():FxnalGroup(RADON,0){}
Radon(const int* Members):FxnalGroup(RADON,0,1,Members){}
};
class Francium:public FxnalGroup{
public:
Francium():FxnalGroup(FRANCIUM,0){}
Francium(const int* Members):FxnalGroup(FRANCIUM,0,1,Members){}
};
class Radium:public FxnalGroup{
public:
Radium():FxnalGroup(RADIUM,0){}
Radium(const int* Members):FxnalGroup(RADIUM,0,1,Members){}
};
class Actinium:public FxnalGroup{
public:
Actinium():FxnalGroup(ACTINIUM,0){}
Actinium(const int* Members):FxnalGroup(ACTINIUM,0,1,Members){}
};
class Thorium:public FxnalGroup{
public:
Thorium():FxnalGroup(THORIUM,0){}
Thorium(const int* Members):FxnalGroup(THORIUM,0,1,Members){}
};
class Protactinium:public FxnalGroup{
public:
Protactinium():FxnalGroup(PROTACTINIUM,0){}
Protactinium(const int* Members):FxnalGroup(PROTACTINIUM,0,1,Members){}
};
class Uranium:public FxnalGroup{
public:
Uranium():FxnalGroup(URANIUM,0){}
Uranium(const int* Members):FxnalGroup(URANIUM,0,1,Members){}
};
class Neptunium:public FxnalGroup{
public:
Neptunium():FxnalGroup(NEPTUNIUM,0){}
Neptunium(const int* Members):FxnalGroup(NEPTUNIUM,0,1,Members){}
};
class Plutonium:public FxnalGroup{
public:
Plutonium():FxnalGroup(PLUTONIUM,0){}
Plutonium(const int* Members):FxnalGroup(PLUTONIUM,0,1,Members){}
};
class Americium:public FxnalGroup{
public:
Americium():FxnalGroup(AMERICIUM,0){}
Americium(const int* Members):FxnalGroup(AMERICIUM,0,1,Members){}
};
class Curium:public FxnalGroup{
public:
Curium():FxnalGroup(CURIUM,0){}
Curium(const int* Members):FxnalGroup(CURIUM,0,1,Members){}
};
class Berkelium:public FxnalGroup{
public:
Berkelium():FxnalGroup(BERKELIUM,0){}
Berkelium(const int* Members):FxnalGroup(BERKELIUM,0,1,Members){}
};
class Californium:public FxnalGroup{
public:
Californium():FxnalGroup(CALIFORNIUM,0){}
Californium(const int* Members):FxnalGroup(CALIFORNIUM,0,1,Members){}
};
class Einsteinium:public FxnalGroup{
public:
Einsteinium():FxnalGroup(EINSTEINIUM,0){}
Einsteinium(const int* Members):FxnalGroup(EINSTEINIUM,0,1,Members){}
};
class Fermium:public FxnalGroup{
public:
Fermium():FxnalGroup(FERMIUM,0){}
Fermium(const int* Members):FxnalGroup(FERMIUM,0,1,Members){}
};
class Mendelevium:public FxnalGroup{
public:
Mendelevium():FxnalGroup(MENDELEVIUM,0){}
Mendelevium(const int* Members):FxnalGroup(MENDELEVIUM,0,1,Members){}
};
class Nobelium:public FxnalGroup{
public:
Nobelium():FxnalGroup(NOBELIUM,0){}
Nobelium(const int* Members):FxnalGroup(NOBELIUM,0,1,Members){}
};
class Lawrencium:public FxnalGroup{
public:
Lawrencium():FxnalGroup(LAWRENCIUM,0){}
Lawrencium(const int* Members):FxnalGroup(LAWRENCIUM,0,1,Members){}
};
class Rutherfordium:public FxnalGroup{
public:
Rutherfordium():FxnalGroup(RUTHERFORDIUM,0){}
Rutherfordium(const int* Members):FxnalGroup(RUTHERFORDIUM,0,1,Members){}
};
class Dubnium:public FxnalGroup{
public:
Dubnium():FxnalGroup(DUBNIUM,0){}
Dubnium(const int* Members):FxnalGroup(DUBNIUM,0,1,Members){}
};
class Seaborgium:public FxnalGroup{
public:
Seaborgium():FxnalGroup(SEABORGIUM,0){}
Seaborgium(const int* Members):FxnalGroup(SEABORGIUM,0,1,Members){}
};
class Bohrium:public FxnalGroup{
public:
Bohrium():FxnalGroup(BOHRIUM,0){}
Bohrium(const int* Members):FxnalGroup(BOHRIUM,0,1,Members){}
};
class Hassium:public FxnalGroup{
public:
Hassium():FxnalGroup(HASSIUM,0){}
Hassium(const int* Members):FxnalGroup(HASSIUM,0,1,Members){}
};
class Meitnerium:public FxnalGroup{
public:
Meitnerium():FxnalGroup(MEITNERIUM,0){}
Meitnerium(const int* Members):FxnalGroup(MEITNERIUM,0,1,Members){}
};
class Darmstadtium:public FxnalGroup{
public:
Darmstadtium():FxnalGroup(DARMSTADTIUM,0){}
Darmstadtium(const int* Members):FxnalGroup(DARMSTADTIUM,0,1,Members){}
};
class Roentgenium:public FxnalGroup{
public:
Roentgenium():FxnalGroup(ROENTGENIUM,0){}
Roentgenium(const int* Members):FxnalGroup(ROENTGENIUM,0,1,Members){}
};
class Copernicium:public FxnalGroup{
public:
Copernicium():FxnalGroup(COPERNICIUM,0){}
Copernicium(const int* Members):FxnalGroup(COPERNICIUM,0,1,Members){}
};
class Ununtrium:public FxnalGroup{
public:
Ununtrium():FxnalGroup(UNUNTRIUM,0){}
Ununtrium(const int* Members):FxnalGroup(UNUNTRIUM,0,1,Members){}
};
class Flerovium:public FxnalGroup{
public:
Flerovium():FxnalGroup(FLEROVIUM,0){}
Flerovium(const int* Members):FxnalGroup(FLEROVIUM,0,1,Members){}
};
class Ununpentium:public FxnalGroup{
public:
Ununpentium():FxnalGroup(UNUNPENTIUM,0){}
Ununpentium(const int* Members):FxnalGroup(UNUNPENTIUM,0,1,Members){}
};
class Livermorium:public FxnalGroup{
public:
Livermorium():FxnalGroup(LIVERMORIUM,0){}
Livermorium(const int* Members):FxnalGroup(LIVERMORIUM,0,1,Members){}
};
class Ununseptium:public FxnalGroup{
public:
Ununseptium():FxnalGroup(UNUNSEPTIUM,0){}
Ununseptium(const int* Members):FxnalGroup(UNUNSEPTIUM,0,1,Members){}
};
class Ununoctium:public FxnalGroup{
public:
Ununoctium():FxnalGroup(UNUNOCTIUM,0){}
Ununoctium(const int* Members):FxnalGroup(UNUNOCTIUM,0,1,Members){}
};



}}//End namespaces




#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ATOMICTYPES_HH_ */
