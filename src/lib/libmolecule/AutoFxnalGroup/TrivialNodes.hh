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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_TRIVIALNODES_HH_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_TRIVIALNODES_HH_

#include "Node.h"
namespace psi{
namespace LibMolecule{
class Hydrogen: public Node {
   public:
      Hydrogen(int i=0):Node(i,"H","HYDROGEN"){}
unsigned Z()const{return 1;}
};
class Helium: public Node {
   public:
      Helium(int i=0):Node(i,"HE","HELIUM"){}
unsigned Z()const{return 2;}
};
class Lithium: public Node {
   public:
      Lithium(int i=0):Node(i,"LI","LITHIUM"){}
unsigned Z()const{return 3;}
};
class Beryllium: public Node {
   public:
      Beryllium(int i=0):Node(i,"BE","BERYLLIUM"){}
unsigned Z()const{return 4;}
};
class Boron: public Node {
   public:
      Boron(int i=0):Node(i,"B","BORON"){}
unsigned Z()const{return 5;}
};
class Carbon: public Node {
   public:
      Carbon(int i=0):Node(i,"C","CARBON"){}
unsigned Z()const{return 6;}
};
class Nitrogen: public Node {
   public:
      Nitrogen(int i=0):Node(i,"N","NITROGEN"){}
unsigned Z()const{return 7;}
};
class Oxygen: public Node {
   public:
      Oxygen(int i=0):Node(i,"O","OXYGEN"){}
unsigned Z()const{return 8;}
};
class Fluorine: public Node {
   public:
      Fluorine(int i=0):Node(i,"F","FLUORINE"){}
unsigned Z()const{return 9;}
};
class Neon: public Node {
   public:
      Neon(int i=0):Node(i,"NE","NEON"){}
unsigned Z()const{return 10;}
};
class Sodium: public Node {
   public:
      Sodium(int i=0):Node(i,"NA","SODIUM"){}
unsigned Z()const{return 11;}
};
class Magnesium: public Node {
   public:
      Magnesium(int i=0):Node(i,"MG","MAGNESIUM"){}
unsigned Z()const{return 12;}
};
class Aluminum: public Node {
   public:
      Aluminum(int i=0):Node(i,"AL","ALUMINUM"){}
unsigned Z()const{return 13;}
};
class Silicon: public Node {
   public:
      Silicon(int i=0):Node(i,"SI","SILICON"){}
unsigned Z()const{return 14;}
};
class Phosphorus: public Node {
   public:
      Phosphorus(int i=0):Node(i,"P","PHOSPHORUS"){}
unsigned Z()const{return 15;}
};
class Sulfur: public Node {
   public:
      Sulfur(int i=0):Node(i,"S","SULFUR"){}
unsigned Z()const{return 16;}
};
class Chlorine: public Node {
   public:
      Chlorine(int i=0):Node(i,"CL","CHLORINE"){}
unsigned Z()const{return 17;}
};
class Argon: public Node {
   public:
      Argon(int i=0):Node(i,"AR","ARGON"){}
unsigned Z()const{return 18;}
};
class Potassium: public Node {
   public:
      Potassium(int i=0):Node(i,"K","POTASSIUM"){}
unsigned Z()const{return 19;}
};
class Calcium: public Node {
   public:
      Calcium(int i=0):Node(i,"CA","CALCIUM"){}
unsigned Z()const{return 20;}
};
class Scandium: public Node {
   public:
      Scandium(int i=0):Node(i,"SC","SCANDIUM"){}
unsigned Z()const{return 21;}
};
class Titanium: public Node {
   public:
      Titanium(int i=0):Node(i,"TI","TITANIUM"){}
unsigned Z()const{return 22;}
};
class Vanadium: public Node {
   public:
      Vanadium(int i=0):Node(i,"V","VANADIUM"){}
unsigned Z()const{return 23;}
};
class Chromium: public Node {
   public:
      Chromium(int i=0):Node(i,"CR","CHROMIUM"){}
unsigned Z()const{return 24;}
};
class Manganese: public Node {
   public:
      Manganese(int i=0):Node(i,"MN","MANGANESE"){}
unsigned Z()const{return 25;}
};
class Iron: public Node {
   public:
      Iron(int i=0):Node(i,"FE","IRON"){}
unsigned Z()const{return 26;}
};
class Cobalt: public Node {
   public:
      Cobalt(int i=0):Node(i,"CO","COBALT"){}
unsigned Z()const{return 27;}
};
class Nickel: public Node {
   public:
      Nickel(int i=0):Node(i,"NI","NICKEL"){}
unsigned Z()const{return 28;}
};
class Copper: public Node {
   public:
      Copper(int i=0):Node(i,"CU","COPPER"){}
unsigned Z()const{return 29;}
};
class Zinc: public Node {
   public:
      Zinc(int i=0):Node(i,"ZN","ZINC"){}
unsigned Z()const{return 30;}
};
class Gallium: public Node {
   public:
      Gallium(int i=0):Node(i,"GA","GALLIUM"){}
unsigned Z()const{return 31;}
};
class Germanium: public Node {
   public:
      Germanium(int i=0):Node(i,"GE","GERMANIUM"){}
unsigned Z()const{return 32;}
};
class Arsenic: public Node {
   public:
      Arsenic(int i=0):Node(i,"AS","ARSENIC"){}
unsigned Z()const{return 33;}
};
class Selenium: public Node {
   public:
      Selenium(int i=0):Node(i,"SE","SELENIUM"){}
unsigned Z()const{return 34;}
};
class Bromine: public Node {
   public:
      Bromine(int i=0):Node(i,"BR","BROMINE"){}
unsigned Z()const{return 35;}
};
class Krypton: public Node {
   public:
      Krypton(int i=0):Node(i,"KR","KRYPTON"){}
unsigned Z()const{return 36;}
};
class Rubidium: public Node {
   public:
      Rubidium(int i=0):Node(i,"RB","RUBIDIUM"){}
unsigned Z()const{return 37;}
};
class Strontium: public Node {
   public:
      Strontium(int i=0):Node(i,"SR","STRONTIUM"){}
unsigned Z()const{return 38;}
};
class Yttrium: public Node {
   public:
      Yttrium(int i=0):Node(i,"Y","YTTRIUM"){}
unsigned Z()const{return 39;}
};
class Zirconium: public Node {
   public:
      Zirconium(int i=0):Node(i,"ZR","ZIRCONIUM"){}
unsigned Z()const{return 40;}
};
class Niobium: public Node {
   public:
      Niobium(int i=0):Node(i,"NB","NIOBIUM"){}
unsigned Z()const{return 41;}
};
class Molybdenum: public Node {
   public:
      Molybdenum(int i=0):Node(i,"MO","MOLYBDENUM"){}
unsigned Z()const{return 42;}
};
class Technetium: public Node {
   public:
      Technetium(int i=0):Node(i,"TC","TECHNETIUM"){}
unsigned Z()const{return 43;}
};
class Ruthenium: public Node {
   public:
      Ruthenium(int i=0):Node(i,"RU","RUTHENIUM"){}
unsigned Z()const{return 44;}
};
class Rhodium: public Node {
   public:
      Rhodium(int i=0):Node(i,"RH","RHODIUM"){}
unsigned Z()const{return 45;}
};
class Palladium: public Node {
   public:
      Palladium(int i=0):Node(i,"PD","PALLADIUM"){}
unsigned Z()const{return 46;}
};
class Silver: public Node {
   public:
      Silver(int i=0):Node(i,"AG","SILVER"){}
unsigned Z()const{return 47;}
};
class Cadmium: public Node {
   public:
      Cadmium(int i=0):Node(i,"CD","CADMIUM"){}
unsigned Z()const{return 48;}
};
class Indium: public Node {
   public:
      Indium(int i=0):Node(i,"IN","INDIUM"){}
unsigned Z()const{return 49;}
};
class Tin: public Node {
   public:
      Tin(int i=0):Node(i,"SN","TIN"){}
unsigned Z()const{return 50;}
};
class Antimony: public Node {
   public:
      Antimony(int i=0):Node(i,"SB","ANTIMONY"){}
unsigned Z()const{return 51;}
};
class Tellurium: public Node {
   public:
      Tellurium(int i=0):Node(i,"TE","TELLURIUM"){}
unsigned Z()const{return 52;}
};
class Iodine: public Node {
   public:
      Iodine(int i=0):Node(i,"I","IODINE"){}
unsigned Z()const{return 53;}
};
class Xenon: public Node {
   public:
      Xenon(int i=0):Node(i,"XE","XENON"){}
unsigned Z()const{return 54;}
};
class Cesium: public Node {
   public:
      Cesium(int i=0):Node(i,"CS","CESIUM"){}
unsigned Z()const{return 55;}
};
class Barium: public Node {
   public:
      Barium(int i=0):Node(i,"BA","BARIUM"){}
unsigned Z()const{return 56;}
};
class Lanthanum: public Node {
   public:
      Lanthanum(int i=0):Node(i,"LA","LANTHANUM"){}
unsigned Z()const{return 57;}
};
class Cerium: public Node {
   public:
      Cerium(int i=0):Node(i,"CE","CERIUM"){}
unsigned Z()const{return 58;}
};
class Praseodymium: public Node {
   public:
      Praseodymium(int i=0):Node(i,"PR","PRASEODYMIUM"){}
unsigned Z()const{return 59;}
};
class Neodymium: public Node {
   public:
      Neodymium(int i=0):Node(i,"ND","NEODYMIUM"){}
unsigned Z()const{return 60;}
};
class Promethium: public Node {
   public:
      Promethium(int i=0):Node(i,"PM","PROMETHIUM"){}
unsigned Z()const{return 61;}
};
class Samarium: public Node {
   public:
      Samarium(int i=0):Node(i,"SM","SAMARIUM"){}
unsigned Z()const{return 62;}
};
class Europium: public Node {
   public:
      Europium(int i=0):Node(i,"EU","EUROPIUM"){}
unsigned Z()const{return 63;}
};
class Gadolinium: public Node {
   public:
      Gadolinium(int i=0):Node(i,"GD","GADOLINIUM"){}
unsigned Z()const{return 64;}
};
class Terbium: public Node {
   public:
      Terbium(int i=0):Node(i,"TB","TERBIUM"){}
unsigned Z()const{return 65;}
};
class Dysprosium: public Node {
   public:
      Dysprosium(int i=0):Node(i,"DY","DYSPROSIUM"){}
unsigned Z()const{return 66;}
};
class Holmium: public Node {
   public:
      Holmium(int i=0):Node(i,"HO","HOLMIUM"){}
unsigned Z()const{return 67;}
};
class Erbium: public Node {
   public:
      Erbium(int i=0):Node(i,"ER","ERBIUM"){}
unsigned Z()const{return 68;}
};
class Thulium: public Node {
   public:
      Thulium(int i=0):Node(i,"TM","THULIUM"){}
unsigned Z()const{return 69;}
};
class Ytterbium: public Node {
   public:
      Ytterbium(int i=0):Node(i,"YB","YTTERBIUM"){}
unsigned Z()const{return 70;}
};
class Lutetium: public Node {
   public:
      Lutetium(int i=0):Node(i,"LU","LUTETIUM"){}
unsigned Z()const{return 71;}
};
class Hafnium: public Node {
   public:
      Hafnium(int i=0):Node(i,"HF","HAFNIUM"){}
unsigned Z()const{return 72;}
};
class Tantalum: public Node {
   public:
      Tantalum(int i=0):Node(i,"TA","TANTALUM"){}
unsigned Z()const{return 73;}
};
class Tungsten: public Node {
   public:
      Tungsten(int i=0):Node(i,"W","TUNGSTEN"){}
unsigned Z()const{return 74;}
};
class Rhenium: public Node {
   public:
      Rhenium(int i=0):Node(i,"RE","RHENIUM"){}
unsigned Z()const{return 75;}
};
class Osmium: public Node {
   public:
      Osmium(int i=0):Node(i,"OS","OSMIUM"){}
unsigned Z()const{return 76;}
};
class Iridium: public Node {
   public:
      Iridium(int i=0):Node(i,"IR","IRIDIUM"){}
unsigned Z()const{return 77;}
};
class Platinum: public Node {
   public:
      Platinum(int i=0):Node(i,"PT","PLATINUM"){}
unsigned Z()const{return 78;}
};
class Gold: public Node {
   public:
      Gold(int i=0):Node(i,"AU","GOLD"){}
unsigned Z()const{return 79;}
};
class Mercury: public Node {
   public:
      Mercury(int i=0):Node(i,"HG","MERCURY"){}
unsigned Z()const{return 80;}
};
class Thallium: public Node {
   public:
      Thallium(int i=0):Node(i,"TL","THALLIUM"){}
unsigned Z()const{return 81;}
};
class Lead: public Node {
   public:
      Lead(int i=0):Node(i,"PB","LEAD"){}
unsigned Z()const{return 82;}
};
class Bismuth: public Node {
   public:
      Bismuth(int i=0):Node(i,"BI","BISMUTH"){}
unsigned Z()const{return 83;}
};
class Polonium: public Node {
   public:
      Polonium(int i=0):Node(i,"PO","POLONIUM"){}
unsigned Z()const{return 84;}
};
class Astatine: public Node {
   public:
      Astatine(int i=0):Node(i,"AT","ASTATINE"){}
unsigned Z()const{return 85;}
};
class Radon: public Node {
   public:
      Radon(int i=0):Node(i,"RN","RADON"){}
unsigned Z()const{return 86;}
};
class Francium: public Node {
   public:
      Francium(int i=0):Node(i,"FR","FRANCIUM"){}
unsigned Z()const{return 87;}
};
class Radium: public Node {
   public:
      Radium(int i=0):Node(i,"RA","RADIUM"){}
unsigned Z()const{return 88;}
};
class Actinium: public Node {
   public:
      Actinium(int i=0):Node(i,"AC","ACTINIUM"){}
unsigned Z()const{return 89;}
};
class Thorium: public Node {
   public:
      Thorium(int i=0):Node(i,"TH","THORIUM"){}
unsigned Z()const{return 90;}
};
class Protactinium: public Node {
   public:
      Protactinium(int i=0):Node(i,"PA","PROTACTINIUM"){}
unsigned Z()const{return 91;}
};
class Uranium: public Node {
   public:
      Uranium(int i=0):Node(i,"U","URANIUM"){}
unsigned Z()const{return 92;}
};
class Neptunium: public Node {
   public:
      Neptunium(int i=0):Node(i,"NP","NEPTUNIUM"){}
unsigned Z()const{return 93;}
};
class Plutonium: public Node {
   public:
      Plutonium(int i=0):Node(i,"PU","PLUTONIUM"){}
unsigned Z()const{return 94;}
};
class Americium: public Node {
   public:
      Americium(int i=0):Node(i,"AM","AMERICIUM"){}
unsigned Z()const{return 95;}
};
class Curium: public Node {
   public:
      Curium(int i=0):Node(i,"CM","CURIUM"){}
unsigned Z()const{return 96;}
};
class Berkelium: public Node {
   public:
      Berkelium(int i=0):Node(i,"BK","BERKELIUM"){}
unsigned Z()const{return 97;}
};
class Californium: public Node {
   public:
      Californium(int i=0):Node(i,"CF","CALIFORNIUM"){}
unsigned Z()const{return 98;}
};
class Einsteinium: public Node {
   public:
      Einsteinium(int i=0):Node(i,"ES","EINSTEINIUM"){}
unsigned Z()const{return 99;}
};
class Fermium: public Node {
   public:
      Fermium(int i=0):Node(i,"FM","FERMIUM"){}
unsigned Z()const{return 100;}
};
class Mendelevium: public Node {
   public:
      Mendelevium(int i=0):Node(i,"MD","MENDELEVIUM"){}
unsigned Z()const{return 101;}
};
class Nobelium: public Node {
   public:
      Nobelium(int i=0):Node(i,"NO","NOBELIUM"){}
unsigned Z()const{return 102;}
};
class Lawrencium: public Node {
   public:
      Lawrencium(int i=0):Node(i,"LR","LAWRENCIUM"){}
unsigned Z()const{return 103;}
};
class Rutherfordium: public Node {
   public:
      Rutherfordium(int i=0):Node(i,"RF","RUTHERFORDIUM"){}
unsigned Z()const{return 104;}
};
class Dubnium: public Node {
   public:
      Dubnium(int i=0):Node(i,"DB","DUBNIUM"){}
unsigned Z()const{return 105;}
};
class Seaborgium: public Node {
   public:
      Seaborgium(int i=0):Node(i,"SG","SEABORGIUM"){}
unsigned Z()const{return 106;}
};
class Bohrium: public Node {
   public:
      Bohrium(int i=0):Node(i,"BH","BOHRIUM"){}
unsigned Z()const{return 107;}
};
class Hassium: public Node {
   public:
      Hassium(int i=0):Node(i,"HS","HASSIUM"){}
unsigned Z()const{return 108;}
};
class Meitnerium: public Node {
   public:
      Meitnerium(int i=0):Node(i,"MT","MEITNERIUM"){}
unsigned Z()const{return 109;}
};
class Darmstadtium: public Node {
   public:
      Darmstadtium(int i=0):Node(i,"DS","DARMSTADTIUM"){}
unsigned Z()const{return 110;}
};
class Roentgenium: public Node {
   public:
      Roentgenium(int i=0):Node(i,"RG","ROENTGENIUM"){}
unsigned Z()const{return 111;}
};
class Copernicium: public Node {
   public:
      Copernicium(int i=0):Node(i,"CN","COPERNICIUM"){}
unsigned Z()const{return 112;}
};
class Ununtrium: public Node {
   public:
      Ununtrium(int i=0):Node(i,"UUT","UNUNTRIUM"){}
unsigned Z()const{return 113;}
};
class Flerovium: public Node {
   public:
      Flerovium(int i=0):Node(i,"FL","FLEROVIUM"){}
unsigned Z()const{return 114;}
};
class Ununpentium: public Node {
   public:
      Ununpentium(int i=0):Node(i,"UUP","UNUNPENTIUM"){}
unsigned Z()const{return 115;}
};
class Livermorium: public Node {
   public:
      Livermorium(int i=0):Node(i,"LV","LIVERMORIUM"){}
unsigned Z()const{return 116;}
};
class Ununseptium: public Node {
   public:
      Ununseptium(int i=0):Node(i,"UUS","UNUNSEPTIUM"){}
unsigned Z()const{return 117;}
};
class Ununoctium: public Node {
   public:
      Ununoctium(int i=0):Node(i,"UUO","UNUNOCTIUM"){}
unsigned Z()const{return 118;}
};


}}//End namespaces
#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_TRIVIALNODES_HH_*/
