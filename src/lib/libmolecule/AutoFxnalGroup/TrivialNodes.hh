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
};
class Helium: public Node {
   public:
      Helium(int i=0):Node(i,"HE","HELIUM"){}
};
class Lithium: public Node {
   public:
      Lithium(int i=0):Node(i,"LI","LITHIUM"){}
};
class Beryllium: public Node {
   public:
      Beryllium(int i=0):Node(i,"BE","BERYLLIUM"){}
};
class Boron: public Node {
   public:
      Boron(int i=0):Node(i,"B","BORON"){}
};
class Carbon: public Node {
   public:
      Carbon(int i=0):Node(i,"C","CARBON"){}
};
class Nitrogen: public Node {
   public:
      Nitrogen(int i=0):Node(i,"N","NITROGEN"){}
};
class Oxygen: public Node {
   public:
      Oxygen(int i=0):Node(i,"O","OXYGEN"){}
};
class Fluorine: public Node {
   public:
      Fluorine(int i=0):Node(i,"F","FLUORINE"){}
};
class Neon: public Node {
   public:
      Neon(int i=0):Node(i,"NE","NEON"){}
};
class Sodium: public Node {
   public:
      Sodium(int i=0):Node(i,"NA","SODIUM"){}
};
class Magnesium: public Node {
   public:
      Magnesium(int i=0):Node(i,"MG","MAGNESIUM"){}
};
class Aluminum: public Node {
   public:
      Aluminum(int i=0):Node(i,"AL","ALUMINUM"){}
};
class Silicon: public Node {
   public:
      Silicon(int i=0):Node(i,"SI","SILICON"){}
};
class Phosphorus: public Node {
   public:
      Phosphorus(int i=0):Node(i,"P","PHOSPHORUS"){}
};
class Sulfur: public Node {
   public:
      Sulfur(int i=0):Node(i,"S","SULFUR"){}
};
class Chlorine: public Node {
   public:
      Chlorine(int i=0):Node(i,"CL","CHLORINE"){}
};
class Argon: public Node {
   public:
      Argon(int i=0):Node(i,"AR","ARGON"){}
};
class Potassium: public Node {
   public:
      Potassium(int i=0):Node(i,"K","POTASSIUM"){}
};
class Calcium: public Node {
   public:
      Calcium(int i=0):Node(i,"CA","CALCIUM"){}
};
class Scandium: public Node {
   public:
      Scandium(int i=0):Node(i,"SC","SCANDIUM"){}
};
class Titanium: public Node {
   public:
      Titanium(int i=0):Node(i,"TI","TITANIUM"){}
};
class Vanadium: public Node {
   public:
      Vanadium(int i=0):Node(i,"V","VANADIUM"){}
};
class Chromium: public Node {
   public:
      Chromium(int i=0):Node(i,"CR","CHROMIUM"){}
};
class Manganese: public Node {
   public:
      Manganese(int i=0):Node(i,"MN","MANGANESE"){}
};
class Iron: public Node {
   public:
      Iron(int i=0):Node(i,"FE","IRON"){}
};
class Cobalt: public Node {
   public:
      Cobalt(int i=0):Node(i,"CO","COBALT"){}
};
class Nickel: public Node {
   public:
      Nickel(int i=0):Node(i,"NI","NICKEL"){}
};
class Copper: public Node {
   public:
      Copper(int i=0):Node(i,"CU","COPPER"){}
};
class Zinc: public Node {
   public:
      Zinc(int i=0):Node(i,"ZN","ZINC"){}
};
class Gallium: public Node {
   public:
      Gallium(int i=0):Node(i,"GA","GALLIUM"){}
};
class Germanium: public Node {
   public:
      Germanium(int i=0):Node(i,"GE","GERMANIUM"){}
};
class Arsenic: public Node {
   public:
      Arsenic(int i=0):Node(i,"AS","ARSENIC"){}
};
class Selenium: public Node {
   public:
      Selenium(int i=0):Node(i,"SE","SELENIUM"){}
};
class Bromine: public Node {
   public:
      Bromine(int i=0):Node(i,"BR","BROMINE"){}
};
class Krypton: public Node {
   public:
      Krypton(int i=0):Node(i,"KR","KRYPTON"){}
};
class Rubidium: public Node {
   public:
      Rubidium(int i=0):Node(i,"RB","RUBIDIUM"){}
};
class Strontium: public Node {
   public:
      Strontium(int i=0):Node(i,"SR","STRONTIUM"){}
};
class Yttrium: public Node {
   public:
      Yttrium(int i=0):Node(i,"Y","YTTRIUM"){}
};
class Zirconium: public Node {
   public:
      Zirconium(int i=0):Node(i,"ZR","ZIRCONIUM"){}
};
class Niobium: public Node {
   public:
      Niobium(int i=0):Node(i,"NB","NIOBIUM"){}
};
class Molybdenum: public Node {
   public:
      Molybdenum(int i=0):Node(i,"MO","MOLYBDENUM"){}
};
class Technetium: public Node {
   public:
      Technetium(int i=0):Node(i,"TC","TECHNETIUM"){}
};
class Ruthenium: public Node {
   public:
      Ruthenium(int i=0):Node(i,"RU","RUTHENIUM"){}
};
class Rhodium: public Node {
   public:
      Rhodium(int i=0):Node(i,"RH","RHODIUM"){}
};
class Palladium: public Node {
   public:
      Palladium(int i=0):Node(i,"PD","PALLADIUM"){}
};
class Silver: public Node {
   public:
      Silver(int i=0):Node(i,"AG","SILVER"){}
};
class Cadmium: public Node {
   public:
      Cadmium(int i=0):Node(i,"CD","CADMIUM"){}
};
class Indium: public Node {
   public:
      Indium(int i=0):Node(i,"IN","INDIUM"){}
};
class Tin: public Node {
   public:
      Tin(int i=0):Node(i,"SN","TIN"){}
};
class Antimony: public Node {
   public:
      Antimony(int i=0):Node(i,"SB","ANTIMONY"){}
};
class Tellurium: public Node {
   public:
      Tellurium(int i=0):Node(i,"TE","TELLURIUM"){}
};
class Iodine: public Node {
   public:
      Iodine(int i=0):Node(i,"I","IODINE"){}
};
class Xenon: public Node {
   public:
      Xenon(int i=0):Node(i,"XE","XENON"){}
};
class Cesium: public Node {
   public:
      Cesium(int i=0):Node(i,"CS","CESIUM"){}
};
class Barium: public Node {
   public:
      Barium(int i=0):Node(i,"BA","BARIUM"){}
};
class Lanthanum: public Node {
   public:
      Lanthanum(int i=0):Node(i,"LA","LANTHANUM"){}
};
class Cerium: public Node {
   public:
      Cerium(int i=0):Node(i,"CE","CERIUM"){}
};
class Praseodymium: public Node {
   public:
      Praseodymium(int i=0):Node(i,"PR","PRASEODYMIUM"){}
};
class Neodymium: public Node {
   public:
      Neodymium(int i=0):Node(i,"ND","NEODYMIUM"){}
};
class Promethium: public Node {
   public:
      Promethium(int i=0):Node(i,"PM","PROMETHIUM"){}
};
class Samarium: public Node {
   public:
      Samarium(int i=0):Node(i,"SM","SAMARIUM"){}
};
class Europium: public Node {
   public:
      Europium(int i=0):Node(i,"EU","EUROPIUM"){}
};
class Gadolinium: public Node {
   public:
      Gadolinium(int i=0):Node(i,"GD","GADOLINIUM"){}
};
class Terbium: public Node {
   public:
      Terbium(int i=0):Node(i,"TB","TERBIUM"){}
};
class Dysprosium: public Node {
   public:
      Dysprosium(int i=0):Node(i,"DY","DYSPROSIUM"){}
};
class Holmium: public Node {
   public:
      Holmium(int i=0):Node(i,"HO","HOLMIUM"){}
};
class Erbium: public Node {
   public:
      Erbium(int i=0):Node(i,"ER","ERBIUM"){}
};
class Thulium: public Node {
   public:
      Thulium(int i=0):Node(i,"TM","THULIUM"){}
};
class Ytterbium: public Node {
   public:
      Ytterbium(int i=0):Node(i,"YB","YTTERBIUM"){}
};
class Lutetium: public Node {
   public:
      Lutetium(int i=0):Node(i,"LU","LUTETIUM"){}
};
class Hafnium: public Node {
   public:
      Hafnium(int i=0):Node(i,"HF","HAFNIUM"){}
};
class Tantalum: public Node {
   public:
      Tantalum(int i=0):Node(i,"TA","TANTALUM"){}
};
class Tungsten: public Node {
   public:
      Tungsten(int i=0):Node(i,"W","TUNGSTEN"){}
};
class Rhenium: public Node {
   public:
      Rhenium(int i=0):Node(i,"RE","RHENIUM"){}
};
class Osmium: public Node {
   public:
      Osmium(int i=0):Node(i,"OS","OSMIUM"){}
};
class Iridium: public Node {
   public:
      Iridium(int i=0):Node(i,"IR","IRIDIUM"){}
};
class Platinum: public Node {
   public:
      Platinum(int i=0):Node(i,"PT","PLATINUM"){}
};
class Gold: public Node {
   public:
      Gold(int i=0):Node(i,"AU","GOLD"){}
};
class Mercury: public Node {
   public:
      Mercury(int i=0):Node(i,"HG","MERCURY"){}
};
class Thallium: public Node {
   public:
      Thallium(int i=0):Node(i,"TL","THALLIUM"){}
};
class Lead: public Node {
   public:
      Lead(int i=0):Node(i,"PB","LEAD"){}
};
class Bismuth: public Node {
   public:
      Bismuth(int i=0):Node(i,"BI","BISMUTH"){}
};
class Polonium: public Node {
   public:
      Polonium(int i=0):Node(i,"PO","POLONIUM"){}
};
class Astatine: public Node {
   public:
      Astatine(int i=0):Node(i,"AT","ASTATINE"){}
};
class Radon: public Node {
   public:
      Radon(int i=0):Node(i,"RN","RADON"){}
};
class Francium: public Node {
   public:
      Francium(int i=0):Node(i,"FR","FRANCIUM"){}
};
class Radium: public Node {
   public:
      Radium(int i=0):Node(i,"RA","RADIUM"){}
};
class Actinium: public Node {
   public:
      Actinium(int i=0):Node(i,"AC","ACTINIUM"){}
};
class Thorium: public Node {
   public:
      Thorium(int i=0):Node(i,"TH","THORIUM"){}
};
class Protactinium: public Node {
   public:
      Protactinium(int i=0):Node(i,"PA","PROTACTINIUM"){}
};
class Uranium: public Node {
   public:
      Uranium(int i=0):Node(i,"U","URANIUM"){}
};
class Neptunium: public Node {
   public:
      Neptunium(int i=0):Node(i,"NP","NEPTUNIUM"){}
};
class Plutonium: public Node {
   public:
      Plutonium(int i=0):Node(i,"PU","PLUTONIUM"){}
};
class Americium: public Node {
   public:
      Americium(int i=0):Node(i,"AM","AMERICIUM"){}
};
class Curium: public Node {
   public:
      Curium(int i=0):Node(i,"CM","CURIUM"){}
};
class Berkelium: public Node {
   public:
      Berkelium(int i=0):Node(i,"BK","BERKELIUM"){}
};
class Californium: public Node {
   public:
      Californium(int i=0):Node(i,"CF","CALIFORNIUM"){}
};
class Einsteinium: public Node {
   public:
      Einsteinium(int i=0):Node(i,"ES","EINSTEINIUM"){}
};
class Fermium: public Node {
   public:
      Fermium(int i=0):Node(i,"FM","FERMIUM"){}
};
class Mendelevium: public Node {
   public:
      Mendelevium(int i=0):Node(i,"MD","MENDELEVIUM"){}
};
class Nobelium: public Node {
   public:
      Nobelium(int i=0):Node(i,"NO","NOBELIUM"){}
};
class Lawrencium: public Node {
   public:
      Lawrencium(int i=0):Node(i,"LR","LAWRENCIUM"){}
};
class Rutherfordium: public Node {
   public:
      Rutherfordium(int i=0):Node(i,"RF","RUTHERFORDIUM"){}
};
class Dubnium: public Node {
   public:
      Dubnium(int i=0):Node(i,"DB","DUBNIUM"){}
};
class Seaborgium: public Node {
   public:
      Seaborgium(int i=0):Node(i,"SG","SEABORGIUM"){}
};
class Bohrium: public Node {
   public:
      Bohrium(int i=0):Node(i,"BH","BOHRIUM"){}
};
class Hassium: public Node {
   public:
      Hassium(int i=0):Node(i,"HS","HASSIUM"){}
};
class Meitnerium: public Node {
   public:
      Meitnerium(int i=0):Node(i,"MT","MEITNERIUM"){}
};
class Darmstadtium: public Node {
   public:
      Darmstadtium(int i=0):Node(i,"DS","DARMSTADTIUM"){}
};
class Roentgenium: public Node {
   public:
      Roentgenium(int i=0):Node(i,"RG","ROENTGENIUM"){}
};
class Copernicium: public Node {
   public:
      Copernicium(int i=0):Node(i,"CN","COPERNICIUM"){}
};
class Ununtrium: public Node {
   public:
      Ununtrium(int i=0):Node(i,"UUT","UNUNTRIUM"){}
};
class Flerovium: public Node {
   public:
      Flerovium(int i=0):Node(i,"FL","FLEROVIUM"){}
};
class Ununpentium: public Node {
   public:
      Ununpentium(int i=0):Node(i,"UUP","UNUNPENTIUM"){}
};
class Livermorium: public Node {
   public:
      Livermorium(int i=0):Node(i,"LV","LIVERMORIUM"){}
};
class Ununseptium: public Node {
   public:
      Ununseptium(int i=0):Node(i,"UUS","UNUNSEPTIUM"){}
};
class Ununoctium: public Node {
   public:
      Ununoctium(int i=0):Node(i,"UUO","UNUNOCTIUM"){}
};

}}//End namespaces
#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_TRIVIALNODES_HH_*/
