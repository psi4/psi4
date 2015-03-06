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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FXNALGROUPTYPES_HH_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FXNALGROUPTYPES_HH_

namespace psi{
namespace LibMolecule{

/** \brief The recognized functional groups
 *
 * - Primitive Groups: Here is what I mean by each term:
 *   - methane (primary-1 carbon) CH4
 *   - methyl (primary carbon) R-CH3
 *   - methene (secondary carbon) R-CH2-R`
 *   - alkenyl1 (primary alkenyl) R=CH2
 *   - methyne (tertiary carbon) R-(CH-R`)-R```
 *   - alkenyl2 (seconary alkenyl) R=CH-R'
 *   - alkynyl1 (primary alkenyl)  R <triple-bond> CH
 *   - carbon4 (quaternary carbon) R-((C-R`)-R``)-R```
 *   - alkenyl3 (tertiary alkenyl) R=(C-R')-R''
 *   - alkynyl2 (secondary alynyl) R <triple-bond>C-R'
 *   - water   (primary-1 oxygen) OH2
 *   - hydroxyl (priarmy oxygen) R-OH
 *   - oxygen2 (secondary oxygen) R-O-R`
 *   - oxydb   (primary oxygen doubl bond) R=O
 *   - ammonia (primary-1 amine) NH3
 *   - amine1 (primary amine) R-NH2
 *   - amine2 (secondary amine) R-HN-R`
 *   - nitrodb1 (primary N double bond) R=N-H
 *   - amine3 (tertiary amine) R-(N-R`)-R``
 *   - nitrodb2 (secondary double N bond) R-N=R`
 *   - nitrotb (primary triple N bond) R <triple bond> N
 * -Derived Groups.  These are functional groups that can be written as
 *   combinations of primitive groups or in some cases other derived groups.
 *    - hydrogen cyanide (alkynyl1)-(nitrotb)
 *    - nitrile R- (alkynyl2) -(nitrotb)
 *    - formaldehyde (alkenyl1)-(oxydb)
 *    - aldehyde R-(alkenyl2)-(oxydb)
 *    - carbonyl R-[(alkenyl3)-(oxydb)]-R`
 *    - carboxyl R-(carbonyl)-(hydroxyl)
 *    - hydroperoxy R-(secondary oxygen)-(hydroxyl)
 *    - peroxide (hydroxyl)-(hydroxyl)
 *    - methoxy R-(secondary oxygen)(methyl)
 *    - methanol (hydroxyl)-(methyl)
 *    - ccdb4 (quatanary double-bond)
 *            R-((alkenyl3)-R')-((alkenyl3)-R'')-R'''
 *    - ccdb3 (ternary double-bond)
 *            R-(alkenyl2)-((alkenyl3)-R')-R''
 *    - ethenyl2 (secondary ethenyl)(alkenyl1)-((alkenyl3)-R)-R''
 *    - ccdb2 (seconardy double bond)
 *            R-(alkenyl2)-(alkenyl2)-R'
 *    - ethenyl1 (primary double bond)
 *            R-(alkenyl2)-(alkenyl1)
 *    - ethene (alkenyl1)-(alkenyl1)
 *    - cctb (triple-bond) R-(alkynyl2)-(alkynyl2)-R'
 *    - ethynyl R-(alkynyl2)-(alkynyl1)
 *    - ethyne (alkynyl1)-(alkynyl1)
 *    - ketimine1 (primary ketimine) R-((alkenyl3)-R')-(nitrodb1)
 *    - ketimine2 (secondary ketimine) R-(alkenyl3)-R')-(nitrodb2)-R''
 *    - aldimine1 (primary aldimine) R-(alkenyl2)-(nitrodb1)
 *    - aldimine2 (secondary aldimine) R-(alkenyl2)-(nitrodb2)-R'
 *    - methanimine (alkenyl1)-(nitrodb1)
 *    - aromaticring
 *
 */
enum FxnGroup_t{NO_GROUP,//The "NULL" Group, useful as a default
  //An individual atom lying by themselves
   HYDROGEN,HELIUM,LITHIUM,BERYLLIUM,BORON,CARBON,NITROGEN,OXYGEN,FLUORINE,
   NEON,SODIUM,MAGNESIUM,ALUMINIUM,SILICON,PHOSPHORUS,SULFUR,CHLORINE,ARGON,
   POTASSIUM,CALCIUM,SCANDIUM,TITANIUM,VANADIUM,CHROMIUM,MANGANESE,IRON,
   COBALT,NICKEL,COPPER,ZINC,GALLIUM,GERMANIUM,ARSENIC,SELENIUM,BROMINE,
   KRYPTON,RUBIDIUM,STRONTIUM,YTTRIUM,ZIRCONIUM,NIOBIUM,MOLYBDENUM,
   TECHNETIUM,RUTHENIUM,RHODIUM,PALLADIUM,SILVER,CADMIUM,INDIUM,TIN,
   ANTIMONY,TELLURIUM,IODINE,XENON,CESIUM,BARIUM,LANTHANUM,CERIUM,
   PRASEODYMIUM,NEODYMIUM,PROMETHIUM,SAMARIUM,EUROPIUM,GADOLINIUM,TERBIUM,
   DYSPROSIUM,HOLMIUM,ERBIUM,THULIUM,YTTERBIUM,LUTETIUM,HAFNIUM,TANTALUM,
   TUNGSTEN,RHENIUM,OSMIUM,IRIDIUM,PLATINUM,GOLD,MERCURY,THALLIUM,LEAD,
   BISMUTH,POLONIUM,ASTATINE,RADON,FRANCIUM,RADIUM,ACTINIUM,THORIUM,
   PROTACTINIUM,URANIUM,NEPTUNIUM,PLUTONIUM,AMERICIUM,CURIUM,BERKELIUM,
   CALIFORNIUM,EINSTEINIUM,FERMIUM,MENDELEVIUM,NOBELIUM,LAWRENCIUM,
   RUTHERFORDIUM,DUBNIUM,SEABORGIUM,BOHRIUM,HASSIUM,MEITNERIUM,
   DARMSTADTIUM,ROENTGENIUM,COPERNICIUM,UNUNTRIUM,FLEROVIUM,UNUNPENTIUM,
   LIVERMORIUM,UNUNSEPTIUM,UNUNOCTIUM,
  //Primitive Groups by atom:
  METHANE,
  METHYL,
  METHENE,ALKENYL1,
  METHYNE,ALKENYL2,ALKYNYL1,
  CARBON4,ALKENYL3,ALKYNYL2,//C
  WATER,
  HYDROXYL,
  OXYGEN2,OXYDB,//O
  AMMONIA,
  AMINE1,
  AMINE2,NITRODB1,
  AMINE3,NITRODB2,NITROTB,//N
  HYDROGENFLUORIDE,FLUORINE1,//F
  HYDROGENCHLORIDE,CHLORINE1,//Cl
  HYDROGENBROMIDE,BROMINE1,//Br
  HYDROGENIODIDE,IODINE1,//I
  //Derived Groups in a relatively random order:
  HYDROGENCYANIDE,NITRILE,//C w/ amine3
  FORMALDEHYDE,ALDEHYDE,CARBONYL, //C w/ oxygen2
  CARBOXYL,PEROXIDE,HYDROPEROXY,METHOXY,
  METHANOL,//Other things w/oxygen2
  CCDB4,CCDB3,CCDB2,ETHENYL1,ETHENYL2,ETHENE,
  CCTB,ETHYNYL,ETHYNE,
  KETIMINE1,KETIMINE2,ALDIMINE1,ALDIMINE2,METHANIMINE,
  AROMATICRING
                /*Still need added
                 *
                 * Haloformyl R(C=O)X
                 * Carbonate ester ROC=OOR'
                 * Carboxylate R(C-O)O
                 * Ester R(C=O)OR'
                 * Peroxy ROOR`
                 * ether ROR`
                 *
                 */
};

}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FXNALGROUPTYPES_HH_ */
