/*! 
  \file element_to_Z.h
  \ingroup
  \brief convert element name or symbol to atomic number
*/

#ifndef _psi_include_element_to_Z_h_
#define _psi_include_element_to_Z_h_

#include <string>
#include <map>

namespace psi {

class Element_to_Z {

  public:
    Element_to_Z() { loaded = false; };
    double operator[](const std::string & elem_sym);

  private:
    bool loaded;
    std::map<std::string,double> element_to_Z;

    double operator[](const std::string & elem_sym) {
      if (!loaded) {
        load_values();
        loaded = true;
      }
      return element_to_Z[elem_sym];
    }

    void load_values(void) {
    loaded = true;
    element_to_Z["G"]       = 0.0;
    element_to_Z["GHOST"]   = 0.0;
    element_to_Z["H"]       = 1.0;
    element_to_Z["HYDROGEN"]= 1.0;
    element_to_Z["HE"]      = 2.0;
    element_to_Z["HELIUM"]  = 2.0;
    element_to_Z["LI"]      = 3.0;
    element_to_Z["LITHIUM"] = 3.0;
    element_to_Z["BE"]      = 4.0;
    element_to_Z["BERYLLIUM"]=4.0;
    element_to_Z["B"]       = 5.0;
    element_to_Z["BORON"]   = 5.0;
    element_to_Z["C"]       = 6.0;
    element_to_Z["CARBON"]  = 6.0;
    element_to_Z["N"]       = 7.0;
    element_to_Z["NITROGEN"]= 7.0;
    element_to_Z["O"]       = 8.0;
    element_to_Z["OXYGEN"]  = 8.0;
    element_to_Z["F"]       = 9.0;
    element_to_Z["FLUORINE"]= 9.0;
    element_to_Z["NE"]      =10.0;
    element_to_Z["NEON"]    =10.0;
/*
"NA") || !strcmp(A,"SODIUM")){ *C = 11.00;
"MG") || !strcmp(A,"MAGNESIUM")){ *C = 12.00;
"AL") || !strcmp(A,"ALUMINUM")){ *C = 13.00;
"SI") || !strcmp(A,"SILICON")){ *C = 14.00;
"P") || !strcmp(A,"PHOSPHORUS")){ *C = 15.00;
"S") || !strcmp(A,"SULPHUR") || !strcmp(A,"SULFUR")){ *C = 16.00;
"CL") || !strcmp(A,"CHLORINE")){ *C = 17.00;
"AR") || !strcmp(A,"ARGON")){ *C = 18.00;
"K") || !strcmp(A,"POTASSIUM")){ *C = 19.00;
"CA") || !strcmp(A,"CALCIUM")){ *C = 20.00;
"SC") || !strcmp(A,"SCANDIUM")){ *C = 21.00;
"TI") || !strcmp(A,"TITANIUM")){ *C = 22.00;
"V") || !strcmp(A,"VANADIUM")){ *C = 23.00;
"CR") || !strcmp(A,"CHROMIUM")){ *C = 24.00;
"MN") || !strcmp(A,"MANGANESE")){ *C = 25.00;
"FE") || !strcmp(A,"IRON")){ *C = 26.00;
"CO") || !strcmp(A,"COBALT")){ *C = 27.00;
"NI") || !strcmp(A,"NICKEL")){ *C = 28.00;
"CU") || !strcmp(A,"COPPER")){ *C = 29.00;
"ZN") || !strcmp(A,"ZINC")){ *C = 30.00;
"GA") || !strcmp(A,"GALLIUM")){ *C = 31.00;
"GE") || !strcmp(A,"GERMANIUM")){ *C = 32.00;
"AS") || !strcmp(A,"ARSENIC")){ *C = 33.00;
"SE") || !strcmp(A,"SELENIUM")){ *C = 34.00;
"BR") || !strcmp(A,"BROMINE")){ *C = 35.00;
"KR") || !strcmp(A,"KRYPTON")){ *C = 36.00;
"RB") || !strcmp(A,"RUBIDIUM")){ *C = 37.00;
"SR") || !strcmp(A,"STRONTIUM")){ *C = 38.00;
"Y") || !strcmp(A,"YTTRIUM")){ *C = 39.00;
"ZR") || !strcmp(A,"ZIRCONIUM")){ *C = 40.00;
"NB") || !strcmp(A,"NIOBIUM")){ *C = 41.00;
"MO") || !strcmp(A,"MOLYBDENUM")){ *C = 42.00;
"TC") || !strcmp(A,"TECHNETIUM")){ *C = 43.00;
"RU") || !strcmp(A,"RUTHENIUM")){ *C = 44.00;
"RH") || !strcmp(A,"RHODIUM")){ *C = 45.00;
"PD") || !strcmp(A,"PALLADIUM")){ *C = 46.00;
"AG") || !strcmp(A,"SILVER")){ *C = 47.00;
"CD") || !strcmp(A,"CADMIUM")){ *C = 48.00;
"IN") || !strcmp(A,"INDIUM")){ *C = 49.00;
"SN") || !strcmp(A,"TIN")){ *C = 50.00;
"SB") || !strcmp(A,"ANTIMONY")){ *C = 51.00;
"TE") || !strcmp(A,"TELLURIUM")){ *C = 52.00;
"I") || !strcmp(A,"IODINE")){ *C = 53.00;
"XE") || !strcmp(A,"XENON")){ *C = 54.00;
"CS") || !strcmp(A,"CESIUM")){ *C = 55.00;
"BA") || !strcmp(A,"BARIUM")){ *C = 56.00;
"LA") || !strcmp(A,"LANTHANUM")){ *C = 57.00;
"CE") || !strcmp(A,"CERIUM")){ *C = 58.00;
"PR") || !strcmp(A,"PRASEODYMIUM")){ *C = 59.00;
"ND") || !strcmp(A,"NEODYMIUM")){ *C = 60.00;
"PM") || !strcmp(A,"PROMETHIUM")){ *C = 61.00;
"SM") || !strcmp(A,"SAMARIUM")){ *C = 62.00;
"EU") || !strcmp(A,"EUROPIUM")){ *C = 63.00;
"GD") || !strcmp(A,"GADOLINIUM")){ *C = 64.00;
"TB") || !strcmp(A,"TERBIUM")){ *C = 65.00;
"DY") || !strcmp(A,"DYSPROSIUM")){ *C = 66.00;
"HO") || !strcmp(A,"HOLMIUM")){ *C = 67.00;
"ER") || !strcmp(A,"ERBIUM")){ *C = 68.00;
"TM") || !strcmp(A,"THULIUM")){ *C = 69.00;
"YB") || !strcmp(A,"YTTERBIUM")){ *C = 70.00;
"LU") || !strcmp(A,"LUTETIUM")){ *C = 71.00;
"HF") || !strcmp(A,"HAFNIUM")){ *C = 72.00;
"TA") || !strcmp(A,"TANTALUM")){ *C = 73.00;
"W") || !strcmp(A,"TUNGSTEN")){ *C = 74.00;
"RE") || !strcmp(A,"RHENIUM")){ *C = 75.00;
"OS") || !strcmp(A,"OSMIUM")){ *C = 76.00;
"IR") || !strcmp(A,"IRIDIUM")){ *C = 77.00;
"PT") || !strcmp(A,"PLATINUM")){ *C = 78.00;
"AU") || !strcmp(A,"GOLD")){ *C = 79.00;
"HG") || !strcmp(A,"MERCURY")){ *C = 80.00;
"TL") || !strcmp(A,"THALLIUM")){ *C = 81.00;
"PB") || !strcmp(A,"LEAD")){ *C = 82.00;
"BI") || !strcmp(A,"BISMUTH")){ *C = 83.00;
"PO") || !strcmp(A,"POLONIUM")){ *C = 84.00;
"AT") || !strcmp(A,"ASTATINE")){ *C = 85.00;
"RN") || !strcmp(A,"RADON")){ *C = 86.00;
"FR") || !strcmp(A,"FRANCIUM")){ *C = 87.00;
"RA") || !strcmp(A,"RADIUM")){ *C = 88.00;
"AC") || !strcmp(A,"ACTINIUM")){ *C = 89.00;
"TH") || !strcmp(A,"THORIUM")){ *C = 90.00;
"PA") || !strcmp(A,"PROTACTINIUM")){ *C = 91.00;
"U") || !strcmp(A,"URANIUM")){ *C = 92.00;
"NP") || !strcmp(A,"NEPTUNIUM")){ *C = 93.00;
"PU") || !strcmp(A,"PLUTONIUM")){ *C = 94.00;
"AM") || !strcmp(A,"AMERICIUM")){ *C = 95.00;
"CM") || !strcmp(A,"CURIUM")){ *C = 96.00;
"BK") || !strcmp(A,"BERKELIUM")){ *C = 97.00;
"CF") || !strcmp(A,"CALIFORNIUM")){ *C = 98.00;
"ES") || !strcmp(A,"EINSTEINIUM")){ *C = 99.00;
"FM") || !strcmp(A,"FERMIUM")){ *C = 100.00;
"MD") || !strcmp(A,"MENDELEVIUM")){ *C = 101.00;
"NO") || !strcmp(A,"NOBELIUM")){ *C = 102.00;
"LR") || !strcmp(A,"LAWRENCIUM")){ *C = 103.00;
"UNQ")){ *C = 104.00;
"UNP")){ *C = 105.00;
"UNH")){ *C = 106.00;
"UNS")){ *C = 107.00;
*/
}

}

}

