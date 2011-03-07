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
  private:
    bool loaded;
    std::map<std::string,double> element_to_Z;

  public:
    Element_to_Z() { loaded = false; };
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
    element_to_Z["NA"]      =11.0;
    element_to_Z["SODIUM"]  =11.00;
    element_to_Z["MG"]      =12.0;
    element_to_Z["MAGNESIUM"]=12.0;
    element_to_Z["AL"]      = 13.0;
    element_to_Z["ALUMINUM"]= 13.00;
    element_to_Z["SI"]      =14.0;
    element_to_Z["SILICON"]= 14.0;
    element_to_Z["P"]      =    15.0;
    element_to_Z["PHOSPHORUS"]= 15.0;
    element_to_Z["S"]      = 16.0;
    element_to_Z["SULPHUR"]= 16.0;
    element_to_Z["SULFUR"] = 16.0;
    element_to_Z["CL"]      = 17.0;
    element_to_Z["CHLORINE"]= 17.0;
    element_to_Z["AR"]   = 18.0;
    element_to_Z["ARGON"]= 18.0;
    element_to_Z["K"]      = 19.0;
    element_to_Z["POTASSIUM"]= 19.0;
    element_to_Z["CA"]      = 20.0;
    element_to_Z["CALCIUM"]= 20.0;
    element_to_Z["SC"]      = 21.0;
    element_to_Z["SCANDIUM"]= 21.0;
    element_to_Z["TI"]      = 22.0;
    element_to_Z["TITANIUM"]= 22.0;
    element_to_Z["V"]       = 23.0;
    element_to_Z["VANADIUM"]= 23.00;
    element_to_Z["CR"]      = 24.0;
    element_to_Z["CHROMIUM"]= 24.0;
    element_to_Z["MN"]      = 25.0;
    element_to_Z["MANGANESE"]=25.0;
    element_to_Z["FE"]  = 26.0;
    element_to_Z["IRON"]= 26.0;
    element_to_Z["CO"]    = 27.0;
    element_to_Z["COBALT"]= 27.0;
    element_to_Z["NI"]    = 28.0;
    element_to_Z["NICKEL"]= 28.0;
    element_to_Z["CU"]    = 29.0;
    element_to_Z["COPPER"]= 29.0;
    element_to_Z["ZN"]    = 30.0;
    element_to_Z["ZINC"]  = 30.0;
    element_to_Z["GA"]    = 31.0;
    element_to_Z["GALLIUM"]=31.0;
    element_to_Z["GE"]       = 32.0;
    element_to_Z["GERMANIUM"]= 32.0;
    element_to_Z["AS"]     = 33.0;
    element_to_Z["ARSENIC"]= 33.0;
    element_to_Z["SE"]      = 34.0;
    element_to_Z["SELENIUM"]= 34.0;
    element_to_Z["BR"]     = 35.0;
    element_to_Z["BROMINE"]= 35.0;
    element_to_Z["KR"]      = 36.0;
    element_to_Z["KRYPTON"]= 36.0;
    element_to_Z["RB"]      = 37.0;
    element_to_Z["RUBIDIUM"]= 37.0;
    element_to_Z["SR"]      =
    element_to_Z["STRONTIUM"]= 38.00;
    element_to_Z["Y"]      =
    element_to_Z["YTTRIUM"]= 39.00;
    element_to_Z["ZR"]      =
    element_to_Z["ZIRCONIUM"]= 40.00;
    element_to_Z["NB"]      =
    element_to_Z["NIOBIUM"]= 41.00;
    element_to_Z["MO"]      =
    element_to_Z["MOLYBDENUM"]= 42.00;
    element_to_Z["TC"]      =
    element_to_Z["TECHNETIUM"]= 43.00;
    element_to_Z["RU"]      =
    element_to_Z["RUTHENIUM"]= 44.00;
    element_to_Z["RH"]      =
    element_to_Z["RHODIUM"]= 45.00;
    element_to_Z["PD"]      =
    element_to_Z["PALLADIUM"]= 46.00;
    element_to_Z["AG"]      =
    element_to_Z["SILVER"]= 47.00;
    element_to_Z["CD"]      =
    element_to_Z["CADMIUM"]= 48.00;
    element_to_Z["IN"]      =
    element_to_Z["INDIUM"]= 49.00;
    element_to_Z["SN"]      =
    element_to_Z["TIN"]= 50.00;
    element_to_Z["SB"]      =
    element_to_Z["ANTIMONY"]= 51.00;
    element_to_Z["TE"]      =
    element_to_Z["TELLURIUM"]= 52.00;
    element_to_Z["I"]      =
    element_to_Z["IODINE"]= 53.00;
    element_to_Z["XE"]      =
    element_to_Z["XENON"]= 54.00;
    element_to_Z["CS"]      =
    element_to_Z["CESIUM"]= 55.00;
    element_to_Z["BA"]      =
    element_to_Z["BARIUM"]= 56.00;
    element_to_Z["LA"]      =
    element_to_Z["LANTHANUM"]= 57.00;
    element_to_Z["CE"]      =
    element_to_Z["CERIUM"]= 58.00;
    element_to_Z["PR"]      =
    element_to_Z["PRASEODYMIUM"]= 59.00;
    element_to_Z["ND"]      =
    element_to_Z["NEODYMIUM"]= 60.00;
    element_to_Z["PM"]      =
    element_to_Z["PROMETHIUM"]= 61.00;
    element_to_Z["SM"]      =
    element_to_Z["SAMARIUM"]= 62.00;
    element_to_Z["EU"]      =
    element_to_Z["EUROPIUM"]= 63.00;
    element_to_Z["GD"]      =
    element_to_Z["GADOLINIUM"]= 64.00;
    element_to_Z["TB"]      =
    element_to_Z["TERBIUM"]= 65.00;
    element_to_Z["DY"]      =
    element_to_Z["DYSPROSIUM"]= 66.00;
    element_to_Z["HO"]      =
    element_to_Z["HOLMIUM"]= 67.00;
    element_to_Z["ER"]      =
    element_to_Z["ERBIUM"]= 68.00;
    element_to_Z["TM"]      =
    element_to_Z["THULIUM"]= 69.00;
    element_to_Z["YB"]      =
    element_to_Z["YTTERBIUM"]= 70.00;
    element_to_Z["LU"]      =
    element_to_Z["LUTETIUM"]= 71.00;
    element_to_Z["HF"]      =
    element_to_Z["HAFNIUM"]= 72.00;
    element_to_Z["TA"]      =
    element_to_Z["TANTALUM"]= 73.00;
    element_to_Z["W"]      =
    element_to_Z["TUNGSTEN"]= 74.00;
    element_to_Z["RE"]      =
    element_to_Z["RHENIUM"]= 75.00;
    element_to_Z["OS"]      =
    element_to_Z["OSMIUM"]= 76.00;
    element_to_Z["IR"]      =
    element_to_Z["IRIDIUM"]= 77.00;
    element_to_Z["PT"]      =
    element_to_Z["PLATINUM"]= 78.00;
    element_to_Z["AU"]      =
    element_to_Z["GOLD"]= 79.00;
    element_to_Z["HG"]      =
    element_to_Z["MERCURY"]= 80.00;
    element_to_Z["TL"]      =
    element_to_Z["THALLIUM"]= 81.00;
    element_to_Z["PB"]      =
    element_to_Z["LEAD"]= 82.00;
    element_to_Z["BI"]      =
    element_to_Z["BISMUTH"]= 83.00;
    element_to_Z["PO"]      =
    element_to_Z["POLONIUM"]= 84.00;
    element_to_Z["AT"]      =
    element_to_Z["ASTATINE"]= 85.00;
    element_to_Z["RN"]      =
    element_to_Z["RADON"]= 86.00;
    element_to_Z["FR"]      =
    element_to_Z["FRANCIUM"]= 87.00;
    element_to_Z["RA"]      =
    element_to_Z["RADIUM"]= 88.00;
    element_to_Z["AC"]      =
    element_to_Z["ACTINIUM"]= 89.00;
    element_to_Z["TH"]      =
    element_to_Z["THORIUM"]= 90.00;
    element_to_Z["PA"]      =
    element_to_Z["PROTACTINIUM"]= 91.00;
    element_to_Z["U"]      =
    element_to_Z["URANIUM"]= 92.00;
    element_to_Z["NP"]      =
    element_to_Z["NEPTUNIUM"]= 93.00;
    element_to_Z["PU"]      =
    element_to_Z["PLUTONIUM"]= 94.00;
    element_to_Z["AM"]      =
    element_to_Z["AMERICIUM"]= 95.00;
    element_to_Z["CM"]      =
    element_to_Z["CURIUM"]= 96.00;
    element_to_Z["BK"]      =
    element_to_Z["BERKELIUM"]= 97.00;
    element_to_Z["CF"]      =
    element_to_Z["CALIFORNIUM"]= 98.00;
    element_to_Z["ES"]      =
    element_to_Z["EINSTEINIUM"]= 99.00;
    element_to_Z["FM"]      =
    element_to_Z["FERMIUM"]= 100.00;
    element_to_Z["MD"]      =
    element_to_Z["MENDELEVIUM"]= 101.00;
    element_to_Z["NO"]      =
    element_to_Z["NOBELIUM"]= 102.00;
    element_to_Z["LR"]      =
    element_to_Z["LAWRENCIUM"]= 103.00;
    element_to_Z["UNQ"]  = 104.00;
    element_to_Z["UNP"] = 105.00;
    element_to_Z["UNH"]  = 106.00;
    element_to_Z["UNS"]  = 107.00;
  }

};

}

#endif
