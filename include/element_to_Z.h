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
    element_to_Z["G"]            =
    element_to_Z["GHOST"]        =   0.0;
    element_to_Z["H"]            =
    element_to_Z["HYDROGEN"]     =   1.0;
    element_to_Z["HE"]           =
    element_to_Z["HELIUM"]       =   2.0;
    element_to_Z["LI"]           =
    element_to_Z["LITHIUM"]      =   3.0;
    element_to_Z["BE"]           =
    element_to_Z["BERYLLIUM"]    =   4.0;
    element_to_Z["B"]            =
    element_to_Z["BORON"]        =   5.0;
    element_to_Z["C"]            =
    element_to_Z["CARBON"]       =   6.0;
    element_to_Z["N"]            =
    element_to_Z["NITROGEN"]     =   7.0;
    element_to_Z["O"]            =
    element_to_Z["OXYGEN"]       =   8.0;
    element_to_Z["F"]            =
    element_to_Z["FLUORINE"]     =   9.0;
    element_to_Z["NE"]           =
    element_to_Z["NEON"]         =  10.0;
    element_to_Z["NA"]           =
    element_to_Z["SODIUM"]       =  11.0;
    element_to_Z["MG"]           =
    element_to_Z["MAGNESIUM"]    =  12.0;
    element_to_Z["AL"]           =
    element_to_Z["ALUMINUM"]     =  13.0;
    element_to_Z["SI"]           =
    element_to_Z["SILICON"]      =  14.0;
    element_to_Z["P"]            =
    element_to_Z["PHOSPHORUS"]   =  15.0;
    element_to_Z["S"]            =
    element_to_Z["SULPHUR"]      =
    element_to_Z["SULFUR"]       =  16.0;
    element_to_Z["CL"]           =
    element_to_Z["CHLORINE"]     =  17.0;
    element_to_Z["AR"]           =
    element_to_Z["ARGON"]        =  18.0;
    element_to_Z["K"]            =
    element_to_Z["POTASSIUM"]    =  19.0;
    element_to_Z["CA"]           =
    element_to_Z["CALCIUM"]      =  20.0;
    element_to_Z["SC"]           =
    element_to_Z["SCANDIUM"]     =  21.0;
    element_to_Z["TI"]           =
    element_to_Z["TITANIUM"]     =  22.0;
    element_to_Z["V"]            =
    element_to_Z["VANADIUM"]     =  23.0;
    element_to_Z["CR"]           =
    element_to_Z["CHROMIUM"]     =  24.0;
    element_to_Z["MN"]           =
    element_to_Z["MANGANESE"]    =  25.0;
    element_to_Z["FE"]           =
    element_to_Z["IRON"]         =  26.0;
    element_to_Z["CO"]           =
    element_to_Z["COBALT"]       =  27.0;
    element_to_Z["NI"]           =
    element_to_Z["NICKEL"]       =  28.0;
    element_to_Z["CU"]           =
    element_to_Z["COPPER"]       =  29.0;
    element_to_Z["ZN"]           =
    element_to_Z["ZINC"]         =  30.0;
    element_to_Z["GA"]           =
    element_to_Z["GALLIUM"]      =  31.0;
    element_to_Z["GE"]           =
    element_to_Z["GERMANIUM"]    =  32.0;
    element_to_Z["AS"]           =
    element_to_Z["ARSENIC"]      =  33.0;
    element_to_Z["SE"]           =
    element_to_Z["SELENIUM"]     =  34.0;
    element_to_Z["BR"]           =
    element_to_Z["BROMINE"]      =  35.0;
    element_to_Z["KR"]           =
    element_to_Z["KRYPTON"]      =  36.0;
    element_to_Z["RB"]           =
    element_to_Z["RUBIDIUM"]     =  37.0;
    element_to_Z["SR"]           =
    element_to_Z["STRONTIUM"]    =  38.0;
    element_to_Z["Y"]            =
    element_to_Z["YTTRIUM"]      =  39.0;
    element_to_Z["ZR"]           =
    element_to_Z["ZIRCONIUM"]    =  40.0;
    element_to_Z["NB"]           =
    element_to_Z["NIOBIUM"]      =  41.0;
    element_to_Z["MO"]           =
    element_to_Z["MOLYBDENUM"]   =  42.0;
    element_to_Z["TC"]           =
    element_to_Z["TECHNETIUM"]   =  43.0;
    element_to_Z["RU"]           =
    element_to_Z["RUTHENIUM"]    =  44.0;
    element_to_Z["RH"]           =
    element_to_Z["RHODIUM"]      =  45.0;
    element_to_Z["PD"]           =
    element_to_Z["PALLADIUM"]    =  46.0;
    element_to_Z["AG"]           =
    element_to_Z["SILVER"]       =  47.0;
    element_to_Z["CD"]           =
    element_to_Z["CADMIUM"]      =  48.0;
    element_to_Z["IN"]           =
    element_to_Z["INDIUM"]       =  49.0;
    element_to_Z["SN"]           = 
    element_to_Z["TIN"]          =  50.0;
    element_to_Z["SB"]           =
    element_to_Z["ANTIMONY"]     =  51.0;
    element_to_Z["TE"]           =
    element_to_Z["TELLURIUM"]    =  52.0;
    element_to_Z["I"]            =
    element_to_Z["IODINE"]       =  53.0;
    element_to_Z["XE"]           =
    element_to_Z["XENON"]        =  54.0;
    element_to_Z["CS"]           =
    element_to_Z["CESIUM"]       =  55.0;
    element_to_Z["BA"]           =
    element_to_Z["BARIUM"]       =  56.0;
    element_to_Z["LA"]           =
    element_to_Z["LANTHANUM"]    =  57.0;
    element_to_Z["CE"]           =
    element_to_Z["CERIUM"]       =  58.0;
    element_to_Z["PR"]           =
    element_to_Z["PRASEODYMIUM"] =  59.0;
    element_to_Z["ND"]           =
    element_to_Z["NEODYMIUM"]    =  60.0;
    element_to_Z["PM"]           =
    element_to_Z["PROMETHIUM"]   =  61.0;
    element_to_Z["SM"]           =
    element_to_Z["SAMARIUM"]     =  62.0;
    element_to_Z["EU"]           =
    element_to_Z["EUROPIUM"]     =  63.0;
    element_to_Z["GD"]           =
    element_to_Z["GADOLINIUM"]   =  64.0;
    element_to_Z["TB"]           =
    element_to_Z["TERBIUM"]      =  65.0;
    element_to_Z["DY"]           =
    element_to_Z["DYSPROSIUM"]   =  66.0;
    element_to_Z["HO"]           =
    element_to_Z["HOLMIUM"]      =  67.0;
    element_to_Z["ER"]           =
    element_to_Z["ERBIUM"]       =  68.0;
    element_to_Z["TM"]           =
    element_to_Z["THULIUM"]      =  69.0;
    element_to_Z["YB"]           =
    element_to_Z["YTTERBIUM"]    =  70.0;
    element_to_Z["LU"]           =
    element_to_Z["LUTETIUM"]     =  71.0;
    element_to_Z["HF"]           =
    element_to_Z["HAFNIUM"]      =  72.0;
    element_to_Z["TA"]           =
    element_to_Z["TANTALUM"]     =  73.0;
    element_to_Z["W"]            =
    element_to_Z["TUNGSTEN"]     =  74.0;
    element_to_Z["RE"]           =
    element_to_Z["RHENIUM"]      =  75.0;
    element_to_Z["OS"]           =
    element_to_Z["OSMIUM"]       =  76.0;
    element_to_Z["IR"]           =
    element_to_Z["IRIDIUM"]      =  77.0;
    element_to_Z["PT"]           =
    element_to_Z["PLATINUM"]     =  78.0;
    element_to_Z["AU"]           =
    element_to_Z["GOLD"]         =  79.0;
    element_to_Z["HG"]           =
    element_to_Z["MERCURY"]      =  80.0;
    element_to_Z["TL"]           =
    element_to_Z["THALLIUM"]     =  81.0;
    element_to_Z["PB"]           =
    element_to_Z["LEAD"]         =  82.0;
    element_to_Z["BI"]           =
    element_to_Z["BISMUTH"]      =  83.0;
    element_to_Z["PO"]           =
    element_to_Z["POLONIUM"]     =  84.0;
    element_to_Z["AT"]           =
    element_to_Z["ASTATINE"]     =  85.0;
    element_to_Z["RN"]           =
    element_to_Z["RADON"]        =  86.0;
    element_to_Z["FR"]           =
    element_to_Z["FRANCIUM"]     =  87.0;
    element_to_Z["RA"]           =
    element_to_Z["RADIUM"]       =  88.0;
    element_to_Z["AC"]           =
    element_to_Z["ACTINIUM"]     =  89.0;
    element_to_Z["TH"]           =
    element_to_Z["THORIUM"]      =  90.0;
    element_to_Z["PA"]           =
    element_to_Z["PROTACTINIUM"] =  91.0;
    element_to_Z["U"]            =
    element_to_Z["URANIUM"]      =  92.0;
    element_to_Z["NP"]           =
    element_to_Z["NEPTUNIUM"]    =  93.0;
    element_to_Z["PU"]           =
    element_to_Z["PLUTONIUM"]    =  94.0;
    element_to_Z["AM"]           =
    element_to_Z["AMERICIUM"]    =  95.0;
    element_to_Z["CM"]           =
    element_to_Z["CURIUM"]       =  96.0;
    element_to_Z["BK"]           =
    element_to_Z["BERKELIUM"]    =  97.0;
    element_to_Z["CF"]           =
    element_to_Z["CALIFORNIUM"]  =  98.0;
    element_to_Z["ES"]           =
    element_to_Z["EINSTEINIUM"]  =  99.0;
    element_to_Z["FM"]           =
    element_to_Z["FERMIUM"]      = 100.0;
    element_to_Z["MD"]           =
    element_to_Z["MENDELEVIUM"]  = 101.0;
    element_to_Z["NO"]           =
    element_to_Z["NOBELIUM"]     = 102.0;
    element_to_Z["LR"]           =
    element_to_Z["LAWRENCIUM"]   = 103.0;
    element_to_Z["UNQ"]          = 104.0;
    element_to_Z["UNP"]          = 105.0;
    element_to_Z["UNH"]          = 106.0;
    element_to_Z["UNS"]          = 107.0;
  }

};

}

#endif
