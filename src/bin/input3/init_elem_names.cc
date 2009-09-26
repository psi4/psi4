/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <cstring>
#include <libchkpt/chkpt.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

void init_elem_names()
{
  elem_name[0] = strdup("GHOST");
  elem_name[1] = strdup("HYDROGEN"); 
  elem_name[2] = strdup("HELIUM");
  elem_name[3] = strdup("LITHIUM");
  elem_name[4] = strdup("BERYLLIUM");
  elem_name[5] = strdup("BORON");
  elem_name[6] = strdup("CARBON");
  elem_name[7] = strdup("NITROGEN");
  elem_name[8] = strdup("OXYGEN");
  elem_name[9] = strdup("FLUORINE");
  elem_name[10] = strdup("NEON");
  elem_name[11] = strdup("SODIUM");
  elem_name[12] = strdup("MAGNESIUM");
  elem_name[13] = strdup("ALUMINUM");
  elem_name[14] = strdup("SILICON");
  elem_name[15] = strdup("PHOSPHORUS");
  elem_name[16] = strdup("SULFUR");
  elem_name[17] = strdup("CHLORINE");
  elem_name[18] = strdup("ARGON");
  elem_name[19] = strdup("POTASSIUM");
  elem_name[20] = strdup("CALCIUM");
  elem_name[21] = strdup("SCANDIUM");
  elem_name[22] = strdup("TITANIUM");
  elem_name[23] = strdup("VANADIUM");
  elem_name[24] = strdup("CHROMIUM");
  elem_name[25] = strdup("MANGANESE");
  elem_name[26] = strdup("IRON");
  elem_name[27] = strdup("COBALT");
  elem_name[28] = strdup("NICKEL");
  elem_name[29] = strdup("COPPER");
  elem_name[30] = strdup("ZINC");
  elem_name[31] = strdup("GALLIUM");
  elem_name[32] = strdup("GERMANIUM");
  elem_name[33] = strdup("ARSENIC");
  elem_name[34] = strdup("SELENIUM");
  elem_name[35] = strdup("BROMINE");
  elem_name[36] = strdup("KRYPTON");
  elem_name[37] = strdup("RUBIDIUM");
  elem_name[38] = strdup("STRONTIUM");
  elem_name[39] = strdup("YTTRIUM");
  elem_name[40] = strdup("ZIRCONIUM");
  elem_name[41] = strdup("NIOBIUM");
  elem_name[42] = strdup("MOLYBDENUM");
  elem_name[43] = strdup("TECHNETIUM");
  elem_name[44] = strdup("RUTHENIUM");
  elem_name[45] = strdup("RHODIUM");
  elem_name[46] = strdup("PALLADIUM");
  elem_name[47] = strdup("SILVER");
  elem_name[48] = strdup("CADMIUM");
  elem_name[49] = strdup("INDIUM");
  elem_name[50] = strdup("TIN");
  elem_name[51] = strdup("ANTIMONY");
  elem_name[52] = strdup("TELLURIUM");
  elem_name[53] = strdup("IODINE");
  elem_name[54] = strdup("XENON");
  elem_name[55] = strdup("CESIUM");
  elem_name[56] = strdup("BARIUM");
  elem_name[57] = strdup("LANTHANUM");
  elem_name[58] = strdup("CERIUM");
  elem_name[59] = strdup("PRASEODYMIUM");
  elem_name[60] = strdup("NEODYMIUM");
  elem_name[61] = strdup("PROMETHIUM");
  elem_name[62] = strdup("SAMARIUM");
  elem_name[63] = strdup("EUROPIUM");
  elem_name[64] = strdup("GADOLINIUM");
  elem_name[65] = strdup("TERBIUM");
  elem_name[66] = strdup("DYSPROSIUM");
  elem_name[67] = strdup("HOLMIUM");
  elem_name[68] = strdup("ERBIUM");
  elem_name[69] = strdup("THULIUM");
  elem_name[70] = strdup("YTTERBIUM");
  elem_name[71] = strdup("LUTETIUM");
  elem_name[72] = strdup("HAFNIUM");
  elem_name[73] = strdup("TANTALUM");
  elem_name[74] = strdup("TUNGSTEN");
  elem_name[75] = strdup("RHENIUM");
  elem_name[76] = strdup("OSMIUM");
  elem_name[77] = strdup("IRIDIUM");
  elem_name[78] = strdup("PLATINUM");
  elem_name[79] = strdup("GOLD");
  elem_name[80] = strdup("MERCURY");
  elem_name[81] = strdup("THALLIUM");
  elem_name[82] = strdup("LEAD");
  elem_name[83] = strdup("BISMUTH");
  elem_name[84] = strdup("POLONIUM");
  elem_name[85] = strdup("ASTATINE");
  elem_name[86] = strdup("RADON");
  elem_name[87] = strdup("FRANCIUM");
  elem_name[88] = strdup("RADIUM");
  elem_name[89] = strdup("ACTINIUM");
  elem_name[90] = strdup("THORIUM");
  elem_name[91] = strdup("PROTACTINIUM");
  elem_name[92] = strdup("URANIUM");
  elem_name[93] = strdup("NEPTUNIUM");
  elem_name[94] = strdup("PLUTONIUM");
  elem_name[95] = strdup("AMERICIUM");
  elem_name[96] = strdup("CURIUM");
  elem_name[97] = strdup("BERKELIUM");
  elem_name[98] = strdup("CALIFORNIUM");
  elem_name[99] = strdup("EINSTEINIUM");
  elem_name[100] = strdup("FERMIUM");
  elem_name[101] = strdup("MENDELEVIUM");
  elem_name[102] = strdup("NOBELIUM");
  elem_name[103] = strdup("LAWRENCIUM");
  elem_name[104] = strdup("UNQ");
  elem_name[105] = strdup("UNP");
  elem_name[106] = strdup("UNH");
  elem_name[107] = strdup("UNS");
}

}} // namespace psi::input
