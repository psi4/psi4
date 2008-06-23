/*! \file read_calculation_options
    \defgroup PSI4
*/

#include <libipv1/ip_lib.h>
#include <liboptions/liboptions.h>

namespace psi {

  int read_options(std::string name, Options & options) {

    ip_cwk_clear();
    ip_cwk_add(":BASIS");
    ip_cwk_add(":DEFAULT");
    ip_cwk_add(":PSI");
    ip_set_uppercase(1);

    if (name == "INPUT") {
      ip_cwk_add(":INPUT");
      options->add_bool_option("NO_REORIENT",false);
      options->add_bool_option("CHKPT_MOS",false);
      options->add_str("LABEL","Default PSI3 Label");
      options->add_bool_option("SHOWNORM",false);
      options->add_bool_option("NORMALIZE",true);
      options->add_bool_option("PUREAM",false);
      options->add_bool_option("EXPERT",false);
      options->add_int_option("PRINT",1);
      options->add_str_option_with_choices("SUBGROUP", NULL, "C1 C2 CS CI C2V C2H D2")
      options->add_str_option_with_choices("UNIQUE_AXIS", NULL, "X Y Z")
      options->add_int_option("NFRAGMENTS",1);
      options->add_bool_option("KEEP_REF_FRAME",false);
      options->add_bool_option("FREEZE_CORE",false);
      options->add_bool_option("FREEZE_VIRT",false);
      options->add_str_option_with_choices("UNITS", "ANGSTROMS", "BOHR AU ANGSTROMS ANGSTROM")
      options->add_bool_option("FRAGMENT_DISTANCE_INVERSE",false);
      // input also does ip calls to read basis
    }
    options->read_options();
    options->print();
  }

} //end ::psi

