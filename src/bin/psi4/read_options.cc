/*! \file read_calculation_options
    \defgroup PSI4
*/

#include <libipv1/ip_lib.h>
#include <liboptions/liboptions.hpp>
#include <psifiles.h>

namespace psi {

int read_options(std::string name, Options & options) {

  ip_cwk_clear();
  ip_cwk_add(":BASIS");
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":PSI");
  ip_set_uppercase(1);

  if (name == "PSI4") {
    options.add_str_option_with_choices("UNITS", "ANGSTROMS", "BOHR AU ANGSTROMS ANGSTROM");
  }
  else if (name == "INPUT") {
    ip_cwk_add(":INPUT");

    /* Keep the chkpt file. */
    options.add_bool_option("KEEP_CHKPT",false);

    /* read MOs from checkpoint file and project onto new basis */
    options.add_bool_option("CHKPT_MOS",false);

  // these may be obseleted in psi4
    /*--- read geometry from checkpoint file (in findif calculations) ---*/
    options.add_bool_option("CHKPT_GEOM",false);

    /*--- don't project MOs but simply keep them ---*/
    options.add_bool_option("NOPROJECT",false);

    /*--- read geometry from geom.dat file (in findif calculations) ---*/
    options.add_bool_option("GEOMDAT",false);

    /* No center of mass shift */
    options.add_bool_option("NO_COM_SHIFT",false);

    /*--- read MOs from checkpoint file and save to a separate file ---*/
    options.add_bool_option("SAVE_MOS",false);

    /* Don't overwrite the output file. */
    options.add_bool_option("KEEP_OUTPUT",false);
  
    options.add_str_option("WFN", NULL);

    options.add_bool_option("NO_REORIENT",false);
    options.add_str_option("LABEL","Default PSI3 Label");
    options.add_bool_option("SHOWNORM",false);
    options.add_bool_option("NORMALIZE",true);
    options.add_bool_option("PUREAM",false);
    options.add_bool_option("EXPERT",false);
    options.add_int_option("PRINT",1);
    options.add_str_option_with_choices("SUBGROUP", NULL, "C1 C2 CS CI C2V C2H D2");
    options.add_str_option_with_choices("UNIQUE_AXIS", NULL, "X Y Z");
    options.add_int_option("NFRAGMENTS",1);
    options.add_bool_option("KEEP_REF_FRAME",false);
    options.add_bool_option("FRAGMENT_DISTANCE_INVERSE",false);
    // input also does ip calls to read basis, geometry
    
    options.add_str_option_with_choices("FREEZE_CORE","FALSE", \
      "FALSE NO TRUE YES SMALL LARGE");
    options.add_int_option("FREEZE_VIRT",0);
  }
  else if (name == "CINTS") {
    ip_cwk_add(":CINTS");
    options.add_str_option("WFN", NULL);
    options.add_int_option("PRINT",1);
    options.add_int_option("NUM_THREADS",1);
    options.add_int_option("CUTOFF",15); // cutoff on integrals
    // scaling for counterfactual fine-structure constant
    options.add_double_option("FINE_STRUCTURE_ALPHA", 1.0);

    options.add_bool_option("MAKE_ERI", 1);
    options.add_bool_option("EMPIRICAL_DISPERSION", 0);
    options.add_bool_option("RESTART", 0);
    options.add_int_option("RESTART_TASK", 0);

    options.add_int_option("S_FILE", PSIF_OEI);
    options.add_int_option("T_FILE", PSIF_OEI);
    options.add_int_option("V_FILE", PSIF_OEI);
    options.add_int_option("ERI_FILE", PSIF_SO_TEI);

  }

  options.read_options();
  options.print(name);
 }

} //end ::psi

