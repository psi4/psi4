/*! \file read_calculation_options
    \defgroup PSI4
*/

#include <libipv1/ip_lib.h>
#include <liboptions/liboptions.h>
#include <psifiles.h>

namespace psi {

int read_options(std::string name, Options & options) {

  ip_cwk_clear();
  ip_cwk_add(":BASIS");
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":PSI");
  ip_set_uppercase(1);

  if (name == "PSI4") {
    options.add_str("UNITS", "ANGSTROMS", "BOHR AU ANGSTROMS ANGSTROM");
  }
  else if (name == "INPUT") {
    ip_cwk_add(":INPUT");

    /* Keep the chkpt file. */
    options.add_bool("KEEP_CHKPT",false);

    /* read MOs from checkpoint file and project onto new basis */
    options.add_bool("CHKPT_MOS",false);

  // these may be obseleted in psi4
    /*--- read geometry from checkpoint file (in findif calculations) ---*/
    options.add_bool("CHKPT_GEOM",false);

    /*--- don't project MOs but simply keep them ---*/
    options.add_bool("NOPROJECT",false);

    /*--- read geometry from geom.dat file (in findif calculations) ---*/
    options.add_bool("GEOMDAT",false);

    /* No center of mass shift */
    options.add_bool("NO_COM_SHIFT",false);

    /*--- read MOs from checkpoint file and save to a separate file ---*/
    options.add_bool("SAVE_MOS",false);

    /* Don't overwrite the output file. */
    options.add_bool("KEEP_OUTPUT",false);
  
    options.add_str("WFN", "");

    options.add_bool("NO_REORIENT",false);
    options.add_str("LABEL","Default PSI3 Label");
    options.add_bool("SHOWNORM",false);
    options.add_bool("NORMALIZE",true);
    options.add_bool("PUREAM",false);
    options.add_bool("EXPERT",false);
    options.add_int("PRINT",1);
    options.add_str("SUBGROUP", "", "C1 C2 CS CI C2V C2H D2");
    options.add_str("UNIQUE_AXIS", "", "X Y Z");
    options.add_int("NFRAGMENTS",1);
    options.add_bool("KEEP_REF_FRAME",false);
    options.add_bool("FRAGMENT_DISTANCE_INVERSE",false);
    // input also does ip calls to read basis, geometry
    
    options.add_str("FREEZE_CORE","FALSE", \
      "FALSE NO TRUE YES SMALL LARGE");
    options.add_int("FREEZE_VIRT",0);
  }
  else if (name == "CINTS") {
    ip_cwk_add(":CINTS");
    options.add_str("WFN", "");
    options.add_int("PRINT",1);
    options.add_int("NUM_THREADS",1);
    options.add_int("CUTOFF",15); // cutoff on integrals
    // scaling for counterfactual fine-structure constant
    options.add_double("FINE_STRUCTURE_ALPHA", 1.0);

    options.add_bool("MAKE_ERI", 1);
    options.add_bool("EMPIRICAL_DISPERSION", 0);
    options.add_bool("RESTART", 0);
    options.add_int("RESTART_TASK", 0);

    options.add_int("S_FILE", PSIF_OEI);
    options.add_int("T_FILE", PSIF_OEI);
    options.add_int("V_FILE", PSIF_OEI);
    options.add_int("ERI_FILE", PSIF_SO_TEI);

  }
  else if (name == "MP2") {
    options.add_str("WFN", "");
    options.add_str("REFERENECE", "RHF");
    options.add_str("JOBTYPE", "SP");
    options.add_str("DERTYPE", "NONE");
    options.add_int("PRINT", 0);
    options.add_int("CACHELEV", 2);
    options.add_int("CACHETYPE", 1);
    options.add_bool("SCS","false");
    options.add_bool("SCS_N", "false");
    options.add_double("SCALE_OS", 6.0/5.0);
    options.add_double("SCALE_SS", 1.0/3.0);
  }
  else if(name == "TRANSQT") {
    options.add_str("WFN", "");
    options.add_str("REFERENCE", "RHF");
    options.add_str("DERTYPE", "NONE");
    options.add_int("PRINT", 1);
    options.add_bool("PRINT_TEI", true);
    options.add_int("TOLERANCE", 14);
    options.add_int("CACHELEV", 2);
    options.add_str("AO_BASIS", "NONE");
    options.add_bool("DELETE_TEI", true);
  }
  else if(name == "CUSP"){
    options.add("FROZEN_DOCC", new ArrayType());
    options.add("FROZEN_UOCC", new ArrayType());
  }   
  else if(name == "CCSORT") {
    options.add_str("WFN", "");
    options.add_str("REFERENCE", "RHF");
    options.add_str("DERTYPE", "NONE");
    options.add_str("PROPERTY", "POLARIZABILITY");
    options.add_bool("LOCAL", false);
    options.add_double("LOCAL_CUTOFF", 0.02);
    options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
    options.add_double("LOCAL_CORE_CUTOFF",0.05);
    options.add_cstr("LOCAL_METHOD","WERNER");
    options.add_cstr("LOCAL_WEAKP");
    options.add_cstr("FREEZE_CORE");
    options.add_cstr("LOCAL_PAIRDEF");
    options.add_bool("LOCAL_DOMAIN_POLAR");
    options.add_bool("LOCAL_DOMAIN_MAG");
    options.add_bool("LOCAL_DOMAIN_SEP");
    options.add_bool("LOCAL_FILTER_SINGLES");
    options.add_cstr("AO_BASIS");
    options.add_cstr("EOM_REFERENCE");
    options.add_int("PRINT");
    options.add_bool("KEEP_TEIFILE");
    options.add_bool("KEEP_OEIFILE");
    options.add_int("TOLERANCE");
    options.add_int("CACHELEV");
    options.add_bool("LOCAL");
    options.add("OMEGA",new ArrayType());
  }
  else if(name == "CCTRIPLES") {
    options.add_str("WFN", "");
    options.add_int("NTHREADS",1);
    options.add_str("REFERENCE","RHF");
    options.add_str("DERTYPE","NONE");
  }
  else if(name == "CLAG") {
    options.add_bool("WRITE_CAS_FILES",0);
    options.add_str("DERTYPE","NONE");
    options.add_str("WFN","NONE");
    options.add_int("ROOT",1);
  }
  options.read_ipv1();
  options.print();
 }

} //end ::psi

