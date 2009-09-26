/*! \file psi3_simulator.cc
*/

//  calls functions according to old PSI3 psi.dat file
//  Rollin King

#include <stdio.h>

#include <psi4-dec.h>
#include "psi4.h"

namespace psi {

  namespace input    { PsiReturnType input(Options &, int argc, char *argv[]); }
  namespace CINTS    { PsiReturnType cints(Options &, int argc, char *argv[]); }
  namespace cscf     { PsiReturnType cscf(Options &, int argc, char *argv[]); }
  namespace transqt2 { PsiReturnType transqt2(Options &, int argc, char *argv[]); }
  namespace psiclean { PsiReturnType psiclean(Options &, int argc, char *argv[]); }

  int read_options(std::string name, Options & options);

// constructs psi.dat form string for execution
int psi3_simulator(Options & options, int argc, char *argv[]) {
  std::string jobtype, wfn, ref, dertype, PsiMethod;

  options.add_str("JOBTYPE", "SP", "SP OPT");
  options.add_str("WFN", "SCF", "SCF MP2 CCSD CCSD_T");
  options.add_str("REF", "RHF", "RHF UHF ROHF TCSCF");
  options.add_str("DERTYPE", "ENERGY", "NONE ENERGY FIRST SECOND");
  options.read_ipv1();

  jobtype = options.get_str("JOBTYPE");
  wfn = options.get_str("WFN");
  ref = options.get_str("REF");
  dertype = options.get_str("DERTYPE");

  if (dertype == "NONE")
    dertype = "ENERGY";

  PsiMethod = wfn + dertype;

  fprintf(outfile,"PSI3 simulator method name: %s\n", PsiMethod.c_str());

  int ndisp = 1;
  int nopt = 40;
  int ncasiter = 20;
  int nbrueckner = 20;

  // make a map of function pointers to the functions
  std::map<std::string, PsiReturnType(*)(Options &, int argc, char *argv[])> dispatch_table;

  dispatch_table["INPUT"] = &(psi::input::input);
  dispatch_table["CSCF"]   = &(psi::cscf::cscf);
  dispatch_table["CINTS"] = &(psi::CINTS::cints);
  dispatch_table["TRANSQT2"] = &(psi::transqt2::transqt2);

  /* the basic programs that were in psi3
  input      = "input"
  uinput     = "input --chkptgeom"
  done       = ("psiclean")
  ints       = "cints"
  scf        = "cscf"
  localize   = "localize"
  deriv      = "cints --deriv1"
  deriv2     = "cints --deriv2"
  propint    = "cints --oeprop"
  mkpt2ints  = "cints --mkpt2"
  transqt    = "transqt2"
  %transqt    = "transqt2"
  cctrans    = "transqt2"  %testing with CC modules; eventually will replace transqt
  backtransqt = "transqt --backtr"
  cphf        = "cphf"
  cphf_X      = "cphf --X_only"
  response    = "response"
  normco      = "normco"
  oeprop      = "oeprop"
  optking     = "optking"
  stable      = "stable"
  cis         = "cis"
*/

  fprintf(outfile,"quitting now\n");

  if (PsiMethod == "SCFENERGY") {

      module.set_prgid("INPUT");
      read_options("INPUT", options);
      dispatch_table["INPUT"](options, argc, argv);

      module.set_prgid("CINTS");
      read_options("CINTS", options);
      dispatch_table["CINTS"](options, argc, argv);

      module.set_prgid("CSCF");
      read_options("CSCF", options);
      dispatch_table["CSCF"](options, argc, argv);
  }
  else {
    fprintf(outfile,"Unknown PSI4 method\n");
    abort();
  }

  return 1;
} 

}
