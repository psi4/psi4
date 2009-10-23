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
  namespace ccsort   { PsiReturnType ccsort(Options &, int argc, char *argv[]); }
  namespace transqt2 { PsiReturnType transqt2(Options &, int argc, char *argv[]); }
  namespace ccsort   { PsiReturnType ccsort(Options &, int argc, char *argv[]); }
  namespace ccenergy { PsiReturnType ccenergy(Options &, int argc, char *argv[]); }
  namespace scf      { PsiReturnType scf(Options&, int, char**); }

  int read_options(std::string name, Options & options);

void launch_module(std::string modname, 
  Options &options, int argc, char* argv[], 
  std::map<std::string, PsiReturnType(*)(Options &, int argc, char *argv[])> 
  dispatch_table);

void execute_sequence(std::string sequence[],
  Options &options, int argc, char* argv[], 
  std::map<std::string, PsiReturnType(*)(Options &, int argc, char *argv[])> 
  dispatch_table);

void psiclean(void);

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

  fprintf(outfile,"\tPSI3 simulator method name: %s\n", PsiMethod.c_str());

  int ndisp = 1;
  int nopt = 40;
  int ncasiter = 20;
  int nbrueckner = 20;

  // make a map of function pointers to the functions
  std::map<std::string, PsiReturnType(*)(Options &, int argc, char *argv[])> dispatch_table;

  dispatch_table["INPUT"]     = &(psi::input::input);
  dispatch_table["SCF"]       = &(psi::scf::scf);
  dispatch_table["CSCF"]      = &(psi::cscf::cscf);
  dispatch_table["CINTS"]     = &(psi::CINTS::cints);
  dispatch_table["TRANSQT2"]  = &(psi::transqt2::transqt2);
  dispatch_table["CCSORT"]    = &(psi::ccsort::ccsort);
  dispatch_table["CCENERGY"]  = &(psi::ccenergy::ccenergy);
//  dispatch_table["CCTRIPLES"] = &(psi::CCTRIPLES::CCTRIPLES);

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

  string seq_scf_energy[]    = {"CINTS", "CSCF", "END"};
  string seq_ccsd_energy[]   = {"CINTS", "CSCF", "TRANSQT2", "CCSORT",
                                "CCENERGY", "END"};
  string seq_ccsd_t_energy[] = {"CINTS", "CSCF", "TRANSQT2", "CCSORT",
                                "CCENERGY", "CCTRIPLES", "END"};


  launch_module("INPUT", options, argc, argv, dispatch_table);

  if (PsiMethod == "SCFENERGY") {
    execute_sequence(seq_scf_energy, options, argc, argv, dispatch_table);
  }
  else if (PsiMethod == "CCSDENERGY") {
    execute_sequence(seq_ccsd_energy, options, argc, argv, dispatch_table);
  }
  else if (PsiMethod == "CCSD_TENERGY") {
    fprintf(outfile, "EXECUTING CCSD_T ENERGY 123\n");
    execute_sequence(seq_ccsd_t_energy, options, argc, argv, dispatch_table);
  }
  else {
    fprintf(outfile,"Unknown PSI4 method\n");
    abort();
  }
 
  psiclean();
  return 1;
} 

void execute_sequence(std::string sequence[],
  Options &options, int argc, char* argv[], 
  std::map<std::string, PsiReturnType(*)(Options &, int argc, char *argv[])> 
  dispatch_table) 
{
  int i=0; 
  while (sequence[i] != "END") {
    launch_module(sequence[i], options, argc, argv, dispatch_table);
    i++;
  } 
}

void launch_module(std::string modname, 
  Options &options, int argc, char* argv[], 
  std::map<std::string, PsiReturnType(*)(Options &, int argc, char *argv[])> 
  dispatch_table) 
{
  module.set_prgid(modname);
  read_options(modname, options); 
  dispatch_table[modname](options, argc, argv);
}

void psiclean(void)
{
  ULI i, nvol;
  int errcod;
  char *vpath;
  char *basename;
  const int MAX_STRING = 300;
  char fileslist[MAX_STRING];
  char cmdstring[MAX_STRING];

  nvol = psio_get_numvols_default();

  errcod = psio_get_filename_default(&basename);

  for (i=0; i<nvol; i++) {
    errcod = psio_get_volpath_default(i, &vpath);

    sprintf(fileslist,"%s%s.*",vpath,basename);
    sprintf(cmdstring,"echo Removing files %s%s",vpath,basename);
    system(cmdstring);
    sprintf(cmdstring,"ls -l %s",fileslist);
    system(cmdstring);
    sprintf(cmdstring,"/bin/rm %s",fileslist);
    system(cmdstring);
    free(vpath);
  }
  free(basename);
}

}
