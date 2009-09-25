//  MCSCF (2008) by Francesco Evangelista - frank@ccc.uga.edu
//  Center for Computational Chemistry, University of Georgia

/**
 *  @defgroup MCSCF MCSCF is a code for SCF/MCSCF computations
 *  @ingroup MCSCF
 *  @brief Contains main() and global variables
*/

// Standard libraries
#include <iostream>
#include <cstdlib>

// PSI libraries
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>
#include <libutil/libutil.h>

// PSI C++ libraries
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>

#include "mcscf.h"
#include "git.h"
#include "scf.h"

// PSI FILES
FILE  *infile, *outfile;
char  *psi_file_prefix;
extern "C" {
  const char* gprgid();
}

namespace psi{
  MOInfoSCF              *moinfo_scf;
  MOInfo                 *moinfo;
  MemoryManager          *_memory_manager_;
  namespace mcscf{
    void add_calculation_options();
  }
} /* End Namespaces */

using namespace std;

/**
 * The main function
 * @param argc
 * @param argv[]
 * @return PSI_RETURN_SUCCESS if the program ran without any problem
 */
int main(int argc, char *argv[])
{
  using namespace psi;
  using namespace psi::mcscf;

  init_psi(argc,argv);



  if(options_get_str("REFERENCE") == "RHF"  ||
     options_get_str("REFERENCE") == "ROHF" ||
     options_get_str("REFERENCE") == "UHF"  ||
     options_get_str("REFERENCE") == "TWOCON"){
    SCF scf;
    scf.compute_energy();
  }else if(options_get_str("REFERENCE") == "MCSCF"){
    fprintf(outfile,"\n\nREFERENCE = MCSCF not implemented yet");
    fflush(outfile);
    return PSI_RETURN_FAILURE;
  }

  close_psi();
  return PSI_RETURN_SUCCESS;
}


namespace psi{ namespace mcscf{

void add_calculation_options()
{

  options_add_int("CONVERGENCE",9);
  options_add_int("LEVELSHIFT",0);
//  options_add_int("DUMPING",0);
  options_add_int("DEBUG",0);
  options_add_int("MAXITER",100);
  options_add_int("NDIIS",7);
  options_add_int("ROOT",1);
  options_add_int("START_FAVG",5);
  options_add_int("TURN_ON_ACTV",0);
  options_add_int("ROTATE_MO_ANGLE",0);
  options_add_int("ROTATE_MO_IRREP",1);  // IRREP is one-based
  options_add_int("ROTATE_MO_P",1);      // P and Q are one-based
  options_add_int("ROTATE_MO_Q",2);

  options_add_bool("CI_DIIS",false);
  options_add_bool("USE_DIIS",true);
  options_add_bool("READ_MOS",true);
  options_add_bool("USE_FAVG",false);
  options_add_bool("CANONICALIZE_ACTIVE_FAVG",false);
  options_add_bool("CANONICALIZE_INACTIVE_FAVG",false);
  options_add_bool("INTERNAL_ROTATIONS",true);
  options_add_bool("FORCE_TWOCON",false);

  options_add_str_with_choices("REFERENCE","RHF","RHF ROHF UHF TWOCON MCSCF GENERAL");
  options_add_str_with_choices("WFN_SYM","1","A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
}

/**
 * Start Psi3, draw the program logo, version, and compilation details
 * @param argc
 * @param argv[]
 */
void init_psi(int argc, char *argv[])
{
  int num_extra_args=0;
  char**  extra_args;

  extra_args = new char*[argc];

  for(int i=1; i<argc; i++) {
    extra_args[num_extra_args++] = argv[i];
/*  Template for argument parsing
    if(strcmp(argv[i], "--opdm") == 0) {
    }
    else {
      extra_args[num_extra_args++] = argv[i];
    }
*/
  }
  psi_start(&infile,&outfile,&psi_file_prefix,num_extra_args,extra_args,0);
  delete[] extra_args;

  _default_psio_lib_ = new PSIO();
  psiopp_ipv1_config(_default_psio_lib_);
  _default_chkpt_lib_ = new Chkpt(_default_psio_lib_, PSIO_OPEN_OLD);

  ip_cwk_clear();
  ip_cwk_add(const_cast<char*>(":PSI"));
  ip_cwk_add(const_cast<char*>(":SCF"));
  ip_cwk_add(const_cast<char*>(":MCSCF"));

  tstart(outfile);

  fprintf(outfile,"\n         ------------------------------------------");
  fprintf(outfile,"\n           MCSCF: a self-consistent field program");
  fprintf(outfile,"\n            written by Francesco A. Evangelista");
  fprintf(outfile,"\n         ------------------------------------------");

  fprintf(outfile,"\n\n\n  Compiled on %s at %s",__DATE__,__TIME__);

  options_init();
  add_calculation_options();
  options_read();
  options_print();

  _default_psio_lib_->open(PSIF_MCSCF,PSIO_OPEN_NEW);

  psi::_memory_manager_   = new MemoryManager();
  moinfo_scf = new MOInfoSCF();
}

/**
 * Close psi by calling psio_done() and psi_stop()
 */
void close_psi()
{
  if(options_get_int("DEBUG") > 0)
    psi::_memory_manager_->MemCheck(outfile);
  delete moinfo_scf;
  delete _memory_manager_;

  fprintf(outfile,"\n\n  MCSCF Execution Completed.\n\n");
  fflush(outfile);

  options_close();

  _default_psio_lib_->close(PSIF_MCSCF,1);

  delete _default_chkpt_lib_;
  delete _default_psio_lib_;

  tstop(outfile);
  psi_stop(infile,outfile,psi_file_prefix);
}

}} /* End Namespaces */

/**
 * @return program ID
 */
const char* gprgid()
{
  return(const_cast<char*>("MCSCF"));
}
