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

namespace psi{  namespace MCSCF{
void add_calculation_options();

using namespace std;

/**
 * The main function
 * @param argc
 * @param argv[]
 * @return PSI_RETURN_SUCCESS if the program ran without any problem
 */
int mcscf(Options& options,int argc, char *argv[])
{
  using namespace psi;
  init_psi(argc,argv);

  if(options.get_str("REFERENCE") == "RHF"  ||
     options.get_str("REFERENCE") == "ROHF" ||
     options.get_str("REFERENCE") == "UHF"  ||
     options.get_str("REFERENCE") == "TWOCON"){
    SCF scf(options);
    scf.compute_energy();
  }else if(options.get_str("REFERENCE") == "MCSCF"){
    fprintf(outfile,"\n\nREFERENCE = MCSCF not implemented yet");
    fflush(outfile);
    return PSI_RETURN_FAILURE;
  }

  close_psi();
  return PSI_RETURN_SUCCESS;
}

void add_calculation_options()
{
}

/**
 * Start Psi3, draw the program logo, version, and compilation details
 * @param argc
 * @param argv[]
 */
void init_psi(Options& options)
{
  tstart();

  fprintf(outfile,"\n         ------------------------------------------");
  fprintf(outfile,"\n           MCSCF: a self-consistent field program");
  fprintf(outfile,"\n            written by Francesco A. Evangelista");
  fprintf(outfile,"\n         ------------------------------------------");

  fprintf(outfile,"\n\n\n  Compiled on %s at %s",__DATE__,__TIME__);

  add_calculation_options();

  _default_psio_lib_->open(PSIF_MCSCF,PSIO_OPEN_NEW);

  psi::_memory_manager_   = new MemoryManager();
  moinfo_scf = new MOInfoSCF(options);
}

/**
 * Close psi by calling psio_done() and psi_stop()
 */
void close_psi(Options& options)
{
  if(options.get_int("DEBUG") > 0)
    psi::_memory_manager_->MemCheck(outfile);
  delete moinfo_scf;
  delete _memory_manager_;

  fprintf(outfile,"\n\n  MCSCF Execution Completed.\n\n");
  fflush(outfile);

  _default_psio_lib_->close(PSIF_MCSCF,1);

  delete _default_chkpt_lib_;
  delete _default_psio_lib_;

  tstop();
}

}} /* End Namespaces */
