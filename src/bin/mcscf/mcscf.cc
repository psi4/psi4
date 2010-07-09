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

namespace psi{
MemoryManager* memory_manager;
MOInfoSCF*     moinfo_scf;

namespace mcscf{

using namespace std;

/**
 * The main function
 * @param argc
 * @param argv[]
 * @return PSI_RETURN_SUCCESS if the program ran without any problem
 */
PsiReturnType mcscf(Options& options)
{
  using namespace psi;
  shared_ptr<PSIO> psio(new PSIO);
  psiopp_ipv1_config(psio);
  shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

  psio->open(PSIF_MCSCF,PSIO_OPEN_NEW);
  init_psi(options);

  memory_manager  = new psi::MemoryManager();
  moinfo_scf      = new psi::MOInfoSCF(options,chkpt);

  if(options.get_str("REFERENCE") == "RHF"  ||
     options.get_str("REFERENCE") == "ROHF" ||
     options.get_str("REFERENCE") == "UHF"  ||
     options.get_str("REFERENCE") == "TWOCON"){
    SCF scf(options,psio,chkpt);
    scf.compute_energy();
  }else if(options.get_str("REFERENCE") == "MCSCF"){
    fprintf(outfile,"\n\nREFERENCE = MCSCF not implemented yet");
    fflush(outfile);
    return Failure;
  }
  delete moinfo_scf;
  delete memory_manager;
  close_psi(options);
  psio->close(PSIF_MCSCF,1);
  return Success;
}

/**
 * Start Psi3, draw the program logo, version, and compilation details
 * @param argc
 * @param argv[]
 */
void init_psi(Options& options_)
{
  fprintf(outfile,"\n         ------------------------------------------");
  fprintf(outfile,"\n           MCSCF: a self-consistent field program");
  fprintf(outfile,"\n            written by Francesco A. Evangelista");
  fprintf(outfile,"\n         ------------------------------------------");
}

/**
 * Close psi by calling psio_done() and psi_stop()
 */
void close_psi(Options& options_)
{
  fprintf(outfile,"\n\n  MCSCF Execution Completed.\n\n");
  fflush(outfile);
  tstop();
}

}} /* End Namespaces */
