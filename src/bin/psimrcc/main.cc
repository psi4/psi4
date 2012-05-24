/*******************************************************************
 *  PSIMRCC (2008) by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ******************************************************************/

/**
 *  @defgroup PSIMRCC PSIMRCC is a code for SR/MRCC computations
 *  @file psimrcc.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Contains main() and global variables
*/

#include <psi4-dec.h>
// Standard libraries
#include <iostream>
#include <complex>
#include <cstdlib>

// PSI libraries
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>
#include <libutil/libutil.h>
#include <libutil/memory_manager.h>
#include <libqt/qt.h>

// PSI C++
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>

#include "blas.h"
#include "git.h"
#include "main.h"
#include "sort.h"
#include "mrcc.h"
#include "debugging.h"
#include "psimrcc.h"
#include "transform.h"

// PSI FILES

using namespace std;
using namespace psi;

namespace psi{
    extern FILE  *outfile;
    namespace psimrcc{
    // Global variables
    Timer               *global_timer;
    CCBLAS              *blas;
    CCSort              *sorter;
    CCTransform         *trans = NULL;
    MOInfo              *moinfo;
    ModelSpace          *model_space;
    Debugging           *debugging;
    MemoryManager       *memory_manager;

PsiReturnType
psimrcc(Options &options)
{
  using namespace psi::psimrcc;
  _default_psio_lib_->open(PSIF_PSIMRCC_INTEGRALS,PSIO_OPEN_NEW);

  fprintf(outfile,"\n  MRCC          MRCC");
  fprintf(outfile,"\n   MRCC  MRCC  MRCC");
  fprintf(outfile,"\n   MRCC  MRCC  MRCC      PSIMRCC Version 0.9.3.3, July 2009");
  fprintf(outfile,"\n   MRCC  MRCC  MRCC      Multireference Coupled Cluster, written by");
  fprintf(outfile,"\n     MRCCMRCCMRCC        Francesco A. Evangelista and Andrew C. Simmonett");
  fprintf(outfile,"\n         MRCC            Compiled on %s at %s",__DATE__,__TIME__);
  fprintf(outfile,"\n         MRCC            id =%s",GIT_ID);
  fprintf(outfile,"\n       MRCCMRCC");

  global_timer = new Timer;
  debugging = new Debugging(options);
  moinfo = new MOInfo(options);

  memory_manager = new MemoryManager(Process::environment.get_memory());
  model_space = new ModelSpace(moinfo);

  moinfo->setup_model_space();  // The is a bug here DELETEME

  int nactmo = moinfo->get_nactv();
  int nactel = moinfo->get_nactive_ael() + moinfo->get_nactive_bel();
  if(nactel > 2 and nactmo > 2){
      fprintf(outfile,"\n   WARNING: PSIMRCC detected that you are not using a CAS(2,n) or CAS(m,2) active space");
      fprintf(outfile,"\n            You requested a CAS(%d,%d) space.  In this case the program will run",nactel,nactmo);
      fprintf(outfile,"\n            but will negled matrix elements of the effective Hamiltonian between");
      fprintf(outfile,"\n            reference determinats that differ by more than two spin orbitals.");
      fprintf(outfile,"\n            The final answer will NOT be the Mk-MRCC energy but only an approximation to it.");
      fprintf(outfile,"\n            If you are going to report this number in a publication make sure that you");
      fprintf(outfile,"\n            understand what is going on and that you document it in your publication.");
  }

  blas   = new CCBLAS(options);
  trans  = new CCTransform();
  if(options.get_str("CORR_WFN")=="PT2"){
    mrpt2(options);
  }else{
    mrccsd(options);
    if(nactel > 2 and nactmo > 2){
        fprintf(outfile,"\n   WARNING: PSIMRCC detected that you are not using a CAS(2,n) or CAS(m,2) active space");
        fprintf(outfile,"\n            You requested a CAS(%d,%d) space.  In this case the program will run",nactel,nactmo);
        fprintf(outfile,"\n            but will negled matrix elements of the effective Hamiltonian between");
        fprintf(outfile,"\n            reference determinats that differ by more than two spin orbitals.");
        fprintf(outfile,"\n            The final answer will NOT be the Mk-MRCC energy but only an approximation to it.");
        fprintf(outfile,"\n            If you are going to report this number in a publication make sure that you");
        fprintf(outfile,"\n            understand what is going on and that you document it in your publication.");
    }
  }

  delete sorter;
  delete trans;
  delete blas;

  fprintf(outfile,"\n\n  PSIMRCC job completed.");
  fprintf(outfile,"\n  Wall Time = %20.6f s",global_timer->get());
  fprintf(outfile,"\n  GEMM Time = %20.6f s",moinfo->get_dgemm_timing());
  fflush(outfile);

  memory_manager->MemCheck(outfile);

  delete model_space;
  delete moinfo;
  delete debugging;
  delete memory_manager;
  delete global_timer;

  _default_psio_lib_->close(PSIF_PSIMRCC_INTEGRALS,1);
  fflush(outfile);
  return Success;
}

}} /* End Namespaces */
