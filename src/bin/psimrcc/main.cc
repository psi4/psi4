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
    MemoryManager       *_memory_manager_;

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

  model_space = new ModelSpace(moinfo);

  moinfo->setup_model_space();  // The is a bug here DELETEME

  blas   = new CCBLAS(options);
  trans  = new CCTransform();
  if(options.get_str("CORR_WFN")=="PT2"){
    mrpt2(options);
  }else if(options.get_str("CORR_WFN")=="MP2-CCSD"){
    mp2_ccsd(options);
  }else{
    mrccsd(options);
  }

  delete sorter;
  delete trans;
  delete blas;


  fprintf(outfile,"\n\n  PSIMRCC job completed.");
  fprintf(outfile,"\n  Wall Time = %20.6f s",global_timer->get());
  fprintf(outfile,"\n  GEMM Time = %20.6f s",moinfo->get_dgemm_timing());
  fflush(outfile);

  _memory_manager_->MemCheck(outfile);

  delete model_space;
  delete moinfo;
  delete debugging;
  delete _memory_manager_;
  delete global_timer;

  _default_psio_lib_->close(PSIF_PSIMRCC_INTEGRALS,1);
  fflush(outfile);
  return Success;
}

}} /* End Namespaces */
