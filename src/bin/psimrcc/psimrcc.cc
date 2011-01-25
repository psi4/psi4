/**
 *  @defgroup PSIMRCC PSIMRCC is a code for SR/MRCC computations
 *  @file psimrcc.cc
 *  @ingroup (PSIMRCC)
 *  @brief Contains psimrcc() and global variables
*/

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
#include <libqt/qt.h>

// PSI C++
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>

//#include "blas.h"
//#include "debugging.h"
//#include "git.h"
//#include "main.h"
//#include "sort.h"
//#include "mrcc.h"
//#include "psimrcc.h"
//#include "transform.h"


//#include <liboptions/liboptions.h>
//#include <libmoinfo/libmoinfo.h>
//#include <libutil/libutil.h>
//
//#include "blas.h"
//#include "sort.h"
//#include "mrcc.h"
//#include "idmrpt2.h"
//#include "mp2_ccsd.h"
//#include "transform.h"
//#include "debugging.h"
//#include "psimrcc.h"
//#include "updater.h"


#include <psi4-dec.h>


using namespace std;
using namespace psi;

namespace psi{ namespace psimrcc{

//namespace psi{
//MOInfo              *moinfo;
//ModelSpace          *model_space;
//MemoryManager       *_memory_manager_;
//
//namespace psimrcc{
//
//// Global variables
//Timer               *global_timer;
//Debugging           *debugging;
//CCBLAS              *blas;
//CCSort              *sorter;
//CCTransform         *trans = NULL;
//
//void read_calculation_options();
//}} /* End Namespaces */

/**
 * The main function
 * @param argc
 * @param argv[]
 * @return PSI_RETURN_SUCCESS if the program ran without any problem
 */
PsiReturnType psimrcc(Options &options)
{
//  using namespace psi::psimrcc;

//  run_psimrcc();

  fprintf(outfile,"\n\n  PSIMRCC job completed.");
  fflush(outfile);
//
//  _memory_manager_->MemCheck(outfile);
//
  return Success;
}

//
//
///**
// * Start Psi3, draw the program logo, version, and compilation details
// * @param argc
// * @param argv[]
// */
//void init_psi(int argc, char *argv[])
//{
//  int num_extra_args=0;
//  char**  extra_args;
//
//  extra_args = new char*[argc];
//
//  for(int i=1; i<argc; i++) {
//    extra_args[num_extra_args++] = argv[i];
///*  Template for argument parsing
//    if(strcmp(argv[i], "--opdm") == 0) {
//    }
//    else {
//      extra_args[num_extra_args++] = argv[i];
//    }
//*/
//  }
//  psi_start(&infile,&outfile,&psi_file_prefix,num_extra_args,extra_args,0);
//  delete[] extra_args;
//
//  _default_psio_lib_ = new PSIO();
//  psiopp_ipv1_config(_default_psio_lib_);
//  _default_chkpt_lib_ = new Chkpt(_default_psio_lib_, PSIO_OPEN_OLD);
//
//  ip_cwk_clear();
//  ip_cwk_add(const_cast<char*>(":PSIMRCC"));
//  ip_cwk_add(const_cast<char*>(":MRCC"));
//  ip_cwk_add(const_cast<char*>(":PSI"));
//  ip_cwk_add(const_cast<char*>(":INPUT"));
//
//  _memory_manager_    = new MemoryManager();
//
//  options_init();
//
//  tstart(outfile);
//
//  _default_psio_lib_->open(PSIF_PSIMRCC_INTEGRALS,PSIO_OPEN_NEW);
//
//  fprintf(outfile,"\n  MRCC          MRCC");
//  fprintf(outfile,"\n   MRCC  MRCC  MRCC");
//  fprintf(outfile,"\n   MRCC  MRCC  MRCC      PSIMRCC Version 0.9.3.3, July 2009");
//  fprintf(outfile,"\n   MRCC  MRCC  MRCC      Multireference Coupled Cluster, written by");
//  fprintf(outfile,"\n     MRCCMRCCMRCC        Francesco A. Evangelista and Andrew C. Simmonett");
//  fprintf(outfile,"\n         MRCC            Compiled on %s at %s",__DATE__,__TIME__);
//  fprintf(outfile,"\n         MRCC            id =%s",GIT_ID);
//  fprintf(outfile,"\n       MRCCMRCC");
//
//  global_timer = new Timer;
//  read_calculation_options();
//
//  ip_cwk_clear();
//  ip_cwk_add(const_cast<char*>(":MRCC"));
//  ip_cwk_add(const_cast<char*>(":PSIMRCC"));
//
//  debugging = new Debugging;
//
//  moinfo = new MOInfo();
//
//  model_space = new ModelSpace(moinfo);
//
//  moinfo->setup_model_space();  // The is a bug here DELETEME
//}
//
///**
// * Close psi by calling psio_done() and psi_stop()
// */
//void close_psi()
//{
//  /***********************
//    Close the checkpoint
//  ***********************/
//
//  delete model_space;
//  delete moinfo;
//  delete debugging;
//  delete _memory_manager_;
//  delete global_timer;
//
//  fflush(outfile);
//  options_close();
//
//  _default_psio_lib_->close(PSIF_PSIMRCC_INTEGRALS,1);
//
//  delete _default_chkpt_lib_;
//  delete _default_psio_lib_;
//
//  tstop(outfile);
//  psi_stop(infile,outfile,psi_file_prefix);
//}
//
//
//
//void run_psimrcc()
//{
//  blas   = new CCBLAS();
//  trans  = new CCTransform();
//
//  if(options_get_str("CORR_WFN")=="PT2"){
//    mrpt2();
//  }else if(options_get_str("CORR_WFN")=="MP2-CCSD"){
//    mp2_ccsd();
//  }else{
//    mrccsd();
//  }
//
//  delete sorter;
//  delete trans;
//  delete blas;
//}
//
///*!
// * Runs a MRPT2 and a MRCCSD computation
// * @todo move this code in the CCMRCC class
// */
//void mrccsd()
//{
//  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
//  CCMRCC        mrcc;
//
//  if(options_get_bool("PERT_CBS") && options_get_bool("PERT_CBS_COUPLING")){
//    mrcc.compute_first_order_amps();
//  }
//
//  // Initialize the appropriate updater
//  Updater* updater;
////  if(options_get_str("CORR_ANSATZ")=="SR")
////    updater = static_cast<Updater*>(new MkUpdater());
//  if(options_get_str("CORR_ANSATZ")=="MK")
//    updater = dynamic_cast<Updater*>(new MkUpdater());
//  if(options_get_str("CORR_ANSATZ")=="BW")
//    updater = dynamic_cast<Updater*>(new BWUpdater());
//
//	// Compute the energy
//  mrcc.compute_energy(updater);
//
//  if(options_get_bool("PERT_CBS")){
//    mrcc.perturbative_cbs();
//  }
//
//  delete updater;
//}
//
///*!
// * Runs a CCSD_MP2 computation
// */
//void mp2_ccsd()
//{
//  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
//  MP2_CCSD        mp2_ccsd;
//
//  // Compute the initial amplitudes and CCSD_MP2 energy
//  mp2_ccsd.compute_mp2_ccsd_energy();
//
//  DEBUGGING(1,
//    blas->print_memory();
//  )
//}
//
///*!
// * Runs a MRPT2 computation
// */
//void mrpt2()
//{
//  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
//  IDMRPT2        idmrpt2;
//
//  Updater* updater = dynamic_cast<Updater*>(new MkUpdater());
//
//  // Compute the initial amplitudes and MP2 energy
//  idmrpt2.compute_mrpt2_energy(updater);
//
//  delete updater;
//
//  DEBUGGING(1,
//    blas->print_memory();
//  )
//}
//
///*!
// * Runs a integral transformation
// * @todo CCTransform is still unused in the code
// */
//void transform_integrals()
//{
////   CCTransform transf;
////   transf.read_so_integrals();
//}

}} /* End Namespaces */


