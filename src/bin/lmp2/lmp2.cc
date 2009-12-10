/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
** \file
** \ingroup LMP2
** \LMP2 evaluation of energy
*/

#include "mpi.h"
//#include <iostream>
//#include <fstream>              // file I/O support
//#include <cstdlib>              // support for exit()
//#include <cstring>
//#include <memory.h>
//#include <assert.h>
//#include <pthread.h>
//#include <libipv1/ip_lib.h>
//#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
//#include <libmints/matrix.h>
#include "globals.h"

//using namespace std;
//using namespace psi;

//std::string to_string(const int val); // In libmints/matrix.cc

namespace psi{

extern int myid;
extern int nprocs;

namespace lmp2{

int myid_lmp2 = myid;
int nprocs_lmp2 = nprocs;

PsiReturnType lmp2(Options &options, int argc, char * argv[]) {

  shared_ptr<PSIO> psio_obj(new PSIO);
  psiopp_ipv1_config(psio_obj);
  shared_ptr<Chkpt> chkpt_obj(new Chkpt(psio_obj, PSIO_OPEN_OLD));

  if(myid == 0) {
    tstart();

    fprintf(outfile, "\t\t\t*************************\n");
    fprintf(outfile, "\t\t\t*                       *\n");
    fprintf(outfile, "\t\t\t*         LMP2          *\n");
    fprintf(outfile, "\t\t\t*                       *\n");
    fprintf(outfile, "\t\t\t*************************\n");
    fprintf(outfile, "\t\t\tRunning on %d processors\n", nprocs);
    fflush(outfile);


  }

  /** Compute the LMP2 energy **/
  {
  LMP2 lmp2_obj(psio_obj, chkpt_obj);
  timer_init();

//  MPIMessageGrp(&argc, &argv);


  lmp2_obj.get_params(options);
  lmp2_obj.get_moinfo();
  lmp2_obj.opdm();
  lmp2_obj.localize();
  lmp2_obj.get_fock();
  lmp2_obj.domains();
  lmp2_obj.projection();

  lmp2_obj.direct_transformation();
  lmp2_obj.allocate_T();

  lmp2_obj.iterate();
  }
  /** LMP2 complete **/

  timer_done();
  if(myid == 0) {
    tstop();
  }
  
//  MPI_Finalize();
  return Success;
}

}}

