/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
** \file
** \ingroup LMP2
** \LMP2 evaluation of energy
*/

//#include "mpi.h"
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
    fprintf(outfile, "\n**************************** Begin LMP2 *********************************\n\n\n");


    fprintf(outfile, "\t\t\t*************************\n");
    fprintf(outfile, "\t\t\t*                       *\n");
  if(options.get_bool("RI_LMP2")){
    fprintf(outfile, "\t\t\t*       DF-LMP2         *\n");
  }
  else {
    fprintf(outfile, "\t\t\t*         LMP2          *\n");
  }
    fprintf(outfile, "\t\t\t*                       *\n");
    fprintf(outfile, "\t\t\t*************************\n");
    fprintf(outfile, "\t\t\tRunning on %d processors\n", nprocs);
    fflush(outfile);


  }

  /** Compute the LMP2 energy **/
  {
  LMP2 lmp2_obj(psio_obj, chkpt_obj);
  timer_init();

  if(myid == 0)
    timer_on("GETPARAMS");
  lmp2_obj.get_params(options);
  if(myid == 0)
    timer_off("GETPARAMS");

  if(myid == 0)
    timer_on("GETMOINFO");
  lmp2_obj.get_moinfo();
  if(myid == 0)
    timer_off("GETMOINFO");

  if(myid == 0)
    timer_on("OPDM");
  lmp2_obj.opdm();
  if(myid == 0)
    timer_off("OPDM");

  if(myid == 0)
    timer_on("LOCALIZE");
  lmp2_obj.localize();
  if(myid == 0)
    timer_off("LOCALIZE");

  if(myid == 0)
    timer_on("GETFOCK");
  lmp2_obj.get_fock();
  if(myid == 0)
    timer_off("GETFOCK");

  if(myid == 0)
    timer_on("DOMAIN");
  lmp2_obj.domains();
  if(myid == 0)
    timer_off("DOMAIN");

  if(myid == 0)
    timer_on("PROJECT");
  lmp2_obj.projection();
  if(myid == 0)
    timer_off("PROJECT");

  if(options.get_bool("RI_LMP2")){
    lmp2_obj.direct_df_transformation();
  }
  else {
    if(myid == 0)
        timer_on("TRANS");
    lmp2_obj.direct_transformation();
    if(myid == 0)
        timer_off("TRANS");
  }


  lmp2_obj.allocate_T();

  if(myid == 0)
    timer_on("ITERATE");
  lmp2_obj.iterate();
  if(myid == 0)
    timer_off("ITERATE");


  }
  /** LMP2 complete **/

  timer_done();
  if(myid == 0) {
    fprintf(outfile, "\n**************************** End of LMP2 *********************************\n");

    tstop();
  }
  
//  MPI_Finalize();
  return Success;
}

}}

