/*! \file
    \ingroup OPTKING
    \brief ENERGY_SAVE.CC Rollin King, 2002 
     function executes if optinfo.mode == MODE_GRAD_SAVE
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <ccfiles.h>
#include <libchkpt/chkpt.h>

namespace psi { //namespace optking {

void energy_save(void) {
  int i,j,dim_carts,total_num_disps;
  double energy, *micro_e, *micro_grad, *grad;

  dim_carts = 3*optinfo.natom;
  fprintf(outfile,"Saving energy.\n");

  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num), sizeof(int));
  psio_read_entry(PSIF_OPTKING, "OPT: Num. of disp.",
      (char *) &(total_num_disps), sizeof(int));

  micro_e = new double[total_num_disps];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(micro_e[0]), total_num_disps*sizeof(double));
  close_PSIF();

  chkpt_init(PSIO_OPEN_OLD);
  energy = chkpt_rd_etot();
  chkpt_close();

  micro_e[optinfo.disp_num] = energy;

  open_PSIF();
  psio_write_entry(PSIF_OPTKING,"OPT: Displaced energies",
      (char *) &(micro_e[0]), total_num_disps*sizeof(double));
  close_PSIF();
  fprintf(outfile,"Energy written: %15.10lf\n", energy);

  delete [] micro_e;

  // increment disp_num
  open_PSIF();
  optinfo.disp_num += 1;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num), sizeof(int));

  /* delete CC temporary files */
  if((strcmp(optinfo.wfn, "MP2")==0)       || (strcmp(optinfo.wfn, "CCSD")==0)  ||
     (strcmp(optinfo.wfn, "CCSD_T")==0)    || (strcmp(optinfo.wfn, "EOM_CCSD")==0)  ||
     (strcmp(optinfo.wfn, "LEOM_CCSD")==0) || (strcmp(optinfo.wfn, "BCCD")==0)  ||
     (strcmp(optinfo.wfn,"BCCD_T")==0)     || (strcmp(optinfo.wfn, "SCF")==0)  ||
     (strcmp(optinfo.wfn,"CIS")==0)        || (strcmp(optinfo.wfn,"RPA")==0)  ||
     (strcmp(optinfo.wfn,"CC2")==0)        || (strcmp(optinfo.wfn,"CC3")==0)  ||
     (strcmp(optinfo.wfn,"EOM_CC3")==0) ) {
      fprintf(outfile, "Deleting CC binary files\n");
      for(i=CC_MIN; i <= CC_MAX; i++) {
        psio_open(i,1); psio_close(i,0);
      }
      psio_open(35,1); psio_close(35,0);
      psio_open(72,1); psio_close(72,0);
  }
  else if ( (strcmp(optinfo.wfn, "DETCI")==0) ) {
    fprintf(outfile, "Deleting CI binary files\n");
    psio_open(35,1); psio_close(35,0);
    psio_open(72,1); psio_close(72,0);
    psio_open(50,1); psio_close(50,0);
    psio_open(51,1); psio_close(51,0);
    psio_open(52,1); psio_close(52,0);
    psio_open(53,1); psio_close(53,0);
  }

  close_PSIF();

  if (optinfo.disp_num == total_num_disps) {
    fprintf(outfile,"Last displacement done, resetting checkpoint prefix.\n");
    chkpt_init(PSIO_OPEN_OLD);
    chkpt_reset_prefix();
    chkpt_commit_prefix();
    chkpt_close();
  }

  return ;
}

}//} /* namespace psi::optking */

