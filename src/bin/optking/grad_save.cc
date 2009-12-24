/*! \file
    \ingroup OPTKING
    \brief ENERGY_SAVE.CC Rollin King, 2002
    function executes if optinfo.mode == MODE_ENERGY_SAVE
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <ccfiles.h>

namespace psi { namespace optking {

void grad_save(const cartesians &carts) {
  int i,j,dim_carts,total_num_disps;
  double energy, *micro_e, *micro_grad, *grad, **ggrad, **refgrad, **rref;

  dim_carts = 3*optinfo.natom;

  fprintf(outfile,"Saving gradient and energy.\n");

  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num), sizeof(int));
  psio_read_entry(PSIF_OPTKING, "OPT: Num. of disp.",
      (char *) &(total_num_disps), sizeof(int));

  micro_e = new double[total_num_disps];
  micro_grad = new double [total_num_disps*3*carts.get_natom()*sizeof(double)];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(micro_grad[0]), total_num_disps*3*carts.get_natom()*sizeof(double));

  /*
  fprintf(outfile,"gradients\n");
  for (i=0; i<total_num_disps*3*carts.get_natom(); ++i)
    fprintf(outfile,"%15.10lf\n",micro_grad[i]);
    */

  psio_read_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(micro_e[0]), total_num_disps*sizeof(double));
  close_PSIF();

  chkpt_init(PSIO_OPEN_OLD);
  grad = chkpt_rd_grad();
  rref = chkpt_rd_rref();
  energy = chkpt_rd_etot();
  chkpt_close();

  // Rotate the gradient back to the reference frame in which all geometries were generated and stores
  int natoms = carts.get_natom();
  ggrad = init_matrix(natoms,3);
  int atomxyz=0;
  for(int atom=0; atom<natoms; atom++)
    for(int xyz=0; xyz<3; xyz++,atomxyz++)
      ggrad[atom][xyz] = grad[atomxyz];
  delete[] grad;
  refgrad = init_matrix(carts.get_natom(),3);
  opt_mmult(ggrad,0,rref,0,refgrad,0,carts.get_natom(),3,3,0);
  free_matrix(rref);
  free_matrix(ggrad);

  micro_e[optinfo.disp_num] = energy;
  for (i=0; i<dim_carts; ++i) {
    int atom = i/3;
    int xyz = i%3;
    micro_grad[3*carts.get_natom()*(optinfo.disp_num)+i] = refgrad[atom][xyz];
  }

  open_PSIF();
  psio_write_entry(PSIF_OPTKING,"OPT: Displaced energies",
      (char *) &(micro_e[0]), total_num_disps*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(micro_grad[0]), total_num_disps*3*carts.get_natom()*sizeof(double));
  close_PSIF();

  delete [] micro_e; delete [] micro_grad;

  // increment disp_num
  open_PSIF();
  optinfo.disp_num += 1;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num), sizeof(int));

    /* delete CC temporary files */
  if(!strcmp(optinfo.wfn, "MP2")       || !strcmp(optinfo.wfn, "CCSD")  ||
     !strcmp(optinfo.wfn, "CCSD_T")    || !strcmp(optinfo.wfn, "EOM_CCSD")  ||
     !strcmp(optinfo.wfn, "LEOM_CCSD") || !strcmp(optinfo.wfn, "BCCD")  ||
     !strcmp(optinfo.wfn,"BCCD_T")     || !strcmp(optinfo.wfn, "SCF")  ||
     !strcmp(optinfo.wfn,"CIS")        || !strcmp(optinfo.wfn,"RPA")  ||
     !strcmp(optinfo.wfn,"CC2")        || !strcmp(optinfo.wfn,"CC3")  ||
     !strcmp(optinfo.wfn,"EOM_CC3") ) {
     fprintf(outfile, "Deleting CC binary files\n");
     for(i=CC_MIN; i <= CC_MAX; i++) {
       psio_open(i,1); psio_close(i,0);
     }
    psio_open(35,1); psio_close(35,0);
    psio_open(72,1); psio_close(72,0);
  }
  else if (!strcmp(optinfo.wfn, "DETCI")) {
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

}} /* namespace psi::optking */

