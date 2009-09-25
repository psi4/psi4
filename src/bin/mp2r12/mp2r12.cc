/*! \defgroup MP2R12 mp2r12: MP2-R12 energy evaluation */

/*! 
** \file
** \ingroup MP2R12
** \brief MP2-R12 energy evaluation
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"
#include "defines.h"

/* First definitions of globals */
extern "C"{
  FILE *infile, *outfile;
  char *psi_file_prefix;
}

namespace psi{
  namespace mp2r12{
  int *ioff;
  struct MOInfo moinfo;
  struct Params params;

  /* Function Prototypes */
  void init_io(int argc, char *argv[]);
  void title();
  void init_ioff();
  void get_parameters();
  void get_moinfo();
  void energy();
  void exit_io();
  void mp2r12_energy();
  } /* Namespace mp2r12 */
} /* Namespace psi */

int main(int argc, char *argv[])
{
  using namespace psi::mp2r12;
  init_io(argc,argv);
  title();
  get_parameters();
  get_moinfo();
  init_ioff();
  mp2r12_energy();
  exit_io();
  exit(0);
}

namespace psi{
  namespace mp2r12{
  void init_io(int argc, char *argv[])
  {
   psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
   tstart(outfile);
   ip_cwk_add(":MP2R12");
   
   psio_init(); psio_ipv1_config();
   chkpt_init(PSIO_OPEN_OLD);
  }

  void title(void)
  {
   fprintf(outfile, "\n");
   fprintf(outfile,"\t---------------------------------------------------\n");
   fprintf(outfile,"\t     MP2R12:  Program to determine the MP2-R12     \n");
   fprintf(outfile,"\t          energy of closed-shell systems.          \n");
   fprintf(outfile,"\t                                                   \n");
   fprintf(outfile,"\t                  Edward Valeev                    \n");
   fprintf(outfile,"\t                    July 1999                      \n");
   fprintf(outfile,"\t                                                   \n");
   fprintf(outfile,"\t        based on Daniel Crawford's MP2 code        \n");
   fprintf(outfile,"\t---------------------------------------------------\n");
   fprintf(outfile, "\n\n");
   fflush(outfile);
  }

  void init_ioff(void)
  {
   int i;
   ioff = (int *) malloc(IOFF * sizeof(int));
   if(ioff == NULL) {
     fprintf(stderr, "(mp2r12): error malloc'ing ioff array\n");
     abort();
   }
   ioff[0] = 0;
   for(i=1; i < IOFF; i++) {
     ioff[i] = ioff[i-1] + i;
   }
  }

  void exit_io(void)
  {
   chkpt_close();
   psio_done();
   tstop(outfile);
   psi_stop(infile,outfile,psi_file_prefix);
  }

  void get_parameters(void)
  {
   int errcod, tol;

   params.print_lvl = 0;
   errcod = ip_data("PRINT_LVL", "%d", &(params.print_lvl),0);

   params.wfn = strdup("MP2-R12/A");

   params.tolerance = 1e-10;
   errcod = ip_data("TOLERANCE","%d",&(tol),0);
   if(errcod == IPE_OK) {
       params.tolerance = 1.0*pow(10.0,(double) -tol);
     }
   params.keep_integrals=1;
   errcod = ip_boolean("KEEP_INTEGRALS", &(params.keep_integrals),0);

   params.c_limit = 0;
   errcod = ip_boolean("C_LIMIT", &(params.c_limit), 0);

   fprintf(outfile,"\tInput Parameters:\n");
   fprintf(outfile,"\t-----------------\n");
   fprintf(outfile,"\tWavefunction           =  %s\n", params.wfn);
   fprintf(outfile,"\tPrint Level            =  %d\n", params.print_lvl);
   fprintf(outfile,"\tTolerance              =  %3.1e\n", params.tolerance);
   fprintf(outfile,"\tKeep Integrals         =  %s\n", (params.keep_integrals ? "Yes": "No"));
   if (params.c_limit)
     fprintf(outfile,"\t*WARNING* : using the asymptotic limit for the C-coefficients!\n\n");
   fflush(outfile);

   return;
  }

  void get_moinfo(void)
  {
   int i,errcod,size;

   moinfo.frdocc = get_frzcpi();
   moinfo.fruocc = get_frzvpi();

   moinfo.nmo = chkpt_rd_nmo();
   moinfo.nirreps = chkpt_rd_nirreps();
   moinfo.iopen = chkpt_rd_iopen();
   moinfo.labels = chkpt_rd_irr_labs();
   moinfo.orbspi = chkpt_rd_orbspi();
   moinfo.clsdpi = chkpt_rd_clsdpi();
   moinfo.openpi = chkpt_rd_openpi();
   moinfo.enuc = chkpt_rd_enuc();
   moinfo.escf = chkpt_rd_escf();
   moinfo.evals = chkpt_rd_evals();
   
   moinfo.virtpi = init_int_array(moinfo.nirreps);
   for(i=0; i < moinfo.nirreps; i++) {
       moinfo.virtpi[i] = moinfo.orbspi[i]-moinfo.clsdpi[i]-moinfo.openpi[i];
     }

   fprintf(outfile,"\n\tFile30 Parameters:\n");
   fprintf(outfile,"\t------------------\n");
   fprintf(outfile,"\tNumber of irreps = %d\n",moinfo.nirreps);
   fprintf(outfile,"\tNumber of MOs    = %d\n\n",moinfo.nmo);
   fprintf(outfile,"\tLabel\t# MOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
   fprintf(outfile,"\t-----\t-----\t------\t------\t------\t------\t------\n");
   for(i=0; i < moinfo.nirreps; i++) {
       fprintf(outfile,
               "\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
               moinfo.labels[i],moinfo.orbspi[i],moinfo.frdocc[i],
               moinfo.clsdpi[i],moinfo.openpi[i],moinfo.virtpi[i],
               moinfo.fruocc[i]);
     }
   fprintf(outfile,"\n\tNuclear Repulsion Energy    = %20.10f\n",moinfo.enuc);
   fprintf(outfile,  "\tTotal SCF Energy            = %20.10f\n",moinfo.escf);
   fflush(outfile);

   /* Set up some other useful parameters */
   moinfo.noeints = moinfo.nmo*(moinfo.nmo+1)/2;
   moinfo.nteints = moinfo.noeints*(moinfo.noeints+1)/2;

   /* Construct first and last indexing arrays for each symblk
      in Pitzer ordering */
   moinfo.first = init_int_array(moinfo.nirreps);
   moinfo.last = init_int_array(moinfo.nirreps);
   moinfo.first[0] = 0;
   moinfo.last[0] = moinfo.orbspi[0] - 1;
   for(i = 1; i < moinfo.nirreps; i++) {
     moinfo.first[i] = moinfo.first[i-1] + moinfo.orbspi[i-1];
     moinfo.last[i] = moinfo.last[i-1] + moinfo.orbspi[i];
   }

   return;
  }

  } /* Namespace mp2r12 */
} /* Namespace psi */
