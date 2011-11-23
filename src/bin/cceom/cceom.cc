/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
/*
**  CCEOM: Program to calculate the EOM CCSD right-hand eigenvector and  energy
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <psi4-dec.h>
#include "Params.h"
#include "MOInfo.h"
#include "Local.h"
#include "globals.h"

namespace psi { namespace cceom {

void init_io(int argc, char *argv[]);
void get_moinfo(void);
void cleanup(void);
void exit_io(void);
void diag(void);
void get_params(Options &);
void get_eom_params(Options &);
void form_dpd_dp(void);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void sort_amps(void);
void hbar_norms(void);

/* local correlation functions */
void local_init(void);
void local_done(void);

}} // namespace psi::cceom

using namespace psi;

PsiReturnType cceom(Options &options)
{
  int i, h, done=0, *cachefiles, **cachelist;
  init_io(argc, argv);
  fprintf(outfile,"\n\t**********************************************************\n");
  fprintf(outfile,"\t*  CCEOM: An Equation of Motion Coupled Cluster Program  *\n");
  fprintf(outfile,"\t**********************************************************\n");

  get_moinfo();
  fflush(outfile);
  get_params();
  get_eom_params();
#ifdef TIME_CCEOM
  timer_on("CCEOM");
#endif

  form_dpd_dp();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if (params.ref == 2) { /* UHF */
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);
    /* cachelist = init_int_matrix(32,32); */

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles,
    cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
    moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  }
  else { /* RHF or ROHF */
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);
    /* cachelist = init_int_matrix(12,12); */

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles,
           cachelist, NULL, 2, moinfo.occpi, moinfo.occ_sym,
           moinfo.virtpi, moinfo.vir_sym);
  }

  if(params.local) local_init();

  diag();

  dpd_close(0);
  if(params.local) local_done();
  cleanup();
#ifdef TIME_CCEOM
  timer_off("CCEOM");
#endif
  exit_io();
  exit(0);
}

namespace psi { namespace cceom {

void init_io(int argc, char *argv[])
{
  int i, num_unparsed;
  char *progid, *argv_unparsed[100];

  eom_params.dot_with_L = 0;
  eom_params.guess = NULL;
  for (i=1, num_unparsed=0; i<argc; ++i) {
    if (!strcmp(argv[i],"--dot_with_L"))
      eom_params.dot_with_L = 1;
    else if(!strcmp(argv[i], "--reuse"))
      eom_params.guess = strdup("DISK");
    else
      argv_unparsed[num_unparsed++] = argv[i];
  }

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(&infile,&outfile,&psi_file_prefix,num_unparsed, argv_unparsed, 0);
  ip_cwk_add(":INPUT");
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);
  psio_init(); psio_ipv1_config();

  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
  /*
  for(i=CC_MIN; i <= CC_MISC; i++) psio_open(i,1);
  for(i=CC_TMP; i <= EOM_D; i++) psio_open(i,0);
  for(i=EOM_Cme; i <= CC_MAX; i++) psio_open(i,0);
  */

  /* it will read it anyway if its there ?
  if(eom_params.guess != NULL) {
    if(!strcmp(eom_params.guess,"DISK"))
      psio_open(EOM_CME,1);
  }
  else
    psio_open(EOM_CME,0);
    */

  if (eom_params.dot_with_L) {
    psio_read_entry(CC_INFO,"EOM L0",(char *) &eom_params.L0, sizeof(double));
    psio_read_entry(CC_INFO,"EOM L Irrep",(char *) &eom_params.L_irr, sizeof(int));
  }
}

void exit_io(void)
{
  int i;
  /* for(i=CC_MIN; i <= CC_MISC; i++) psio_close(i,1);
  for(i=CC_TMP; i <= CC_MAX; i++) psio_close(i,1);
  */
  for(i=CC_MIN; i <= CC_DIIS_AMP; i++) psio_close(i,1);
  for(i=CC_TMP; i <= CC_TMP11; i++) psio_close(i,0);
  for(i=CC_TMP11+1; i <= CC_MAX; i++) psio_close(i,1);

  psio_done();
  tstop(outfile);
  psi_stop(infile,outfile,psi_file_prefix);
}

void form_dpd_dp(void) {
  int h, h0, h1, cnt, nirreps;
  nirreps = moinfo.nirreps;

  dpd_dp = (int ***) malloc(nirreps * sizeof(int **));
  for(h=0; h < nirreps; h++) {
      dpd_dp[h] = init_int_matrix(nirreps,2);
      cnt=0;
      for(h0=0; h0 < nirreps; h0++) {
          for(h1=0; h1 < nirreps; h1++) {
              if((h0^h1)==h) {
                  dpd_dp[h][cnt][0] = h0;
                  dpd_dp[h][cnt++][1] = h1;
                }
            }
        }
    }
}

}} // namespace psi::cceom
