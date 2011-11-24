/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
/*
**  CCEOM: Program to calculate the EOM CCSD right-hand eigenvector and  energy
*/
#include <cstdio>
#include <cstdlib>
#include <string>
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

void init_io(void);
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

namespace psi { namespace cceom {

PsiReturnType cceom(Options &options)
{
  int i, h, done=0, *cachefiles, **cachelist;
  init_io();
  fprintf(outfile,"\n\t**********************************************************\n");
  fprintf(outfile,"\t*  CCEOM: An Equation of Motion Coupled Cluster Program  *\n");
  fprintf(outfile,"\t**********************************************************\n");

  get_moinfo();
  fflush(outfile);
  get_params(options);
  get_eom_params(options);
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

void init_io(void)
{
  tstart();
  for(int i = CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
}

void exit_io(void)
{
  int i;
  for(i=CC_MIN; i <= CC_DIIS_AMP; i++) psio_close(i,1);
  for(i=CC_TMP; i <= CC_TMP11; i++) psio_close(i,0);
  for(i=CC_TMP11+1; i <= CC_MAX; i++) psio_close(i,1);
  tstop();
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
