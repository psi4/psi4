/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
/*
**  ccresponse: Program to compute CC linear response properties.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <physconst.h>
#include <psifiles.h>
#include "Params.h"
#include "MOInfo.h"
#include "Local.h"
#include "globals.h"

namespace psi { namespace ccresponse {

/* Max length of ioff array */
#define IOFF_MAX 32641

/* Function prototypes */
void init_io(void);
void init_ioff(void);
void title(void);
void get_moinfo(void);
void get_params(Options &);
void cleanup(void);
void exit_io(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);
void hbar_extra(void);
void cc2_hbar_extra(void);
void sort_lamps(void);
void lambda_residuals(void);

void local_init(void);
void local_done(void);

void polar(void);
void optrot(void);
void roa(void);

void preppert(void);

int ccresponse(Options &options)
{
  int **cachelist, *cachefiles;

  init_io();
  init_ioff();
  title();
  get_moinfo();
  get_params(options);

  timer_on("ccresponse");

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /*** UHF references ***/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist,
         NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi, moinfo.avir_sym,
         moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  }
  else { /*** RHF/ROHF references ***/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
         2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);
  }

  if(params.local) local_init();

  if (params.wfn == "CC2") {
    cc2_hbar_extra();
  }
  else {
    hbar_extra();
  }

  sort_lamps(); /* should be removed sometime - provided by cclambda */
  if(params.wfn != "CC2") lambda_residuals(); /* don't do this for CC2 */

  preppert();

  if(params.prop == "POLARIZABILITY") polar();
  if(params.prop == "ROTATION") optrot();
  if(params.prop == "ROA") roa();

  if(params.local) local_done();

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup();

  timer_off("ccresponse");

  exit_io();
  return 0;
}

void init_io(void)
{
  int i;

  tstart();

  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i, 1);

  /* Clear out DIIS TOC Entries */
  psio_close(CC_DIIS_AMP, 0);
  psio_close(CC_DIIS_ERR, 0);

  psio_open(CC_DIIS_AMP, 0);
  psio_open(CC_DIIS_ERR, 0);
}

void title(void)
{
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*       ccresponse       *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;

  /* Close all dpd data files here */
  for(i=CC_MIN; i < CC_TMP; i++) psio_close(i,1);
  for(i=CC_TMP; i <= CC_TMP11; i++) psio_close(i,0);  /* get rid of TMP files */
  for(i=CC_TMP11+1; i <= CC_MAX; i++) psio_close(i,1);

  tstop();
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}

}} // namespace psi::ccresponse
