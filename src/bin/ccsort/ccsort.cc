/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here
*/
/*
**  CCSORT: Program to reorganize integrals for CC and MBPT calculations.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <liboptions/liboptions.h>
#include <psifiles.h>
#include <psi4-dec.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#include "globals.h"

namespace psi { namespace ccsort {

#define IOFF_MAX 32641

void init_io();
void init_ioff(void);
void title(void);
void get_params(Options & options);
void get_moinfo(void);
void sort_oei(void);
void sort_tei(void);
void b_sort(void);
void c_sort(void);
void d_sort(void);
void e_sort(void);
void f_sort(void);
void d_spinad(void);
void e_spinad(void);
void f_spinad(void);
void scf_check(void);
void fock(void);
void denom(void);
void exit_io(void);
void cleanup(void);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);
void local_init(Options & options);
void local_done(void);
void cc_memcheck(void);

int ccsort(Options &options)
{
  int i;
  int **cachelist, *cachefiles;
  int bamount,famount; /* Truncated theoretical number of B/F-type ints that
                          could be stored in cache at once                  */

  unsigned long int ia_size, ab_size, ij_size, f_size, t2_size, b_size;

  init_io();
  init_ioff();
  title();

  get_params(options);
  get_moinfo();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /*** UHF references ***/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist,
         NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi, moinfo.avir_sym,
         moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  }
  else { /*** RHF/ROHF references ***/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist,
         NULL, 2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi,
         moinfo.vir_sym);
  }

  /* run a small computation of memory and disk requirements */
  cc_memcheck();

  fprintf(outfile, "\n");

  sort_oei();
  sort_tei();
  c_sort();
  d_sort();
  e_sort();
  f_sort();
  if(params.ref == 0) {
    d_spinad();
    e_spinad();
/*     f_spinad(); */
  }
  scf_check();
  fock();
  denom();

  /* CPHF stuff for local correlation tests */
  if(params.local) {
    local_init(options);
    local_done();
  }

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup();
  exit_io();
  return (Success);
}

void init_io()
{
//  int i, num_unparsed;
//  char **argv_unparsed;

//  argv_unparsed = (char **) malloc(argc * sizeof(char *));
  params.reset = 0;
//  for(i=1, num_unparsed=0; i < argc; i++) {
//    if(!strcmp(argv[i], "--reset")) params.reset = 1;
//    else argv_unparsed[num_unparsed++] = argv[i];
//  }

//  free(argv_unparsed);

  tstart();

  int i;
  if(params.reset) for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,0);
  else for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*         CCSORT         *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

void exit_io(void)
{
  int i;
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

void cleanup(void)
{
  int i;

  psio_write_entry(CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
           sizeof(double));

  if(params.ref == 2) {

    free(moinfo.pitz2qt_A);
    free(moinfo.pitz2qt_B);
    free(moinfo.qt2pitz_A);
    free(moinfo.qt2pitz_B);

    free(moinfo.aocc);
    free(moinfo.bocc);
    free(moinfo.avir);
    free(moinfo.bvir);
    free(moinfo.all_aocc);
    free(moinfo.all_bocc);
    free(moinfo.all_avir);
    free(moinfo.all_bvir);
    free(moinfo.aoccpi);
    free(moinfo.boccpi);
    free(moinfo.avirtpi);
    free(moinfo.bvirtpi);
    free(moinfo.all_aoccpi);
    free(moinfo.all_boccpi);
    free(moinfo.all_avirtpi);
    free(moinfo.all_bvirtpi);

    free(moinfo.cc_aocc);
    free(moinfo.cc_bocc);
    free(moinfo.cc_avir);
    free(moinfo.cc_bvir);
    free(moinfo.qt_aocc);
    free(moinfo.qt_bocc);
    free(moinfo.qt_avir);
    free(moinfo.qt_bvir);
    free(moinfo.aocc_sym);
    free(moinfo.bocc_sym);
    free(moinfo.avir_sym);
    free(moinfo.bvir_sym);

    free(moinfo.cc_allaocc);
    free(moinfo.cc_allbocc);
    free(moinfo.cc_allavir);
    free(moinfo.cc_allbvir);
    free(moinfo.qt_allaocc);
    free(moinfo.qt_allbocc);
    free(moinfo.qt_allavir);
    free(moinfo.qt_allbvir);
    free(moinfo.allaocc_sym);
    free(moinfo.allbocc_sym);
    free(moinfo.allavir_sym);
    free(moinfo.allbvir_sym);

    free(moinfo.aocc_off);
    free(moinfo.bocc_off);
    free(moinfo.avir_off);
    free(moinfo.bvir_off);
    free(moinfo.all_aocc_off);
    free(moinfo.all_bocc_off);
    free(moinfo.all_avir_off);
    free(moinfo.all_bvir_off);
  }
  else {

    free(moinfo.pitz2qt);
    free(moinfo.qt2pitz);

    free(moinfo.occ);
    free(moinfo.vir);
    free(moinfo.all_occ);
    free(moinfo.all_vir);
    free(moinfo.socc);
    free(moinfo.all_socc);
    free(moinfo.occpi);
    free(moinfo.virtpi);
    free(moinfo.all_occpi);
    free(moinfo.all_virtpi);
    free(moinfo.orbsym);

    free(moinfo.cc_occ);
    free(moinfo.cc_vir);
    free(moinfo.qt_occ);
    free(moinfo.qt_vir);
    free(moinfo.occ_sym);
    free(moinfo.vir_sym);

    free(moinfo.cc_allocc);
    free(moinfo.cc_allvir);
    free(moinfo.qt_allocc);
    free(moinfo.qt_allvir);
    free(moinfo.allocc_sym);
    free(moinfo.allvir_sym);

    free(moinfo.occ_off);
    free(moinfo.vir_off);
    free(moinfo.all_occ_off);
    free(moinfo.all_vir_off);

  }

  free(moinfo.pitzer2qt);
  free(moinfo.qt2pitzer);

  free(moinfo.sopi);
  free(moinfo.orbspi);
  free(moinfo.clsdpi);
  free(moinfo.openpi);
  free(moinfo.uoccpi);
  // Wavefunction owns these arrays
//  free(moinfo.fruocc);
//  free(moinfo.frdocc);
  for(i=0; i < moinfo.nirreps; i++)
    free(moinfo.labels[i]);
  free(moinfo.labels);
  free(moinfo.frozen);
  free(ioff);
}

}} //namespace psi::ccsort

