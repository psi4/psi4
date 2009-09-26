/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
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

namespace psi { namespace cis {

void init_io(int argc, char *argv[]);
void title(void);
void exit_io(void);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void get_moinfo(void);
void get_params(Options & options);
void cleanup(void);
void build_A(void);
void d_corr(void);
void local_init(Options & options);
void local_done(void);
void amp_write_T1(dpdfile2 *T1, int length, FILE *outfile);
void diag(void);

PsiReturnType cis(Options & options, int argc, char *argv[])
{
  char lbl[32];
  int **cachelist, *cachefiles;
  int h, i, j, jj, k, nroot, a, max_j, max_a;
  double *weakp, value, d_value, **B_PAO;
  double *Bt, max_val, test_val, sum_val;
  int *spin, count, ivalue;
  dpdfile2 B;

  init_io(argc, argv);
  title();

  get_moinfo();
  get_params(options);

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    cachelist = cacheprep_rhf(0, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, 
        cachelist, NULL, 2, moinfo.occpi, moinfo.occ_sym, 
        moinfo.virtpi, moinfo.vir_sym);
  }
  else if(params.ref == 2) { /** UHF **/
    cachelist = cacheprep_uhf(0, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, 
        cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
        moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  }

  /* build components of CIS response matrix, A */
  build_A();

  /* diagonalize A within each irrep */
  diag();

  /* moved up with impunity??? */
  if(params.local) local_init(options);

  /* print eigenvalues and largest components of eigenvectors in ascending order */
  if(params.ref == 0) {

    fprintf(outfile, "\tRHF-CIS Singlet Excitation Energies:\n");
    fprintf(outfile, "\t------------------------------------\n\n");
    fprintf(outfile, "\t  Root Irrep       Hartree          eV          cm-1  \n");
    fprintf(outfile, "\t  ---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.singlet_evals[h][i];

        fprintf(outfile,   "CIS State %4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);

        fprintf(outfile, "\nLargest components of singlet excited wave function #%d/#%d:\n",
            h, i);
        sprintf(lbl, "BIA(%d)[%d] singlet", i, h);
        dpd_file2_init(&B, CC_OEI, h, 0, 1, lbl);
        amp_write_T1(&B, 5, outfile);
        dpd_file2_close(&B);

      }
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "\tRHF-CIS Triplet Excitation Energies:\n");
    fprintf(outfile, "\t------------------------------------\n\n");
    fprintf(outfile, "\t  Root Irrep       Hartree          eV          cm-1  \n");
    fprintf(outfile, "\t  ---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.triplet_evals[h][i];

        fprintf(outfile,   "CIS State %4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);

        fprintf(outfile, "\nLargest components of triplet excited wave function #%d/#%d:\n", h, i);
        sprintf(lbl, "BIA(%d)[%d] triplet", i, h);
        dpd_file2_init(&B, CC_OEI, h, 0, 1, lbl);
        amp_write_T1(&B, 5, outfile);
        dpd_file2_close(&B);
        fprintf(outfile, "\n");

      }
    }
    fprintf(outfile, "\n");

  }

  fflush(outfile);

  /* compute the (D) correction to each CIS singlet excitation energy */
  d_corr();

  fprintf(outfile, "\n");

  if(params.local) local_done();

  /* print corrected eigenvalues in ascending order */
  if(params.ref == 0) {

    fprintf(outfile, "\tRHF-CIS(D) Singlet Corrections:\n");
    fprintf(outfile, "\t-------------------------------\n");
    fprintf(outfile, "\t  Root Irrep       Hartree          eV          cm-1  \n");
    fprintf(outfile, "\t  ---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.singlet_d[h][i];

        fprintf(outfile,   "CIS State %4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);
      }
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "\tRHF-CIS(D) Singlet Excitation Energies:\n");
    fprintf(outfile, "\t---------------------------------------\n\n");
    fprintf(outfile, "\t  Root Irrep       Hartree          eV          cm-1  \n");
    fprintf(outfile, "\t  ---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.singlet_evals[h][i];
        value += moinfo.singlet_d[h][i];

        fprintf(outfile,   "CIS State %4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);
      }
    }
    fprintf(outfile, "\n");

    if(params.local) {
      fprintf(outfile, "\tRHF-CIS(D) Weak Pair Corrections:\n");
      fprintf(outfile, "\t---------------------------------\n");
      fprintf(outfile, "\t  Root Irrep       Hartree          eV          cm-1  \n");
      fprintf(outfile, "\t  ---- ----- ------------------  ---------  ----------\n");
      for(h=0; h < moinfo.nirreps; h++) {
        for(i=0; i < params.rpi[h]; i++) {
          value = moinfo.singlet_weakp[h][i];

          fprintf(outfile,   "CIS State %4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
              value, value*_hartree2ev, value*_hartree2wavenumbers);
        }
      }
      fprintf(outfile, "\n");
    }
  }

  else if(params.ref == 2) {

    fprintf(outfile, "\tUHF-CIS Excitation Energies:\n");
    fprintf(outfile, "\t----------------------------\n\n");
    fprintf(outfile, "\t  Root Irrep       Hartree          eV          cm-1\n");
    fprintf(outfile, "\t  ---- ----- ------------------  ---------  -----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.uhf_evals[h][i];

        fprintf(outfile,   "CIS State %4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);
      }
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "\tUHF-CIS(D) Corrections:\n");
    fprintf(outfile, "\t-----------------------\n\n");

    fprintf(outfile, "\t  Root Irrep       Hartree          eV         cm-1\n");
    fprintf(outfile, "\t  ---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.uhf_d[h][i];

        fprintf(outfile,   "CIS State %4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);
      }
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "\tUHF-CIS(D) Excitation Energies:\n");
    fprintf(outfile, "\t-------------------------------\n\n");

    fprintf(outfile, "\t  Root Irrep       Hartree          eV         cm-1\n");
    fprintf(outfile, "\t  ---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.uhf_evals[h][i];
        value += moinfo.uhf_d[h][i];

        fprintf(outfile,   "CIS State %4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);
      }
    }
    fprintf(outfile, "\n");

  }

  dpd_close(0);
  if(params.ref == 0 || params.ref == 1) cachedone_rhf(cachelist);
  else if(params.ref == 2) cachedone_uhf(cachelist);
  free(cachefiles);

  cleanup();

  exit_io();
  return(Success);
}

void init_io(int argc, char *argv[])
{
  tstart();

  psio_open(CC_INFO, 1);
  psio_open(CC_OEI, 1);
  psio_open(CC_CINTS, 1);
  psio_open(CC_DINTS, 1);
  psio_open(CC_EINTS, 1);
  psio_open(CC_FINTS, 1);
  psio_open(CC_DENOM, 1);
  psio_open(CC_MISC, 0);
  psio_open(CC_TMP0, 0);
  psio_open(CC_TMP1, 0);
}

void title(void)
{
  fprintf(outfile, "\t\t\t*************************\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*          CIS          *\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*************************\n");
}

void exit_io(void)
{
  psio_close(CC_INFO, 1);
  psio_close(CC_OEI, 1);
  psio_close(CC_CINTS, 1);
  psio_close(CC_DINTS, 1);
  psio_close(CC_EINTS, 1);
  psio_close(CC_FINTS, 1);
  psio_close(CC_DENOM, 1);
  psio_close(CC_MISC, 1);
  psio_close(CC_TMP0, 0);
  psio_close(CC_TMP1, 0);

  tstop();
}

}} // namespace psi::cis
