/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \defgroup STABLE stable: Perform MO Stability Analysis */

/*! 
** \file
** \ingroup STABLE
** \brief Module to calculate eigenvalues and eigenvectors of the molecular
**   orbital Hessian
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include "Params.h"
#include "MOInfo.h"
#include "globals.h"
#include "liboptions/liboptions.h"

namespace psi { namespace stable {

/* Max length of ioff array */
#define IOFF_MAX 32641

/* Function prototypes */
void init_io();
void init_ioff(void);
void get_moinfo(void);
void get_params(Options &options);
void cleanup(void);
void exit_io(void);
void build_A_RHF(void);
void build_A_ROHF(void);
void build_A_UHF(void);
void diag_A_RHF(void);
void diag_A_ROHF(void);
void diag_A_UHF(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void print_evals(double **evals, int *rank);



PsiReturnType stability(Options & options)
{
  int **cachelist, *cachefiles;

  init_io();
  init_ioff();
  get_params(options);
  get_moinfo();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /*** UHF references ***/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    std::vector<int*> spaces;
    spaces.push_back(moinfo.aoccpi);
    spaces.push_back(moinfo.aocc_sym);
    spaces.push_back(moinfo.avirtpi);
    spaces.push_back(moinfo.avir_sym);
    spaces.push_back(moinfo.boccpi);
    spaces.push_back(moinfo.bocc_sym);
    spaces.push_back(moinfo.bvirtpi);
    spaces.push_back(moinfo.bvir_sym);
    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 4, spaces);
  } 
  else { /*** RHF/ROHF references ***/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    std::vector<int*> spaces;
    spaces.push_back(moinfo.occpi);
    spaces.push_back(moinfo.occ_sym);
    spaces.push_back(moinfo.virtpi);
    spaces.push_back(moinfo.vir_sym);
    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 2, spaces);
  }

  /* Final list of evals for pretty output */
  moinfo.A_evals = block_matrix(moinfo.nirreps,5);
  moinfo.rank = init_int_array(moinfo.nirreps);
  if(params.ref == 0)
    moinfo.A_triplet_evals = block_matrix(moinfo.nirreps, 5);

  if(params.ref == 0) {
    build_A_RHF();
    diag_A_RHF();
  }
  else if(params.ref == 1) {
    build_A_ROHF();
    diag_A_ROHF();
  }
  else if(params.ref == 2) {
    build_A_UHF();
    diag_A_UHF();
  }

  /* Print the eigenvalues in a nice format */
  if(params.ref == 0) {
    outfile->Printf( "\n\t     RHF->RHF stability eigenvalues:\n");
    outfile->Printf(   "\t     -------------------------------\n");
    print_evals(moinfo.A_evals, moinfo.rank);

    outfile->Printf( "\n\t     RHF->UHF stability eigenvalues:\n");
    outfile->Printf(   "\t     -------------------------------\n");
    print_evals(moinfo.A_triplet_evals, moinfo.rank);
  }
  else if(params.ref == 1) {
    outfile->Printf( "\n\t     ROHF->ROHF stability eigenvalues:\n");
    outfile->Printf(   "\t     ---------------------------------\n");
    print_evals(moinfo.A_evals, moinfo.rank);
  }
  else if(params.ref == 2) {
    outfile->Printf( "\n\t     UHF->UHF stability eigenvalues:\n");
    outfile->Printf(   "\t     -------------------------------\n");
    print_evals(moinfo.A_evals, moinfo.rank);
  }

  if(params.ref == 0) free_block(moinfo.A_triplet_evals);
  free_block(moinfo.A_evals);
  free(moinfo.rank);

  dpd_close(0);
  cleanup();
  exit_io();
  exit(PSI_RETURN_SUCCESS);
}



void init_io()
{
  int i;


  tstart();

  psio_open(PSIF_MO_HESS,0);
  psio_open(PSIF_CC_INFO, PSIO_OPEN_OLD);
  psio_open(PSIF_CC_OEI, PSIO_OPEN_OLD);
  psio_open(PSIF_CC_CINTS, PSIO_OPEN_OLD);
  psio_open(PSIF_CC_DINTS, PSIO_OPEN_OLD);
  psio_open(PSIF_CC_TMP0, PSIO_OPEN_NEW);
}

void exit_io(void)
{
  int i;
 
  psio_close(PSIF_MO_HESS,1);
  psio_close(PSIF_CC_INFO, 1);
  psio_close(PSIF_CC_OEI, 1);
  psio_close(PSIF_CC_CINTS, 1);
  psio_close(PSIF_CC_DINTS, 1);
  psio_close(PSIF_CC_TMP0, 1);

  tstop();
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}

}} // namespace psi::stable
