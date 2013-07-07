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

/*! \defgroup RESPONSE response: Compute various response properties */

/*!
** \file
** \ingroup RESPONSE
** \brief Module to compute various response properties
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
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"

namespace psi { namespace response {

/* Max length of ioff array */
#define IOFF_MAX 32641

/* Function prototypes */
void init_io(int argc, char *argv[]);
void init_ioff(void);
void get_moinfo(void);
void get_params(void);
void cleanup(void);
void exit_io(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void build_A_RHF(void);
void build_B_RHF(void);
void invert_RPA_RHF(double omega);
void polar(void);
void optrot(void);
void dipquad(void);


int response(int argc, char *argv[])
{
  int **cachelist, *cachefiles;

  init_io(argc, argv);
  init_ioff();
  get_moinfo();
  get_params();

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

  if(params.ref == 0) {
    build_A_RHF();
    build_B_RHF();
    invert_RPA_RHF(params.omega[0]);
    polar();
    dipquad();
    optrot();
  }

  dpd_close(0);
  cleanup();
  exit_io();
  exit(0);
}


void init_io(int argc, char *argv[])
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
  psio_done();
  tstop();
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}

}} // namespace psi::response
