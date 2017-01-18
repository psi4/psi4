/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

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
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/physconst.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/wavefunction.h"
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
void get_moinfo(std::shared_ptr<Wavefunction>);
void get_params(std::shared_ptr<Wavefunction>, Options&);
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
void optrot(std::shared_ptr<Molecule> molecule);
void roa(void);

void preppert(std::shared_ptr<BasisSet> primary);

int ccresponse(std::shared_ptr<Wavefunction> ref_wfn, Options &options)
{
  int **cachelist, *cachefiles;

  init_io();
  init_ioff();
  title();
  get_moinfo(ref_wfn);
  get_params(ref_wfn, options);

  timer_on("ccresponse");

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

  if(params.local) local_init();

  if (params.wfn == "CC2") {
    cc2_hbar_extra();
  }
  else {
    hbar_extra();
 }

  sort_lamps(); /* should be removed sometime - provided by cclambda */
  if(params.wfn != "CC2") lambda_residuals(); /* don't do this for CC2 */

  preppert(ref_wfn->basisset());

  if(params.prop == "POLARIZABILITY") polar();
  if(params.prop == "ROTATION") optrot(ref_wfn->molecule());
  if(params.prop == "ROA_TENSOR") roa();

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

  for(i=PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i, 1);

  /* Clear out DIIS TOC Entries */
  psio_close(PSIF_CC_DIIS_AMP, 0);
  psio_close(PSIF_CC_DIIS_ERR, 0);

  psio_open(PSIF_CC_DIIS_AMP, 0);
  psio_open(PSIF_CC_DIIS_ERR, 0);
}

void title(void)
{
  outfile->Printf( "\t\t\t**************************\n");
  outfile->Printf( "\t\t\t*                        *\n");
  outfile->Printf( "\t\t\t*       ccresponse       *\n");
  outfile->Printf( "\t\t\t*                        *\n");
  outfile->Printf( "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;

  /* Close all dpd data files here */
  for(i=PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio_close(i,1);
  for(i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio_close(i,0);  /* get rid of TMP files */
  for(i=PSIF_CC_TMP11+1; i <= PSIF_CC_MAX; i++) psio_close(i,1);

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
