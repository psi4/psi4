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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
/*
**  CCLAMBDA: Program to calculate the coupled-cluster lambda vector.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/wavefunction.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#include "globals.h"
#include "cclambda.h"

namespace psi { namespace cclambda {

void init_io(void);
void title(void);
void get_moinfo(std::shared_ptr<Wavefunction> wfn);
void get_params(Options& options);
void cleanup(void);
void init_amps(struct L_Params L_params);
double pseudoenergy(struct L_Params L_params);
void exit_io(void);
void G_build(int L_irr);
void L1_build(struct L_Params L_params);
void L2_build(struct L_Params L_params);
void sort_amps(int L_irr);
void Lsave(int L_irr);
void Lnorm(struct L_Params L_params);
void Lmag(void);
void update(void);
int converged(int L_irr);
void diis(int iter, int L_irr);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void denom(struct L_Params);
void overlap(int L_irr);
void overlap_LAMPS(struct L_Params L_params);
void Lsave_index(struct L_Params L_params);
void Lamp_write(struct L_Params L_params);
void check_ortho(struct L_Params *pL_params);
void projections(struct L_Params *pL_params);
void L_zero(int irrep);
void c_clean(dpdfile2 *LIA, dpdfile2 *Lia, dpdbuf4 *LIJAB, dpdbuf4 *Lijab, dpdbuf4 *LIjAb);
void L_clean(struct L_Params pL_params);
void zeta_norm(struct L_Params pL_params);
void spinad_amps(void);
void status(const char *, std::string);
void hbar_extra(void);
void ortho_Rs(struct L_Params *pL_params, int current_L);

void cc2_L1_build(struct L_Params L_params);
void cc2_L2_build(struct L_Params L_params);
void cc2_Gai_build(int L_irr);
void cc2_hbar_extra(void);

void cc3_t3z(void);
void cc3_t3x(void);
void cc3_l3l2(void);
void cc3_l3l1(void);

void local_init(void);
void local_done(void);

}} //namespace psi::cclambda

// Forward declaration to call cctriples
namespace psi { namespace cctriples {
PsiReturnType cctriples(std::shared_ptr<Wavefunction> ref_wfn, Options &options);
}}

namespace psi { namespace cclambda {

CCLambdaWavefunction::CCLambdaWavefunction(std::shared_ptr<Wavefunction>
reference_wavefunction, Options &options)
    : Wavefunction(options)
{
    set_reference_wavefunction(reference_wavefunction);
    psio_ = _default_psio_lib_;
    init();
}

CCLambdaWavefunction::~CCLambdaWavefunction()
{

}

void CCLambdaWavefunction::init()
{
    shallow_copy(reference_wavefunction_);
}

double CCLambdaWavefunction::compute_energy()
{
    energy_ = 0.0;
    int done=0, i, root_L_irr;
    int **cachelist, *cachefiles;
    dpdfile2 L1;

    init_io();
    title();
    moinfo.iter=0;
    get_moinfo(reference_wavefunction_);
    get_params(options_);

    /* throw any existing CC_LAMBDA, CC_DENOM away */
    /* Do this only if we're not running an analytic gradient on the
       ground state. Keeping the files around should allow us to
       restart from old Lambda amplitudes. -TDC, 11/2007 */
    if(!(params.dertype==1 && !cc_excited(params.wfn))) {
      outfile->Printf( "\tDeleting old CC_LAMBDA data.\n");
      psio_close(PSIF_CC_LAMBDA,0);
      psio_open(PSIF_CC_LAMBDA,PSIO_OPEN_NEW);
      psio_close(PSIF_CC_DENOM,0);
      psio_open(PSIF_CC_DENOM,PSIO_OPEN_NEW);
    }

    cachefiles = init_int_array(PSIO_MAXUNIT);

    if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

      cachelist = cacheprep_rhf(params.cachelev, cachefiles);

      std::vector<int*> spaces;
      spaces.push_back(moinfo.occpi);
      spaces.push_back(moinfo.occ_sym);
      spaces.push_back(moinfo.virtpi);
      spaces.push_back(moinfo.vir_sym);
      dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 2, spaces);

      if(params.aobasis) { /* Set up new DPD for AO-basis algorithm */
          std::vector<int*> aospaces;
          aospaces.push_back(moinfo.occpi);
          aospaces.push_back(moinfo.occ_sym);
          aospaces.push_back(moinfo.sopi);
          aospaces.push_back(moinfo.sosym);
          dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 2, aospaces);
          dpd_set_default(0);
      }

    }
    else if(params.ref == 2) { /** UHF **/

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

      if(params.aobasis) { /* Set up new DPD's for AO-basis algorithm */
          std::vector<int*> aospaces;
          aospaces.push_back(moinfo.aoccpi);
          aospaces.push_back(moinfo.aocc_sym);
          aospaces.push_back(moinfo.sopi);
          aospaces.push_back(moinfo.sosym);
          aospaces.push_back(moinfo.boccpi);
          aospaces.push_back(moinfo.bocc_sym);
          aospaces.push_back(moinfo.sopi);
          aospaces.push_back(moinfo.sosym);
          dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 4, aospaces);
          dpd_set_default(0);
      }
    }

    if(params.local) local_init();

    if(params.ref == 0) {
      if (params.wfn == "CC2" || params.wfn == "EOM_CC2")
        cc2_hbar_extra();
      else
        hbar_extra();
    }

    /* CC3: Z-build */
    if(params.wfn == "CC3") cc3_t3z();

    for (i=0; i<params.nstates; ++i) {

      /* delete and reopen intermediate files */
      psio_close(PSIF_CC_TMP,0); psio_close(PSIF_CC_TMP0,0);
      psio_close(PSIF_CC_TMP1,0); psio_close(PSIF_CC_TMP2,0);
      psio_open(PSIF_CC_TMP,0); psio_open(PSIF_CC_TMP0,0);
      psio_open(PSIF_CC_TMP1,0); psio_open(PSIF_CC_TMP2,0);
      /* Keep the old lambda amps if this is a ground-state geomopt */
      if(!(params.dertype==1 && !cc_excited(params.wfn))) {
        psio_close(PSIF_CC_LAMBDA,0);
        psio_open(PSIF_CC_LAMBDA,PSIO_OPEN_NEW);
        psio_close(PSIF_CC_DENOM,0); /* aren't these recomputed anyway - perhaps should always delete? */
        psio_open(PSIF_CC_DENOM,PSIO_OPEN_NEW);
      }

      outfile->Printf("\tSymmetry of left-hand state: %s\n",
              moinfo.labels[ moinfo.sym^(pL_params[i].irrep) ]);
      outfile->Printf("\tSymmetry of left-hand eigenvector: %s\n",
              moinfo.labels[(pL_params[i].irrep)]);

      denom(pL_params[i]); /* uses L_params.cceom_energy for excited states */
      init_amps(pL_params[i]); /* uses denominators for initial zeta guess */

      outfile->Printf( "\n\t          Solving Lambda Equations\n");
      outfile->Printf( "\t          ------------------------\n");
      outfile->Printf( "\tIter     PseudoEnergy or Norm         RMS  \n");
      outfile->Printf( "\t----     ---------------------     --------\n");

      moinfo.lcc = pseudoenergy(pL_params[i]);
      update();

      for(moinfo.iter=1 ; moinfo.iter <= params.maxiter; moinfo.iter++) {
        sort_amps(pL_params[i].irrep);

        /* must zero New L before adding RHS */
        L_zero(pL_params[i].irrep);

        if(params.wfn == "CC3") cc3_t3x();

        if(params.wfn == "CC2" || params.wfn == "EOM_CC2") {

      cc2_Gai_build(pL_params[i].irrep);
      cc2_L1_build(pL_params[i]);
      if(params.print & 2) status("L1 amplitudes", "outfile");
      cc2_L2_build(pL_params[i]);

        }
        else {
      G_build(pL_params[i].irrep);
      L1_build(pL_params[i]);
      if(params.print & 2) status("L1 amplitudes", "outfile");
      L2_build(pL_params[i]);

      if(params.wfn == "CC3") {
        cc3_l3l2();
        cc3_l3l1();
      }
        }

        if (params.ref == 1) L_clean(pL_params[i]);
        if (params.nstates > 2) ortho_Rs(pL_params, i);

        if(converged(pL_params[i].irrep)) {
          done = 1;  /* Boolean for convergence */
          Lsave(pL_params[i].irrep); /* copy "New L" to "L" */
          moinfo.lcc = pseudoenergy(pL_params[i]);
          update();
          if (!pL_params[i].ground && !params.zeta) {
            Lnorm(pL_params[i]); /* normalize against R */
          }
          Lsave_index(pL_params[i]); /* save Ls with indices in LAMPS */
          Lamp_write(pL_params[i]); /* write out largest  Ls */
      /* sort_amps(); to be done by later functions */
          outfile->Printf( "\n\tIterations converged.\n");

          moinfo.iter = 0;
          break;
        }

        if(params.diis) diis(moinfo.iter, pL_params[i].irrep);
        Lsave(pL_params[i].irrep);
        moinfo.lcc = pseudoenergy(pL_params[i]);
        update();
      }
      outfile->Printf( "\n");
      if(!done) {
        outfile->Printf( "\t ** Lambda not converged to %2.1e ** \n",
            params.convergence);

        dpd_close(0);
        cleanup();
        exit_io();
        throw PsiException("cclambda: error", __FILE__, __LINE__);
      }
      if (pL_params[i].ground)
        overlap(pL_params[i].irrep);
    }

    if (params.zeta) {
      zeta_norm(pL_params[0]);
    }
    else if (params.nstates > 1) { /* some excited states are present */
      check_ortho(pL_params);
      projections(pL_params);
    }

    if(params.local) local_done();

    dpd_close(0);

    if(params.ref == 2) cachedone_uhf(cachelist);
    else cachedone_rhf(cachelist);
    free(cachefiles);

    cleanup();
    exit_io();

    if ((options_.get_str("WFN") == "CCSD_AT")) {

      // Run cctriples
      if (psi::cctriples::cctriples(reference_wavefunction_, options_) == Success)
          energy_ = Process::environment.globals["CURRENT ENERGY"];
      else
          energy_ = 0.0;
    }

    return energy_;
}

// must be fixed with options later for excited states
void init_io(void)
{
  int i, num_unparsed;
  char *lbl, *argv_unparsed[100];

  params.all=0;    /* do all Ls including ground state */
  params.zeta=0; /* only do ground-state L */
/*
  for (i=1, num_unparsed=0; i<argc; ++i) {
    if (!strcmp(argv[i],"--all")) {
      params.all = 1;
    }
    else if (!strcmp(argv[i],"--zeta")) {
      params.zeta = 1;
    }
    else {
      argv_unparsed[num_unparsed++] = argv[i];
    }
  }
*/
  tstart();

  for(i=PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i,1);
}

void title(void)
{
  outfile->Printf( "\n");
  outfile->Printf( "\t\t\t**************************\n");
  outfile->Printf( "\t\t\t*        CCLAMBDA        *\n");
  outfile->Printf( "\t\t\t**************************\n");
  outfile->Printf( "\n");
}

void exit_io(void)
{
  int i;

  for(i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) {
    psio_close(i,0);
    psio_open(i,PSIO_OPEN_NEW);
  }
  psio_close(PSIF_CC_DENOM,0);
  psio_open(PSIF_CC_DENOM,PSIO_OPEN_NEW);

  /* Close all dpd data files here */
  for(i=PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio_close(i,1);
  for(i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio_close(i,0); /* delete CC_TMP files */
  for(i=PSIF_CC_TMP11+1; i <= PSIF_CC_MAX; i++) psio_close(i,1);

  tstop();
}

/* put copies of L for excited states in LAMPS with irrep and index label */
void Lsave_index(struct L_Params L_params) {
  int L_irr;
  dpdfile2 L1;
  dpdbuf4 L2, LIjAb, LIjbA;
  char *L1A_lbl, *L1B_lbl, *L2AA_lbl, *L2BB_lbl, *L2AB_lbl, *L2RHF_lbl, lbl[32];
  L1A_lbl = L_params.L1A_lbl;
  L1B_lbl = L_params.L1B_lbl;
  L2AA_lbl = L_params.L2AA_lbl;
  L2BB_lbl = L_params.L2BB_lbl;
  L2AB_lbl = L_params.L2AB_lbl;
  L2RHF_lbl = L_params.L2RHF_lbl;
  L_irr = L_params.irrep;

  if(params.ref == 0 || params.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1, PSIF_CC_LAMPS, L1A_lbl);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
    global_dpd_->file2_copy(&L1, PSIF_CC_LAMPS, L1B_lbl);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMPS, L2AA_lbl);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMPS, L2BB_lbl);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMPS, L2AB_lbl);
    global_dpd_->buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1, PSIF_CC_LAMPS, L1A_lbl);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
    global_dpd_->file2_copy(&L1, PSIF_CC_LAMPS, L1B_lbl);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMPS, L2AA_lbl);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMPS, L2BB_lbl);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMPS, L2AB_lbl);
    global_dpd_->buf4_close(&L2);
  }

  if (params.ref == 0) { /** RHF for those codes that can use them **/
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMPS, L_irr, 0, 5, 0, 5, 0, L2AB_lbl);
    global_dpd_->buf4_sort(&LIjAb, PSIF_CC_TMP, pqsr, 0, 5, "LIjbA");
    global_dpd_->buf4_copy(&LIjAb, PSIF_CC_LAMPS, L2RHF_lbl);
    global_dpd_->buf4_close(&LIjAb);

    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMPS, L_irr, 0, 5, 0, 5, 0, L2RHF_lbl);
    global_dpd_->buf4_scm(&LIjAb, 2.0);
    global_dpd_->buf4_init(&LIjbA, PSIF_CC_TMP, L_irr, 0, 5, 0, 5, 0, "LIjbA");
    global_dpd_->buf4_axpy(&LIjbA, &LIjAb, -1.0);
    global_dpd_->buf4_close(&LIjbA);
    global_dpd_->buf4_close(&LIjAb);
  }
  return;
}

void L_zero(int L_irr) {
  dpdfile2 LIA, Lia;
  dpdbuf4 LIJAB, Lijab, LIjAb;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_scm(&LIA, 0.0);
    global_dpd_->file2_close(&LIA);
  }
  else if(params.ref == 1) { /** RHF/ROHF **/
    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
    global_dpd_->file2_scm(&LIA, 0.0);
    global_dpd_->file2_scm(&Lia, 0.0);
    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_close(&Lia);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");
    global_dpd_->file2_scm(&LIA, 0.0);
    global_dpd_->file2_scm(&Lia, 0.0);
    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_close(&Lia);
  }

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_scm(&LIjAb, 0.0);
    global_dpd_->buf4_close(&LIjAb);
  }
  else if (params.ref == 1 ) { /** ROHF **/
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_scm(&LIJAB, 0.0);
    global_dpd_->buf4_close(&LIJAB);
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_scm(&Lijab, 0.0);
    global_dpd_->buf4_close(&Lijab);
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_scm(&LIjAb, 0.0);
    global_dpd_->buf4_close(&LIjAb);
  }
  else { /** UHF **/
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_scm(&LIJAB, 0.0);
    global_dpd_->buf4_close(&LIJAB);
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    global_dpd_->buf4_scm(&Lijab, 0.0);
    global_dpd_->buf4_close(&Lijab);
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    global_dpd_->buf4_scm(&LIjAb, 0.0);
    global_dpd_->buf4_close(&LIjAb);
  }
}


/* Cleaning out L vectors for open-shell cases  */
void L_clean(struct L_Params L_params) {
  int L_irr, i;
  dpdfile2 LIA, Lia;
  dpdbuf4 LIJAB, Lijab, LIjAb;
  char lbl[80];

  L_irr = L_params.irrep;

  global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
  global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
  global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
  global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
  global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

  c_clean(&LIA, &Lia, &LIJAB, &Lijab, &LIjAb);

  global_dpd_->file2_close(&LIA);
  global_dpd_->file2_close(&Lia);
  global_dpd_->buf4_close(&LIJAB);
  global_dpd_->buf4_close(&Lijab);
  global_dpd_->buf4_close(&LIjAb);
}

void zeta_norm(struct L_Params L_params) {
  int Z_irr, i;
  dpdfile2 ZIA, Zia;
  dpdbuf4 ZIJAB, Zijab, ZIjAb;
  double tval;
  Z_irr = L_params.irrep;

  if (params.ref == 0 || params.ref == 1) {
    global_dpd_->file2_init(&ZIA, PSIF_CC_LAMPS, Z_irr, 0, 1, "ZIA");
    tval = global_dpd_->file2_dot_self(&ZIA);
    global_dpd_->file2_close(&ZIA);
    global_dpd_->file2_init(&Zia, PSIF_CC_LAMPS, Z_irr, 0, 1, "Zia");
    tval += global_dpd_->file2_dot_self(&Zia);
    global_dpd_->file2_close(&Zia);
    global_dpd_->buf4_init(&ZIJAB, PSIF_CC_LAMPS, Z_irr, 2, 7, 2, 7, 0, "ZIJAB");
    tval += global_dpd_->buf4_dot_self(&ZIJAB);
    global_dpd_->buf4_close(&ZIJAB);
    global_dpd_->buf4_init(&Zijab, PSIF_CC_LAMPS, Z_irr, 2, 7, 2, 7, 0, "Zijab");
    tval += global_dpd_->buf4_dot_self(&Zijab);
    global_dpd_->buf4_close(&Zijab);
    global_dpd_->buf4_init(&ZIjAb, PSIF_CC_LAMPS, Z_irr, 0, 5, 0, 5, 0, "ZIjAb");
    tval += global_dpd_->buf4_dot_self(&ZIjAb);
    global_dpd_->buf4_close(&ZIjAb);
  }
  else { /* UHF */
    global_dpd_->file2_init(&ZIA, PSIF_CC_LAMPS, Z_irr, 0, 1, "ZIA");
    tval = global_dpd_->file2_dot_self(&ZIA);
    global_dpd_->file2_close(&ZIA);
    global_dpd_->file2_init(&Zia, PSIF_CC_LAMPS, Z_irr, 2, 3, "Zia");
    tval += global_dpd_->file2_dot_self(&Zia);
    global_dpd_->file2_close(&Zia);
    global_dpd_->buf4_init(&ZIJAB, PSIF_CC_LAMPS, Z_irr, 2, 7, 2, 7, 0, "ZIJAB");
    tval += global_dpd_->buf4_dot_self(&ZIJAB);
    global_dpd_->buf4_close(&ZIJAB);
    global_dpd_->buf4_init(&Zijab, PSIF_CC_LAMPS, Z_irr, 12, 17, 12, 17, 0, "Zijab");
    tval += global_dpd_->buf4_dot_self(&Zijab);
    global_dpd_->buf4_close(&Zijab);
    global_dpd_->buf4_init(&ZIjAb, PSIF_CC_LAMPS, Z_irr, 22, 28, 22, 28, 0, "ZIjAb");
    tval += global_dpd_->buf4_dot_self(&ZIjAb);
    global_dpd_->buf4_close(&ZIjAb);
  }
  outfile->Printf("Norm of Zeta: %20.15lf\n", sqrt(tval) );
  return;
}

}} // namespace psi::cclambda
