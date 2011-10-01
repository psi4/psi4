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
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#include "globals.h"

namespace psi { namespace cclambda {

void init_io(void);
void title(void);
void get_moinfo(void);
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
void status(const char *, FILE *);
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

PsiReturnType cclambda(Options& options)
{
  int done=0, i, root_L_irr;
  int **cachelist, *cachefiles;
  dpdfile2 L1;

  init_io();
  title();
  moinfo.iter=0;
  get_moinfo();
  get_params(options);

  /* throw any existing CC_LAMBDA, CC_DENOM away */
  /* Do this only if we're not running an analytic gradient on the
     ground state. Keeping the files around should allow us to
     restart from old Lambda amplitudes. -TDC, 11/2007 */
  if(!(params.dertype==1 && !cc_excited(params.wfn))) {
    fprintf(outfile, "\tDeleting old CC_LAMBDA data.\n");
    psio_close(CC_LAMBDA,0);
    psio_open(CC_LAMBDA,PSIO_OPEN_NEW);
    psio_close(CC_DENOM,0);
    psio_open(CC_DENOM,PSIO_OPEN_NEW);
  }

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
         2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);

    if(params.aobasis) { /* Set up new DPD for AO-basis algorithm */
      dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
           2, moinfo.occpi, moinfo.occ_sym, moinfo.sopi, moinfo.sosym);
      dpd_set_default(0);
    }

  }
  else if(params.ref == 2) { /** UHF **/

    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles,
         cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
         moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);

    if(params.aobasis) { /* Set up new DPD's for AO-basis algorithm */
      dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
               4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.sopi, moinfo.sosym,
               moinfo.boccpi, moinfo.bocc_sym, moinfo.sopi, moinfo.sosym);
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
    psio_close(CC_TMP,0); psio_close(CC_TMP0,0);
    psio_close(CC_TMP1,0); psio_close(CC_TMP2,0);
    psio_open(CC_TMP,0); psio_open(CC_TMP0,0);
    psio_open(CC_TMP1,0); psio_open(CC_TMP2,0);
    /* Keep the old lambda amps if this is a ground-state geomopt */
    if(!(params.dertype==1 && !cc_excited(params.wfn))) {
      psio_close(CC_LAMBDA,0);
      psio_open(CC_LAMBDA,PSIO_OPEN_NEW);
      psio_close(CC_DENOM,0); /* aren't these recomputed anyway - perhaps should always delete? */
      psio_open(CC_DENOM,PSIO_OPEN_NEW);
    }

    fprintf(outfile,"\tSymmetry of left-hand state: %s\n",
            moinfo.labels[ moinfo.sym^(pL_params[i].irrep) ]);
    fprintf(outfile,"\tSymmetry of left-hand eigenvector: %s\n",
            moinfo.labels[(pL_params[i].irrep)]);

    denom(pL_params[i]); /* uses L_params.cceom_energy for excited states */
    init_amps(pL_params[i]); /* uses denominators for initial zeta guess */

    fprintf(outfile, "\n\t          Solving Lambda Equations\n");
    fprintf(outfile, "\t          ------------------------\n");
    fprintf(outfile, "\tIter     PseudoEnergy or Norm         RMS  \n");
    fprintf(outfile, "\t----     ---------------------     --------\n");

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
    if(params.print & 2) status("L1 amplitudes", outfile);
    cc2_L2_build(pL_params[i]);

      }
      else {
    G_build(pL_params[i].irrep);
    L1_build(pL_params[i]);
    if(params.print & 2) status("L1 amplitudes", outfile);
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
        fprintf(outfile, "\n\tIterations converged.\n");
        fflush(outfile);
        moinfo.iter = 0;
        break;
      }

      if(params.diis) diis(moinfo.iter, pL_params[i].irrep);
      Lsave(pL_params[i].irrep);
      moinfo.lcc = pseudoenergy(pL_params[i]);
      update();
    }
    fprintf(outfile, "\n");
    if(!done) {
      fprintf(outfile, "\t ** Lambda not converged to %2.1e ** \n",
          params.convergence);
      fflush(outfile);
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
  return Success;
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

  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*        CCLAMBDA        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

void exit_io(void)
{
  int i;

  for(i=CC_TMP; i <= CC_TMP11; i++) {
    psio_close(i,0);
    psio_open(i,PSIO_OPEN_NEW);
  }
  psio_close(CC_DENOM,0);
  psio_open(CC_DENOM,PSIO_OPEN_NEW);

  /* Close all dpd data files here */
  for(i=CC_MIN; i < CC_TMP; i++) psio_close(i,1);
  for(i=CC_TMP; i <= CC_TMP11; i++) psio_close(i,0); /* delete CC_TMP files */
  for(i=CC_TMP11+1; i <= CC_MAX; i++) psio_close(i,1);

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
    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_LAMPS, L1A_lbl);
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 0, 1, "Lia");
    dpd_file2_copy(&L1, CC_LAMPS, L1B_lbl);
    dpd_file2_close(&L1);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_copy(&L2, CC_LAMPS, L2AA_lbl);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_copy(&L2, CC_LAMPS, L2BB_lbl);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_copy(&L2, CC_LAMPS, L2AB_lbl);
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_LAMPS, L1A_lbl);
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 2, 3, "Lia");
    dpd_file2_copy(&L1, CC_LAMPS, L1B_lbl);
    dpd_file2_close(&L1);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_copy(&L2, CC_LAMPS, L2AA_lbl);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_copy(&L2, CC_LAMPS, L2BB_lbl);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_copy(&L2, CC_LAMPS, L2AB_lbl);
    dpd_buf4_close(&L2);
  }

  if (params.ref == 0) { /** RHF for those codes that can use them **/
    dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, L2AB_lbl);
    dpd_buf4_sort(&LIjAb, CC_TMP, pqsr, 0, 5, "LIjbA");
    dpd_buf4_copy(&LIjAb, CC_LAMPS, L2RHF_lbl);
    dpd_buf4_close(&LIjAb);

    dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, L2RHF_lbl);
    dpd_buf4_scm(&LIjAb, 2.0);
    dpd_buf4_init(&LIjbA, CC_TMP, L_irr, 0, 5, 0, 5, 0, "LIjbA");
    dpd_buf4_axpy(&LIjbA, &LIjAb, -1.0);
    dpd_buf4_close(&LIjbA);
    dpd_buf4_close(&LIjAb);
  }
  return;
}

void L_zero(int L_irr) {
  dpdfile2 LIA, Lia;
  dpdbuf4 LIJAB, Lijab, LIjAb;

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_scm(&LIA, 0.0);
    dpd_file2_close(&LIA);
  }
  else if(params.ref == 1) { /** RHF/ROHF **/
    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "New Lia");
    dpd_file2_scm(&LIA, 0.0);
    dpd_file2_scm(&Lia, 0.0);
    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "New Lia");
    dpd_file2_scm(&LIA, 0.0);
    dpd_file2_scm(&Lia, 0.0);
    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);
  }

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_scm(&LIjAb, 0.0);
    dpd_buf4_close(&LIjAb);
  }
  else if (params.ref == 1 ) { /** ROHF **/
    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_scm(&LIJAB, 0.0);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_scm(&Lijab, 0.0);
    dpd_buf4_close(&Lijab);
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_scm(&LIjAb, 0.0);
    dpd_buf4_close(&LIjAb);
  }
  else { /** UHF **/
    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_scm(&LIJAB, 0.0);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_scm(&Lijab, 0.0);
    dpd_buf4_close(&Lijab);
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_scm(&LIjAb, 0.0);
    dpd_buf4_close(&LIjAb);
  }
}


/* Cleaning out L vectors for open-shell cases  */
void L_clean(struct L_Params L_params) {
  int L_irr, i;
  dpdfile2 LIA, Lia;
  dpdbuf4 LIJAB, Lijab, LIjAb;
  char lbl[80];

  L_irr = L_params.irrep;

  dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
  dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "New Lia");
  dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
  dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
  dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

  c_clean(&LIA, &Lia, &LIJAB, &Lijab, &LIjAb);

  dpd_file2_close(&LIA);
  dpd_file2_close(&Lia);
  dpd_buf4_close(&LIJAB);
  dpd_buf4_close(&Lijab);
  dpd_buf4_close(&LIjAb);
}

void zeta_norm(struct L_Params L_params) {
  int Z_irr, i;
  dpdfile2 ZIA, Zia;
  dpdbuf4 ZIJAB, Zijab, ZIjAb;
  double tval;
  Z_irr = L_params.irrep;

  if (params.ref == 0 || params.ref == 1) {
    dpd_file2_init(&ZIA, CC_LAMPS, Z_irr, 0, 1, "ZIA");
    tval = dpd_file2_dot_self(&ZIA);
    dpd_file2_close(&ZIA);
    dpd_file2_init(&Zia, CC_LAMPS, Z_irr, 0, 1, "Zia");
    tval += dpd_file2_dot_self(&Zia);
    dpd_file2_close(&Zia);
    dpd_buf4_init(&ZIJAB, CC_LAMPS, Z_irr, 2, 7, 2, 7, 0, "ZIJAB");
    tval += dpd_buf4_dot_self(&ZIJAB);
    dpd_buf4_close(&ZIJAB);
    dpd_buf4_init(&Zijab, CC_LAMPS, Z_irr, 2, 7, 2, 7, 0, "Zijab");
    tval += dpd_buf4_dot_self(&Zijab);
    dpd_buf4_close(&Zijab);
    dpd_buf4_init(&ZIjAb, CC_LAMPS, Z_irr, 0, 5, 0, 5, 0, "ZIjAb");
    tval += dpd_buf4_dot_self(&ZIjAb);
    dpd_buf4_close(&ZIjAb);
  }
  else { /* UHF */
    dpd_file2_init(&ZIA, CC_LAMPS, Z_irr, 0, 1, "ZIA");
    tval = dpd_file2_dot_self(&ZIA);
    dpd_file2_close(&ZIA);
    dpd_file2_init(&Zia, CC_LAMPS, Z_irr, 2, 3, "Zia");
    tval += dpd_file2_dot_self(&Zia);
    dpd_file2_close(&Zia);
    dpd_buf4_init(&ZIJAB, CC_LAMPS, Z_irr, 2, 7, 2, 7, 0, "ZIJAB");
    tval += dpd_buf4_dot_self(&ZIJAB);
    dpd_buf4_close(&ZIJAB);
    dpd_buf4_init(&Zijab, CC_LAMPS, Z_irr, 12, 17, 12, 17, 0, "Zijab");
    tval += dpd_buf4_dot_self(&Zijab);
    dpd_buf4_close(&Zijab);
    dpd_buf4_init(&ZIjAb, CC_LAMPS, Z_irr, 22, 28, 22, 28, 0, "ZIjAb");
    tval += dpd_buf4_dot_self(&ZIjAb);
    dpd_buf4_close(&ZIjAb);
  }
  fprintf(outfile,"Norm of Zeta: %20.15lf\n", sqrt(tval) );
  return;
}

}} // namespace psi::cclambda
