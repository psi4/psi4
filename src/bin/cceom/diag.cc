/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
/*
 *   diag() diagonalized Hbar with Davidson-Liu algorithm to obtain
 *    right-hand eigenvector and eigenvalue
 */

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"
#include <psi4-dec.h>

namespace psi { namespace cceom {
#include <physconst.h>
#include <psifiles.h>

extern void test_dpd();
extern void rzero(int C_irr, int *converged);
extern void rzero_rhf(int C_irr, int *converged);
void init_S1(int index, int irrep);
void init_S2(int index, int irrep);
void init_C1(int i, int C_irr);
void init_C0(int i);
void init_S0(int i);
void init_C2(int index, int irrep);
extern void write_Rs(int C_irr, double *energies, int *converged);
extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
extern double norm_C_full(double C0, dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
extern double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);
extern double norm_C_rhf_full(double C0, dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);
extern void scm_C1(dpdfile2 *CME, dpdfile2 *Cme, double a);
extern void scm_C1_full(double C0, dpdfile2 *CME, dpdfile2 *Cme, double a);
extern void scm_C(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF,
    dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a);
extern void scm_C_full(double C0, dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF,
    dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a);
extern void scm_C2(dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a);
extern void restart(double **alpha, int L, int num, int irrep, int ortho, double **alpha_old, int L_old, int use_alpha_old);
extern void precondition(dpdfile2 *RIA, dpdfile2 *Ria,
    dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb, double eval);
extern void precondition_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, double eval);
void form_diagonal(int irrep);
extern void schmidt_add(dpdfile2 *RIA, dpdfile2 *Ria,
    dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int *numCs, int irrep);
extern void schmidt_add_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, int *numCs, int irrep);
void c_clean(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
void sigmaSS(int index, int irrep);
void sigmaSD(int index, int irrep);
void sigmaDS(int index, int irrep);
void sigmaDD(int index, int irrep);
void sigma00(int index, int irrep);
void sigma0S(int index, int irrep);
void sigma0D(int index, int irrep);
void sigmaS0(int index, int irrep);
void sigmaSS_full(int index, int irrep);
void sigmaD0(int index, int irrep);
void sigmaDS_full(int index, int irrep);
void sigmaDD_full(int index, int irrep);
void sigmaCC3(int i, int C_irr, double omega);
void diagSS(int irrep);
void hbar_extra(void);
void hbar_norms(void);
extern void sort_C(int index, int irrep);

void dgeev_eom(int L, double **G, double *evals, double **alpha);

double local_G1_dot(dpdfile2 *, dpdfile2 *);
double local_G2_dot(dpdbuf4 *, dpdbuf4 *);
void local_filter_T1_nodenom(dpdfile2 *);
void local_filter_T2_nodenom(dpdbuf4 *);
void local_guess(void);

void read_guess(int);
void read_guess_init(void);
extern double norm_C1_rhf(dpdfile2 *C1A);
void cc3_HC1(int i, int C_irr); /* compute [H,C1] */
void cc3_HC1ET1 (int i, int C_irr);
void norm_HC1(int i, int C_irr); /* compute norms for [H,C1] */
extern void restart_with_root(int i, int C_irr);
extern void save_C_ccsd(int i, int C_irr);
extern int follow_root(int L, double **alpha, int C_irr);

void cc2_hbar_extra(void);
void cc2_sigma(int index, int irrep);
void amp_write_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, int length);
void amp_write_UHF(dpdfile2 *, dpdfile2 *, dpdbuf4 *, dpdbuf4 *, dpdbuf4 *, int length);
void amp_write_ROHF(dpdfile2 *, dpdfile2 *, dpdbuf4 *, dpdbuf4 *, dpdbuf4 *, int length);

void diag(void) {
  dpdfile2 CME, CME2, Cme, SIA, Sia, RIA, Ria, DIA, Dia, tIA, tia, LIA, Lia;
  dpdbuf4 CMNEF, Cmnef, CMnEf, SIJAB, Sijab, SIjAb, RIJAB, Rijab, RIjAb, RIjbA;
  dpdbuf4 CMnEf1, CMnfE1, CMnfE, CMneF, C2;
  char lbl[32];
  int num_converged, num_converged_index=0, *converged, keep_going, already_sigma;
  int irrep, numCs, iter, lwork, info, vectors_per_root, nsigma_evaluations=0;
  int get_right_ev = 1, get_left_ev = 0, first_irrep=1;
  int L,h,i,j,k,a,nirreps,errcod,C_irr;
  double norm, tval, **G, *work, *evals_complex, **alpha, **evectors_left;
  double *lambda, *lambda_old, totalE, **G_old, **alpha_old; 
  int num_vecs, cc3_index, num_cc3_restarts = 0, ignore_G_old=0;
  double ra, rb, r2aa, r2bb, r2ab, cc3_eval, cc3_last_converged_eval=0.0, C0, S0, R0;
  int cc3_stage; /* 0=eom_ccsd; 1=eom_cc3 (reuse sigmas), 2=recompute sigma */
  int L_start_iter, L_old;
  char *keyw;


#ifdef TIME_CCEOM
timer_on("HBAR_EXTRA");
#endif
  if (params.wfn == "EOM_CC2") cc2_hbar_extra();  
  else hbar_extra(); /* sort hbar matrix elements for sigma equations */
#ifdef TIME_CCEOM
timer_off("HBAR_EXTRA");
#endif

#ifdef EOM_DEBUG
  hbar_norms();
#endif

//  if(eom_params.guess == "INPUT")
//    read_guess_init();

  if (params.wfn == "EOM_CC3")
    cc3_stage = 0; /* do EOM_CCSD first */

  fprintf(outfile,"Symmetry of ground state: %s\n", moinfo.irr_labs[moinfo.sym]);
  /* loop over symmetry of C's */
  for (C_irr=0; C_irr<moinfo.nirreps; ++C_irr) {

    ignore_G_old = 1;
    already_sigma = 0;
    iter = 0;
    keep_going = 1;
    num_converged = 0;
    if (eom_params.cs_per_irrep[C_irr] == 0) continue;
#ifdef TIME_CCEOM
timer_on("INIT GUESS");
#endif
    fprintf(outfile,"Symmetry of excited state: %s\n", moinfo.irr_labs[moinfo.sym ^ C_irr]);
    fprintf(outfile,"Symmetry of right eigenvector: %s\n",moinfo.irr_labs[C_irr]);
    if (params.eom_ref == 0)
      fprintf(outfile,"Seeking states with multiplicity of %d\n", eom_params.mult);

    /* zero out files between irreps */
    if (!first_irrep) {
      dpd_file2_cache_close();
      dpd_file2_cache_init();
      dpd_file4_cache_close();
      dpd_file4_cache_init();
    }
    for (i=EOM_D; i<=EOM_R; ++i) {
      if (eom_params.restart_eom_cc3 && (i >= EOM_CME) && (i<= EOM_CMnEf)) continue;
      psio_close(i,0);
      psio_open(i,0);
    }
    psio_close(EOM_TMP,0);
    psio_open(EOM_TMP,0);
    first_irrep = 0; /* used to zero out files for new irreps above */

    /* Store approximate diagonal elements of Hbar */
    form_diagonal(C_irr);

    if (eom_params.restart_eom_cc3) {
    }
    else if(params.local) {
      if(eom_params.guess == "DISK") { /* only do this if we don't already have guesses on disk */
	fprintf(outfile, "\n\tUsing C1 vectors on disk as initial guesses.\n");
      }
      else {
	local_guess();
	if(local.do_singles) diagSS(C_irr);
      }
    }
    else {
      if(eom_params.guess == "SINGLES") {
	/* Diagonalize Hbar-SS to obtain initial CME and Cme guess */
	fprintf(outfile,"Obtaining initial guess from singles-singles block of Hbar...");
	diagSS(C_irr);
        if (!eom_params.print_singles) fprintf(outfile,"Done.\n\n");
      }
//      else if(eom_params.guess == "INPUT") {
//        read_guess(C_irr);
//      }
      else if(eom_params.guess == "DISK" && params.ref == 0) {
        fprintf(outfile, "Using C1 vectors on disk as initial guesses.\n");
        /* normalize first guess */
        sprintf(lbl, "%s %d", "CME", 0);
        dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
        norm = norm_C1_rhf(&CME);
        dpd_file2_scm(&CME, 1.0/norm);
        dpd_file2_close(&CME);
        /* reorthoganalize and normalize other guesses */
        for(i=1; i < eom_params.cs_per_irrep[C_irr]; i++) {
          sprintf(lbl, "%s %d", "CME", i);
          dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
          for(j=0; j < i; j++) {
            sprintf(lbl, "%s %d", "CME", j);
            dpd_file2_init(&CME2, EOM_CME, C_irr, 0, 1, lbl);
            norm = 2.0 * dpd_file2_dot(&CME, &CME2);
            dpd_file2_axpy(&CME2, &CME, -1.0*norm, 0);
            dpd_file2_close(&CME2);
          }
          norm = norm_C1_rhf(&CME);
          dpd_file2_scm(&CME, 1.0/norm);
          dpd_file2_close(&CME);
        }
#ifdef EOM_DEBUG
        /* check initial guesses - overlap matrix */
        fprintf(outfile,"Checking overlap of orthogonalized initial guesses\n");
        for(i=0; i < eom_params.cs_per_irrep[C_irr]; i++) {
          sprintf(lbl, "%s %d", "CME", i);
          dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
          for(j=0; j < eom_params.cs_per_irrep[C_irr]; j++) {
            sprintf(lbl, "%s %d", "CME", j);
            dpd_file2_init(&CME2, EOM_CME, C_irr, 0, 1, lbl);
            fprintf(outfile,"C[%d][%d] = %15.10lf\n",i,j, 2.0 * dpd_file2_dot(&CME, &CME2));
            dpd_file2_close(&CME2);
          }
          dpd_file2_close(&CME);
        }
#endif
      }
      else {
        fprintf(outfile, "Invalid initial guess method.\n");
        exit(PSI_RETURN_FAILURE);
      }
    }

#ifdef TIME_CCEOM
timer_off("INIT GUESS");
#endif

#ifdef EOM_DEBUG
    /* printout initial guesses */
    for (i=0;i<eom_params.cs_per_irrep[C_irr]; ++i) {
      sprintf(lbl, "%s %d", "CME", i);
      dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
      dpd_file2_print(&CME,outfile);
      dpd_file2_close(&CME);
      if (params.eom_ref == 1) {
        sprintf(lbl, "%s %d", "Cme", i);
        dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
        dpd_file2_print(&Cme,outfile);
        dpd_file2_close(&Cme);
      }
      else if (params.eom_ref == 2) {
        sprintf(lbl, "%s %d", "Cme", i);
        dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, lbl);
        dpd_file2_print(&Cme,outfile);
        dpd_file2_close(&Cme);
      }
    }
    fflush(outfile);
#endif

    /* Setup and zero initial C2 and S2 vector to go with Hbar_SS */
    for (i=0;i<eom_params.cs_per_irrep[C_irr];++i) {
      /* init_S1(i, C_irr); gets done at first iteration anyway */
      init_S1(i, C_irr);
      if (!eom_params.restart_eom_cc3) init_C2(i, C_irr);
      if (params.full_matrix) { init_C0(i); init_S0(i); }
      /* init_S2(i, C_irr); */
    }

#ifdef EOM_DEBUG
    check_sum("reset", 0, 0); /* reset checksum */ 
#endif

    converged = init_int_array(eom_params.cs_per_irrep[C_irr]);
    lambda_old = init_array(eom_params.cs_per_irrep[C_irr]);
    L = eom_params.cs_per_irrep[C_irr];
    /* allocate G_old just once */
    i = (eom_params.vectors_per_root+1)*eom_params.cs_per_irrep[C_irr];
    G_old = block_matrix(i,i);

    vectors_per_root = eom_params.vectors_per_root; /* used for CCSD */

    while ((keep_going == 1) && (iter < eom_params.max_iter)) {
      fprintf(outfile,"Iter=%-4d L=%-4d", iter+1, L); fflush(outfile);
      keep_going = 0;
      numCs = L_start_iter = L;
      num_converged = 0;

      for (i=already_sigma;i<L;++i) {
        /* Form a zeroed S vector for each C vector
	   SIA and Sia do get overwritten by sigmaSS
	   so this may only be necessary for debugging */
        ++nsigma_evaluations;
        if (params.full_matrix) init_S0(i);
        init_S1(i, C_irr);
        init_S2(i, C_irr);

        sort_C(i, C_irr);

        /* Computing sigma vectors */
#ifdef EOM_DEBUG
        check_sum("reset",0,0);
#endif
#ifdef TIME_CCEOM
        timer_on("SIGMA ALL");
        if (params.wfn == "EOM_CC2") {
          timer_on("sigmacc2"); cc2_sigma(i,C_irr); timer_off("sigmacc2");
        }
        else {
          timer_on("sigmaSS"); sigmaSS(i,C_irr); timer_off("sigmaSS");
          timer_on("sigmaSD"); sigmaSD(i,C_irr); timer_off("sigmaSD");
          timer_on("sigmaDS"); sigmaDS(i,C_irr); timer_off("sigmaDS");
          timer_on("sigmaDD"); sigmaDD(i,C_irr); timer_off("sigmaDD");
          if ( ((params.wfn == "EOM_CC3") && (cc3_stage>0)) || eom_params.restart_eom_cc3) {
            timer_on("cc3_HC1"); cc3_HC1(i,C_irr); timer_off("cc3_HC1");
            timer_on("cc3_HC1ET1"); cc3_HC1ET1(i,C_irr); timer_off("cc3_HC1ET1");
            timer_on("sigmaCC3"); sigmaCC3(i,C_irr,cc3_eval); timer_off("sigmaCC3");
          }
        }
        timer_off("SIGMA ALL");
#else 
        if (params.wfn == "EOM_CC2") cc2_sigma(i,C_irr);  
        else {
          sigmaSS(i,C_irr);
          sigmaSD(i,C_irr);
          sigmaDS(i,C_irr);
          sigmaDD(i,C_irr);
        }
        if ( ((params.wfn == "EOM_CC3") && (cc3_stage>0)) || eom_params.restart_eom_cc3) {
          cc3_HC1(i,C_irr);
          cc3_HC1ET1(i,C_irr);
          sigmaCC3(i,C_irr,cc3_eval);
        }
#endif
        if (params.full_matrix) {
          sigma00(i,C_irr);
          sigma0S(i,C_irr);
          sigma0D(i,C_irr); 
          sigmaS0(i,C_irr);
          sigmaSS_full(i,C_irr);
          sigmaD0(i,C_irr);
          sigmaDS_full(i,C_irr);
          sigmaDD_full(i,C_irr);
        }

#ifdef EOM_DEBUG
        check_sum("reset",0,0);
#endif

        /* Cleaning out sigma vectors for open-shell cases  */
        if (params.eom_ref == 1) {
          sprintf(lbl, "%s %d", "SIA", i);
          dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
          sprintf(lbl, "%s %d", "Sia", i);
          dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
          sprintf(lbl, "%s %d", "SIJAB", i);
          dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
          sprintf(lbl, "%s %d", "Sijab", i);
          dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
          sprintf(lbl, "%s %d", "SIjAb", i);
          dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);

          c_clean(&SIA, &Sia, &SIJAB, &Sijab, &SIjAb);

#ifdef EOM_DEBUG
          check_sum("reset",0,0);
          sprintf(lbl, "Total sigma%d norm af clean", i);
          check_sum(lbl, i, C_irr);
#endif

          dpd_file2_close(&SIA);
          dpd_file2_close(&Sia);
          dpd_buf4_close(&SIJAB);
          dpd_buf4_close(&Sijab);
          dpd_buf4_close(&SIjAb);
          fflush(outfile);
        }
      }

#ifdef TIME_CCEOM
      timer_on("BUILD G");
#endif /*timing*/
      /* Form G = C'*S matrix */
      G = block_matrix(L,L);

      /* reuse values from old G matrix */
      /* if last step was restart, sigma is OK but recompute full G matrix */
      if (ignore_G_old) already_sigma = 0;
      for (i=0; i<already_sigma; ++i)
        for (j=0; j<already_sigma; ++j)
          G[i][j] = G_old[i][j];

      for (i=0;i<L;++i) {

        if(params.eom_ref == 0) {
          /* Spin-adapt C */
          sprintf(lbl, "%s %d", "CME", i);
          dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
					if (params.full_matrix) {
             sprintf(lbl, "%s %d", "C0", i);
             psio_read_entry(EOM_CME, lbl, (char *) &C0, sizeof(double));
					}

          sprintf(lbl, "%s %d", "CMnEf", i);
          dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
          dpd_buf4_copy(&CMnEf, EOM_TMP, "CMnEf");
          dpd_buf4_sort(&CMnEf, EOM_TMP, pqsr, 0, 5, "CMnfE");
          dpd_buf4_close(&CMnEf);

          dpd_buf4_init(&CMnEf1, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMnEf");
          dpd_buf4_scm(&CMnEf1, 2.0);
          dpd_buf4_init(&CMnfE1, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMnfE");
          dpd_buf4_axpy(&CMnfE1, &CMnEf1, -1.0);
          dpd_buf4_close(&CMnfE1);
          dpd_buf4_close(&CMnEf1);

          /* dpd_file2_init(&CME, EOM_TMP, C_irr, 0, 1, "CME");*/
          dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMnEf");
        }
        else if (params.eom_ref == 1) {
          sprintf(lbl, "%s %d", "CME", i);
          dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
          sprintf(lbl, "%s %d", "Cme", i);
          dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
          sprintf(lbl, "%s %d", "CMNEF", i);
          dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
          sprintf(lbl, "%s %d", "Cmnef", i);
          dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
          sprintf(lbl, "%s %d", "CMnEf", i);
          dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
        }
        else if (params.eom_ref == 2) {
          sprintf(lbl, "%s %d", "CME", i);
          dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
          sprintf(lbl, "%s %d", "Cme", i);
          dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, lbl);
          sprintf(lbl, "%s %d", "CMNEF", i);
          dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
          sprintf(lbl, "%s %d", "Cmnef", i);
          dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, lbl);
          sprintf(lbl, "%s %d", "CMnEf", i);
          dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
        }

        /* Dot C's and sigma vectors together to form G matrix */
        for (j=0;j<L;++j) {
				  if (i<already_sigma && j<already_sigma)
					  continue;
           /* fprintf(outfile,"Computing G[%d][%d].\n",i,j); */

          if(params.eom_ref == 0) {
            sprintf(lbl, "%s %d", "SIA", j);
            dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
            tval = 2.0 * dpd_file2_dot(&CME, &SIA);
            sprintf(lbl, "%s %d", "SIjAb", j);
            dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
            tval += dpd_buf4_dot(&CMnEf, &SIjAb);
            dpd_file2_close(&SIA);
						if (params.full_matrix) {
              sprintf(lbl, "%s %d", "S0", j);
							psio_read_entry(EOM_SIA, lbl, (char *) &S0, sizeof(double));
							tval += C0 * S0;
						}
            dpd_buf4_close(&SIjAb);
          }
          else if (params.eom_ref == 1) {
            sprintf(lbl, "%s %d", "SIA", j);
            dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
            tval = dpd_file2_dot(&CME, &SIA);
            dpd_file2_close(&SIA);
            sprintf(lbl, "%s %d", "Sia", j);
            dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
            tval += dpd_file2_dot(&Cme, &Sia);
            dpd_file2_close(&Sia);
            sprintf(lbl, "%s %d", "SIJAB", j);
            dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
            tval += r2aa = dpd_buf4_dot(&CMNEF, &SIJAB);
            dpd_buf4_close(&SIJAB);
            sprintf(lbl, "%s %d", "Sijab", j);
            dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
            tval += r2bb = dpd_buf4_dot(&Cmnef, &Sijab);
            dpd_buf4_close(&Sijab);
            sprintf(lbl, "%s %d", "SIjAb", j);
            dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
            tval += r2ab = dpd_buf4_dot(&CMnEf, &SIjAb);
            dpd_buf4_close(&SIjAb);
            /* fprintf(outfile,"r2aa %12.7lf r2bb %12.7lf r2ab %12.7lf\n", r2aa, r2bb, r2ab); */
          }
          else if (params.eom_ref == 2) {
            sprintf(lbl, "%s %d", "SIA", j);
            dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
            tval = dpd_file2_dot(&CME, &SIA);
            dpd_file2_close(&SIA);
            sprintf(lbl, "%s %d", "Sia", j);
            dpd_file2_init(&Sia, EOM_Sia, C_irr, 2, 3, lbl);
            tval += dpd_file2_dot(&Cme, &Sia);
            dpd_file2_close(&Sia);
            sprintf(lbl, "%s %d", "SIJAB", j);
            dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
            tval += dpd_buf4_dot(&CMNEF, &SIJAB);
            dpd_buf4_close(&SIJAB);
            sprintf(lbl, "%s %d", "Sijab", j);
            dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 12, 17, 12, 17, 0, lbl);
            tval += dpd_buf4_dot(&Cmnef, &Sijab);
            dpd_buf4_close(&Sijab);
            sprintf(lbl, "%s %d", "SIjAb", j);
            dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, lbl);
            tval += dpd_buf4_dot(&CMnEf, &SIjAb);
            dpd_buf4_close(&SIjAb);
          }
          G[i][j] = tval;
        }

        dpd_file2_close(&CME);
        dpd_buf4_close(&CMnEf);
        if (params.eom_ref > 0) {
          dpd_file2_close(&Cme);
          dpd_buf4_close(&CMNEF);
          dpd_buf4_close(&Cmnef);
        }
      } /* end build of G */

      ignore_G_old = 0;
      already_sigma = L;

#ifdef TIME_CCEOM
      timer_off("BUILD G");
#endif /* timing */
#ifdef EOM_DEBUG
      fprintf(outfile,"The G Matrix\n");
      mat_print(G, L, L, outfile);
#endif
      for (i=0;i<L;++i) {
        for (j=0;j<L;++j)
          G_old[i][j] = G[i][j];
      }

      /* Diagonalize G Matrix */ 
      lambda = init_array(L);        /* holds real part of eigenvalues of G */
      alpha = block_matrix(L,L);     /* will hold eigenvectors of G */
      dgeev_eom(L, G, lambda, alpha);
      eigsort(lambda, alpha, L);
      /* eivout(alpha, lambda, L, L, outfile);*/
      free_block(G);

      /* Open up residual vector files */
      if (params.eom_ref == 0) {
        dpd_file2_init(&RIA, EOM_R, C_irr, 0, 1, "RIA");
        dpd_buf4_init(&RIjAb, EOM_R, C_irr, 0, 5, 0, 5, 0, "RIjAb");
      }
      else if (params.eom_ref == 1) {
        dpd_file2_init(&RIA, EOM_R, C_irr, 0, 1, "RIA");
        dpd_file2_init(&Ria, EOM_R, C_irr, 0, 1, "Ria");
        dpd_buf4_init(&RIJAB, EOM_R, C_irr, 2, 7, 2, 7, 0, "RIJAB");
        dpd_buf4_init(&Rijab, EOM_R, C_irr, 2, 7, 2, 7, 0, "Rijab");
        dpd_buf4_init(&RIjAb, EOM_R, C_irr, 0, 5, 0, 5, 0, "RIjAb");
      }
      else if (params.eom_ref == 2) {
        dpd_file2_init(&RIA, EOM_R, C_irr, 0, 1, "RIA");
        dpd_file2_init(&Ria, EOM_R, C_irr, 2, 3, "Ria");
        dpd_buf4_init(&RIJAB, EOM_R, C_irr, 2, 7, 2, 7, 0, "RIJAB");
        dpd_buf4_init(&Rijab, EOM_R, C_irr, 12, 17, 12, 17, 0, "Rijab");
        dpd_buf4_init(&RIjAb, EOM_R, C_irr, 22, 28, 22, 28, 0, "RIjAb");
      }
      fprintf(outfile,"  Root    EOM Energy     Delta E   Res. Norm    Conv?\n");
      for (k=0;k<eom_params.cs_per_irrep[C_irr];++k) {
#ifdef TIME_CCEOM
      timer_on("CALC RES");
#endif /* timing */

      /* rezero residual vector for each root */
      if (params.full_matrix) {
        R0 = 0.0;
        psio_write_entry(EOM_R, "R0", (char *) &R0, sizeof(double));
      }
      dpd_file2_scm(&RIA, 0.0);
      dpd_buf4_scm(&RIjAb, 0.0);
      if (params.eom_ref > 0) {
        dpd_file2_scm(&Ria, 0.0);
        dpd_buf4_scm(&RIJAB, 0.0);
        dpd_buf4_scm(&Rijab, 0.0);
      }

        /* only one cc3 root can be sought */
        if ( (params.wfn == "EOM_CC3") && (cc3_stage>0) )  {
          if (eom_params.follow_root)
            k = follow_root(L, alpha, C_irr);
          else
            k = eom_params.prop_root;
        }

        converged[k] = 0;
        for (i=0;i<L;++i) {
          if (params.eom_ref == 0) { /* RHF residual */
            sprintf(lbl, "%s %d", "SIA", i);
            dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
            sprintf(lbl, "%s %d", "CME", i);
            dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
            dpd_file2_axpbycz(&CME, &SIA, &RIA, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_file2_close(&CME);
            dpd_file2_close(&SIA);

            sprintf(lbl, "%s %d", "CMnEf", i);
            dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
            sprintf(lbl, "%s %d", "SIjAb", i);
            dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
            dpd_buf4_axpbycz(&CMnEf, &SIjAb, &RIjAb, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_buf4_close(&CMnEf);
            dpd_buf4_close(&SIjAb);

            if (params.full_matrix) {
              sprintf(lbl, "%s %d", "S0", i);
              psio_read_entry(EOM_SIA, lbl, (char *) &(S0), sizeof(double));
              sprintf(lbl, "%s %d", "C0", i);
              psio_read_entry(EOM_CME, lbl, (char *) &(C0), sizeof(double));
              psio_read_entry(EOM_R, "R0", (char *) &(R0), sizeof(double));
              R0 += -1.0*lambda[k]*alpha[i][k]*C0 + alpha[i][k]*S0;
              psio_write_entry(EOM_R, "R0", (char *) &(R0), sizeof(double));
            }
          }
          else if (params.eom_ref == 1) { /* ROHF residual */
            sprintf(lbl, "%s %d", "SIA", i);
            dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
            sprintf(lbl, "%s %d", "CME", i);
            dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
            dpd_file2_axpbycz(&CME, &SIA, &RIA, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_file2_close(&CME);
            dpd_file2_close(&SIA);
  
            sprintf(lbl, "%s %d", "CMnEf", i);
            dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
            sprintf(lbl, "%s %d", "SIjAb", i);
            dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
            dpd_buf4_axpbycz(&CMnEf, &SIjAb, &RIjAb, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_buf4_close(&CMnEf);
            dpd_buf4_close(&SIjAb);

            sprintf(lbl, "%s %d", "Sia", i);
            dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
            sprintf(lbl, "%s %d", "Cme", i);
            dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
            dpd_file2_axpbycz(&Cme, &Sia, &Ria, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_file2_close(&Cme);
            dpd_file2_close(&Sia);

            sprintf(lbl, "%s %d", "CMNEF", i);
            dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
            sprintf(lbl, "%s %d", "SIJAB", i);
            dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
            dpd_buf4_axpbycz(&CMNEF, &SIJAB, &RIJAB, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_buf4_close(&CMNEF);
            dpd_buf4_close(&SIJAB);

            sprintf(lbl, "%s %d", "Cmnef", i);
            dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
            sprintf(lbl, "%s %d", "Sijab", i);
            dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
            dpd_buf4_axpbycz(&Cmnef, &Sijab, &Rijab, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_buf4_close(&Cmnef);
            dpd_buf4_close(&Sijab);
          }
          else if (params.eom_ref == 2) { /* UHF residual */
            sprintf(lbl, "%s %d", "SIA", i);
            dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
            sprintf(lbl, "%s %d", "CME", i);
            dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
            dpd_file2_axpbycz(&CME, &SIA, &RIA, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_file2_close(&CME);
            dpd_file2_close(&SIA);
  
            sprintf(lbl, "%s %d", "CMnEf", i);
            dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
            sprintf(lbl, "%s %d", "SIjAb", i);
            dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, lbl);
            dpd_buf4_axpbycz(&CMnEf, &SIjAb, &RIjAb, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_buf4_close(&CMnEf);
            dpd_buf4_close(&SIjAb);

            sprintf(lbl, "%s %d", "Sia", i);
            dpd_file2_init(&Sia, EOM_Sia, C_irr, 2, 3, lbl);
            sprintf(lbl, "%s %d", "Cme", i);
            dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, lbl);
            dpd_file2_axpbycz(&Cme, &Sia, &Ria, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_file2_close(&Cme);
            dpd_file2_close(&Sia);

            sprintf(lbl, "%s %d", "CMNEF", i);
            dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
            sprintf(lbl, "%s %d", "SIJAB", i);
            dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
            dpd_buf4_axpbycz(&CMNEF, &SIJAB, &RIJAB, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_buf4_close(&CMNEF);
            dpd_buf4_close(&SIJAB);

            sprintf(lbl, "%s %d", "Cmnef", i);
            dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, lbl);
            sprintf(lbl, "%s %d", "Sijab", i);
            dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 12, 17, 12, 17, 0, lbl);
            dpd_buf4_axpbycz(&Cmnef, &Sijab, &Rijab, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
            dpd_buf4_close(&Cmnef);
            dpd_buf4_close(&Sijab);
          }
        }

#ifdef EOM_DEBUG
        if (params.eom_ref == 0) {
          dpd_buf4_sort(&RIjAb, EOM_TMP, pqsr, 0, 5, "RIjbA");
          dpd_buf4_init(&RIjbA, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "RIjbA");
          norm = norm_C_rhf(&RIA, &RIjAb, &RIjbA);
          dpd_buf4_close(&RIjbA);
        }
        else norm = norm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb);
        fprintf(outfile,"Norm of residual vector %d  before precondition %18.13lf\n",k,norm);
#endif
#ifdef TIME_CCEOM
      timer_off("CALC RES");
#endif /* timing */
    /* moved back down 8-06 why was this ever before residual norm? 
	if(params.eom_ref == 0) precondition_RHF(&RIA, &RIjAb, lambda[k]);
	else precondition(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, lambda[k]); */

        if (params.eom_ref == 0) {
          dpd_buf4_sort(&RIjAb, EOM_TMP, pqsr, 0, 5, "RIjbA");
          dpd_buf4_init(&RIjbA, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "RIjbA");
          if (!params.full_matrix) norm = norm_C_rhf(&RIA, &RIjAb, &RIjbA);
          else {
            psio_read_entry(EOM_R, "R0", (char *) &R0, sizeof(double));
            norm = norm_C_rhf_full(R0, &RIA, &RIjAb, &RIjbA);
          }
        }
        else norm = norm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb);

#ifdef EOM_DEBUG
        fprintf(outfile,"Norm of residual vector %d  after precondition %18.13lf\n",k,norm);
#endif

        fprintf(outfile,"%22d%15.10lf%11.2e%12.2e",k+1,lambda[k], lambda[k]-lambda_old[k], norm);

        /* Check for convergence and add new vector if not converged */
        if ( (norm > eom_params.residual_tol) || (fabs(lambda[k]-lambda_old[k]) > eom_params.eval_tol) ) {
          fprintf(outfile,"%7s\n","N");

          if(params.eom_ref == 0) precondition_RHF(&RIA, &RIjAb, lambda[k]);
          else precondition(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, lambda[k]);

	  if(params.eom_ref == 0) {

            /* Normalize R */
            dpd_buf4_sort(&RIjAb, EOM_TMP, pqsr, 0, 5, "RIjbA");
            dpd_buf4_init(&RIjbA, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "RIjbA");

            if (!params.full_matrix) norm = norm_C_rhf(&RIA, &RIjAb, &RIjbA);
            else {
              psio_read_entry(EOM_R, "R0", (char *) &R0, sizeof(double));
              norm = norm_C_rhf_full(R0, &RIA, &RIjAb, &RIjbA);
            }
            dpd_buf4_close(&RIjbA);

            if (params.full_matrix) {
              R0 *= 1.0/norm;
              psio_write_entry(EOM_R, "R0", (char *) &R0, sizeof(double));
            }
            dpd_file2_scm(&RIA, 1.0/norm);
            dpd_buf4_scm(&RIjAb, 1.0/norm);
          }
          else {
            norm = norm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb);
            scm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, 1.0/norm);
          }

#ifdef EOM_DEBUG
          fprintf(outfile,"Norm of residual vector af preconditioning %18.13lf\n",norm);
#endif

          if(params.eom_ref == 0) schmidt_add_RHF(&RIA, &RIjAb, &numCs, C_irr);
          else schmidt_add(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, &numCs, C_irr);
        }
        else {
          fprintf(outfile,"%7s\n","Y");
          ++num_converged;
          converged[k] = 1;
        }

        /* only one cc3 root can be sought */
        if ( (params.wfn == "EOM_CC3") && (cc3_stage>0) ) {
          cc3_index = k;
          cc3_eval = lambda[k];
          /* fprintf(outfile,"Setting CC3 eigenvalue to %15.10lf\n",cc3_eval); */
        }
      }

      dpd_file2_close(&RIA);
      dpd_buf4_close(&RIjAb);
      if (params.eom_ref > 0) {
        dpd_file2_close(&Ria);
        dpd_buf4_close(&RIJAB);
        dpd_buf4_close(&Rijab);
      }

      for (i=0;i<eom_params.cs_per_irrep[C_irr];++i) lambda_old[i] = lambda[i];
      free(lambda);
      if ( (params.wfn == "EOM_CC3") && (cc3_stage>0) ) {
          lambda_old[cc3_index] = cc3_eval; /* a hack to make Delta E work next iteration */
      }

      /* restart with new B vectors if there are too many */
      if (L >= vectors_per_root * eom_params.cs_per_irrep[C_irr]) {
        if ( (params.wfn == "EOM_CC3") && (cc3_stage>0) ) {
          restart(alpha, L, cc3_index+1, C_irr, 1, alpha_old, L_old, eom_params.collapse_with_last);
          L_old = L;
          L = cc3_index+1;
          if (eom_params.collapse_with_last) L *= 2;
          fprintf(outfile,"Collapsing to %d vector(s).\n",L);
          already_sigma = 0;
          ignore_G_old = 1;
        }
        else {
          restart(alpha, L, eom_params.cs_per_irrep[C_irr], C_irr, 1,alpha_old,L_old,eom_params.collapse_with_last);
          L_old = L;
          if (eom_params.collapse_with_last)
            L = 2 * eom_params.cs_per_irrep[C_irr];
          else
            L = eom_params.cs_per_irrep[C_irr];
          already_sigma = L;
          ignore_G_old = 1;
        }
        keep_going = 1;
        /* keep track of number of triples restarts */
        if ( (params.wfn == "EOM_CC3") && (cc3_stage>0) ) {
          fprintf(outfile,"Change in CC3 energy from last iterated value %15.10lf\n",
              cc3_eval - cc3_last_converged_eval);
          cc3_last_converged_eval = cc3_eval;
          ++num_cc3_restarts;
          cc3_stage = 2;
        }
      }
      else {
        /* If any new vectors were added, then continue */
        if (numCs > L) {
          keep_going = 1;
          L_old = L;
          L = numCs;
        }
      }

      ++iter;
      if (eom_params.restart_eom_cc3 && (cc3_stage==0)) {
        fprintf(outfile,"Jumping to EOM CC3 iterations\n");
        keep_going = 0;
      }
      if ((keep_going == 0) && (iter < eom_params.max_iter)) {

        /* for CC3: done with EOM CCSD - now do EOM CC3 */
        if ( (params.wfn == "EOM_CC3") && (cc3_stage == 0) ) {
          fprintf(outfile, "Completed EOM_CCSD\n");
          fprintf(outfile,"Collapsing to only %d vector(s).\n", eom_params.cs_per_irrep[C_irr]);
          if (!eom_params.restart_eom_cc3) {
            restart(alpha, L, eom_params.cs_per_irrep[C_irr], C_irr, 1,alpha_old,L_old,0);
            save_C_ccsd(eom_params.prop_root, C_irr);
          }

          cc3_last_converged_eval = cc3_eval = lambda_old[eom_params.prop_root];
          fprintf(outfile,"Setting initial CC3 eigenvalue to %15.10lf\n",cc3_eval);

          L_old = L;
          L = eom_params.cs_per_irrep[C_irr];
          eom_params.cs_per_irrep[C_irr] = 1; /* only get 1 CC3 solution */
          keep_going = 1;
          already_sigma = 0;
          ignore_G_old = 1; 
          iter = 0;
          cc3_stage = 1;
          vectors_per_root = eom_params.vectors_cc3;
        }
        else if( (params.wfn == "EOM_CC3") && /* can't trust sigmas yet */
               ( (cc3_stage == 1) || fabs(cc3_eval-cc3_last_converged_eval)>eom_params.eval_tol)) {
          /* for CC3: restart one time if no cc3_restarts have yet been done */
          if (cc3_stage == 1) fprintf(outfile, "Forcing one restart with sigma recomputation.\n");
          else fprintf(outfile,"Forcing restart to make sure new sigma vectors give same eigenvalue.\n");
          fprintf(outfile,"Collapsing to only %d vector(s).\n", cc3_index+1);
          restart(alpha, L, cc3_index+1, C_irr, 1, alpha_old, L_old, 0);
          cc3_eval = lambda_old[cc3_index];
          if (cc3_stage == 1)
            fprintf(outfile,"Change in CC3 energy from last iterated value %15.10lf\n", cc3_eval - 0.0);
          else
            fprintf(outfile,"Change in CC3 energy from last iterated value %15.10lf\n",
              cc3_eval-cc3_last_converged_eval);
          cc3_last_converged_eval = cc3_eval;
          fprintf(outfile,"Setting old CC3 eigenvalue to %15.10lf\n",cc3_eval);
          keep_going = 1;
          already_sigma = 0;
          ignore_G_old = 1;
          L_old = L;
          L = cc3_index+1;
          cc3_stage = 2;
        }
        else if (params.wfn == "EOM_CC3") {
          /* for CC3: collapse to one final root and place in location 0 */
          fprintf(outfile,"Collapsing to only %d vector(s).\n", cc3_index+1);
          restart(alpha, L, cc3_index+1, C_irr, 0, alpha_old, L_old, 0);
          if (cc3_index > 0) restart_with_root(cc3_index, C_irr);
          converged[0] = 1;
          cc3_eval = lambda_old[0] = lambda_old[cc3_index];
          fprintf(outfile,"Change in CC3 energy from last iterated value %15.10lf\n",
              cc3_eval - cc3_last_converged_eval);
        }
        else {
          fprintf(outfile,"Collapsing to only %d vector(s).\n", eom_params.cs_per_irrep[C_irr]);
          restart(alpha, L, eom_params.cs_per_irrep[C_irr], C_irr, 0, alpha_old, L_old, 0);
        }
      }
      alpha_old = block_matrix(L_start_iter,L_start_iter);
      for (k=0;k<L_start_iter;++k)
        for (i=0;i<L_start_iter;++i)
          alpha_old[i][k] = alpha[i][k];
      free_block(alpha);
    }
    free_block(G_old);

    fprintf(outfile,"\nProcedure converged for %d root(s).\n",num_converged);
    if (num_converged == eom_params.cs_per_irrep[C_irr]) { }
    else if (iter == eom_params.max_iter) {
      fprintf(outfile,"\nMaximum number of iterations exceeded, ");
      fprintf(outfile,"so not all roots converged!\n\n");
    }
    else if ( (params.wfn == "EOM_CC3") && (num_converged == 1) ) { }
    else {
      fprintf(outfile,"\nAlgorithm failure: No vectors could be added, ");
      fprintf(outfile,"though not all roots converged!\n\n");
    }

    /* write Cs and energies to RAMPS file */
    write_Rs(C_irr, lambda_old, converged);
    /* compute R0 and normalize - also do any orthogonality checks */
    if (params.eom_ref == 0)
      rzero_rhf(C_irr, converged);
    else
      rzero(C_irr, converged);

    if (num_converged > 0) {
      fprintf(outfile,"\nFinal Energetic Summary for Converged Roots of Irrep %s\n",
	      moinfo.irr_labs[moinfo.sym^C_irr]);
      fprintf(outfile,"                     Excitation Energy              Total Energy\n");
      fprintf(outfile,"                (eV)     (cm^-1)     (au)             (au)\n");
      for (i=0;i<eom_params.cs_per_irrep[C_irr];++i) {
        if (converged[i] == 1) {
          if (!params.full_matrix) totalE =lambda_old[i]+moinfo.eref+moinfo.ecc; 
          else totalE =lambda_old[i]+moinfo.eref;

          // save a list of all converged energies in order by irrep and then energy
          eom_params.state_energies[num_converged_index] = totalE;
          // Put this list in environment for testing - I don't like it much because it mixes all
          // the irreps together
          /*- strings so that variable-name psi variables get parsed in docs -*/
          /*- Process::environment.globals["CC ROOT $N$ TOTAL ENERGY"] -*/

          std::stringstream s;
          s << "CC ROOT " << (num_converged_index+1) << " TOTAL ENERGY";
          Process::environment.globals[s.str()] = totalE;
					
          fprintf(outfile,"EOM State %d %10.3lf %10.1lf %14.10lf  %17.12lf\n", ++num_converged_index,
             lambda_old[i]* _hartree2ev, lambda_old[i]* _hartree2wavenumbers, lambda_old[i], totalE);

          /* print out largest components of wavefunction */
          fprintf(outfile, "\nLargest components of excited wave function #%d:\n", num_converged_index);
          if(params.eom_ref == 0) {
            sprintf(lbl, "%s %d %d", "RIA", C_irr, i);
            dpd_file2_init(&CME, CC_RAMPS, C_irr, 0, 1, lbl);
            sprintf(lbl, "%s %d %d", "RIjAb", C_irr, i);
            dpd_buf4_init(&CMnEf, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);

            amp_write_RHF(&CME, &CMnEf, eom_params.amps_to_print);

            dpd_file2_close(&CME);
            dpd_buf4_close(&CMnEf);
          }
          else if (params.eom_ref == 1) {
            sprintf(lbl, "%s %d %d", "RIA", C_irr, i);
            dpd_file2_init(&CME, CC_RAMPS, C_irr, 0, 1, lbl);
            sprintf(lbl, "%s %d %d", "Ria", C_irr, i);
            dpd_file2_init(&Cme, CC_RAMPS, C_irr, 0, 1, lbl);
            sprintf(lbl, "%s %d %d", "RIJAB", C_irr, i);
            dpd_buf4_init(&CMNEF, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
            sprintf(lbl, "%s %d %d", "Rijab", C_irr, i);
            dpd_buf4_init(&Cmnef, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
            sprintf(lbl, "%s %d %d", "RIjAb", C_irr, i);
            dpd_buf4_init(&CMnEf, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);

            amp_write_ROHF(&CME, &Cme, &CMNEF, &Cmnef, &CMnEf, eom_params.amps_to_print);

            dpd_file2_close(&CME);
            dpd_file2_close(&Cme);
            dpd_buf4_close(&CMNEF);
            dpd_buf4_close(&Cmnef);
            dpd_buf4_close(&CMnEf);
          }
          else if (params.eom_ref == 2) {
            sprintf(lbl, "%s %d %d", "RIA", C_irr, i);
            dpd_file2_init(&CME, CC_RAMPS, C_irr, 0, 1, lbl);
            sprintf(lbl, "%s %d %d", "Ria", C_irr, i);
            dpd_file2_init(&Cme, CC_RAMPS, C_irr, 2, 3, lbl);
            sprintf(lbl, "%s %d %d", "RIJAB", C_irr, i);
            dpd_buf4_init(&CMNEF, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
            sprintf(lbl, "%s %d %d", "Rijab", C_irr, i);
            dpd_buf4_init(&Cmnef, CC_RAMPS, C_irr, 12, 17, 12, 17, 0, lbl);
            sprintf(lbl, "%s %d %d", "RIjAb", C_irr, i);
            dpd_buf4_init(&CMnEf, CC_RAMPS, C_irr, 22, 28, 22, 28, 0, lbl);

            amp_write_UHF(&CME, &Cme, &CMNEF, &Cmnef, &CMnEf, eom_params.amps_to_print);

            dpd_file2_close(&CME);
            dpd_file2_close(&Cme);
            dpd_buf4_close(&CMNEF);
            dpd_buf4_close(&Cmnef);
            dpd_buf4_close(&CMnEf);
          } // UHF

          // The 'key' or 'property' root is stored in eom_params.prop_sym and prop_root
          // by default it is the uppermost state but not necesarily
          if (C_irr == eom_params.prop_sym && i == eom_params.prop_root) {
            fprintf(outfile,"\n\tPutting into environment energy for root of R irrep %d and root %d.\n", C_irr+1, i+1);
            Process::environment.globals["CURRENT ENERGY"] = totalE;
            fprintf(outfile,"\tPutting into environment CURRENT ENERGY:             %15.10lf\n", totalE);
//            Process::environment.globals["CURRENT CORRELATION ENERGY"] = lambda_old[i]+moinfo.ecc;
//            fprintf(outfile,"\tPutting into environment CURRENT CORRELATION ENERGY: %15.10lf\n", lambda_old[i]+moinfo.ecc);
          }

        } // converged[i] == 1
      } // i
    } // if num_converged > 0
    fprintf(outfile,"\n");

    free(lambda_old);
    free_block(alpha_old);
    free(converged);
  }

  /* Save the energy vector to checkpoint */
  chkpt_init(PSIO_OPEN_OLD);
  keyw = chkpt_build_keyword("EOM Number of States");
  psio_write_entry(PSIF_CHKPT, keyw, (char *) &(eom_params.number_of_states), sizeof(int));
  free(keyw);
  keyw = chkpt_build_keyword("EOM State Energies");
  psio_write_entry(PSIF_CHKPT, keyw, (char *) eom_params.state_energies, eom_params.number_of_states * sizeof(double));
  free(keyw);
  chkpt_close();

  fprintf(outfile,"\tTotal # of sigma evaluations: %d\n",nsigma_evaluations);
  return;
}

void init_C0(int i) {
  char lbl[32];
  double zip = 0.0;
  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "C0", i);
    psio_write_entry(EOM_CME, lbl, (char *) &(zip), sizeof(double));
  }
}
void init_S0(int i) {
  char lbl[32];
  double zip = 0.0;
  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "S0", i);
    psio_write_entry(EOM_SIA, lbl, (char *) &(zip), sizeof(double));
  }
}

/* zeroes ith CME (and Cme) vectors on disk */
void init_C1(int i, int C_irr ){
  dpdfile2 CME, Cme;
  char lbl[32];
  double zip = 0.0;
  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "CME", i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_scm(&CME, 0.0);
    dpd_file2_close(&CME);
    if (params.full_matrix) {
      sprintf(lbl, "%s %d", "C0", i);
      psio_write_entry(EOM_CME, lbl, (char *) &(zip), sizeof(double));
    }
  }
  else {
    sprintf(lbl, "%s %d", "CME", i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Cme", i);
    if (params.eom_ref == 1) dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
    else if (params.eom_ref == 2) dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, lbl);
    scm_C1(&CME, &Cme, 0.0);
    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }
}

/* zeroes ith SIA (and Sia) vectors on disk */
void init_S1(int i, int C_irr) {
  dpdfile2 SIA, Sia;
  char lbl[32];
  double zip = 0.0;
  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    dpd_file2_scm(&SIA, 0.0);
    dpd_file2_close(&SIA);
    if (params.full_matrix) {
      sprintf(lbl, "%s %d", "S0", i);
      psio_write_entry(EOM_SIA, lbl, (char *) &(zip), sizeof(double));
    }
  }
  else {
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", i);
    if (params.eom_ref == 1) dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
    else if (params.eom_ref == 2) dpd_file2_init(&Sia, EOM_Sia, C_irr, 2, 3, lbl);
    scm_C1(&SIA, &Sia, 0.0);
    dpd_file2_close(&SIA);
    dpd_file2_close(&Sia);
  }
}

/* zeroes ith CMnEf (+ CMNEF + Cmnef) on disk */
void init_C2(int i, int C_irr) {
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  char lbl[32];
  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_scm(&CMnEf, 0.0);
    dpd_buf4_close(&CMnEf);
  }
  else {
    sprintf(lbl, "%s %d", "CMNEF", i);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);

    sprintf(lbl, "%s %d", "Cmnef", i);
    if (params.eom_ref == 1)
      dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
    else if (params.eom_ref == 2)
      dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, lbl);

    sprintf(lbl, "%s %d", "CMnEf", i);
    if (params.eom_ref == 1)
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    else if (params.eom_ref == 2)
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);

    /* scm_C2(&CMNEF, &Cmnef, &CMnEf, 0.0); */
    dpd_buf4_scm(&CMNEF, 0.0);
    dpd_buf4_scm(&Cmnef, 0.0);
    dpd_buf4_scm(&CMnEf, 0.0);

    dpd_buf4_close(&CMNEF);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_close(&CMnEf);
  }
}


void init_S2(int i, int C_irr) {
  dpdbuf4 SIJAB, Sijab, SIjAb;
  char lbl[32];
  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "SIjAb", i);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_scm(&SIjAb, 0.0);
    dpd_buf4_close(&SIjAb);
  }
  else {
    sprintf(lbl, "%s %d", "SIJAB", i);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "Sijab", i);
    if (params.eom_ref == 1) 
      dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
    else if (params.eom_ref == 2) 
      dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 12, 17, 12, 17, 0, lbl);

    sprintf(lbl, "%s %d", "SIjAb", i);
    if (params.eom_ref == 1)
      dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
    else if (params.eom_ref == 2)
      dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, lbl);
    scm_C2(&SIJAB, &Sijab, &SIjAb, 0.0);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_close(&Sijab);
    dpd_buf4_close(&SIjAb);
  }
}

}} // namespace psi::cceom
