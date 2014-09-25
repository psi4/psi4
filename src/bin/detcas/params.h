/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
/*
** PARAMS.H
** 
** C. David Sherrill
** University of California, Berkeley
** 1998
*/

#ifndef _psi_src_bin_detcas_params_h
#define _psi_src_bin_detcas_params_h

namespace psi { namespace detcas {

/*
** parameters structure: holds user-specified parameters
*/
struct params {
   char *dertype;           /* derivative level: none, first, etc.          */
   char *wfn;               /* wavefunction, CASSCF, RASSCF, ..             */
   int print_lvl;           /* print verbosity level                        */ 
   int print_mos;           /* print the molecular orbitals ?               */
   int rms_grad_convergence;/* convergence, 10^-n, on RMS of orbital grad   */
   int energy_convergence;  /* convergence, 10^-n, on CI energy             */
   int oei_file;            /* file number for one-electron integrals       */
   int oei_erase;           /* erase onel ints after reading them?          */
   int tei_file;            /* file number for two-electron integrals       */
   int tei_erase;           /* erase twoel ints after reading them?         */
   int opdm_file;           /* file number for one-particle density matrix  */
   int opdm_erase;          /* erase onepdm ints after reading?             */
   int tpdm_file;           /* file number for two-particle density matrix  */
   int tpdm_erase;          /* erase twopdm after reading?                  */
   int lag_file;            /* file number for lagrangian                   */
   int lag_erase;           /* erase lagrangian after reading?              */
   int ignore_ras_ras;      /* ignore RAS/RAS rotations in independ pairs?  */
   int ignore_fz;           /* ignore FZC/FZV in independent pair list?     */
   int filter_ints;         /* filter out the frozen orbital integrals?     */
   int scale_grad;          /* scale the orbital gradient by the appx Hess? */
   int diis_start;          /* how many diis vectors built up before start  */
   int diis_freq;           /* how many iters to go before a diis step      */
   int diis_min_vecs;       /* how many vectors required before do diis?    */
   int diis_max_vecs;       /* how many vectors maximum to hold?            */
   double scale_step;       /* stepsize scaling factor                      */
   char *hessian;           /* string describing type of MO Hessian         */
                            /* DIAG, APPROX_DIAG, or FULL                   */
   int use_fzc_h;           /* Use frozen-core operator h?(1) Or bare h?(0) */
                            /* this determines which onel ints are read     */ 
   int level_shift;         /* Allow for level shifting of the hessian?     */
   double shift;            /* How much do I level shift the hessian        */
   double determ_min;       /* Min det of MO Hessian before levelshift      */
   double step_max;         /* Biggest single allowed theta step            */
   int invert_hessian;      /* If=1, directly invert the Hessian, 
                               if=0, solve linear equations H delta = -g    */
   int use_thetas;          /* If=1, use Givens matrix formalism,
                               if=0, use YY 2nd-order expansion U=e^R       */
   int force_step;          /* Ignore usual updating and force a user
                               specified step?  (For debugging)             */
   int force_pair;          /* If force_step=1, which indep pair to step?   */
   double force_value;      /* If force_step=1, how far to step?            */
   double scale_act_act;    /* Scale the active-active Hessian by this      */
   int bfgs;                /* Do BFGS update of Hessian?                   */
   int ds_hessian;          /* Do a DS Hessian update?                      */
  };

}} // end namespace psi::detcas

#endif // header guard
