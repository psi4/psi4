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

#include <string>

namespace psi { namespace detcas {

/*
** parameters structure: holds user-specified parameters
*/
struct params {
   std::string *dertype;        /* derivative level: none, first, etc.          */
   std::string *wfn;            /* wavefunction, CASSCF, RASSCF, ..             */
   int print_lvl;               /* print verbosity level                        */ 
   bool print_mos;              /* print the molecular orbitals ?               */
   double rms_grad_convergence; /* convergence on RMS of orbital grad           */
   double energy_convergence;   /* convergence on CI energy                     */
   int oei_file;                /* file number for one-electron integrals       */
   bool oei_erase;              /* erase onel ints after reading them?          */
   int tei_file;                /* file number for two-electron integrals       */
   bool tei_erase;              /* erase twoel ints after reading them?         */
   int opdm_file;               /* file number for one-particle density matrix  */
   bool opdm_erase;             /* erase onepdm ints after reading?             */
   int tpdm_file;               /* file number for two-particle density matrix  */
   bool tpdm_erase;             /* erase twopdm after reading?                  */
   int lag_file;                /* file number for lagrangian                   */
   bool lag_erase;              /* erase lagrangian after reading?              */
   bool ignore_ras_ras;         /* ignore RAS/RAS rotations in independ pairs?  */
   bool ignore_fz;              /* ignore FZC/FZV in independent pair list?     */
   int filter_ints;             /* filter out the frozen orbital integrals?     */
   bool scale_grad;             /* scale the orbital gradient by the appx Hess? */
   int diis_start;              /* how many diis vectors built up before start  */
   int diis_freq;               /* how many iters to go before a diis step      */
   int diis_min_vecs;           /* how many vectors required before do diis?    */
   int diis_max_vecs;           /* how many vectors maximum to hold?            */
   double scale_step;           /* stepsize scaling factor                      */
   std::string *hessian;        /* string describing type of MO Hessian         */
                                /* DIAG, APPROX_DIAG, or FULL                   */
   bool use_fzc_h;              /* Use frozen-core operator h?(1) Or bare h?(0) */
                                /* this determines which onel ints are read     */ 
   bool level_shift;            /* Allow for level shifting of the hessian?     */
   double shift;                /* How much do I level shift the hessian        */
   double determ_min;           /* Min det of MO Hessian before levelshift      */
   double step_max;             /* Biggest single allowed theta step            */
   bool invert_hessian;         /* If=True, directly invert the Hessian, 
                                   if=False, solve linear equations H delta = -g*/
   bool use_thetas;             /* If=True, use Givens matrix formalism,
                                   if=False, use YY 2nd-order expansion U=e^R   */
   bool force_step;             /* Ignore usual updating and force a user
                                   specified step?  (For debugging)             */
   int force_pair;              /* If force_step=True, which indep pair to step?*/
   double force_value;          /* If force_step=True, how far to step?         */
   double scale_act_act;        /* Scale the active-active Hessian by this      */
   bool bfgs;                   /* Do BFGS update of Hessian?                   */
   bool ds_hessian;             /* Do a DS Hessian update?                      */
  };

}} // end namespace psi::detcas

#endif // header guard
