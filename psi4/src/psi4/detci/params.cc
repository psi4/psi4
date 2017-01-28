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
    \ingroup DETCI
    \brief Enter brief description of file here
*/
/*
** PARAMS.CC: File contains functions which get or print the running
**    parameters for the CI calculation.
**
** David Sherrill, 16 November 1994
**
*/

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/psifiles.h"
#include "psi4/detci/ciwave.h"
#include "psi4/detci/structs.h"
#include "psi4/psi4-dec.h"

namespace psi { namespace detci {

/*
** get_parameters(): Function gets the program running parameters such
**   as convergence.  These are stored in the Parameters data structure.
*/
void CIWavefunction::get_parameters(Options &options)
{
  int i, j, k, errcod;
  int iopen=0, tval;
  char line1[133];
  double junk;

  /* need to figure out wheter to filter tei's */
  Parameters_->dertype = options.get_str("DERTYPE");
  Parameters_->wfn = options.get_str("WFN");


  // CDS-TODO: Check these
  print_ = 1;

  Parameters_->ex_lvl = options.get_int("EX_LEVEL");
  Parameters_->cc_ex_lvl = options.get_int("CC_EX_LEVEL");
  Parameters_->val_ex_lvl = options.get_int("VAL_EX_LEVEL");
  Parameters_->cc_val_ex_lvl = options.get_int("CC_VAL_EX_LEVEL");

  Parameters_->cc_a_val_ex_lvl = -1;
  Parameters_->cc_b_val_ex_lvl = -1;

  Parameters_->num_roots = options.get_int("NUM_ROOTS");

  Parameters_->istop = options["ISTOP"].to_integer();
  Parameters_->print_ciblks = options["CIBLKS_PRINT"].to_integer();

  if (options["PRINT"].has_changed()) {
    print_ = options.get_int("PRINT");
  }

  Parameters_->opentype = PARM_OPENTYPE_UNKNOWN;

  Parameters_->ref_sym = options.get_int("REFERENCE_SYM");

  Parameters_->hd_ave = EVANGELISTI;

  if (Parameters_->wfn == "CASSCF") {
    Parameters_->fci = 1;
    Parameters_->mcscf = 1;
  }
  else if (Parameters_->wfn == "RASSCF"){
    Parameters_->fci = 0;
    Parameters_->mcscf = 1;
  }
  else {
    Parameters_->fci = 0;
    Parameters_->mcscf = 0;
  }

  if (Parameters_->wfn == "ZAPTN") {
    Parameters_->mpn = 1;
    Parameters_->zaptn = 1;
  }
  else {
    Parameters_->mpn = 0;
    Parameters_->zaptn = 0;
  }

  Parameters_->z_scale_H = 0;
  Parameters_->ras1_lvl = -1;
  Parameters_->ras1_min = -1;
  Parameters_->a_ras1_lvl = -1;
  Parameters_->a_ras1_min = -1;
  Parameters_->b_ras1_lvl = -1;
  Parameters_->b_ras1_min = -1;
  Parameters_->ras3_lvl = -1;
  Parameters_->ras4_lvl = -1;

  Parameters_->guess_vector = PARM_GUESS_VEC_H0_BLOCK;

  Parameters_->neg_only = 1;
  Parameters_->nunits = 1;
  Parameters_->hd_filenum = options["CI_FILE_START"].to_integer();
  Parameters_->c_filenum = options["CI_FILE_START"].to_integer() + 1;
  Parameters_->s_filenum = options["CI_FILE_START"].to_integer() + 2;
  Parameters_->d_filenum = options["CI_FILE_START"].to_integer() + 3;

  Parameters_->cc = options["CC"].to_integer();

  if (Parameters_->dertype == "FIRST" ||
      Parameters_->wfn == "DETCAS" ||
      Parameters_->wfn == "CASSCF" ||
      Parameters_->wfn == "RASSCF")
  {
     Parameters_->convergence = 1e-7;
     Parameters_->energy_convergence = 1e-8;
     Parameters_->opdm = 1;
     Parameters_->tpdm = 1;
     Parameters_->maxiter = 12;
     Parameters_->die_if_not_converged = false;
  }

  else {
    if (Parameters_->cc) {
      Parameters_->convergence = 1e-5;
      Parameters_->energy_convergence = 1e-7;
      Parameters_->maxiter = 24;
      Parameters_->die_if_not_converged = true;
    }
    else {
      Parameters_->convergence = 1e-4;
      Parameters_->energy_convergence = 1e-6;
      Parameters_->maxiter = 24;
      Parameters_->die_if_not_converged = true;
    }
    Parameters_->opdm = 0;
    Parameters_->tpdm = 0;
  }

  Parameters_->opdm_diag = 0;

  if (options["CI_MAXITER"].has_changed()) {
    Parameters_->maxiter = options.get_int("CI_MAXITER");
  }
  if (options["R_CONVERGENCE"].has_changed()) {
    Parameters_->convergence = options.get_double("R_CONVERGENCE");
  }
  if (options["E_CONVERGENCE"].has_changed()) {
    Parameters_->energy_convergence = options.get_double("E_CONVERGENCE");
  }
  if (options["DIE_IF_NOT_CONVERGED"].has_changed()) {
    Parameters_->die_if_not_converged = options.get_bool("DIE_IF_NOT_CONVERGED");
  }

  Parameters_->multp = molecule_->multiplicity();
  Parameters_->S = (((double) Parameters_->multp) - 1.0) / 2.0;
  if (options["S"].has_changed())
    Parameters_->S = options.get_double("S");

  /* specify OPENTYPE (SINGLET if open-shell singlet) */
  Parameters_->ref = options.get_str("REFERENCE");
  if (Parameters_->ref == "RHF") {
    Parameters_->opentype = PARM_OPENTYPE_NONE;
    Parameters_->Ms0 = 1;
  }
  else if (Parameters_->ref == "ROHF") {
    if (Parameters_->multp == 1) {
      Parameters_->opentype = PARM_OPENTYPE_SINGLET;
      Parameters_->Ms0 = 1;
    }
    else {
      Parameters_->opentype = PARM_OPENTYPE_HIGHSPIN;
      Parameters_->Ms0 = 0;
    }
  }  /* end ROHF parsing */
  else {
    throw InputException("Invalid DETCI reference " + Parameters_->ref,
      "REFERENCE", __FILE__, __LINE__);
  }

  if (options["MS0"].has_changed())
    Parameters_->Ms0 = options["MS0"].to_integer();

  Parameters_->h0blocksize = options.get_int("H0_BLOCKSIZE");
  Parameters_->h0guess_size = Parameters_->h0blocksize;
  if (options["H0_GUESS_SIZE"].has_changed())
    Parameters_->h0guess_size = options.get_int("H0_GUESS_SIZE");
  if (Parameters_->h0guess_size > Parameters_->h0blocksize)
    Parameters_->h0guess_size = Parameters_->h0blocksize;

  Parameters_->h0block_coupling_size = options.get_int("H0_BLOCK_COUPLING_SIZE");
  Parameters_->h0block_coupling = options["H0_BLOCK_COUPLING"].to_integer();

  Parameters_->nprint = options.get_int("NUM_DETS_PRINT");
  Parameters_->cc_nprint = options.get_int("NUM_AMPS_PRINT");

  if (options["FCI"].has_changed())
    Parameters_->fci = options["FCI"].to_integer();

  Parameters_->fci_strings = 0;
  if (Parameters_->fci) Parameters_->fci_strings = 1;
  if (options["FCI_STRINGS"].has_changed())
    Parameters_->fci_strings = options["FCI_STRINGS"].to_integer();

  Parameters_->mixed = options["MIXED"].to_integer();
  Parameters_->mixed4 = options["MIXED4"].to_integer();
  Parameters_->r4s = options["R4S"].to_integer();
  Parameters_->repl_otf = options["REPL_OTF"].to_integer();
  Parameters_->calc_ssq = options["CALC_S_SQUARED"].to_integer();

  if (options["MPN"].has_changed())
    Parameters_->mpn = options["MPN"].to_integer();

  // NOTE: In cases like this where the default value changes based
  // on other parameters, we need to make sure there's always a default
  // set here in this file, because we can't get the default from input
  // parsing --- the input parser is only called if the keyword is
  // present in the input file, because we're checking
  // options["XX"].has_changed() --CDS 8/2011
  if (Parameters_->mpn) {
    Parameters_->fci = 1;
    Parameters_->mpn_schmidt = FALSE;
    Parameters_->wigner = TRUE;
    Parameters_->guess_vector = PARM_GUESS_VEC_UNIT;
    Parameters_->hd_ave = ORB_ENER;
    Parameters_->update = UPDATE_DAVIDSON;
    Parameters_->hd_otf = TRUE;
    Parameters_->nodfile = TRUE;
  }
  else {
    Parameters_->update = UPDATE_DAVIDSON;
    Parameters_->mpn_schmidt = FALSE;
    Parameters_->wigner = FALSE;
    Parameters_->hd_otf = TRUE;
    Parameters_->nodfile = FALSE;
  }

  Parameters_->save_mpn2 = options.get_int("MPN_ORDER_SAVE");
  Parameters_->perturbation_parameter =
    options.get_double("PERTURB_MAGNITUDE");

  if (Parameters_->perturbation_parameter <= 1.0 &&
      Parameters_->perturbation_parameter >= -1.0) Parameters_->z_scale_H = 1;
/*
   else { outfile->Printf( "Parameters_->perturbation_parameters beyond the"
                 "bounds of -1.0 >= z <= 1.0\n");
         exit(0);
        }
*/
  if (options["MPN_SCHMIDT"].has_changed())
    Parameters_->mpn_schmidt = options["MPN_SCHMIDT"].to_integer();

  if (options["MPN_WIGNER"].has_changed())
    Parameters_->wigner = options["MPN_WIGNER"].to_integer();

  Parameters_->a_ras3_max = options.get_int("A_RAS3_MAX");
  Parameters_->b_ras3_max = options.get_int("B_RAS3_MAX");
  Parameters_->ras3_max = options.get_int("RAS3_MAX");
  Parameters_->ras4_max = options.get_int("RAS4_MAX");
  Parameters_->ras34_max = options.get_int("RAS34_MAX");

  Parameters_->cc_a_ras3_max = options.get_int("CC_A_RAS3_MAX");
  Parameters_->cc_b_ras3_max = options.get_int("CC_B_RAS3_MAX");
  Parameters_->cc_ras3_max = options.get_int("CC_RAS3_MAX");
  Parameters_->cc_ras4_max = options.get_int("CC_RAS4_MAX");
  Parameters_->cc_ras34_max = options.get_int("CC_RAS34_MAX");

  if (options["GUESS_VECTOR"].has_changed()) {
    std::string line1 = options.get_str("GUESS_VECTOR");
    if (line1 == "UNIT")
      Parameters_->guess_vector = PARM_GUESS_VEC_UNIT;
    else if (line1 == "H0_BLOCK")
      Parameters_->guess_vector = PARM_GUESS_VEC_H0_BLOCK;
    else if (line1 == "DFILE")
      Parameters_->guess_vector = PARM_GUESS_VEC_DFILE;
    else Parameters_->guess_vector = PARM_GUESS_VEC_UNIT;
  }

  Parameters_->icore = options.get_int("ICORE");

  if (options["HD_AVG"].has_changed()) {
    std::string line1 = options.get_str("HD_AVG");
    if (line1 == "HD_EXACT")    Parameters_->hd_ave = HD_EXACT;
    if (line1 == "HD_KAVE")     Parameters_->hd_ave = HD_KAVE;
    if (line1 == "ORB_ENER")    Parameters_->hd_ave = ORB_ENER;
    if (line1 == "EVANGELISTI") Parameters_->hd_ave = EVANGELISTI;
    if (line1 == "LEININGER")   Parameters_->hd_ave = LEININGER;
    if (line1 == "Z_KAVE")      Parameters_->hd_ave = Z_HD_KAVE;
    /* if (Parameters_->mpn) Parameters_->hd_ave = ORB_ENER; */
  }

  if (options["HD_OTF"].has_changed())
    Parameters_->hd_otf = options["HD_OTF"].to_integer();

  if (Parameters_->hd_otf == 0) {
    outfile->Printf( "Warning: HD_OTF FALSE has not been tested recently\n");
  }

  if (options["NO_DFILE"].has_changed())
    Parameters_->nodfile = options["NO_DFILE"].to_integer();
  if (Parameters_->num_roots > 1) Parameters_->nodfile = FALSE;

  Parameters_->diag_method = METHOD_DAVIDSON_LIU_SEM;
  if (options["DIAG_METHOD"].has_changed()) {
      std::string line1 = options.get_str("DIAG_METHOD");
      if (line1 == "RSP") {
          Parameters_->diag_method = METHOD_RSP;
          Parameters_->neg_only = 0;
      }
      else if (line1 == "OLSEN") {
          Parameters_->diag_method = METHOD_OLSEN;
      }
      else if (line1 == "MITRUSHENKOV") {
          Parameters_->diag_method = METHOD_MITRUSHENKOV;
      }
      else if (line1 == "DAVIDSON") {
          Parameters_->diag_method = METHOD_DAVIDSON_LIU_SEM;
      }
      else if (line1 == "SEM") {
          Parameters_->diag_method = METHOD_DAVIDSON_LIU_SEM;
      }
      else if (line1 == "SEMTEST") {
          Parameters_->diag_method = METHOD_RSPTEST_OF_SEM;
      }
  }

  if ((Parameters_->diag_method == METHOD_RSP) & (Parameters_->icore != 1)){
    outfile->Printf("RSP only works with icore = 1, switching.");
    Parameters_->icore = 1;
  }

  Parameters_->precon = PRECON_DAVIDSON;
  if (options["PRECONDITIONER"].has_changed()) {
    std::string line1 = options.get_str("PRECONDITIONER");
    if (line1 == "LANCZOS") Parameters_->precon = PRECON_LANCZOS;
    if (line1 == "DAVIDSON") Parameters_->precon = PRECON_DAVIDSON;
    if (line1 == "GEN_DAVIDSON")
      Parameters_->precon = PRECON_GEN_DAVIDSON;
    if (line1 == "H0BLOCK") Parameters_->precon = PRECON_GEN_DAVIDSON;
    if (line1 == "H0BLOCK_INV")
      Parameters_->precon = PRECON_H0BLOCK_INVERT;
    if (line1 == "ITER_INV")
      Parameters_->precon = PRECON_H0BLOCK_ITER_INVERT;
    if (line1 == "H0BLOCK_COUPLING")
      Parameters_->precon = PRECON_H0BLOCK_COUPLING;
    if (line1 == "EVANGELISTI")
      Parameters_->precon = PRECON_EVANGELISTI;
  }

  if (options["UPDATE"].has_changed()) {
    std::string line1 = options.get_str("UPDATE");
    if (line1 == "DAVIDSON") Parameters_->update = UPDATE_DAVIDSON;
    if (line1 == "OLSEN") Parameters_->update = UPDATE_OLSEN;
  }

  if (Parameters_->diag_method < METHOD_DAVIDSON_LIU_SEM &&
      Parameters_->update==UPDATE_DAVIDSON) {
    outfile->Printf("DAVIDSON update not available for OLSEN or MITRUSH"
            " iterators\n");
    Parameters_->update = UPDATE_OLSEN;
  }
  if (Parameters_->precon==PRECON_EVANGELISTI &&
      (Parameters_->update!=UPDATE_DAVIDSON
        || Parameters_->diag_method!=METHOD_DAVIDSON_LIU_SEM)) {
    outfile->Printf("EVANGELISTI preconditioner not available for OLSEN or"
                    " MITRUSH iterators or updates.\n");
    Parameters_->update = UPDATE_DAVIDSON;
  }

  // errcod = ip_boolean("ZERO_BLOCKS",&(Parameters_->zero_blocks),0);
  // if (Parameters_->icore || !Parameters_->mpn) Parameters_->zero_blocks = 0;
  Parameters_->num_init_vecs = Parameters_->num_roots;
  if (options["NUM_INIT_VECS"].has_changed())
    Parameters_->num_init_vecs = options.get_int("NUM_INIT_VECS");

  Parameters_->collapse_size = options.get_int("COLLAPSE_SIZE");
  if (Parameters_->collapse_size < 1) Parameters_->collapse_size = 1;

  Parameters_->lse = options["LSE"].to_integer();

  Parameters_->lse_collapse = options.get_int("LSE_COLLAPSE");
  if (Parameters_->lse_collapse < 1) Parameters_->lse_collapse = 3;

  Parameters_->lse_tolerance = options.get_double("LSE_TOLERANCE");

  Parameters_->maxnvect = options.get_int("MAX_NUM_VECS");

  if (Parameters_->maxnvect == 0 &&
      Parameters_->diag_method == METHOD_DAVIDSON_LIU_SEM) {
    Parameters_->maxnvect = Parameters_->maxiter * Parameters_->num_roots
      + Parameters_->num_init_vecs;
  }
  else if (Parameters_->maxnvect == 0 &&
           Parameters_->diag_method == METHOD_RSPTEST_OF_SEM) {
    Parameters_->maxnvect = Parameters_->maxiter * Parameters_->num_roots
      + Parameters_->num_init_vecs;
  }
  else if (Parameters_->maxnvect == 0 &&
           Parameters_->diag_method == METHOD_MITRUSHENKOV) {
    Parameters_->maxnvect = 2;
  }
  else if (Parameters_->maxnvect == 0 &&
           Parameters_->diag_method == METHOD_OLSEN) {
    Parameters_->maxnvect = 1;
  }
  else { /* the user tried to specify a value for maxnvect...check it */
  /*    if (Parameters_->maxnvect / (Parameters_->collapse_size *
        Parameters_->num_roots) < 2) {
        outfile->Printf( "maxnvect must be at least twice collapse_size *");
        outfile->Printf( " num_roots.\n");
        exit(0);
        }
  */
  }


  Parameters_->restart = options["RESTART"].to_integer();


  Parameters_->bendazzoli = options["BENDAZZOLI"].to_integer();
  if (Parameters_->bendazzoli && !Parameters_->fci) Parameters_->bendazzoli=0;

  /* Parse the OPDM stuff.  It is possible to give incompatible options,
   * but we will try to eliminate some of those.  Parameters_opdm will
   * function as the master switch for all other OPDM parameters.
   */
  Parameters_->opdm_diag = options["NAT_ORBS"].to_integer();

  if (Parameters_->opdm_diag || Parameters_->opdm_ave)
    Parameters_->opdm = 1;
  if (options["OPDM"].has_changed())
    Parameters_->opdm = options["OPDM"].to_integer();

  /* transition density matrices */
  Parameters_->transdens = 0;

  if (Parameters_->num_roots > 1)
    Parameters_->transdens = 1;
  else
    Parameters_->transdens = 0;

  if (options["TDM"].has_changed())
    Parameters_->transdens = options["TDM"].to_integer();

  /* dipole or transition dipole moment? */
  if (Parameters_->opdm) Parameters_->dipmom = 1;
  else Parameters_->dipmom = 0;

  if (Parameters_->transdens) Parameters_->dipmom = 1;

  if (Parameters_->wfn == "RASSCF" ||
      Parameters_->wfn == "CASSCF" ||
      Parameters_->wfn == "DETCAS")
    Parameters_->dipmom = 0;

  if (options["DIPMOM"].has_changed())
    Parameters_->dipmom = options["DIPMOM"].to_integer();

  if (Parameters_->dipmom == 1) Parameters_->opdm = 1;

  Parameters_->root = options.get_int("FOLLOW_ROOT");
  if (Parameters_->root < 0) Parameters_->root = 0;

  if (options["TPDM"].has_changed())
    Parameters_->tpdm = options["TPDM"].to_integer();

  if (Parameters_->num_init_vecs < Parameters_->num_roots)
    Parameters_->num_init_vecs = Parameters_->num_roots;
  if (Parameters_->guess_vector == PARM_GUESS_VEC_UNIT &&
      Parameters_->num_init_vecs > 1) {
    Parameters_->guess_vector = PARM_GUESS_VEC_H0_BLOCK;
    outfile->Printf("Warning: Unit vec option not available for more than"
            " one root\n");
  }
  if (Parameters_->guess_vector == PARM_GUESS_VEC_UNIT)
    Parameters_->h0blocksize = Parameters_->h0guess_size = 1;


  Parameters_->nthreads = Process::environment.get_n_threads();
  if (options["CI_NUM_THREADS"].has_changed()){
     Parameters_->nthreads = options.get_int("CI_NUM_THREADS");
  }
  if (Parameters_->nthreads < 1) Parameters_->nthreads = 1;

  Parameters_->sf_restrict = options["SF_RESTRICT"].to_integer();
  Parameters_->print_sigma_overlap = options["SIGMA_OVERLAP"].to_integer();

  if (Parameters_->cc) Parameters_->ex_lvl = Parameters_->cc_ex_lvl + 2;

  Parameters_->ex_allow.resize(Parameters_->ex_lvl);
  if (options["EX_ALLOW"].has_changed()) {
    i = options["EX_ALLOW"].size(); // CDS-TODO: Check that this really works
    if (i != Parameters_->ex_lvl) {
      std::string str = "Dim. of EX_ALLOW must be";
      str += std::to_string(Parameters_->ex_lvl);
      throw PsiException(str,__FILE__,__LINE__);
    }
    options.fill_int_array("EX_ALLOW", Parameters_->ex_allow.data());
  }
  else {
    for (i=0;i<Parameters_->ex_lvl;i++) {
      Parameters_->ex_allow[i] = 1;
    }
  }

  /* The filter_guess options are used to filter out some trial
     vectors which may not have the appropriate phase convention
     between two determinants.  This is useful to remove, e.g.,
     delta states when a sigma state is desired.  The user
     inputs two determinants (by giving the absolute alpha string
     number and beta string number for each), and also the
     desired phase between these two determinants for guesses
     which are to be kept.
   */

  Parameters_->filter_guess = options["FILTER_GUESS"].to_integer();

  if (Parameters_->filter_guess == 1) {
    Parameters_->filter_guess_sign = options.get_int("FILTER_GUESS_SIGN");
    if (Parameters_->filter_guess_sign != 1 &&
         Parameters_->filter_guess_sign != -1) {
      outfile->Printf( "FILTER_GUESS_SIGN should be 1 or -1 !\n");
      abort();
    }

    if (options["FILTER_GUESS_DET1"].size() != 2) {
      outfile->Printf( "Need to specify FILTER_GUESS_DET1 = "
                        "(alphastr betastr)\n");
      abort();
    }
    Parameters_->filter_guess_Ia = options["FILTER_GUESS_DET1"][0].to_integer();
    Parameters_->filter_guess_Ib = options["FILTER_GUESS_DET1"][1].to_integer();

    if (options["FILTER_GUESS_DET2"].size() != 2) {
      outfile->Printf( "Need to specify FILTER_GUESS_DET2 = "
                        "(alphastr betastr)\n");
      abort();
    }
    Parameters_->filter_guess_Ja = options["FILTER_GUESS_DET2"][0].to_integer();
    Parameters_->filter_guess_Jb = options["FILTER_GUESS_DET2"][1].to_integer();

  } /* end the filter_guess stuff */

  /* sometimes the filter guess stuff is not sufficient in that
     some states come in that we don't want.  We can help exclude
     them by explicitly zeroing out certain determinants, if that
     is correct for the desired state.  This stuff will allow the
     user to select a determinant which should always have a zero
     coefficient in the desired target state
   */
   Parameters_->filter_zero_det = 0;
  if (options["FILTER_ZERO_DET"].has_changed()) {
    if (options["FILTER_ZERO_DET"].size() != 2) {
      outfile->Printf( "Need to specify FILTER_ZERO_DET = "
                       "(alphastr betastr)\n");
      abort();
    }
    Parameters_->filter_zero_det = 1;
    Parameters_->filter_zero_det_Ia = options["FILTER_ZERO_DET"][0].to_integer();
    Parameters_->filter_zero_det_Ib = options["FILTER_ZERO_DET"][1].to_integer();
  }

  /* Does the user request a state-averaged calculation? */
  if (options["AVG_STATES"].has_changed()) {
      i = options["AVG_STATES"].size();
      if (i < 1 || i > Parameters_->num_roots) {
          std::string str = "Invalid number of states to average (";
          str += std::to_string(i);
          str += ")";
          throw PsiException(str, __FILE__, __LINE__);
      }

      Parameters_->average_states.resize(i);
      Parameters_->average_weights.resize(i);
      Parameters_->average_num = i;
      for (i = 0; i < Parameters_->average_num; i++) {
          Parameters_->average_states[i] = options["AVG_STATES"][i].to_integer();
          if (Parameters_->average_states[i] >= Parameters_->num_roots) {
              throw PsiException("Average state greater than the number of roots!", __FILE__, __LINE__);
          }
          if (Parameters_->average_states[i] < 0) {
              std::string str = "AVG_STATES start numbering from 0.\n";
              str += "Invalid state number ";
              str += std::to_string(Parameters_->average_states[i]);
              throw PsiException(str, __FILE__, __LINE__);
          }
          Parameters_->average_weights[i] = 1.0 / ((double)Parameters_->average_num);
      }

      if (options["AVG_WEIGHTS"].has_changed()) {
          if (options["AVG_WEIGHTS"].size() != Parameters_->average_num) {
              std::string str = "Mismatched number of average weights (";
              str += std::to_string(i);
              str += ")";
              throw PsiException(str, __FILE__, __LINE__);
          }
          for (i = 0; i < Parameters_->average_num; i++) {
              Parameters_->average_weights[i] = options["AVG_WEIGHTS"][i].to_double();
          }
      }

    if (Parameters_->average_num > 1) Parameters_->opdm_ave = 1;

    if ((!options["FOLLOW_ROOT"].has_changed()) && (Parameters_->average_num==1)) {
      Parameters_->root = Parameters_->average_states[0];
    }
  }

  else {
    Parameters_->average_num = 1;
    Parameters_->average_states.resize(1);
    Parameters_->average_weights.resize(1);
    Parameters_->average_states[0] = Parameters_->root;
    Parameters_->average_weights[0] = 1.0;
  } /* end state-average parsing */


  Parameters_->follow_vec_num = 0;
  if (options["FOLLOW_VECTOR"].has_changed()) {
    i = options["FOLLOW_VECTOR"].size();
    Parameters_->follow_vec_num = i;
    Parameters_->follow_vec_coef.resize(i);
    Parameters_->follow_vec_Ia.resize(i);
    Parameters_->follow_vec_Ib.resize(i);
    Parameters_->follow_vec_Iac.resize(i);
    Parameters_->follow_vec_Ibc.resize(i);
    Parameters_->follow_vec_Iaridx.resize(i);
    Parameters_->follow_vec_Ibridx.resize(i);

    /* now parse each piece */
    for (i=0; i<Parameters_->follow_vec_num; i++) {

      int isize = options["FOLLOW_VECTOR"][i].size();
      if (isize != 2) {
        outfile->Printf( "Need format FOLLOW_VECTOR = \n");
        outfile->Printf( "  [ [[alphastr_i, betastr_i], coeff_i], ... ] \n");
        abort();
      }

      int iisize = options["FOLLOW_VECTOR"][i][0].size();
      if (iisize != 2) {
        outfile->Printf( "Need format FOLLOW_VECTOR = \n");
        outfile->Printf( "  [ [[alphastr_i, betastr_i], coeff_i], ... ] \n");
        abort();
      }

      Parameters_->follow_vec_Ia[i] =
        options["FOLLOW_VECTOR"][i][0][0].to_integer();

      Parameters_->follow_vec_Ib[i] =
        options["FOLLOW_VECTOR"][i][0][1].to_integer();

      Parameters_->follow_vec_coef[i] =
        options["FOLLOW_VECTOR"][i][1].to_double();

    } /* end loop over parsing */
  } /* end follow vector stuff */


  /* make sure SA weights add up to 1.0 */
  for (i=0,junk=0.0; i<Parameters_->average_num; i++) {
    junk += Parameters_->average_weights[i];
  }
  if (junk <= 0.0) {
    std::string str = "Error: AVERAGE WEIGHTS add up to ";
    char*str2 = new char[25];
    sprintf(str2,"%20.15lf",junk);
    str += str2;
    delete[] str2;
    throw PsiException(str,__FILE__,__LINE__);
  }
  for (i=0; i<Parameters_->average_num; i++) {
    Parameters_->average_weights[i] /= junk;
  }

  Parameters_->cc_export = options["CC_VECS_WRITE"].to_integer();
  Parameters_->cc_import = options["CC_VECS_READ"].to_integer();
  Parameters_->cc_fix_external = options["CC_FIX_EXTERNAL"].to_integer();
  Parameters_->cc_fix_external_min = options.get_int("CC_FIX_EXTERNAL_MIN");
  Parameters_->cc_variational = options["CC_VARIATIONAL"].to_integer();
  Parameters_->cc_mixed = options["CC_MIXED"].to_integer();
  /* update using orb eigvals or not? */
  Parameters_->cc_update_eps = options["CC_UPDATE_EPS"].to_integer();

  /* DIIS only kicks in for CC anyway, no need to prefix with CC_ */
  Parameters_->diis  = options["DIIS"].to_integer();
  /* Iteration to turn on DIIS */
  Parameters_->diis_start  = options.get_int("DIIS_START_ITER");
  /* Do DIIS every n iterations */
  Parameters_->diis_freq  = options.get_int("DIIS_FREQ");
  Parameters_->diis_min_vecs  = options.get_int("DIIS_MIN_VECS");
  Parameters_->diis_max_vecs  = options.get_int("DIIS_MAX_VECS");
  Parameters_->mcscf_type  = options_.get_str("MCSCF_TYPE");
}


/*
** print_parameters(): Function prints the program's running parameters
**   found in the Parameters structure.
*/
void CIWavefunction::print_parameters(void) {
    int i;

    outfile->Printf("\n");
    outfile->Printf("   ==> Parameters <==\n\n");
    outfile->Printf("    EX LEVEL       =   %6d      H0 BLOCKSIZE  =   %6d\n",
                    Parameters_->ex_lvl, Parameters_->h0blocksize);
    outfile->Printf("    VAL EX LEVEL   =   %6d      H0 GUESS SIZE =   %6d\n",
                    Parameters_->val_ex_lvl, Parameters_->h0guess_size);
    outfile->Printf("    H0COUPLINGSIZE =   %6d      H0 COUPLING   =   %6s\n",
                    Parameters_->h0block_coupling_size,
                    Parameters_->h0block_coupling ? "YES" : "NO");
    outfile->Printf("    MAXITER        =   %6d      NUM PRINT     =   %6d\n",
                    Parameters_->maxiter, Parameters_->nprint);
    outfile->Printf("    NUM ROOTS      =   %6d      ICORE         =   %6d\n",
                    Parameters_->num_roots, Parameters_->icore);
    outfile->Printf("    PRINT LVL      =   %6d      FCI           =   %6s\n",
                    print_, Parameters_->fci ? "YES" : "NO");
    outfile->Printf("    R CONV         = %6.2e      MIXED         =   %6s\n",
                    Parameters_->convergence,
                    Parameters_->mixed ? "YES" : "NO");
    outfile->Printf("    E CONV         = %6.2e      MIXED4        =   %6s\n",
                    Parameters_->energy_convergence,
                    Parameters_->mixed4 ? "YES" : "NO");
    outfile->Printf("    R4S            =   %6s      REPL OTF      =   %6s\n",
                    Parameters_->r4s ? "YES" : "NO",
                    Parameters_->repl_otf ? "YES" : "NO");
    outfile->Printf("    DIAG METHOD    =   ");

    switch (Parameters_->diag_method) {
        case 0:
            outfile->Printf("%6s", "RSP");
            break;
        case 1:
            outfile->Printf("%6s", "OLSEN");
            break;
        case 2:
            outfile->Printf("%6s", "MITRUS");
            break;
        case 3:
            outfile->Printf("%6s", "SEM");
            break;
        case 4:
            outfile->Printf("%6s", "SEMTEST");
            break;
        default:
            outfile->Printf("%6s", "???");
            break;
    }
    outfile->Printf("      FOLLOW ROOT   =   %6d\n", Parameters_->root);

    outfile->Printf("    PRECONDITIONER = ");
    switch (Parameters_->precon) {
        case PRECON_LANCZOS:
            outfile->Printf("%6s", " LANCZOS    ");
            break;
        case PRECON_DAVIDSON:
            outfile->Printf("%6s", "DAVIDSON    ");
            break;
        case PRECON_GEN_DAVIDSON:
            outfile->Printf("%6s", "GEN_DAVIDSON");
            break;
        case PRECON_H0BLOCK_INVERT:
            outfile->Printf("%6s", "H0BLOCK_INV ");
            break;
        case PRECON_H0BLOCK_ITER_INVERT:
            outfile->Printf("%6s", "ITER_INV    ");
            break;
        case PRECON_H0BLOCK_COUPLING:
            outfile->Printf("%6s", "H0_COUPLING ");
            break;
        case PRECON_EVANGELISTI:
            outfile->Printf("%6s", "EVANGELISTI ");
            break;
        default:
            outfile->Printf("%6s", "???         ");
            break;
    }

    outfile->Printf("  UPDATE        = ");
    switch (Parameters_->update) {
        case 1:
            outfile->Printf("DAVIDSON\n");
            break;
        case 2:
            outfile->Printf("OLSEN\n");
            break;
        default:
            outfile->Printf("???\n");
            break;
    }

    outfile->Printf("    S              =   %.4lf      Ms0           =   %6s\n",
                    Parameters_->S, Parameters_->Ms0 ? "YES" : "NO");
    outfile->Printf("    GUESS VECTOR   =  ");
    switch (Parameters_->guess_vector) {
        case PARM_GUESS_VEC_UNIT:
            outfile->Printf("%7s", "UNIT");
            break;
        case PARM_GUESS_VEC_H0_BLOCK:
            outfile->Printf("%7s", "H0BLOCK");
            break;
        case PARM_GUESS_VEC_DFILE:
            outfile->Printf("%7s", "D FILE");
            break;
        default:
            outfile->Printf("%7s", "???");
            break;
    }
    outfile->Printf("      OPENTYPE      = ");
    switch (Parameters_->opentype) {
        case PARM_OPENTYPE_NONE:
            outfile->Printf("%8s\n", "NONE");
            break;
        case PARM_OPENTYPE_HIGHSPIN:
            outfile->Printf("%8s\n", "HIGHSPIN");
            break;
        case PARM_OPENTYPE_SINGLET:
            outfile->Printf("%8s\n", "SINGLET");
            break;
        default:
            outfile->Printf("%8s\n", "???");
            break;
    }

    outfile->Printf("    COLLAPSE SIZE  =   %6d", Parameters_->collapse_size);
    outfile->Printf("      HD AVG        = ");
    switch (Parameters_->hd_ave) {
        case HD_EXACT:
            outfile->Printf("HD_EXACT\n");
            break;
        case HD_KAVE:
            outfile->Printf("HD_KAVE\n");
            break;
        case ORB_ENER:
            outfile->Printf("ORB_ENER\n");
            break;
        case EVANGELISTI:
            outfile->Printf("EVANGELISTI\n");
            break;
        case LEININGER:
            outfile->Printf("LEININGER\n");
            break;
        default:
            outfile->Printf("???\n");
            break;
    }

    // HD OTF
    // if (Parameters_->hd_otf){
    //   outfile->Printf( "   HD OTF         =   %6s      NO DFILE      = %6s\n",
    //           Parameters_->hd_otf ? "YES" : "NO", Parameters_->nodfile ?
    //           "YES":"NO");
    // }

    if (Parameters_->mpn) {
        outfile->Printf("MPN options:\n");
        outfile->Printf("    MPN            =   %6s      MPN SCHMIDT   =   %6s\n",
                        Parameters_->mpn ? "YES" : "NO",
                        Parameters_->mpn_schmidt ? "YES" : "NO");
        outfile->Printf("    ZAPTN          =   %6s      MPN WIGNER    =   %6s\n",
                        Parameters_->zaptn ? "YES" : "NO",
                        Parameters_->wigner ? "YES" : "NO");
        outfile->Printf("    PERT Z         =   %1.4f\n",
                        Parameters_->perturbation_parameter);
    }

    if (Parameters_->filter_guess || Parameters_->sf_restrict) {
        outfile->Printf("    FILTER GUESS   =   %6s      SF RESTRICT   =   %6s\n",
                        Parameters_->filter_guess ? "YES" : "NO",
                        Parameters_->sf_restrict ? "YES" : "NO");
    }
    if (Parameters_->cc && Parameters_->diis) {
        outfile->Printf("CC options:\n");
        outfile->Printf("    DIIS START     =   %6d      DIIS FREQ     =   %6d\n",
                        Parameters_->diis_start, Parameters_->diis_freq);
        outfile->Printf("    DIIS MIN VECS  =   %6d      DIIS MAX VECS=   %6d\n",
                        Parameters_->diis_min_vecs, Parameters_->diis_max_vecs);
    }
    if (Parameters_->restart){
      outfile->Printf("    RESTART        =   %6s\n",
                      Parameters_->restart ? "YES" : "NO");
    }
    outfile->Printf("    MAX NUM VECS   =   %6d   ", Parameters_->maxnvect);
    if (Parameters_->ref_sym == -1){
        outfile->Printf("   REF SYM       =   %6s\n", "AUTO");
    } else {
        outfile->Printf("   REF SYM       =   %6d\n", Parameters_->ref_sym);
    }
    outfile->Printf("    IOPEN        =   %6s\n", CalcInfo_->iopen ? "YES" : "NO");

    // LSE
    if (Parameters_->lse) {
        outfile->Printf("    LSE            =   %6s      LSE ITER      =   %6d\n",
                        Parameters_->lse ? "YES" : "NO", Parameters_->lse_iter);
    }
    // outfile->Printf( "   OPDM           =   %6s      TRANS DENSITY=   %6s\n",
    //         Parameters_->opdm ?  "YES":"NO",
    //         Parameters_->transdens ? "YES":"NO");
    // outfile->Printf( "\n   FILES          = %3d %2d %2d %2d\n",
    //    Parameters_->hd_filenum, Parameters_->c_filenum,
    //    Parameters_->s_filenum, Parameters_->d_filenum);

    outfile->Printf("\n    EX ALLOW       = ");
    for (i  = 0; i < Parameters_->ex_lvl; i++) {
        outfile->Printf("%2d ", Parameters_->ex_allow[i]);
    }

    outfile->Printf("\n    STATE AVERAGE  = ");
    for (i  = 0; i < Parameters_->average_num; i++) {
        if (i % 5 == 0 && i != 0) outfile->Printf("\n");
        outfile->Printf("%2d(%4.2lf) ", Parameters_->average_states[i],
                        Parameters_->average_weights[i]);
    }

    if (Parameters_->follow_vec_num > 0) {
        outfile->Printf("\nDensity matrices will follow vector like:\n");
        for (i  = 0; i < Parameters_->follow_vec_num; i++)
            outfile->Printf("(%d %d) %12.6lf\n", Parameters_->follow_vec_Ia[i],
                            Parameters_->follow_vec_Ib[i],
                            Parameters_->follow_vec_coef[i]);
    }

    outfile->Printf("\n\n");
}

/*
** set_ras_parms(): Set the RAS parameters or their conventional equivalents
**   (i.e. fermi level, etc).
**
*/
void CIWavefunction::set_ras_parameters(void) {
    int i, j, cnt;
    int errcod;
    int tot_expl_el, nras2alp, nras2bet, betsocc;
    int *ras1, *ras2, *ras3;
    int *orbsym;

    /* If the user asked for FCI=true, then override the other keywords
       if necessary to ensure that it's really a FCI
     */
    if (Parameters_->fci == 1 &&
        (CalcInfo_->num_alp_expl + CalcInfo_->num_bet_expl) >
            Parameters_->ex_lvl) {
        Parameters_->val_ex_lvl = 0;
        Parameters_->ex_lvl = CalcInfo_->num_alp_expl + CalcInfo_->num_bet_expl;
        Parameters_->ex_allow.clear();
        Parameters_->ex_allow.resize(Parameters_->ex_lvl);
        for (i = 0; i < Parameters_->ex_lvl; i++) Parameters_->ex_allow[i] = 1;

        if (print_ > 2) {
            outfile->Printf("Note: Calculation requested is a full CI.\n");
            outfile->Printf(
                "Resetting EX_LEVEL to %d and turning on all excitations\n\n",
                Parameters_->ex_lvl);
        }

    } /* end FCI override */

    /* reset ex_lvl if incompatible with number of electrons */
    if (Parameters_->cc &&
        (Parameters_->cc_ex_lvl >
         CalcInfo_->num_alp_expl + CalcInfo_->num_bet_expl)) {
        Parameters_->cc_ex_lvl =
            CalcInfo_->num_alp_expl + CalcInfo_->num_bet_expl;
    }
    if (Parameters_->ex_lvl >
        CalcInfo_->num_alp_expl + CalcInfo_->num_bet_expl) {
        Parameters_->ex_lvl = CalcInfo_->num_alp_expl + CalcInfo_->num_bet_expl;
    }

    for (i = 0, j = 0; i < CalcInfo_->nirreps; i++)
        j += CalcInfo_->ras_opi[0][i];
    Parameters_->a_ras1_lvl = Parameters_->b_ras1_lvl = Parameters_->ras1_lvl =
        j - 1;

    /* figure out how many electrons are in RAS II */
    /* alpha electrons */
    for (i = 0, nras2alp = 0, betsocc = 0; i < CalcInfo_->nirreps; i++) {
        j = CalcInfo_->docc[i] - CalcInfo_->dropped_docc[i] -
            CalcInfo_->ras_opi[0][i];
        if (Parameters_->opentype == PARM_OPENTYPE_HIGHSPIN) {
            j += CalcInfo_->socc[i];
        } else if (Parameters_->opentype == PARM_OPENTYPE_SINGLET) {
            if (betsocc + CalcInfo_->socc[i] <= CalcInfo_->spab)
                betsocc += CalcInfo_->socc[i];
            else {
                j += CalcInfo_->socc[i] - (CalcInfo_->spab - betsocc);
                betsocc = CalcInfo_->spab;
            }
        }
        if (j > 0) nras2alp += j;
        if (j > CalcInfo_->ras_opi[1][i]) {
            outfile->Printf("(set_ras_parms): detecting %d electrons ",
                            j - CalcInfo_->ras_opi[1][i]);
            outfile->Printf("in RAS III for irrep %d.\n", i);
            outfile->Printf(
                "Some parts of DETCI assume all elec in I and II\n");
        }
    }
    /* beta electrons */
    for (i = 0, nras2bet = 0, betsocc = 0; i < CalcInfo_->nirreps; i++) {
        j = CalcInfo_->docc[i] - CalcInfo_->dropped_docc[i] -
            CalcInfo_->ras_opi[0][i];
        if (Parameters_->opentype == PARM_OPENTYPE_SINGLET &&
            CalcInfo_->socc[i]) {
            if (betsocc + CalcInfo_->socc[i] <= CalcInfo_->spab)
                j += CalcInfo_->socc[i];
            else {
                j += CalcInfo_->spab - betsocc;
                betsocc = CalcInfo_->spab;
            }
        }
        if (j > 0) nras2bet += j;
        if (j > CalcInfo_->ras_opi[1][i]) {
            outfile->Printf("(set_ras_parms): detecting %d electrons ",
                            j - CalcInfo_->ras_opi[1][i]);
            outfile->Printf("in RAS III for irrep %d.\n", i);
            outfile->Printf(
                "Some parts of DETCI assume all elec in I and II\n");
        }
    }

    Parameters_->a_ras1_max =
        (CalcInfo_->num_alp_expl > Parameters_->a_ras1_lvl + 1)
            ? Parameters_->a_ras1_lvl + 1
            : (CalcInfo_->num_alp_expl);
    // CDS 4/15
    // if (Parameters_->fzc) Parameters_->a_ras1_max += CalcInfo_->num_fzc_orbs;

    Parameters_->b_ras1_max =
        (CalcInfo_->num_bet_expl > Parameters_->b_ras1_lvl + 1)
            ? Parameters_->b_ras1_lvl + 1
            : (CalcInfo_->num_bet_expl);
    // CDS 4/15
    // if (Parameters_->fzc) Parameters_->b_ras1_max += CalcInfo_->num_fzc_orbs;

    for (i = 0, j = 0; i < CalcInfo_->nirreps; i++)
        j += CalcInfo_->ras_opi[1][i];
    Parameters_->ras3_lvl = Parameters_->ras1_lvl + j + 1;

    for (i = 0, j = 0; i < CalcInfo_->nirreps; i++)
        j += CalcInfo_->ras_opi[2][i];
    Parameters_->ras4_lvl = Parameters_->ras3_lvl + j;

    /* check Parameters to make sure everything consistent */

    if (Parameters_->cc) {
        if (Parameters_->cc_a_val_ex_lvl == -1)
            Parameters_->cc_a_val_ex_lvl = Parameters_->cc_val_ex_lvl;
        if (Parameters_->cc_b_val_ex_lvl == -1)
            Parameters_->cc_b_val_ex_lvl = Parameters_->cc_val_ex_lvl;
        if (Parameters_->cc_a_val_ex_lvl >
            Parameters_->ras3_lvl - Parameters_->ras1_lvl - 1)
            Parameters_->cc_a_val_ex_lvl =
                Parameters_->ras3_lvl - Parameters_->ras1_lvl - 1;
        if (Parameters_->cc_b_val_ex_lvl >
            Parameters_->ras3_lvl - Parameters_->ras1_lvl - 1)
            Parameters_->cc_b_val_ex_lvl =
                Parameters_->ras3_lvl - Parameters_->ras1_lvl - 1;
        if (Parameters_->cc_val_ex_lvl >
            Parameters_->cc_a_val_ex_lvl + Parameters_->cc_b_val_ex_lvl)
            Parameters_->cc_val_ex_lvl =
                Parameters_->cc_a_val_ex_lvl + Parameters_->cc_b_val_ex_lvl;
    }

    /* deduce Parameters_->cc_a_ras3_max and Parameters_->cc_b_ras3_max if
     * needed */
    if (Parameters_->cc & (Parameters_->cc_a_ras3_max == -1 ||
                           Parameters_->cc_b_ras3_max == -1)) {
        if (Parameters_->cc_ras3_max != -1) { /* have parsed cc_ras3_max */
            Parameters_->cc_a_ras3_max =
                (Parameters_->cc_ras3_max <= CalcInfo_->num_alp_expl)
                    ? Parameters_->cc_ras3_max
                    : CalcInfo_->num_alp_expl;
            Parameters_->cc_b_ras3_max =
                (Parameters_->cc_ras3_max <= CalcInfo_->num_bet_expl)
                    ? Parameters_->cc_ras3_max
                    : CalcInfo_->num_bet_expl;
        } else {
            Parameters_->cc_a_ras3_max =
                (Parameters_->cc_ex_lvl <= CalcInfo_->num_alp_expl)
                    ? Parameters_->cc_ex_lvl
                    : CalcInfo_->num_alp_expl;
            Parameters_->cc_b_ras3_max =
                (Parameters_->cc_ex_lvl <= CalcInfo_->num_bet_expl)
                    ? Parameters_->cc_ex_lvl
                    : CalcInfo_->num_bet_expl;
        }
    }

    if (Parameters_->cc) {
        Parameters_->a_ras3_max =
            (Parameters_->cc_a_ras3_max + 2 <= CalcInfo_->num_alp_expl)
                ? Parameters_->cc_a_ras3_max + 2
                : CalcInfo_->num_alp_expl;
        Parameters_->b_ras3_max =
            (Parameters_->cc_b_ras3_max + 2 <= CalcInfo_->num_bet_expl)
                ? Parameters_->cc_b_ras3_max + 2
                : CalcInfo_->num_bet_expl;
    }

    if (Parameters_->a_ras3_max == -1 || Parameters_->b_ras3_max == -1) {
        if (Parameters_->ras3_max != -1) { /* have parsed ras3_max */
            Parameters_->a_ras3_max =
                (Parameters_->ras3_max <= CalcInfo_->num_alp_expl)
                    ? Parameters_->ras3_max
                    : CalcInfo_->num_alp_expl;
            Parameters_->b_ras3_max =
                (Parameters_->ras3_max <= CalcInfo_->num_bet_expl)
                    ? Parameters_->ras3_max
                    : CalcInfo_->num_bet_expl;
        } else {
            Parameters_->a_ras3_max =
                (Parameters_->ex_lvl <= CalcInfo_->num_alp_expl)
                    ? Parameters_->ex_lvl
                    : CalcInfo_->num_alp_expl;
            Parameters_->b_ras3_max =
                (Parameters_->ex_lvl <= CalcInfo_->num_bet_expl)
                    ? Parameters_->ex_lvl
                    : CalcInfo_->num_bet_expl;
        }
    }

    if (Parameters_->cc) {
        if (Parameters_->cc_ras4_max != -1) { /* have parsed */
            Parameters_->cc_a_ras4_max =
                (Parameters_->cc_ras4_max <= CalcInfo_->num_alp_expl)
                    ? Parameters_->cc_ras4_max
                    : CalcInfo_->num_alp_expl;
            Parameters_->cc_b_ras4_max =
                (Parameters_->cc_ras4_max <= CalcInfo_->num_bet_expl)
                    ? Parameters_->cc_ras4_max
                    : CalcInfo_->num_bet_expl;
        } else {
            Parameters_->cc_a_ras4_max = Parameters_->cc_a_ras3_max;
            Parameters_->cc_b_ras4_max = Parameters_->cc_b_ras3_max;
        }
        Parameters_->a_ras4_max =
            (Parameters_->cc_a_ras4_max + 2 <= CalcInfo_->num_alp_expl)
                ? Parameters_->cc_a_ras4_max + 2
                : CalcInfo_->num_alp_expl;
        Parameters_->b_ras4_max =
            (Parameters_->cc_b_ras4_max + 2 <= CalcInfo_->num_bet_expl)
                ? Parameters_->cc_b_ras4_max + 2
                : CalcInfo_->num_bet_expl;
    } else {
        if (Parameters_->ras4_max != -1) { /* have parsed */
            Parameters_->a_ras4_max =
                (Parameters_->ras4_max <= CalcInfo_->num_alp_expl)
                    ? Parameters_->ras4_max
                    : CalcInfo_->num_alp_expl;
            Parameters_->b_ras4_max =
                (Parameters_->ras4_max <= CalcInfo_->num_bet_expl)
                    ? Parameters_->ras4_max
                    : CalcInfo_->num_bet_expl;
        } else {
            Parameters_->a_ras4_max = Parameters_->a_ras3_max;
            Parameters_->b_ras4_max = Parameters_->b_ras3_max;
        }
    }

    if (Parameters_->cc) {
        if (Parameters_->cc_ras34_max != -1) { /* have parsed */
            Parameters_->cc_a_ras34_max = Parameters_->cc_ras34_max;
            Parameters_->cc_b_ras34_max = Parameters_->cc_ras34_max;
        } else {
            Parameters_->cc_a_ras34_max =
                Parameters_->cc_a_ras3_max + Parameters_->cc_a_ras4_max;
            Parameters_->cc_b_ras34_max =
                Parameters_->cc_b_ras3_max + Parameters_->cc_b_ras4_max;
        }
        if (Parameters_->ras34_max != -1) { /* have parsed */
            Parameters_->a_ras34_max = Parameters_->ras34_max;
            Parameters_->b_ras34_max = Parameters_->ras34_max;
        } else {
            Parameters_->a_ras34_max = Parameters_->cc_a_ras34_max + 2;
            if (Parameters_->a_ras34_max > CalcInfo_->num_alp_expl)
                Parameters_->a_ras34_max = CalcInfo_->num_alp_expl;
            Parameters_->b_ras34_max = Parameters_->cc_b_ras34_max + 2;
            if (Parameters_->b_ras34_max > CalcInfo_->num_bet_expl)
                Parameters_->b_ras34_max = CalcInfo_->num_bet_expl;
            Parameters_->ras34_max = Parameters_->cc_ras34_max + 2;
            if (Parameters_->ras34_max >
                Parameters_->a_ras34_max + Parameters_->b_ras34_max)
                Parameters_->ras34_max =
                    Parameters_->a_ras34_max + Parameters_->b_ras34_max;
        }
    } else {                                /* non-CC */
        if (Parameters_->ras34_max != -1) { /* have parsed */
            Parameters_->a_ras34_max = Parameters_->ras34_max;
            Parameters_->b_ras34_max = Parameters_->ras34_max;
        } else {
            Parameters_->a_ras34_max = Parameters_->a_ras3_max;
            Parameters_->b_ras34_max = Parameters_->b_ras3_max;
        }
    }

    i = Parameters_->ras4_lvl - Parameters_->ras3_lvl;
    if (Parameters_->a_ras3_max > i) Parameters_->a_ras3_max = i;
    if (Parameters_->b_ras3_max > i) Parameters_->b_ras3_max = i;
    if (Parameters_->cc) {
        if (Parameters_->cc_a_ras3_max > i) Parameters_->cc_a_ras3_max = i;
        if (Parameters_->cc_b_ras3_max > i) Parameters_->cc_b_ras3_max = i;
    }

    i = CalcInfo_->num_ci_orbs - Parameters_->ras4_lvl;
    if (Parameters_->a_ras4_max > i) Parameters_->a_ras4_max = i;
    if (Parameters_->b_ras4_max > i) Parameters_->b_ras4_max = i;
    if (Parameters_->cc) {
        if (Parameters_->cc_a_ras4_max > i) Parameters_->cc_a_ras4_max = i;
        if (Parameters_->cc_b_ras4_max > i) Parameters_->cc_b_ras4_max = i;
    }

    i = CalcInfo_->num_ci_orbs - Parameters_->ras3_lvl;
    if (Parameters_->a_ras34_max > i) Parameters_->a_ras34_max = i;
    if (Parameters_->b_ras34_max > i) Parameters_->b_ras34_max = i;
    if (Parameters_->cc) {
        if (Parameters_->cc_a_ras34_max > i) Parameters_->cc_a_ras34_max = i;
        if (Parameters_->cc_b_ras34_max > i) Parameters_->cc_b_ras34_max = i;
    }

    i = (CalcInfo_->num_alp_expl <= Parameters_->a_ras1_lvl + 1)
            ? CalcInfo_->num_alp_expl
            : Parameters_->a_ras1_lvl + 1;
    Parameters_->a_ras1_min = i - Parameters_->ex_lvl - Parameters_->val_ex_lvl;
    if (Parameters_->a_ras1_min < 0) Parameters_->a_ras1_min = 0;
    // CDS 4/15 no longer include dropped core in ras1_min
    // Parameters_->a_ras1_min += CalcInfo_->num_fzc_orbs;
    Parameters_->a_ras1_min += CalcInfo_->num_expl_cor_orbs;

    i = (CalcInfo_->num_bet_expl <= Parameters_->b_ras1_lvl + 1)
            ? CalcInfo_->num_bet_expl
            : Parameters_->b_ras1_lvl + 1;
    Parameters_->b_ras1_min = i - Parameters_->ex_lvl - Parameters_->val_ex_lvl;
    if (Parameters_->b_ras1_min < 0) Parameters_->b_ras1_min = 0;
    // CDS 4/15 no longer include dropped core in ras1_min
    // Parameters_->b_ras1_min += CalcInfo_->num_fzc_orbs;
    Parameters_->b_ras1_min += CalcInfo_->num_expl_cor_orbs;

    tot_expl_el = CalcInfo_->num_alp_expl + CalcInfo_->num_bet_expl;
    if (Parameters_->cc) {
        if (Parameters_->cc_val_ex_lvl != 0)
            i = Parameters_->cc_val_ex_lvl;
        else
            i = Parameters_->cc_ex_lvl;
        if (Parameters_->cc_ras3_max == -1) {
            Parameters_->cc_ras3_max = (i <= tot_expl_el) ? i : tot_expl_el;
        } else {
            if (Parameters_->cc_ras3_max > tot_expl_el)
                Parameters_->cc_ras3_max = tot_expl_el;
        }
        if (Parameters_->ras3_max == -1)
            Parameters_->ras3_max = Parameters_->cc_ras3_max + 2;
    }
    if (Parameters_->ras3_max == -1) {
        Parameters_->ras3_max = (Parameters_->ex_lvl <= tot_expl_el)
                                    ? Parameters_->ex_lvl
                                    : tot_expl_el;
    } else {
        if (Parameters_->ras3_max > tot_expl_el)
            Parameters_->ras3_max = tot_expl_el;
    }

    i = 2 * (Parameters_->ras4_lvl - Parameters_->ras3_lvl);
    if (i < Parameters_->ras3_max) Parameters_->ras3_max = i;
    if (Parameters_->cc) {
        if (i < Parameters_->cc_ras3_max) Parameters_->cc_ras3_max = i;
    }

    i = (tot_expl_el < 2 * (Parameters_->ras1_lvl + 1))
            ? tot_expl_el
            : 2 * (Parameters_->ras1_lvl + 1);

    // CDS 4/15 no longer include dropped core in ras1_min
    // Parameters_->ras1_min = i - Parameters_->ex_lvl -
    //   Parameters_->val_ex_lvl + 2 * CalcInfo_->num_fzc_orbs;
    Parameters_->ras1_min = i - Parameters_->ex_lvl - Parameters_->val_ex_lvl;

    if (Parameters_->a_ras1_min + Parameters_->b_ras1_min >
        Parameters_->ras1_min)
        Parameters_->ras1_min =
            Parameters_->a_ras1_min + Parameters_->b_ras1_min;

    if (Parameters_->cc && Parameters_->cc_ras4_max == -1) {
        Parameters_->cc_ras4_max = (Parameters_->cc_ex_lvl <= tot_expl_el)
                                       ? Parameters_->cc_ex_lvl
                                       : tot_expl_el;
    }

    if (Parameters_->ras4_max == -1) {
        if (Parameters_->cc) {
            Parameters_->ras4_max =
                (Parameters_->cc_ras4_max + 2 <= tot_expl_el)
                    ? Parameters_->cc_ras4_max + 2
                    : tot_expl_el;
        } else
            Parameters_->ras4_max = (Parameters_->ex_lvl <= tot_expl_el)
                                        ? Parameters_->ex_lvl
                                        : tot_expl_el;
    }

    i = 2 * (CalcInfo_->num_ci_orbs - Parameters_->ras4_lvl);
    if (i < Parameters_->ras4_max) Parameters_->ras4_max = i;
    if (Parameters_->cc) {
        if (i < Parameters_->cc_ras4_max) Parameters_->cc_ras4_max = i;
    }

    if (Parameters_->cc && Parameters_->cc_ras34_max == -1)
        Parameters_->cc_ras34_max =
            Parameters_->cc_ras3_max + Parameters_->cc_ras4_max;
    i = 2 * (CalcInfo_->num_ci_orbs - Parameters_->ras3_lvl);
    if (Parameters_->cc) {
        if (i < Parameters_->cc_ras34_max) Parameters_->cc_ras34_max = i;
    }

    if (Parameters_->ras34_max == -1 && !Parameters_->cc)
        Parameters_->ras34_max = Parameters_->ras3_max;
    else
        Parameters_->ras34_max = Parameters_->cc_ras34_max + 2;

    i = 2 * (CalcInfo_->num_ci_orbs - Parameters_->ras3_lvl);
    if (i < Parameters_->ras34_max) Parameters_->ras34_max = i;
    if (Parameters_->ras34_max > tot_expl_el)
        Parameters_->ras34_max = tot_expl_el;

    if (Parameters_->a_ras34_max >
        Parameters_->a_ras3_max + Parameters_->a_ras4_max)
        Parameters_->a_ras34_max =
            Parameters_->a_ras3_max + Parameters_->a_ras4_max;
    if (Parameters_->cc &&
        (Parameters_->cc_a_ras34_max >
         Parameters_->cc_a_ras3_max + Parameters_->cc_b_ras4_max))
        Parameters_->cc_a_ras34_max =
            Parameters_->cc_a_ras3_max + Parameters_->cc_a_ras4_max;

    if (Parameters_->b_ras34_max >
        Parameters_->b_ras3_max + Parameters_->b_ras4_max)
        Parameters_->b_ras34_max =
            Parameters_->b_ras3_max + Parameters_->b_ras4_max;
    if (Parameters_->cc &&
        (Parameters_->cc_b_ras34_max >
         Parameters_->cc_b_ras3_max + Parameters_->cc_b_ras4_max))
        Parameters_->cc_b_ras34_max =
            Parameters_->cc_b_ras3_max + Parameters_->cc_b_ras4_max;

    /* now just re-check some basic things */
    if (Parameters_->a_ras34_max > CalcInfo_->num_alp_expl)
        Parameters_->a_ras34_max = CalcInfo_->num_alp_expl;
    if (Parameters_->b_ras34_max > CalcInfo_->num_bet_expl)
        Parameters_->b_ras34_max = CalcInfo_->num_bet_expl;
    if (Parameters_->cc) {
        if (Parameters_->cc_a_ras34_max > CalcInfo_->num_alp_expl)
            Parameters_->cc_a_ras34_max = CalcInfo_->num_alp_expl;
        if (Parameters_->cc_b_ras34_max > CalcInfo_->num_bet_expl)
            Parameters_->cc_b_ras34_max = CalcInfo_->num_bet_expl;
    }

    if (Parameters_->ras34_max >
        Parameters_->a_ras34_max + Parameters_->b_ras34_max)
        Parameters_->ras34_max =
            Parameters_->a_ras34_max + Parameters_->b_ras34_max;
}

std::string _concat_dim(std::string label, size_t spacer1, Dimension dim, size_t spacer2, size_t spacer3){
  std::stringstream ret;
  ret << std::setw(spacer1) << label;
  ret << std::setw(spacer2) << dim.sum();
  for (size_t h=0; h < dim.n(); h++ ) {
    ret << std::setw(spacer3) << dim[h];
  };
  ret << std::endl;
  return ret.str();
}
/*
** print_ras_parms(): Set the RAS parameters or their conventional equivalents
**   (i.e. fermi level, etc).
**
*/
void CIWavefunction::print_ras_parameters(void) {

  outfile->Printf("   ==> CI Orbital and Space information <==\n\n");

  // Print out any specific space info
  if (!Parameters_->fci & !Parameters_->cc) {
      outfile->Printf("    RAS1 LVL      =   %6d      A RAS3 MAX    =   %6d\n",
                      Parameters_->ras1_lvl, Parameters_->a_ras3_max);
      outfile->Printf("    RAS1 MIN      =   %6d      B RAS3 MAX    =   %6d\n",
                      Parameters_->ras1_min, Parameters_->b_ras3_max);
      outfile->Printf("    A RAS1 LVL    =   %6d      RAS4 LVL      =   %6d\n",
                      Parameters_->a_ras1_lvl, Parameters_->ras4_lvl);
      outfile->Printf("    A RAS1 MIN    =   %6d      A RAS4 MAX    =   %6d\n",
                      Parameters_->a_ras1_min, Parameters_->a_ras4_max);
      outfile->Printf("    A RAS1 MAX    =   %6d      B RAS4 MAX    =   %6d\n",
                      Parameters_->a_ras1_max, Parameters_->b_ras4_max);
      outfile->Printf("    B RAS1 LVL    =   %6d      RAS4 MAX      =   %6d\n",
                      Parameters_->b_ras1_lvl, Parameters_->ras4_max);
      outfile->Printf("    B RAS1 MIN    =   %6d      A RAS34 MAX   =   %6d\n",
                      Parameters_->b_ras1_min, Parameters_->a_ras34_max);
      outfile->Printf("    B RAS1 MAX    =   %6d      B RAS34 MAX   =   %6d\n",
                      Parameters_->b_ras1_max, Parameters_->b_ras34_max);
      outfile->Printf("    RAS3 LVL      =   %6d      RAS34 MAX     =   %6d\n",
                      Parameters_->ras3_lvl, Parameters_->ras34_max);
      outfile->Printf("    RAS3 MAX      =   %6d\n", Parameters_->ras3_max);
  }
  if (Parameters_->cc) {
      outfile->Printf("    CC RAS3 MAX   =   %6d      CC RAS4 MAX   =   %6d\n",
                      Parameters_->cc_ras3_max, Parameters_->cc_ras4_max);
      outfile->Printf("    CC A RAS3 MAX =   %6d      CC B RAS3 MAX =   %6d\n",
                      Parameters_->cc_a_ras3_max, Parameters_->cc_b_ras3_max);
      outfile->Printf("    CC A RAS4 MAX =   %6d      CC B RAS4 MAX =   %6d\n",
                      Parameters_->cc_a_ras4_max, Parameters_->cc_b_ras4_max);
      outfile->Printf("    CC RAS34 MAX  =   %6d\n", Parameters_->cc_ras34_max);
      outfile->Printf("    CC A RAS34 MAX=  %6d      CC B RAS34 MAX =  %6d\n",
                      Parameters_->cc_a_ras34_max,
                      Parameters_->cc_b_ras34_max);
      outfile->Printf("    CC MIXED      =   %6s      CC FIX EXTERN  =  %6s\n",
                      Parameters_->cc_mixed ? "YES"  : "NO",
                      Parameters_->cc_fix_external ? "YES"  : "NO");
      outfile->Printf("   CC VARIATIONAL =  %6s\n",
                      Parameters_->cc_variational ? "YES"  : "NO");
  }

  CharacterTable ct = molecule_->point_group()->char_table();
  size_t sdist = 20;
  size_t tdist = 9;
  size_t hdist = 6;

  // Build header
  std::stringstream headers;
  headers << std::setw(sdist) << "Space";
  headers << std::setw(tdist) << "Total";
  for (size_t h=0; h<nirrep_; h++){
    headers << std::setw(hdist) << ct.gamma(h).symbol();
  }
  headers << std::endl;
  std::string header = headers.str();

  // Figure out seperator size
  std::string sep_template =
    "-------------------------------------------------------------------------------------";
  std::string seperator = "   " + sep_template.substr(0, header.size()) + "\n";

  // Start building orbital info
  std::stringstream orbital_info;
  orbital_info << seperator;
  orbital_info << header;

  // General info
  orbital_info << seperator;
  orbital_info << _concat_dim("Nso", sdist, nsopi_, tdist, hdist);
  orbital_info << _concat_dim("Nmo", sdist, nmopi_, tdist, hdist);
  orbital_info << _concat_dim("Ndocc", sdist, doccpi_, tdist, hdist);
  orbital_info << _concat_dim("Nsocc", sdist, soccpi_, tdist, hdist);
  orbital_info << seperator;

  // Occupied spaces
  if (Parameters_->mcscf) {
    orbital_info << std::setw(header.size()/ 2 + 10) << "MCSCF Spaces\n";
    orbital_info << seperator;
    orbital_info << _concat_dim("Frozen DOCC", sdist, get_dimension("FZC"), tdist, hdist);
    orbital_info << _concat_dim("Restricted DOCC", sdist, get_dimension("DOCC"), tdist, hdist);
  } else {
    orbital_info << std::setw(header.size()/ 2 + 6) << "CI Spaces\n";
    orbital_info << seperator;
    orbital_info << _concat_dim("Dropped DOCC", sdist, get_dimension("DRC"), tdist, hdist);
  }

  // Active spaces
  if (Parameters_->fci) {
    orbital_info << _concat_dim("Active", sdist, get_dimension("ACT"), tdist, hdist);
  } else {
    orbital_info << _concat_dim("RAS1", sdist, get_dimension("RAS1"), tdist, hdist);
    orbital_info << _concat_dim("RAS2", sdist, get_dimension("RAS2"), tdist, hdist);
    orbital_info << _concat_dim("RAS3", sdist, get_dimension("RAS3"), tdist, hdist);
    orbital_info << _concat_dim("RAS4", sdist, get_dimension("RAS4"), tdist, hdist);
    orbital_info << _concat_dim("Active (total)", sdist, get_dimension("ACT"), tdist, hdist);
  }

  // Virtual spaces
  if (Parameters_->mcscf) {
    orbital_info << _concat_dim("Restricted UOCC", sdist, get_dimension("VIR"), tdist, hdist);
    orbital_info << _concat_dim("Frozen UOCC", sdist, get_dimension("FZV"), tdist, hdist);
  } else {
    orbital_info << _concat_dim("Dropped UOCC", sdist, get_dimension("DRV"), tdist, hdist);
  }

  orbital_info << seperator;

  // Print out original info
  outfile->Printf("%s", orbital_info.str().c_str());

}
}} // namespace psi::detci
