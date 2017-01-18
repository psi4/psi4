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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "Params.h"
#include "Local.h"
#include "MOInfo.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::get_params(Options &options)
{
  int errcod, iconv, forceit;
  std::string cachetype = "";
  std::string junk;

  params_.newtrips = options.get_bool("NEW_TRIPLES");

  params_.wfn = options.get_str("WFN");

  if(params_.wfn == "NONE")
     throw PsiException("Invalid value of input keyword WFN", __FILE__, __LINE__);

  if(params_.wfn == "BCCD" || params_.wfn == "BCCD_T")
    params_.brueckner = 1;
  else params_.brueckner = 0;

  params_.df = options.get_str("CC_TYPE") == "DF";

  params_.semicanonical = 0;
  junk = options.get_str("REFERENCE");
  /* if no reference is given, assume rhf */
  if(junk == "RHF") params_.ref = 0;
  else if(junk == "ROHF" &&
    (params_.wfn == "MP2" || params_.wfn == "CCSD_T" || params_.wfn == "CCSD_AT" ||
    params_.wfn == "CC3" || params_.wfn == "EOM_CC3" ||
    params_.wfn == "CC2" || params_.wfn == "EOM_CC2" ||
    params_.wfn == "BCCD" || params_.wfn == "BCCD_T")) {
    params_.ref = 2;
    params_.semicanonical = 1;
  }
  else if(junk == "ROHF") params_.ref = 1;
  else if(junk == "UHF" ) params_.ref = 2;
  else
   throw PsiException("Invalid value of input keyword REFERENCE", __FILE__, __LINE__);

  // Allow user to force semicanonical
  if(options["SEMICANONICAL"].has_changed()) {
   params_.semicanonical = options.get_bool("SEMICANONICAL");
   params_.ref = 2;
  }

  params_.analyze = options.get_bool("ANALYZE");

  params_.dertype = 0;
  junk = options.get_str("DERTYPE");
  if(junk == "NONE") params_.dertype = 0;
  else if(junk == "FIRST") params_.dertype = 1;
  else if(junk == "RESPONSE") params_.dertype = 3; /* linear response */
  else
   throw PsiException("Invalid value of input keyword DERTYPE", __FILE__, __LINE__);

  params_.print = options.get_int("PRINT");
  params_.maxiter = options.get_int("MAXITER");
  params_.convergence = options.get_double("R_CONVERGENCE");
  params_.e_convergence = options.get_double("E_CONVERGENCE");
  params_.restart = options.get_bool("RESTART");

  params_.memory = Process::environment.get_memory();

  params_.aobasis = options.get_str("AO_BASIS");
  params_.cachelev = options.get_int("CACHELEVEL");

  params_.cachetype = 1;
  cachetype = options.get_str("CACHETYPE");
  if(cachetype == "LOW") params_.cachetype = 1;
  else if(cachetype == "LRU") params_.cachetype = 0;
  else
    throw PsiException("Error in input: invalid CACHETYPE", __FILE__, __LINE__);


 if(params_.ref == 2) /* No LOW cacheing yet for UHF references */
    params_.cachetype = 0;

  params_.nthreads = Process::environment.get_n_threads();
  if (options["CC_NUM_THREADS"].has_changed()){
     params_.nthreads = options.get_int("CC_NUM_THREADS");
  }
  params_.diis = options.get_bool("DIIS");
  params_.t2_coupled = options.get_bool("T2_COUPLED");
  params_.prop = options.get_str("PROPERTY");
  params_.abcd = options.get_str("ABCD");
  params_.local = options.get_bool("LOCAL");
  local_.cutoff = options.get_double("LOCAL_CUTOFF");
  local_.method = options.get_str("LOCAL_METHOD");
  local_.weakp = options.get_str("LOCAL_WEAKP");

  //local.filter_singles = options.get_bool("LOCAL_FILTER_SINGLES");
  //if(params.dertype == 3) local.filter_singles = 0;

  local_.cphf_cutoff = options.get_double("LOCAL_CPHF_CUTOFF");
  std::string freeze_docc = options.get_str("FREEZE_CORE");
  local_.freeze_core = (freeze_docc != "FALSE");

  local_.pairdef = options.get_str("LOCAL_PAIRDEF");
  if(params_.local && params_.dertype == 3)
    local_.pairdef = "RESPONSE";
  else if(params_.local)
    local_.pairdef = "BP";

  params_.num_amps = options.get_int("NUM_AMPS_PRINT");
  params_.bconv = options.get_double("BRUECKNER_ORBS_R_CONVERGENCE");

  // Tying orbital convergence to the desired e_conv,
  //   particularly important for sane numerical frequencies by energy
  if (options["BRUECKNER_ORBS_R_CONVERGENCE"].has_changed())
      params_.bconv = options.get_double("BRUECKNER_ORBS_R_CONVERGENCE");
  else
      params_.bconv = 100.0 * params_.e_convergence;

  params_.print_mp2_amps = options.get_bool("MP2_AMPS_PRINT");
  params_.print_pair_energies = options.get_bool("PAIR_ENERGIES_PRINT");
  params_.spinadapt_energies = options.get_bool("SPINADAPT_ENERGIES");
  params_.t3_Ws_incore = options.get_bool("T3_WS_INCORE");

  /* get parameters related to SCS-MP2 or SCS-N-MP2 */
  /* see papers by S. Grimme or J. Platz */
  params_.scsn = options.get_bool("SCSN_MP2");
  params_.scs = options.get_bool("SCS_MP2");
  params_.scscc = options.get_bool("SCS_CCSD");
  params_.scsmp2_scale_os = options.get_double("MP2_OS_SCALE");
  params_.scsmp2_scale_ss = options.get_double("MP2_SS_SCALE");
  /* see paper by T. Takatani*/
  params_.scscc_scale_os = options.get_double("CC_OS_SCALE");
  params_.scscc_scale_ss = options.get_double("CC_SS_SCALE");

  if (options["MP2_OS_SCALE"].has_changed() || options["MP2_SS_SCALE"].has_changed()) {
    params_.scs = 1;
    }

  if (options["CC_OS_SCALE"].has_changed() || options["CC_SS_SCALE"].has_changed()) {
    params_.scscc = 1;
    }


  outfile->Printf( "\n    Input parameters:\n");
  outfile->Printf( "    -----------------\n");
  outfile->Printf( "    Wave function   =     %s\n", params_.wfn.c_str());

  if(params_.semicanonical) {
    outfile->Printf( "    Reference wfn   =     ROHF changed to UHF for Semicanonical Orbitals\n");
  }
  else {
    outfile->Printf( "    Reference wfn   =     %s\n",
        (params_.ref == 0) ? "RHF" : ((params_.ref == 1) ? "ROHF" : "UHF"));
  }
  outfile->Printf("    Brueckner       =     %s\n", params_.brueckner ? "Yes" : "No");
  if(params_.brueckner)
    outfile->Printf( "    Brueckner conv. =     %3.1e\n", params_.bconv);
  outfile->Printf( "    Memory (Mbytes) =     %5.1f\n",params_.memory/1e6);
  outfile->Printf( "    Maxiter         =   %4d\n", params_.maxiter);
  outfile->Printf( "    R_Convergence   =     %3.1e\n", params_.convergence);
  outfile->Printf( "    E_Convergence   =     %3.1e\n", params_.e_convergence);
  outfile->Printf( "    Restart         =     %s\n",
      params_.restart ? "Yes" : "No");
  outfile->Printf( "    DIIS            =     %s\n", params_.diis ? "Yes" : "No");
  outfile->Printf( "    AO Basis        =     %s\n", params_.aobasis.c_str());
  outfile->Printf( "    ABCD            =     %s\n", params_.abcd.c_str());
  outfile->Printf( "    Cache Level     =     %1d\n", params_.cachelev);
  outfile->Printf( "    Cache Type      =    %4s\n",
      params_.cachetype ? "LOW" : "LRU");
  outfile->Printf( "    Print Level     =     %1d\n",  params_.print);
  outfile->Printf( "    Num. of threads =     %d\n",  params_.nthreads);
  outfile->Printf( "    # Amps to Print =     %1d\n",  params_.num_amps);
  outfile->Printf( "    Print MP2 Amps? =     %s\n",  params_.print_mp2_amps ?
      "Yes" : "No" );
  outfile->Printf( "    Analyze T2 Amps =     %s\n",  params_.analyze ? "Yes" : "No" );
  outfile->Printf( "    Print Pair Ener =     %s\n",  params_.print_pair_energies ? "Yes" : "No" );

  if (params_.print_pair_energies)
    outfile->Printf( "    Spinadapt Ener. =     %s\n",  params_.spinadapt_energies ? "Yes" : "No" );
  outfile->Printf( "    Local CC        =     %s\n", params_.local ? "Yes" : "No");

  if ( params_.wfn == "CC3" || params_.wfn == "EOM_CC3")
    outfile->Printf( "    T3 Ws incore    =     %s\n", params_.t3_Ws_incore ? "Yes" : "No");

  if(params_.local) {
    outfile->Printf( "    Local Cutoff       =     %3.1e\n", local_.cutoff);
    outfile->Printf( "    Local Method      =     %s\n", local_.method.c_str());
    outfile->Printf( "    Weak pairs        =     %s\n", local_.weakp.c_str());
    outfile->Printf( "    Filter singles    =     %s\n", local_.filter_singles ? "Yes" : "No");
    outfile->Printf( "    Local pairs       =     %s\n", local_.pairdef.c_str());
    outfile->Printf( "    Local CPHF cutoff =     %3.1e\n", local_.cphf_cutoff);
  }
  outfile->Printf( "    SCS-MP2         =     %s\n", (params_.scs == 1) ? "True" : "False");
  outfile->Printf( "    SCSN-MP2        =     %s\n", (params_.scsn == 1) ? "True" : "False");
  outfile->Printf( "    SCS-CCSD        =     %s\n", (params_.scscc == 1) ? "True" : "False");
  if (params_.scs) {
    outfile->Printf( "    SCS_MP2_OS_SCALE =     %.2f\n",params_.scsmp2_scale_os);
    outfile->Printf( "    SCS_MP2_SS_SCALE =     %.2f\n",params_.scsmp2_scale_ss);
  }
  if (params_.scsn) {
    outfile->Printf( "    SCSN_MP2_OS_SCALE =     %.2f\n",0.0);
    outfile->Printf( "    SCSN_MP2_SS_SCALE =     %.2f\n",1.76);
  }
  if (params_.scscc) {
    outfile->Printf( "    CC_OS_SCALE     =     %.2f\n",params_.scscc_scale_os);
    outfile->Printf( "    CC_SS_SCALE     =     %.2f\n",params_.scscc_scale_ss);
  }

  outfile->Printf( "\n");

}
}} // namespace psi::ccenergy
