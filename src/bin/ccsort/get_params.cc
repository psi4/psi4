/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
    \ingroup CCSORT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <libciomr/libciomr.h>
#include <libpsi4util/libpsi4util.h>
#include <liboptions/liboptions.h>
#include <psi4-dec.h>
#include <psifiles.h>
#include <physconst.h>
#include "Params.h"
#include "MOInfo.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void get_params(Options & options)
{
  int errcod, tol, count, i;
  int *mu_irreps;
  std::string junk;

  params.wfn = options.get_str("WFN");
  if(params.wfn != "MP2" && params.wfn!="CCSD" &&
     params.wfn!="CCSD_T" && params.wfn!="EOM_CCSD" &&
     params.wfn!="LEOM_CCSD" && params.wfn!= "BCCD" &&
     params.wfn!="BCCD_T" && params.wfn!="SCF" &&
     params.wfn!="CIS" && params.wfn!="RPA" &&
     params.wfn!="CC2" && params.wfn!="CC3" &&
     params.wfn!="EOM_CC3" && params.wfn!="EOM_CC2" &&
     params.wfn!="CCSD_MVD" && params.wfn!="CCSD_AT") {
    outfile->Printf( "Invalid value of input keyword WFN: %s\n", params.wfn.c_str());
    throw PsiException("ccsort failure", __FILE__, __LINE__);
  }

  /* NB: SCF wfns are allowed because, at present, ccsort is needed for
     RPA calculations */

  params.semicanonical = 0;
  junk = options.get_str("REFERENCE");
  if(junk=="RHF") params.ref = 0;
  else if(junk=="ROHF" &&
          (params.wfn=="MP2" || params.wfn=="CCSD_T" || params.wfn=="CCSD_AT" ||
           params.wfn=="CC3" || params.wfn=="EOM_CC3" ||
           params.wfn=="CC2" || params.wfn=="EOM_CC2")) {
    params.ref = 2;
    params.semicanonical = 1;
  }
  else if(junk=="ROHF") params.ref = 1;
  else if(junk=="UHF") params.ref = 2;
  else {
    printf("Invalid value of input keyword REFERENCE: %s\n", junk.c_str());
    throw PsiException("ccsort failure", __FILE__, __LINE__);
  }

  // Allow user to force semicanonical
  if(options["SEMICANONICAL"].has_changed()) {
   params.semicanonical = options.get_bool("SEMICANONICAL");
   params.ref = 2;
  }

  junk = options.get_str("DERTYPE");
  if(junk=="NONE") params.dertype = 0;
  else if(junk=="FIRST") params.dertype = 1;
  else if(junk=="RESPONSE") params.dertype = 3; /* linear response */
  else {
    printf("Invalid value of input keyword DERTYPE: %s\n", junk.c_str());
    throw PsiException("ccsort failure", __FILE__, __LINE__);
  }

  if(params.wfn == "CCSD_T" && params.dertype == 1 && params.semicanonical)
    throw PsiException("ROHF-CCSD(T) gradients not yet available.", __FILE__, __LINE__);

  params.prop = options.get_str("PROPERTY");
  if(params.prop!="POLARIZABILITY" && params.prop!="ROTATION" &&
     params.prop!="ALL" && params.prop!="MAGNETIZABILITY"
     && params.prop!="ROA") {
     outfile->Printf( "Invalid choice of response property: %s\n", params.prop.c_str());
     throw PsiException("ccsort error", __FILE__, __LINE__);
  }

  params.local = options.get_bool("LOCAL");
  local.cutoff = options.get_double("LOCAL_CUTOFF");
  local.cphf_cutoff = options.get_double("LOCAL_CPHF_CUTOFF");
  local.core_cutoff = options.get_double("LOCAL_CORE_CUTOFF");

  local.method = options.get_str("LOCAL_METHOD");
  if(local.method!="AOBASIS" && local.method!="WERNER") {
    outfile->Printf( "Invalid local correlation method: %s\n", local.method.c_str());
    throw PsiException("ccsort error", __FILE__, __LINE__);
  }

  local.weakp = options.get_str("LOCAL_WEAKP");
  if(local.weakp!="MP2" && local.weakp!="NEGLECT" && local.weakp!="NONE") {
    outfile->Printf( "Invalid method for treating local pairs: %s\n", local.weakp.c_str());
    throw PsiException("ccsort failure", __FILE__, __LINE__);
  }

  local.freeze_core = options.get_str("FREEZE_CORE");

  local.pairdef = options.get_str("LOCAL_PAIRDEF");
  if(local.pairdef!="BP" && local.pairdef!="RESPONSE") {
    outfile->Printf( "Invalid keyword for strong/weak pair definition: %s\n", local.pairdef.c_str());
    throw PsiException("ccsort failure", __FILE__, __LINE__);
  }

  if(params.local && local.pairdef!="BP") {
    if(params.dertype == 3) {
      if(params.prop=="POLARIZABILITY") {
        local.domain_polar = 1;
        local.domain_mag = 0;
      }
      else if(params.prop=="MAGNETIZABILITY") {
        local.domain_polar = 0;
        local.domain_mag = 1;
      }
      else if(params.prop=="ROTATION") {
        local.domain_polar = 1;
        local.domain_mag = 1;
      }
    }
    else {
      local.domain_polar = 0;
      local.domain_mag = 0;
    }
  }
  local.domain_polar = options.get_bool("LOCAL_DOMAIN_POLAR");
  local.domain_mag = options.get_bool("LOCAL_DOMAIN_MAG");
  local.domain_sep = options.get_bool("LOCAL_DOMAIN_SEP");
  local.filter_singles = options.get_bool("LOCAL_FILTER_SINGLES");

  params.aobasis = options.get_str("AO_BASIS");

  /* Do we need MO-basis <ab|cd> integrals? */
  if(params.wfn=="MP2" || params.aobasis=="DISK" ||
     params.aobasis=="DIRECT") {
    params.make_abcd = 0;
  }
  else {
    params.make_abcd = 1;
  }

  /* If we need the MO-basis <ab|cd> integrals, do we need the fully unpacked list? */
  params.make_unpacked_abcd = 0;
  if(params.make_abcd) {
    if(params.ref != 0 || params.dertype == 1 || params.wfn=="EOM_CC2" ||
       params.wfn=="CC3" || params.wfn=="EOM_CC3" || params.wfn=="CCSD_MVD") {
      params.make_unpacked_abcd = 1;
    }
   junk = options.get_str("EOM_REFERENCE");
   if(junk=="ROHF") params.make_unpacked_abcd = 1;
  }

  /* for now, generate <ai|bc> ordering if CC gradient, ROHF-CC, CC2, or CC3 */
  if(params.dertype == 1 || params.ref == 1 || params.wfn=="CC2" ||
     params.wfn=="CC3" || params.wfn=="EOM_CC3" ||
     params.wfn=="EOM_CC2")
    params.make_aibc = 1;
  else params.make_aibc = 0;
  junk = options.get_str("EOM_REFERENCE");
  if(junk=="ROHF") params.make_aibc = 1;

  params.print_lvl = options.get_int("PRINT");
  params.keep_TEIFile = options.get_bool("KEEP_TEIFILE");
  params.keep_OEIFile = options.get_bool("KEEP_OEIFILE");

  params.tolerance = options.get_double("INTS_TOLERANCE");
  params.memory = Process::environment.get_memory();
  params.cachelev = options.get_int("CACHELEVEL");
  params.local = options.get_bool("LOCAL");

  /* grab the field frequencies from input -- a few different units are converted to E_h */
#if 0
  count = options["OMEGA"].size();
  if(count == 1) { /* assume Hartrees and only one frequency */
    params.nomega = 1;
    params.omega = init_array(1);
    params.omega[0] = options["OMEGA"][0].to_double();
  }
  else if(count >= 2) {
    params.nomega = count-1;
    params.omega = init_array(params.nomega);

    std::string units = options["OMEGA"][count-1].to_string();
    to_upper(units);

    for(i=0; i < count-1; i++) {
      params.omega[i] = options["OMEGA"][i].to_double();

      if(units=="HZ") params.omega[i] *= _h / _hartree2J;
      else if(units=="AU") 1; /* do nothing */
      else if(units=="NM") params.omega[i] = (_c*_h*1e9)/(params.omega[i]*_hartree2J);
      else if(units=="EV") params.omega[i] /= _hartree2ev;
      else {
        outfile->Printf( "\n\tError in unit for input field frequencies.  Must use one of:\n");
        outfile->Printf(   "\tau, hz, nm, or ev.\n");
        throw PsiException("ccsort error", __FILE__, __LINE__);
      }
    }
  }
  else {
    outfile->Printf( "\n\tError reading input field frequencies.  Please use the format:\n");
    outfile->Printf(   "\t  omega = (value1 value2 ... units)\n");
    outfile->Printf(   "\twhere units = hartrees, hz, nm, or ev.\n");
    throw PsiException("Failure in ccsort.", __FILE__, __LINE__);
  }
  mu_irreps = init_int_array(3);
  moinfo.irrep_x = options["MU_IRREPS"][0].to_integer();
  moinfo.irrep_y = options["MU_IRREPS"][1].to_integer();
  moinfo.irrep_z = options["MU_IRREPS"][2].to_integer();
#endif

  outfile->Printf( "\n\tInput parameters:\n");
  outfile->Printf( "\t-----------------\n");
  outfile->Printf( "\tWave function   =\t%s\n", params.wfn.c_str());
  if(params.semicanonical) {
    outfile->Printf( "\tReference wfn   =\tROHF changed to UHF for Semicanonical Orbitals\n");
  }
  else {
    outfile->Printf( "\tReference wfn   =\t%s\n",
            (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  }
  if(params.dertype == 0) outfile->Printf( "\tDerivative      =\tNone\n");
  else if(params.dertype == 1) outfile->Printf( "\tDerivative      =\tFirst\n");
  else if(params.dertype == 3) outfile->Printf( "\tDerivative      =\tResponse\n");
  outfile->Printf( "\tMemory (Mbytes) =\t%.1f\n", params.memory/1e6);
  outfile->Printf( "\tAO Basis        =\t%s\n", params.aobasis.c_str());
  outfile->Printf( "\tMake (ab|cd)    =\t%s\n",
          (params.make_abcd == 1) ? "True" : "False");
  outfile->Printf( "\tMake unpacked (ab|cd) =\t%s\n",
          (params.make_unpacked_abcd == 1) ? "True" : "False");
  outfile->Printf( "\tCache Level     =\t%d\n", params.cachelev);
  outfile->Printf( "\tCache Type      =\t%s\n", "LRU");
  outfile->Printf( "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    outfile->Printf( "\tLocal Cutoff       = %3.1e\n", local.cutoff);
    outfile->Printf( "\tLocal Core Cutoff  = %3.1e\n", local.core_cutoff);
    outfile->Printf( "\tLocal Method      =    %s\n", local.method.c_str());
    outfile->Printf( "\tWeak pairs        =    %s\n", local.weakp.c_str());
    outfile->Printf( "\tFilter singles    =    %s\n", local.filter_singles ? "Yes" : "No");
    outfile->Printf( "\tLocal pairs       =    %s\n", local.pairdef.c_str());
    outfile->Printf( "\tLocal CPHF cutoff =  %3.1e\n", local.cphf_cutoff);
  }
  outfile->Printf( "\n");
  
}

}} // namespace psi::ccsort