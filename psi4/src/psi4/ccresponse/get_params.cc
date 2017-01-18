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
#include <cstdio>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/physconst.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/multipolesymmetry.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void get_params(std::shared_ptr<Wavefunction> wfn, Options &options)
{
  int i, errcod, ref, count, iconv, *tmpi;
  std::string units;
  std::string junk;

  params.wfn = options.get_str("WFN");
  if(params.wfn != "CCSD" && params.wfn != "CC2") {
    throw PsiException("Invalid value of input keyword WFN",__FILE__,__LINE__);
  }

  params.print = options.get_int("PRINT");

  params.memory = Process::environment.get_memory();

  params.cachelev = options.get_int("CACHELEVEL");
  params.cachelev = 0;

  junk = options.get_str("REFERENCE");
  /* if no reference is given, assume rhf */
  if(junk == "RHF") ref = 0;
  else if(junk == "ROHF") ref = 1;
  else if(junk == "UHF") ref = 2;
  else {
    throw PsiException("Invalid value of input keyword REFERENCE",__FILE__,__LINE__);
  }

  /* Make sure the value of ref matches that from CC_INFO */
  if(params.ref != ref) {
    outfile->Printf( "Value of REFERENCE from input.dat (%1d) and CC_INFO (%1d) do not match!\n",
            ref, params.ref);
    outfile->Printf( "Is this what you want to do?\n");
    params.ref = ref;
  }

  junk = options.get_str("DERTYPE");
  if(junk == "NONE") params.dertype = 0;
  else if(junk == "FIRST") params.dertype = 1;
  else if(junk == "RESPONSE") params.dertype = 3; /* linear response */
  else {
    throw PsiException("Invalid value of input keyword DERTYPE",__FILE__,__LINE__);
  }

  params.gauge = options.get_str("GAUGE");
  if(params.gauge != "LENGTH" && params.gauge != "VELOCITY" && params.gauge != "BOTH") {
    throw PsiException("Invalid choice of gauge",__FILE__,__LINE__);
  }

  // grab the field freqs from input -- a few units are converted to E_h
  count = options["OMEGA"].size();
  if(count == 0) { // Assume 0.0 E_h for field energy
    params.nomega = 1;
    params.omega = init_array(1);
    params.omega[0] = 0.0;
  }
  else if(count == 1) { // Assume E_h for field energy and read value
    params.nomega = 1;
    params.omega = init_array(1);
    params.omega[0] = options["OMEGA"][0].to_double();
  }
  else if(count >= 2) {
    params.nomega = count-1;
    params.omega = init_array(params.nomega);
    units = options["OMEGA"][count-1].to_string();
    for(i=0; i < (count-1); i++) {
      params.omega[i] = options["OMEGA"][i].to_double();
      if(units == "HZ" || units == "Hz" || units == "hz")
        params.omega[i] *= pc_h / pc_hartree2J;
      else if(units == "AU" || units == "Au" || units == "au") continue; // do nothing
      else if(units == "NM" || units == "nm")
        params.omega[i] = (pc_c*pc_h*1e9)/(params.omega[i]*pc_hartree2J);
      else if(units == "EV" || units == "ev" || units == "eV")
        params.omega[i] /= pc_hartree2ev;
      else
        throw PsiException("Error in unit for input field frequencies, should be au, Hz, nm, or eV", __FILE__,__LINE__);
    }
  }

  std::shared_ptr<Molecule> mol = wfn->molecule();
  std::shared_ptr<IntegralFactory> intfact = wfn->integral();
  std::shared_ptr<MatrixFactory> matfact = wfn->matrix_factory();

  OperatorSymmetry dipsym(1, mol, intfact, matfact);
  moinfo.mu_irreps = init_int_array(3);
  moinfo.mu_irreps[0] = dipsym.component_symmetry(0);
  moinfo.mu_irreps[1] = dipsym.component_symmetry(1);
  moinfo.mu_irreps[2] = dipsym.component_symmetry(2);

  /* compute the irreps of the angular momentum operator while we're here */
  moinfo.l_irreps = init_int_array(3);
  for(i=0; i < 3; i++)
    moinfo.l_irreps[i] = moinfo.mu_irreps[(int) (i+1)%3] ^ moinfo.mu_irreps[(int) (i+2)%3];

  params.maxiter = options.get_int("MAXITER");
  params.convergence = options.get_double("R_CONVERGENCE");
  params.diis = options.get_bool("DIIS");

  params.prop = options.get_str("PROPERTY");
  if(params.prop != "POLARIZABILITY" && params.prop != "ROTATION"
     && params.prop != "ROA" && params.prop != "ROA_TENSOR"
     && params.prop != "ALL") {
    throw PsiException("Invalid choice of resp. property",__FILE__,__LINE__);
  }

  params.abcd = options.get_str("ABCD");
  if(params.abcd != "NEW" && params.abcd != "OLD") {
    throw PsiException("Invalid ABCD algorith",__FILE__,__LINE__);
  }


  params.restart = options.get_bool("RESTART");

  params.local = options.get_bool("LOCAL");
  local.cutoff = options.get_double("LOCAL_CUTOFF");

  local.method = options.get_str("LOCAL_METHOD");
  if(local.method != "AOBASIS" && local.method != "WERNER") {
      throw PsiException("Invalid local correlation method",__FILE__,__LINE__);
  }

  local.weakp = options.get_str("LOCAL_WEAKP");
  if(local.weakp != "MP2" && local.weakp != "NEGLECT" && local.weakp != "NONE") {
    throw PsiException("Invalid method for treating local pairs",__FILE__,__LINE__);
  }

  if(params.dertype == 3)
    local.filter_singles = 0;
  else
    local.filter_singles = 1;
  local.filter_singles = options.get_bool("LOCAL_FILTER_SINGLES");

  local.cphf_cutoff = options.get_double("LOCAL_CPHF_CUTOFF");
  local.freeze_core = options.get_str("FREEZE_CORE");

  if(options["LOCAL_PAIRDEF"].has_changed()){
    local.pairdef = options.get_str("LOCAL_PAIRDEF");
    if(local.pairdef != "BP" && local.pairdef != "RESPONSE") {
      throw PsiException("Invalid keyword for strong/weak pair definating", __FILE__,__LINE__);
    }
  }
  else if(params.local && params.dertype == 3)
    local.pairdef = strdup("RESPONSE");
  else if(params.local)
    local.pairdef = strdup("BP");

  params.analyze = options.get_bool("ANALYZE");
  params.num_amps = options.get_int("NUM_AMPS_PRINT");
  params.sekino = options.get_bool("SEKINO");
  params.linear = options.get_bool("LINEAR");

  outfile->Printf( "\n\tInput parameters:\n");
  outfile->Printf( "\t-----------------\n");
  if(params.prop == "ALL")
    outfile->Printf( "\tProperty         =    POLARIZABILITY + ROTATION\n");
  else
    outfile->Printf( "\tProperty         =    %s\n", params.prop.c_str());
  outfile->Printf( "\tReference wfn    =    %s\n",
          (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  outfile->Printf( "\tMemory (Mbytes)  =    %5.1f\n",params.memory/1e6);
  outfile->Printf( "\tCache Level      =    %d\n", params.cachelev);
  outfile->Printf( "\tPrint Level      =    %d\n",  params.print);
  outfile->Printf( "\tMaxiter          =    %d\n",  params.maxiter);
  outfile->Printf( "\tConvergence      =    %3.1e\n", params.convergence);
  outfile->Printf( "\tRestart          =    %s\n", params.restart ? "Allowed" : "Not Allowed");
  outfile->Printf( "\tDIIS             =    %s\n", params.diis ? "Yes" : "No");
  outfile->Printf( "\tModel III        =    %s\n", params.sekino ? "Yes" : "No");
  outfile->Printf( "\tLinear Model     =    %s\n", params.linear ? "Yes" : "No");
  outfile->Printf( "\tABCD             =    %s\n", params.abcd.c_str());
  outfile->Printf( "\tIrrep X          =    %s\n", moinfo.labels[moinfo.mu_irreps[0]]);
  outfile->Printf( "\tIrrep Y          =    %s\n", moinfo.labels[moinfo.mu_irreps[1]]);
  outfile->Printf( "\tIrrep Z          =    %s\n", moinfo.labels[moinfo.mu_irreps[2]]);
  outfile->Printf( "\tIrrep RX         =    %s\n", moinfo.labels[moinfo.l_irreps[0]]);
  outfile->Printf( "\tIrrep RY         =    %s\n", moinfo.labels[moinfo.l_irreps[1]]);
  outfile->Printf( "\tIrrep RZ         =    %s\n", moinfo.labels[moinfo.l_irreps[2]]);
  outfile->Printf( "\tGauge            =    %s\n", params.gauge.c_str());
  for(i=0; i < params.nomega; i++) {
    if(params.omega[i] == 0.0)
      outfile->Printf( "\tApplied field %2d =  0.000\n", i);
    else
      outfile->Printf( "\tApplied field %2d =    %5.3f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", i, params.omega[i],
              (pc_c*pc_h*1e9)/(pc_hartree2J*params.omega[i]), pc_hartree2ev*params.omega[i],
              pc_hartree2wavenumbers*params.omega[i]);
  }
  outfile->Printf( "\tAnalyze X2 Amps  =    %s\n", params.analyze ? "Yes" : "No");
  outfile->Printf( "\tLocal CC         =    %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    outfile->Printf( "\tLocal Cutoff      =    %3.1e\n", local.cutoff);
    outfile->Printf( "\tLocal Method      =    %s\n", local.method.c_str());
    outfile->Printf( "\tWeak pairs        =    %s\n", local.weakp.c_str());
    outfile->Printf( "\tFilter singles    =    %s\n", local.filter_singles ? "Yes" : "No");
    outfile->Printf( "\tLocal pairs       =    %s\n", local.pairdef.c_str());
    outfile->Printf( "\tLocal CPHF cutoff =    %3.1e\n", local.cphf_cutoff);
  }
  outfile->Printf( "\n");
}


}} // namespace psi::ccresponse
