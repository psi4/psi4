/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <psi4-dec.h>
#include <psifiles.h>
#include <physconst.h>
#include <libmints/wavefunction.h>
#include <libmints/molecule.h>
#include <libmints/integral.h>
#include <libmints/multipolesymmetry.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void get_params(Options &options)
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

  params.cachelev = options.get_int("CACHELEV");
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
    fprintf(outfile, "Value of REFERENCE from input.dat (%1d) and CC_INFO (%1d) do not match!\n",
            ref, params.ref);
    fprintf(outfile, "Is this what you want to do?\n");
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

  /* grab the field frequencies from input -- a few different units are converted to E_h */
  units = options.get_str("OMEGA_UNITS");
  count = options["OMEGA"].size();
  if(count == 0) {
    params.nomega = 1;
    params.omega = init_array(1);
    params.omega[0] = 0.0;
  }
  else {
    params.nomega = count;
    params.omega = init_array(params.nomega);

    for(i=0; i < count; i++) {
      params.omega[i] = options["OMEGA"][i].to_double();

      if(units == "HZ") params.omega[i] *= _h / _hartree2J;
      else if(units == "AU") 1; /* do nothing */
      else if(units == "NM") params.omega[i] = (_c*_h*1e9)/(params.omega[i]*_hartree2J);
      else if(units == "EV") params.omega[i] /= _hartree2ev;
      else
        throw PsiException("Error in unit for input field frequencies, should be au, Hz, nm, or eV", __FILE__,__LINE__);
    }
  }

  boost::shared_ptr<Wavefunction> wfn = Process::environment.reference_wavefunction();
  boost::shared_ptr<Molecule> mol = wfn->molecule();
  boost::shared_ptr<IntegralFactory> fact = wfn->integral();

  OperatorSymmetry dipsym(1, mol, fact);
  moinfo.mu_irreps = init_int_array(3);
  moinfo.mu_irreps[0] = dipsym.component_symmetry(0);
  moinfo.mu_irreps[1] = dipsym.component_symmetry(1);
  moinfo.mu_irreps[2] = dipsym.component_symmetry(2);

  /* compute the irreps of the angular momentum operator while we're here */
  moinfo.l_irreps = init_int_array(3);
  for(i=0; i < 3; i++)
    moinfo.l_irreps[i] = moinfo.mu_irreps[(int) (i+1)%3] ^ moinfo.mu_irreps[(int) (i+2)%3];

  params.maxiter = options.get_int("MAXITER");
  params.convergence = 1.0*pow(10.0, (double) -options.get_int("CONVERGENCE"));
  params.diis = options.get_bool("DIIS");

  params.prop = options.get_str("PROPERTY");
  if(params.prop != "POLARIZABILITY" && params.prop != "ROTATION"
     && params.prop != "ROA" && params.prop != "ALL") {
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
  local.filter_singles = options.get_bool("LOCAL_FILER_SINGLES");

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
  params.num_amps = options.get_int("NUM_AMPS");
  params.sekino = options.get_bool("SEKINO");
  params.linear = options.get_bool("LINEAR");

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  if(params.prop == "ALL")
    fprintf(outfile, "\tProperty         =    POLARIZABILITY + ROTATION\n");
  else
    fprintf(outfile, "\tProperty         =    %s\n", params.prop.c_str());
  fprintf(outfile, "\tReference wfn    =    %5s\n",
          (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  fprintf(outfile, "\tMemory (Mbytes)  =  %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tCache Level      =    %1d\n", params.cachelev);
  fprintf(outfile, "\tPrint Level      =    %1d\n",  params.print);
  fprintf(outfile, "\tMaxiter          =    %3d\n",  params.maxiter);
  fprintf(outfile, "\tConvergence      = %3.1e\n", params.convergence);
  fprintf(outfile, "\tRestart          =     %s\n", params.restart ? "Allowed" : "Not Allowed");
  fprintf(outfile, "\tDIIS             =     %s\n", params.diis ? "Yes" : "No");
  fprintf(outfile, "\tModel III        =     %s\n", params.sekino ? "Yes" : "No");
  fprintf(outfile, "\tLinear Model     =     %s\n", params.linear ? "Yes" : "No");
  fprintf(outfile, "\tABCD             =     %s\n", params.abcd.c_str());
  fprintf(outfile, "\tIrrep X          =    %3s\n", moinfo.labels[moinfo.mu_irreps[0]]);
  fprintf(outfile, "\tIrrep Y          =    %3s\n", moinfo.labels[moinfo.mu_irreps[1]]);
  fprintf(outfile, "\tIrrep Z          =    %3s\n", moinfo.labels[moinfo.mu_irreps[2]]);
  fprintf(outfile, "\tIrrep RX         =    %3s\n", moinfo.labels[moinfo.l_irreps[0]]);
  fprintf(outfile, "\tIrrep RY         =    %3s\n", moinfo.labels[moinfo.l_irreps[1]]);
  fprintf(outfile, "\tIrrep RZ         =    %3s\n", moinfo.labels[moinfo.l_irreps[2]]);
  fprintf(outfile, "\tGauge            =    %s\n", params.gauge.c_str());
  for(i=0; i < params.nomega; i++) {
    if(params.omega[i] == 0.0)
      fprintf(outfile, "\tApplied field %2d =  0.000\n", i);
    else
      fprintf(outfile, "\tApplied field %2d =    %5.3f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", i, params.omega[i],
              (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i],
              _hartree2wavenumbers*params.omega[i]);
  }
  fprintf(outfile, "\tAnalyze X2 Amps  =    %s\n", params.analyze ? "Yes" : "No");
  fprintf(outfile, "\tLocal CC         =    %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff      = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method      =    %s\n", local.method.c_str());
    fprintf(outfile, "\tWeak pairs        =    %s\n", local.weakp.c_str());
    fprintf(outfile, "\tFilter singles    =    %s\n", local.filter_singles ? "Yes" : "No");
    fprintf(outfile, "\tLocal pairs       =    %s\n", local.pairdef.c_str());
    fprintf(outfile, "\tLocal CPHF cutoff =  %3.1e\n", local.cphf_cutoff);
  }
  fprintf(outfile, "\n");
}


}} // namespace psi::ccresponse
