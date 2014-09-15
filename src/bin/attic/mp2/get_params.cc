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
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void get_params(Options &options)
{
  int errcod;
  std::string cachetype = "NULL";
  std::string ref;
  
  params.wfn = options.get_str("WFN");
  ref = options.get_str("REFERENCE");

  /* Default reference is RHF */
  params.ref = 0;
  params.semicanonical = 0;
  if(ref == "RHF") params.ref = 0;
  else if(ref == "ROHF" && params.wfn == "MP2") {
    params.ref = 2;
    params.semicanonical = 1;
  }
  else if(ref == "ROHF") params.ref = 1;
  else if(ref == "UHF") params.ref = 2;
  else {
    throw PsiException("Invalid Reference", __FILE__, __LINE__);
  }
  
  /* Default Jobtype */
  params.jobtype = options.get_str("JOBTYPE");

  /* Default Dertype */
  params.dertype = options.get_str("DERTYPE");

  params.gradient = options.get_str("DERTYPE") == "FIRST";
  params.relax_opdm = options.get_bool("OPDM_RELAX");
  params.opdm = options.get_bool("OPDM");
//  if(params.jobtype == "SP") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 0;
//  }
//  else if(params.jobtype == "OEPROP" && params.dertype == "NONE") {
//    params.opdm = 1;
//    params.relax_opdm = 0;
//    params.gradient = 0;
//  }
//  else if(params.jobtype == "OEPROP" && params.dertype == "FIRST") {
//    params.opdm = 1;
//    params.relax_opdm = 1;
//    params.gradient = 0;
//  }
//  else if(params.jobtype == "OPT" && params.dertype == "NONE") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 0;
//  }
//  else if(params.jobtype == "OPT" && params.dertype == "FIRST") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 1;
//  }
//  else if(params.jobtype == "OPT_FC" && params.dertype == "FIRST") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 1;
//  }
//  else if(params.jobtype == "SYMM_FC" && params.dertype == "FIRST") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 1;
//  }
//  else if(params.jobtype == "FREQ" && params.dertype == "NONE") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 0;
//  }
//  else if(params.jobtype == "FREQ" && params.dertype == "FIRST") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 1;
//  }
//  else {
//    throw PsiException("Invalid combination of JOBTYPE and DERTYPE", __FILE__, __LINE__);
//  }

  if((params.relax_opdm || params.gradient) && 
     (mo.nfzdocc != 0 || mo.nfzvirt != 0)) {
    throw PsiException("The Z-vector equations DO NOT work with frozen orbitals ... yet", __FILE__, __LINE__);
  }

  params.print = options.get_int("PRINT");
  
  params.cachelev = options.get_int("CACHELEVEL");
  
  cachetype = options.get_str("CACHETYPE");
  if(cachetype != "NULL") {
    if(cachetype == "LOW")
      params.cachetype = 1;
    else if(cachetype == "LRU")
      params.cachetype = 0;
    else {
      outfile->Printf( "Invalide CACHETYPE = %s\n",cachetype.c_str());
      abort();
    }
  }
  
  /* get parameters related to SCS-MP2 or SCS-N-MP2 */
  /* see papers by S. Grimme or J. Platz */
  params.scs = options.get_int("SCS_N");
  if (params.scs == 1) {
    params.scs_scale_os = 0.0;
    params.scs_scale_ss = 1.76;
  }
  if(options.get_int("SCS") == 1) 
    params.scs = 1;
  if (params.scs == 1) { 
    params.scs_scale_os = options.get_double("MP2_OS_SCALE");
    params.scs_scale_ss = options.get_double("MP2_SS_SCALE");
  }

  params.memory = Process::environment.get_memory();
 
  outfile->Printf( "\n");
  outfile->Printf( "\tInput parameters:\n");
  outfile->Printf( "\t-----------------\n");
  outfile->Printf( "\tWave function \t=\t%s\n", params.wfn.c_str());
  if(params.semicanonical) {
  outfile->Printf( "\tReference WFN \t=\tROHF changed to UHF for Semicanonical Orbitals\n");
  }
  else {
  outfile->Printf( "\tReference WFN \t=\t%s\n", (params.ref==0)?"RHF":((params.ref==1)?"ROHF":"UHF"));
  } 
  outfile->Printf( "\tDerivative    \t=\t%s\n", params.dertype.c_str());
  outfile->Printf( "\tCache Level   \t=\t%d\n", params.cachelev);
  outfile->Printf( "\tCache Type    \t=\t%s\n", params.cachetype ? "LOW":"LRU");
  outfile->Printf( "\tMemory (MB)   \t=\t%.1f\n",params.memory/1e6);
  outfile->Printf( "\tPrint Level   \t=\t%d\n", params.print);
  outfile->Printf( "\tOPDM          \t=\t%s\n", params.opdm ? "YES":"NO");
  outfile->Printf( "\tSCS           \t=\t%s\n", (params.scs == 1) ? "True" : "False");
  outfile->Printf( "\tMP2_OS_SCALE  \t=\t%.6f\n",params.scs_scale_os);
  outfile->Printf( "\tMP2_SS_SCALE  \t=\t%.6f\n",params.scs_scale_ss);

  if (params.scs && params.dertype != "NONE") {
    outfile->Printf("\nWarning: SCS-MP2 computation requested but\n");
    outfile->Printf("derivative will be evaluated for standard MP2 energy.\n");
  }

}

}} /* End namespaces */
