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
** PARAMS.C
**
** Gets/prints user specified input
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
*/

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include "globals.h"
#include "psi4-dec.h"

namespace psi { namespace detcas {

/*
** get_parameters(): Function gets the program running parameters such
**   as convergence.  These are stored in the Parameters data structure.
*/
void get_parameters(Options &options)
{
  int i, errcod;
  char line1[133];
   
  Params.wfn = options.get_str("WFN");
  Params.dertype = options.get_str("DERTYPE");
  if (Params.dertype == "NONE") {
    Params.rms_grad_convergence = 1e-4;
    Params.energy_convergence = 1e-7;
  }
  else {
    Params.rms_grad_convergence = 1e-7;
    Params.energy_convergence = 1e-10;
  }

  if (options["R_CONVERGENCE"].has_changed()) {
    Params.rms_grad_convergence = options.get_double("R_CONVERGENCE");
  }
  if (options["E_CONVERGENCE"].has_changed()) {
    Params.energy_convergence = options.get_double("E_CONVERGENCE");
  }


  Params.filter_ints = 0;  /* assume we need all for MCSCF */
  Params.oei_file = PSIF_OEI;  /* contains frozen core operator */
  Params.tei_file = PSIF_MO_TEI;
  Params.opdm_file = PSIF_MO_OPDM;
  Params.tpdm_file = PSIF_MO_TPDM;
  Params.lag_file = PSIF_MO_LAG;

  Params.ignore_fz = true;     /* ignore frozen orbitals for ind pairs? */
  
  if ((Params.wfn == "CASSCF") || (Params.wfn == "DETCAS"))
    Params.ignore_ras_ras = true;   /* ignore RAS/RAS independent pairs? */
  else
    Params.ignore_ras_ras = false;

  Params.print_lvl = options.get_int("PRINT");
  Params.print_mos = options.get_bool("PRINT_MOS");

  Params.oei_erase = options.get_bool("OEI_ERASE");
  Params.tei_erase = options.get_bool("TEI_ERASE");
  Params.ignore_fz = options.get_bool("IGNORE_FZ");
  Params.scale_grad = true; /* scale orb grad by inverse Hessian? */

  Params.diis_start = options.get_int("DIIS_START");
  Params.diis_freq = options.get_int("DIIS_FREQ");
  Params.diis_min_vecs = options.get_int("DIIS_MIN_VECS");
  Params.diis_max_vecs = options.get_int("DIIS_MAX_VECS");

  Params.scale_step = options.get_double("SCALE_STEP");
  Params.use_fzc_h = options.get_bool("USE_FZC_H");
  Params.invert_hessian = options.get_bool("INVERT_HESSIAN");
  Params.hessian = options.get_str("HESSIAN");
 
  Params.level_shift = options.get_bool("DO_LEVEL_SHIFT");
  Params.shift = options.get_double("SHIFT");
  Params.determ_min = options.get_double("DETERM_MIN");
  Params.step_max = options.get_double("STEP_MAX");
  Params.use_thetas = options.get_bool("USE_THETAS");
  Params.force_step = options.get_bool("FORCE_STEP");
  Params.force_pair = options.get_int("FORCE_PAIR");
  Params.force_value = options.get_double("FORCE_VALUE");
  Params.scale_act_act = options.get_double("SCALE_ACT_ACT");
  Params.bfgs = options.get_bool("BFGS");
  Params.ds_hessian = options.get_bool("DS_HESSIAN");

}


/*
** print_parameters(): Function prints the program's running parameters
**   found in the Parameters structure.
*/
void print_parameters(void)
{
  outfile->Printf("\n") ;
  outfile->Printf("PARAMETERS: \n") ;
  outfile->Printf("   PRINT         =   %6d      PRINT_MOS     =   %6s\n", 
      Params.print_lvl, Params.print_mos ? "yes" : "no");
  outfile->Printf("   R_CONVERGENCE   =   %6.2e  E CONVERG     =   %6.2e\n",
      Params.rms_grad_convergence, Params.energy_convergence);
  outfile->Printf("   IGNORE_RAS_RAS=   %6s      IGNORE_FZ     =   %6s\n", 
      Params.ignore_ras_ras ? "yes" : "no", Params.ignore_fz ? "yes" : "no") ;
  outfile->Printf("   OEI FILE      =   %6d      OEI ERASE     =   %6s\n", 
      Params.oei_file, Params.oei_erase ? "yes" : "no");
  outfile->Printf("   TEI FILE      =   %6d      TEI ERASE     =   %6s\n", 
      Params.tei_file, Params.tei_erase ? "yes" : "no");
  outfile->Printf("   OPDM FILE     =   %6d      OPDM ERASE    =   %6s\n", 
      Params.lag_file, Params.opdm_erase ? "yes" : "no");
  outfile->Printf("   TPDM FILE     =   %6d      TPDM ERASE    =   %6s\n", 
      Params.tpdm_file, Params.tpdm_erase ? "yes" : "no");
  outfile->Printf("   LAG FILE      =   %6d      LAG ERASE     =   %6s\n", 
      Params.lag_file, Params.lag_erase ? "yes" : "no");
  outfile->Printf("   DIIS START    =   %6d      DIIS FREQ     =   %6d\n", 
      Params.diis_start, Params.diis_freq);
  outfile->Printf("   DIIS MIN VECS =   %6d      DIIS MAX VECS =   %6d\n", 
      Params.diis_min_vecs, Params.diis_max_vecs);
  outfile->Printf("   SCALE STEP    =   %6.2E    MAX STEP      =   %6.2lf\n",
      Params.scale_step, Params.step_max);
  outfile->Printf("   DO_LEVEL SHIFT   =   %6s   SHIFT         =   %6.2lf\n",
      Params.level_shift ? "yes" : "no", Params.shift);
  outfile->Printf("   USE FZC H     =   %6s      HESSIAN       = %-12s\n",
      Params.use_fzc_h ? "yes" : "no", Params.hessian.c_str());
  outfile->Printf("\n") ;
  //fflush(outfile);
}

}} // end namespace psi::detcas

