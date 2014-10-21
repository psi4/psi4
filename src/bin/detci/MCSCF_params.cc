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
**   as convergence.  These are stored in the MCSCF_Parameters data structure.
*/
void get_parameters(Options &options)
{
  int i, errcod;
  char line1[133];
   
  MCSCF_Parameters.wfn = options.get_str("WFN");
  MCSCF_Parameters.dertype = options.get_str("DERTYPE");
  if (MCSCF_Parameters.dertype == "NONE") {
    MCSCF_Parameters.rms_grad_convergence = 1e-4;
    MCSCF_Parameters.energy_convergence = 1e-7;
  }
  else {
    MCSCF_Parameters.rms_grad_convergence = 1e-7;
    MCSCF_Parameters.energy_convergence = 1e-10;
  }

  if (options["R_CONVERGENCE"].has_changed()) {
    MCSCF_Parameters.rms_grad_convergence = options.get_double("R_CONVERGENCE");
  }
  if (options["E_CONVERGENCE"].has_changed()) {
    MCSCF_Parameters.energy_convergence = options.get_double("E_CONVERGENCE");
  }


  MCSCF_Parameters.filter_ints = 0;  /* assume we need all for MCSCF */
  MCSCF_Parameters.oei_file = PSIF_OEI;  /* contains frozen core operator */
  MCSCF_Parameters.tei_file = PSIF_MO_TEI;
  MCSCF_Parameters.opdm_file = PSIF_MO_OPDM;
  MCSCF_Parameters.tpdm_file = PSIF_MO_TPDM;
  MCSCF_Parameters.lag_file = PSIF_MO_LAG;

  MCSCF_Parameters.ignore_fz = true;     /* ignore frozen orbitals for ind pairs? */
  
  if ((MCSCF_Parameters.wfn == "CASSCF") || (MCSCF_Parameters.wfn == "DETCAS"))
    MCSCF_Parameters.ignore_ras_ras = true;   /* ignore RAS/RAS independent pairs? */
  else
    MCSCF_Parameters.ignore_ras_ras = false;

  MCSCF_Parameters.print_lvl = options.get_int("PRINT");
  MCSCF_Parameters.print_mos = options.get_bool("PRINT_MOS");

  MCSCF_Parameters.oei_erase = options.get_bool("OEI_ERASE");
  MCSCF_Parameters.tei_erase = options.get_bool("TEI_ERASE");
  MCSCF_Parameters.ignore_fz = options.get_bool("IGNORE_FZ");
  MCSCF_Parameters.scale_grad = true; /* scale orb grad by inverse Hessian? */

  MCSCF_Parameters.diis_start = options.get_int("DIIS_START");
  MCSCF_Parameters.diis_freq = options.get_int("DIIS_FREQ");
  MCSCF_Parameters.diis_min_vecs = options.get_int("DIIS_MIN_VECS");
  MCSCF_Parameters.diis_max_vecs = options.get_int("DIIS_MAX_VECS");

  MCSCF_Parameters.scale_step = options.get_double("SCALE_STEP");
  MCSCF_Parameters.use_fzc_h = options.get_bool("USE_FZC_H");
  MCSCF_Parameters.invert_hessian = options.get_bool("INVERT_HESSIAN");
  MCSCF_Parameters.hessian = options.get_str("HESSIAN");
 
  MCSCF_Parameters.level_shift = options.get_bool("DO_LEVEL_SHIFT");
  MCSCF_Parameters.shift = options.get_double("SHIFT");
  MCSCF_Parameters.determ_min = options.get_double("DETERM_MIN");
  MCSCF_Parameters.step_max = options.get_double("STEP_MAX");
  MCSCF_Parameters.use_thetas = options.get_bool("USE_THETAS");
  MCSCF_Parameters.force_step = options.get_bool("FORCE_STEP");
  MCSCF_Parameters.force_pair = options.get_int("FORCE_PAIR");
  MCSCF_Parameters.force_value = options.get_double("FORCE_VALUE");
  MCSCF_Parameters.scale_act_act = options.get_double("SCALE_ACT_ACT");
  MCSCF_Parameters.bfgs = options.get_bool("BFGS");
  MCSCF_Parameters.ds_hessian = options.get_bool("DS_HESSIAN");

}


/*
** print_parameters(): Function prints the program's running parameters
**   found in the MCSCF_Parameters structure.
*/
void print_parameters(void)
{
  outfile->Printf("\n") ;
  outfile->Printf("PARAMETERS: \n") ;
  outfile->Printf("   PRINT          =   %6d      PRINT_MOS     =   %6s\n", 
      MCSCF_Parameters.print_lvl, MCSCF_Parameters.print_mos ? "yes" : "no");
  outfile->Printf("   R_CONVERGENCE  =   %6.2e    E CONVERG     =   %6.2e\n",
      MCSCF_Parameters.rms_grad_convergence, MCSCF_Parameters.energy_convergence);
  outfile->Printf("   IGNORE_RAS_RAS =   %6s      IGNORE_FZ     =   %6s\n", 
      MCSCF_Parameters.ignore_ras_ras ? "yes" : "no", MCSCF_Parameters.ignore_fz ? "yes" : "no") ;
  outfile->Printf("   OEI FILE       =   %6d      OEI ERASE     =   %6s\n", 
      MCSCF_Parameters.oei_file, MCSCF_Parameters.oei_erase ? "yes" : "no");
  outfile->Printf("   TEI FILE       =   %6d      TEI ERASE     =   %6s\n", 
      MCSCF_Parameters.tei_file, MCSCF_Parameters.tei_erase ? "yes" : "no");
  outfile->Printf("   OPDM FILE      =   %6d      OPDM ERASE    =   %6s\n", 
      MCSCF_Parameters.lag_file, MCSCF_Parameters.opdm_erase ? "yes" : "no");
  outfile->Printf("   TPDM FILE      =   %6d      TPDM ERASE    =   %6s\n", 
      MCSCF_Parameters.tpdm_file, MCSCF_Parameters.tpdm_erase ? "yes" : "no");
  outfile->Printf("   LAG FILE       =   %6d      LAG ERASE     =   %6s\n", 
      MCSCF_Parameters.lag_file, MCSCF_Parameters.lag_erase ? "yes" : "no");
  outfile->Printf("   DIIS START     =   %6d      DIIS FREQ     =   %6d\n", 
      MCSCF_Parameters.diis_start, MCSCF_Parameters.diis_freq);
  outfile->Printf("   DIIS MIN VECS  =   %6d      DIIS MAX VECS =   %6d\n", 
      MCSCF_Parameters.diis_min_vecs, MCSCF_Parameters.diis_max_vecs);
  outfile->Printf("   SCALE STEP     =   %6.2E    MAX STEP      =   %6.2lf\n",
      MCSCF_Parameters.scale_step, MCSCF_Parameters.step_max);
  outfile->Printf("   DO_LEVEL SHIFT =   %6s      SHIFT         =   %6.2lf\n",
      MCSCF_Parameters.level_shift ? "yes" : "no", MCSCF_Parameters.shift);
  outfile->Printf("   USE FZC H      =   %6s      HESSIAN       = %-12s\n",
      MCSCF_Parameters.use_fzc_h ? "yes" : "no", MCSCF_Parameters.hessian.c_str());
  outfile->Printf("\n") ;
  //fflush(outfile);
}

}} // end namespace psi::detcas

