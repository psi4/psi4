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
#include "MCSCF_globals.h"
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
   
  Parameters.wfn = options.get_str("WFN");
  Parameters.dertype = options.get_str("DERTYPE");
  if (Parameters.dertype == "NONE") {
    Parameters.rms_grad_convergence = 1e-4;
    Parameters.energy_convergence = 1e-7;
  }
  else {
    Parameters.rms_grad_convergence = 1e-7;
    Parameters.energy_convergence = 1e-10;
  }

  if (options["R_CONVERGENCE"].has_changed()) {
    Parameters.rms_grad_convergence = options.get_double("R_CONVERGENCE");
  }
  if (options["E_CONVERGENCE"].has_changed()) {
    Parameters.energy_convergence = options.get_double("E_CONVERGENCE");
  }


  Parameters.filter_ints = 0;  /* assume we need all for MCSCF */
  Parameters.oei_file = PSIF_OEI;  /* contains frozen core operator */
  Parameters.tei_file = PSIF_MO_TEI;
  Parameters.opdm_file = PSIF_MO_OPDM;
  Parameters.tpdm_file = PSIF_MO_TPDM;
  Parameters.lag_file = PSIF_MO_LAG;

  Parameters.ignore_fz = true;     /* ignore frozen orbitals for ind pairs? */
  
  if ((Parameters.wfn == "CASSCF") || (Parameters.wfn == "DETCAS"))
    Parameters.ignore_ras_ras = true;   /* ignore RAS/RAS independent pairs? */
  else
    Parameters.ignore_ras_ras = false;

  Parameters.print_lvl = options.get_int("PRINT");
  Parameters.print_mos = options.get_bool("PRINT_MOS");

  Parameters.oei_erase = options.get_bool("OEI_ERASE");
  Parameters.tei_erase = options.get_bool("TEI_ERASE");
  Parameters.ignore_fz = options.get_bool("IGNORE_FZ");
  Parameters.scale_grad = true; /* scale orb grad by inverse Hessian? */

  Parameters.diis_start = options.get_int("DIIS_START");
  Parameters.diis_freq = options.get_int("DIIS_FREQ");
  Parameters.diis_min_vecs = options.get_int("DIIS_MIN_VECS");
  Parameters.diis_max_vecs = options.get_int("DIIS_MAX_VECS");

  Parameters.scale_step = options.get_double("SCALE_STEP");
  Parameters.use_fzc_h = options.get_bool("USE_FZC_H");
  Parameters.invert_hessian = options.get_bool("INVERT_HESSIAN");
  Parameters.hessian = options.get_str("HESSIAN");
 
  Parameters.level_shift = options.get_bool("DO_LEVEL_SHIFT");
  Parameters.shift = options.get_double("SHIFT");
  Parameters.determ_min = options.get_double("DETERM_MIN");
  Parameters.step_max = options.get_double("STEP_MAX");
  Parameters.use_thetas = options.get_bool("USE_THETAS");
  Parameters.force_step = options.get_bool("FORCE_STEP");
  Parameters.force_pair = options.get_int("FORCE_PAIR");
  Parameters.force_value = options.get_double("FORCE_VALUE");
  Parameters.scale_act_act = options.get_double("SCALE_ACT_ACT");
  Parameters.bfgs = options.get_bool("BFGS");
  Parameters.ds_hessian = options.get_bool("DS_HESSIAN");

}


/*
** print_parameters(): Function prints the program's running parameters
**   found in the Parameters structure.
*/
void print_parameters(void)
{
  outfile->Printf("\n") ;
  outfile->Printf("PARAMETERS: \n") ;
  outfile->Printf("   PRINT          =   %6d      PRINT_MOS     =   %6s\n", 
      Parameters.print_lvl, Parameters.print_mos ? "yes" : "no");
  outfile->Printf("   R_CONVERGENCE  =   %6.2e    E CONVERG     =   %6.2e\n",
      Parameters.rms_grad_convergence, Parameters.energy_convergence);
  outfile->Printf("   IGNORE_RAS_RAS =   %6s      IGNORE_FZ     =   %6s\n", 
      Parameters.ignore_ras_ras ? "yes" : "no", Parameters.ignore_fz ? "yes" : "no") ;
  outfile->Printf("   OEI FILE       =   %6d      OEI ERASE     =   %6s\n", 
      Parameters.oei_file, Parameters.oei_erase ? "yes" : "no");
  outfile->Printf("   TEI FILE       =   %6d      TEI ERASE     =   %6s\n", 
      Parameters.tei_file, Parameters.tei_erase ? "yes" : "no");
  outfile->Printf("   OPDM FILE      =   %6d      OPDM ERASE    =   %6s\n", 
      Parameters.lag_file, Parameters.opdm_erase ? "yes" : "no");
  outfile->Printf("   TPDM FILE      =   %6d      TPDM ERASE    =   %6s\n", 
      Parameters.tpdm_file, Parameters.tpdm_erase ? "yes" : "no");
  outfile->Printf("   LAG FILE       =   %6d      LAG ERASE     =   %6s\n", 
      Parameters.lag_file, Parameters.lag_erase ? "yes" : "no");
  outfile->Printf("   DIIS START     =   %6d      DIIS FREQ     =   %6d\n", 
      Parameters.diis_start, Parameters.diis_freq);
  outfile->Printf("   DIIS MIN VECS  =   %6d      DIIS MAX VECS =   %6d\n", 
      Parameters.diis_min_vecs, Parameters.diis_max_vecs);
  outfile->Printf("   SCALE STEP     =   %6.2E    MAX STEP      =   %6.2lf\n",
      Parameters.scale_step, Parameters.step_max);
  outfile->Printf("   DO_LEVEL SHIFT =   %6s      SHIFT         =   %6.2lf\n",
      Parameters.level_shift ? "yes" : "no", Parameters.shift);
  outfile->Printf("   USE FZC H      =   %6s      HESSIAN       = %-12s\n",
      Parameters.use_fzc_h ? "yes" : "no", Parameters.hessian.c_str());
  outfile->Printf("\n") ;
  //fflush(outfile);
}

}} // end namespace psi::detcas

