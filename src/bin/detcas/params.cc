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
    Params.convergence = options.get_double("R_CONVERGENCE");
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
  Params.scale_grad = options.get_bool("SCALE_GRAD");

  Params.diis_start = options.get_int("DIIS_START");
  Params.diis_freq = options.get_int("DIIS_FREQ");
  Params.diis_min_vecs = options.get_int("DIIS_MIN_VECS");
  Params.diis_max_vecs = options.get_int("DIIS_MAX_VECS");

  Params.scale_step = options.get_double("SCALE_STEP");
  Params.use_fzc_h = options.get_bool("USE_FZC_H");
  Params.invert_hessian = options.get_bool("INVERT_HESSIAN");
  Params.hessian = option.get_str("HESSIAN");
 
  Params.level_shift = option.get_bool("LEVEL_SHIFT");
  Params.shift = option.get_double("SHIFT");
  Params.determ_min = option.get_double("DETERM_MIN");
  Params.step_max = option.get_double("MAX_STEP");
  Params.use_thetas = option.get_bool("USE_THETAS");
  Params.force_step = option.get_bool("FORCE_STEP");
  Params.force_pair = option.get_int("FORCE_PAIR");
  Params.force_value = option.get_double("FORCE_VALUE");
  Params.scale_act_act = option.get_double("SCALE_ACT_ACT");
  Params.bfgs = option.get_bool("BFGS");
  Params.ds_hessian = option.get_bool("DS_HESSIAN");

  Params.level_shift = options.get_bool("LEVEL_SHIFT");
  Params.shift = options.get_double("SHIFT");
  Params.determ_min = options.get_double("DETERM_MIN");
  Params.step_max = options.get_double("MAX_STEP");
  Params.use_thetas = options.get_bool("USE_THETAS");
  Params.force_step = options.get_bool("FORCE_STEP");
  Params.force_pair = options.get_int("FORCE_PAIR");
  Params.force_value = options.get_double("FORCE_VALUE");
  /* at the moment, the following only work for diagonal (non-YY) Hessians */
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
  fprintf(outfile, "\n") ;
  fprintf(outfile, "PARAMETERS: \n") ;
  fprintf(outfile, "   PRINT         =   %6d      PRINT_MOS     =   %6s\n", 
      Params.print_lvl, Params.print_mos ? "yes" : "no");
  fprintf(outfile, "   CONVERGENCE   =   %6d      E CONVERG     =   %6d\n",
      Params.rms_grad_convergence, Params.energy_convergence);
  fprintf(outfile, "   IGNORE_RAS_RAS=   %6s      IGNORE_FZ     =   %6s\n", 
      Params.ignore_ras_ras ? "yes" : "no", Params.ignore_fz ? "yes" : "no") ;
  fprintf(outfile, "   OEI FILE      =   %6d      OEI ERASE     =   %6s\n", 
      Params.oei_file, Params.oei_erase ? "yes" : "no");
  fprintf(outfile, "   TEI FILE      =   %6d      TEI ERASE     =   %6s\n", 
      Params.tei_file, Params.tei_erase ? "yes" : "no");
  fprintf(outfile, "   OPDM FILE     =   %6d      OPDM ERASE    =   %6s\n", 
      Params.lag_file, Params.opdm_erase ? "yes" : "no");
  fprintf(outfile, "   TPDM FILE     =   %6d      TPDM ERASE    =   %6s\n", 
      Params.tpdm_file, Params.tpdm_erase ? "yes" : "no");
  fprintf(outfile, "   LAG FILE      =   %6d      LAG ERASE     =   %6s\n", 
      Params.lag_file, Params.lag_erase ? "yes" : "no");
  fprintf(outfile, "   DIIS START    =   %6d      DIIS FREQ     =   %6d\n", 
      Params.diis_start, Params.diis_freq);
  fprintf(outfile, "   DIIS MIN VECS =   %6d      DIIS MAX VECS =   %6d\n", 
      Params.diis_min_vecs, Params.diis_max_vecs);
  fprintf(outfile, "   SCALE STEP    =   %6.2E    MAX STEP      =   %6.2lf\n",
      Params.scale_step, Params.step_max);
  fprintf(outfile, "   LEVEL SHIFT   =   %6s      SHIFT         =   %6.2lf\n",
      Params.level_shift ? "yes" : "no", Params.shift);
  fprintf(outfile, "   USE FZC H     =   %6s      HESSIAN       = %-12s\n",
      Params.use_fzc_h ? "yes" : "no", Params.hessian);
  fprintf(outfile, "\n") ;
  fflush(outfile) ;
}

}} // end namespace psi::detcas

