/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "defines.h"
#include "params.h"
#include "mo_overlap.h"

namespace psi { namespace DBOC {

extern Params_t Params;
extern "C" FILE *outfile;

extern void done(const char *);
extern double eval_rhf_derwfn_overlap(DisplacementIndex LDisp, DisplacementIndex RDisp);
extern double eval_rohf_derwfn_overlap(DisplacementIndex LDisp, DisplacementIndex RDisp);
extern double eval_uhf_derwfn_overlap(DisplacementIndex LDisp, DisplacementIndex RDisp);
extern double eval_rci_derwfn_overlap(DisplacementIndex LDisp, DisplacementIndex RDisp);
extern double eval_roci_derwfn_overlap(DisplacementIndex LDisp, DisplacementIndex RDisp);

double eval_derwfn_overlap(bool symm)
{
  // Pointer to the function that we need to use
  double (*eval_overlap)(DisplacementIndex, DisplacementIndex);

  double S;

  if (!strcmp(Params.wfn,"SCF")) {
    if (Params.reftype == Params_t::rhf) {
      eval_overlap = eval_rhf_derwfn_overlap;
    }
    else if (Params.reftype == Params_t::rohf) {
      eval_overlap = eval_rohf_derwfn_overlap;
    }
    else if (Params.reftype == Params_t::uhf) {
      eval_overlap = eval_uhf_derwfn_overlap;
    }
    else
      done("This HF SCF method is not supported at the moment");
  }
  else if (!strcmp(Params.wfn,"CCSD")) {
    done("CCSD method with this reference is not supported at the moment");
  }
  else if (!strcmp(Params.wfn,"DETCI") || !strcmp(Params.wfn,"DETCAS")) {
    if (Params.reftype == Params_t::rhf) {
      eval_overlap = eval_rci_derwfn_overlap;
    }
    else if (Params.reftype == Params_t::rohf) {
      eval_overlap = eval_roci_derwfn_overlap;
    }
    else
      done("CI method with this reference is not supported at the moment");
  }

  //
  // Only need (+Delta|-Delta) overlap if using a 2-point formula
  //
  if (Params.disp_per_coord == 2) {
    double S_P1_M1 = eval_overlap(PlusDelta,MinusDelta);
    if (Params.print_lvl >= PrintLevels::print_contrib) {
      fprintf(outfile,"  +1 -1 wave function overlap = %25.15lf\n",S_P1_M1);
    }
    S = (1.0-S_P1_M1)/(2.0*Params.delta*Params.delta);
  }
  //
  // Need up to 6 different overlaps if using a 4-point formula
  //
  else if (Params.disp_per_coord == 4) {
    double S_P1_M1 = eval_overlap(PlusDelta,MinusDelta);
    double S_P2_M2 = eval_overlap(Plus2Delta,Minus2Delta);
    double S_P2_M1 = eval_overlap(Plus2Delta,MinusDelta);
    double S_P2_P1 = eval_overlap(Plus2Delta,PlusDelta);
    double S_M2_P1, S_M2_M1;
    if (!symm) {
      S_M2_P1 = eval_overlap(Minus2Delta,PlusDelta);
      S_M2_M1 = eval_overlap(Minus2Delta,MinusDelta);
    }
    else {
      S_M2_P1 = S_P2_M1;
      S_M2_M1 = S_P2_P1;
    }

    if (Params.print_lvl >= PrintLevels::print_contrib) {
      fprintf(outfile,"  +1 -1 wave function overlap = %25.15lf\n",S_P1_M1);
      fprintf(outfile,"  +2 -2 wave function overlap = %25.15lf\n",S_P2_M2);
      fprintf(outfile,"  +2 -1 wave function overlap = %25.15lf\n",S_P2_M1);
      fprintf(outfile,"  -2 +1 wave function overlap = %25.15lf\n",S_M2_P1);
      fprintf(outfile,"  +2 +1 wave function overlap = %25.15lf\n",S_P2_P1);
      fprintf(outfile,"  -2 -1 wave function overlap = %25.15lf\n",S_M2_M1);
    }

    S = (128.0*(1.0-S_P1_M1) + 2.0*(1.0-S_P2_M2) + 16.0*((S_P2_M1 - S_P2_P1) + (S_M2_P1 - S_M2_M1)))/
        (144.0 * Params.delta*Params.delta);
  }

  return S;
}  

}} // namespace psi::DBOC
