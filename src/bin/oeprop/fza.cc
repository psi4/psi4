/*! \file
    \ingroup OEPROP
    \brief fza returns extraplated values of QED function F(Z*alpha), tables for which have been tabulated 
        by Mohr and co-workers.  These are the n=1 values (there is a small n dependence).  See
        Pyykko et al Phys. Rev. A 63, 024502.
*/

#define EXTERN
#include "includes.h"
#include "globals.h"

namespace psi { namespace oeprop {

/* F(Z alpha) values for 1 electron tabulated by Mohr for n=1*/
#define NENTRIES 35
static double FZA[NENTRIES][2] = {
{ 1.0, 10.3168},
{ 2.0, 8.5283},
{ 3.0, 7.5045},
{ 4.0, 6.7928},
{ 5.0, 6.2516},
{10.0, 4.6542},
{15.0, 3.8014},
{20.0, 3.2463},
{25.0, 2.8501},
{26.0, 2.7839},
{30.0, 2.5520},
{35.0, 2.3199},
{36.0, 2.2796},
{40.0, 2.1352},
{45.0, 1.9859},
{50.0, 1.8643},
{54.0, 1.7831},
{55.0, 1.7648},
{60.0, 1.6838},
{65.0, 1.6186},
{66.0, 1.6073},
{70.0, 1.5674},
{75.0, 1.5289},
{79.0, 1.5070},
{80.0, 1.5027},
{82.0, 1.4956},
{83.0, 1.4928},
{85.0, 1.4887},
{90.0, 1.4875},
{92.0, 1.4909},
{95.0, 1.5005},
{100.0, 1.5301},
{105.0, 1.5809},
{110.0, 1.6600},
{200.0, 1.6600} /* assume no more change */
};

double F_Zalpha(double Z_alpha) {
  int i,j;
  double x1, x2, y1, y2, tval;

  if (Z_alpha < 1.0) {
    throw PsiException("Error: Z_alpha must be at least 1", __FILE__, __LINE__);
//    fprintf(outfile,"Error: Z_alpha must be at least 1.\n");
//    exit(PSI_RETURN_FAILURE);
  }

  /* check for exact match first */
  for (i=0; i<NENTRIES; ++i) {
    if (Z_alpha == FZA[i][0])
      return FZA[i][1];
  }

  for (i=0; i<NENTRIES; ++i) {
    if (Z_alpha < FZA[i][0]) {
      x1 = FZA[i-1][0];
      x2 = FZA[i][0];
      y1 = FZA[i-1][1];
      y2 = FZA[i][1];
      tval = (y2-y1)/(x2-x1) * (Z_alpha-x1) + y1;
      return tval;
    }
  }
  return FZA[33][1];
}

}}

