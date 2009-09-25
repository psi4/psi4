#include <iostream>
#include <cmath>
#include <cstdio>

#include "scf.h"

extern FILE* outfile;

namespace psi{ namespace MCSCF{

void SCF::check_orthonormality()
{
  SBlockMatrix CSC("CSC",nirreps,sopi,sopi);
  transform(S,CSC,C);

  double    diagonal = 0.0;
  for(int h = 0; h < nirreps; ++h)
    for(int i = 0; i < sopi[h]; ++i)
      diagonal += fabs(CSC->get(h,i,i));

  double offdiagonal = 0.0;
  for(int h = 0; h < nirreps; ++h)
    for(int i = 0; i < sopi[h]; ++i)
      for(int j = i + 1; j < sopi[h]; ++j)
        offdiagonal += fabs(CSC->get(h,i,j));

  if((offdiagonal > 1.0e-8) || ((diagonal-double(nso)) > 1.0e-8)){
    fprintf(outfile,"\n\n  Warning: CSC has an orthonormality index of %lf",offdiagonal);
    fprintf(outfile,"\n  Trace(CSC) - nso = %lf",diagonal-nso);
    fprintf(outfile,"      Sum_i>j (CSC)ij  = %lf",offdiagonal);
  }else{
    fprintf(outfile,"\n  MOs orthonormality check passed.");
  }
}

}} // End namespace
