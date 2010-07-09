/* STEP_LIMIT limits maximum change in primitive to value in au or rad */

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"

namespace psi { //namespace optking {

// no internal coordinate value can change more than optinfo.step_limit
// dq is passed in in Angstroms or radians

void step_limit(const simples_class & simples, const salc_set &symm, double *dq) {

  int i, j, dim, max_i, simple, sub_index, sub_index2;
  double min_scale = 1.0, scale, coeff, prefactor, tval, dq_simple;
  double inv_limit, tval2, R, dR, inv_R_min, inv_R_max, step_limit;
  Intco_type intco_type;
  
  dim = symm.get_num();
  step_limit = optinfo.step_limit;

  for (i=0; i<dim; ++i) { // loop over symm SALCS vectors
    prefactor = symm.get_prefactor(i);
    for (j=0;j<symm.get_length(i);++j) { // loop over simples in each salc
      coeff = symm.get_coeff(i,j);
      dq_simple = prefactor * coeff * dq[i]; // change in value of primitive
      scale = 1.0;

      // dq is in Angstroms or radians
      // If coordinate is a strech, and we need to convert Angstroms -> bohr
      simple = symm.get_simple(i,j);
      simples.locate_id(simple,&intco_type,&sub_index,&sub_index2);

      if (intco_type != FRAG) { // regular, intrafragment coordinates
        if (intco_type == STRE) dq_simple /= _bohr2angstroms;

        if (fabs(dq_simple) > step_limit)
          scale = step_limit / fabs(dq_simple);
      }
      else if (sub_index2 > 0) { // interfragment angles
        if (fabs(dq_simple) > step_limit)
          scale = step_limit / fabs(dq_simple);
      }
      else if (!optinfo.frag_dist_rho) { // interfragment R(A-B)
        dq_simple /= _bohr2angstroms;
        if (fabs(dq_simple) > step_limit)
          scale = step_limit / fabs(dq_simple);
      }
      else {  // interfragment 1/R(A-B)
        dq_simple *= _bohr2angstroms;
        // fprintf(outfile, "dq_simple (1/au) %15.10lf\n", dq_simple);
        R = 1.0 / simples.get_val(FRAG, sub_index, 0) * _bohr2angstroms;
        // fprintf(outfile, "R in au %15.10lf\n", R);
        inv_R_min = - step_limit / (R * (R + step_limit));
        // fprintf(outfile, "1/R min in au %15.10lf\n", inv_R_min);
        inv_R_max = step_limit / (R * (R - step_limit));
        // fprintf(outfile, "1/R max in au %15.10lf\n", inv_R_max);

        if (dq_simple < inv_R_min) 
          scale = inv_R_min / dq_simple;
        else if (dq_simple > inv_R_max)
          scale = inv_R_max / dq_simple;
      }

      if (scale < min_scale) {
        min_scale = scale;
        max_i = i;
      }
    }
  }

  // avoid piddly symmetry breaking
  for (i=0;i<dim;++i)
    if (fabs(dq[i]) < MIN_DQ_STEP) dq[i] = 0.0;

  if (min_scale < 1) {
    fprintf(outfile,"\nMaximum change in SALC %d exceeds STEP_LIMIT\n", max_i+1);
    fprintf(outfile,"Scaling displacements by %lf\n",min_scale);
    for (i=0;i<dim;++i)
      dq[i] *= min_scale;   
  }

  return;
}

}//}
