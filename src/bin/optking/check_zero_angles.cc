/* CHECK_ZERO_ANGLES checks to see if interfragment angles (theta-A and theta_B) are
     passing through 180 or 0 - which spells big problems
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

namespace psi { namespace optking {

// angle bends cannot go < 0
// dq is in radians

void check_zero_angles(const simples_class & simples, const salc_set & symm, double *dq) {
  int i, j, dim, simple, sub_index, sub_index2;
  double coeff, prefactor, dq_simple, val;
  Intco_type intco_type;
  char error[100];
  
  dim = symm.get_num();

  // find bend angles with range {0,180}
  for (i=0; i<dim; ++i) { // loop over symm SALCS vectors
    prefactor = symm.get_prefactor(i);
    for (j=0;j<symm.get_length(i);++j) { // loop over simples in each salc
      coeff = symm.get_coeff(i,j);
      dq_simple = prefactor * coeff * dq[i]; // change in value of primitive

      simple = symm.get_simple(i,j);
      simples.locate_id(simple,&intco_type,&sub_index,&sub_index2);

      if (intco_type == BEND) {
        val = simples.get_val_A_or_rad(BEND, sub_index); // in rad
        if (val + dq_simple < 0.0)
          throw("Bond angle passing through zero. Try new bonds, angles, or considering higher symmetry");
      }
      else if ((intco_type == FRAG) && (sub_index2 == 1)) {
        val = simples.get_val_A_or_rad(FRAG, sub_index, sub_index2); // in rad
        if (val + dq_simple < 0.0)
          throw("Theta-A angle passing through zero. Try adding 180 to tau and chi-A or higher symmetry");
      }
      else if ((intco_type == FRAG) && (sub_index2 == 2)) {
        val = simples.get_val_A_or_rad(FRAG, sub_index, sub_index2);
        if (val + dq_simple < 0.0)
          throw("Theta-B angle passing through zero. Try adding 180 to tau and chi-B or higher symmetry");
      }
    }
  }
  return;
}

}}
