/*! \file
    \ingroup OPTKING
    \brief Returns the values of salcs given the simple internals
    and salc_set the value of the simple internals must already be computed.
    value is in Angstroms or radians
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"

namespace psi { //namespace optking {

double *compute_q(const simples_class & simples, const salc_set & symm) {
  int i, j, simple, sub_index, sub_index2;
  double *q, coeff, prefactor;
  Intco_type itype;

  q = init_array(symm.get_num());

  // q is built and returned in angstroms or radians
  for (i=0;i<symm.get_num();++i) {
    prefactor = symm.get_prefactor(i);

    for (j=0;j<symm.get_length(i);++j) {
      simple = symm.get_simple(i,j);
      coeff = symm.get_coeff(i,j);
      simples.locate_id(simple,&itype,&sub_index,&sub_index2);
      // sub_index2 is only used for multidimensional (interfragment coordinates)
      // fprintf(outfile,"itype %d sub_index %d sub_index2 %d\n", itype, sub_index, sub_index2);
      q[i] += prefactor * coeff * simples.get_val_A_or_rad(itype, sub_index, sub_index2);
    }

  }

  return q;
}

}//} /* namespace psi::optking */

