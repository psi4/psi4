/*! \file
    \ingroup OPTKING
    \brief This function constructs the B matrix for a set of salcs
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

namespace psi { //namespace optking {

double **compute_B(const simples_class & simples, const salc_set & symm) {
  int i,j,a,b,c,d, simple, sub_index, sub_index2, atom, xyz;
  double **B, coeff, prefactor, weight;
  Intco_type intco_type;

  B = block_matrix(symm.get_num(),3*optinfo.natom);

  for (i=0;i<symm.get_num();++i) {
    prefactor = symm.get_prefactor(i);

    for (j=0;j<symm.get_length(i);++j) {
      simple = symm.get_simple(i,j);
      coeff = symm.get_coeff(i,j);
      simples.locate_id(simple, &intco_type, &sub_index, &sub_index2);
      // fprintf(outfile,"intcotype %d, sub_index %d\n", intco_type, sub_index);

      if (intco_type == STRE) {
        a = simples.get_atom(STRE, sub_index, 0);
        b = simples.get_atom(STRE, sub_index, 1);
        for (xyz=0;xyz<3;++xyz) {
          B[i][3*a+xyz] += prefactor * coeff * simples.get_s(STRE, sub_index, 0, xyz);
          B[i][3*b+xyz] += prefactor * coeff * simples.get_s(STRE, sub_index, 1, xyz);
        }
      }
      else if (intco_type == BEND) {
        a = simples.get_atom(BEND, sub_index, 0);
        b = simples.get_atom(BEND, sub_index, 1);
        c = simples.get_atom(BEND, sub_index, 2);
        for (xyz=0;xyz<3;++xyz) {
          B[i][3*a+xyz] += prefactor * coeff * simples.get_s(BEND, sub_index, 0, xyz);
          B[i][3*b+xyz] += prefactor * coeff * simples.get_s(BEND, sub_index, 1, xyz);
          B[i][3*c+xyz] += prefactor * coeff * simples.get_s(BEND, sub_index, 2, xyz);
        }
      }
      else if (intco_type == TORS) {
        a = simples.get_atom(TORS, sub_index, 0);
        b = simples.get_atom(TORS, sub_index, 1);
        c = simples.get_atom(TORS, sub_index, 2);
        d = simples.get_atom(TORS, sub_index, 3);
        for (xyz=0;xyz<3;++xyz) {
          B[i][3*a+xyz] += prefactor * coeff * simples.get_s(TORS, sub_index, 0, xyz);
          B[i][3*b+xyz] += prefactor * coeff * simples.get_s(TORS, sub_index, 1, xyz);
          B[i][3*c+xyz] += prefactor * coeff * simples.get_s(TORS, sub_index, 2, xyz);
          B[i][3*d+xyz] += prefactor * coeff * simples.get_s(TORS, sub_index, 3, xyz);
        }
      }
      else if (intco_type == OUT) {
        a = simples.get_atom(OUT, sub_index, 0);
        b = simples.get_atom(OUT, sub_index, 1);
        c = simples.get_atom(OUT, sub_index, 2);
        d = simples.get_atom(OUT, sub_index, 3);
        for (xyz=0;xyz<3;++xyz) {
          B[i][3*a+xyz] += prefactor * coeff * simples.get_s(OUT, sub_index, 0, xyz);
          B[i][3*b+xyz] += prefactor * coeff * simples.get_s(OUT, sub_index, 1, xyz);
          B[i][3*c+xyz] += prefactor * coeff * simples.get_s(OUT, sub_index, 2, xyz);
          B[i][3*d+xyz] += prefactor * coeff * simples.get_s(OUT, sub_index, 3, xyz);
        }
      }
      else if (intco_type == LINB) {
        a = simples.get_atom(LINB, sub_index, 0);
        b = simples.get_atom(LINB, sub_index, 1);
        c = simples.get_atom(LINB, sub_index, 2);
        for (xyz=0;xyz<3;++xyz) {
          B[i][3*a+xyz] += prefactor * coeff * simples.get_s(LINB, sub_index, 0, xyz);
          B[i][3*b+xyz] += prefactor * coeff * simples.get_s(LINB, sub_index, 1, xyz);
          B[i][3*c+xyz] += prefactor * coeff * simples.get_s(LINB, sub_index, 2, xyz);
        }
      }
      else if (intco_type == FRAG) {
        for (atom=0; atom<simples.get_natom(FRAG, sub_index, FRAG_A); ++atom) {
          a = simples.get_atom(FRAG, sub_index, atom, FRAG_A);                
          for (xyz=0;xyz<3;++xyz)
            B[i][3*a+xyz] += prefactor * coeff * simples.get_s(FRAG, sub_index, atom, xyz, sub_index2, FRAG_A);
        }
        for (atom=0; atom<simples.get_natom(FRAG, sub_index, FRAG_B); ++atom) {
          b = simples.get_atom(FRAG, sub_index, atom, FRAG_B);
          for (xyz=0;xyz<3;++xyz)
            B[i][3*b+xyz] += prefactor * coeff * simples.get_s(FRAG, sub_index, atom, xyz, sub_index2, FRAG_B);
        }
      }
    }
  }

  //fprintf(outfile, "B matrix\n");
  //print_mat2(B, symm.get_num(),3*optinfo.natom, outfile);

  return B;
}

}//} /* namespace psi::optking */
