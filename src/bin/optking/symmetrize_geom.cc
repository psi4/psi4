/*! \file 
    \ingroup OPTKING
    \brief symmetrize_geom - symmetrizes a cartesian geometry
    \param double *geom
*/

#define EXTERN
#include "def.h"
#undef EXTERN

#include <libchkpt/chkpt.h>

namespace psi { namespace optking {

void symmetrize_geom(double *x) {
  int ua, op, xyz, num_uniques, i;
  int atom1, atom2, *ua2a, **ict, diag_ind, natom;
  double **cartrep, *x_temp;
  int stab_order, nirreps;

  chkpt_init(PSIO_OPEN_OLD);
  num_uniques = chkpt_rd_num_unique_atom();
  ua2a = chkpt_rd_ua2a();
  nirreps = chkpt_rd_nirreps();
  ict = chkpt_rd_ict();
  cartrep = chkpt_rd_cartrep();
  natom = chkpt_rd_natom();
  chkpt_close();

  x_temp = init_array(natom*3);

  for(ua=0; ua<num_uniques; ua++) {
    atom1 = ua2a[ua];
    stab_order = 0;
    for(op=0; op < nirreps; op++) {
      atom2 = ict[op][atom1] - 1;
      if (atom1 == atom2)
        stab_order++;
    }
    for(op=0; op < nirreps; op++) {
      atom2 = ict[op][atom1] - 1;
      for(xyz=0; xyz<3; xyz++) {
        diag_ind = xyz*3 + xyz;
        x_temp[3*atom2+xyz] += cartrep[op][diag_ind] * x[3*atom1+xyz] / stab_order;
      }
    }
  }

  for (i=0;i<3*natom;++i)
    x[i] = x_temp[i];

  free_array(x_temp);
  free_int_array(ua2a);
  free_int_matrix(ict);
  free_matrix(cartrep);
  return;
}

}} /* namespace psi::optking */
