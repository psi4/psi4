/*! \file symm.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<stdlib.h>
#include<libipv1/ip_lib.h>
#include<libciomr/libciomr.h>
#include<libchkpt/chkpt.h>

#include<libint/libint.h>
#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace CINTS {

/*-------------------------------
  Explicit function declarations
 -------------------------------*/
static void init_dp_table(void);

void init_symmetry()
{
  int i, j, count;

  Symmetry.symlabel = chkpt_rd_sym_label();
  Symmetry.nirreps = chkpt_rd_nirreps();
  Symmetry.num_so = chkpt_rd_nso();
  Symmetry.num_unique_atoms = chkpt_rd_num_unique_atom();
  Symmetry.num_unique_shells = chkpt_rd_num_unique_shell();

  Symmetry.atom_positions = chkpt_rd_atom_position();
  Symmetry.ua2a = chkpt_rd_ua2a();
  Symmetry.us2s = chkpt_rd_us2s();
  Symmetry.sopi = chkpt_rd_sopi();
  Symmetry.sym_oper = chkpt_rd_symoper();
  Symmetry.irr_labels = chkpt_rd_irr_labs();
  Symmetry.ict = chkpt_rd_ict();
  Symmetry.cartrep = chkpt_rd_cartrep();
  Symmetry.cdsalcpi = chkpt_rd_cdsalcpi();
  Symmetry.cdsalc2cd = chkpt_rd_cdsalc2cd();
  Symmetry.cdsalc_ioffset = init_int_array(Symmetry.nirreps);
  Symmetry.cdsalc_ioffset[0] = 0;
  for(i=1;i<Symmetry.nirreps;i++)
    Symmetry.cdsalc_ioffset[i] = Symmetry.cdsalc_ioffset[i-1] + Symmetry.cdsalcpi[i-1];

  if (Symmetry.nirreps) {
  /* Symmetry.dp_table = */ init_dp_table();
#if SCF_ONLY
    if (Symmetry.nirreps < 4)
      UserOptions.scf_only = 0;
    if (UserOptions.scf_only) {
      Symmetry.so2symblk = init_int_array(Symmetry.num_so);
      count = 0;
      for(i=0;i<Symmetry.nirreps;i++)
	for(j=0;j<Symmetry.sopi[i];j++)
	  Symmetry.so2symblk[count++] = i;
    }
#endif
  }

  return;
}


void cleanup_symmetry()
{
  if (Symmetry.nirreps > 1)
    free_int_matrix(Symmetry.dp_table);
  free(Symmetry.sopi);
  free_block(Symmetry.usotao);
  free(Symmetry.us2s);
  free(Symmetry.atom_positions);
  free(Symmetry.cartrep);

  return;
}


/*!----------------------------------------------------------------------
  Compute direct product multiplication table for the given point group
  NOTE: This matrix is really pointless at the moment, I left it here
  just in case
 ----------------------------------------------------------------------*/
void init_dp_table(void)
{
  int i,j;

  Symmetry.dp_table = init_int_matrix(Symmetry.nirreps,
				      Symmetry.nirreps);
  for(i=0;i<Symmetry.nirreps;i++)
    for(j=0;j<i;j++) {
      /*------------------------------------------------------------
	The line below works only in a case of Abelian point group!
	Have to do honest multiplication of rows of character table
	if non-Abelian groups to be used
       ------------------------------------------------------------*/
      Symmetry.dp_table[i][j] = Symmetry.dp_table[j][i] = i ^ j;
    }

  return;
}
};};
