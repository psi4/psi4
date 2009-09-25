/*! \file
    \ingroup MVO
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include "MOInfo.h"
#include "globals.h"
#include "params.h"

namespace psi { namespace mvo {

/*
** transform_density
** 
** transforms the density from the MO basis to the SO basis for a given
** spin.  spin==0 is alpha, spin==1 is beta
*/
void transform_density(double **onepdm, double **psq_so, int spin)
{ 
  int i,j,k,l,dim_i,count,ibf;
  double **tmp_mat;

  /* CDS 2/02 */
  int irrep, mo_offset, so_offset, i_ci, j_ci, max_opi, errcod;
  int populated_orbs;
  int *docc, *socc, *frozen_docc, *frozen_uocc, *reorder; 
  double **scfvec, **opdm_blk;

  int nmo, *orbspi, *sopi, nirreps, nbfso;

  nirreps = moinfo.nirreps;
  sopi = moinfo.sopi;

  max_opi = 0;
  for (irrep=0; irrep<nirreps; irrep++)
    if (sopi[irrep] > max_opi) max_opi = sopi[irrep];
  opdm_blk = block_matrix(max_opi, max_opi); 
  tmp_mat = block_matrix(max_opi, max_opi);

  docc = moinfo.clsdpi;
  socc = moinfo.openpi;
  frozen_docc = moinfo.frdocc;
  frozen_uocc = moinfo.fruocc;
  reorder = moinfo.order;

  nmo = moinfo.nmo;
  nbfso = moinfo.nso;
  orbspi = moinfo.orbspi;

  /*
  if (params.print_lvl >= 4) {
    fprintf(outfile, "Reorder array:\n");
    for (i=0; i<nmo; i++) fprintf(outfile, "%d ", reorder[i]);
    fprintf(outfile, "\n");
  }
  */

  populated_orbs = nmo;
  for (irrep=0; irrep<nirreps; irrep++)
    populated_orbs -= frozen_uocc[irrep];

  mo_offset = 0;
  so_offset = 0;
  for (irrep=0; irrep<nirreps; irrep++) {
    if (orbspi[irrep] == 0) continue;
    for (i=0; i<orbspi[irrep]-frozen_uocc[irrep]; i++) {
      for (j=0; j<orbspi[irrep]-frozen_uocc[irrep]; j++) {
        i_ci = reorder[i+mo_offset];
        j_ci = reorder[j+mo_offset];
        opdm_blk[i][j] = onepdm[i_ci][j_ci];
      }
    }

    if (params.print_lvl >= 4) {
      fprintf(outfile, "Irrep %d (MO basis)\n", irrep);
      print_mat(opdm_blk,orbspi[irrep]-frozen_uocc[irrep], 
                orbspi[irrep]-frozen_uocc[irrep],outfile);
    }

    if (spin==0)
      scfvec = chkpt_rd_alpha_scf_irrep(irrep);
    else
      scfvec = chkpt_rd_beta_scf_irrep(irrep);

    if (params.print_lvl >= 4) {
      fprintf(outfile, "SCF Vector spin=%d\n", spin);
      print_mat(scfvec, sopi[irrep], orbspi[irrep], outfile);
    }

    mmult(opdm_blk,0,scfvec,1,tmp_mat,0,orbspi[irrep]-frozen_uocc[irrep],
          orbspi[irrep]-frozen_uocc[irrep],sopi[irrep],0);
    mmult(scfvec,0,tmp_mat,0,opdm_blk,0,sopi[irrep],
          orbspi[irrep]-frozen_uocc[irrep],sopi[irrep],0);
    if (params.print_lvl >= 4) {
      fprintf(outfile,"Irrep %d (SO basis)\n", irrep);
      print_mat(opdm_blk,sopi[irrep], sopi[irrep],outfile);
    }
    for (i=0; i<sopi[irrep]; i++)
      for (j=0; j<sopi[irrep]; j++)
        psq_so[i+so_offset][j+so_offset] = opdm_blk[i][j];
    mo_offset += orbspi[irrep];
    so_offset += sopi[irrep];

    free_block(scfvec);
  }

  if (params.print_lvl >= 4) {
    fprintf(outfile,"  Total %s density matrix in SO basis :\n",
      spin ? "beta" : "alpha");
    print_mat(psq_so,nbfso,nbfso,outfile);
    fprintf(outfile,"\n");
  }
    
  free_block(onepdm);
  free_block(opdm_blk);
  free_block(tmp_mat);
}

}} // end namespace psi::mvo

