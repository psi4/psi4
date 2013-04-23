/*! \file
    \ingroup RESPONSE
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include <psifiles.h>
#include "globals.h"

namespace psi { namespace response {

void trans_mu(void)
{
  int stat, nao, noei_ao, nso, nmo;
  int i, j, ij;
  double *scratch, **TMP, **X;
  double **MUX, **MUY, **MUZ;  /* MO-basis dipole integrals */

  nao = moinfo.nao;
  nso = moinfo.nso;
  nmo = moinfo.nmo;
  noei_ao = moinfo.noei_ao;

  /**** Transform the dipole integrals to the MO basis ****/

  MUX = block_matrix(nmo,nmo);
  MUY = block_matrix(nmo,nmo);
  MUZ = block_matrix(nmo,nmo);

  TMP = block_matrix(nao, nao);
  X = block_matrix(nao, nao);
  scratch = init_array(noei_ao);

  stat = iwl_rdone(PSIF_OEI, PSIF_AO_MX, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++, ij++) {
      TMP[i][j] = TMP[j][i] = scratch[ij];
    }

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(moinfo.usotao[0][0]),nao,
          0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(moinfo.usotao[0][0]),nao,&(X[0][0]),nao,
          0,&(TMP[0][0]),nao);

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(moinfo.scf[0][0]),nmo,
          0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(moinfo.scf[0][0]),nmo,&(X[0][0]),nao,
          0,&(MUX[0][0]),nmo);

  zero_arr(scratch,noei_ao);

  stat = iwl_rdone(PSIF_OEI, PSIF_AO_MY, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++, ij++) {
      TMP[i][j] = TMP[j][i] = scratch[ij];
    }

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(moinfo.usotao[0][0]),nao,
          0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(moinfo.usotao[0][0]),nao,&(X[0][0]),nao,
          0,&(TMP[0][0]),nao);

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(moinfo.scf[0][0]),nmo,
          0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(moinfo.scf[0][0]),nmo,&(X[0][0]),nao,
          0,&(MUY[0][0]),nmo);

  zero_arr(scratch,noei_ao);

  stat = iwl_rdone(PSIF_OEI, PSIF_AO_MZ, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++, ij++) {
      TMP[i][j] = TMP[j][i] = scratch[ij];
    }

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(moinfo.usotao[0][0]),nao,
          0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(moinfo.usotao[0][0]),nao,&(X[0][0]),nao,
          0,&(TMP[0][0]),nao);

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(moinfo.scf[0][0]),nmo,
          0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(moinfo.scf[0][0]),nmo,&(X[0][0]),nao,
          0,&(MUZ[0][0]),nmo);

  free_block(TMP);
  free_block(X);
  free(scratch);

  moinfo.MUX = MUX;
  moinfo.MUY = MUY;
  moinfo.MUZ = MUZ;

  /*
  fprintf(outfile, "MO-Basis MuX Integrals:\n");
  mat_print(moinfo.MUX, nmo, nmo, outfile);
  fprintf(outfile, "MO-Basis MuY Integrals:\n");
  mat_print(moinfo.MUY, nmo, nmo, outfile);
  fprintf(outfile, "MO-Basis MuZ Integrals:\n");
  mat_print(moinfo.MUZ, nmo, nmo, outfile);
  */

}

}} // namespace psi::response
