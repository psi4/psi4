/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <libmints/mints.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

using namespace std;

namespace psi { namespace ccdensity {

void transL(double sign)
{
  int nmo, nso;
  double **scf_qt, **X;
  double **LX_MO, **LY_MO, **LZ_MO;
  double **LX_SO, **LY_SO, **LZ_SO;

  nso = moinfo.nso;
  nmo = moinfo.nmo;
  scf_qt = moinfo.scf_qt;

  /*** Transform the SO nabla integrals to the MO basis ***/
  MintsHelper mints(Process::environment.options, 0);
  vector<SharedMatrix> dipole = mints.so_angular_momentum();
  for(int i=0; i < 3; i++) dipole[i]->scale(-0.5 * sign);
  LX_SO = dipole[0]->to_block_matrix();
  LY_SO = dipole[1]->to_block_matrix();
  LZ_SO = dipole[2]->to_block_matrix();

  X = block_matrix(nmo,nso); /* just a temporary matrix */

  LX_MO=block_matrix(nmo,nmo);
  LY_MO=block_matrix(nmo,nmo);
  LZ_MO=block_matrix(nmo,nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(LX_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(LX_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(LY_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(LY_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(LZ_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(LZ_MO[0][0]),nmo);

  free_block(X);

  moinfo.L = (double ***) malloc(3 * sizeof(double **));
  moinfo.L[0] = LX_MO;
  moinfo.L[1] = LY_MO;
  moinfo.L[2] = LZ_MO;

  free_block(LX_SO); 
  free_block(LY_SO); 
  free_block(LZ_SO);

  return;
}

}} // namespace psi::ccdensity
