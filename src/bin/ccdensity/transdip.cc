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

void transdip(void)
{
  int nmo, nso;
  double **scf_qt, **X;
  double **MUX_MO, **MUY_MO, **MUZ_MO;
  double **MUX_SO, **MUY_SO, **MUZ_SO;

  nso = moinfo.nso;
  nmo = moinfo.nmo;
  scf_qt = moinfo.scf_qt;

  /*** Transform the SO dipole integrals to the MO basis ***/
  MintsHelper mints;
  vector<SharedMatrix> dipole = mints.so_dipole();
  MUX_SO = dipole[0]->to_block_matrix();
  MUY_SO = dipole[1]->to_block_matrix();
  MUZ_SO = dipole[2]->to_block_matrix();

  X = block_matrix(nmo,nso); /* just a temporary matrix */

  MUX_MO=block_matrix(nmo,nmo);
  MUY_MO=block_matrix(nmo,nmo);
  MUZ_MO=block_matrix(nmo,nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(MUX_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(MUX_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(MUY_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(MUY_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(MUZ_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(MUZ_MO[0][0]),nmo);

  free_block(X);

  moinfo.dip = (double ***) malloc(3 * sizeof(double **));
  moinfo.dip[0] = MUX_MO;
  moinfo.dip[1] = MUY_MO;
  moinfo.dip[2] = MUZ_MO;

  free_block(MUX_SO); 
  free_block(MUY_SO); 
  free_block(MUZ_SO);

  return;
}

}} // namespace psi::ccdensity
