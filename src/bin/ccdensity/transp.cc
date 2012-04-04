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

void transp(double sign)
{
  int nmo, nso;
  double **scf_qt, **X;
  double **NX_MO, **NY_MO, **NZ_MO;
  double **NX_SO, **NY_SO, **NZ_SO;

  nso = moinfo.nso;
  nmo = moinfo.nmo;
  scf_qt = moinfo.scf_qt;

  /*** Transform the SO nabla integrals to the MO basis ***/
  MintsHelper mints(Process::environment.options, 0);
  vector<SharedMatrix> dipole = mints.so_nabla();
  for(int i=0; i < 3; i++) dipole[i]->scale(-1.0 * sign);
  NX_SO = dipole[0]->to_block_matrix();
  NY_SO = dipole[1]->to_block_matrix();
  NZ_SO = dipole[2]->to_block_matrix();

  X = block_matrix(nmo,nso); /* just a temporary matrix */

  NX_MO=block_matrix(nmo,nmo);
  NY_MO=block_matrix(nmo,nmo);
  NZ_MO=block_matrix(nmo,nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(NX_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(NX_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(NY_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(NY_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(NZ_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(NZ_MO[0][0]),nmo);

  free_block(X);

  moinfo.nabla = (double ***) malloc(3 * sizeof(double **));
  moinfo.nabla[0] = NX_MO;
  moinfo.nabla[1] = NY_MO;
  moinfo.nabla[2] = NZ_MO;

  free_block(NX_SO); 
  free_block(NY_SO); 
  free_block(NZ_SO);

  return;
}

}} // namespace psi::ccdensity
