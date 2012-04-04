/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libmints/mints.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

using namespace std;

namespace psi { namespace ccresponse {

void sort_pert(const char *pert, double **pertints, int irrep);

/* preppert(): Prepare DPD structures for all currently known one-electron 
** property integrals in the MO basis.
**
** -TDC, 6/11
*/

void preppert()
{
  int i, j, ij;

  char **cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");
  char lbl[32];
 
  MintsHelper mints(Process::environment.options, 0);
  vector<SharedMatrix> dipole = mints.so_dipole();
  vector<SharedMatrix> nabla = mints.so_nabla();
  vector<SharedMatrix> angmom = mints.so_angular_momentum();
  vector<SharedMatrix> trquad = mints.so_traceless_quadrupole();

  int nso = moinfo.nso;
  int nmo = moinfo.nmo;

  double **TMP2 = block_matrix(nso,nso);

  // Electric dipole integrals
  for(i=0; i < 3; i++) {
    double **TMP1 = dipole[i]->to_block_matrix();
    C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
    C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP1[0],nmo);
    moinfo.MU[i] = TMP1;
    sprintf(lbl, "Mu_%1s", cartcomp[i]);
    sort_pert(lbl, moinfo.MU[i], moinfo.mu_irreps[i]);
  }

  // Velocity-gauge electric dipole integrals
  for(i=0; i < 3; i++) {
    double **TMP1 = nabla[i]->to_block_matrix();
    C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
    C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP1[0],nmo);
    moinfo.P[i] = TMP1;
    sprintf(lbl, "P_%1s", cartcomp[i]);
    sort_pert(lbl, moinfo.P[i], moinfo.mu_irreps[i]);
  }

  // Complex conjugate of velocity-gauge electric dipole integrals
  for(i=0; i < 3; i++) nabla[i]->scale(-1.0);
  for(i=0; i < 3; i++) {
    double **TMP1 = nabla[i]->to_block_matrix();
    C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
    C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP1[0],nmo);
    moinfo.Pcc[i] = TMP1;
    sprintf(lbl, "P*_%1s", cartcomp[i]);
    sort_pert(lbl, moinfo.Pcc[i], moinfo.mu_irreps[i]);
  }

  // Magnetic dipole integrals (these require a -1/2 prefactor)
  for(i=0; i < 3; i++) {
    angmom[i]->scale(-0.5);
    double **TMP1 = angmom[i]->to_block_matrix();
    sprintf(lbl, "L_%1s", cartcomp[i]);
  //  fprintf(outfile, "%s Angular Momentum Integrals (SO)\n",lbl);
//    mat_print(TMP1,nmo, nmo, outfile);
    C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
    C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP1[0],nmo);
    moinfo.L[i] = TMP1;
    sort_pert(lbl, moinfo.L[i], moinfo.l_irreps[i]);
  }

  // Complex conjugate of magnetic dipole integrals
  for(i=0; i < 3; i++) angmom[i]->scale(-1.0);
  for(i=0; i < 3; i++) {
    double **TMP1 = angmom[i]->to_block_matrix();
    C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
    C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP1[0],nmo);
    moinfo.Lcc[i] = TMP1;
    sprintf(lbl, "L*_%1s", cartcomp[i]);
    sort_pert(lbl, moinfo.Lcc[i], moinfo.l_irreps[i]);
  }

  // Traceless quadrupole integrals
  for(i=0,ij=0; i < 3; i++) {
    for(j=i; j < 3; j++,ij++) {
      double **TMP1 = trquad[ij]->to_block_matrix();
      C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
      C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP1[0],nmo);
      moinfo.Q[i][j] = TMP1;
      sprintf(lbl, "Q_%1s%1s", cartcomp[i], cartcomp[j]);
      sort_pert(lbl, moinfo.Q[i][j], moinfo.mu_irreps[i]^moinfo.mu_irreps[j]);
      if(i!=j) {
        moinfo.Q[j][i] = TMP1;
        sprintf(lbl, "Q_%1s%1s", cartcomp[j], cartcomp[i]);
        sort_pert(lbl, moinfo.Q[j][i], moinfo.mu_irreps[j]^moinfo.mu_irreps[i]);
      }
    }
  }

  free_block(TMP2);
}

}} // namespace psi::ccresponse
