/*! \file
    \ingroup response
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace response {

/* transpert(): Transform various one-electron property integrals from
** the AO to the MO basis.  In some cases, we must also add
** appropriate prefactors and signs.  The only argument is a
** character string indicating the type of integral we want: "Mu",
** "L", "L*", "P", or "P*".  The cints code produces only lower
** triangles, so we must unpack the integrals and keep up with
** symmetric vs. antisymmetric cases.
**
** Notes on specific integrals (all produced by "cints --oeprop") used
** here:
**
** (1) Mu: Length-gauge electric dipole moment integrals = -r.  These
** already include the electronic charge, and they are symmetric wrt
** index permutation.

** (2) L: Magnetic dipole integrals = -i/2 (r x p).  The
** input integrals are angular momentum integrals (r x p) = - i (r x Del),
** but we multiply these by -0.5 to account for both the sign of the
** electronic charge and the definition of the magnetic dipole.  These
** integrals are antisymmetric wrt index permutation.  Use "L*" to select
** the complex conjugate of the operator (i.e., multiply by -1).
**
** (3) P: Velocity-gauge electric dipole moment integrals = -p.  The input
** integrals are +del integrals, but the definition of
** the linear momentum operator involves a -1 and the electron charge is
**  -1, so no special multiplcation is necessary (i.e., p = -i del, so
** e*p = i del).  These integrals are antisymmetric wrt to index permutation.
** Use "P*" to select the complex conjugate of the operator (i.e., multiply
** by -1).
**
** (4) Q: Traceless quadrupole integrals: Q_ab = -1/2 (3 r_a r_b - r^2).
** They are real, already include the electronic charge and the leading 1/2
** and are symmetric wrt index permutation.
**
** NB: The magnetic-dipole and velocity-gauge electric-dipole operators are
** both pure-imaginary, which means that one must take care to account for
** the factor of i included implicity in their definition.  This matters,
** for example, in the computation of optical rotations.  See the notes in
** optrot.cc for specifics.
**
** -TDC, 11/05
** Updated by TDC 4/09
*/

void transpert(const char *pert)
{
  int nao, nso, nmo, noei_ao;
  int alpha;
  int i, j, ij;
  double *scratch, **TMP, **X, **target;
  const char *name;
  double prefactor, anti, sign;

  nao = moinfo.nao;
  nso = moinfo.nso;
  nmo = moinfo.nmo;
  noei_ao = moinfo.noei_ao;

  TMP = block_matrix(nao, nao);
  X = block_matrix(nao, nao);
  scratch = init_array(noei_ao);

  if(!strcmp(pert,"Mu_X") || !strcmp(pert,"Mu_Y") || !strcmp(pert,"Mu_Z")) {
    prefactor = 1.0; anti = 1.0; sign = 1.0;
  }
  if(!strcmp(pert,"L_X") || !strcmp(pert,"L_Y") || !strcmp(pert,"L_Z")) {
    prefactor = -0.5; anti = -1.0; sign = 1.0;
  }
  if(!strcmp(pert,"L*_X") || !strcmp(pert,"L*_Y") || !strcmp(pert,"L*_Z")) {
    prefactor = -0.5; anti = -1.0; sign = -1.0;
  }
  if(!strcmp(pert,"P_X") || !strcmp(pert,"P_Y") || !strcmp(pert,"P_Z")) {
    prefactor = 1.0; anti = -1.0; sign = 1.0;
  }
  if(!strcmp(pert,"P*_X") || !strcmp(pert,"P*_Y") || !strcmp(pert,"P*_Z")) {
    prefactor = 1.0; anti = -1.0; sign = -1.0;
  }
  if(!strcmp(pert,"Q_XX") || !strcmp(pert,"Q_XY") || !strcmp(pert,"Q_XZ") ||
     !strcmp(pert,"Q_YX") || !strcmp(pert,"Q_YY") || !strcmp(pert,"Q_YZ") ||
     !strcmp(pert,"Q_ZX") || !strcmp(pert,"Q_ZY") || !strcmp(pert,"Q_ZZ")) {
    prefactor = 1.0; anti = 1.0; sign = 1.0;
  }

  target = block_matrix(nmo,nmo);

  if(!strcmp(pert,"Mu_X")) { name = PSIF_AO_MX; moinfo.MU[0] = target; }
  if(!strcmp(pert,"Mu_Y")) { name = PSIF_AO_MY; moinfo.MU[1] = target; }
  if(!strcmp(pert,"Mu_Z")) { name = PSIF_AO_MZ; moinfo.MU[2] = target; }

  if(!strcmp(pert,"L_X") || !strcmp(pert, "L*_X")) {
    name = PSIF_AO_LX; moinfo.L[0] = target;
  }
  if(!strcmp(pert,"L_Y") || !strcmp(pert, "L*_Y")) {
    name = PSIF_AO_LY; moinfo.L[1] = target;
  }
  if(!strcmp(pert,"L_Z") || !strcmp(pert, "L*_Z")) {
    name = PSIF_AO_LZ; moinfo.L[2] = target;
  }

  if(!strcmp(pert,"P_X") || !strcmp(pert, "P*_X")) {
    name = PSIF_AO_NablaX; moinfo.P[0] = target;
  }
  if(!strcmp(pert,"P_Y") || !strcmp(pert, "P*_Y")) {
    name = PSIF_AO_NablaY; moinfo.P[1] = target;
  }
  if(!strcmp(pert,"P_Z") || !strcmp(pert, "P*_Z")) {
    name = PSIF_AO_NablaZ; moinfo.P[2] = target;
  }

  if(!strcmp(pert,"Q_XX")) { name = PSIF_AO_TXX; moinfo.Q[0][0] = target; }
  if(!strcmp(pert,"Q_XY")) { name = PSIF_AO_TXY; moinfo.Q[0][1] = target; }
  if(!strcmp(pert,"Q_XZ")) { name = PSIF_AO_TXZ; moinfo.Q[0][2] = target; }
  if(!strcmp(pert,"Q_YX")) { name = PSIF_AO_TXY; moinfo.Q[1][0] = target; }
  if(!strcmp(pert,"Q_YY")) { name = PSIF_AO_TYY; moinfo.Q[1][1] = target; }
  if(!strcmp(pert,"Q_YZ")) { name = PSIF_AO_TYZ; moinfo.Q[1][2] = target; }
  if(!strcmp(pert,"Q_ZX")) { name = PSIF_AO_TXZ; moinfo.Q[2][0] = target; }
  if(!strcmp(pert,"Q_ZY")) { name = PSIF_AO_TYZ; moinfo.Q[2][1] = target; }
  if(!strcmp(pert,"Q_ZZ")) { name = PSIF_AO_TZZ; moinfo.Q[2][2] = target; }

  iwl_rdone(PSIF_OEI, const_cast<char*>(name), scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = prefactor * sign * scratch[ij];
      TMP[j][i] = anti * prefactor * sign * scratch[ij];
    }

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(moinfo.usotao[0][0]),nao,
          0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(moinfo.usotao[0][0]),nao,&(X[0][0]),nao,
          0,&(TMP[0][0]),nao);

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(moinfo.scf[0][0]),nmo,
          0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(moinfo.scf[0][0]),nmo,&(X[0][0]),nao,
          0,&(target[0][0]),nmo);

  zero_arr(scratch,noei_ao);

  free(scratch);
  free_block(TMP);
  free_block(X);
}

}} // namespace psi::response
