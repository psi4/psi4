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

namespace psi { namespace ccresponse {

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

void transpert(const char *pert, double **target)
{
  int nso, nmo;
  int i, j;
  double **scratch, **TMP, **X, **target;
  char *name;
  double prefactor, sign;

  nso = moinfo.nso;
  nmo = moinfo.nmo;

  TMP = block_matrix(nao, nao);
  X = block_matrix(nao, nao);
  scratch = block_matrix(nso, nso);

  if(pert == "Mu_X" || pert == "Mu_Y" || pert == "Mu_Z") { 
    prefactor = 1.0; sign = 1.0; 
  }

  if(pert == "L_X" || pert == "L_Y" || pert == "L_Z") { 
    prefactor = -0.5; sign = 1.0; 
  }

  if(pert == "L*_X" || pert == "L*_Y" || pert == "L*_Z") { 
    prefactor = -0.5; sign = -1.0; 
  }

  if(pert == "P_X" || pert == "P_Y" || pert == "P_Z") { 
    prefactor = 1.0; sign = 1.0; 
  }

  if(pert == "P*_X" || pert == "P*_Y" || pert == "P*_Z") { 
    prefactor = 1.0; sign = -1.0; 
  }

  if(pert == "Q_XX" || pert == "Q_XY" || pert == "Q_XZ" || 
     pert == "Q_YX" || pert == "Q_YY" || pert == "Q_YZ" ||
     pert == "Q_ZX" || pert == "Q_ZY" || pert == "Q_ZZ") { 
    prefactor = 1.0; sign = 1.0; 
  }

  target = block_matrix(nmo,nmo);

  if(pert == "Mu_X") { name = PSIF_AO_MX; moinfo.MU[0] = target; }
  if(pert == "Mu_Y") { name = PSIF_AO_MY; moinfo.MU[1] = target; }
  if(pert == "Mu_Z") { name = PSIF_AO_MZ; moinfo.MU[2] = target; }

  if(pert == "L_X" || pert == "L*_X") {
    name = PSIF_AO_LX; moinfo.L[0] = target; 
  }
  if(pert == "L_Y" || pert == "L*_Y") {
    name = PSIF_AO_LY; moinfo.L[1] = target; 
  }
  if(pert == "L_Z" || pert == "L*_Z") {
    name = PSIF_AO_LZ; moinfo.L[2] = target;
  }

  if(pert == "P_X" || pert == "P*_X") {
    name = PSIF_AO_NablaX; moinfo.P[0] = target; 
  }
  if(pert == "P_Y" || pert == "P*_Y") {
    name = PSIF_AO_NablaY; moinfo.P[1] = target;
  }
  if(pert == "P_Z" || pert == "P*_Z") {
    name = PSIF_AO_NablaZ; moinfo.P[2] = target;
  }

  if(pert == "Q_XX") { name = PSIF_AO_TXX; moinfo.Q[0][0] = target; }
  if(pert == "Q_XY") { name = PSIF_AO_TXY; moinfo.Q[0][1] = target; }
  if(pert == "Q_XZ") { name = PSIF_AO_TXZ; moinfo.Q[0][2] = target; }
  if(pert == "Q_YX") { name = PSIF_AO_TXY; moinfo.Q[1][0] = target; }
  if(pert == "Q_YY") { name = PSIF_AO_TYY; moinfo.Q[1][1] = target; }
  if(pert == "Q_YZ") { name = PSIF_AO_TYZ; moinfo.Q[1][2] = target; }
  if(pert == "Q_ZX") { name = PSIF_AO_TXZ; moinfo.Q[2][0] = target; }
  if(pert == "Q_ZY") { name = PSIF_AO_TYZ; moinfo.Q[2][1] = target; }
  if(pert == "Q_ZZ") { name = PSIF_AO_TZZ; moinfo.Q[2][2] = target; }

  MintsHelper mints;
  vector<SharedMatrix> dipole = mints.so_dipole();
  double **dipole_x = dipole[0]->to_block_matrix();

  for(i=0,ij=0; i < nso; i++)
    for(j=0; j < nso; j++,ij++) {
      TMP[i][j] = prefactor * sign * scratch[ij];
    }

  free_block(dipole_x);

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(moinfo.scf[0][0]),nmo,
	  0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(moinfo.scf[0][0]),nmo,&(X[0][0]),nao,
	  0,&(target[0][0]),nmo);

  free(scratch);
  free_block(TMP);
  free_block(X);
}

}} // namespace psi::ccresponse
