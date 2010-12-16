/*! \file v3d.cc
    \ingroup OPT10
    \brief v3d functions
*/

#include "v3d.h"

#define EXTERN
#include "globals.h"

namespace opt { namespace v3d {

// Compute angle in radians A-B-C (between vector B->A and vector B->C)
// if points are absurdly close or far apart, returns false
bool v3d_angle(const double *A, const double *B, const double *C, double & phi) {
  double dotprod, eBA[3], eBC[3];

  // eBA
  if (! v3d_eAB(B, A, eBA) ) {
fprintf(outfile, "could not normalize eBA, B:");
for (int i=0; i<3; ++i) fprintf(outfile,"%15.10lf", B[i]);
fprintf(outfile,"\n A:");
for (int i=0; i<3; ++i) fprintf(outfile,"%15.10lf", A[i]);
    return false;
  }

  // eBC
  if (! v3d_eAB(B, C, eBC) ) {
fprintf(outfile, "could not normalize eBC, B:");
for (int i=0; i<3; ++i) fprintf(outfile,"%15.10lf", B[i]);
fprintf(outfile,"\n A:");
for (int i=0; i<3; ++i) fprintf(outfile,"%15.10lf", A[i]);
    return false;
  }

  dotprod = v3d_dot(eBA,eBC);

  if (dotprod >= -1.0 && dotprod <= 1.0)
    phi = acos(dotprod);
  else if (dotprod > 1.0)
    phi = 0.0;
  else
    phi = acos(-1);

  return true;
}

// returns false if bond angles are too large for good torsion definition
bool v3d_tors(const double *A, const double *B, const double *C, const double *D,
  double & tau) {
  double tval, phi_123, phi_234, ulim, llim;
  double eAB[3], eBC[3], eCD[3], tmp[3], tmp2[3];

  tau = 0.0;

  // form e unit vectors 
  if ( !v3d_eAB(A,B,eAB) || !v3d_eAB(B,C,eBC) || !v3d_eAB(C,D,eCD) )
    throw ("v3d_tors: distances are not reasonably normalized for e vectors.");

  //printf("v3d_eAB : %15.10lf %15.10lf %15.10lf \n", eAB[0], eAB[1], eAB[2]);
  //printf("v3d_eBC : %15.10lf %15.10lf %15.10lf \n", eBC[0], eBC[1], eBC[2]);
  //printf("v3d_eCD : %15.10lf %15.10lf %15.10lf \n", eCD[0], eCD[1], eCD[2]);

  // compute bond angles
  if ( !v3d_angle(A, B, C, phi_123) || !v3d_angle(B, C, D, phi_234) )
    throw ("v3d_tors: distances are not reasonably normalized for angles.");

  //printf("v3d_tors : phi123 = %15.10lf\n", phi_123);
  //printf("v3d_tors : phi234 = %15.10lf\n", phi_234);

  // check bond angles; don't allow torsions with angles between 1.8 degrees of singularity
  const double pi = acos(-1);
  llim = Opt_params.tors_angle_lim * pi;
  ulim = (1-Opt_params.tors_angle_lim) * pi;
  if (phi_123 < llim || phi_123 > ulim)
    return false;
  else if (phi_234 < llim || phi_234 > ulim)
    return false;

  v3d_cross_product(eAB,eBC,tmp);
  v3d_cross_product(eBC,eCD,tmp2);
  tval = v3d_dot(tmp,tmp2) / (sin(phi_123) * sin(phi_234) );

  if (tval >= 1.0 - Opt_params.tors_cos_tol)
    tau = 0.0;
  else if (tval <= -1.0 + Opt_params.tors_cos_tol)
    tau = pi;
  else
    tau = acos(tval);

  // determine sign of torsion ; this convention matches Wilson, Decius and Cross
  v3d_cross_product(eBC,eCD,tmp);
  tval = v3d_dot(eAB, tmp);
  if (tval < 0) tau *= -1;

  return true;
}


}}

