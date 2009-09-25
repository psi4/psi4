/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#include "intco.h"
#include "transform.h"
#include "params.h"
#include "atom.h"
#include "bmat.h"

#define EXTERN
#include "globals.h"
#include "3dmatrix.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>

using namespace psi::intder;

namespace psi { namespace intder {
extern InternalCoordinates gIntCo;
extern Params gParams;
extern BMat gBMat;
}}
extern void open_PSIF();
extern void close_PSIF();
  
extern void secondDerivativeOut(int, double **);
extern void secondDerivativeIn(int, double **);

Transform::Transform()
{
/* Numerical testing stuff that is untested 

  if(gParams.derlvl + gParams.atEquilibrium > 1) {
    gParams.intco2ndDer = NULL;
    SR2Matrix = NULL;
    gParams.intco2ndDer = init_matrix(gParams.ncartesians, gParams.ncartesians);
    SR2Matrix = init_matrix(gParams.ncartesians, gParams.ncartesians);
  }
  if(gParams.derlvl + gParams.atEquilibrium > 2) {
    gParams.intco3rdDer = new C3DMatrix(gParams.ncartesians, gParams.ncartesians, gParams.ncartesians);
    gParams.intco3rdDer->array = init_array(gParams.intco3rdDer->Size(gParams.ncartesians,gParams.ncartesians,gParams.ncartesians));
    SR3Matrix = new C3DMatrix(gParams.ncartesians, gParams.ncartesians, gParams.ncartesians);
    SR3Matrix->array = init_array(SR3Matrix->Size(gParams.ncartesians,gParams.ncartesians,gParams.ncartesians));
  }
  if(gParams.derlvl + gParams.atEquilibrium > 3) {
    gParams.intco4thDer = new C4DMatrix(gParams.ncartesians, gParams.ncartesians, gParams.ncartesians, gParams.ncartesians);
    gParams.intco4thDer->array = init_array(gParams.intco4thDer->Size(gParams.ncartesians,gParams.ncartesians,gParams.ncartesians, gParams.ncartesians));
    gParams.intco4thDerBlock = new C3DMatrix(gParams.ncartesians, gParams.ncartesians, gParams.ncartesians);
    gParams.intco4thDerBlock->array = init_array(gParams.intco4thDerBlock->Size(gParams.ncartesians,gParams.ncartesians, gParams.ncartesians));
    SR4Matrix = new C4DMatrix(gParams.ncartesians, gParams.ncartesians, gParams.ncartesians, gParams.ncartesians);
    SR4Matrix->array = init_array(SR4Matrix->Size(gParams.ncartesians,gParams.ncartesians,gParams.ncartesians, gParams.ncartesians));
  }  */

}

Transform::~Transform()
{
/*  if(gParams.derlvl + gParams.atEquilibrium > 1) {
    //~SecondDerivative();
  }
  if(gParams.derlvl + gParams.atEquilibrium > 2) {
    //~ThirdDerivative();
  }
  if(gParams.derlvl + gParams.atEquilibrium > 3) {
    //FourthDerDelete();
  } 
*/
}

/* linear_transform -- Performs a linear transformation from internal to cartesians *
 * or vice versa                                                                    *
 * Passed in variables are dependent on whether or not atomic masses are            *
 * included, i.e. if transformType is greater than or less than zero                *
 *                                                                                  *
 * if transformType > 0, coord1 = gParams.ncartesians, coord2 = gParams.nintco      *
 * transformMatrix = BMatrix                                                        *

 * if transformType < 0, coord1 = gParams.nintco, coord2 = gParams.ncartesians      *
 * transformMatrix = AMatrix                                                        *
 * transformType < 0 is not yet implemented!                                        */

void Transform::der_transform(int coord1, int coord2, double **transformMatrix)
{
  int i,j,m;
  double cf1, cf2, cf3, cf4; //conversion factors, which will be appropriate cf's or equal to 1
  double *temp_mass_array;

  temp_mass_array = init_array(gParams.nintco);

  if(abs(gParams.transformType != 3)) {
    cf1 = _hartree2J / (_bohr2angstroms * 1.0E-18);
    cf2 = cf1 / _bohr2angstroms;
    cf3 = cf2 / _bohr2angstroms;
    cf4 = cf3 / _bohr2angstroms;
  }
  else {
    cf1 = 1.0;
    cf2 = 1.0;
    cf3 = 1.0;
    cf4 = 1.0;
  }
  
  if(gParams.transformType >= 1) {
    cf1 = 1.0;
    cf2 = 1.0;
    cf3 = 1.0;
    cf4 = 1.0;
  }  //this seems dumb
  
  if(gParams.atEquilibrium) {
    if(gParams.transformType >= 1)
      gParams.FConst1 = init_array(coord1);
      }
  for(i = 0; i < coord2; i++) {
    temp_mass_array[i] = gParams.mass_array[i] * cf1;
    for(m = 0; m < coord1; m++)
      gParams.FConst1[m] += temp_mass_array[i] * gBMat.AMatrix[i][m];
  }
}

/*
SecondDerivative::SecondDerivative()
{
  H11 = init_matrix(3,3);
  H21 = init_matrix(3,3);
  H22 = init_matrix(3,3);
  H31 = init_matrix(3,3);
  H32 = init_matrix(3,3);
  H33 = init_matrix(3,3);
  H41 = init_matrix(3,3);
  H42 = init_matrix(3,3);
  H43 = init_matrix(3,3);
  H44 = init_matrix(3,3);
}

SecondDerivative::~SecondDerivative()
{
  free_matrix(H11,3);
  free_matrix(H21,3);
  free_matrix(H22,3);
  free_matrix(H31,3);
  free_matrix(H32,3);
  free_matrix(H33,3);
  free_matrix(H41,3);
  free_matrix(H42,3);
  free_matrix(H43,3);
  free_matrix(H44,3);

  free_matrix(SR2Matrix, gParams.ncartesians);
}

ThirdDerivative::ThirdDerivative()
{
  H111 = new C3DMatrix(3,3,3);
  H112 = new C3DMatrix(3,3,3); 
  H221 = new C3DMatrix(3,3,3); 
  H222 = new C3DMatrix(3,3,3); 
  H113 = new C3DMatrix(3,3,3); 
  H123 = new C3DMatrix(3,3,3); 
  H223 = new C3DMatrix(3,3,3); 
  H331 = new C3DMatrix(3,3,3); 
  H332 = new C3DMatrix(3,3,3); 
  H333 = new C3DMatrix(3,3,3); 
  H411 = new C3DMatrix(3,3,3); 
  H421 = new C3DMatrix(3,3,3); 
  H422 = new C3DMatrix(3,3,3);
  H431 = new C3DMatrix(3,3,3); 
  H432 = new C3DMatrix(3,3,3); 
  H433 = new C3DMatrix(3,3,3); 
  H441 = new C3DMatrix(3,3,3); 
  H442 = new C3DMatrix(3,3,3); 
  H443 = new C3DMatrix(3,3,3); 
  H444 = new C3DMatrix(3,3,3); 

}

ThirdDerivative::~ThirdDerivative()
{
  delete H111;
  delete H112;
  delete H221;
  delete H222;
  delete H113;
  delete H123;
  delete H223;
  delete H331;
  delete H332;
  delete H333;
  delete H411;
  delete H421;
  delete H431;
  delete H432;
  delete H433;
  delete H441;
  delete H442;
  delete H443;
  delete H444;

  delete SR3Matrix;
}

FourthDerivative::FourthDerivative()
{

}

FourthDerivative::~FourthDerivative()
{

}

//This will be MACHX
void SecondDerivative::i_to_c()
{
  int i,j,r;
  double s = 0.0;
  
  XRow(SR2Matrix, 0, r); 
  
  //Does this command ever get called in intder2000.f?!
  if(gParams.intcoInclude[r]) {
    secondDerivativeIn(-r, SR2Matrix);
    secondDerivativeIn(r, gParams.intco2ndDer);
    fprintf(outfile, "Numerical SR(I,J) and X(M,N) matrices used for simple internal coordinate %i\n",r);
  } 
  else
    fprintf(outfile, "SR(I,J) and X(M,N) Matrices set to zero for simple internal coordinate %i\n",r);
  
  secondDerivativeOut(-r, SR2Matrix);
  secondDerivativeOut(r, gParams.intco2ndDer);
  //Again, matrix dimensions are not right. Also need to add lines dealing with symmetric internals
}

//MACHY
void ThirdDerivative::i_to_c()
{
  int i,j,r;
  
  YRow(SR3Matrix, 0, r);
  
  if(gParams.intcoInclude[r]) {
    thirdDerivativeIn(-r, SR3Matrix);
    thirdDerivativeIn(r, gParams.intco3rdDer);
    fprintf(outfile, "Numerical SR(I,J,K) and Y(M,N,P) matrices used for simple internal coordinate %i\n",r);
  } 
  else
    fprintf(outfile, "SR(I,J,K) and Y(M,N,P) Matrices set to zero for simple internal coordinate %i\n",r);
  
  fill3a(gParams.ncartesians, gParams.ncartesians, gParams.intco3rdDer);
  thirdDerivativeOut(-r, SR3Matrix);
  thirdDerivativeOut(r, gParams.intco3rdDer);
}

void FourthDerivative::i_to_c()
{


} */

