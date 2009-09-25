/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#include "intco.h"
#include "params.h"
#include "displacements.h"
#include "bmat.h"
//#include "transform.h"

#define EXTERN
#include "globals.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>
#include <libpsio/psio.h>
#include <psifiles.h>

using namespace psi::intder;

namespace psi { namespace intder {
extern InternalCoordinates gIntCo;
extern Params gParams;
//extern Transform gTransform;
}}

BMat::BMat()
{
  BMatrix = NULL;
  BMatSave = NULL;
  SVectArray = NULL;
  AMatrix = NULL;
}

BMat::~BMat()
{
  if (BMatrix) {
    free_matrix(BMatrix, gParams.nintco);  //Frees number of rows. Hope this is right!
    free_matrix(BMatSave, gParams.nintco);
  }
  if (SVectArray) {
    delete[] SVectArray;
  }
  if (AMatrix)
    free_matrix(AMatrix, gParams.ncartesians);
}

void BMat::init()
{
  fprintf(outfile, "\nNumber of nintcos %i number of cartesians %i\n", gParams.nintco, gParams.ncartesians);
  BMatrix = init_matrix(gParams.nintco, gParams.ncartesians);
  BMatSave = init_matrix(gParams.nintco, gParams.ncartesians);
  SVectArray = new double[gParams.nintco];
  AMatrix = init_matrix(gParams.ncartesians, gParams.ncartesians);

  fprintf(outfile, "BMatrix created with (%d, %d) dimensions\n", gParams.nintco, gParams.ncartesians);
  fprintf(outfile, "SVectArray created with (%d) dimensions\n", gParams.nintco);
}

void BMat::make_BMat()
{
  int i, j, k;
  double *BR;
  BR = init_array(gParams.ncartesians);
  
  for(i = 0; i < gParams.nintco; i++) {
    BRow(BMatrix[i], 0.0, i, SVectArray);
    //gIntCo.setValue(SVectArray[i]);
    //for(j = 0; j < gParams.ncartesians; j++)
    // BMatrix[i][j] = BR[j]; 
  }
}

void BMat::BRow(double *BMatrow, double disp, int r, double *BValue)
{
  
  int i, j, nIntcoAtoms;
  int atomA, atomB, atomC, atomD;
  
  double *s1, *s2, *s3, *s4;
  s1 = new double[3];
  s2 = new double[3];
  s3 = new double[3];
  s4 = new double[3];
  
  nIntcoAtoms = gIntCo.intcoSwitch(disp, r, &atomA, &atomB, &atomC, &atomD, s1, s2, s3, s4, &BValue[r]);
  
  if(nIntcoAtoms == 2) {
    if(gIntCo.vectorInternals[r]->getType() == STRE) 
      StoreElement(BMatrow, atomA, atomB, s1, s2);
    fprintf(outfile, "Returned value of nintco %i is %lf, atomA = %i, atomB = %i\n", r, BValue[r], atomA+1, atomB+1); 
  }
  
  if(nIntcoAtoms == 3) {
    StoreElement(BMatrow, atomA, atomB, atomC, s1, s2, s3);
    fprintf(outfile, "Returned value of nintco %i is %lf, atomA = %i, atomB = %i, atomC = %i\n", r, BValue[r], atomA+1, 
	    atomB+1, atomC+1); 
  }
  
  if(nIntcoAtoms == 4) {
    StoreElement(BMatrow, atomA, atomB, atomC, atomD, s1, s2, s3, s4);
    fprintf(outfile, "Returned value of nintco %i is %lf, atomA = %i, atomB = %i, atomC = %i, atomD = %i\n", r, BValue[r], 
	    atomA+1, atomB+1, atomC+1, atomD+1); 
  }
}

//For two-element intco's
void BMat::StoreElement(double *row, int aa, int ab, double *s1, double *s2)
{
  int j;
  
  for(j = 0; j < 3; j++) {
    row[(aa*3)+j] = s1[j];
    
    row[(ab*3)+j] = s2[j];
  }
}

//For three-element intco's
void BMat::StoreElement(double *row, int aa, int ab, int ac, double *s1, double *s2, double *s3)
{
  int j;

  for(j = 0; j < 3; j++) {
    row[(aa*3)+j] = s1[j];

    row[(ab*3)+j] = s2[j];

    row[(ac*3)+j] = s3[j];
  }
}

//For four-element intco's (lin1)
void BMat::StoreElement(double *row, int aa, int ab, int ac, int ad, double *s1, double *s2, double *s3)
{
  int j;

  for(j = 0; j < 3; j++) {
    row[(aa*3)+j] = s1[j];

    row[(ab*3)+j] = s2[j];

    row[(ac*3)+j] = s3[j];

    /* We might have some work to do to make this work. Which one is dummy atom? we have 4 atoms
    and only 3 s-vectors. Is fourth s-vector zero? Or can we just eliminate this and pass into
    3-element StoreElement?  */

    //BMatrix[i][(ad*3)+j] = s1[j];
    //BMatrix[i][(ad*3)+j] = s2[j];
    //BMatrix[i][(ad*3)+j] = s3[j];
  }
}

//For four-element intco's
void BMat::StoreElement(double *row, int aa, int ab, int ac, int ad, double *s1, double *s2, double *s3, double*s4)
{
  int j;

  for(j = 0; j < 3; j++) {
    row[(aa*3)+j] = s1[j];

    row[(ab*3)+j] = s2[j];

    row[(ac*3)+j] = s3[j];

    row[(ad*3)+j] = s4[j];
  }
}

/*Hey compadre, this is commented out in the intder.cc driver. It doesn't give the correct determinant
  So something's wrong with the code. */
void BMat::invert_BMat()
{
  int i,j,k;
  int ntri = (gParams.nintco * (gParams.nintco + 1)) / 2;
  
  double tol_check = 0.0;

  double **BuBT;   //We're calling BuBT the equivalent of D in intder2000.f
  double **invBuBT;
  double *BuBTevals;
  double **BuBTevects;
  double **tol_matrix;
  double determinant;
  double *Barray;

  BuBT = block_matrix(gParams.nintco,gParams.nintco);
  invBuBT = block_matrix(gParams.nintco, gParams.nintco);
  BuBTevals = new double[gParams.nintco];
  BuBTevects = block_matrix(gParams.nintco,gParams.nintco);
  tol_matrix = block_matrix(gParams.nintco,gParams.nintco);
  Barray = init_array(ntri);
   
  for(k = 0; k < gParams.ncartesians; k++) 
    fprintf(outfile, "\nChecking for correct mass of atom %i = %lf", (k / 3) + 1, gParams.mass_array[k / 3]);
  
  for(j = 0; j < gParams.nintco; j++) {
    for(i = 0; i < gParams.nintco; i++) {
      for(k = 0; k < gParams.ncartesians; k++) {
	BuBT[i][j] += BMatrix[i][k]*BMatrix[j][k] / gParams.mass_array[k / 3];
      }
    }
  }

  fprintf(outfile, "\n\nB*BT Matrix for Internal Coordinates\n");
  print_mat(BuBT, gParams.nintco, gParams.nintco, outfile);

  //Here, INTDER2000.f makes the second half of the BuBT matrix a unit symmetric matrix, and passes both to FLIN
  //Instead of FLIN, we're going to make it a symmetric square matrix (BuB^T) and use Rollin's Symm_matrix_invert

  invBuBT = symm_matrix_invert(BuBT, gParams.nintco,1,0);
  
  fprintf(outfile, "\nPrinting inverse of BuBt:\n");
  print_mat(invBuBT, gParams.nintco, gParams.nintco, outfile);

  for(j = 0; j < gParams.nintco; j++) {
    for(i = 0; i < gParams.ncartesians; i++) {
      for(k = 0; k < gParams.nintco; k++) {
	AMatrix[i][j] += BMatrix[k][i] * invBuBT[j][k] / gParams.mass_array[i / 3];
      }
    }
  }

  if(gParams.print_array[1] > 3) {
    fprintf(outfile, "\n A Matrix for Internal Coordinates\n");
    print_mat(AMatrix,  gParams.ncartesians, gParams.nintco, outfile);
  }
  
  if(gParams.matrixTest == 0)
    return; //Memory leak here
  else if(gParams.matrixTest == 1) {
    for(j = 0; j < gParams.nintco; j++) {
      for(i = 0; i < gParams.nintco; i++) {
	for(k = 0; k < gParams.ncartesians; k++) {
	  tol_check += BMatrix[i][k] * AMatrix[k][j];
	}
	tol_matrix[i][j] == tol_check;
	if(i == j)
	  tol_check--;
	if(fabs(tol_check) > gParams.invtol) {
	  fprintf(outfile, "\n\nB Matrix inversion for given tolerance was unsuccessful!\n");
	  print_mat(BuBT, gParams.nintco, gParams.nintco, outfile);
	  exit(2); 
	}
      }
    }
  }
  else if(gParams.matrixTest > 1) {

    for(j = 0; j < gParams.nintco; j++) {
      for(i = 0; i < gParams.nintco; i++) {
	for(k = 0; k < gParams.ncartesians; k++) {
	  tol_matrix[i][j] += BMatrix[i][k]*BMatrix[j][k] / gParams.mass_array[k / 3];
	}
      }
    }
  
    fprintf(outfile, "\n\ntol_matrix debug for internal coordinates\n");
    print_mat(tol_matrix, gParams.nintco, gParams.nintco, outfile);
    
    for(i = 0; i < gParams.nintco; i++)
      for(j = 0; j < gParams.nintco; j++) 
	if(i != j)
	  tol_matrix[i][j] /= (sqrt(tol_matrix[i][i] * tol_matrix[j][j]));
    
    for(i = 0; i < gParams.nintco; i++)
      tol_matrix[i][i] = 1.0;
    
    fprintf(outfile, "\n\nNormalized overlap matrix for internal coordinates\n");
    print_mat(tol_matrix, gParams.nintco, gParams.nintco, outfile);
    
    k = 0;
    for(i = 0; i < gParams.nintco; i++)
      for(j = 0; j < i; j++) {
	Barray[k] = tol_matrix[i][j];
	k++;
      }
    
    rsp(gParams.nintco, gParams.nintco, ntri, Barray, BuBTevals, 1, BuBTevects, 1.0E-10);
    
    fprintf(outfile, "\n\nEigenvectors and eigenvalues of normalized overlap matrix\n");
    
    print_mat(BuBTevects, gParams.nintco, gParams.nintco, outfile);
    
    determinant = 1.0;
    fprintf(outfile, "\n\n");
    for(i = 0; i < gParams.nintco; i++) {
      fprintf(outfile, "%lf\t", BuBTevals[i]);
      determinant *= BuBTevals[i];
    }
    
    fprintf(outfile, "\nDeterminant of overlap matrix = %lf\n", determinant);
  }
  
  free_block(BuBT);
  free_block(invBuBT);
  free_block(BuBTevects);
  free_block(tol_matrix);
  delete[] BuBTevals;
}

void BMat::print_BMat()
{
  int x;
  
  fprintf(outfile, "\nValues of simple internal coordinates (ANG. OR DEG.): \n");
  for(x = 0; x < gParams.nintco; x++)
    if(gIntCo.vectorInternals[x]->getType() == STRE) 
      fprintf(outfile, "\n%i = %10.7lf", x, gParams.UGF[x][0]);
    else {
      gParams.UGF[x][0] = (gParams.UGF[x][0] * 180.0 / _pi);
      fprintf(outfile, "\n%i = %10.7lf", x, gParams.UGF[x][0]);
    }
  
  fprintf(outfile, "\n\nB Matrix for Internal Coordinates:\n");
  print_mat(BMatSave, gParams.nintco,  gParams.ncartesians, outfile) ;
}

