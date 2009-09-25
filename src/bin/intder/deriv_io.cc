/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "globals.h"
//#include "transform.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "params.h"
#include "3dmatrix.h"

#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>
#include <libpsio/psio.h>
#include <psifiles.h>

namespace psi { namespace intder {

extern Params gParams;

/* This is INPFKM and it's used when NINV (transformType) > 2, which tells intder   *
 * to read the internal coordinate 1st/2nd/3rd/4th derivatives from the ider.dat    *
 * file                                                                             */   

void readIntCoDerivatives()
{
  /* file11 - first derivatives
     file15 - second derivatives
     file20 - third derivatives
     file24 - fourth derivatives
  */

  int i,j;
  char *buffer;
  double value;

  if(gParams.derlvl)
    gParams.FConst1 = init_array(gParams.nintco);
  if(gParams.derlvl > 1)
    gParams.FConst2 = block_matrix(gParams.nintco, gParams.nintco);
  if(gParams.derlvl > 2)
    gParams.FConst3->Init(gParams.nintco, gParams.nintco, gParams.nintco);
  if(gParams.derlvl > 3)
    gParams.FConst4->Init(gParams.nintco, gParams.nintco, gParams.nintco, gParams.nintco);
  
  if(!gParams.atEquilibrium) {
    ffile(&fp_ider, "ider.dat", 2);
  }	
  
  if(gParams.derlvl == 2) {

    int tempi;
    int tempj;
    int ntri = (gParams.nintco * gParams.nintco + 1);

    ffile(&fp_ider, "ider.dat", 2);
    rewind(fp_ider);
    buffer = new char[MAX_LINE];
        
    for (i=0; i < ntri;i++) {
      fgets(buffer,MAX_LINE,fp_ider);
      sscanf(buffer,"%i %i %lf",&tempi, &tempj, &value);
      if((tempi == 0) || (tempj == 0))
	break;
      gParams.FConst2[tempi - 1][tempj - 1] = value;
    } //The logic here is sketchy and might not work for the general case. 
    delete [] buffer;
    fclose(fp_ider);  
  }

  for(i = 0; i < gParams.nintco; i++)
    for(j = 0; j < i; j++)
      gParams.FConst2[j][i] = gParams.FConst2[i][j];
  
  fprintf(outfile, "\nPrinting the values of the internal coordinate second derivative matrix\n");
  print_mat(gParams.FConst2, gParams.nintco, gParams.nintco, outfile);
}

/*
void SecondDerivative::secondDerivativeOut(int SRFlag, double **SR2Matrix)
{
  long unsigned int bufsize = 0;
  long unsigned int ioff = 0;
  psio_address next;

  bufsize = 2 * (long unsigned int) pow(gParams.ncartesians,2);
  ioff = abs(SRFlag) * bufsize;
  if(SRFlag > 0)
    ioff += (long unsigned int) pow(gParams.ncartesians,2);

  open_PSIF();
  
  next = psio_get_address(PSIO_ZERO, ioff);

  psio_write(PSIF_OPTKING, "Second Derivative / SR Matrix", (char *) (SR2Matrix), (long unsigned int) pow(gParams.ncartesians,2) * sizeof(double), next, &next);
  
  if(SRFlag > 0)
    fprintf(outfile, "X Matrix for Coordinate %i written to File 1 TOC\n", abs(SRFlag));
  else
    fprintf(outfile, "SR Matrix for Coordinate %i written to File 1 TOC\n", abs(SRFlag));
  close_PSIF();
}

void SecondDerivative::secondDerivativeIn(int SRFlag, double **matrix)
{
  long unsigned int bufsize = 0;
  long unsigned int ioff = 0;
  psio_address next;

  bufsize = 2 * (long unsigned int) pow(gParams.ncartesians,2);
  ioff = abs(SRFlag) * bufsize;
  if(SRFlag > 0)
    ioff += (long unsigned int) pow(gParams.ncartesians,2);

  open_PSIF();
  
  next = psio_get_address(PSIO_ZERO, ioff);

  psio_read(PSIF_OPTKING, "Second Derivative / SR Matrix", (char *) (matrix), (long unsigned int) pow(gParams.ncartesians,2) * sizeof(double), next, &next);
  
  if(SRFlag > 0)
    fprintf(outfile, "X Matrix for Coordinate %i read from File 1 TOC\n", abs(SRFlag));
  else
    fprintf(outfile, "SR Matrix for Coordinate %i read from File 1 TOC\n", abs(SRFlag));
  close_PSIF();
  return;
}

void ThirdDerivative::thirdDerivativeOut(int SRFlag, C3DMatrix *SR3Matrix)
{
  long unsigned int bufsize = 0;
  long unsigned int ioff = 0;
  psio_address next;

  bufsize = 2 * (long unsigned int) pow(gParams.ncartesians,3);
  ioff = abs(SRFlag) * bufsize;
  if(SRFlag > 0)
    ioff += (long unsigned int) pow(gParams.ncartesians,3);

  open_PSIF();
  
  next = psio_get_address(PSIO_ZERO, ioff);

  psio_write(PSIF_OPTKING, "Third Derivative / SR Matrix", (char *) (SR3Matrix), (long unsigned int) pow(gParams.ncartesians,3) * sizeof(double), next, &next);
  
  if(SRFlag > 0)
    fprintf(outfile, "Y Matrix for Coordinate %i written to File 1 TOC\n", abs(SRFlag));
  else
    fprintf(outfile, "SR Matrix for Coordinate %i written to File 1 TOC\n", abs(SRFlag));
  close_PSIF();
}

void FourthDerivative::thirdDerivativeBlockOut(int SRFlag, int LFlag, C3DMatrix *YRMatrix)
{
  long unsigned int bufsize = 0;
  long unsigned int ioff = 0;
  long int LR = 0; //What is this?
  psio_address next;

  if(SRFlag < 0)
    LR = -(2 * (SRFlag+1) * gParams.ncartesians * LFlag);
  else if(SRFlag > 0)
    LR = (2 * SRFlag - 1) * gParams.ncartesians * LFlag;
  else {
    fprintf(outfile,"\nR = 0 encountered in subroutine YOUT2");
    exit(0);
  }
  
  bufsize = (2 * (long unsigned int) pow(gParams.ncartesians,3) * (LR - 1)) + 1;
  ioff = abs(SRFlag) * bufsize;
  
  open_PSIF();

  next = psio_get_address(PSIO_ZERO, ioff);

  psio_write(PSIF_OPTKING, "Fourth Derivative / YR Block", (char *) (YRMatrix), (long unsigned int) pow(gParams.ncartesians,3) * sizeof(double), next, &next);
  
  //Need print check statement here
  close_PSIF();
}

void ThirdDerivative::thirdDerivativeIn(int SRFlag, C3DMatrix *YWrite)
{
  long unsigned int bufsize = 0;
  long unsigned int ioff = 0;
  psio_address next;
  
  bufsize = 2 * (long unsigned int) pow(gParams.ncartesians,3);
  ioff = abs(SRFlag) * bufsize;
  if(SRFlag > 0)
    ioff += (long unsigned int) pow(gParams.ncartesians,3);

  open_PSIF();
  
  next = psio_get_address(PSIO_ZERO, ioff);

  psio_read(PSIF_OPTKING, "Third Derivative / SR Matrix", (char *) (YWrite), (long unsigned int) pow(gParams.ncartesians,3) * sizeof(double), next, &next);
  
  if(SRFlag > 0)
    fprintf(outfile, "Y Matrix for Coordinate %i read from File 1 TOC\n", abs(SRFlag));
  else
    fprintf(outfile, "SR Matrix for Coordinate %i read from File 1 TOC\n", abs(SRFlag));
  close_PSIF();
  return;
}

void FourthDerivative::thirdDerivativeBlockIn(int SRFlag, int LFlag, C3DMatrix *YWrite)
{
  long unsigned int bufsize = 0;
  long unsigned int ioff = 0;
  psio_address next;
  
  bufsize = 2 * (long unsigned int) pow(gParams.ncartesians,3);
  ioff = abs(SRFlag) * bufsize;
  if(SRFlag > 0)
    ioff += (long unsigned int) pow(gParams.ncartesians,3);

  open_PSIF();
  
  next = psio_get_address(PSIO_ZERO, ioff);

  psio_read(PSIF_OPTKING, "Third Derivative / SR Matrix", (char *) (YWrite), (long unsigned int) pow(gParams.ncartesians,3) * sizeof(double), next, &next);
  
  if(SRFlag > 0)
    fprintf(outfile, "Y Matrix for Coordinate %i read from File 1 TOC\n", abs(SRFlag));
  else
    fprintf(outfile, "SR Matrix for Coordinate %i read from File 1 TOC\n", abs(SRFlag));
  close_PSIF();
  return;
}

void Transform::open_PSIF() {
  psio_open(PSIF_OPTKING, PSIO_OPEN_OLD);
  return;
}

void Transform::close_PSIF() {
  psio_close(PSIF_OPTKING, 1);
  return;
}

*/

}} // namespace psi::intder
