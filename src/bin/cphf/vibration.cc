/*! \file
    \ingroup CPHF
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include <physconst.h>
#include <masses.h>
#define EXTERN
#include "globals.h"

#define _D2esucm 1e-18

namespace psi { namespace cphf {

/* vibration(): Computes the harmonic vibrational frequencies and IR
** integrated absorption coefficients (intensities) using the
** cartesian hessian from build_hessian() and dipole derivatives from
** build_dipder().  
**
** The procedure for computing the normal modes is described in:
**
**  E.B. Wilson, J.C Decius, and P.C. Cross, "Molecular Vibrations:
** The Theory of Infrared and Raman Vibrational Spectra", Dover, New
** York, 1955.
**
** The relationship between IR intensities and dipole derivatives is
** described in the same text (Ch. 7) as well as in:
**
** I. Mills, T. Cvitas, K. Homann, N. Kallay, and K. Kuchitsu,
** "Quantities, Units, and Symbols in Physical Chemistry", Oxford,
** London, 1993, pp. 33-35.
**
** B.S. Galabov and T. Dudev, "Vibrational Intensities", v.22 of
** Vibrational Spectra and Structure, J.R. Durig, ed., Elsevier, 1996,
** pp. 2-12.
**
** G. Zerbi, "Introduction to the Theory of Vibrational Frequencies
** and Vibrational Intensities," in Vibrational Intensities in
** Infrared and Raman Spectroscopy, v.20 of Studies in Physical and
** Theoretical Chemistry, W.B. Person and G. Zerbi, eds., Elsevier
** 1982, pp. 45-51.
**
** TDC, October 2002.
*/

void vibration(double **hessian, double **lx)
{
  int i, j, k;
  double **M, *irint;
  double **TMP;
  double *km, k_convert, cm_convert;
  double *work;
  int stat;
  double dipder_conv, ir_prefactor;
  double ds;
  double freq;
  double *vib_temp, *vib_energy, T, *q_vib, *s_vib;
  double ZPVE, evib;
  double etrans, erot, ethermal;
  double molecmass, qnuc;
  double P, qtrans;
  double qelec, qvib, qtotal;
  double Selec, Snuc, Strans, Srot, Svib, Stotal;

  /* Print out the force constants for subsequent optimizations */
  TMP = block_matrix(natom*3,natom*3);
  for(i=0; i<3*natom; i++) {
    for(j=0; j<3*natom; j++)
      TMP[i][j] = hessian[i][j] * _hartree2J * 1.0E18 /
        (_bohr2angstroms * _bohr2angstroms);
  }
  psio_open(PSIF_OPTKING,PSIO_OPEN_OLD);
  fprintf(outfile,"\n\tWriting Cartesian Force Constants to PSIF_OPTKING file\n");
  psio_write_entry(PSIF_OPTKING, "Cartesian Force Constants",
        (char *) &(TMP[0][0]),3*natom*3*natom*sizeof(double));
  psio_close(PSIF_OPTKING,1);
  // fprintf(outfile, "\n\tHessian matrix in aJ/Ang^2:\n");
  // print_mat(TMP, natom*3, natom*3, outfile);
  free_block(TMP);

  /* mass-weight the hessian */
  M = block_matrix(natom*3, natom*3);
  for(i=0; i < natom; i++) {
    for(j=0; j < 3; j++)  {
      M[i*3+j][i*3+j] = 1/sqrt(an2masses[(int) zvals[i]]);
    }
  }

  //fprintf(outfile, "\n\tM^-1/2 matrix:\n");
  //print_mat(M, natom*3, natom*3, outfile);

  TMP = block_matrix(natom*3,natom*3);
  C_DGEMM('n','n', natom*3, natom*3, natom*3, 1.0, &(M[0][0]), natom*3,
          &(hessian[0][0]), natom*3, 0.0, &(TMP[0][0]), natom*3);
  C_DGEMM('n','n', natom*3, natom*3, natom*3, 1.0, &(TMP[0][0]), natom*3,
	  &(M[0][0]), natom*3, 0.0, &(hessian[0][0]), natom*3);
  free_block(TMP);

  //fprintf(outfile, "\n\tMass-Weighted Hessian matrix:\n");
  //print_mat(hessian, natom*3, natom*3, outfile);

  /* diagonalize mass-weighted hessian */
  km = init_array(natom*3);  /* mass-weighted force constants */
  work = init_array(natom*3*3); /* scratch array */

  if(stat = C_DSYEV('v','u',natom*3,&(hessian[0][0]),natom*3,&(km[0]),&(work[0]),natom*3*3)) {
    fprintf(outfile, "vibration(): Error in hessian diagonalization. stat = %d\n", stat);
    throw PsiException("Error", __FILE__, __LINE__);
  }

  /*
  fprintf(outfile, "\n\tEigenvalues of Diagonalized Hessian Matrix\n");
  for(i=0; i<natom*3; i++) fprintf(outfile,"\t%d\t%12.8lf\n",i,km[i]);
  fprintf(outfile,"\n");
  */

  /* Construct the mass-weighted normal coordinates, lx */
  /* (note, after C_DSYEV, hessian contains the eigenvectors) */
  C_DGEMM('n','t',natom*3,natom*3,natom*3,1,&(M[0][0]),natom*3,
          &(hessian[0][0]),natom*3,0,&(lx[0][0]),natom*3);

  if(print_lvl & 1) {
    fprintf(outfile, "\n\tLX matrix:\n");
    print_mat(lx, natom*3, natom*3, outfile);
  }

  /* Transform dipole derivatives to normal coordinates */
  C_DGEMM('n','n',3,natom*3,natom*3,1,&(dipder[0][0]),natom*3,
          &(lx[0][0]),natom*3,0,&(dipder_q[0][0]),natom*3);

  /*
  fprintf(outfile, "\n\tDipole Derivatives W.R.T Normal Coordinates (debye/A-amu+1/2):\n");
  print_mat(dipder_q, 3, natom*3, outfile);
  */

  /* Compute the IR intensities in D^2/(A^2 amu) */
  irint = init_array(natom*3);
  for(i=0; i < natom*3; i++)
    for(j=0; j < 3; j++)
      irint[i] += dipder_q[j][i] * dipder_q[j][i];

  /*
  fprintf(outfile,"\nIR intensities in D^2/(A^2 amu)\n");
  for(i=natom*3-1; i >= (natom*3-nnc); i--) {
    fprintf(outfile,"%12.6lf\n",irint[i]);
  }
  */

  k_convert = _hartree2J/(_bohr2m * _bohr2m * _amu2kg);
  cm_convert = 1.0/(2.0 * _pi * _c * 100.0);

  /*
  fprintf(outfile,"\nDipole Strength (10^-40 esu^2 cm^2):\n");
  for(i=natom*3-1; i >= (natom*3-nnc); i--) {
    freq = cm_convert*sqrt(k_convert*km[i]);
    ds = (_h*_D2esucm*_D2esucm*irint[i]*1e40)/(_amu2kg*2*_c*4*_pi*_pi*100*freq*1e-20);
    fprintf(outfile,"%10.4lf\n",ds);
  }
  */
  
  /* conversion factor from D^2/(A^2 amu) to C^2/kg */
  dipder_conv = _dipmom_debye2si*_dipmom_debye2si/(1e-20 * _amu2kg);
  for(i=0; i < natom*3; i++)
    irint[i] *= dipder_conv;

  /* IR integrated absorption coefficient prefactor */
  ir_prefactor = _na * _pi/(3.0 * _c * _c * 4.0 * _pi * _e0 * 1000.0);

  /* fprintf(outfile, "\n\tIR conversion = %20.10f\n", dipder_conv * ir_prefactor); */

  /* compute the frequencies and spit them out in a nice table */
  fprintf(outfile, "\n\t        Harmonic Frequency   Infrared Intensity\n");
  fprintf(outfile,   "\t              (cm-1)               (km/mol)    \n");
  fprintf(outfile,   "\t-----------------------------------------------\n");
  k_convert = _hartree2J/(_bohr2m * _bohr2m * _amu2kg);
  cm_convert = 1.0/(2.0 * _pi * _c * 100.0);
  for(i=natom*3-1; i >= 0; i--) {
    if(km[i] < 0.0)
      fprintf(outfile, "\t  %3d   %17.3fi       %10.4f\n", (natom*3-i), 
              cm_convert * sqrt(-k_convert * km[i]), irint[i]*ir_prefactor);
    else
      fprintf(outfile, "\t  %3d   %17.3f        %10.4f\n", (natom*3-i), 
              cm_convert * sqrt(k_convert * km[i]), irint[i]*ir_prefactor);
  }
  fprintf(outfile,   "\t-----------------------------------------------\n");

  fprintf(outfile, "\nNormal Modes (mass-weighted)\n");

  for(i=0; i < 3*natom; i++) {
    if (fabs(cm_convert * sqrt(k_convert * fabs(km[i]))) < 5.0) continue;
    fprintf(outfile,"\n");
    if(km[i] < 0.0)
      fprintf(outfile, "   Frequency:    %8.2fi\n", cm_convert * sqrt(-k_convert * km[i]));
    else
      fprintf(outfile, "   Frequency:    %8.2f\n", cm_convert * sqrt(k_convert * km[i]));
    fprintf(outfile,   "   IR Intensity: %8.2f\n", irint[i]*ir_prefactor);
    fprintf(outfile, "\t     X       Y       Z \t\n");
    for(j=0; j < natom; j++) {
      fprintf(outfile, "  %s \t", asymbol[3*j]);
      for(k=0; k < 3; k++)  { 
        fprintf(outfile, "%8.3f", lx[3*j+k][i]);
      }
      fprintf(outfile, "\n");
    }  
  }

  double *vibs = init_array(3*natom);
  for(i=0; i<3*natom; ++i) {
    if (km[3*natom-1-i] < 0.0)
      vibs[i] = -1.0 * cm_convert * sqrt(-k_convert * km[3*natom-1-i]);
    else
      vibs[i] = cm_convert * sqrt(k_convert * km[3*natom-1-i]);
  }
 
  // chkpt_wt_vib_freqs will only write the first 3n-5/6 frequencies
  chkpt_init(PSIO_OPEN_OLD);
 // chkpt_wt_vib_freqs(vibs);
  chkpt_close();
  free(vibs);


/*
  // compute G, H, and S 
  fprintf(outfile, "\n Thermochemistry \n");
  // compute Vib Temp (K)
  vib_temp = init_array(natom*3);
  T = 298.15;
  P = 101325;
  vib_energy = init_array(natom*3);
  q_vib = init_array(natom*3);
  s_vib = init_array(natom*3);
  ZPVE = 0.0;
  evib = 0.0;
  for(i=0; i < 3*natom; i++)  {
    vib_temp[i] = _c * 100 * _h * cm_convert * sqrt(k_convert * km[i]) / _kb;
    //fprintf(outfile, "%16.8f\n", vib_temp[i]);
    vib_energy[i] = vib_temp[i] / T * (0.5 + (1 / ((exp(vib_temp[i] / T) - 1))));
    //fprintf(outfile, "%17.6f\n", vib_energy[i]);
    q_vib[i] = (exp(-vib_temp[i] / (2*T))) / (1 - exp(-vib_temp[i] / T));
    //fprintf(outfile, "%7.6f\n", q_vib[i]);
    s_vib[i] = vib_temp[i] / (T * (exp(vib_temp[i] / T) - 1)) - log(1 - exp(-vib_temp[i] / T));
    //fprintf(outfile, "%7.6f\n", s_vib[i]);
    ZPVE += (cm_convert * sqrt(k_convert * km[i]) / 2);
    //fprintf(outfile, "%17.6f\n", ZPVE);
    evib += vib_energy[i];
    //fprintf(outfile, "%17.6f\n", evib);
  }

  etrans = 1.5;
  if(rottype == 3) 
    erot = 1;
  else
    erot = 1.5; 
  ethermal = evib + etrans + erot;
  //fprintf(outfile, "%17.6f\n", ethermal);

  // calculating the partition functions
  qelec = 1 ;//degen; 
  for(i=0; i < natom; i++)
    molecmass += an2masses[(int) zvals[i]]; 
  //fprintf(outfile, "%17.6f\n", molecmass);
  fprintf(outfile, "%17.6f\n", qelec);
  qnuc = 1;
  qtrans = pow((_twopi * molecmass * _amu2kg * _kb * T / (_h * _h)), 1.5) * (_na * _kb * T / (P * _na));
  //fprintf(outfile, "%17.6f\n", qtrans);
//  qrot = pow(_pi, 0.5) / 1 * pow((_kb * T / _h), 1.5) * pow((A * B * C), -0.5);
//  fprintf(outfile, "%17.6f\n", qrot);
  qvib = 1;
  for(i=0; i < natom*3; i++)
    qvib = qvib * q_vib[i];
  //fprintf(outfile, "%17.6f\n", qvib);
  qtotal = qelec * qnuc * qtrans * qvib; // qrot
  //fprintf(outfile, "%17.6f\n", qtotal);
  
  // calculating the contributions to entropy
  Selec = _psi3_R * log(qelec) / _cal2J;
  //fprintf(outfile, "%17.6f\n", Selec);
  Snuc = 0;
  Strans = _psi3_R * (log(qtrans) + 2.5) / _cal2J;
  //Srot = _psi3_R * (log(qrot) + 1.5) / _cal2J;
  for(i=0; i < natom*3; i++)
    Svib += s_vib[i];
  Svib = Svib * _psi3_R / _cal2J;
  //fprintf(outfile, "%17.6f\n", Svib);
  Stotal = Selec + Snuc + Strans + Svib; // + Srot
  //fprintf(outfile, "%17.6f\n", Stotal);
*/

  free(work);
  free(irint);
  free_block(M);
  free(km);
}

}} // namespace psi::cphf
