/*! \file
    \ingroup CPHF
    \brief Evaluate total dipole moment derivatives
*/

/*! \defgroup CPHF cphf: Solve the Coupled-Perturbed Hartree-Fock Equations */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <physconst.h>
#define EXTERN
#include "globals.h"

/* build_dipder(): Compute the total dipole moment derivatives from
** the derivative dipole integrals computed by "cints --oeprop" and
** the nuclear-perturbation CPHF coefficients computed in cphf_X().
** The equation for the electronic contribution to the dipole
** derivatives is (MO basis):
**
** (d mu/d X)_elec = 4 UX_ai * d_ai - 2 SX_ij d_ij + 2 (d^X)_ii
**
** where X is a nuclear coordinate, d_pq is a dipole moment integral,
** (d^X)_pq is a derivative dipole moment integral, UX_pq is a CPHF
** coefficient, and SX_pq is an overlap derivative.  Summation over
** all indices is implied. The nuclear contribution to the dipole
** derivatives is simple:
**
** (d mu_g/d X)_nuc = delta(g,X) Z_X
**
** where mu_g is the g-th component of the dipole moment (x, y, or z),
** and Z_X is the nuclear Z-value for the atom associated with
** coordinate X.
**
** For more details (and clearer notation) see:
**
** Y. Yamaguchi et al., "A New Dimension to Quantum Chemistry:
** Analytic Derivative Methods in Ab Initio Molecular Electronic
** Structure Theory", Oxford Press, New York, 1994.  Ch.17,
** esp. pp. 326-328.
**
** TDC, October 2002
*/

namespace psi { namespace cphf {

void build_dipder(double ***UX)
{
  int coord, ij, stat;
  int i, isym, ifirst, ilast;
  int j, jsym, jfirst, jlast;
  int a, asym, afirst, alast;
  double ***MU, **S;
  double *scratch, *scratch2;
  char *label;
  double **TMP1, **TMP2;
  double **MUX, **MUY, **MUZ; /* derivative dipole integrals */
  FILE *file17;

  MU = (double ***) malloc(3 * sizeof(double **));
  for(i=0; i < 3; i++) MU[i] = block_matrix(nmo,nmo);

  /* Grab the MO-basis dipole integrals from disk */
  scratch = init_array(ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_MO_MX, scratch, ntri, 0, 0, outfile);
  for(i=0,ij=0; i < nmo; i++)
    for(j=0; j <= i; j++,ij++)
      MU[0][i][j] = MU[0][j][i] = scratch[ij];
  zero_arr(scratch,ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_MO_MY, scratch, ntri, 0, 0, outfile);
  for(i=0,ij=0; i < nmo; i++)
    for(j=0; j <= i; j++,ij++)
      MU[1][i][j] = MU[1][j][i] = scratch[ij];
  zero_arr(scratch,ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_MO_MZ, scratch, ntri, 0, 0, outfile);
  for(i=0,ij=0; i < nmo; i++)
    for(j=0; j <= i; j++,ij++)
      MU[2][i][j] = MU[2][j][i] = scratch[ij];
  zero_arr(scratch,ntri);

  S = block_matrix(nmo,nmo);
  label = (char *) malloc(PSIO_KEYLEN * sizeof(char));
  scratch2 = init_array(noei_ao);
  TMP1 = block_matrix(nao, nao);
  TMP2 = block_matrix(nao, nao);

  MUX = block_matrix(nmo,nmo);
  MUY = block_matrix(nmo,nmo);
  MUZ = block_matrix(nmo,nmo);

  /* Build the total dipole moment derivatives */
  for(coord=0; coord < natom*3; coord++) {
    /* Grab the MO-basis overlap derivatives from disk */
    sprintf(label, "MO-basis Overlap Derivs (%d)", coord);
    stat = iwl_rdone(PSIF_OEI, label, scratch, ntri, 0, 0, NULL);
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';

    for(i=0,ij=0; i < nmo; i++)
      for(j=0; j <= i; j++,ij++)
	S[i][j] = S[j][i] = scratch[ij];

    for(isym=0; isym < nirreps; isym++) {
      ifirst = ofirst[isym];
      ilast = olast[isym];
      for(i=ifirst; i <= ilast; i++) {
	for(jsym=0; jsym < nirreps; jsym++) {
	  jfirst = ofirst[jsym];
	  jlast = olast[jsym];
	  for(j=jfirst; j <= jlast; j++) {
	    dipder[0][coord] -= 2.0 * S[i][j] * MU[0][i][j];
	    dipder[1][coord] -= 2.0 * S[i][j] * MU[1][i][j];
	    dipder[2][coord] -= 2.0 * S[i][j] * MU[2][i][j];
	  }
	}
      }
    }

    for(asym=0; asym < nirreps; asym++) {
      afirst = vfirst[asym];
      alast = vlast[asym];
      for(a=afirst; a <= alast; a++) {
	for(isym=0; isym < nirreps; isym++) {
	  ifirst = ofirst[isym];
	  ilast = olast[isym];
	  for(i=ifirst; i <= ilast; i++) {
	    dipder[0][coord] += 4.0 * UX[coord][a][i] * MU[0][a][i];
	    dipder[1][coord] += 4.0 * UX[coord][a][i] * MU[1][a][i];
	    dipder[2][coord] += 4.0 * UX[coord][a][i] * MU[2][a][i];
	  }
	}
      }
    }

    /* Read AO-basis dipole derivative integrals from disk and transform to MO basis*/
    sprintf(label, "AO-basis MUX Derivs (%d)", coord);
    stat = iwl_rdone(PSIF_OEI, label, scratch2, noei_ao, 0, 0, NULL);
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';
    for(i=0,ij=0; i < nao; i++)
      for(j=0; j <= i; j++,ij++)
	TMP1[i][j] = TMP1[j][i] = scratch2[ij];

    /*
      fprintf(outfile, "AO-basis MUX (%d)", coord);
      print_mat(TMP1, nao, nao, outfile);
    */

    C_DGEMM('n','t',nao,nso,nao,1,&(TMP1[0][0]),nao,&(usotao[0][0]),nao,
	    0,&(TMP2[0][0]),nao);
    C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(TMP2[0][0]),nao,
	    0,&(TMP1[0][0]),nao);

    C_DGEMM('n','n',nso,nmo,nso,1,&(TMP1[0][0]),nao,&(scf[0][0]),nmo,
	    0,&(TMP2[0][0]),nao);
    C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(TMP2[0][0]),nao,
	    0,&(MUX[0][0]),nmo);

    /*
      fprintf(outfile, "MO-basis MUX (%d)", coord);
      print_mat(MUX, nmo, nmo, outfile);
    */

    sprintf(label, "AO-basis MUY Derivs (%d)", coord);
    stat = iwl_rdone(PSIF_OEI, label, scratch2, noei_ao, 0, 0, NULL);
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';
    for(i=0,ij=0; i < nao; i++)
      for(j=0; j <= i; j++,ij++)
	TMP1[i][j] = TMP1[j][i] = scratch2[ij];

    /*
      fprintf(outfile, "AO-basis MUY (%d)", coord);
      print_mat(TMP1, nao, nao, outfile);
    */

    C_DGEMM('n','t',nao,nso,nao,1,&(TMP1[0][0]),nao,&(usotao[0][0]),nao,
	    0,&(TMP2[0][0]),nao);
    C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(TMP2[0][0]),nao,
	    0,&(TMP1[0][0]),nao);

    C_DGEMM('n','n',nso,nmo,nso,1,&(TMP1[0][0]),nao,&(scf[0][0]),nmo,
	    0,&(TMP2[0][0]),nao);
    C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(TMP2[0][0]),nao,
	    0,&(MUY[0][0]),nmo);

    /*
      fprintf(outfile, "MO-basis MUY (%d)", coord);
      print_mat(MUY, nmo, nmo, outfile);
    */

    sprintf(label, "AO-basis MUZ Derivs (%d)", coord);
    stat = iwl_rdone(PSIF_OEI, label, scratch2, noei_ao, 0, 0, NULL);
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';
    for(i=0,ij=0; i < nao; i++)
      for(j=0; j <= i; j++,ij++)
	TMP1[i][j] = TMP1[j][i] = scratch2[ij];

    /*
      fprintf(outfile, "AO-basis MUZ (%d)", coord);
      print_mat(TMP1, nao, nao, outfile);
    */

    C_DGEMM('n','t',nao,nso,nao,1,&(TMP1[0][0]),nao,&(usotao[0][0]),nao,
	    0,&(TMP2[0][0]),nao);
    C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(TMP2[0][0]),nao,
	    0,&(TMP1[0][0]),nao);

    C_DGEMM('n','n',nso,nmo,nso,1,&(TMP1[0][0]),nao,&(scf[0][0]),nmo,
	    0,&(TMP2[0][0]),nao);
    C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(TMP2[0][0]),nao,
	    0,&(MUZ[0][0]),nmo);

    /*
      fprintf(outfile, "MO-basis MUZ (%d)", coord);
      print_mat(MUZ, nmo, nmo, outfile);
    */

    for(isym=0; isym < nirreps; isym++) {
      ifirst = ofirst[isym];
      ilast = olast[isym];
      for(i=ifirst; i <= ilast; i++) {
	dipder[0][coord] += 2.0 * MUX[i][i];
	dipder[1][coord] += 2.0 * MUY[i][i];
	dipder[2][coord] += 2.0 * MUZ[i][i];
      }
    }
  }

  fprintf(outfile,"\n\tAtomic Polar Tensor:\n");
  fprintf(outfile,"\n\tUnits: au\n");
  fprintf(outfile,"\n\tTerms: Electronic\n");
  fprintf(outfile,"\n\t     Ex\t\t     Ey\t\t     Ez\n");
  for(i=0; i<natom*3; i++) {
    if((i%3)==0)
      fprintf(outfile,"%3sx\t%10.6lf\t%10.6lf\t%10.6lf\n",asymbol[i],dipder[0][i],dipder[1][i],dipder[2][i]);
    if((i%3)==1)
      fprintf(outfile,"%3sy\t%10.6lf\t%10.6lf\t%10.6lf\n",asymbol[i],dipder[0][i],dipder[1][i],dipder[2][i]);
    if((i%3)==2)
      fprintf(outfile,"%3sz\t%10.6lf\t%10.6lf\t%10.6lf\n",asymbol[i],dipder[0][i],dipder[1][i],dipder[2][i]);
    if((i+1)%3==0) fprintf(outfile,"\n");
  }


  /* Add nuclear contributions */
  for(coord=0; coord < natom; coord++) {
    dipder[0][coord*3] += zvals[coord];
    dipder[1][coord*3+1] += zvals[coord];
    dipder[2][coord*3+2] += zvals[coord];
  }

  fprintf(outfile,"\n\tAtomic Polar Tensor:\n");
  fprintf(outfile,"\n\tUnits: au\n");
  fprintf(outfile,"\n\tTerms: Electronic + Nuclear\n");
  fprintf(outfile,"\n\t     Ex\t\t     Ey\t\t     Ez\n");
  for(i=0; i<natom*3; i++) {
    if((i%3)==0)
      fprintf(outfile,"%3sx\t%10.6lf\t%10.6lf\t%10.6lf\n",asymbol[i],dipder[0][i],dipder[1][i],dipder[2][i]);
    if((i%3)==1)
      fprintf(outfile,"%3sy\t%10.6lf\t%10.6lf\t%10.6lf\n",asymbol[i],dipder[0][i],dipder[1][i],dipder[2][i]);
    if((i%3)==2)
      fprintf(outfile,"%3sz\t%10.6lf\t%10.6lf\t%10.6lf\n",asymbol[i],dipder[0][i],dipder[1][i],dipder[2][i]);
    if((i+1)%3==0) fprintf(outfile,"\n");
  }

  /* Convert to debye/bohr */
  for(i=0; i < 3; i++)
    for(j=0; j < natom*3; j++)
      dipder[i][j] *= _dipmom_au2debye;
  
  if (print_lvl > 4) {
    fprintf(outfile, "\n\tDipole Derivatives W.R.T. Cartesian Coordinates (debye/a0):\n");
    print_mat(dipder, 3, natom*3, outfile);
  }

  /* Convert to debye/A */
  for(i=0; i < 3; i++)
    for(j=0; j < natom*3; j++)
      dipder[i][j] /= _bohr2angstroms;

  /*
  fprintf(outfile, "\n\tDipole Derivatives W.R.T. Cartesian Coordinates (debye/A):\n");
  print_mat(dipder, 3, natom*3, outfile);
  */

  /* write the dipole derivatives to file17 in the PSI2 standard format */
  ffile(&file17, "file17.dat", 0);
  fprintf(file17, "%5d%5d\n", natom, natom*3);
  for(i=0; i < 3; i++) {
    for(j=0; j < natom; j++) {
      fprintf(file17, "%20.10f%20.10f%20.10f\n", dipder[i][j*3], dipder[i][j*3+1], dipder[i][j*3+2]);
    }
  }
  fclose(file17);

  free_block(MUX);
  free_block(MUY);
  free_block(MUZ);
  free_block(TMP1);
  free_block(TMP2);
  free(scratch2);
  free_block(S);
  free(label);
  for(i=0; i < 3; i++) free_block(MU[i]);
  free(MU);
  free(scratch);
}

}} // namespace psi::cphf
