/*! \file
    \ingroup CUSP
    \brief Enter brief description of file here 
*/
/*
** phi(): Compute the values of the atomic orbitals at a given point.
** -TDC, June 2001
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libint/libint.h>  /* for the maximum angluar momentum, LIBINT_MAX_AM */
#include <psiconfig.h>
#include <psifiles.h>

namespace psi { namespace cusp {

#define MAXFACT 100

extern int nao;
int natom, nprim, nshell;
int *stype, *snuc, *snumg, *sprim;
double *a, *c, **geom, df[MAXFACT], **norm;
int **xexp, **yexp, **zexp, *l_length;

void setup_phi(void)
{
  static int use_cca_integrals_standard = 1;
  static int done=0;
  int i,l,j,ao;

  if(done) return;

  chkpt_init(PSIO_OPEN_OLD);
  nao = chkpt_rd_nao();
  natom = chkpt_rd_natom();
  nprim = chkpt_rd_nprim();
  nshell = chkpt_rd_nshell();
  stype = chkpt_rd_stype();
  snuc = chkpt_rd_snuc();
  snumg = chkpt_rd_snumg();
  sprim = chkpt_rd_sprim();
  a = chkpt_rd_exps();
  c = chkpt_rd_contr();
  geom = chkpt_rd_geom();
  chkpt_close();

  /* compute double factorial --- df[n] = (n-1)!! */
  df[0] = 1.0;  df[1] = 1.0;  df[2] = 1.0;
  for(i=3; i < MAXFACT; i++) df[i] = (i-1) * df[i-2];
  
  /* Compute the AO ordering within a shell and the angle-independent
     normalization contants */
  xexp = (int **) malloc(LIBINT_MAX_AM * sizeof(int *));
  yexp = (int **) malloc(LIBINT_MAX_AM * sizeof(int *));
  zexp = (int **) malloc(LIBINT_MAX_AM * sizeof(int *));
  norm = (double **) malloc(LIBINT_MAX_AM * sizeof(double *));
  l_length = init_int_array(LIBINT_MAX_AM);
  l_length[0] = 1;
  for(l=0; l < (LIBINT_MAX_AM); l++) {

      if(l) l_length[l] = l_length[l-1] + l + 1;
      
      xexp[l] = init_int_array(l_length[l]);
      yexp[l] = init_int_array(l_length[l]);
      zexp[l] = init_int_array(l_length[l]);
      norm[l] = init_array(l_length[l]);

      for(i=0,ao=0; i <= l; i++) {

	  for(j=0; j <= i; j++) {

	      xexp[l][ao] = l - i;
	      yexp[l][ao] = i - j;
	      zexp[l][ao] = j;

	      if (use_cca_integrals_standard)
		norm[l][ao] = 1.0;
	      else
		norm[l][ao] = sqrt(df[2*l]/(df[2*(l-i)]*df[2*(i-j)]*df[2*j]));

/*	      printf("%d %d %20.10f\n", l, ao, norm[l][ao]); */

	      ao++;
	    }
	}
    }

  done = 1;

}

void compute_phi(double *phi, double x, double y, double z)
{
  int i, j, l, firstp, lastp;
  int am, atom, ao;
  double cexpr, dx, dy, dz, rr;

  setup_phi();

  /* Loop over the shells */
  for(i=0,ao=0; i < nshell; i++) {

      firstp = sprim[i]-1;
      lastp = firstp + snumg[i];

      atom = snuc[i] - 1;
      dx = x - geom[atom][0];
      dy = y - geom[atom][1];
      dz = z - geom[atom][2];
      
      /* Compute r^2 for this center */
      rr = dx*dx + dy*dy + dz*dz;

      /* Angular momentum level for this shell */
      am = stype[i] - 1;

      if(am > LIBINT_MAX_AM) {
	  printf("Angular momentum max. exceeded.\n");
	  exit(PSI_RETURN_FAILURE);
	}

      /* Loop over the primitive Gaussians in this shell */
      cexpr = 0;
      for(j=firstp; j < lastp; j++) {
          cexpr += c[j] * exp(-a[j] * rr);
        }

      for(l=0; l < l_length[am]; l++) {
          phi[ao+l] += pow(dx,(double) xexp[am][l]) *
		       pow(dy,(double) yexp[am][l]) *
		       pow(dz,(double) zexp[am][l]) *
		       cexpr * norm[am][l];
	}

      ao += l_length[am];

    }
}


}} // namespace psi::cusp
