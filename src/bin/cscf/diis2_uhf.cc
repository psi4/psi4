/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* diis_uhf(): DIIS extrapolation for the UHF SCF procedure.  This version
** uses Pulay's method [see P. Pulay, J. Comp. Chem. 3, 556-560 (1982)].
**
** (1) The elements of the B-matrix are taken to be the sums of the dot products
** involving alpha and beta error vectors: B_ij = <e_i|e_j>_alpha + <e_i|e_j>_beta.
** I tried independent extrapolation of alpha and beta, but this failed horribly.
**
** (2) The overlap matrix used for the expression e = FDS - SDF is the
** P^-1 built in shalf.c rather than the true overlap.  This is
** important for linearly dependent basis sets, I think.
**
** The original version of this code was written by S. Brown,
** ca. 11/99.  This version was written by TDC, 4/04.
*/

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

extern double delta;
static double *btemp, **bold, **bmat;

static struct diis_data {
    double ****fock;
    double ****error;
} *diism, dtemp;

void diis_uhf(void)
{
  int i,j,k,ij, kk;
  int errcod;
  int a,b,c,e;
  int m,n,nn,mm,num_mo;
  int last, col, _try;
  double etemp, dotp, norm, determ;
  double scale;
  struct symm *s;
  struct diis_data *d;
  int diis_print;
  double **S, **F, **D, **X, **Y, **Z;
  double sum, maximum;
    
  if(diism == NULL){ /* first call, so allocate memory */

    /* use one DIIS b-matrix for alpha+beta error vectors */
    bmat = block_matrix(ndiis+1,ndiis+1);
    bold = block_matrix(ndiis,ndiis);
    btemp = init_array(ndiis+1);

    diism = (struct diis_data *) malloc(sizeof(struct diis_data)*ndiis);
    for(m=0; m < ndiis ; m++) {
      d = &diism[m];
      d->fock = (double ****) malloc(2 * sizeof(double ***));
      d->error = (double ****) malloc(2 * sizeof(double ***));
      for(n=0; n < 2; n++) {
	d->fock[n] = (double ***) malloc(num_ir * sizeof(double **));
	d->error[n] = (double ***) malloc(num_ir * sizeof(double **));
	for(j=0;j < num_ir; j++){
	  if(nn=scf_info[j].num_so) {
	    d->fock[n][j] = block_matrix(nn,nn);
	    d->error[n][j] = block_matrix(nn,nn);
	  }
	}
      }
    }
  }	

  scale = 1.0 + dampsv;

  S = block_matrix(nsfmax,nsfmax);  /* to store overlap matrix for DIIS */
  D = block_matrix(nsfmax,nsfmax);  /* to store density matrix for DIIS */
  X = block_matrix(nsfmax,nsfmax);  /* temp matrix for DIIS */
  Z = block_matrix(nsfmax,nsfmax);  /* temp matrix for DIIS */

  last = iter-1;
  col = iter+1;
  if(iter > ndiis) {
    last = ndiis-1;
    col = ndiis+1;
  }

  if (iter > ndiis) {
    /* shift DIIS structs up */
    dtemp = diism[0];
    for (i=0; i < last ; i++) {
      diism[i] = diism[i+1];
    }
    diism[last] = dtemp;
  } 
	
  diis_print = 0;
  errcod = ip_boolean("DIIS_PRINT",&diis_print,0);
  if(iter == 1 && diis_print)
    ffile(&diis_out,"diis_out.dat",0);

  d = &diism[last];
  for(n=0; n < 2; n++){
    for(m=0; m < num_ir; m++) {
      s = &scf_info[m];
      if(nn=s->num_so) {
	num_mo = s->num_mo;

	/*
	zero_mat(S, nsfmax, nsfmax);
	tri_to_sq(scf_info[m].smat, S, nn);
	*/

	/* Generate the density for this irrep */
	zero_mat(D, nsfmax, nsfmax);
	C_DGEMM('n', 't', nn, nn, spin_info[n].scf_spin[m].noccup, 1.0,
		spin_info[n].scf_spin[m].cmat[0], nn,
		spin_info[n].scf_spin[m].cmat[0], nn, 0.0, D[0], nsfmax);

	zero_mat(X, nsfmax, nsfmax);
	zero_mat(Z, nsfmax, nsfmax);

	tri_to_sq(spin_info[n].scf_spin[m].fock_pac, d->fock[n][m], nn);

	/* SDF */
	C_DGEMM('n', 'n', nn, nn, nn, 1.0, s->pinv[0], nn, D[0], nsfmax, 0.0, X[0], nsfmax);
	C_DGEMM('n', 'n', nn, nn, nn, 1.0, X[0], nsfmax, d->fock[n][m][0], nn, 0.0, Z[0], nsfmax);
	zero_mat(X, nsfmax, nsfmax);

	/* FDS - SDF */
	C_DGEMM('n', 'n', nn, nn, nn, 1.0, d->fock[n][m][0], nn, D[0], nsfmax, 0.0, X[0], nsfmax);
	C_DGEMM('n', 'n', nn, nn, nn, 1.0, X[0], nsfmax, s->pinv[0], nn, -1.0, Z[0], nsfmax);

	/* Transform the error vector to the orthogonal AO basis */
	C_DGEMM('t', 'n', num_mo, nn, nn, 1.0, s->sahalf[0], nn, Z[0], nsfmax, 0.0, X[0], nsfmax);
	C_DGEMM('n', 'n', num_mo, num_mo, nn, 1.0, X[0], nsfmax, s->sahalf[0], nn, 0.0, 
		d->error[n][m][0], nn);

	/*
	for(i=0; i < nn; i++)
	  for(j=0; j < nn; j++)
	    d->error[n][m][i][j] = Z[i][j];
	*/

	for(i=0; i < num_mo; i++)
	  for(j=0; j <= i; j++) {
	    etemp = fabs(d->error[n][m][i][j]);
	    diiser = MAX0(diiser,etemp);
	  }

      }
    }
  }

  free_block(S);
  free_block(D);
  free_block(X);
  free_block(Z);

  if(iter > ndiis) {
    for(i=0; i < last; i++) {
      for(j=0; j <= i; j++) {
	bold[i][j] = bold[j][i] = bold[i+1][j+1];
      }
    }
  }

  for(i=0; i <= last; i++) {
    etemp = 0.0;
    for(n=0; n < 2; n++) {
      for(m=0; m < num_ir; m++) {
	s = &scf_info[m];
	if(nn=s->num_so) {
	  sdot(diism[i].error[n][m], diism[last].error[n][m], s->num_mo, &dotp);
	  etemp += dotp;
	}
      }
    }
    bold[i][last] = bold[last][i] = etemp;
  }

  if(diis_print) {
    fprintf(diis_out, "\nRaw BMAT for iter %d.\n", iter);
    print_mat(bold, last+1, last+1, diis_out);
  }

  bmat[0][0] = 0.0;
  btemp[0] = -1.0;
  /*    norm = 1.0/bold[n][0][0]; */
  for(i=1; i <= last+1; i++) {
    bmat[i][0] = bmat[0][i] = -1.0;
    btemp[i] = 0.0;
    for(j=1; j <= i; j++) {
      bmat[i][j] = bmat[j][i] = bold[i-1][j-1];
      if(i==j) bmat[i][j] *= scale;
    }
  }

  /* find the maximum in B and scale all the elements */
  maximum = fabs(bmat[1][1]);
  for(i=1; i <= last+1; i++) {
    for(j=1; j <= i; j++) {
      if(fabs(bmat[i][j]) > maximum) maximum = fabs(bmat[i][j]);
    }
  }
  for(i=1; i <= last+1; i++) {
    for(j=1; j <= last+1; j++) {
      bmat[i][j] /= maximum;
    }
  }

  if(diis_print){
    fprintf(diis_out,"\nBMAT for iter %d\n", iter);
    print_mat(bmat,col,col,diis_out);
  }

  if(iter-1) {
    flin(bmat,btemp,col,1,&determ);

    if(diis_print) {
      fprintf(diis_out, "BMAT determinant for iter %d = %20.12f\n", iter, determ);
      fprintf(diis_out, "DIIS coeffs for iter %d\n", iter);
      sum = 0;
      for(i=0; i < col; i++) {
	fprintf(diis_out, "%d %20.12f\n", i, btemp[i]);
	if(i) sum += btemp[i];
      }
      fprintf(diis_out, "sum of DIIS coeffs for iter %d = %20.12f\n", iter, sum);
    }

    _try = 0;
    while(fabs(determ) < 1.0e-19 && _try < last) {
      _try++;
      col--;
      bmat[0][0] = 0.0;
      btemp[0] = -1.0;
      norm = 1.0/bold[_try][_try];
      for(i=1; i <= ndiis-_try; i++) {
	bmat[i][0] = bmat[0][i] = -1.0;
	btemp[i] = 0.0;
	for(j=1; j <= i; j++) {
	  bmat[i][j] = bmat[j][i] = bold[i+_try-1][j+_try-1]*norm;
	  if(i==j) bmat[i][j] *= scale;
	}
      }

      if(diis_print) {
	fprintf(diis_out, "\nCorrected BMAT for iter %d\n", iter);
	print_mat(bmat, col, col, diis_out);
      }

      flin(bmat, btemp, col, 1, &determ);
    }

    if(fabs(determ) < 1e-20) {
      printf("\nDIIS extrapolation failed in iter %d\n", iter);
      return;
    }

    if((iter >= it_diis)) {
      for(n=0; n < 2; n++) {
	for(m=0; m < num_ir; m++) {
	  s = &scf_info[m];
	  if(nn=s->num_so) {
	    num_mo = s->num_mo;
	    for(i=0,ij=0; i < nn; i++) {
	      for(j=0; j <= i; j++,ij++) {
		etemp =0.0;
		for(k=_try,kk=1; k < last+1; k++,kk++)
		  etemp += btemp[kk] * diism[k].fock[n][m][i][j];

		spin_info[n].scf_spin[m].fock_pac[ij] = etemp;
	      }
	    }

	  }

	}
      }
    }

  } /* if(iter-1) */
}

}} // namespace psi::cscf
