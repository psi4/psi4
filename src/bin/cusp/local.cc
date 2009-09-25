/*! \file
    \ingroup CUSP
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libint/libint.h>
#include <libchkpt/chkpt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libqt/qt.h>

extern FILE *outfile;

namespace psi { namespace cusp {

void local(void)
{
  int iter, s, t, A, k, l, m, inew, iold, max;
  int i, j, ij, am, atom, shell_length, offset, stat;
  int nirreps, nao, nmo, nso, natom, nshell, noei, nocc;
  int *stype, *snuc, *aostart, *aostop, *ao2atom, *l_length;
  int *clsdpi, *openpi, *orbspi, *dummy, *order;
  double *ss, **S, **scf, **u, **Ctmp, **C, *evals_tmp, *evals;
  double P, PiiA, Pst, Pss, Ptt, Ast, Bst, AB;
  double Uss, Utt, Ust, Uts, Cks, Ckt, **U, **V, **VV, **F;
  double cos4a, alpha, alphamax, alphalast, conv;
  double **R, **D;
  double norm, *charge, tmp;
  double *SR, **X, *Z;
  int row, col, junk;
  int *rank, *boolean, **domain, next_atom;
  int *orb_order, *orb_boolean;
  double fR, cutoff, det;
  double **TMP, **MuX, **MuY, **MuZ, *scratch;

  chkpt_init(PSIO_OPEN_OLD);
  nao = chkpt_rd_nao();
  nmo = chkpt_rd_nmo();
  nso = nmo;
  natom = chkpt_rd_natom();
  nshell = chkpt_rd_nshell();
  stype = chkpt_rd_stype();
  snuc = chkpt_rd_snuc();
  u = chkpt_rd_usotao();
  nirreps = chkpt_rd_nirreps();
  clsdpi = chkpt_rd_clsdpi();
  openpi = chkpt_rd_openpi();
  orbspi = chkpt_rd_orbspi();
  scf = chkpt_rd_scf();
  evals_tmp = chkpt_rd_evals();
  chkpt_close();

  /* Compute the length of each AM block */
  junk = LIBINT_MAX_AM;
  l_length = init_int_array(LIBINT_MAX_AM);
  l_length[0] = 1;
  for(l=0; l < (LIBINT_MAX_AM); l++) if(l) l_length[l] = l_length[l-1] + l + 1;

  /* Set up the atom->AO and AO->atom lookup arrays */
  aostart = init_int_array(natom);
  aostop = init_int_array(natom);
  for(i=0,atom=-1,offset=0; i < nshell; i++) {
    am = stype[i] - 1;
    shell_length = l_length[am];

    if(atom != snuc[i]-1) {
      if(atom != -1) aostop[atom] = offset-1;
      atom = snuc[i]-1;
      aostart[atom] = offset;
    }

    offset += shell_length;
  }
  aostop[atom] = offset-1;

  ao2atom = init_int_array(nao);
  for(i=0; i < natom; i++)
    for(j=aostart[i]; j <= aostop[i]; j++) ao2atom[j] = i;

  /* Get the overlap integrals */
  noei = nmo*(nmo+1)/2;
  ss = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_S,ss,noei,0,0,outfile);
  S = block_matrix(nao,nao);
  for(i=0,ij=0; i < nmo; i++)
    for(j=0; j <= i; j++,ij++) {
      S[i][j] = S[j][i] = ss[ij];
    }
  free(ss);


  /* transform the MO coefficients to the AO basis */
  Ctmp = block_matrix(nao,nmo);
  C_DGEMM('t','n',nao,nmo,nso,1,&(u[0][0]),nao,&(scf[0][0]),nmo,
	  0,&(Ctmp[0][0]),nmo);

  /* Sort the MOs to an occupation ordering */
  dummy = init_int_array(nmo);
  order = init_int_array(nmo);
  C = block_matrix(nao,nmo);
  evals = init_array(nmo);
  reorder_qt(clsdpi, openpi, dummy, dummy, order, orbspi, nirreps);
  for(i=0; i < nmo; i++) {
    inew = order[i];
    for(j=0; j < nao; j++) C[j][inew] = Ctmp[j][i];
    evals[inew] = evals_tmp[i];
  }
  free(dummy);
  free(order);
  free_block(Ctmp);
  free(evals_tmp);

  /* Grab the dipole moment integrals and transform them to the MO basis */
  /*
  TMP = block_matrix(nao,nao);
  X = block_matrix(nao,nmo);

  scratch = init_array(noei);

  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MX,scratch,noei,0,0,outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = TMP[j][i] = scratch[ij];
    }
  zero_arr(scratch,noei);

  MuX = block_matrix(nmo,nmo);
  C_DGEMM('n','n',nao,nmo,nao,1,&(TMP[0][0]),nao,&(C[0][0]),nmo,
	  0,&(X[0][0]),nmo);
  C_DGEMM('t','n',nmo,nmo,nao,1,&(C[0][0]),nmo,&(X[0][0]),nmo,
	  0,&(MuX[0][0]),nmo);

  fprintf(outfile, "\tMO-basis MuX:\n");
  print_mat(MuX,nmo,nmo,outfile);
  free_block(MuX);

  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MY,scratch,noei,0,0,outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = TMP[j][i] = scratch[ij];
    }
  zero_arr(scratch,noei);

  MuY = block_matrix(nmo,nmo);
  C_DGEMM('n','n',nao,nmo,nao,1,&(TMP[0][0]),nao,&(C[0][0]),nmo,
          0,&(X[0][0]),nmo);
  C_DGEMM('t','n',nmo,nmo,nao,1,&(C[0][0]),nmo,&(X[0][0]),nmo,
          0,&(MuY[0][0]),nmo);
  
  fprintf(outfile, "\tMO-basis MuY:\n");
  print_mat(MuY,nmo,nmo,outfile);
  free_block(MuY);

  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MZ,scratch,noei,0,0,outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = TMP[j][i] = scratch[ij];
    }
  zero_arr(scratch,noei);

  MuZ = block_matrix(nmo,nmo);
  C_DGEMM('n','n',nao,nmo,nao,1,&(TMP[0][0]),nao,&(C[0][0]),nmo,
          0,&(X[0][0]),nmo);
  C_DGEMM('t','n',nmo,nmo,nao,1,&(C[0][0]),nmo,&(X[0][0]),nmo,
          0,&(MuZ[0][0]),nmo);
  
  fprintf(outfile, "\tMO-basis MuZ:\n");
  print_mat(MuZ,nmo,nmo,outfile);
  free_block(MuZ);

  free_block(TMP);
  free_block(X);
  free(scratch);
  
  exit(PSI_RETURN_FAILURE);
  */

  /*
    fprintf(outfile, "\tCanonical MO's:\n");
    print_mat(C,nao,nmo,outfile);
  */


  /* Compute nocc --- closed-shells only */
  for(i=0,nocc=0; i < nirreps; i++) {
    if(openpi[i]) exit(PSI_RETURN_FAILURE);
    nocc += clsdpi[i];
  }

  /*
  for(i=0; i < nocc; i++) fprintf(outfile, "%d %20.15f\n", i, evals[i]);
  */

  fprintf(outfile, "Number of doubly occupied orbitals: %d\n", nocc);

  fprintf(outfile, "\tIter     Pop. Localization   Max. Rotation Angle       Conv\n");
  fprintf(outfile, "\t------------------------------------------------------------\n");

  U = block_matrix(nocc, nocc);
  V = block_matrix(nocc, nocc);
  for(i=0; i < nocc; i++) V[i][i] = 1.0;
  VV = block_matrix(nocc, nocc);

  for(iter=0; iter < 100; iter++) {

    P = 0.0;
    for(i=0; i < nocc; i++) {
      for(A=0; A < natom; A++) {
	PiiA = 0.0;

	for(l=aostart[A]; l <= aostop[A]; l++)
	  for(k=0; k < nao; k++) 
	    PiiA += C[k][i] * C[l][i] * S[k][l];

	P += PiiA * PiiA;
      }
    }

    /* Compute 2x2 rotations to Pipek-Mezey localization */
    alphamax = 0.0;
    for(s=0; s < nocc; s++) {
      for(t=0; t < s; t++) {

	Ast = Bst = 0.0;

	for(A=0; A < natom; A++) {

	  Pst = Pss = Ptt = 0.0;
	      
	  for(l=aostart[A]; l <= aostop[A]; l++) {
	    for(k=0; k < nao; k++) {
	      Pst += 0.5 * (C[k][s] * C[l][t] +
			    C[l][s] * C[k][t]) * S[k][l];

	      Pss += C[k][s] * C[l][s] * S[k][l];

	      Ptt += C[k][t] * C[l][t] * S[k][l];
	    }
	  }

	  Ast += Pst * Pst - 0.25 * (Pss - Ptt) * (Pss - Ptt);
	  Bst += Pst * (Pss - Ptt);

	} /* A-loop */

	/* Compute the rotation angle */
	AB = Ast * Ast + Bst * Bst;
	if(fabs(AB) > 0.0) {
	  cos4a = -Ast/sqrt(AB);
	  alpha = 0.25 * acos(cos4a) * (Bst > 0 ? 1 : -1);
	}
	else alpha = 0.0;

	/* Keep up with the maximum 2x2 rotation angle */
	alphamax = (alpha > alphamax ? alpha : alphamax);

	Uss = cos(alpha);
	Utt = cos(alpha);
	Ust = sin(alpha);
	Uts = -Ust;

	/* Now do the rotation */
	for(k=0; k < nao; k++) {
	  Cks = C[k][s];
	  Ckt = C[k][t];
	  C[k][s] = Uss * Cks + Ust * Ckt;
	  C[k][t] = Uts * Cks + Utt * Ckt;
	}

	zero_mat(U, nocc, nocc);
	for(i=0; i < nocc; i++) U[i][i] = 1.0;

	U[s][s] = Uss;
	U[t][t] = Utt;
	U[s][t] = Ust;
	U[t][s] = Uts;

	zero_mat(VV, nocc, nocc);
	for(i=0; i < nocc; i++) {
	  for(j=0; j < nocc; j++) {
	    for(k=0; k < nocc; k++) {
	      VV[i][j] += V[i][k] * U[j][k];
	    }
	  }
	}

	for(i=0; i < nocc; i++)
	  for(j=0; j < nocc; j++)
	    V[i][j] = VV[i][j];
	      
      } /* t-loop */
    } /* s-loop */

    conv = fabs(alphamax - alphalast);
    fprintf(outfile, "\t%4d  %20.10f  %20.10f  %4.3e\n", iter, P, alphamax, conv);
    if(iter && (conv < 1e-12)) break;
    alphalast = alphamax;

    fflush(outfile);
      
  } /* iter-loop */

  /*  print_mat(V, nocc, nocc, outfile); */

  /* Transform occupied orbital eigenvalues */
  F = block_matrix(nocc, nocc);
  for(i=0; i < nocc; i++)
    for(j=0; j < nocc; j++)
      for(k=0; k < nocc; k++) 
	F[i][j] += V[k][i] * evals[k] * V[k][j];

  fprintf(outfile, "\nTransformed Orbital Energies:\n");
  print_mat(F, nocc, nocc, outfile);

  /* Compute a reordering array based on the diagonal elements of F */
  orb_order = init_int_array(nocc);
  orb_boolean = init_int_array(nocc);
  for(i=0; i < nocc; i++) { orb_order[i] = 0;  orb_boolean[i] = 0; }

  for(i=0,max=0; i < nocc; i++) /* First, find the overall maximum */
    if(fabs(F[i][i]) > fabs(F[max][max])) max = i;

  orb_order[0] = max;  orb_boolean[max] = 1;

  for(i=1; i < nocc; i++) {
    max = 0;
    while(orb_boolean[max]) max++; /* Find an unused max */
    for(j=0; j < nocc; j++) 
      if((fabs(F[j][j]) >= fabs(F[max][max])) && !orb_boolean[j]) max = j;
    orb_order[i] = max; orb_boolean[max] = 1;
  }

  /*
  for(i=0; i < nocc; i++) fprintf(outfile, "%d %d\n", i, orb_order[i]);

  fprintf(outfile, "\n\tPipek-Mezey Localized MO's (before sort):\n");
  print_mat(C, nao, nmo, outfile);
  */

  /* Now reorder the localized MO's according to F */
  Ctmp = block_matrix(nao,nocc);
  for(i=0; i < nocc; i++)
    for(j=0; j < nao; j++) Ctmp[j][i] = C[j][i];

  for(i=0; i < nocc; i++) {
    iold = orb_order[i];
    for(j=0; j < nao; j++) C[j][i] = Ctmp[j][iold];
    evals[i] = F[iold][iold];
  }
  free_block(Ctmp);


  /*
    fprintf(outfile, "\n\tAO Overlap Integrals:\n");
    print_mat(S, nao, nao, outfile);
  */

  /*
  fprintf(outfile, "\n\tPipek-Mezey Localized MO's (after sort):\n");
  print_mat(C, nao, nmo, outfile);
  */

  /* Check MO normalization */
  for(i=0; i < nmo; i++) {
    norm = 0.0;
    for(j=0; j < nao; j++) 
      for(k=0; k < nao; k++) {
	norm += C[j][i] * C[k][i] * S[j][k];
      }

    /*    fprintf(outfile, "norm[%d] = %20.10f\n", i, norm); */
  }

  /* Construct orbital domains and check for completeness */
  charge = init_array(natom);
  rank = init_int_array(natom);
  domain = init_int_matrix(nocc,natom);
  boolean = init_int_array(natom);
  SR = init_array(nao);
  X = block_matrix(nao,nao);
  Z = init_array(nao);

  cutoff = 0.001;  /* Completness cutoff value */
  for(i=0; i < nocc; i++) {

    zero_arr(charge,natom);
    zero_int_array(rank,natom);
    zero_int_array(boolean,natom);

    /* Compute the contribution of each atom to this orbital's charge */
    for(j=0; j < natom; j++) {
      charge[j] = 0.0;
      for(k=aostart[j]; k <= aostop[j]; k++) {
	tmp = 0.0;
	for(l=0; l < nao; l++) tmp += S[k][l] * C[l][i];
	tmp *= C[k][i];
	charge[j] += tmp;
      }
      /*
	fprintf(outfile, "orbital %d charge[%d] = %20.10f\n", i,j,charge[j]);
      */
    }

    /* Now rank the atoms' contributions */
    for(j=0; j < natom; j++) { rank[j] = 0;  boolean[j] = 0; }

    for(j=0,max=0; j < natom; j++) /* First, find the overall maximum */
      if(fabs(charge[j]) > fabs(charge[max])) max = j;

    rank[0] = max;    boolean[max] = 1;

    for(j=1; j < natom; j++) {
      max = 0;
      while(boolean[max]) max++; /* Find an unused max */
      for(k=0; k < natom; k++) 
	if((fabs(charge[k]) >= fabs(charge[max])) && !boolean[k]) max = k;
      rank[j] = max; boolean[max] = 1;
    }

    /*
      for(j=0; j < natom; j++) fprintf(outfile, "%d charge = %20.10f rank = %d\n", j, charge[j], rank[j]);
    */

    /* Now build the domain */
    zero_arr(SR,nao);
    zero_mat(X,nao,nao);
    zero_arr(Z,nao);

    for(j=0; j < nao; j++) 
      for(k=0; k < nao; k++)
	SR[j] += S[j][k] * C[k][i];

    domain[i][rank[0]] = 1; /* by default */

    fR = 1.0;
    next_atom = 1;
    while(fabs(fR) > cutoff) {

      for(j=0,row=0; j < natom; j++) {
	if(domain[i][j]) { /* If this atom is in the domain... */
	  for(k=aostart[j]; k <= aostop[j]; k++,row++) {
	    Z[row] = SR[k];

	    for(l=0,col=0; l < natom; l++) {
	      if(domain[i][l]) {
		for(m=aostart[l]; m <= aostop[l]; m++,col++)
		  X[row][col] = S[k][m];
	      }
	    }
	  }
	}
      }

      /*
	print_mat(X, row, col, outfile);
      */

      /* Solve X * Y = Z */
      /*      stat = pople(X, Z, row, 1, 1e-8, outfile, 0); */  /* row should equal col */
      flin(X, Z, row, 1, &det);

      stat = 0;
      if(!stat) {

	/*
	  fprintf(outfile, "\n\tR' coefficients:\n");
	  for(j=0; j < nao; j++) fprintf(outfile, "%d %20.10f\n", j, Z[j]);
	*/

	/* Now check the completeness of the chosen domain */
	fR = 1.0;
	for(j=0,row=0; j < natom; j++) {
	  if(domain[i][j]) {
	    for(k=aostart[j]; k <= aostop[j]; k++,row++)
	      for(l=0; l < nao; l++) fR -= Z[row] * S[k][l] * C[l][i];
	  }
	}

	/*
	  fprintf(outfile, "\n\tCompleteness check %d = %20.10f\n", i, fR);
	*/

	/* Augment the domain if necessary */
	if(fabs(fR) > cutoff) domain[i][rank[next_atom++]] = 1;
      }
      else fR = 0;
    }

    /* Print out info on the domain */
    if(!stat) {
      fprintf(outfile, "Domain of Occupied Orbital %d:", i);
      for(j=0; j < natom; j++) if(domain[i][j]) fprintf(outfile, "%2d ", j+1);
      fprintf(outfile, "   Completeness = %20.10f\n", fR);
    }
    else fprintf(outfile, "Problem with R' for orbital %d.\n", i);

    fflush(outfile);

  }

  /* Build the SCF closed-shell density matrix/2 */
  /*
    D = block_matrix(nao,nao);
    for(i=0; i < nao; i++) 
    for(j=0; j < nao; j++)
    for(k=0; k < nocc; k++)
    D[i][j] += C[i][k] * C[j][k];

    fprintf(outfile, "\n\tAO-basis SCF Density:\n");
    print_mat(D, nao, nao, outfile);
  */

  /* Compute the virtual space projector */
  /*
    R = block_matrix(nao,nao);
    for(i=0; i < nao; i++) R[i][i] = 1.0;

    C_DGEMM('n','n',nao,nao,nao,-1.0,&(D[0][0]),nao,&(S[0][0]),nao,
    1.0,&(R[0][0]),nao);

    fprintf(outfile, "\n\tVirtual-Space Projector:\n");
    print_mat(R, nao, nao, outfile);
  */

}

}} // namespace psi::cusp
