/*! \defgroup LOCALIZE localize: Localize the orbitals */

/*! 
** \file
** \ingroup LOCALIZE
** \brief Localize the orbitals
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
//#include <libint/libint.h>
#include <libchkpt/chkpt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libqt/qt.h>
#include <psi4-def.h>
#include <psi4-dec.h>

//Default
int LIBINT_MAX_AM = 0;

namespace psi { namespace LOCALIZE {
  using namespace psi;
  void title(void);

int localize(Options & options, int argc, char *argv[])
{
  using namespace psi::LOCALIZE;
  int iter, s, t, A, k, l, m, p, q, inew, iold, max, max_col, phase_ok, phase_chk;
  int i, j, ij, am, atom, shell_length, offset, stat;
  int nirreps, nao, nmo, nso, natom, nshell, noei, nocc, errcod, nfzc;
  int *stype, *snuc, *aostart, *aostop, *ao2atom, *l_length;
  int *clsdpi, *openpi, *orbspi, *dummy, *order, *frdocc;
  double *ss, **S, **scf, **u, **Ctmp, **C, *evals, **MO_S, **scf_old, **X;

  int *orb_order, *orb_boolean, puream;
  double P, PiiA, Pst, Pss, Ptt, Ast, Bst, AB;
  double Uss, Utt, Ust, Uts, Cks, Ckt, **U, **V, **VV, **F;
  double cos4a, alpha, alphamax, alphalast, conv, norm;
  int print;

  alphalast = 1.0;
  title();

  chkpt_init(PSIO_OPEN_OLD);
  puream = chkpt_rd_puream();
  nao = chkpt_rd_nao();
  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();
  natom = chkpt_rd_natom();
  nshell = chkpt_rd_nshell();
  stype = chkpt_rd_stype();
  snuc = chkpt_rd_snuc();
  u = chkpt_rd_usotao();
  nirreps = chkpt_rd_nirreps();
  clsdpi = chkpt_rd_clsdpi();
  openpi = chkpt_rd_openpi();
  orbspi = chkpt_rd_orbspi();
  C = chkpt_rd_scf();
  evals = chkpt_rd_evals();
  chkpt_close();

  /* A couple of error traps */
  if(nirreps != 1) {
    throw PsiException("Error: localization is only valid in C1 symmetry!", __FILE__, __LINE__);
  }
  if(openpi[0]) {
    throw PsiException("Error: localization available for closed-shells only!", __FILE__, __LINE__);
  }

  /* Frozen orbital info */
  frdocc = get_frzcpi();
  nfzc = frdocc[0];
  free(frdocc);

  /* Compute the length of each AM block */
  l_length = init_int_array(LIBINT_MAX_AM);
  l_length[0] = 1;
  for(l=1; l < (LIBINT_MAX_AM); l++) {
    if(puream) l_length[l] = 2 * l + 1;
    else l_length[l] = l_length[l-1] + l + 1;
  }

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

  ao2atom = init_int_array(nso);
  for(i=0; i < natom; i++)
    for(j=aostart[i]; j <= aostop[i]; j++) ao2atom[j] = i;

  /* Get the overlap integrals -- these should be identical to AO S */
  noei = nso*(nso+1)/2;
  ss = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_SO_S,ss,noei,0,0,outfile);
  S = block_matrix(nso,nso);
  for(i=0,ij=0; i < nso; i++)
    for(j=0; j <= i; j++,ij++) {
      S[i][j] = S[j][i] = ss[ij];
    }
  free(ss);

  /* Compute nocc --- closed-shells only */
  for(i=0,nocc=0; i < nirreps; i++) nocc += clsdpi[i];

  fprintf(outfile, "\tNumber of doubly occupied orbitals: %d\n\n", nocc);

  fprintf(outfile, "\tIter     Pop. Localization   Max. Rotation Angle       Conv\n");
  fprintf(outfile, "\t------------------------------------------------------------\n");

  U = block_matrix(nocc, nocc);
  V = block_matrix(nocc, nocc);
  for(i=0; i < nocc; i++) V[i][i] = 1.0;
  VV = block_matrix(nocc, nocc);

  for(iter=0; iter < 100; iter++) {

    P = 0.0;
    for(i=nfzc; i < nocc; i++) {
      for(A=0; A < natom; A++) {
	PiiA = 0.0;

	for(l=aostart[A]; l <= aostop[A]; l++)
	  for(k=0; k < nso; k++) 
	    PiiA += C[k][i] * C[l][i] * S[k][l];

	P += PiiA * PiiA;
      }
    }

    /* Compute 2x2 rotations for Pipek-Mezey localization */
    alphamax = 0.0;
    for(s=nfzc; s < nocc; s++) {
      for(t=nfzc; t < s; t++) {

	Ast = Bst = 0.0;

	for(A=0; A < natom; A++) {

	  Pst = Pss = Ptt = 0.0;
	      
	  for(l=aostart[A]; l <= aostop[A]; l++) {
	    for(k=0; k < nso; k++) {
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
	alphamax = (fabs(alpha) > alphamax ? alpha : alphamax);

	Uss = cos(alpha);
	Utt = cos(alpha);
	Ust = sin(alpha);
	Uts = -Ust;

	/* Now do the rotation */
	for(k=0; k < nso; k++) {
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

    conv = fabs(alphamax) - fabs(alphalast);
    fprintf(outfile, "\t%4d  %20.10f  %20.10f  %4.3e\n", iter, P, alphamax, conv);
    if((iter > 2) && ((fabs(conv) < 1e-12) || alphamax == 0.0)) break;
    alphalast = alphamax;

    fflush(outfile);
      
  } /* iter-loop */

  /*  print_mat(V, nocc, nocc, outfile);  */

  /* Transform occupied orbital eigenvalues */
  F = block_matrix(nocc, nocc);
  for(i=0; i < nocc; i++)
    for(j=0; j < nocc; j++)
      for(k=0; k < nocc; k++) 
	F[i][j] += V[k][i] * evals[k] * V[k][j];

  /*
    fprintf(outfile, "\nTransformed Orbital Energies:\n");
    print_mat(F, nocc, nocc, outfile);
  */

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
  */

  /*
    fprintf(outfile, "\n\tPipek-Mezey Localized MO's (before sort):\n");
    print_mat(C, nso, nmo, outfile);
  */

  /* Now reorder the localized MO's according to F */
  Ctmp = block_matrix(nso,nocc);
  for(i=0; i < nocc; i++)
    for(j=0; j < nso; j++) Ctmp[j][i] = C[j][i];

  for(i=0; i < nocc; i++) {
    iold = orb_order[i];
    for(j=0; j < nso; j++) C[j][i] = Ctmp[j][iold];
    evals[i] = F[iold][iold];
  }
  free_block(Ctmp);

  print = 0;
  print = options.get_bool("PRINT_MOS");
  if(print) {
    fprintf(outfile, "\n\tPipek-Mezey Localized MO's (after sort):\n");
    print_mat(C, nso, nmo, outfile);
  }

  /* Check MO normalization */
  /*
    for(i=0; i < nmo; i++) {
    norm = 0.0;
    for(j=0; j < nso; j++) 
    for(k=0; k < nso; k++) {
    norm += C[j][i] * C[k][i] * S[j][k];
    }

    fprintf(outfile, "norm[%d] = %20.10f\n", i, norm);
    }
  */

  /* correct orbital phases for amplitude restarts */
  chkpt_init(PSIO_OPEN_OLD);
  scf_old = chkpt_rd_local_scf();
  chkpt_close();
  if (scf_old != NULL) {
    MO_S = block_matrix(nmo, nmo);
    X = block_matrix(nso, nmo);
    C_DGEMM('n','n',nso, nmo, nso, 1, &(S[0][0]), nso, &(C[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('t','n',nmo, nmo, nso, 1, &(scf_old[0][0]), nmo, &(X[0][0]), nmo,
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    /*
    fprintf(outfile, "Approximate Overlap Matrix\n");
    print_mat(MO_S, nmo, nmo, outfile);
    */

    for(p=0; p < nmo; p++) {
      max = 0.0;
      for(q=0; q < nmo; q++) {
	if(fabs(MO_S[p][q]) > max) {
	  max = fabs(MO_S[p][q]); max_col = q;
	}
      }
      if(max_col != p) phase_ok = 0;
    }

    chkpt_init(PSIO_OPEN_OLD);
    chkpt_wt_phase_check(phase_ok);
    chkpt_close();
    if(phase_ok) {
      for(p=0; p < nmo; p++) {
	if(MO_S[p][p] < 0.0) {
	  for(q=0; q < nso; q++)
	    C[q][p] *= -1.0;
	}
      }
    }

    free_block(MO_S);
    free_block(scf_old);
  }
  free_block(S);

  /* Write the new MO's to chkpt */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_scf(C);
  chkpt_wt_local_scf(C);
  chkpt_close();

  free_block(C);
  free(evals);

  fprintf(outfile, "\n\tLocalization of occupied orbitals complete.\n");

  exit(PSI_RETURN_SUCCESS);
}
void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*          LOCAL         *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}
}// namespace LOCALIZE
}// namespace psi 


