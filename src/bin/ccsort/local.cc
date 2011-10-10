/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libint/libint.h>
#include <libchkpt/chkpt.h>
#include <liboptions/liboptions.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "Local.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

/*!
** local_init(): Set up parameters of local excitation domains.
**
** The orbital domains constructed here are based on those described
** in Boughton and Pulay, J. Comp. Chem. 14, 736-740 (1993).  The
** localization of the occupied orbitals is done elsewhere (see the
** program "local").  Pair domains are defined as the union of pairs
** of single occupied orbital domains.  "Weak pairs", which are
** defined as pair domains whose individual occupied orbital domains
** have no atoms in common, are identified (cf. int *weak_pairs).
**
** TDC, Jan-June 2002
*/

void domain_print(int, int, int *, int **, double *);
void transpert(const char *);
void sort_pert(const char *, double **, double **, double **, int, int, int);
void build_F_RHF(double);
void build_B_RHF(double);
void cphf_F(const char *);
void cphf_B(const char *);
void local_polar(const char*, int **, int *, int, int *, int *);
void local_magnetic(const char*, int **, int *, int, int *, int *);

void local_init(Options & options)
{
  int i, j, k, ij, stat, a, b, l, I, L;
  int nmo, nso, nocc, nocc_all, nvir, noei, nirreps, nfzc;
  double domain_tot, domain_ave;
  double **C;  /* AO -> localized MO transformation matrix */
  double **Ci; /* localized MO -> AO transformation matrix */
  double **D;  /* 1/2 SCF closed-shell density matrix (AO) */
  double **Rt, **Rt_full; /* Projected, redundant virtual transform (R-tilde) */
  double **S;  /* AO overlap */
  double **St; /* Projected virtual overlap */
  double **Xt; /* Projected, non-redundant virtual transform (X-tilde) */
  double ***V;  /* MO -> projected, redundant virtual transform */
  double **Fmo;/* MO basis Fock matrix */
  double **F;  /* AO basis Fock matrix */
  double **Ft; /* Projected, redundant virtual Fock matrix */
  double **Fbar; /* Projected, non-redundant virtual Fock matrix */
  double ***W;  /* Transformation matrix from tilde -> bar for each ij pair*/
  double *eps_occ; /* occupied orbital energies for local denominators */
  double **eps_vir; /* virtual orbital energies for local denominators */
  double **X, **Y;
  double *evals, **evecs;
  double *eps_all; /* All MO energies */
  dpdfile2 fock;

  int natom, atom, am, offset, nshell, shell_length, next_atom;
  int row, col, max, m, errcod, cnt;
  int *rank, *boolean, *ipiv;
  int *l_length, *aostart, *aostop, *ao2atom;
  int *stype, *snuc;
  int **domain, *domain_len, **pairdomain, *pairdom_len, *pairdom_nrlen;
  int **domain_bp, *domain_len_bp;
  int *weak_pairs;
  double *fR, *charge, *SR, *Z, tmp, *ss;

  int print_test, num_entries, entry_len, orbital;
  int t1_length, t2_length, puream, weak;
  double norm;

  int num_zero;
  double **RS;

  chkpt_init(PSIO_OPEN_OLD);
  C = chkpt_rd_scf();
  natom = chkpt_rd_natom();
  nshell = chkpt_rd_nshell();
  puream = chkpt_rd_puream();
  eps_all = chkpt_rd_evals();
  stype = chkpt_rd_stype();
  snuc = chkpt_rd_snuc();
  chkpt_close();

  timer_on("Local");

  /* C1 symmetry only */
  nirreps = moinfo.nirreps;
  if(nirreps != 1) {
    fprintf(outfile, "\nError: localization must use C1 symmetry.\n");
    exit(PSI_RETURN_FAILURE);
  }

  nso = moinfo.nso;
  nmo = moinfo.nmo; /* should be the same as nso */
  if(nmo != nso) {
    fprintf(outfile, "\nError: NMO != NSO!  %d != %d\n", nmo, nso);
    exit(PSI_RETURN_FAILURE);
  }

  nocc = moinfo.occpi[0]; /* active doubly occupied orbitals */
  nfzc = moinfo.frdocc[0];  /* frozen doubly occupied orbitals */
  nocc_all = nocc + nfzc; /* all doubly occupied orbitals */
  nvir = moinfo.virtpi[0]; /* active virtual orbitals */

  local.nso = nso;
  local.natom = natom;
  local.nocc = nocc;
  local.nvir = nvir;

  /* A couple of scratch arrays */
  X = block_matrix(nso, nso);
  Y = block_matrix(nso, nso);

  /* Invert C */
  Ci = block_matrix(nso, nso);
  for(i=0; i < nso; i++)
    for(j=0; j < nso; j++)
      Y[i][j] = C[i][j];

  invert_matrix(C, Ci, nso, outfile);

  for(i=0; i < nso; i++)
    for(j=0; j < nso; j++)
      C[i][j] = Y[i][j];

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

  /*
    fprintf(outfile, "\n\tAO Overlap (S)\n");
    print_mat(S, nso, nso, outfile);
  */

  /* Build the SCF closed-shell density matrix/2 */
  D = block_matrix(nso,nso);
  for(i=0; i < nso; i++)
    for(j=0; j < nso; j++)
      for(k=0; k < nocc_all; k++)
        D[i][j] += C[i][k] * C[j][k];

  /*
    fprintf(outfile, "\n\tAO-basis SCF Density (D):\n");
    print_mat(D, nso, nso, outfile);
  */


  /* Compute the length of each AM block */
  l_length = init_int_array(LIBINT_MAX_AM);
  l_length[0] = 1;
  for(l=1; l < LIBINT_MAX_AM; l++) {
    if(puream) l_length[l] = 2 * l + 1;
    else l_length[l] = l_length[l-1] + l + 1;
  }

  /* Set up the atom->AO and AO->atom lookups */
  aostart = init_int_array(natom);
  aostop = init_int_array(natom);
  for(i=0,atom=-1,offset=0; i<nshell; i++) {
    am = stype[i] - 1;
    shell_length = l_length[am];

    if(atom != snuc[i] - 1) {
      if(atom != -1) aostop[atom] = offset-1;
      atom = snuc[i]-1;
      aostart[atom] = offset;
    }
    offset += shell_length;
  }
  aostop[atom] = offset-1;

  ao2atom = init_int_array(nso);
  for(i=0; i < natom; i++)
    for(j=aostart[i]; j < aostop[i]; j++) ao2atom[j] = i;

  free(stype);
  free(snuc);

  /************* Build the orbital domains ************/

  domain = init_int_matrix(nocc, natom);
  domain_len = init_int_array(nocc);
  domain_bp = init_int_matrix(nocc, natom);
  domain_len_bp = init_int_array(nocc);
  charge = init_array(natom);
  rank = init_int_array(natom);
  boolean = init_int_array(natom);
  SR = init_array(nso);
  Z = init_array(nso);
  ipiv = init_int_array(nso);
  fR = init_array(nocc);

  for(i=nfzc; i < nocc_all; i++) {

    /* Compute the contribution of each atom to this orbital's charge/population */
    for(j=0; j < natom; j++) {
      charge[j] = 0.0;
      for(k=aostart[j]; k <= aostop[j]; k++) {
        tmp = 0.0;
        for(l=0; l < nso; l++) tmp += S[k][l] * C[l][i];
        tmp *= C[k][i];
        charge[j] += tmp;
      }
    }

    /* Rank the atomic contributions to the orbital's charge */
    for(j=0; j < natom; j++) { rank[j] = 0; boolean[j] = 0; }
    for(j=0,max=0; j < natom; j++) /* find the overall maximum */
      if(fabs(charge[j]) >= fabs(charge[max])) max = j;
    rank[0] = max; boolean[max] = 1;
    for(j=1; j < natom; j++) {
      max = 0;
      while(boolean[max]) max++; /* find an unused max */
      for(k=0; k < natom; k++)
        if((fabs(charge[k]) >= fabs(charge[max])) && !boolean[k]) max = k;
      rank[j] = max; boolean[max] = 1;
    }

    /* Build the orbital's domain starting in order of decreasing charge contribution */
    for(j=0; j < nso; j++) {
      SR[j] = 0.0;
      for(k=0; k < nso; k++)
        SR[j] += S[j][k] * C[k][i];
    }

    domain[i-nfzc][rank[0]] = 1; /* at least one atom must be in the domain */
    domain_len[i-nfzc] = 1;

    fR[i-nfzc] = 1.0;
    next_atom = 1;
    while(fabs(fR[i-nfzc]) > local.cutoff) {

      /* Completeness check */
      for(j=0,row=0; j < natom; j++) {
        if(domain[i-nfzc][j]) {
          for(k=aostart[j]; k <= aostop[j]; k++,row++) {

            Z[row] = SR[k];

            for(l=0,col=0; l < natom; l++) {
              if(domain[i-nfzc][l]) {

                for(m=aostart[l]; m <= aostop[l]; m++,col++)
                  X[row][col] = S[k][m];

              }
            } /* l */

          } /* k */
        }
      } /* j */

      errcod = C_DGESV(row, 1, &(X[0][0]), nso, &(ipiv[0]), &(Z[0]), nso);
      if(errcod) {
        fprintf(outfile, "\nError in DGESV return in orbital domain construction.\n");
        exit(PSI_RETURN_FAILURE);
      }

      fR[i-nfzc] = 1.0;
      for(j=0,row=0; j < natom; j++) {
        if(domain[i-nfzc][j]) {
          for(k=aostart[j]; k <= aostop[j]; k++,row++) {
            for(l=0; l < nso; l++) fR[i-nfzc] -= Z[row] * S[k][l] * C[l][i];
          }
        }
      }

      /* Augment the domain if necessary */
      if(fabs(fR[i-nfzc]) > local.cutoff) {
        domain[i-nfzc][rank[next_atom++]] = 1;
        domain_len[i-nfzc]++;
      }
    } /* cutoff check */
  } /* i */

  for(i=0; i<nocc; i++) {
    domain_len_bp[i] = domain_len[i];
    for(k=0; k<natom; k++)
      domain_bp[i][k] = domain[i][k];
  }
  /* Print the orbital domains */
  fprintf(outfile, "\n   ****** Boughton-Pulay Occupied Orbital Domains ******\n");
  domain_print(nocc, natom, domain_len_bp, domain_bp, fR);

  /* Identify and/or remove weak pairs -- using Bougton-Pulay domains */
  weak_pairs = init_int_array(nocc*nocc);
  if(local.pairdef=="BP") {
    fprintf(outfile, "\n");
    for(i=0,ij=0; i < nocc; i++)
      for(j=0; j < nocc; j++,ij++) {
        weak = 1;
        for(k=0; k < natom; k++)
          if(domain[i][k] && domain[j][k]) weak = 0;

        if(weak && local.weakp!="NONE") {
          weak_pairs[ij] = 1;

          if(local.weakp=="MP2")
            fprintf(outfile, "\tPair %d %d [%d] is weak and will be treated with MP2.\n", i, j, ij);
          else if(local.weakp=="NEGLECT") {
            fprintf(outfile, "\tPair %d %d = [%d] is weak and will be deleted.\n", i, j, ij);
          }
        }
        else weak_pairs[ij] = 0;
      }
  }

  /* If this is a response calculation, augment domains using polarized orbitals */
  fprintf(outfile, "\n");
  if(local.domain_polar) {
    fprintf(outfile, "\tGenerating electric-field CPHF solutions for local-CC.\n");
    fflush(outfile);
    transpert("Mu");
    sort_pert("Mu", moinfo.MUX, moinfo.MUY, moinfo.MUZ,
              moinfo.irrep_x, moinfo.irrep_y, moinfo.irrep_z);

    if(local.domain_sep) {
      /* Zero omega */
      build_F_RHF(0);
      cphf_F("X");
      local_polar("X", domain, domain_len, natom, aostart, aostop);
      fprintf(outfile, "\nMu_X (%d)\n", 0);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_F("Y");
      local_polar("Y", domain, domain_len, natom, aostart, aostop);
      fprintf(outfile, "\nMu_Y (%d)\n", 0);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_F("Z");
      local_polar("Z", domain, domain_len, natom, aostart, aostop);
      fprintf(outfile, "\nMu_Z (%d)\n", 0);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      /* Positive omega */
      build_F_RHF(params.omega[0]);
      cphf_F("X");
      local_polar("X", domain, domain_len, natom,	aostart, aostop);
      fprintf(outfile, "\nMu_X (%lf)\n", params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_F("Y");
      local_polar("Y", domain, domain_len, natom,
                  aostart, aostop);
      fprintf(outfile, "\nMu_Y (%lf)\n", params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_F("Z");
      local_polar("Z", domain, domain_len, natom,	aostart, aostop);
      fprintf(outfile, "\nMu_Z (%lf)\n", params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      /* Negative omega */
      build_F_RHF(-params.omega[0]);
      cphf_F("X");
      local_polar("X", domain, domain_len, natom,	aostart, aostop);
      fprintf(outfile, "\nMu_X (%lf)\n", -params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_F("Y");
      local_polar("Y", domain, domain_len, natom,
                  aostart, aostop);
      fprintf(outfile, "\nMu_Y (%lf)\n", -params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_F("Z");
      local_polar("Z", domain, domain_len, natom,	aostart, aostop);
      fprintf(outfile, "\nMu_Z (%lf)\n", -params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }
    }
    else {
      build_F_RHF(0);
      cphf_F("X");
      local_polar("X", domain, domain_len, natom,	aostart, aostop);
      cphf_F("Y");
      local_polar("Y", domain, domain_len, natom,	aostart, aostop);
      cphf_F("Z");
      local_polar("Z", domain, domain_len, natom,	aostart, aostop);

    }
  }
  if(local.domain_mag) {
    fprintf(outfile, "\tGenerating magnetic-field CPHF solutions for local-CC.\n");
    fflush(outfile);
    transpert("L");
    sort_pert("L", moinfo.LX, moinfo.LY, moinfo.LZ,
              moinfo.irrep_x, moinfo.irrep_y, moinfo.irrep_z);

    if(local.domain_sep) {
      /* Zero omega */
      build_B_RHF(0);
      cphf_B("X");
      local_magnetic("X", domain, domain_len, natom, aostart, aostop);
      fprintf(outfile, "\nL_X (%d)\n", 0);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_B("Y");
      local_magnetic("Y", domain, domain_len, natom,
                     aostart, aostop);
      fprintf(outfile, "\nL_Y (%d)\n", 0);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_B("Z");
      local_magnetic("Z", domain, domain_len, natom, aostart, aostop);
      fprintf(outfile, "\nL_Z (%d)\n", 0);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      /* Positive omega */
      build_B_RHF(params.omega[0]);
      cphf_B("X");
      local_magnetic("X", domain, domain_len, natom, aostart, aostop);
      fprintf(outfile, "\nL_X (%lf)\n", params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_B("Y");
      local_magnetic("Y", domain, domain_len, natom, aostart, aostop);
      fprintf(outfile, "\nL_Y (%lf)\n", params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_B("Z");
      local_magnetic("Z", domain, domain_len, natom, aostart, aostop);
      fprintf(outfile, "\nL_Z (%lf)\n", params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      /* Negative omega */
      build_B_RHF(-params.omega[0]);
      cphf_B("X");
      local_magnetic("X", domain, domain_len, natom, aostart, aostop);
      fprintf(outfile, "\nL_X (%lf)\n", -params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_B("Y");
      local_magnetic("Y", domain, domain_len, natom, aostart, aostop);
      fprintf(outfile, "\nL_Y (%lf)\n", -params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }

      cphf_B("Z");
      local_magnetic("Z", domain, domain_len, natom, aostart, aostop);
      fprintf(outfile, "\nL_Z (%lf)\n", -params.omega[0]);
      domain_print(nocc, natom, domain_len, domain, fR);
      for(i=0; i<nocc; i++) {
        domain_len[i] = domain_len_bp[i];
        for(k=0; k<natom; k++)
          domain[i][k] = domain_bp[i][k];
      }
    }
    else {
      build_B_RHF(0);
      cphf_B("X");
      local_magnetic("X", domain, domain_len, natom, aostart, aostop);
      cphf_B("Y");
      local_magnetic("Y", domain, domain_len, natom, aostart, aostop);
      cphf_B("Z");
      local_magnetic("Z", domain, domain_len, natom, aostart, aostop);
    }
  }

  /* Allow user input of selected domains */
  num_entries = options["DOMAINS"].size();
  for(i=0; i < num_entries; i++) {
    entry_len = options["DOMAINS"][i].size();
    orbital = options["DOMAINS"][i][0].to_integer();

    /* Clear out the current domain for this orbital */
    for(j=0; j < natom; j++) domain[orbital][j] = 0;
    domain_len[orbital] = 0;

    for(j=1; j < entry_len; j++) {
      atom = options["DOMAINS"][i][j].to_integer();
      domain[orbital][atom] = 1;
      domain_len[orbital]++;
    }
  }

  /* Recheck Completeness */
  for(i=nfzc; i < nocc_all; i++) {

    /* Build the orbital's domain starting in order of decreasing charge contribution */
    for(j=0; j < nso; j++) {
      SR[j] = 0.0;
      for(k=0; k < nso; k++)
        SR[j] += S[j][k] * C[k][i];
    }

    for(j=0,row=0; j < natom; j++) {
      if(domain[i-nfzc][j]) {
        for(k=aostart[j]; k <= aostop[j]; k++,row++) {

          Z[row] = SR[k];

          for(l=0,col=0; l < natom; l++) {
            if(domain[i-nfzc][l]) {

              for(m=aostart[l]; m <= aostop[l]; m++,col++)
                X[row][col] = S[k][m];

            }
          } /* l */

        } /* k */
      }
    } /* j */

    errcod = C_DGESV(row, 1, &(X[0][0]), nso, &(ipiv[0]), &(Z[0]), nso);
    if(errcod) {
      fprintf(outfile, "\nError in DGESV return in orbital domain construction.\n");
      exit(PSI_RETURN_FAILURE);
    }

    fR[i-nfzc] = 1.0;
    for(j=0,row=0; j < natom; j++) {
      if(domain[i-nfzc][j]) {
        for(k=aostart[j]; k <= aostop[j]; k++,row++) {
          for(l=0; l < nso; l++) fR[i-nfzc] -= Z[row] * S[k][l] * C[l][i];
        }
      }
    }

  } /* i */


  /* Print the orbital domains */
  fprintf(outfile, "\n   ****** Final Occupied Orbital Domains ******\n");
  if(!local.domain_sep)
    domain_print(nocc, natom, domain_len, domain, fR);

  /* Build the pair domains */
  pairdomain = init_int_matrix(nocc*nocc,natom);
  pairdom_len = init_int_array(nocc*nocc);
  for(i=0,ij=0; i < nocc; i++)
    for(j=0; j < nocc; j++,ij++)
      for(k=0; k < natom; k++) {
        if(domain[i][k] || domain[j][k]) {
          pairdomain[ij][k] = 1;
          pairdom_len[ij] += aostop[k] - aostart[k] + 1;
        }
      }

  /* Identify and/or remove weak pairs -- for CPHF "response" domains */
  if(local.domain_polar || local.domain_mag) {
    fprintf(outfile, "\n");
    for(i=0,ij=0; i < nocc; i++)
      for(j=0; j < nocc; j++,ij++) {
        weak = 1;
        for(k=0; k < natom; k++)
          if(domain[i][k] && domain[j][k]) weak = 0;

        if(weak && local.weakp!="NONE") {
          weak_pairs[ij] = 1;

          if(local.weakp=="MP2")
            fprintf(outfile, "\tPair %d %d [%d] is weak and will be treated with MP2.\n", i, j, ij);
          else if(local.weakp=="NEGLECT") {
            fprintf(outfile, "\tPair %d %d = [%d] is weak and will be deleted.\n", i, j, ij);
          }
        }
        else weak_pairs[ij] = 0;
      }
  }

  /* Compute the total number of singles and doubles */
  /* replacement code from TDC on 11-5-02 */
  t1_length = t2_length = 0;
  for(i=0,ij=0; i < nocc; i++) {
    for(k=0; k < natom; k++) {
      if(domain[i][k])
        for(a=aostart[k]; a <= aostop[k]; a++) t1_length++;
    }
    for(j=0; j < nocc; j++,ij++) {
      for(k=0; k < natom; k++) {
        for(l=0; l < natom; l++) {
          if(pairdomain[ij][k] && pairdomain[ij][l] && !weak_pairs[ij]) {
            for(a=aostart[k]; a <= aostop[k]; a++)
              for(b=aostart[l]; b <= aostop[l]; b++)
                t2_length++;
          }
        }
      }
    }
  }

  /* Print excitation space reduction info */
  fprintf(outfile, "\n\tT1 Length = %d (local), %d (canonical)\n",
          t1_length, nocc*nvir);
  fprintf(outfile, "\tT2 Length = %d (local), %d (canonical)\n\n",
          t2_length, nocc*nocc*nvir*nvir);
  fflush(outfile);

  local.domain = domain;
  local.domain_len = domain_len;
  local.pairdomain = pairdomain;
  local.pairdom_len = pairdom_len;
  local.weak_pairs = weak_pairs;
  local.aostart = aostart;
  local.aostop = aostop;

  free_int_matrix(domain_bp);
  free(domain_len_bp);

  free(ao2atom);
  free(l_length);
  free(charge);
  free(rank);
  free(boolean);
  free(SR);
  free(Z);
  free(ipiv);
  free(fR);

  print_test = 0;
  print_test = options.get_bool("DOMAIN_PRINT");
  if(print_test) {
    fprintf(outfile, "Printing of orbital domains requested...exiting.\n\n");
    exit(PSI_RETURN_FAILURE);
  }

  /************* Orbital Domains Complete ***************/

  /* Compute the complete virtual space projector */
  Rt_full = block_matrix(nso,nso);
  for(i=0; i < nso; i++) Rt_full[i][i] = 1.0;

  C_DGEMM('n','n',nso,nso,nso,-1.0,&(D[0][0]),nso,&(S[0][0]),nso,
          1.0,&(Rt_full[0][0]),nso);

  /*
    fprintf(outfile, "\n\tVirtual-Space Projector (R-tilde):\n");
    print_mat(Rt_full, nso, nso, stdout);
  */

  /* Compute the norm of each PAO */
  for(i=0; i < nso; i++) {
    norm = 0.0;
    for(j=0; j < nso; j++) {
      norm += Rt_full[j][i] * Rt_full[j][i];
    }
    norm = sqrt(norm);
    if(norm < local.core_cutoff && local.freeze_core!="FALSE") {
      fprintf(outfile, "\tNorm of orbital %4d = %20.12f...deleteing\n", i, norm);
      for(j=0; j < nso; j++) Rt_full[j][i] = 0.0;
    }
  }
  fprintf(outfile, "\n");
  fflush(outfile);

  /* Grab the MO-basis Fock matrix */
  Fmo = block_matrix(nso, nso);
  for(i=0; i < nfzc; i++) Fmo[i][i] = eps_all[i];
  dpd_file2_init(&fock, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&fock);
  dpd_file2_mat_rd(&fock);
  for(i=0; i < nocc; i++)
    for(j=0; j < nocc; j++)
      Fmo[i+nfzc][j+nfzc] = fock.matrix[0][i][j];
  dpd_file2_mat_close(&fock);
  dpd_file2_close(&fock);

  dpd_file2_init(&fock, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fock);
  dpd_file2_mat_rd(&fock);
  for(i=0; i < nvir; i++)
    for(j=0; j < nvir; j++)
      Fmo[i+nfzc+nocc][j+nfzc+nocc] = fock.matrix[0][i][j];
  dpd_file2_mat_close(&fock);
  dpd_file2_close(&fock);

  /*
    fprintf(outfile, "\n\tMO Basis Fock matrix:\n");
    print_mat(Fmo, nso, nso, outfile);
  */

  /* Build the AO-basis Fock matrix */
  F = block_matrix(nso,nso);
  C_DGEMM('t','n',nso,nso,nso,1.0,&(Ci[0][0]),nso,&(Fmo[0][0]),nso,
          0.0,&(X[0][0]),nso);
  C_DGEMM('n','n',nso,nso,nso,1.0,&(X[0][0]),nso,&(Ci[0][0]),nso,
          0.0,&(F[0][0]),nso);

  /* Build the occupied orbital energy list */
  eps_occ = init_array(nocc);
  for(i=0;i < nocc; i++) eps_occ[i] = Fmo[i+nfzc][i+nfzc];

  /*
    fprintf(outfile, "\n\tAO-Basis Fock Matrix:\n");
    print_mat(Fmo, nso, nso, outfile);
  */

  /* Compute R^+ S for virtual orbitals */
  RS = block_matrix(nvir,nso);
  for(a=0; a < nvir; a++)
    for(i=0; i < nso; i++)
      X[i][a] = C[i][a+nocc_all];
  C_DGEMM('t','n',nvir,nso,nso,1.0,&(X[0][0]),nso,&(S[0][0]),nso,0.0,&(RS[0][0]),nso);

  /* Build the virtual metric and W transforms for each pair domain */
  Rt = block_matrix(nso, nso);
  W = (double ***) malloc(nocc * nocc * sizeof(double **));
  V = (double ***) malloc(nocc * nocc * sizeof(double **));
  eps_vir = (double **) malloc(nocc * nocc * sizeof(double *));
  pairdom_nrlen = init_int_array(nocc * nocc); /* dimension of non-redundant basis */
  num_zero = 0;

  for(ij=0; ij < nocc * nocc; ij++) {

    zero_mat(Rt, nso, nso);

    /* Build the virtual space projector for this pair */
    for(k=0,L=0; k < natom; k++) {
      if(pairdomain[ij][k]) {
        for(l=aostart[k]; l <= aostop[k]; l++,L++) {
          for(m=0; m < nso; m++) {
            Rt[m][L] = Rt_full[m][l];
          }
        }
      }
    }

    /* Compute the MO -> projected virtual transformation matrix */
    V[ij] = block_matrix(nvir,pairdom_len[ij]);
    C_DGEMM('n','n',nvir,pairdom_len[ij],nso,1.0,&(RS[0][0]),nso,&(Rt[0][0]),nso,0.0,
            &(V[ij][0][0]),pairdom_len[ij]);

    /*
      fprintf(outfile, "\nV[%d]:\n", ij);
      fprintf(outfile,   "======\n");
      print_mat(V[ij], nvir, pairdom_len[ij], outfile);
    */

    /* Virtual space metric */
    St = block_matrix(pairdom_len[ij],pairdom_len[ij]);
    C_DGEMM('n','n',nso,pairdom_len[ij],nso,1.0,&(S[0][0]),nso,&(Rt[0][0]),nso,
            0.0,&(X[0][0]),nso);
    C_DGEMM('t','n',pairdom_len[ij],pairdom_len[ij],nso,1.0,&(Rt[0][0]),nso,&(X[0][0]),nso,
            0.0,&(St[0][0]),pairdom_len[ij]);

    /*
      fprintf(outfile, "\n\tVirtual-Space Metric (S-tilde) for ij = %d:\n", ij);
      print_mat(St, pairdom_len[ij], pairdom_len[ij], outfile);
    */

    /* Diagonalize metric */
    evals = init_array(pairdom_len[ij]);
    evecs = block_matrix(pairdom_len[ij],pairdom_len[ij]);
    sq_rsp(pairdom_len[ij],pairdom_len[ij],St,evals,1,evecs,1e-12);

    /* Count the number of zero eigenvalues */
    for(i=0,cnt=0; i < pairdom_len[ij]; i++) if(evals[i] <= 1e-6) cnt++;

    pairdom_nrlen[ij] = pairdom_len[ij]-cnt;

    /*
      fprintf(outfile, "\n\tS-tilde eigenvalues for ij = %d:\n", ij);
      for(i=0; i < pairdom_len[ij]; i++) fprintf(outfile, "\t%d %20.12f\n", i, evals[i]);

      fprintf(outfile, "\n\tS-tilde eigenvectors for ij = %d:\n", ij);
      print_mat(evecs,pairdom_len[ij],pairdom_len[ij],outfile);
    */

    /* Build the projected, non-redundant transform (X-tilde) */
    Xt = block_matrix(pairdom_len[ij],pairdom_nrlen[ij]);
    for(i=0,I=0; i < pairdom_len[ij]; i++) {
      if(evals[i] > 1e-6) {
        for(j=0; j < pairdom_len[ij]; j++)
          Xt[j][I] = evecs[j][i]/sqrt(evals[i]);
        I++;
      }
      else num_zero++;
    }


    /*
      fprintf(outfile, "\n\tTransform to non-redundant, projected virtuals (X-tilde) for ij = %d:\n", ij);
      print_mat(Xt, pairdom_len[ij], pairdom_nrlen[ij], outfile);
    */

    free_block(evecs);
    free(evals);

    /* Build the projected (redundant) virtual Fock matrix */
    Ft = block_matrix(pairdom_len[ij], pairdom_len[ij]);
    C_DGEMM('t','n',pairdom_len[ij],nso,nso,1.0,&(Rt[0][0]),nso,&(F[0][0]),nso,
            0.0,&(X[0][0]),nso);
    C_DGEMM('n','n',pairdom_len[ij],pairdom_len[ij],nso,1.0,&(X[0][0]),nso,&(Rt[0][0]),nso,
            0.0,&(Ft[0][0]),pairdom_len[ij]);

    /* Project the Fock matrix into the non-redundant virtual space */
    Fbar = block_matrix(pairdom_nrlen[ij],pairdom_nrlen[ij]);
    C_DGEMM('t','n',pairdom_nrlen[ij],pairdom_len[ij],pairdom_len[ij],1.0,
            &(Xt[0][0]),pairdom_nrlen[ij],&(Ft[0][0]),pairdom_len[ij],0.0,&(X[0][0]),nso);
    C_DGEMM('n','n',pairdom_nrlen[ij],pairdom_nrlen[ij],pairdom_len[ij],1.0,
            &(X[0][0]),nso,&(Xt[0][0]),pairdom_nrlen[ij],0.0,&(Fbar[0][0]),pairdom_nrlen[ij]);

    /*
      fprintf(outfile, "\n\tFbar matrix for ij = %d:\n", ij);
      print_mat(Fbar,pairdom_nrlen[ij],pairdom_nrlen[ij],outfile);
    */

    /* Diagonalize Fbar */
    evals = init_array(pairdom_nrlen[ij]);
    evecs = block_matrix(pairdom_nrlen[ij],pairdom_nrlen[ij]);
    sq_rsp(pairdom_nrlen[ij],pairdom_nrlen[ij],Fbar,evals,1,evecs,1e-12);

    /*
      fprintf(stdout, "\n\tFbar eigenvectors for ij = %d:\n", ij);
      print_mat(evecs,pairdom_nrlen[ij],pairdom_nrlen[ij],outfile);
    */

    /* Finally, build the W matrix */
    W[ij] = block_matrix(pairdom_len[ij],pairdom_nrlen[ij]);
    C_DGEMM('n','n',pairdom_len[ij],pairdom_nrlen[ij],pairdom_nrlen[ij],1.0,
            &(Xt[0][0]),pairdom_nrlen[ij],&(evecs[0][0]),pairdom_nrlen[ij],
            0.0,&(W[ij][0][0]),pairdom_nrlen[ij]);

    /*
      fprintf(outfile, "\n\tW Transformation Matrix for ij = %d:\n", ij);
      print_mat(W[ij],pairdom_len[ij],pairdom_nrlen[ij],outfile);
    */

    /* build the orbital energy list */
    eps_vir[ij] = init_array(pairdom_nrlen[ij]);
    for(i=0; i < pairdom_nrlen[ij]; i++)
      eps_vir[ij][i] = evals[i]; /* virtual orbital energies */

    /*
      fprintf(outfile, "\n\tVirtual orbital Energies for ij = %d:\n", ij);
      for(i=0; i < pairdom_nrlen[ij]; i++)
      fprintf(outfile, "%d %20.12f\n", i, eps_vir[ij][i]);
    */

    free(evals);
    free_block(evecs);

    free_block(St);
    free_block(Xt);
    free_block(Fbar);
    free(Ft);

  } /* ij loop */

  free_block(RS);
  free_block(F);
  free_block(Fmo);
  free_block(S);
  free_block(Rt_full);
  free_block(D);
  free_block(Ci);
  free_block(C);

  free_block(X);
  free_block(Y);

  local.W = W;
  local.V = V;
  local.eps_occ = eps_occ;
  local.eps_vir = eps_vir;
  local.pairdom_nrlen = pairdom_nrlen;

  local.weak_pair_energy = 0.0;

  fflush(outfile);
  timer_off("Local");
}

void local_done(void)
{
  int i, ij, h, nocc, nvir, natom;
  psio_address next;

  nocc = local.nocc;
  nvir = local.nvir;
  natom = local.natom;

  psio_write_entry(CC_INFO, "Local Cutoff", (char *) &local.cutoff,
                   sizeof(double));
  psio_write_entry(CC_INFO, "Local Domain Length", (char *) local.domain_len,
                   nocc*sizeof(int));
  psio_write_entry(CC_INFO, "Local Pair Domain Length", (char *) local.pairdom_len,
                   nocc*nocc*sizeof(int));
  psio_write_entry(CC_INFO, "Local Pair Domain NR Length", (char *) local.pairdom_nrlen,
                   nocc*nocc*sizeof(int));
  psio_write_entry(CC_INFO, "Local Weak Pairs", (char *) local.weak_pairs,
                   nocc*nocc*sizeof(int));
  psio_write_entry(CC_INFO, "Local Occupied Orbital Energies", (char *) local.eps_occ,
                   nocc*sizeof(double));

  next = PSIO_ZERO;
  for(i=0; i<nocc; i++)
    psio_write(CC_INFO, "Local Domains", (char *) local.domain[i],
               natom*sizeof(int), next, &next);
  next = PSIO_ZERO;
  for(ij=0; ij<nocc*nocc; ij++)
    psio_write(CC_INFO, "Local Pair Domains", (char *) local.pairdomain[ij],
               natom*sizeof(int), next, &next);
  next = PSIO_ZERO;
  for(ij=0; ij < nocc*nocc; ij++)
      psio_write(CC_INFO, "Local Virtual Orbital Energies", (char *) local.eps_vir[ij],
                 local.pairdom_nrlen[ij]*sizeof(double), next, &next);
  next = PSIO_ZERO;
  for(ij=0; ij < nocc*nocc; ij++)
      psio_write(CC_INFO, "Local Transformation Matrix (W)", (char *) local.W[ij][0],
                 local.pairdom_len[ij]*local.pairdom_nrlen[ij]*sizeof(double), next, &next);
  next = PSIO_ZERO;
  for(ij=0; ij < nocc*nocc; ij++)
      psio_write(CC_INFO, "Local Residual Vector (V)", (char *) local.V[ij][0],
                 nvir*local.pairdom_len[ij]*sizeof(double), next, &next);

  if(params.ref == 0 || params.ref == 1) {
    for(h=0; h < moinfo.nirreps; h++)
      if(moinfo.sopi[h] && moinfo.virtpi[h]) free_block(moinfo.C[h]);
    free(moinfo.C);
  }

  free_int_matrix(local.pairdomain);
  free_int_matrix(local.domain);
  for(i=0; i < nocc*nocc; i++) {
    free_block(local.W[i]);
    free_block(local.V[i]);
    free(local.eps_vir[i]);
  }
  free(local.W);
  free(local.V);
  free(local.eps_vir);

  free(local.aostart);
  free(local.aostop);

  free(local.eps_occ);
  free(local.domain_len);
  free(local.pairdom_len);
  free(local.pairdom_nrlen);
  free(local.weak_pairs);
}

}} // namespace psi::ccsort
