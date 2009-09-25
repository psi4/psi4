/*! \file
    \ingroup MVO
    \brief Enter brief description of file here 
*/
/*
**  UMP2NOs
**  
**  Code to generate Unrestricted MP2 One-PDMs and Natural Orbitals
**
**  C. David Sherrill and Gerry O. Hyde
**  Georgia Institute of Technology
**  March 24, 2002
**
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <libiwl/iwl.h>
#include "MOInfo.h"
#include "params.h"
#include "globals.h"


/* First definitions of globals */
extern "C" {
  extern FILE *infile, *outfile;
}

namespace psi { namespace mvo {

extern int *ioff;
extern struct MOInfo moinfo;
extern struct Params params;

/* Max length of ioff array */
#define IOFF_MAX 32641
#define TOL 1E-14

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

extern void transform_density(double **onepdm, double **psq_so, int spin);

void get_mp2nos(void)
{
  double *evals_alpha, *evals_beta;
  double *Iaaaa, *Ibbbb, **Iaabb;
  double *taaaa, *tbbbb, *taabb;
  double **Paa, **Pbb;
  double **P_so_aa, **P_so_bb, **P_so_tot, **P_mo_tot, **tmat;
  double *P_eigvals, **P_eigvecs, **P_mo_block, **P_so_block, *Stri, **Smat;
  double **scfvec, **scfvec_irrep;
  double energy, E_scf, tval;

  int nmo, occ, vir, ntri, irrep, nirreps, irrep_dim;
  int printflg = 0, errcod;
  int i,j,i_ci,j_ci,mo_offset;

  double *evals_unsrt;
  struct iwlbuf Buf;  
  void get_ampl(double *ampl, double *ints, 
                double *evals1, double *evals2, int mixed);
  void get_dens_alpha(double **Paa, double *taaaa, double *taabb);
  void get_dens_beta(double **Pbb, double *tbbbb, double *taabb);
  double compute_energy(double *taaaa, double *tbbbb, double *taabb, 
                        double *Iaaaa, double *Ibbbb, double **Iaabb);
  double check_density_trace(double **P);


  nmo = moinfo.nmo;
  occ = moinfo.ndocc;
  vir = nmo - occ;
  ntri = (nmo * (nmo+1))/2; 

  /* read the eigenvalues */
  evals_alpha = init_array(nmo); 
  evals_beta = init_array(nmo); 
  chkpt_init(PSIO_OPEN_OLD);
  E_scf = chkpt_rd_escf();
  evals_unsrt = chkpt_rd_alpha_evals();

  /* we'll need "restricted" eigenvalues in the checkpoint file later on.
   * let's just put in the alpha ones, we don't use them in DETCI for
   * anything but preconditioners and/or guessers, but we need something
   * there 
   */
  chkpt_wt_evals(evals_unsrt);

  for (i=0; i<nmo; i++) {
    j = moinfo.order[i];
    evals_alpha[j] = evals_unsrt[i];
    }
  free(evals_unsrt);

  evals_unsrt = chkpt_rd_beta_evals();
  for (i=0; i<nmo; i++) {
    j = moinfo.order[i];
    evals_beta[j] = evals_unsrt[i];
    }
  free(evals_unsrt);


  /* allocate memory for integrals */
  Iaaaa = init_array((ntri*(ntri+1))/2);
  Ibbbb = init_array((ntri*(ntri+1))/2);
  Iaabb = block_matrix(ntri, ntri);

  /* allocate memory for amplitudes */
  taaaa = init_array(((occ*(occ+1))/2)*((vir*(vir+1))/2));
  tbbbb = init_array(((occ*(occ+1))/2)*((vir*(vir+1))/2));
  taabb = init_array((occ*occ)*(vir*vir));

  /* read the integrals */
  iwl_rdtwo(PSIF_MO_AA_TEI, Iaaaa, ioff, nmo, 0, 0, printflg, outfile); 
  iwl_rdtwo(PSIF_MO_BB_TEI, Ibbbb, ioff, nmo, 0, 0, printflg, outfile);

  /* mixed spin case requires special read routine since no pq/rs symm */ 
  iwl_buf_init(&Buf, PSIF_MO_AB_TEI, 0.0, 1, 1);
  iwl_buf_rd_all2(&Buf, Iaabb, ioff, ioff, 0, ioff, printflg, outfile);
  iwl_buf_close(&Buf, 1);

  /* get the amplitudes */
  get_ampl(taaaa, Iaaaa, evals_alpha, evals_alpha, 0);
  get_ampl(tbbbb, Ibbbb, evals_beta, evals_beta, 0);
  get_ampl(taabb, Iaabb[0], evals_alpha, evals_beta, 1);

  /* test the amplitudes by computing the energy */
  energy = compute_energy(taaaa, tbbbb, taabb, Iaaaa, Ibbbb, Iaabb);
  fprintf(outfile, "SCF energy              : %14.9lf\n", E_scf);
  fprintf(outfile, "Total Correlation energy: %14.9lf\n", energy);
  fprintf(outfile, "MP2 energy              : %14.9lf\n", energy+E_scf);

  /* allocate memory for opdms */
  Paa = block_matrix(nmo,nmo);
  Pbb = block_matrix(nmo,nmo);

  /* get the densities */
  get_dens_alpha(Paa, taaaa, taabb);
  get_dens_beta(Pbb, tbbbb, taabb);

  if (params.print_lvl > 3) {
    fprintf(outfile, "Alpha Density Matrix:\n");
    print_mat(Paa, nmo, nmo, outfile); 
    fprintf(outfile, "Beta Density Matrix:\n");
    print_mat(Pbb, nmo, nmo, outfile); 
  }

  tval = check_density_trace(Paa);
  if (params.print_lvl > 0) 
    fprintf(outfile, "Trace of alpha density matrix = %lf\n", tval);
  tval = check_density_trace(Pbb);
  if (params.print_lvl > 0) 
    fprintf(outfile, "Trace of beta density matrix  = %lf\n", tval);

  P_so_aa = block_matrix(nmo,nmo);
  P_so_bb = block_matrix(nmo,nmo);

  transform_density(Paa, P_so_aa, 0);
  transform_density(Pbb, P_so_bb, 1);

  P_so_tot = block_matrix(nmo,nmo);

  for (i=0;i<nmo; i++) {
    for (j=0;j<nmo; j++) {
      P_so_tot[i][j] = P_so_aa[i][j] + P_so_bb[i][j];
    }
  }

  /* we need to transform the total density matrix to an orthonormal
     basis so we can diagonalize it */
  /* the alpha MO's make a convenient set of orthonormal orbitals ... */
  P_mo_tot = block_matrix(nmo, nmo);
  tmat = block_matrix(nmo, nmo);
  scfvec = chkpt_rd_alpha_scf();
  
  Stri = init_array(ntri);
  Smat = block_matrix(nmo, nmo);
  errcod = iwl_rdone(PSIF_OEI, PSIF_SO_S, Stri, ntri, 0, 0, outfile);
  tri_to_sq(Stri, Smat, nmo);
  free(Stri);

  /* I ignored the difference between nmo and nso many places including 
     below.  This will cause problems if any linear dependencies. */
  mmult(P_so_tot, 0, Smat, 0, tmat, 0, nmo, nmo, nmo, 0);
  mmult(tmat, 0, scfvec, 0, P_so_tot, 0, nmo, nmo, nmo, 0);
  mmult(Smat, 0, P_so_tot, 0, tmat, 0, nmo, nmo, nmo, 0);
  mmult(scfvec, 1, tmat, 0, P_mo_tot, 0, nmo, nmo, nmo, 0);
   
  if (params.print_lvl > 3) {
    fprintf(outfile, "Total Density (alpha MO basis): \n");
    print_mat(P_mo_tot, nmo, nmo, outfile); 
  }

  /* I believe the total density matrix is correctly symmetry-blocked
   * (Pitzer order) at this stage.  However, if we diagonalize it as
   * one giant block, the orbitals will rearrange and we will lose the
   * symmetry blocking.  Therefore, we need to diagonalize a symmetry 
   * block at a time.  CDS 11/02.
   */

  /* this will be bigger than we need */
  P_eigvals = init_array(nmo);
  P_eigvecs = block_matrix(nmo,nmo);
  P_mo_block = block_matrix(nmo,nmo);
  P_so_block = block_matrix(nmo,nmo);
  nirreps = moinfo.nirreps;

  /* zero out the MO coefficient file.  block_matrix zeroes the matrix */
  chkpt_wt_scf(P_mo_block);

  for (irrep=0,mo_offset=0; irrep<nirreps; irrep++) {
    irrep_dim = moinfo.orbspi[irrep];
    if (irrep_dim==0) continue;
    scfvec_irrep = chkpt_rd_alpha_scf_irrep(irrep);
    for (i=0; i<irrep_dim; i++) {
      for (j=0; j<irrep_dim; j++) {
        /* seems already in right order */      
        /*
        i_ci = moinfo.order[i+mo_offset];
	j_ci = moinfo.order[j+mo_offset];
	*/
        i_ci = i+mo_offset;
	j_ci = j+mo_offset;
	P_mo_block[i][j] = P_mo_tot[i_ci][j_ci];
      }
    }
    
    /* Print the onepdm block */
    if (params.print_lvl > 2) {
      fprintf(outfile, "%3s block of one-particle density matrix\n",
              moinfo.labels[irrep]);
      print_mat(P_mo_block, irrep_dim, irrep_dim, outfile);
      fprintf(outfile, "\n");
    }

    /* Diagonalize the OPDM in MO basis */
    sq_rsp(irrep_dim, irrep_dim, P_mo_block, P_eigvals, 3, P_eigvecs, TOL);

    /* Print the orbitals */
    if (params.print_lvl > 2) {
      fprintf(outfile, "%3s MP2 NOs in terms of molecular orbitals\n",
              moinfo.labels[irrep]);
      eivout(P_eigvecs, P_eigvals, irrep_dim, irrep_dim, outfile);
    }

    /* Now we have to transform from MO's to SO's */
    mmult(scfvec_irrep, 0, P_eigvecs, 0, P_so_block, 0, irrep_dim, 
          irrep_dim, irrep_dim, 0);

    /* Print the orbitals */
    if (params.print_lvl > 2) {
      fprintf(outfile, "%3s MP2 NOs in terms of symmetry orbitals\n",
              moinfo.labels[irrep]);
      eivout(P_so_block, P_eigvals, irrep_dim, irrep_dim, outfile);
    }

    /* Write Natural Orbitals to PSIF_CHKPT */
    chkpt_wt_scf_irrep(P_so_block,irrep);
    free_block(scfvec_irrep);
    mo_offset += moinfo.orbspi[irrep];
  }

  fprintf(outfile, "UMP2 NOs written to checkpoint file\n\n");
  chkpt_close();

  free_block(tmat);
  free_block(Smat);
  free_block(scfvec);
  free_block(P_mo_tot);
  free_block(P_mo_block);
  free_block(P_so_tot);
  free_block(P_so_block);
  free_block(P_eigvecs);
  free(P_eigvals);
}

void get_ampl(double *ampl, double *ints, 
              double *evals1, double *evals2, int mixed)
{

  int i,j,a,b, jmax, bmax;
  int ia,jb,iajb,ib,ja,ibja,ij,ab,ijab;
  int nmo, occ, vir;
  double denom, numerator, value;

  nmo = moinfo.nmo;
  occ = moinfo.ndocc;
  vir = nmo - occ;

  for (i=0; i<occ; i++) {
    if (mixed) jmax = occ;
    else jmax = i;
    for (j=0; j<jmax; j++) {
      for (a=0; a<vir; a++) {
        ia = ioff[a+occ]+i; 
        if (mixed) bmax = vir;
        else bmax = a;
        for (b=0; b<bmax; b++) {
          denom = evals1[i] + evals2[j] - evals1[a+occ] - evals2[b+occ];
          jb = ioff[b+occ]+j;
          /* if mixed, no pq/rs symmetry, and we didn't make the ints
             sparse to take advantage of one occ, one virt index */
          if (mixed)
            iajb = ia*(nmo*(nmo+1)/2)+jb;
          else
            iajb = ioff[ia]+jb;
          numerator = ints[iajb];

          /* in the pure-spin case, we also need to subtract (ib|ja) */
          if (!mixed) {
            ib = ioff[b+occ]+i;
            ja = ioff[a+occ]+j;
            ibja = INDEX(ib,ja);
            numerator -= ints[ibja];
          }

          value = numerator/denom;
          if (mixed) { 
            ij = i*occ+j;
            ab = a*vir+b;
            ijab = ij*vir*vir+ab;  
          }
          else {
            ij = ioff[i]+j;
            ab = ioff[a]+b;
            ijab = ij*((vir*(vir+1))/2)+ab;
          }
          
          ampl[ijab] = value;

          if (fabs(ampl[ijab]) > 10.0) {
            fprintf(outfile, "Unrealistically large amplitude %lf\n", 
                    ampl[ijab]);
            fprintf(outfile, "found for %s i=%d,j=%d,a=%d,b=%d\n",
                    mixed ? "mixed" : "unmixed", i, j, a, b);
            fprintf(outfile, "numerator = %lf, denominator = %lf\n", 
                    numerator, denom);
          }

        } /* end loop over b */
      } /* end loop over a */
    } /* end loop over j */
  } /* end loop over i */ 
 

}


void get_dens_alpha(double **Paa, double *taaaa, double *taabb)
{

  int occ, vir;
  int i,j,k,a,b,c;
  int ik_aa, jk_aa, ab_aa, bc_aa, ac_aa;
  int ik_ab, jk_ab, ab_ab, bc_ab, ac_ab;
  int ij_aa, ij_ab;
  int ikab_aa, jkab_aa, ikab_ab, jkab_ab;
  int ijac_aa, ijbc_aa, ijac_ab, ijbc_ab;
  double del;
  double value, value1, value2;
  double ki_sign, kj_sign, ac_sign, bc_sign, tot_sign; /* sign on aa terms */

  occ = moinfo.ndocc;
  vir = moinfo.nmo - occ;

  /* first we get the ij elements of the density matrix */
  for (i=0; i<occ; i++) {  /* begin loop over i */
    for (j=0; j<occ; j++){ /* begin loop over j */
      del = 0.0;
      if (i==j) {
        del = 1.0;
      }
      value = 0.0;
      /* We simultaneously sum over pure-spin (bb) and over mixed-spin (ba).
         Note that for the bb amplitudes ij (occ,occ) and ab (vir,vir)
         are both lower triangles, whereas for the ab amplitudes ij and ab
         are both squares.*/ 
      for (k=0; k<occ; k++) { /* begin loop over dummy variable k */
        if (k > i) ki_sign = -1.0;
        else ki_sign = 1.0;
        if (k > j) kj_sign = -1.0;
        else kj_sign = 1.0;
        ik_aa = INDEX(i,k);
        ik_ab = i*occ+k;
        jk_aa = INDEX(j,k);
        jk_ab = j*occ+k;
        tot_sign = ki_sign * kj_sign;
        for (a=0; a<vir; a++) { 
          for (b=0; b<vir; b++) {
            ab_aa = INDEX(a,b);
            ab_ab = a*vir+b;
            /* sign on ab doesn't matter since that factor is squared */
            ikab_aa = ik_aa*((vir*(vir+1))/2)+ab_aa;
            ikab_ab = ik_ab*(vir*vir)+ab_ab;
            jkab_aa = jk_aa*((vir*(vir+1))/2)+ab_aa;
            jkab_ab = jk_ab*(vir*vir)+ab_ab;
            value += 0.5 * (tot_sign * taaaa[ikab_aa] * taaaa[jkab_aa]) + 
	             (taabb[ikab_ab] * taabb[jkab_ab]);
          } 
        } 
      } 
      Paa[i][j] = del - value;
    }
  }

  /* now we get the ab elements of the density matrix */
  for (a=0; a<vir; a++) { /* begin loop over a */
    for (b=0; b<vir; b++) { /* begin loop over b */
      value = 0.0;
      /* We simultaneously sum over pure-spin (aa) and over mixed-spin (ab).
         Note that for the aa amplitudes ij (occ,occ) and ab (vir,vir)
         are both lower triangles, whereas for the ab amplitudes ij and ab
         are both squares.*/
      for (c=0; c<vir; c++){ /* begin loop over dummy variable c */
        ac_aa = INDEX(a,c);
        ac_ab = a*vir+c;
        bc_aa = INDEX(b,c);
        bc_ab = b*vir+c;
        if (c > a) ac_sign = -1.0;
        else ac_sign = 1.0;
        if (c > b) bc_sign = -1.0;
        else bc_sign = 1.0;
        tot_sign = ac_sign * bc_sign;
        for (i=0; i<occ; i++){ /* begin loop over i */
          for (j=0; j<occ; j++){ /* begin loop over j */
            ij_aa = INDEX(i,j);
            ij_ab = i*occ+j;
            ijac_aa = ij_aa*((vir*(vir+1))/2)+ac_aa;
            ijac_ab = ij_ab*(vir*vir)+ac_ab;
            ijbc_aa = ij_aa*((vir*(vir+1))/2)+bc_aa;
            ijbc_ab = ij_ab*(vir*vir)+bc_ab;
            value += 0.5 * (tot_sign * taaaa[ijac_aa] * taaaa[ijbc_aa]) + 
                     (taabb[ijac_ab] * taabb[ijbc_ab]);
          } /* end loop over j */
        } /* end loop over i */
      } /* end loop over c */
      Paa[a+occ][b+occ] = value;
    }
  }

}


void get_dens_beta(double **Pbb, double *tbbbb, double *taabb)
{

  int occ, vir;
  int i,j,k,a,b,c;
  int ki_ab, kj_ab, ba_ab, cb_ab, ca_ab;
  int ik_bb, jk_bb;
  int ab_bb, ac_bb, bc_bb, ij_bb, ji_ab;
  int ikab_bb, jkab_bb, kiba_ab, kjba_ab;
  int ijac_bb, ijbc_bb, jica_ab, jicb_ab;
  double del;
  double value;
  double ki_sign, kj_sign, ac_sign, bc_sign, tot_sign; /* sign on bb terms */

  occ = moinfo.ndocc;
  vir = moinfo.nmo - occ;

  /* first we get the ij elements of the density matrix */
  for (i=0; i<occ; i++){ /* begin loop over i */
    for (j=0; j<occ; j++){ /* begin loop over j */
      del = 0.0;
      if (i==j) {
        del = 1.0;
      }
      value = 0.0;
      /* We simultaneously sum over pure-spin (bb) and over mixed-spin (ba).
         Note that for the bb amplitudes ij (occ,occ) and ab (vir,vir)
         are both lower triangles, whereas for the ab amplitudes ij and ab
         are both squares.*/
      for (k=0; k<occ; k++){ /* begin loop over dummy variable k */
        if (k > i) ki_sign = -1.0;
        else ki_sign = 1.0;
        if (k > j) kj_sign = -1.0;
        else kj_sign = 1.0;
        ik_bb = INDEX(i,k);
        ki_ab = k*occ+i;
        jk_bb = INDEX(j,k);
        kj_ab = k*occ+j;
        tot_sign = ki_sign * kj_sign;
        for (a=0; a<vir; a++){ /* begin loop over a */
          for (b=0; b<vir; b++){ /* begin loop over b */
            ab_bb = INDEX(a,b);
            ba_ab = b*vir+a;
            /* sign on ab doesn't matter since that factor is squared */
            ikab_bb = ik_bb*((vir*(vir+1))/2)+ab_bb;
            kiba_ab = ki_ab*(vir*vir)+ba_ab;
            jkab_bb = jk_bb*((vir*(vir+1))/2)+ab_bb;
            kjba_ab = kj_ab*(vir*vir)+ba_ab;
            value += 0.5 * (tot_sign * tbbbb[ikab_bb]*tbbbb[jkab_bb]) + 
                     (taabb[kiba_ab]*taabb[kjba_ab]);
          } /* end loop over b */
        } /* end loop over a */
      } /* end loop over k */
      Pbb[i][j] = del - value;
    }
  }

  /* now we get the ab elements of the density matrix */
  for (a=0; a<vir; a++){ /* begin loop over a */
    for (b=0; b<vir; b++){ /* begin loop over b */
      value = 0.0;
      /* We simultaneously sum over pure-spin (bb) and over mixed-spin (ab).
         Note that for the bb amplitudes ij (occ,occ) and ab (vir,vir)
         are both lower triangles, whereas for the ab amplitudes ij and ab
         are both squares.*/
      for (c=0; c<vir; c++){ /* begin loop over dummy variable c */
        if (c > a) ac_sign = -1.0;
        else ac_sign = 1.0;
        if (c > b) bc_sign = -1.0;
        else bc_sign = 1.0;
        tot_sign = ac_sign * bc_sign;
        ac_bb = INDEX(a,c);
        ca_ab = c*vir+a;
        bc_bb = INDEX(b,c);
        cb_ab = c*vir+b;
        for (i=0; i<occ; i++){ /* begin loop over i */
          for (j=0; j<occ; j++){ /* begin loop over j */
            ij_bb = INDEX(i,j);
            ji_ab = j*occ+i;
            ijac_bb = ij_bb*((vir*(vir+1))/2)+ac_bb;
            ijbc_bb = ij_bb*((vir*(vir+1))/2)+bc_bb;
            jica_ab = ji_ab*(vir*vir)+ca_ab;
		    ijbc_bb = ij_bb*((vir*(vir+1))/2)+bc_bb;
		    jicb_ab = ji_ab*(vir*vir)+cb_ab;
		    value += 0.5 * (tot_sign * tbbbb[ijac_bb]*tbbbb[ijbc_bb]) + 
			     (taabb[jica_ab]*taabb[jicb_ab]);
		  } /* end loop over j */
		} /* end loop over i */
	      } /* end loop over c */
	      Pbb[a+occ][b+occ] = value;
    }
  }

}


double compute_energy(double *taaaa, double *tbbbb, double *taabb, 
                      double *Iaaaa, double *Ibbbb, double **Iaabb)
{
  int nmo, occ, vir, ntri;
  int i, j, a, b, ij, ab, ijab, ia, jb, iajb, ib, ja, ibja;
  double Eaa, Ebb, Eab;
  
  nmo = moinfo.nmo;
  occ = moinfo.ndocc;
  vir = nmo - occ;
  ntri = (nmo * (nmo+1))/2;
 
  /* alpha-alpha component */
  Eaa = 0.0;
  for (i=0, ij=0; i<occ; i++) {
    for (j=0; j<=i; j++, ij++) {
      for (a=0, ab=0; a<vir; a++) {
        ia = ioff[a+occ]+i;
        ja = ioff[a+occ]+j;
        for (b=0; b<=a; b++, ab++) {
          ijab = ij*((vir*(vir+1))/2)+ab;
          ib = ioff[b+occ]+i;
          jb = ioff[b+occ]+j;
          iajb = ioff[ia]+jb; 
          ibja = INDEX(ib,ja);
          Eaa += taaaa[ijab]*(Iaaaa[iajb] - Iaaaa[ibja]);
        }
      }
    }
  }     

  fprintf(outfile, "\nContributions to the correlation energy:\n");
  fprintf(outfile, "     Eaa = %14.6lf\n", Eaa);

  /* beta-beta component */
  Ebb = 0.0;
  for (i=0, ij=0; i<occ; i++) {
    for (j=0; j<=i; j++, ij++) {
      for (a=0, ab=0; a<vir; a++) {
        ia = ioff[a+occ]+i;
        ja = ioff[a+occ]+j;
        for (b=0; b<=a; b++, ab++) {
          ijab = ij*((vir*(vir+1))/2)+ab;
          ib = ioff[b+occ]+i;
          jb = ioff[b+occ]+j;
          iajb = ioff[ia]+jb; 
          ibja = INDEX(ib,ja);
          Ebb += tbbbb[ijab]*(Ibbbb[iajb] - Ibbbb[ibja]);
        }
      }
    }
  }  
  fprintf(outfile, "     Ebb = %14.6lf\n", Ebb);

  /* mixed (alpha-beta) component */
  Eab = 0.0;
  for (i=0; i<occ; i++) {
    for (j=0; j<occ; j++) {
      ij = i*occ+j;
      for (a=0; a<vir; a++) {
        ia = ioff[a+occ]+i;
        for (b=0; b<vir; b++) {
          jb = ioff[b+occ]+j;
          ab = a*vir+b;
          ijab = ij*vir*vir+ab;
          Eab += taabb[ijab]*Iaabb[ia][jb];
        }
      }
    }
  }
  fprintf(outfile, "     Eab = %14.6lf\n", Eab);

  return(Eaa+Ebb+Eab);
}


double check_density_trace(double **P) 
{
  int nmo, p;
  double tval;

  nmo = moinfo.nmo;

  tval = 0.0;
  for (p=0; p<nmo; p++) 
    tval += P[p][p];
  
  return(tval); 
}

}} // end namespace psi::mvo
