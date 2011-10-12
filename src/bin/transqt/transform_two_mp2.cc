/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"
#include "yoshimine.h"

namespace psi {
extern FILE* outfile;
namespace transqt {

#define MAXIOFF3 255
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void exit_io();
static void make_arrays(double ****Cdocc, double ****Cvirt,
                 int **active_docc, int **active_virt,
                 int **focact, int **locact, int **fvract, int **lvract,
                 int **occ, int **vir, int **ioff3, int nvirt);

void transform_two_mp2(void)
{
  int p,q,r,s,pq,rs;
  int psym, qsym, rsym, ssym,pqsym;
  int pfirst, qfirst, rfirst, sfirst;
  int plast, qlast, rlast, slast;
  int k,l,kl;
  int ksym, lsym, klsym;
  int kfirst, lfirst, ifirst, jfirst;
  int klast, llast, ilast, jlast;
  int R,S,K,L,P,Q;
  int h;
  int ntri,nmo,nso;
  int ndocc, nvirt;
  int nirreps;
  int *clsdpi,*orbspi, *sopi;
  int *first, *last, *first_so, *last_so;
  int *active_virt, *active_docc;
  int *focact, *locact, *fvract, *lvract;
  int *occ, *vir;
  int *ioff3;
  int print_integrals;
  long int maxcor, maxcord;
  double **A, **B;
  double ***Cdocc, ***Cvirt;
  double *P_block;
  double *J_block;
  struct iwlbuf PBuff;
  struct iwlbuf JBuff;
  struct iwlbuf MBuff;
  struct yoshimine YBuffP, YBuffJ;
  int max_buckets, print_lvl, first_tmp_file;
  int presort_file, keep_presort, jfile, keep_half_tf, mfile;
  double tolerance;


  if(moinfo.iopen) {
      fprintf(outfile,"\n\tThis system seems to be open-shell.\n");
      fprintf(outfile,"MP2 restricted transformations are available only ");
      fprintf(outfile,"for closed-shell molecules.\n");
      exit_io();
      exit(PSI_RETURN_FAILURE);
    }

  if (params.backtr) {
      fprintf(outfile,"Warning! You have set backtr = true.\n");
      fprintf(outfile,"This is not yet implemented for MP2. \n");
      fprintf(outfile,"Ignoring backtr option.\n");
    }

  ntri = moinfo.noeints;
  nso = moinfo.nso;
  nmo = moinfo.nmo;
  nirreps = moinfo.nirreps;
  sopi = moinfo.sopi;
  orbspi = moinfo.orbspi;
  clsdpi = moinfo.clsdpi;
  first_so = moinfo.first_so;
  last_so = moinfo.last_so;
  first = moinfo.first;
  last = moinfo.last;
  max_buckets = params.max_buckets;
  print_lvl = params.print_lvl;
  first_tmp_file = params.first_tmp_file;
  presort_file = params.presort_file;
  keep_presort = params.keep_presort;
  jfile = params.jfile;
  keep_half_tf = params.keep_half_tf;
  mfile = params.mfile;
  tolerance = params.tolerance;
  print_integrals = params.print_te_ints;
  maxcor = params.maxcor;
  maxcord = params.maxcord;


  /** Pre-sort the two-electron Integrals **/
  if (params.print_lvl) {
    fprintf(outfile, "\n\tPre-sorting two-electron ints...\n\n");
    fflush(outfile);
  }

  yosh_init(&YBuffP, ntri, ntri, maxcor, maxcord, max_buckets,
            first_tmp_file, tolerance, outfile);

  if (print_lvl > 1) {
    fprintf(outfile, "\tPresort");
    yosh_print(&YBuffP, outfile);
    fprintf(outfile, "\n");
    fflush(outfile);
  }

  yosh_init_buckets(&YBuffP);

  yosh_rdtwo(&YBuffP, params.src_tei_file, params.delete_src_tei,
             sopi, nirreps, ioff, 0,
             params.fzc && moinfo.nfzc, moinfo.fzc_density,
             moinfo.fzc_operator, 1, (print_lvl > 5), outfile);

  yosh_close_buckets(&YBuffP, 0);

  yosh_sort(&YBuffP, presort_file, 0, ioff, NULL, nso, ntri,
            0, 1, 0, 0, 1, (print_lvl > 5), outfile);

  yosh_done(&YBuffP);  /* Pre-transform complete */

  /* no the frozen core here */

  ndocc = 0;
  nvirt = 0;
  for(h=0; h < nirreps; h++) {
      ndocc += clsdpi[h];
      nvirt += orbspi[h] - clsdpi[h];
    }

  fprintf(outfile,
          "\tEntering MP2 two-electron integral transformation...\n\n");

  iwl_buf_init(&PBuff, presort_file, tolerance, 1, 1);
  yosh_init(&YBuffJ, ndocc*nvirt, ntri, maxcor, maxcord, max_buckets,
            first_tmp_file, tolerance, outfile);

  yosh_init_buckets(&YBuffJ);

  fprintf(outfile, "\tHalf-transform");
  yosh_print(&YBuffJ, outfile);
  fprintf(outfile, "\n");
  fflush(outfile);

  A = block_matrix(nso,nso);
  B = block_matrix(nso,nso);
  P_block = init_array(ntri);

  make_arrays(&Cdocc, &Cvirt, &active_docc, &active_virt,
              &focact, &locact, &fvract, &lvract, &occ, &vir, &ioff3,nvirt);

  for(psym=0; psym < nirreps; psym++) {
      pfirst = first_so[psym];
      plast = last_so[psym];
      for(p=pfirst; p <= plast; p++) {
          for(qsym=0; qsym <= psym; qsym++) {
              qfirst = first_so[qsym];
              qlast = last_so[qsym];
              pqsym = psym^qsym;
              for(q=qfirst; (q<=qlast) && (q <= p); q++) {
                  pq = ioff[p] + q;

                  zero_arr(P_block,ntri);
                  iwl_buf_rd(&PBuff, pq, P_block, ioff, ioff, 0, 0, outfile);

                 for(rsym=0; rsym < nirreps; rsym++) {
                      rfirst = first_so[rsym];
                      rlast = last_so[rsym];
                      kfirst = focact[rsym];
                      klast = locact[rsym];
                      ssym = pqsym^rsym;
                      sfirst = first_so[ssym];
                      slast = last_so[ssym];
                      lfirst = fvract[ssym];
                      llast = lvract[ssym];
                      for(r=rfirst,R=0; r <= rlast; r++,R++) {
                          for(s=sfirst,S=0; s <= slast; s++,S++) {
                              rs = INDEX(r,s);
                              A[R][S] = P_block[rs];
                            }
                        }
                      mmult(A,0,Cvirt[ssym],0,B,0,sopi[rsym],sopi[ssym],
                            active_virt[ssym],0);
                      mmult(Cdocc[rsym],1,B,0,A,0,active_docc[rsym],
                            sopi[rsym],active_virt[ssym],0);
                      yosh_wrt_arr_mp2(&YBuffJ, p, q, pq, pqsym, A,
                                       rsym, focact, locact, fvract, lvract,
                                       1, ndocc, nvirt, occ, vir, ioff3,
                                       (print_lvl >4), outfile);
                   }
                }
            }
        }
    }

  fprintf(outfile, "\tSorting half-transformed integrals...\n");

  iwl_buf_close(&PBuff, keep_presort);
  yosh_flush(&YBuffJ);
  yosh_close_buckets(&YBuffJ,0);
  yosh_sort(&YBuffJ, jfile, 0, ioff3, ioff, nso, ntri, 0, 1, 1, nvirt,
            0, (print_lvl > 5), outfile);
  yosh_done(&YBuffJ);

  fprintf(outfile, "\tFinished half-transform...\n");
  fprintf(outfile, "\tWorking on second half...\n");
  fflush(outfile);

  iwl_buf_init(&JBuff, jfile, tolerance, 1, 1);
  iwl_buf_init(&MBuff, mfile, tolerance, 0, 0);

  J_block = P_block;

  for(ksym=0; ksym < nirreps; ksym++) {
      kfirst = (focact[ksym] < 0) ? -1 : occ[focact[ksym]];
      klast = (locact[ksym] < 0) ? -2 : occ[locact[ksym]];
      for(k=kfirst; k <= klast; k++) {
          for(lsym=0; lsym < nirreps; lsym++) {
              lfirst = (fvract[lsym] < 0) ? -1 : vir[fvract[lsym]];
              llast = (lvract[lsym] < 0) ? -2 : vir[lvract[lsym]];
              klsym = ksym^lsym;
              for(l=lfirst; l <= llast; l++) {
                  kl = ioff3[k] + l;

                  zero_arr(J_block, ntri);
                  iwl_buf_rd(&JBuff, kl, J_block, ioff3, ioff, 1, 0, outfile);

                  for(psym=0; psym < nirreps; psym++) {
                      pfirst = first_so[psym];
                      plast = last_so[psym];
                      ifirst = focact[psym];
                      ilast = locact[psym];
                      qsym = klsym^psym;
                      qfirst = first_so[qsym];
                      qlast = last_so[qsym];
                      jfirst = fvract[qsym];
                      jlast = lvract[qsym];
                      for(p=pfirst,P=0; p <= plast; p++,P++) {
                          for(q=qfirst,Q=0; q <= qlast; q++,Q++) {
                              pq = INDEX(p,q);
                              A[P][Q] = J_block[pq];
                            }
                        }
                      mmult(A,0,Cvirt[qsym],0,B,0,sopi[psym],sopi[qsym],
                            active_virt[qsym],0);
                      mmult(Cdocc[psym],1,B,0,A,0,active_docc[psym],
                            sopi[psym],active_virt[qsym],0);
                      iwl_buf_wrt_mp2(&MBuff, k, l, kl, klsym, A, psym,
                                      focact, locact, fvract, lvract,
                                      occ,vir,ioff3,print_integrals,outfile);
                    }
                }
            }
        }
    }

  free(J_block);
  free_block(A);
  free_block(B);
  iwl_buf_close(&JBuff, keep_half_tf);
  /* Need to flush last buffer, else it's not written to disk */
  iwl_buf_flush(&MBuff, 1);
  iwl_buf_close(&MBuff, 1);

  fprintf(outfile, "\n\tTransformation finished.\n");

  fprintf(outfile, "\tTwo-electron integrals written to file%d.\n",mfile);
  fflush(outfile);

  return;
}

void make_arrays(double ****Cdocc, double ****Cvirt,
                 int **active_docc, int **active_virt,
                 int **focact, int **locact, int **fvract, int **lvract,
                 int **occ, int **vir, int **ioff3, int nvirt)
{
  int h, i, offset, count;
  int first_offset, last_offset;
  int p,q,row,col;

  /* active_docc[] and active_uocc[] are the number of active doubly-occupied
     and unoccupied orbitals, respectively, in each irrep. */
  *active_docc = init_int_array(moinfo.nirreps);
  *active_virt = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) {
      (*active_docc)[h] = moinfo.clsdpi[h] - moinfo.frdocc[h];
      (*active_virt)[h] = moinfo.virtpi[h] - moinfo.fruocc[h];
    }

  /* focact[] and locact[] supply the first and last orbital indices (Pitzer
     ordering) for the occupied orbitals in each irrep. */
  *focact = init_int_array(moinfo.nirreps);
  *locact = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) {
      (*focact)[h] = -1;
      (*locact)[h] = -2;
    }
  first_offset = moinfo.frdocc[0];
  last_offset = moinfo.clsdpi[0] - 1;
  (*focact)[0] = first_offset;
  (*locact)[0] = last_offset;
  for(h=1; h < moinfo.nirreps; h++) {
      first_offset += moinfo.orbspi[h-1]+moinfo.frdocc[h]-moinfo.frdocc[h-1];
      last_offset += moinfo.virtpi[h-1]+moinfo.orbspi[h]-moinfo.virtpi[h];
      if((*active_docc)[h]) {
          (*focact)[h] = first_offset;
          (*locact)[h] = last_offset;
        }
    }

  /* frvact[] and lrvact[] supply the first and last orbital indices (Pitzer
     ordering) for the virtual orbitals in each irrep. */
  *fvract = init_int_array(moinfo.nirreps);
  *lvract = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) {
      (*fvract)[h] = -1;
      (*lvract)[h] = -2;
    }
  first_offset = moinfo.clsdpi[0];
  last_offset = moinfo.orbspi[0] - moinfo.fruocc[0] - 1;
  (*fvract)[0] = first_offset;
  (*lvract)[0] = last_offset;
  for(h=1; h < moinfo.nirreps; h++) {
      first_offset += moinfo.virtpi[h-1] + moinfo.clsdpi[h];
      last_offset += moinfo.fruocc[h-1] + moinfo.orbspi[h]-moinfo.fruocc[h];
      if((*active_virt)[h]) {
          (*fvract)[h] = first_offset;
          (*lvract)[h] = last_offset;
        }
    }

  /* Construct occupied and virtual Pitzer -> QTS ordering arrays for
   occupied (occ[]) and virtual (vir[]) orbitals */
  *occ = init_int_array(moinfo.nmo);
  *vir = init_int_array(moinfo.nmo);
  for(i=0; i< moinfo.nmo; i++) {
      (*occ)[i] = -1;
      (*vir)[i] = -1;
    }

  offset = 0;
  count=0;
  for(h=0; h < moinfo.nirreps; h++) {
      if(h)
          offset += moinfo.orbspi[h-1];
      for(i=offset; i < (offset+moinfo.clsdpi[h]); i++) {
          (*occ)[i] = count++;
        }
    }

  count=0;
  offset=0;
  for(h=0; h < moinfo.nirreps; h++) {
      if(h)
          offset += moinfo.orbspi[h-1];
      for(i=offset+moinfo.clsdpi[h]; i < (offset+moinfo.orbspi[h]); i++) {
          (*vir)[i] = count++;
        }
    }

  /* Generate ioff3 array.  This array gives the row offset for an
     ndocc x nvirt matrix */
  *ioff3 = init_int_array(MAXIOFF3);
  for(i=0; i < MAXIOFF3; i++) {
      (*ioff3)[i] = i*nvirt;
    }

  /* Organize MOs for docc and virt sets only */
  if(params.print_mos) {
      fprintf(outfile,"\n\tSCF Eigenvectors (Occupied and Virtual Sets):\n");
      fprintf(outfile,  "\t---------------------------------------------\n");
    }
  *Cdocc = (double ***) malloc(moinfo.nirreps * sizeof(double **));
  *Cvirt = (double ***) malloc(moinfo.nirreps * sizeof(double **));
  for(h=0; h < moinfo.nirreps; h++) {
      if((*active_docc)[h]) {
          (*Cdocc)[h] = block_matrix(moinfo.sopi[h],(*active_docc)[h]);
          row = -1;
          for(p=moinfo.first_so[h]; p <= moinfo.last_so[h]; p++) {
              row++; col=-1;
              for(q=(*focact)[h]; q <= (*locact)[h]; q++) {
                  col++;
                  (*Cdocc)[h][row][col] = moinfo.scf_vector[p][q];
                }
            }
          if(params.print_mos) {
              fprintf(outfile,"\n\tDoubly Occupied Orbitals for Irrep %s\n",
                      moinfo.labels[h]);
              print_mat((*Cdocc)[h],moinfo.sopi[h],(*active_docc)[h],
                        outfile);
            }
        }
      if((*active_virt)[h]) {
          (*Cvirt)[h] = block_matrix(moinfo.sopi[h],(*active_virt)[h]);
          row = -1;
          for(p=moinfo.first_so[h]; p <= moinfo.last_so[h]; p++) {
              row++; col=-1;
              for(q=(*fvract)[h]; q <= (*lvract)[h]; q++) {
                  col++;
                  (*Cvirt)[h][row][col] = moinfo.scf_vector[p][q];
                }
            }
          if(params.print_mos) {
              fprintf(outfile,"\n\tVirtual Orbitals for Irrep %s\n",
                      moinfo.labels[h]);
              print_mat((*Cvirt)[h],moinfo.sopi[h],(*active_virt)[h],
                        outfile);
            }
        }
    }

}

}} // end namespace psi::transqt

