/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<memory.h>
#include<pthread.h>
#include<libipv1/ip_lib.h>
#include<libciomr/libciomr.h>
#include<libpsio/psio.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

#include"quartet_data.h"      /* From Default_Ints */
#include"norm_quartet.h"
#include"hash.h"
#include"transmat.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
#endif
#include"schwartz.h"
#include"shell_block_matrix.h"

using namespace std;

extern pthread_mutex_t fock_mutex;

namespace psi { 
  namespace cints {
extern void *hf_fock_thread(void *);

/*--- To be accessed by all HF Fock threads ---*/
double ****Gskel, ****Gskel_o;       /* Shell-blocked skeleton G matrices */

void hf_fock()
{
  pthread_attr_t thread_attr;
  pthread_t *thread_id;
  
  int nstri;
  int count ;
  int dum;
  int g, i, j, k, l, m, ii, jj, kk, ll;
  int si, sj, ni, nj, li, lj, si_g, sj_g;

  double temp;
  double **tmpmat1;
  double ****Gfull, ****Gfull_o;  /* Shell-blocked G matrices in AO basis*/
  double ****Gsym, ****Gsym_o;    /* Shell-blocked symmetrized (Gskel + Gskel(transp.)) G matrices */
  double ***ao_type_transmat;
  double *Gtri, *Gtri_o;          /* Total G matrices in lower*/
                                  /* triagonal form*/
  /*---------------
    Initialization
   ---------------*/
#ifdef USE_TAYLOR_FM
/*  init_Taylor_Fm_Eval(BasisSet.max_am*4-4,UserOptions.cutoff);*/
  init_Taylor_Fm_Eval(BasisSet.max_am*4-4,1.0E-20);
#else
  init_fjt(BasisSet.max_am*4);
#endif
  init_libint_base();

  /*------------------------------------------
    Compute integrals for Schwartz inequality
   ------------------------------------------*/
  schwartz_eri();

  /*------------------------------------
    Allocate shell-blocked skeleton G's
   ------------------------------------*/
  Gskel = init_shell_block_matrix();
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf)
    Gskel_o = init_shell_block_matrix();

  thread_id = (pthread_t *) malloc(UserOptions.num_threads*sizeof(pthread_t));
  pthread_attr_init(&thread_attr);
  pthread_attr_setscope(&thread_attr,
			PTHREAD_SCOPE_SYSTEM);
  pthread_mutex_init(&fock_mutex,NULL);
  for(long int i=0;i<UserOptions.num_threads-1;i++)
    pthread_create(&(thread_id[i]),&thread_attr,
		   hf_fock_thread,(void *)i);
  hf_fock_thread( (void *) (UserOptions.num_threads - 1) );
  for(i=0;i<UserOptions.num_threads-1;i++)
    pthread_join(thread_id[i], NULL);
  free(thread_id);
  pthread_mutex_destroy(&fock_mutex);
  
  /*-------------------------------
    Gskel = Gskel + Gskel(transp.)
   -------------------------------*/
  Gsym = init_shell_block_matrix();
  GplusGt(Gskel,Gsym);
  free_shell_block_matrix(Gskel);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    Gsym_o = init_shell_block_matrix();
    GplusGt(Gskel_o,Gsym_o);
    free_shell_block_matrix(Gskel_o);
  }
  
  /*-----------------
    Symmetrize Gskel
   -----------------*/
  if (Symmetry.nirreps > 1) {
    ao_type_transmat = build_transmat(Symmetry.sym_oper, Symmetry.nirreps, BasisSet.max_am);
    Gfull = init_shell_block_matrix();
    for(g=0;g<Symmetry.nirreps;g++) {
      for(si=0;si<BasisSet.num_shells;si++) {
	ni = ioff[BasisSet.shells[si].am];
	li = BasisSet.shells[si].am-1;
	si_g = BasisSet.shells[si].trans_vec[g] - 1;
	for(sj=0;sj<BasisSet.num_shells;sj++) {
	  sj_g = BasisSet.shells[sj].trans_vec[g] - 1;
	  nj = ioff[BasisSet.shells[sj].am];
	  lj = BasisSet.shells[sj].am-1;

	  for(i=0;i<ni;i++)
	    for(j=0;j<nj;j++)
	      Gfull[si_g][sj_g][i][j] += ao_type_transmat[li][g][i]*
					   ao_type_transmat[lj][g][j]*
					   Gsym[si][sj][i][j];
	}
      }
    }
  }
  else
    Gfull = Gsym;

  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    if (Symmetry.nirreps > 1) {
    Gfull_o = init_shell_block_matrix();
    for(g=0;g<Symmetry.nirreps;g++) {
      for(si=0;si<BasisSet.num_shells;si++) {
	ni = ioff[BasisSet.shells[si].am];
	li = BasisSet.shells[si].am-1;
	si_g = BasisSet.shells[si].trans_vec[g] - 1;
	for(sj=0;sj<BasisSet.num_shells;sj++) {
	  sj_g = BasisSet.shells[sj].trans_vec[g] - 1;
	  nj = ioff[BasisSet.shells[sj].am];
	  lj = BasisSet.shells[sj].am-1;

	  for(i=0;i<ni;i++)
	    for(j=0;j<nj;j++)
	      Gfull_o[si_g][sj_g][i][j] += ao_type_transmat[li][g][i]*
					   ao_type_transmat[lj][g][j]*
					   Gsym_o[si][sj][i][j];
	}
      }
    }
    }
    else
      Gfull_o = Gsym_o;
  }
  /*--------------------
    Print out G for now
   --------------------*/
  G = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  shell_block_to_block(Gfull,G);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    Go =  block_matrix(BasisSet.num_ao,BasisSet.num_ao);
    shell_block_to_block(Gfull_o,Go);
  }
   /*----------------------
    Transform to SO basis
   ----------------------*/
  if (Symmetry.nirreps > 1 || BasisSet.puream) {
    tmpmat1 = block_matrix(Symmetry.num_so,BasisSet.num_ao);
    mmult(Symmetry.usotao,0,G,0,tmpmat1,0,Symmetry.num_so,BasisSet.num_ao,BasisSet.num_ao,0);
    mmult(tmpmat1,0,Symmetry.usotao,1,G,0,Symmetry.num_so,BasisSet.num_ao,Symmetry.num_so,0);
    if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
      mmult(Symmetry.usotao,0,Go,0,tmpmat1,0,Symmetry.num_so,BasisSet.num_ao,BasisSet.num_ao,0);
      mmult(tmpmat1,0,Symmetry.usotao,1,Go,0,Symmetry.num_so,BasisSet.num_ao,Symmetry.num_so,0);
    }
    free_block(tmpmat1);
  }
  /*
  fprintf(outfile,"  Closed-shell Fock matrix in SO basis:\n");
  print_mat(G,Symmetry.num_so,Symmetry.num_so,outfile);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
      fprintf(outfile,"  Open-shell Fock matrix in SO basis:\n");
      print_mat(Go,Symmetry.num_so,Symmetry.num_so,outfile);
      } */
  
  /*-------------------------
    Write G-matrices to disk
   -------------------------*/
  nstri = ioff[Symmetry.num_so];
  Gtri = init_array(nstri);
  sq_to_tri(G,Gtri,Symmetry.num_so);
  free_block(G);
  psio_open(IOUnits.itapDSCF, PSIO_OPEN_OLD);
  switch (UserOptions.reftype) {
  case rohf:
      Gtri_o = init_array(nstri);
      sq_to_tri(Go,Gtri_o,Symmetry.num_so);
      free_block(Go);
      psio_write_entry(IOUnits.itapDSCF, "Open-shell JX G-matrix", (char *) Gtri_o, sizeof(double)*nstri);
      free(Gtri_o);
  case rhf:
      psio_write_entry(IOUnits.itapDSCF, "Total JX G-matrix", (char *) Gtri, sizeof(double)*nstri);
      free(Gtri);
      break;

  case uhf:
      Gtri_o = init_array(nstri);
      sq_to_tri(Go,Gtri_o,Symmetry.num_so);
      free_block(Go);
      /*--- Form alpha and beta Fock matrices first and then write them out ---*/
      for(i=0;i<nstri;i++) {
	temp = Gtri[i] + Gtri_o[i];
	Gtri[i] = Gtri[i] - Gtri_o[i];
	Gtri_o[i] = temp;
      }
      
      psio_write_entry(IOUnits.itapDSCF, "Alpha JX G-matrix", (char *) Gtri, sizeof(double)*nstri);
      psio_write_entry(IOUnits.itapDSCF, "Beta JX G-matrix", (char *) Gtri_o, sizeof(double)*nstri);
      free(Gtri);
      free(Gtri_o);
      break;
  }
  psio_close(IOUnits.itapDSCF, 1);

  /*---------
    Clean-up
   ---------*/
  free_shell_block_matrix(Gsym);
  if (Symmetry.nirreps > 1) {
    free_shell_block_matrix(Gfull);
  }
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    free_shell_block_matrix(Gsym_o);
    if (Symmetry.nirreps > 1)
      free_shell_block_matrix(Gfull_o);
  }
#ifdef USE_TAYLOR_FM
  free_Taylor_Fm_Eval();
#else
  free_fjt();
#endif

  return;
}

}}
