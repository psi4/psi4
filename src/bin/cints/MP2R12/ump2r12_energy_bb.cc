#include<cmath>
#include<cstdio>
#include<cstring>
#include<memory.h>
#include<cstdlib>
#include<libciomr/libciomr.h>
#include<libqt/qt.h>
#include<libiwl/iwl.h>
#include<libint/libint.h>
#include<libr12/libr12.h>
#include<pthread.h>
#include <stdexcept>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"r12_quartet_data.h"
#include"norm_quartet.h"
#include"int_fjt.h"
#include"quartet_permutations.h"
#include"rmp2r12_energy.h"

namespace psi{ namespace cints{

extern void make_transqt_arrays_uhf(int **first, int **last, int **fstocc_alpha, int **fstocc_beta, int **lstocc_alpha, 
                         int **lstocc_beta, int **occ_alpha, int **occ_beta, int **act2fullQTS_alpha,
                         int **act2fullQTS_beta, int **ioff3);

namespace ump2r12_bb{

/*-------------------------------------------------------
  Algorithm

  ***Split into threads. Each thread do the following:
  
  Loop over I batches (batch size num_i_per_batch) of active DOCC

    Loop over all symmetry-unique shells UR, US<=UR
      ***if this UR, US is to be handled by this thread - proceed, else skip to next one
      Find all symmetry-distinct shell doublets resulting from (UR US|
      Loop over all resulting shell doublets R, S

        Loop over all symmetry-unique shells UP, UQ<=UP
	  Find all symmetry-distinct shell quartets resulting from (R S|UP UQ)
	  Loop over the resulting set of P, Q doublets
            Evaluate (RS|PQ), (RS|r12|PQ), and (RS|[r12,T1]|PQ)
	    Loop over p in P, q in Q, r in R, s in S, i in I
	      (rs|iq) += Cpi * (rs|pq)
	      (rs|ip) += Cqi * (rs|pq)
	      same for (rs|r12|pq) and (rs|[r12,T1]|pq)
	    End p, q, r, s, i loop
          End P,Q loop
        End UP, UQ loop

        Loop over r in R, s in S
          Loop over q < nao, x < num_mo, i in I
	    (rs|ix) += Cxs * (rs|is)
	    same for (rs|r12|is) and (rs|[r12,T1]|is)
	  End q, x, i
	End r, s loop
	***Lock (js| and (jr| (either individual or shell blocks depending on LOCK_RS_SHELL in rmp2r12_energy.h)
	Loop over r in R, s in S
	  Loop over i in I, x < num_mo, j <= i
	    (js|ix) += Cjr * (rs|ix)
	    (jr|ix) += Cjs * (rs|ix)
	    same for (rs|r12|ix), but
	    (js|[r12,T1]|ix) += Cjr * (rs|[r12,T1]|ix)
	    (jr|[r12,T1]|ix) -= Cjs * (rs|[r12,T1]|ix)                   <---- Note the minus sign here!!!
	  End i, x, j loop
        End r, s loop
	***Unlock (js| and (jr|

      End R, S loop
    End UR, US loop

    ***Barrier: threads wait until everyone is done
    ***Do the following in one thread only
    Loop over i in I, j <= i
      Loop over r < nao, x < num_mo, y < num_mo
        (jy|ix) += Cys * (js|ix)
	same for (js|r12|ix) and (js|[r12,T1]|ix)
      End r, x, y loop
    End i, j loop

  End I loop

  ***Merge all threads
  
 -------------------------------------------------------*/

/*-------------------------------
  Explicit function declarations
 -------------------------------*/
extern void *ump2r12_energy_thread_bb(void *);

/*--------------------------------
  Variables common to all threads
 --------------------------------*/
pthread_mutex_t rmp2r12_energy_mutex;
pthread_mutex_t *rmp2r12_sindex_mutex;
pthread_cond_t rmp2r12_energy_cond;
RMP2R12_Status_t RMP2R12_Status;
double *jsix_buf[NUM_TE_TYPES];             /* buffer for (js|ia) integrals, where j runs over all d.-o. MOs,
 					 s runs over all AOs, i - over I-batch, x - over all MOs */
double *jyix_buf[NUM_TE_TYPES];             /* buffer contains all MP2-R12/A-type integrals */
double *jsix_buf2;             /* buffer for (js|ia) integrals, where j runs over all d.-o. MOs,
 					 s runs over all AOs, i - over I-batch, x - over all MOs */
double *jyix_buf2;             /* buffer contains all MP2-R12/A-type integrals */

int *first, *last;                          /* first and last absolute (Pitzer) orbital indices in symblk */
int *fstocc_alpha, *lstocc_alpha;           /* first and last occupied indices in Pitzer ordering for each symblk */
int *fstocc_beta, *lstocc_beta;             /* first and last occupied indices in Pitzer ordering for each symblk */
int *occ_alpha, *occ_beta;                  /* Pitzer to "full"(no frozen core) QTS index mapping */
int *act2fullQTS_alpha, *act2fullQTS_beta;  /* Maps "active"(taking frozen core into account) QTS into "frozen" QTS index */
int *ioff3;                                 /* returns pointers to the beginning of rows in a rectangular matrix */

void ump2r12_energy_bb()
{
  pthread_attr_t thread_attr;
  pthread_t *thread_id;

  Libr12_t Libr12;
  long int libr12_memory;
  int max_bf_per_shell;
  int max_num_prim_comb;
  int te_type;
  long int i;

  int num_ibatch, num_i_per_ibatch, ibatch, ibatch_first, ibatch_length;
  int imin, imax, jmin;

  /*---------------
    Initialization
   ---------------*/
  init_fjt(BasisSet.max_am*4+1);
  init_libr12_base();
// make_transqt_arrays_uhf(&first, &last, &fstocc_alpha, &fstocc_beta, &lstocc_alpha, &lstocc_beta, 
//                     &occ_alpha, &occ_beta, &act2fullQTS_alpha, &act2fullQTS_beta, &ioff3);
  timer_init();
  RMP2R12_Status.num_arrived = 0;
  fprintf(outfile,"  Performing direct AO->MO integral tranformation for the UHF MP2-R12/A energy\n");
  fprintf(outfile,"  (BB|BB) and (BB|AA) spin cases\n");
  
  /*-------------------------
    Allocate data structures
   -------------------------*/
  max_bf_per_shell = ioff[BasisSet.max_am];
  /*--- Use this dirty trick to get how much memory integrals library needs ---*/
  max_num_prim_comb = (BasisSet.max_num_prims*
                       BasisSet.max_num_prims)*
                      (BasisSet.max_num_prims*
                       BasisSet.max_num_prims);
  libr12_memory = libr12_storage_required(BasisSet.max_am-1,max_num_prim_comb);
  UserOptions.memory -= libr12_memory*UserOptions.num_threads;

  /*---
    Minimum number of I-batches - 
    take sizes of rsiq_buf, rsix_buf, jsix_buf, and
    jyix_buf into account
    plus enough for the (BB|AA) spin case over the
    [r12,T1] operator only ACS (01/06)
   ---*/
  num_i_per_ibatch = UserOptions.memory / ((NUM_TE_TYPES-1)*
					   (UserOptions.num_threads*(BasisSet.num_ao*max_bf_per_shell*max_bf_per_shell +
								    MOInfo.num_mo*max_bf_per_shell*max_bf_per_shell) +
					    MOInfo.num_mo*MOInfo.beta_act_occ*BasisSet.num_ao +
					    MOInfo.num_mo*MOInfo.beta_act_occ*MOInfo.num_mo) +
                                            UserOptions.num_threads*(MOInfo.num_mo*MOInfo.alpha_act_occ*BasisSet.num_ao +
                                            MOInfo.num_mo*MOInfo.alpha_act_occ*MOInfo.num_mo) );
  if (num_i_per_ibatch > MOInfo.beta_act_occ)
    num_i_per_ibatch = MOInfo.beta_act_occ;
  if (num_i_per_ibatch < 1)
    throw std::domain_error("Not enough memory for direct MP2-R12/A transformation");
  num_ibatch = (MOInfo.beta_act_occ + num_i_per_ibatch - 1) / num_i_per_ibatch;
  /*--- Recompute number of MOs per I-batch ---*/
  num_i_per_ibatch = (MOInfo.beta_act_occ + num_ibatch - 1) / num_ibatch;
  RMP2R12_Status.num_ibatch = num_ibatch;
  RMP2R12_Status.num_i_per_ibatch = num_i_per_ibatch;
  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
    jsix_buf[te_type] = init_array(MOInfo.beta_act_occ*BasisSet.num_ao*
				   num_i_per_ibatch*MOInfo.num_mo);
    jyix_buf[te_type] = init_array(MOInfo.beta_act_occ*MOInfo.num_mo*
				   num_i_per_ibatch*MOInfo.num_mo);
  }
  /* (BB|AA) ints are needed for [r12,T1], due to non-hermiticity -- ACS (01/06)*/
  jsix_buf2 = init_array(MOInfo.alpha_act_occ*BasisSet.num_ao*
                	     num_i_per_ibatch*MOInfo.num_mo);
  jyix_buf2 = init_array(MOInfo.alpha_act_occ*MOInfo.num_mo*
			     num_i_per_ibatch*MOInfo.num_mo);

  fprintf(outfile,"  Using total of %d %s\n",num_ibatch, (num_ibatch == 1) ? "pass" : "passes");
  if (UserOptions.restart) {
    fprintf(outfile,"  (Re)starting at pass %d\n\n",UserOptions.restart_task);
    RMP2R12_Status.ibatch_first = UserOptions.restart_task;
  }
  else {
    fprintf(outfile,"\n");
    RMP2R12_Status.ibatch_first = 0;
  }

  /*--------------------------
    Start compute threads now
   --------------------------*/
  thread_id = (pthread_t *) malloc(UserOptions.num_threads*sizeof(pthread_t));
  pthread_attr_init(&thread_attr);
  pthread_attr_setscope(&thread_attr,
			PTHREAD_SCOPE_SYSTEM);
  pthread_mutex_init(&rmp2r12_energy_mutex,NULL);
  pthread_cond_init(&rmp2r12_energy_cond,NULL);
#if LOCK_RS_SHELL
  rmp2r12_sindex_mutex = (pthread_mutex_t *) malloc(ioff[Symmetry.num_unique_shells]*sizeof(pthread_mutex_t));
  for(i=0;i<ioff[Symmetry.num_unique_shells];i++)
#else
  rmp2r12_sindex_mutex = (pthread_mutex_t *) malloc(BasisSet.num_ao*sizeof(pthread_mutex_t));
  for(i=0;i<BasisSet.num_ao;i++)
#endif
      pthread_mutex_init(&(rmp2r12_sindex_mutex[i]),NULL);
  for(i=0;i<UserOptions.num_threads-1;i++)
    pthread_create(&(thread_id[i]),&thread_attr,
		   ump2r12_energy_thread_bb,(void *)i);
  ump2r12_energy_thread_bb( (void *) (UserOptions.num_threads - 1) );
  for(i=0;i<UserOptions.num_threads-1;i++)
    pthread_join(thread_id[i], NULL);
  free(thread_id);
  pthread_mutex_destroy(&rmp2r12_energy_mutex);
#if LOCK_RS_SHELL
  for(i=0;i<ioff[Symmetry.num_unique_shells];i++)
#else
  for(i=0;i<BasisSet.num_ao;i++)
#endif
      pthread_mutex_destroy(&(rmp2r12_sindex_mutex[i]));

  
  /*---------
    Clean-up
   ---------*/
  free(rmp2r12_sindex_mutex);
  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
    free(jyix_buf[te_type]);
    free(jsix_buf[te_type]);
  }
  free(jyix_buf2);
  free(jsix_buf2);
  free_fjt();
  timer_done();

  return;
}

} /* End ump2r12_bb namespace */

void make_transqt_arrays_uhf(int **first, int **last, int **fstocc_alpha, int **fstocc_beta, int **lstocc_alpha, 
                         int **lstocc_beta, int **occ_alpha, int **occ_beta, int **act2fullQTS_alpha,
                         int **act2fullQTS_beta, int **ioff3)
{
  int h, i, offset, count, irrep;
  int first_offset, last_offset;
  int p,q,row,col,mo,mo_max,qts_offset,pitz_offset;
  double **scf_evec_so;



  /*
     Construct first and last index arrays: this defines the first
     absolute orbital index (Pitzer ordering) and last absolute orbital
     index for each irrep.  When there are no orbitals for an irrep, the
     value is -1 for first[] and -2 for last[].  Note that there must be
     orbitals in the first irrep (i.e. totally symmetric) for this to work.
  */
  *first = init_int_array(Symmetry.nirreps);
  *last = init_int_array(Symmetry.nirreps);
  for(h=0; h < Symmetry.nirreps; h++) {
      (*first)[h] = -1;
      (*last)[h] = -2;
  }
  first_offset = 0;
  last_offset = MOInfo.orbspi[0] - 1; 
  (*first)[0] = first_offset;
  (*last)[0] = last_offset;
  for(h=1; h < Symmetry.nirreps; h++) {
      first_offset += MOInfo.orbspi[h-1];
      last_offset += MOInfo.orbspi[h];
      if(MOInfo.orbspi[h]) {
          (*first)[h] = first_offset;
          (*last)[h] = last_offset;
        }
    }
  
  /* fstocc[] and lstocc[] supply the first and last orbital indices (Pitzer
     ordering) for the occupied orbitals in each irrep. */
  /* The alpha orbitals */
  *fstocc_alpha = init_int_array(Symmetry.nirreps);
  *lstocc_alpha = init_int_array(Symmetry.nirreps);
  for(h=0; h < Symmetry.nirreps; h++) {
      (*fstocc_alpha)[h] = -1;
      (*lstocc_alpha)[h] = -2;
    }
  first_offset = 0;
  last_offset = MOInfo.clsdpi[0] + MOInfo.openpi[0] - 1;
  (*fstocc_alpha)[0] = first_offset;
  (*lstocc_alpha)[0] = last_offset;
  for(h=1; h < Symmetry.nirreps; h++) {
      first_offset += MOInfo.orbspi[h-1];
      last_offset += MOInfo.virtpi_alpha[h-1]+MOInfo.orbspi[h]-MOInfo.virtpi_alpha[h];
      if(MOInfo.clsdpi[h] + MOInfo.openpi[h]) {
          (*fstocc_alpha)[h] = first_offset;
          (*lstocc_alpha)[h] = last_offset;
        }
    }

  /* The beta orbitals */
  *fstocc_beta = init_int_array(Symmetry.nirreps);
  *lstocc_beta = init_int_array(Symmetry.nirreps);
  for(h=0; h < Symmetry.nirreps; h++) {
      (*fstocc_beta)[h] = -1;
      (*lstocc_beta)[h] = -2;
    }
  first_offset = 0;
  last_offset = MOInfo.clsdpi[0] + MOInfo.openpi[0] - 1;
  (*fstocc_beta)[0] = first_offset;
  (*lstocc_beta)[0] = last_offset;
  for(h=1; h < Symmetry.nirreps; h++) {
      first_offset += MOInfo.orbspi[h-1];
      last_offset += MOInfo.virtpi_beta[h-1]+MOInfo.orbspi[h]-MOInfo.virtpi_beta[h];
      if(MOInfo.clsdpi[h] + MOInfo.openpi[h]) {
          (*fstocc_beta)[h] = first_offset;
          (*lstocc_beta)[h] = last_offset;
        }
    }
  /* Construct occupied Pitzer -> QTS ordering arrays for
   occupied (occ[]) orbitals */
  /* The alpha orbitals */
  *occ_alpha = init_int_array(MOInfo.num_mo);
  for(i=0; i< MOInfo.num_mo; i++) {
      (*occ_alpha)[i] = -1;
    }
  
  offset = 0;
  count=0;
  for(h=0; h < Symmetry.nirreps; h++) {
      if(h)
          offset += MOInfo.orbspi[h-1];
      for(i=offset; i < (offset+MOInfo.clsdpi[h]+MOInfo.openpi[h]); i++) {
          (*occ_alpha)[i] = count++;
        }
    }

  /* The beta orbitals */
  *occ_beta = init_int_array(MOInfo.num_mo);
  for(i=0; i< MOInfo.num_mo; i++) {
      (*occ_beta)[i] = -1;
    }
  
  offset = 0;
  count=0;
  for(h=0; h < Symmetry.nirreps; h++) {
      if(h)
          offset += MOInfo.orbspi[h-1];
      for(i=offset; i < (offset+MOInfo.clsdpi[h]); i++) {
          (*occ_beta)[i] = count++;
        }
    }

  /* Construct active -> full QTS ordering array for occupied alpha orbitals */
  *act2fullQTS_alpha = init_int_array(MOInfo.ndocc+MOInfo.nsocc);
  offset = 0;
  count=0;
  for(h=0; h < Symmetry.nirreps; h++) {
    for(i=0; i < MOInfo.frozen_docc[h]; i++) {
      (*act2fullQTS_alpha)[count] = offset+i;
      count++;
    }
    offset += MOInfo.clsdpi[h] + MOInfo.openpi[h];
  }
  offset = 0;
  for(h=0; h < Symmetry.nirreps; h++) {
    for(i=MOInfo.frozen_docc[h]; i < MOInfo.clsdpi[h] + MOInfo.openpi[h]; i++) {
      (*act2fullQTS_alpha)[count] = offset+i;
      count++;
    }
    offset += MOInfo.clsdpi[h] + MOInfo.openpi[h];
  }
 
  /* Construct active -> full QTS ordering array for occupied beta orbitals */
  *act2fullQTS_beta = init_int_array(MOInfo.ndocc);
  offset = 0;
  count=0;
  for(h=0; h < Symmetry.nirreps; h++) {
    for(i=0; i < MOInfo.frozen_docc[h]; i++) {
      (*act2fullQTS_beta)[count] = offset+i;
      count++;
    }
    offset += MOInfo.clsdpi[h];
  }
  offset = 0;
  for(h=0; h < Symmetry.nirreps; h++) {
    for(i=MOInfo.frozen_docc[h]; i < MOInfo.clsdpi[h]; i++) {
      (*act2fullQTS_beta)[count] = offset+i;
      count++;
    }
    offset += MOInfo.clsdpi[h];
  }
  
  
  /* Generate ioff3 array.  This array gives the row offset for an
     ndocc x nmo matrix */
  *ioff3 = init_int_array(MOInfo.ndocc);
  for(i=0; i < MOInfo.ndocc; i++) {
      (*ioff3)[i] = i*MOInfo.num_mo;
    }

  return;
}


}} /* End psi::cints namespaces */
