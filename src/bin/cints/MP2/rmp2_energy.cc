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
#include<libciomr/libciomr.h>
#include<libchkpt/chkpt.h>
#include<libqt/qt.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"schwartz.h"
#include"quartet_data.h"
#include"norm_quartet.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
  #include"fjt.h"
#endif
#include"quartet_permutations.h"
#include"rmp2_energy.h"

namespace psi { 
  namespace CINTS {

    /*-------------------------------
      Explicit function declarations
      -------------------------------*/
    extern void *rmp2_energy_thread(void *);
    
    /*--------------------------------
      Variables common to all threads
      --------------------------------*/
    pthread_mutex_t rmp2_energy_mutex;
    pthread_mutex_t *rmp2_sindex_mutex;
    pthread_cond_t rmp2_energy_cond;
    RMP2_Status_t RMP2_Status;
    double *jsia_buf;             /* buffer for (js|ia) integrals, where j runs over all d.-o. MOs,
				     s runs over all AOs, i - over I-batch, a - over all virtuals */
    double *jbia_buf;             /* buffer contains all MP2-type integrals */
    
    
    /*!-------------------------------------------------------
      Algorithm (for the original one see Ida Nielsen's Ph.D. thesis)
      
      ***Split into threads. Each thread do the following:
      
      Loop over I batches (batch size num_i_per_batch)
      
      Loop over all symmetry-unique shells UR, US<=UR
      ***if this UR, US is to be handled by this thread - proceed, else skip to next one
      Find all symmetry-distinct shell doublets resulting from (UR US|
      Loop over all resulting shell doublets R, S
      
      Loop over all symmetry-unique shells UP, UQ<=UP
      Find all symmetry-distinct shell quartets resulting from (R S|UP UQ)
	  Loop over the resulting set of P, Q doublets
            Evaluate (RS|PQ)
	    Loop over p in P, q in Q, r in R, s in S, i in I
	      (rs|iq) += Cpi * (rs|pq)
	      (rs|ip) += Cqi * (rs|pq)
	    End p, q, r, s, i loop
          End P,Q loop
        End UP, UQ loop

        Loop over r in R, s in S
          Loop over q < nao, a < nuocc, i in I
	    (rs|ia) += Cas * (rs|is)
	  End q, a, i
        End r, s loop
	***Lock (js| and (jr| (either individual or shell blocks depending on LOCK_RS_SHELL in rmp2_energy.h)
	Loop over r in R, s in S
	  Loop over i in I, a < nuocc, j <= i
	    (js|ia) += Cjr * (rs|ia)
	    (jr|ia) += Cjs * (rs|ia)
	  End i, a, j loop
        End r, s loop
	***Unlock (js| and (jr|

      End R, S loop
    End UR, US loop

    ***Barrier: threads wait until everyone is done
    ***Do the following in one thread only
    Loop over i in I, j <= i
      Loop over r < nao, a < nuocc, b < nuocc
        (jb|ia) += Cbs * (js|ia)
      End r, a, b loop
    End i, j loop

  End I loop

  ***Merge all threads
  
  -------------------------------------------------------*/
    void rmp2_energy()
    {
      pthread_attr_t thread_attr;
      pthread_t *thread_id;
  
      /*--- Various data structures ---*/
      Libint_t Libint;
      long int libint_memory;
      int max_bf_per_shell;
      int max_num_prim_comb;
      
      int i;
      
      
      int num_ibatch, num_i_per_ibatch, ibatch, ibatch_length;
      int imin, imax, jmin;
      int mo_i, mo_j, mo_a, mo_b, mo_ij;
      int rs_offset, rsi_offset, rsp_offset;
      
      double AB2, CD2;
      
      double *raw_data;             /* pointer to the unnormalized taregt quartet of integrals */
      double *data;                 /* pointer to the transformed normalized target quartet of integrals */
#ifdef NONDOUBLE_INTS
      REALTYPE *target_ints;            /* Pointer to the location of the target quartet on the stack of
					   integrals quartets if libint.a is using other than regular doubles */
#endif
      
      double *rspq_ptr;
      double temp;
      double *mo_vec;
      double *rsiq_buf;             /* buffer for (rs|iq) integrals, where r,s run over shell sets,
				       i runs over I-batch, q runs over all AOs */
      double *rsi_row;
      double **ia_buf;              /* buffer for one |ia) ket */
      double *i_row;
      double *jsi_row;
      double *jbi_row;
      double iajb, ibja, pfac;
      
      double temp1,temp2,*iq_row,*ip_row;
      int rs,qrs;
      
      /*---------------
	Initialization
	---------------*/
#ifdef USE_TAYLOR_FM
      init_Taylor_Fm_Eval(BasisSet.max_am*4-4,UserOptions.cutoff);
#else
      init_fjt(BasisSet.max_am*4);
#endif
      init_libint_base();
      timer_init();
      timer_on("Schwartz");
      schwartz_eri();
      timer_off("Schwartz");
      fprintf(outfile,"  Computing RHF MP2 energy via direct algorithm\n");
      RMP2_Status.Emp2_0 = 0.0;
      RMP2_Status.Emp2_1 = 0.0;
      RMP2_Status.num_arrived = 0;
      
      
      /*-------------------------
	Allocate data structures
	-------------------------*/
      max_bf_per_shell = ioff[BasisSet.max_am];
      /*--- Use this dirty trick to get how much memory integrals library needs ---*/
      max_num_prim_comb = (BasisSet.max_num_prims*
			   BasisSet.max_num_prims)*
	(BasisSet.max_num_prims*
	 BasisSet.max_num_prims);
      libint_memory = libint_storage_required(BasisSet.max_am-1,max_num_prim_comb);
      UserOptions.memory -= libint_memory*UserOptions.num_threads;
      
      /*---
	Minimum number of I-batches - 
	take sizes of rsiq_buf, rsia_buf, jsia_buf, and
	jbia_buf into account
   ---*/
      num_i_per_ibatch = UserOptions.memory / (UserOptions.num_threads*(BasisSet.num_ao*max_bf_per_shell*max_bf_per_shell +
									MOInfo.nactuocc*max_bf_per_shell*max_bf_per_shell) +
					       MOInfo.nactuocc*MOInfo.nactdocc*BasisSet.num_ao +
					       MOInfo.nactuocc*MOInfo.nactdocc*MOInfo.nactuocc);
      if (num_i_per_ibatch > MOInfo.nactdocc)
	num_i_per_ibatch = MOInfo.nactdocc;
      if (num_i_per_ibatch < 1)
	throw std::domain_error("Not enough memory for direct MP2");
      num_ibatch = (MOInfo.nactdocc + num_i_per_ibatch - 1) / num_i_per_ibatch;
      /*--- Recompute number of MOs per I-batch ---*/
      num_i_per_ibatch = (MOInfo.nactdocc + num_ibatch - 1) / num_ibatch;
      RMP2_Status.num_ibatch = num_ibatch;
      RMP2_Status.num_i_per_ibatch = num_i_per_ibatch;
      RMP2_Status.emp2_0 = init_array(ioff[MOInfo.ndocc]);
      RMP2_Status.emp2_1 = init_array(ioff[MOInfo.ndocc]);
      jsia_buf = init_array(MOInfo.nactdocc*BasisSet.num_ao*
			    num_i_per_ibatch*MOInfo.nactuocc);
      jbia_buf = init_array(MOInfo.nactdocc*MOInfo.nactuocc*
			    num_i_per_ibatch*MOInfo.nactuocc);
      fprintf(outfile,"  Using %d %s\n\n",num_ibatch, (num_ibatch == 1) ? "pass" : "passes");
      
      /*--------------------------
	Start compute threads now
	--------------------------*/
      thread_id = (pthread_t *) malloc(UserOptions.num_threads*sizeof(pthread_t));
      pthread_attr_init(&thread_attr);
      pthread_attr_setscope(&thread_attr,
			    PTHREAD_SCOPE_SYSTEM);
      pthread_mutex_init(&rmp2_energy_mutex,NULL);
      pthread_cond_init(&rmp2_energy_cond,NULL);
#if LOCK_RS_SHELL
      rmp2_sindex_mutex = (pthread_mutex_t *) malloc(ioff[BasisSet.num_shells]*sizeof(pthread_mutex_t));
      for(i=0;i<ioff[BasisSet.num_shells];i++)
#else
	rmp2_sindex_mutex = (pthread_mutex_t *) malloc(BasisSet.num_ao*sizeof(pthread_mutex_t));
      for(i=0;i<BasisSet.num_ao;i++)
#endif
	pthread_mutex_init(&(rmp2_sindex_mutex[i]),NULL);
    for(long int i=0;i<UserOptions.num_threads-1;i++)
      pthread_create(&(thread_id[i]),&thread_attr,
		       rmp2_energy_thread,(void *)i);
    rmp2_energy_thread( (void *) (UserOptions.num_threads - 1) );
    for(i=0;i<UserOptions.num_threads-1;i++)
      pthread_join(thread_id[i], NULL);
    free(thread_id);
    pthread_mutex_destroy(&rmp2_energy_mutex);
#if LOCK_RS_SHELL
      for(i=0;i<ioff[BasisSet.num_shells];i++)
#else
      for(i=0;i<BasisSet.num_ao;i++)
#endif
	pthread_mutex_destroy(&(rmp2_sindex_mutex[i]));
      
      fprintf(outfile,"\n");
      fprintf(outfile,"  Singlet pair energies:\n");
      fprintf(outfile,"    i       j         e(ij)\n");
      fprintf(outfile,"  -----   -----   ------------\n");
      for(mo_i=MOInfo.nfrdocc;mo_i<MOInfo.ndocc;mo_i++)
	for(mo_j=MOInfo.nfrdocc;mo_j<=mo_i;mo_j++) {
	  mo_ij = INDEX(mo_i, mo_j);
	  fprintf(outfile,"  %3d     %3d     %12.9lf\n",mo_i+1,mo_j+1,RMP2_Status.emp2_0[mo_ij]);
	}
      fprintf(outfile,"\n");
      fprintf(outfile,"  Triplet pair energies:\n");
      fprintf(outfile,"    i       j         e(ij)\n");
      fprintf(outfile,"  -----   -----   ------------\n");
      for(mo_i=MOInfo.nfrdocc;mo_i<MOInfo.ndocc;mo_i++)
	for(mo_j=MOInfo.nfrdocc;mo_j<mo_i;mo_j++) {
	  mo_ij = INDEX(mo_i, mo_j);
	  fprintf(outfile,"  %3d     %3d     %12.9lf\n",mo_i+1,mo_j+1,RMP2_Status.emp2_1[mo_ij]);
	}
      fprintf(outfile,"\n");
      
      fprintf(outfile,"  MBPT[2] Correlation Energy (singlet) = %20.10lf\n",RMP2_Status.Emp2_0);
      fprintf(outfile,"  MBPT[2] Correlation Energy (triplet) = %20.10lf\n",RMP2_Status.Emp2_1);
      fprintf(outfile,"  MBPT[2] Correlation Energy (total)   = %20.10lf\n",RMP2_Status.Emp2_0+RMP2_Status.Emp2_1);
      fprintf(outfile,"  Total MBPT[2] Energy                 = %20.10lf\n\n",
	      RMP2_Status.Emp2_0+RMP2_Status.Emp2_1+MOInfo.Escf);
      
      chkpt_wt_etot(RMP2_Status.Emp2_0+RMP2_Status.Emp2_1+MOInfo.Escf);
      
      /*---------
	Clean-up
	---------*/
      free(rmp2_sindex_mutex);
      free(jbia_buf);
      free(jsia_buf);
#ifdef USE_TAYLOR_FM
      free_Taylor_Fm_Eval();
#else
      free_fjt();
#endif
      
      timer_done();
      
      return;
    }
  }
}

