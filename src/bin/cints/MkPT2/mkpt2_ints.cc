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
#include<libiwl/iwl.h>
#include<libint/libint.h>

#include"defines.h"
#include"psifiles.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"schwartz.h"
#include"quartet_data.h"
#include"norm_quartet.h"
#include "data_structs.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
  #include"fjt.h"
#endif
#include"quartet_permutations.h"
#include"mkpt2_ints.h"


namespace psi { 
  namespace cints {
    namespace mkpt2 {
    
    /*-------------------------------
      Explicit function declarations
      -------------------------------*/
    extern void *mkpt2_ints_thread(void *);
    
    /*--------------------------------
      Varixbles common to all threads
      --------------------------------*/
    pthread_mutex_t mkpt2_energy_mutex;
    pthread_mutex_t *mkpt2_sindex_mutex;
    pthread_cond_t mkpt2_energy_cond;
    MkPT2_Status_t MkPT2_Status;
    double *jsix_buf;
    double *jyix_buf;
    double *asij_buf;
    double *abij_buf;
    unsigned long int n_ij;
    unsigned long int n_xy;
    unsigned long int n_ab;
#if MkPT2_USE_IWL
    tebuf *iwl_buf;
    iwlbuf ERIOUT;
    int iwl_count;
#else
    double * xy_buf;
#endif
    
    
    void mkpt2_ints()
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
      double **ix_buf;              /* buffer for one |ix) ket */
      double *i_row;
      double *jsi_row;
      double *jbi_row;
      double ixjy, iyjx, pfac;
      
      double temp1,temp2,*iq_row,*ip_row;
      int rs,qrs;
      
      /*---------------
	Initixlization
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
      MkPT2_Status.num_arrived = 0;
      
      
      /*-------------------------
	Allocate data structures
	-------------------------*/
      n_ij = MOInfo.ndocc*(MOInfo.ndocc+1)/2;
      n_ab = MOInfo.nuocc*(MOInfo.nuocc+1)/2;
      n_xy = MOInfo.num_mo*(MOInfo.num_mo+1)/2;
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
	take sizes of rsiq_buf, rsix_buf, jsix_buf,
	jyix_buf ,abij_buf, asij_buf, rsij_buf, xy_buf into account
      ---*/
      fprintf(outfile,"\n  Computing MkPT2 integrals\n");
      num_i_per_ibatch = UserOptions.memory / (UserOptions.num_threads*
                              (
                                BasisSet.num_ao*max_bf_per_shell*max_bf_per_shell + /*rsiq*/
                                MOInfo.num_mo*max_bf_per_shell*max_bf_per_shell + /*rsix*/
                                MOInfo.ndocc*max_bf_per_shell*max_bf_per_shell /*rsij*/
                              ) +
			      MOInfo.num_mo*MOInfo.ndocc*BasisSet.num_ao + /*jsix*/
			      MOInfo.nuocc*MOInfo.ndocc*BasisSet.num_ao + /*asij*/
			      n_ab*MOInfo.ndocc + /*abij*/
			      MOInfo.num_mo*MOInfo.ndocc*MOInfo.num_mo + /*jyix*/
#if MkPT2_USE_IWL
                              MOInfo.num_mo*MOInfo.num_mo*sizeof(struct tebuf)/sizeof(double) /*iwl_buf*/
#else
                              MOInfo.num_mo*MOInfo.num_mo /*xy_buf*/
#endif
                             );
      if (num_i_per_ibatch > MOInfo.ndocc)
	num_i_per_ibatch = MOInfo.ndocc;
      if (num_i_per_ibatch < 1)
	throw std::domain_error("Not enough memory for direct MkPT2 integrals");
      num_ibatch = (MOInfo.ndocc + num_i_per_ibatch - 1) / num_i_per_ibatch;
      /*--- Recompute number of MOs per I-batch ---*/
      num_i_per_ibatch = (MOInfo.ndocc + num_ibatch - 1) / num_ibatch;
      MkPT2_Status.num_ibatch = num_ibatch;
      MkPT2_Status.num_i_per_ibatch = num_i_per_ibatch;
      jsix_buf = init_array(MOInfo.ndocc*BasisSet.num_ao* num_i_per_ibatch*MOInfo.num_mo);
      jyix_buf = init_array(MOInfo.ndocc*MOInfo.num_mo* num_i_per_ibatch*MOInfo.num_mo);
      asij_buf = init_array(MOInfo.nuocc*BasisSet.num_ao*num_i_per_ibatch*MOInfo.ndocc);
      abij_buf = init_array(n_ab* num_i_per_ibatch*MOInfo.ndocc);
#if MkPT2_USE_IWL
      iwl_buf = (struct tebuf*) malloc(MOInfo.num_mo*MOInfo.num_mo*sizeof(struct tebuf));
#else
      xy_buf = init_array(MOInfo.num_mo*MOInfo.num_mo);
#endif
      fprintf(outfile,"  Using %d %s\n\n",num_ibatch, (num_ibatch == 1) ? "pass" : "passes");
      
      /*--------------------------
	Start compute threads now
	--------------------------*/
      thread_id = (pthread_t *) malloc(UserOptions.num_threads*sizeof(pthread_t));
      pthread_attr_init(&thread_attr);
      pthread_attr_setscope(&thread_attr,
			    PTHREAD_SCOPE_SYSTEM);
      pthread_mutex_init(&mkpt2_energy_mutex,NULL);
      pthread_cond_init(&mkpt2_energy_cond,NULL);
#if LOCK_RS_SHELL
      mkpt2_sindex_mutex = (pthread_mutex_t *) malloc(ioff[BasisSet.num_shells]*sizeof(pthread_mutex_t));
      for(i=0;i<ioff[BasisSet.num_shells];i++)
#else
	mkpt2_sindex_mutex = (pthread_mutex_t *) malloc(BasisSet.num_ao*sizeof(pthread_mutex_t));
      for(i=0;i<BasisSet.num_ao;i++)
#endif
	pthread_mutex_init(&(mkpt2_sindex_mutex[i]),NULL);
    for(long int i=0;i<UserOptions.num_threads-1;i++)
      pthread_create(&(thread_id[i]),&thread_attr,
		       mkpt2_ints_thread,(void *)i);
    mkpt2_ints_thread( (void *) (UserOptions.num_threads - 1) );
    for(i=0;i<UserOptions.num_threads-1;i++)
      pthread_join(thread_id[i], NULL);
    free(thread_id);
    pthread_mutex_destroy(&mkpt2_energy_mutex);
#if LOCK_RS_SHELL
      for(i=0;i<ioff[BasisSet.num_shells];i++)
#else
      for(i=0;i<BasisSet.num_ao;i++)
#endif
	pthread_mutex_destroy(&(mkpt2_sindex_mutex[i]));
      
#if MkPT2_TEST
     double *temp_arr = init_array(MOInfo.num_mo*MOInfo.num_mo*MOInfo.num_mo*MOInfo.num_mo); 
     iwl_buf_init(&ERIOUT,PSIF_MO_TEI,UserOptions.cutoff,1,1);
     iwl_buf_rd_all(&ERIOUT,temp_arr,ioff,ioff,1,ioff,1,outfile);
     iwl_buf_close(&ERIOUT, 1);  
     free(temp_arr);
#endif

      /*---------
	Clean-up
	---------*/
      free(mkpt2_sindex_mutex);
#if MkPT2_USE_IWL
      free(iwl_buf);
#else
      free(xy_buf);
#endif
      free(jyix_buf);
      free(abij_buf);
      free(jsix_buf);
      free(asij_buf);
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
}

