/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstring>
#include <cstdio>
#include <memory.h>
#include <cstdlib>
#include <pthread.h>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libint/libint.h>
#include <libderiv/libderiv.h>
#include "defines.h"
#define EXTERN
#include "global.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
#endif
#include "small_fns.h"
#include <Tools/prints.h>

#define USE_SYM 1

namespace psi { namespace CINTS {
extern void *te_deriv1_scf_thread(void *);
extern void *te_deriv1_scf_thread_symm(void *);
extern pthread_mutex_t deriv1_mutex;
extern double **grad_te;
extern void te_deriv1_print(void);

void te_deriv1_scf()
{
  pthread_attr_t thread_attr;
  pthread_t *thread_id;
  
  int i;

  /*---------------
    Initialization
   ---------------*/
#ifdef USE_TAYLOR_FM
  init_Taylor_Fm_Eval(BasisSet.max_am*4-4+DERIV_LVL,UserOptions.cutoff);
#else
  init_fjt(BasisSet.max_am*4+DERIV_LVL);
#endif
  init_libderiv_base();
  grad_te = block_matrix(Molecule.num_atoms,3);

  thread_id = (pthread_t *) malloc(UserOptions.num_threads*sizeof(pthread_t));
  pthread_attr_init(&thread_attr);
  pthread_attr_setscope(&thread_attr,
			PTHREAD_SCOPE_SYSTEM);
  pthread_mutex_init(&deriv1_mutex,NULL);
#if USE_SYM
  for(long int i=0;i<UserOptions.num_threads-1;i++)
    pthread_create(&(thread_id[i]),&thread_attr,
		   te_deriv1_scf_thread_symm,(void *)i);
  te_deriv1_scf_thread_symm( (void *) (UserOptions.num_threads - 1) );
#else
  for(long int i=0;i<UserOptions.num_threads-1;i++)
    pthread_create(&(thread_id[i]),&thread_attr,
		   te_deriv1_scf_thread,(void *)i);
  te_deriv1_scf_thread( (void *) (UserOptions.num_threads - 1) );
#endif

#if PRINT_DERIV1
  te_deriv1_print();
#endif

  for(i=0;i<UserOptions.num_threads-1;i++)
    pthread_join(thread_id[i], NULL);
  free(thread_id);
  
  if (UserOptions.print_lvl >= PRINT_TEDERIV)
    print_atomvec("Two-electron contribution to the forces (a.u.)",grad_te);

  for(i=0;i<Molecule.num_atoms;i++) {
    Grad[i][0] += grad_te[i][0];
    Grad[i][1] += grad_te[i][1];
    Grad[i][2] += grad_te[i][2];
  }
  
  /*---------
    Clean-up
   ---------*/
  free_block(grad_te);
#ifdef USE_TAYLOR_FM
  free_Taylor_Fm_Eval();
#else
  free_fjt();
#endif

  return;
}
}}
