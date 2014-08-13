#ifndef _GLOBALP_H_
#define _GLOBALP_H_

#include <stdio.h>

#include "gaconfig.h"

#ifdef __crayx1
#undef CRAY
#endif

#ifdef FALSE
#undef FALSE
#endif
#ifdef TRUE
#undef TRUE
#endif
#ifdef CRAY_YMP
#define FALSE _btol(0)
#define TRUE  _btol(1)
#else
#define FALSE (logical) 0
#define TRUE  (logical) 1
#endif

#if HAVE_WINDOWS_H
#   include <windows.h>
#   define sleep(x) Sleep(1000*(x))
#endif
#include "macdecls.h"

#define GA_OFFSET   1000           /* offset for handle numbering */

/* types/tags of messages used internally by GA */
#define     GA_TYPE_SYN   GA_MSG_OFFSET + 1
#define     GA_TYPE_GSM   GA_MSG_OFFSET + 5
#define     GA_TYPE_GOP   GA_MSG_OFFSET + 15
#define     GA_TYPE_BRD   GA_MSG_OFFSET + 16

/* GA operation ids */
#define     GA_OP_GET 1          /* Get                         */
#define     GA_OP_END 2          /* Terminate                   */
#define     GA_OP_CRE 3          /* Create                      */
#define     GA_OP_PUT 4          /* Put                         */
#define     GA_OP_ACC 5          /* Accumulate                  */
#define     GA_OP_DES 6          /* Destroy                     */
#define     GA_OP_DUP 7          /* Duplicate                   */
#define     GA_OP_ZER 8          /* Zero                        */
#define     GA_OP_DDT 9          /* dot product                 */
#define     GA_OP_SCT 10         /* scatter                     */
#define     GA_OP_GAT 11         /* gather                      */
#define     GA_OP_RDI 15         /* Integer read and increment  */
#define     GA_OP_ACK 16         /* acknowledgment              */
#define     GA_OP_LCK 17         /* acquire lock                */
#define     GA_OP_UNL 18         /* release lock                */


#ifdef ENABLE_TRACE
  static Integer     op_code;
#endif


#define GA_MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define GA_MIN(a,b) (((a) <= (b)) ? (a) : (b))
#define GA_ABS(a)   (((a) >= 0) ? (a) : (-(a)))

typedef struct ga_typeinfo_t {
  int active;
  size_t size;
} ga_typeinfo_t;

extern ga_typeinfo_t ga_types[];

#define GA_TYPES_MAX 256
#define GA_TYPES_RESERVED 17 /**Should match num lines initialized in ga_types struct*/

#define GAtypebuiltinM(_type) ((_type)>=MT_BASE && (_type)<(MT_BASE+GA_TYPES_RESERVED))
#define GAsizeofM(_type)   ga_types[(_type)-MT_BASE].size
#define GAvalidtypeM(_type) ((_type)>=MT_BASE && (_type)<(MT_BASE+GA_TYPES_MAX) && ga_types[(_type)-MT_BASE].active!=0)

#define NAME_STACK_LEN 10
#define PAGE_SIZE  4096

struct ga_stat_t {
         long   numcre; 
         long   numdes;
         long   numget;
         long   numput;
         long   numacc;
         long   numsca;
         long   numgat;
         long   numrdi;
         long   numser;
         long   curmem; 
         long   maxmem; 
         long   numget_procs;
         long   numput_procs;
         long   numacc_procs;
         long   numsca_procs;
         long   numgat_procs;
};

struct ga_bytes_t{ 
         double acctot;
         double accloc;
         double gettot;
         double getloc;
         double puttot;
         double putloc;
         double rditot;
         double rdiloc;
         double gattot;
         double gatloc;
         double scatot;
         double scaloc;
};

#define STAT_AR_SZ sizeof(ga_stat_t)/sizeof(long)

extern long *GAstat_arr;  
extern struct ga_stat_t GAstat;
extern struct ga_bytes_t GAbytes;
extern char *GA_name_stack[NAME_STACK_LEN];    /* stack for names of GA ops */ 
extern int GA_stack_size;
extern int _ga_sync_begin;
extern int _ga_sync_end;
extern int *_ga_argc;
extern char ***_ga_argv;

#define  GA_PUSH_NAME(name) (GA_name_stack[GA_stack_size++] = (name)) 
#define  GA_POP_NAME        (GA_stack_size--)

/* periodic operations */
#define PERIODIC_GET 1
#define PERIODIC_PUT 2
#define PERIODIC_ACC 3

#define FLUSH_CACHE
#ifdef  CRAY_T3D
#       define ALLIGN_SIZE      32
#else
#       define ALLIGN_SIZE      128
#endif

#define allign__(n, SIZE) (((n)%SIZE) ? (n)+SIZE - (n)%SIZE: (n))
#define allign_size(n) allign__((long)(n), ALLIGN_SIZE)
#define allign_page(n) allign__((long)(n), PAGE_SIZE)

extern void    ga_free(void *ptr);
extern void*   ga_malloc(Integer nelem, int type, char *name);
extern void    gai_init_onesided();
extern void    gai_finalize_onesided();
extern void    gai_print_subscript(char *pre,int ndim, Integer subscript[], char* post);
extern Integer GAsizeof(Integer type);
extern void    ga_sort_gath(Integer *pn, Integer *i, Integer *j, Integer *base);
extern void    ga_sort_permutation(Integer *pn, Integer *index, Integer *base);
extern void    ga_sort_scat(Integer *pn, void *v, Integer *i, Integer *j, Integer *base, Integer type);
extern void    gai_hsort(Integer *list, int num);
extern void    ga_init_nbhandle(Integer *nbhandle);
extern int     nga_test_internal(Integer *nbhandle);
extern int     nga_wait_internal(Integer *nbhandle);
extern int     ga_icheckpoint_init(Integer *gas, int num);
extern int     ga_icheckpoint(Integer *gas, int num);
extern int     ga_irecover(int rid);
extern int     ga_icheckpoint_finalize(int g_a);
extern void    ga_checkpoint_arrays(Integer *gas,int num);
extern int     ga_recover_arrays(Integer *gas, int num);
extern void    set_ga_group_is_for_ft(int val);
extern void    ga_set_spare_procs(int *spare);

#endif /* _GLOBALP_H_ */
