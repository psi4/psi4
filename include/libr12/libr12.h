#ifndef _psi3_libr12_h
#define _psi3_libr12_h

#include <libint/libint.h>
/* Maximum angular momentum of functions in a basis set plus 1 */
#define LIBR12_MAX_AM 4
#define LIBR12_OPT_AM 4
#define NUM_TE_TYPES 4

typedef struct {
  REALTYPE AB[3];
  REALTYPE CD[3];
  REALTYPE AC[3];
  REALTYPE ABdotAC, CDdotCA;
  } contr_data;

typedef struct {
  REALTYPE *int_stack;
  prim_data *PrimQuartet;
  contr_data ShellQuartet;
  REALTYPE *te_ptr[NUM_TE_TYPES];
  REALTYPE *t1vrr_classes[7][7];
  REALTYPE *t2vrr_classes[7][7];
  REALTYPE *rvrr_classes[7][7];
  REALTYPE *gvrr_classes[8][8];
  REALTYPE *r12vrr_stack;

  } Libr12_t;

#ifdef __cplusplus
extern "C" {
#endif
extern void (*build_r12_gr[4][4][4][4])(Libr12_t *, int);
extern void (*build_r12_grt[4][4][4][4])(Libr12_t *, int);
void init_libr12_base();

int  init_libr12(Libr12_t *, int max_am, int max_num_prim_quartets);
void free_libr12(Libr12_t *);
int  libr12_storage_required(int max_am, int max_num_prim_quartets);

#ifdef __cplusplus
}
#endif

#endif
