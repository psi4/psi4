#ifndef _psi3_libderiv_h
#define _psi3_libderiv_h

#include <libint/libint.h>

/* Maximum angular momentum of functions in a basis set plus 1 */
#define LIBDERIV_MAX_AM1 5
#define LIBDERIV_MAX_AM12 4
#ifdef DERIV_LVL
 #undef DERIV_LVL
#endif
#define DERIV_LVL 2

typedef struct {
  double *int_stack;
  prim_data *PrimQuartet;
  double *zero_stack;
  double *ABCD[12+144];
  double AB[3];
  double CD[3];
  double *deriv_classes[9][9][12];
  double *deriv2_classes[9][9][144];
  double *dvrr_classes[9][9];
  double *dvrr_stack;
  } Libderiv_t;

#ifdef __cplusplus
extern "C" {
#endif
extern void (*build_deriv1_eri[5][5][5][5])(Libderiv_t *, int);
extern void (*build_deriv12_eri[4][4][4][4])(Libderiv_t *, int);
void init_libderiv_base();

int  init_libderiv1(Libderiv_t *, int max_am, int max_num_prim_quartets, int max_cart_class_size);
int  init_libderiv12(Libderiv_t *, int max_am, int max_num_prim_quartets, int max_cart_class_size);
void free_libderiv(Libderiv_t *);

int  libderiv1_storage_required(int max_am, int max_num_prim_quartets, int max_cart_class_size);
int  libderiv12_storage_required(int max_am, int max_num_prim_quartets, int max_cart_class_size);
#ifdef __cplusplus
}
#endif

#endif
