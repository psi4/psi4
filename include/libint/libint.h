#ifndef _psi3_libint_h
#define _psi3_libint_h

/* Maximum angular momentum of functions in a basis set plus 1 */
#define REALTYPE double
#define LIBINT_MAX_AM 5
#define LIBINT_OPT_AM 3
typedef struct pdata{
  REALTYPE F[17];
  REALTYPE U[6][3];
  REALTYPE twozeta_a;
  REALTYPE twozeta_b;
  REALTYPE twozeta_c;
  REALTYPE twozeta_d;
  REALTYPE oo2z;
  REALTYPE oo2n;
  REALTYPE oo2zn;
  REALTYPE poz;
  REALTYPE pon;
  REALTYPE oo2p;
  REALTYPE ss_r12_ss;
  } prim_data;

typedef struct {
  REALTYPE *int_stack;
  prim_data *PrimQuartet;
  REALTYPE AB[3];
  REALTYPE CD[3];
  REALTYPE *vrr_classes[9][9];
  REALTYPE *vrr_stack;
  } Libint_t;

#ifdef __cplusplus
extern "C" {
#endif
extern REALTYPE *(*build_eri[5][5][5][5])(Libint_t *, int);
void init_libint_base();
int  init_libint(Libint_t *, int max_am, int max_num_prim_comb);
void free_libint(Libint_t *);
int  libint_storage_required(int max_am, int max_num_prim_comb);
#ifdef __cplusplus
}
#endif

#endif
